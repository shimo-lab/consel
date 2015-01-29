/*
  opt.c

  Time-stamp: <2011-01-25 16:57:16 shimo>

  shimo@ism.ac.jp 
  Hidetoshi Shimodaira

  Molphy: luinverse (J.Adachi)
  Ref: Numerical Recipes in C (Press et al)

*/

#include <stdio.h>
#include <math.h>
#include "misc.h"

static const char rcsid[] = "$Id: opt.c,v 1.3 2011/05/12 07:21:41 shimo Exp $";

/*
  INVERSION OF MATRIX ON LU DECOMPOSITION
*/
void luinverse(double **omtrx, double **imtrx, int size) /* From Molphy */
{
double eps = 1.0e-20; /* ! */
	int i, j, k, l, maxi, idx, ix, jx;
	double sum, tmp, maxb, aw;
	int *index;
	double *wk;

	index = (int *) MALLOC((unsigned)size * sizeof(int));
	wk = (double *) MALLOC((unsigned)size * sizeof(double));

	aw = 1.0;
	for (i = 0; i < size; i++) {
		maxb = 0.0;
		for (j = 0; j < size; j++) {
			if (fabs(omtrx[i][j]) > maxb)
				maxb = fabs(omtrx[i][j]);
		}
		if (maxb == 0.0) {
		  error("luinverse: singular matrix");
		}
		wk[i] = 1.0 / maxb;
	}
	for (j = 0; j < size; j++) {
		for (i = 0; i < j; i++) {
			sum = omtrx[i][j];
			for (k = 0; k < i; k++)
				sum -= omtrx[i][k] * omtrx[k][j];
			omtrx[i][j] = sum;
		}
		maxb = 0.0;
		maxi=0;
		for (i = j; i < size; i++) {
			sum = omtrx[i][j];
			for (k = 0; k < j; k++)
				sum -= omtrx[i][k] * omtrx[k][j];
			omtrx[i][j] = sum;
			tmp = wk[i] * fabs(sum);
			if (tmp >= maxb) {
				maxb = tmp;
				maxi = i;
			}
		}
		if (j != maxi) {
			for (k = 0; k < size; k++) {
				tmp = omtrx[maxi][k];
				omtrx[maxi][k] = omtrx[j][k];
				omtrx[j][k] = tmp;
			}
			aw = -aw;
			wk[maxi] = wk[j];
		}
		index[j] = maxi;
		if (omtrx[j][j] == 0.0)
			omtrx[j][j] = eps;
		if (j != size - 1) {
			tmp = 1.0 / omtrx[j][j];
			for (i = j + 1; i < size; i++)
				omtrx[i][j] *= tmp;
		}
	}
	for (jx = 0; jx < size; jx++) {
		for (ix = 0; ix < size; ix++)
			wk[ix] = 0.0;
		wk[jx] = 1.0;
		l = -1;
		for (i = 0; i < size; i++) {
			idx = index[i];
			sum = wk[idx];
			wk[idx] = wk[i];
			if (l != -1) {
				for (j = l; j < i; j++)
					sum -= omtrx[i][j] * wk[j];
			} else if (sum != 0.0)
				l = i;
			wk[i] = sum;
		}
		for (i = size - 1; i >= 0; i--) {
			sum = wk[i];
			for (j = i + 1; j < size; j++)
				sum -= omtrx[i][j] * wk[j];
			wk[i] = sum / omtrx[i][i];
		}
		for (ix = 0; ix < size; ix++)
			imtrx[ix][jx] = wk[ix];
	}
	FREE(wk);
	FREE(index);
} /*_ luinverse */


/*
  symmetrize a matrix
*/
double sym_mat(double **mat, int m)
{
  double x,sum;
  int i,j;

  sum=0.0;
  for(i=0;i<m;i++)
    for(j=0;j<i;j++) {
      sum += fabs(mat[i][j]-mat[j][i]);
      x = mat[i][j]+mat[j][i];
      mat[i][j]=mat[j][i]=0.5*x;
    }
  return sum;
}

/* weighted least squares fitting

   X: m x n matrix of predictors
   Y: n vector of response
   W: n vector of weight

   min RSS(beta), where
   RSS = sum_{k=0,...,n} W[k]*(Y[k] - sum_{i=0,...,m} X[i][k]*beta[i])

   returns beta and rss
   *acmat may refer to the accuracy matrix

 */

double *lsfit(double **X, double *Y, double *W,
	      int m, int n,
	      double *beta, double *rss, double ***acmat)
{
  int i,j,k;
  double x,y;
  static double **covmat=NULL, **invmat=NULL, *xyvec=NULL;
  static int m0=0;

  /* memory allocation */
  if(m0!=m){
    covmat=renew_mat(covmat,m,m);
    invmat=renew_mat(invmat,m,m);
    xyvec=renew_vec(xyvec,m);
    m0=m;
  }
  if(!beta) beta=new_vec(m);
  if(acmat) *acmat=invmat; /* reference only */

  /* getting covariances */
  for(i=0;i<m;i++) {
    for(x=0.0,k=0;k<n;k++) x+=X[i][k]*Y[k]*W[k];
    xyvec[i]=x;
    for(j=0;j<=i;j++) {
      for(x=0.0,k=0;k<n;k++) x+=X[i][k]*X[j][k]*W[k];
      covmat[i][j]=covmat[j][i]=x;
    }
  }

  /* invmat = inverse matrix of covmat */
  luinverse(covmat,invmat,m);
  x=sym_mat(invmat,m);
  if(x>1e-5) warning("lsfit: covmat singularity %g",x);
  if(x>1e-3) warning("lsfit: COVARIANCE MATRIX IS SINGULAR");

  /* calculate the beta */
  for(i=0;i<m;i++) {
    for(x=0.0,j=0;j<m;j++) x+=invmat[i][j]*xyvec[j];
    beta[i]=x;
  }

  /* obtain the rss */
  for(x=0.0,k=0;k<n;k++) {
    for(y=0.0,i=0;i<m;i++) y+=X[i][k]*beta[i];
    x+=W[k]*(Y[k]-y)*(Y[k]-y);
  }
  *rss=x;
  
  return beta;
}

/**
   recipe
**/ 

static double maxarg1,maxarg2;
#define DMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
(maxarg1) : (maxarg2))

void dfsimple(double parm[], double diff[], double xh[],
	      int n, double(*func)(double []))
{
  int i,j;
  double x1,x2,y1,y2;
  static int n0=0; static double *ptmp=NULL;
  if(n>n0) { ptmp=renew_vec(ptmp,n); n0=n; }
  for(i=0;i<n;i++) {
    for(j=0;j<n;j++) ptmp[j]=parm[j];
    x1=ptmp[i]=parm[i]+xh[i];
    y1=(*func)(ptmp);
    x2=ptmp[i]=parm[i]-xh[i];
    y2=(*func)(ptmp);
    diff[i]=(y2-y1)/(x2-x1);
  }
}

#define CON 1.4
#define CON2 (CON*CON)
#define BIG 1.0e30
#define NTAB 10
#define SAFE 5.0
double dfridr(double (*func)(double), double x, double h, double *err)
{
  int i,j;
  double errt,fac,hh,ans;
  static double **a=0;
  
  if(h==0.0) error("h must be nonzero in dfridr");
  if(!a) a=new_mat(NTAB,NTAB);
  hh=h;
  a[0][0]=((*func)(x+hh)-(*func)(x-hh))/(2.0*hh);
  *err=BIG;
  for(i=1;i<NTAB;i++) {
    hh /= CON;
    a[0][i]=((*func)(x+hh)-(*func)(x-hh))/(2.0*hh);
    fac=CON2;
    for(j=1;j<=i;j++) {
      a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
      fac=CON2*fac;
      errt=DMAX(fabs(a[j][i]-a[j-1][i]),fabs(a[j][i]-a[j-1][i-1]));
      if(errt<=*err) {
	*err=errt; ans=a[j][i];
      }
    }
    if(fabs(a[i][i]-a[i-1][i-1]) >= SAFE*(*err)) break;
  }
  return ans;
}
#undef CON
#undef CON2
#undef BIG
#undef NTAB
#undef SAFE


static double (*dfmpridrfunc)(double []);
static int dfmpridrindex=0;
static int dfmpridrn0=0;
static double *dfmpridrpar=0;
static double dfmpridr1(double x)
{
  dfmpridrpar[dfmpridrindex]=x;
  return (*dfmpridrfunc)(dfmpridrpar);
}

void dfmpridr(double parm[], double diff[], double xh[],
	      int n, double(*func)(double []))
{
  int j; double err;
  if(n>dfmpridrn0) {
    dfmpridrpar=renew_vec(dfmpridrpar,n); dfmpridrn0=n;
  }
  dfmpridrfunc=func;
  for(dfmpridrindex=0;dfmpridrindex<n;dfmpridrindex++) {
    for(j=0;j<n;j++) dfmpridrpar[j]=parm[j];
    diff[dfmpridrindex]=
      dfridr(dfmpridr1,parm[dfmpridrindex],xh[dfmpridrindex],&err);
  }
}


#define SETPARM {for(k=0;k<n;k++) ptmp[k]=parm[k];}
void ddfsimple(double parm[], double **diff, double xh[],
	       int n, double(*func)(double []))
{
  int i,j,k;
  double y11,y12,y21,y22,y0;
  static int n0=0; static double *ptmp=NULL;
  if(n>n0) { ptmp=renew_vec(ptmp,n); n0=n; }
  SETPARM; y0=(*func)(ptmp);
  for(i=0;i<n;i++) {
    SETPARM; ptmp[i]+=xh[i]; y22=(*func)(ptmp);
    SETPARM; ptmp[i]-=xh[i]; y11=(*func)(ptmp);
    diff[i][i]=(y11+y22-2.0*y0)/(xh[i]*xh[i]);
    for(j=0;j<i;j++) {
      SETPARM; ptmp[i]+=xh[i]; ptmp[j]+=xh[j]; y22=(*func)(ptmp);
      SETPARM; ptmp[i]-=xh[i]; ptmp[j]-=xh[j]; y11=(*func)(ptmp);
      SETPARM; ptmp[i]+=xh[i]; ptmp[j]-=xh[j]; y21=(*func)(ptmp);
      SETPARM; ptmp[i]-=xh[i]; ptmp[j]+=xh[j]; y12=(*func)(ptmp);
      diff[i][j]=diff[j][i]=(y11+y22-y12-y21)/(xh[i]*xh[j]*4.0);
    }
  }
}
#undef SETPARM

void dfnmin(double p[], int n, double gtol, int itmax, int maxback,
	    int *iter, double *fret,
	    double ***hesinv,
	    double(*func)(double []), void (*dfunc)(double [], double []),
	    void (*ddfunc)(double [], double **))
{
  double *g, **A, **Ainv, *xi, *pnew;
  double sum,fp,fnew,lam,x;
  int loop,i,j,k;

  g=new_vec(n); xi=new_vec(n); pnew=new_vec(n);  
  A=new_mat(n,n); Ainv=new_mat(n,n);
  fp=(*func)(p); /* function */      
  for(loop=1;loop<=itmax;loop++) {
    (*dfunc)(p,g); /* derivative */
    (*ddfunc)(p,A); /* second derivative */
    luinverse(A,Ainv,n);
    x=sym_mat(Ainv,n);

    sum=0.0;
    for(i=0;i<n;i++) {
      x=0.0; for(j=0;j<n;j++) x+=Ainv[i][j]*g[j];
      sum+=g[i]*x; xi[i] = -x;
    }
    if(sum>=0.0) lam=1.0; else lam=-1.0;
    for(k=0;k<maxback;k++) {
      for(i=0;i<n;i++) pnew[i]=p[i]+lam*xi[i];
      fnew=(*func)(pnew); /* function */      
      mydprintf(3,"\n### dfnmin: lam=%g fnew=%g fp=%g",lam,fnew,fp);
      if(fnew < fp) break;
      lam *= 0.1;
    }
    if(k==maxback) break;
    fp=fnew;
    for(i=0;i<n;i++) p[i]=pnew[i];
    mydprintf(3,"\n### dfnmin: loop=%d sum=%g fp=%g",loop,sum,fp);
    if(sum>=0 && sum<gtol) break;
  }

  *fret=fp; /* function */
  *iter=loop;
  *hesinv=Ainv;

  free_mat(A); free_vec(g); free_vec(xi); free_vec(pnew);
}

