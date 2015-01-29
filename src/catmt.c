/*

  catmt.c : concanate mt files and selection by KH-test

  Time-stamp: <2001-06-23 11:41:24 shimo>

  shimo@ism.ac.jp 
  Hidetoshi Shimodaira

  typical usage:
  # foo1.mt foo2.mt -> foo.mt
  catmt foo1 foo2 foo
  # selecting by KH-test (pvalue >= 0.01)
  # generates goo.mt and goo.vt
  # goo.mt is the selected mt, goo.vt is the index file
  catmt --kht 0.01 foo goo

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "misc.h"
#include "rand.h"

static const char rcsid[] = "$Id: catmt.c,v 1.3 2001/08/10 06:04:55 shimo Exp $";

typedef struct {
  int m;
  int n;
  double **me; /* mm x nn matrix */
} mat;

mat *fcatread_mat(FILE *fp, mat *x)
{
  int i,j,m,n,n0;

  m=fread_i(fp); /* number of items (rows) */
  n=fread_i(fp); /* number of samples (columns) */
  
  if(x){
    if(m != x->m) error("matrix row size mismatch");
    n0=x->n; x->n+=n;
    for(i=0;i<m;i++) x->me[i]=renew_vec(x->me[i],x->n);
  } else {
    n0=0; x=NEW(mat); x->m=m; x->n=n;
    x->me=NEW_A(m,double*);
    for(i=0;i<m;i++) x->me[i]=new_vec(x->n);
  }
  for(i=0;i<m;i++)
    for(j=n0; j<x->n; j++)
      x->me[i][j]=fread_d(fp);

  return x;
}

void putdot() {putchar('.'); fflush(STDOUT);}
void byebye() {error("error in command line");}

double varadd=1.0; /* adding to wtmat */
#define VARADDWARN 100.0 /* warning if var < varadd*VARADDWARN */
#define HUGENUM 1.0e30;
double khtalpha=0.0; /* alpha of kh-test */
int sw_kht=0;


char *fname_vt = NULL;
char **fnamev_in = NULL;
int nfile;
char *fname_out = NULL;
char *fext_vt = ".vt";
char *fext_mt = ".mt";

mat *datmat=NULL;

int do_kht();

int main(int argc, char** argv)
{
  /* working variables */
  int i,ifile,n;
  FILE *fp;
  char *cbuf;

  printf("# %s",rcsid);

  fnamev_in=NEW_A(argc-1,char*);
  nfile=0;
  /* args */
  for(i=1;i<argc;i++) {
    if(argv[i][0] != '-') {
      fnamev_in[nfile]=argv[i];
      nfile++;
    } else if(streq(argv[i],"-d")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%d",&debugmode) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"--kht")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%lf",&khtalpha) != 1)
	byebye();
      sw_kht=1;
      i+=1;
    } else byebye();
  }

  if(nfile<2) {
    error("must specify input-file and output-file");
    byebye();
  } else {
    nfile--;
    fname_out=fnamev_in[nfile];
    fname_vt=rmvaxt(fname_out);
  }

  /* reading the mat */
  for(n=ifile=0;ifile<nfile;ifile++){
    fp=openfp(fnamev_in[ifile],fext_mt,"r",&cbuf);
    printf("\n# reading %s",cbuf);
    datmat=fcatread_mat(fp,datmat);
    printf(" %d",datmat->n-n); n=datmat->n;
    FREE(cbuf); fclose(fp);
  }
  printf("\n# M:%d N:%d",datmat->m,datmat->n);

  if(sw_kht) do_kht();

  /* writing the mat */
  fp=openfp(fname_out,fext_mt,"w",&cbuf);
  printf("\n# writing %s",cbuf);
  fwrite_mat(fp,datmat->me,datmat->m,datmat->n);
  FREE(cbuf); fclose(fp);

  printf("\n# exit normally\n");
  exit(0);
}

int do_kht()
{
  int i,j,i0,m,n,m0;
  double v,s,x,y,x0,mn;
  double *datvec,*khtpv,**X;
  mat *datnew; int *orderv;
  FILE *fp;
  char *cbuf;

  m=datmat->m; n=datmat->n; X=datmat->me;
  datvec=new_vec(m); khtpv=new_vec(m); orderv=new_ivec(m);
  x0=-HUGENUM;
  for(i=0;i<m;i++) {
    x=0.0; for(j=0;j<n;j++) x+=X[i][j];
    datvec[i]=x;
    if(x>x0) {x0=x; i0=i;}
  }

  printf("\n# %4s %6s %6s %5s %6s","item","obs","se","T","kh");
  for(i=0;i<m;i++) {
    if(i==i0) { 
      khtpv[i]=1.0; s=x=0.0;
    } else {
      mn=(datvec[i0]-datvec[i])/n;
      x=0.0; for(j=0;j<n;j++) {y=X[i0][j]-X[i][j]-mn; x+=y*y;}
      x=(double)n*x; 
      if(x<varadd*VARADDWARN)
	warning("small variance [%d,%d]: %g<%g*%g",i0+1,i+1,
		x,varadd,VARADDWARN);
      v=(x+varadd)/(n-1) ; /* variance */
      s=sqrt(v); /* std err */
      x=(datvec[i0]-datvec[i])/s; /* test statistic */
      khtpv[i]=1.0-poz(x);
    }
    printf("\n# %4d %6.2f %6.2f %5.2f %6.3f",
	   i+1,datvec[i0]-datvec[i],s,x,khtpv[i]);
  }

  /* screening */
  if(khtalpha>0.0) {
    datnew=NEW(mat); datnew->m=m; datnew->n=n;
    datnew->me=NEW_A(m,double*);
    m0=0;
    for(i=0;i<m;i++) if(khtpv[i]>=khtalpha) {
      datnew->me[m0]=datmat->me[i];
      orderv[m0]=i;
      m0++;
    } 
    datnew->m=m0;
    printf("\n# SCREENING BY KH > %6.3f\n# M:%d -> %d",
	   khtalpha,m,m0);
	   
    /* write order vector */
    fp=openfp(fname_vt,fext_vt,"w",&cbuf);
    printf("\n# writing %s",cbuf);
    fwrite_ivec(fp,orderv,m0);
    FREE(cbuf); fclose(fp);

    datmat=datnew;
  }

  return 0;
}
