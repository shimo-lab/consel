/*

  makerep.c : make rep-file by the RELL method

  Time-stamp: <2011-01-25 16:56:43 shimo>

  shimo@ism.ac.jp 
  Hidetoshi Shimodaira

  typical usage:
  # sample.pa, hyp.ass and foo.mt -> foo.rep
  makerep -p sample -a hyp foo 

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rand.h"
#include "misc.h"
#include "freadmat.h"

static const char rcsid[] = "$Id: makerep.c,v 1.5 2011/05/12 07:20:53 shimo Exp $";

typedef struct {
  double *ve;
  int len;
} dvec;

void putdot() {putchar('.'); fflush(STDOUT);}
void byebye() {error("error in command line");}

unsigned long seed=0;
char *fname_pa = NULL; char *fext_pa = ".pa";
char *fname_mt = NULL;
char *fname_vt = NULL; char *fext_vt = ".vt";
char *fname_rep = NULL; char *fext_rep = ".rep";
char *fname_ass = NULL; char *fext_ass = ".ass";

/* msboot parameter */
int kk,*bb,*bn;
double *rr;
int kk0=10;
double rr0[]={0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4};
int bb0[]={10000,10000,10000,10000,10000,10000,10000,10000,10000,10000};

/* for associations */
int mm=0; /* no. of items */
int cm0=0; /* no. of associations (read) */
int **assvec0=NULL; /* association vectors (read) */
int *asslen0=NULL; /* lengths of the associations (read) */
int **cassvec0=NULL; /* complement of assvec0 */
int cm=0; /* truncate small items and use only top cm itmes */
int **assvec=NULL; /* sorted assvec */
int *asslen=NULL; /* sorted asslen */
int **cassvec=NULL; /* complement of assvec used in mctest */
int *orderv=NULL; /* cm-vector of the item id's */

/* for mt file */
int mm; int nn;
double **datmat=NULL; /* mm x nn matrix of data */
double *datvec=NULL;  /* mm-vector of data */

/* for rep file */
double *obsvec=NULL; /* cm-vector of observed statistics */

/* work */
double *buf1; /* mm-vector */
double *buf2; /* cm0-vector */

/* switches */
int sw_nosort=0; /* dont sort the items */
/* misc */
#define HUGENUM 1.0e30

/* routines */
double *calcmaxs(double *xx, int m);
double *calcmaxass(double *xx, double **wt, int m, 
		   int *ass, int *cass, int n, double *out);
double *calcassmins(double *xx, 
		    int **assvec, int *asslen, int cm, double *yy);
int compdvec(dvec *xp, dvec *yp);
int *getcass(int *ass, int n, int m, int *cass);
double *scaleboot1(double **datmat, double *repvec, int m,
		   int n, int bn, double *rw, double *rv);
void do_makerep(FILE *fp);


int main(int argc, char** argv)
{
  /* working variables */
  int i,j;
  double x,t0,t1;
  FILE *fp;
  char *cbuf,*fext;
  dvec **dvbuf;

  printf("# %s",rcsid);

  /* args */
  for(i=j=1;i<argc;i++) {
    if(argv[i][0] != '-') {
      switch(j) {
      case 1: fname_mt=argv[i]; break;
      case 2: fname_rep=argv[i]; break;
      default: byebye();
      }
      j++;
    } else if(streq(argv[i],"-s")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%lu",&seed) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"-p")) {
      if(i+1>=argc) byebye();
      fname_pa=argv[i+1];
      i+=1;
    } else if(streq(argv[i],"-a")) {
      if(i+1>=argc) byebye();
      fname_ass=argv[i+1];
      i+=1;
    } else if(streq(argv[i],"--molphy")) {
      seqmode=SEQ_MOLPHY;
    } else if(streq(argv[i],"--paml")) {
      seqmode=SEQ_PAML;
    } else if(streq(argv[i],"--paup")) {
      seqmode=SEQ_PAUP;
    } else if(streq(argv[i],"--puzzle")) {
      seqmode=SEQ_PUZZLE;
    } else byebye();
  }

  /* random seed */
  smrand(seed);

  /* reading parameters */
  if(fname_pa!=NULL) {
    fp=openfp(fname_pa,fext_pa,"r",&cbuf);
    printf("\n# reading %s",cbuf);
    kk=0;
    rr=fread_vec(fp,&kk); putdot();
    bb=fread_ivec(fp,&kk); putdot();
    fclose(fp); FREE(cbuf);
  } else {
    kk=kk0; rr=rr0; bb=bb0; 
  }
  printf("\n# seed:%lu",seed);
  printf("\n# K:%d",kk);
  printf("\n# R:"); for(i=0;i<kk;i++) printf("%g ",rr[i]);
  printf("\n# B:"); for(i=0;i<kk;i++) printf("%d ",bb[i]);

  /* open file */
  switch(seqmode) {
  case SEQ_MOLPHY: fext=fext_molphy; break;
  case SEQ_PAML: fext=fext_paml; break;
  case SEQ_PAUP: fext=fext_paup; break;
  case SEQ_PUZZLE: fext=fext_puzzle; break;
  case SEQ_MT: 
  default: fext=fext_mt; break;
  }
  if(fname_mt) {
    fp=openfp(fname_mt,fext,"r",&cbuf);
    printf("\n# reading %s",cbuf);
  } else {
    fp=STDIN;
    printf("\n# reading from stdin");
  }

  /* read file */
  mm=nn=0;
  switch(seqmode) {
  case SEQ_MOLPHY: 
    datmat = fread_mat_lls(fp, &mm, &nn); break;
  case SEQ_PAML: 
    datmat = fread_mat_lfh(fp, &mm, &nn); break;
  case SEQ_PAUP: 
    datmat = fread_mat_paup(fp, &mm, &nn); break;
  case SEQ_PUZZLE: 
    datmat = fread_mat_puzzle(fp, &mm, &nn); break;
  case SEQ_MT: 
  default: 
    datmat = fread_mat(fp, &mm, &nn); break;  
  }
  if(fname_mt) {fclose(fp);  FREE(cbuf);}
  printf("\n# M:%d N:%d",mm,nn);


  /* calculate bn */
  bn=new_ivec(kk);
  for(i=0;i<kk;i++) {
    bn[i]=(int)(rr[i]*nn); /* sample size for bootstrap */
    rr[i]=(double)bn[i]/nn; /* recalculate rr for integer adjustment */  
  }

  /* calculate the log-likelihoods */
  datvec=new_vec(mm);
  for(i=0;i<mm;i++) {
    x=0; for(j=0;j<nn;j++) x+=datmat[i][j];
    datvec[i]=x;
  }

  /* reading association vector */
  if(fname_ass){
    fp=openfp(fname_ass,fext_ass,"r",&cbuf);
    printf("\n# reading %s",cbuf);    
    cm0=fread_i(fp);
    assvec0=NEW_A(cm0,int*); asslen0=NEW_A(cm0,int);
    cassvec0=NEW_A(cm0,int*);
    for(i=0;i<cm0;i++) {
      asslen0[i]=0; assvec0[i]=fread_ivec(fp,&(asslen0[i]));
      for(j=0;j<asslen0[i];j++)
	if(assvec0[i][j]<0 || assvec0[i][j]>=mm)
	  error("ass out of range: ass[%d][%d]=%d",
		i+1,j+1,assvec0[i][j]+1);
    }
    fclose(fp); FREE(cbuf);
  } else { /* generate the identity association */
    printf("\n# generate the identity association");
    cm0=mm;
    assvec0=NEW_A(cm0,int*); asslen0=NEW_A(cm0,int);
    cassvec0=NEW_A(cm0,int*);
    for(i=0;i<cm0;i++) {
      asslen0[i]=1;
      assvec0[i]=NEW_A(1,int);
      assvec0[i][0]=i;
    }
  }
  for(i=0;i<cm0;i++) cassvec0[i]=getcass(assvec0[i],asslen0[i],mm,NULL);
  printf("\n# CM:%d",cm0);
  if(cm<=0 || cm>cm0) cm=cm0;
  if(cm<cm0) printf(" -> %d",cm);

  /* allocating the memory */
  buf1=new_vec(mm); buf2=new_vec(cm0); orderv=new_ivec(cm0);
  dvbuf=NEW_A(cm0,dvec*);
  obsvec=new_vec(cm);
  assvec=NEW_A(cm,int*); asslen=NEW_A(cm,int);
  cassvec=NEW_A(cm,int*);

  /* sorting the association by the test statistics */
  for(i=0;i<cm0;i++) {
    dvbuf[i]=NEW(dvec); dvbuf[i]->len=asslen0[i];
    dvbuf[i]->ve=new_vec(asslen0[i]);
    for(j=0;j<asslen0[i];j++) 
      calcmaxass(datvec,NULL,mm,assvec0[i],cassvec0[i],
		 asslen0[i],dvbuf[i]->ve);
    sort(dvbuf[i]->ve,NULL,dvbuf[i]->len);
  }
  if(sw_nosort) {
    for(i=0;i<cm0;i++) orderv[i]=i;
  } else {
    mypsort((void**)dvbuf,orderv,cm0,(int (*)(void *, void *))&compdvec);
  }
  mydprintf(1,"\n# observed statistics for the associations");
  mydprintf(1,"\n# %4s %4s ","rank","item");
  for(i=0;i<cm0;i++) {
    mydprintf(1,"\n# %4d %4d ",i+1,orderv[i]+1);
    for(j=0;j<dvbuf[i]->len;j++) mydprintf(1," %6.2f",dvbuf[i]->ve[j]);
  }

  /* copying and truncating the sorted associations */
  for(i=0;i<cm;i++) {
    obsvec[i]=dvbuf[i]->ve[0]; /* take the min */
    asslen[i]=asslen0[orderv[i]];
    assvec[i]=NEW_A(asslen[i],int);
    for(j=0;j<asslen[i];j++) assvec[i][j]=assvec0[orderv[i]][j];
    cassvec[i]=NEW_A(mm-asslen[i],int);
    for(j=0;j<mm-asslen[i];j++) cassvec[i][j]=cassvec0[orderv[i]][j];
  }

  /* freeing unnecessary associations */
  for(i=0;i<cm0;i++) {
    free_vec(dvbuf[i]->ve); FREE(dvbuf[i]);
    free_ivec(assvec0[i]); free_ivec(cassvec0[i]);
  }
  FREE(dvbuf); FREE(assvec0); FREE(cassvec0); free_ivec(asslen0); 

  /* make rep and write rep */
  if(fname_mt && !fname_rep) fname_rep=rmvaxt(fname_mt);  
  if(fname_rep) { /* binary write to file */
    fp=openfp(fname_rep,fext_rep,"wb",&cbuf);
    printf("\n# writing %s",cbuf);
    fwrite_bivec(fp,orderv,cm); fwrite_bvec(fp,obsvec,cm);
    fwrite_bvec(fp,rr,kk); fwrite_bivec(fp,bb,kk);
    fwrite_bi(fp,kk);
  } else { /* ascii write to stdout */
    fp=NULL;
    printf("\n# writing to stdout\n");
    write_ivec(orderv,cm); write_vec(obsvec,cm);
    write_vec(rr,kk); write_ivec(bb,kk);
    write_i(kk);
  }
  fflush(STDOUT);
  t0=get_time();
  do_makerep(fp);
  t1=get_time();
  fclose(fp); FREE(cbuf);
  printf("\n# time elapsed for bootstrap t=%g sec",t1-t0);

  printf("\n# exit normally\n");
  exit(0);
}

void do_makerep(FILE *fp)
{
  double *repvec; /* m-vector of a replicate of dat */
  double *statvec; /* cm-vector of a replicate of statistic */
  double **statmat; /* cm x bb[k] matrix of replicates */
  double *rw; /* n-vector */
  double *rv; /* bn-vector */
  int i,j,k;


  repvec=new_vec(mm); statvec=new_vec(cm); rw=new_vec(nn); 
  statmat=NULL; rv=NULL;

  for(k=0;k<kk;k++) {
    rv=renew_vec(rv,bn[k]);
    statmat=renew_mat(statmat,cm,bb[k]);
    for(i=0;i<bb[k];i++) {
      scaleboot1(datmat,repvec,mm,nn,bn[k],rw,rv);
      calcmaxs(repvec,mm);
      calcassmins(repvec,assvec,asslen,cm,statvec);
      for(j=0;j<cm;j++) statmat[j][i]=statvec[j];
    }
    if(fp) {fwrite_bmat(fp,statmat,cm,bb[k]); putdot();}
    else write_mat(statmat,cm,bb[k]);
  }

  free_vec(rv); free_mat(statmat);
  free_vec(repvec); free_vec(statvec); free_vec(rw);
}


/*
  Scaleboot1 simply generates a replicate of sample size bn. The
  original data have n samples of m items. The replicate is
  calculated as the sum of the bn vectors of size n.
*/
double *scaleboot1(double **datmat, /* m x n data matrix */
		   double *repvec, /* m-vector output */
		   int m,  /* number of items */
		   int n,  /* original sample size */
		   int bn,  /* replicate sample size */
		   double *wv, /* n-vector */
		   double *rv /* bn-vector */
		   ) {
  int j,k;
  double *xp,r;

  if(repvec==NULL) repvec=new_vec(m);
  r=(double)bn/(double)n;

  /* first, get wv=weight vector */
  mrandlist(rv,bn);
  for(j=0;j<n;j++) wv[j]=0.0;
  for(j=0;j<bn;j++) {
    k = (int) (rv[j]*n);
    if(k<0) {  /* for safety */
      k=0; warning("negative? in rand");
    } else if(k>=n) {
      k=n-1; warning("unity? in rand");
    }
    wv[k] += 1.0;
  }

  /* then, summing up datmat with the weights */
  for(k=0;k<m;k++) {
    register int j; register double x;
    xp=datmat[k]; x=0.0;
    for(j=0;j<n;j++) x+=wv[j]*xp[j]; /* <- time consuming */
    repvec[k]=x/r; /* rescaling with 1/r */
  }

  return repvec;
}


int compdvec(dvec *xp, dvec *yp)
{
  int i,len,out;
  
  if(xp->len < yp->len) { out=-1; len=xp->len;}
  else if(xp->len > yp->len) { out=1; len=yp->len; }
  else { out=0; len=yp->len; }
  for(i=0;i<len;i++) {
    if(xp->ve[i] < yp->ve[i]) {out=-1; break;}
    else if(xp->ve[i] > yp->ve[i]) {out=1; break;}
  }
  return out;
}


/*
  acc = size n
  cacc = size >= m-n
  acc U cacc = 0,...,m-1, disjoint.
 */
int *getcass(int *ass, int n, int m, int *cass)
{
  int i,j;
  static int m0=0;
  static int *buf=NULL;

  if(m>m0) { buf=renew_ivec(buf,m); m0=m; }

  for(i=0;i<m;i++) buf[i]=0;
  for(i=0;i<n;i++) buf[ass[i]]++;
  for(i=j=0;i<m;i++) if(buf[i]==0) j++;
  if(!cass) cass=new_ivec(j); /* j >= m-n */
  if(j != (m-n)) error("\n# duplicate associations");
  for(i=j=0;i<m;i++) if(buf[i]==0) cass[j++]=i;

  return cass;
}


/*
  out[i] = max_{j ne i} xx[j] - xx[i], then xx <- out.
  compute it in two pass for insitu
 */
double *calcmaxs(double *xx, int m)
{
  int i;
  double x;
  static int m0=0;
  static double *yy=NULL;
  if(m>m0) { yy=renew_vec(yy,m); m0=m; }

  x=-HUGENUM;
  for(i=0;i<m;i++) {
    yy[i]=x;
    if(xx[i]>x) x=xx[i];
  }
  x=-HUGENUM;
  for(i=m-1;i>=0;i--) {
    if(x>yy[i]) yy[i]=x;
    if(xx[i]>x) x=xx[i];
  }
  /* now yy[i] = max(xx[0],x[1],...,x[i-1],x[i+1],...,x[m-1]) */
  for(i=0;i<m;i++) xx[i]=yy[i]-xx[i];

  return xx;
}

/*
  xx = m-vector of input
  ass = n-vector of association
  cass = (m-n)-vector of complementary of ass
  out = n-vector of output
  wt = m x m matrix for weight

  out[i] = max_{j in iass} wt[i,j](xx[j] - xx[i]), for i in ass
 */
double *calcmaxass(double *xx, double **wt, int m, 
		   int *ass, int *cass, int n, double *out)
{
  int i,j,ia,ja;
  double x,x0,y,*wp;

  if(wt) {
    for(i=0;i<n;i++) {
      ia=ass[i]; x0=xx[ia]; wp=wt[ia]; x=-HUGENUM;
      for(j=0;j<m-n;j++) {
	ja=cass[j]; y=wp[ja]*(xx[ja]-x0);
	if(y>x) x=y;
      }
      out[i]=x;
    } 
  } else {
    x=-HUGENUM;
    for(j=0;j<m-n;j++) {
      y=xx[cass[j]]; if(y>x) x=y;
    }
    for(i=0;i<n;i++) out[i]=x-xx[ass[i]];
  }
  return out;
}

/*
  yy[i] = min_{j in assvec[i]} xx[j]
  for i=1,...,cm, where
  assvec[i] is of size asslen[i].
 */
double *calcassmins(double *xx, 
		    int **assvec, int *asslen, int cm,
		    double *yy)
{
  int i,j; double x,y;
  for(i=0;i<cm;i++) {
    y=HUGENUM;
    for(j=0;j<asslen[i];j++) {
      x=xx[assvec[i][j]]; if(x<y) y=x;
    }
    yy[i]=y;
  }

  return yy;
}

