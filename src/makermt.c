/*

  makermt.c : make rmt-file by the RELL method

  Time-stamp: <2001-04-19 09:22:33 shimo>

  shimo@ism.ac.jp 
  Hidetoshi Shimodaira

  typical usage:
  # sample.pa and foo.mt -> foo.rmt
  makermt -p sample foo
  # sample.pa and foo.mt -> foo.rmt and aho.vt 
  makermt -p sample foo aho

*/

#include <stdio.h>
#include <math.h>
#include "rand.h"
#include "misc.h"

static const char rcsid[] = "$Id: makermt.c,v 1.3 2001/04/16 07:01:51 shimo Exp shimo $";


/*
  Scaleboot simply generates the bb replicates of sample size bn. The
  original data have n samples of m items. Each replicate is
  calculated as the sum of the bn vectors of size n.
*/
double **scaleboot(double **datmat, /* m x n data matrix */
		   double **repmat, /* m x bb output */
		   int m,  /* number of items */
		   int n,  /* original sample size */
		   int bn, /* replicate sample size */
		   int bb  /* number of bootstrap replicates */
		   ) {
  int i,j,k;
  double x,*wv,*rv;
  double r;

  wv=new_vec(n); rv=new_vec(bn);
  if(repmat==NULL) repmat=new_mat(m,bb);
  r=(double)bn/(double)n;

  for(i=0;i<bb;i++) {
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
      x=0.0;
      for(j=0;j<n;j++) x+=wv[j]*datmat[k][j];
      repmat[k][i]=x/r;
    }
  }

  free_vec(wv); free_vec(rv);
  return repmat;
}



void putdot() {putchar('.'); fflush(STDOUT);}
void byebye() {error("error in command line");}

int seed=123;
char *fname_pa = NULL;
char *fname_mt = NULL;
char *fname_rmt = NULL;
char *fname_vt = NULL;
char *fext_pa = ".pa";
char *fext_mt = ".mt";
char *fext_rmt = ".rmt";
char *fext_vt = ".vt";

/* msboot parameter */
int kk,*bb,*bn;
double *rr;
int kk0=10;
double rr0[]={0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4};
int bb0[]={1000,1000,1000,1000,1000,1000,1000,1000,1000,1000};

double **datmat=NULL;
int mm; int nn;
double ***repmats=NULL;
double *datvec=NULL;

int main(int argc, char** argv)
{
  /* working variables */
  int i,j;
  double x,t0,t1;
  FILE *fp;

  printf("# %s",rcsid);

  /* args */
  for(i=j=1;i<argc;i++) {
    if(argv[i][0] != '-') {
      switch(j) {
      case 1: fname_mt=fname_rmt=argv[i]; break;
      case 2: fname_rmt=argv[i]; break;
      default: byebye();
      }
      j++;
    } else if(streq(argv[i],"-s")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%d",&seed) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"-p")) {
      if(i+1>=argc) byebye();
      fname_pa=argv[i+1];
      i+=1;
    } else if(streq(argv[i],"-l")) {
      if(i+1>=argc) byebye();
      fname_vt=argv[i+1];
      i+=1;
    } else byebye();
  }

  /* random seed */
  smrand(seed);

  /* reading parameters */
  if(fname_pa!=NULL) {
    fname_pa = mstrcat(fname_pa,fext_pa);
    if((fp=fopen(fname_pa,"r"))==NULL) error("cant open %s",fname_pa);
    printf("\n# reading %s.",fname_pa);
    kk=0;
    rr=fread_vec(fp,&kk); putdot();
    bb=fread_ivec(fp,&kk); putdot();
    fclose(fp);
  } else {
    kk=kk0; rr=rr0; bb=bb0; 
  }
  printf("\n# seed:%d",seed);
  printf("\n# K:%d",kk);
  printf("\n# R:"); for(i=0;i<kk;i++) printf("%g ",rr[i]);
  printf("\n# B:"); for(i=0;i<kk;i++) printf("%d ",bb[i]);

  /* reading datmat */
  mm=nn=0;
  if(fname_mt!=NULL) {
    fname_mt = mstrcat(fname_mt,fext_mt);
    if((fp=fopen(fname_mt,"r"))==NULL) error("cant open %s",fname_mt);
    printf("\n# reading %s.",fname_mt);
    datmat=fread_mat(fp,&mm,&nn); putdot();
    fclose(fp);
  } else {
    printf("\n# reading from stdin.");
    datmat=fread_mat(STDIN,&mm,&nn); putdot();
  }
  printf("\n# M:%d N:%d",mm,nn);

  /* allocating buffers */
  datvec=new_vec(mm);
  repmats=(double ***)MALLOC((sizeof(double **))*kk);
  for(i=0;i<kk;i++) repmats[i]=new_mat(mm,bb[i]);
  bn=new_ivec(kk);

  /* calculate the log-likelihoods */
  for(i=0;i<mm;i++) {
    x=0; for(j=0;j<nn;j++) x+=datmat[i][j];
    datvec[i]=x;
  }

  /* generating the replicates by resampling*/
  for(i=j=0;i<kk;i++) j+=bb[i];
  printf("\n# start generating total %d replicates for %d items",j,mm);
  fflush(STDOUT);
  t0=get_time();
  for(i=0;i<kk;i++) {
    bn[i]=(int)(rr[i]*nn); /* sample size for bootstrap */
    rr[i]=(double)bn[i]/nn; /* recalculate rr for integer adjustment */
    scaleboot(datmat,repmats[i],mm,nn,bn[i],bb[i]);
    putdot();
  }
  t1=get_time();
  printf("\n# time elapsed for bootstrap t=%g sec",t1-t0);

  /* saving the results */
  if(fname_vt!=NULL) {
    fname_vt = mstrcat(fname_vt,fext_vt);
    if((fp=fopen(fname_vt,"w"))==NULL) error("cant open %s",fname_vt);
    printf("\n# writing %s.",fname_vt);
    fwrite_vec(fp,datvec,mm); putdot();
    fclose(fp);
  }
  if(fname_rmt!=NULL) {
    fname_rmt = mstrcat(fname_rmt,fext_rmt);
    if((fp=fopen(fname_rmt,"w"))==NULL) error("cant open %s",fname_rmt);
    printf("\n# writing %s.",fname_rmt);
    fwrite_bvec(fp,datvec,mm);
    fwrite_bvec(fp,rr,kk);
    fwrite_bivec(fp,bb,kk);
    fwrite_bi(fp,kk);
    for(i=0;i<kk;i++) {
      fwrite_bmat(fp,repmats[i],mm,bb[i]);
      putdot();
    }
    fclose(fp);
  } else {
    printf("\n# writing to stdout.");
    printf("\n# L:\n"); write_vec(datvec,mm);
    printf("\n# R:\n"); write_vec(rr,kk);
    printf("\n# B:\n"); write_ivec(bb,kk);
    printf("\n# RMAT:\n");
    printf("%d\n",kk);
    for(i=0;i<kk;i++) {
      printf("\n## RMAT[%d]:\n",i); write_mat(repmats[i],mm,bb[i]);
    }
  }

  printf("\n# exit normally\n");
  exit(0);
}
