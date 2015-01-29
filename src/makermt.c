/*

  makermt.c : make rmt-file by the RELL method

  Time-stamp: <2010-01-30 01:43:55 shimo>

  shimo@is.titech.ac.jp
  Hidetoshi Shimodaira

  typical usage:
  # sample.pa and foo.mt -> foo.rmt
  makermt -p sample foo
  # sample.pa and foo.mt -> foo.rmt and aho.vt 
  makermt -p sample foo aho
  # foo.svt -> aho1.rmt aho2.rmt (foo.svt = aho1 aho2)
  makermt -g foo

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rand.h"
#include "misc.h"
#include "freadmat.h"

static const char rcsid[] = "$Id: makermt.c,v 1.16 2010/01/29 16:45:36 shimo Exp $";


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
  double *wv,*rv,*xp;
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
      register int j; register double x;
      xp=datmat[k]; x=0.0;
      for(j=0;j<n;j++) x+=wv[j]*xp[j]; /* <- time consuming */
      repmat[k][i]=x/r; /* rescaling with 1/r */
    }
  }

  free_vec(wv); free_vec(rv);
  return repmat;
}

void putdot() {putchar('.'); fflush(STDOUT);}
void byebye() {error("error in command line");}

unsigned long seed=0;
char *fname_pa = NULL;
char *fname_mt = NULL;
char *fname_rmt = NULL;
char *fname_vt = NULL;
char *fname_svt = NULL;
char *fext_pa = ".pa";
char *fext_rmt = ".rmt";
char *fext_vt = ".vt";
char *fext_svt = ".svt";

double bbfact=1.0;

/* msboot parameter */
int kk,*bb,*bn;
double *rr,*rr1;
int kk0=10;
double rr0[]={0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4};
int bb0[]={10000,10000,10000,10000,10000,10000,10000,10000,10000,10000};
double rr00[]={1.0}; int bb00[]={10000}; int kk00=1; 
int sw_fastrep=0;
int sw_multi=0;

double **datmat=NULL;
int mm=0; int nn=0;
double **repmat=NULL;
double *datvec=NULL;


int genrmt(char *infile, char *outfile)
{
  int i,j;
  FILE *fp;
  double x,t0,t1;
  char *cbuf,*fext;

  /* open file */
  switch(seqmode) {
  case SEQ_MOLPHY: fext=fext_molphy; break;
  case SEQ_PAML: fext=fext_paml; break;
  case SEQ_PAUP: fext=fext_paup; break;
  case SEQ_PUZZLE: fext=fext_puzzle; break;
  case SEQ_PHYML: fext=fext_phyml; break;
  case SEQ_MT: 
  default: fext=fext_mt; break;
  }
  if(infile) {
    fp=openfp(infile,fext,"r",&cbuf);
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
  case SEQ_PHYML: 
    datmat = fread_mat_phyml(fp, &mm, &nn); break;
  case SEQ_MT: 
  default: 
    datmat = fread_mat(fp, &mm, &nn); break;  
  }
  if(infile) {fclose(fp);  FREE(cbuf);}
  printf("\n# M:%d N:%d",mm,nn);

  /* allocating buffers */
  datvec=new_vec(mm);
  bn=new_ivec(kk); rr1=new_vec(kk);

  /* calculate the log-likelihoods */
  for(i=0;i<mm;i++) {
    x=0; for(j=0;j<nn;j++) x+=datmat[i][j];
    datvec[i]=x;
  }
  
  /* calculate scales */
  for(i=0;i<kk;i++) {
    bn[i]=(int)(rr[i]*nn); /* sample size for bootstrap */
    rr1[i]=(double)bn[i]/nn; /* recalculate rr for integer adjustment */
  }

  /* open out file */
  if(outfile) {
    /* vt ascii write to file */
    fp=openfp(outfile,fext_vt,"w",&cbuf);
    printf("\n# writing %s",cbuf);
    fwrite_vec(fp,datvec,mm);
    fclose(fp); FREE(cbuf);
    /* rmt binary write to file */
    fp=openfp(outfile,fext_rmt,"wb",&cbuf);
    printf("\n# writing %s",cbuf);
    fwrite_bvec(fp,datvec,mm);
    fwrite_bvec(fp,rr1,kk);
    fwrite_bivec(fp,bb,kk);
    fwrite_bi(fp,kk);
  } else {
    /* rmt ascii write to stdout */
    printf("\n# writing to stdout");
    printf("\n# OBS:\n"); write_vec(datvec,mm);
    printf("\n# R:\n"); write_vec(rr1,kk);
    printf("\n# B:\n"); write_ivec(bb,kk);
    printf("\n# RMAT:\n");
    printf("%d\n",kk);
  }


  /* generating the replicates by resampling*/
  for(i=j=0;i<kk;i++) j+=bb[i];
  printf("\n# start generating total %d replicates for %d items",j,mm);
  fflush(STDOUT);
  t0=get_time();

  for(i=0;i<kk;i++) {
    repmat=new_lmat(mm,bb[i]);
    scaleboot(datmat,repmat,mm,nn,bn[i],bb[i]);
    if(outfile) {
      fwrite_bmat(fp,repmat,mm,bb[i]);
      putdot();
    } else {
      printf("\n## RMAT[%d]:\n",i); write_mat(repmat,mm,bb[i]);
    }
    free_lmat(repmat,mm);
  }

  t1=get_time();
  printf("\n# time elapsed for bootstrap t=%g sec",t1-t0);

  if(outfile) {
    fclose(fp); FREE(cbuf);
  }

  /* freeing buffers */
  free_vec(bn); free_vec(rr1); free_vec(datvec); free_mat(datmat);

  return 0;
}


int main(int argc, char** argv)
{
  /* working variables */
  int i,j,mf;
  FILE *fp;
  char *cbuf,**infiles,**outfiles;

  printf("# %s",rcsid);

  /* args */
  for(i=j=1;i<argc;i++) {
    if(argv[i][0] != '-') {
      switch(j) {
      case 1: fname_mt=argv[i]; break;
      case 2: fname_rmt=argv[i]; break;
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
    } else if(streq(argv[i],"-b")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%lf",&bbfact) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"-f")) {
      sw_fastrep=1;
    } else if(streq(argv[i],"-g")) {
      sw_multi=1;
    } else if(streq(argv[i],"--molphy")) {
      seqmode=SEQ_MOLPHY;
    } else if(streq(argv[i],"--paml")) {
      seqmode=SEQ_PAML;
    } else if(streq(argv[i],"--paup")) {
      seqmode=SEQ_PAUP;
    } else if(streq(argv[i],"--puzzle")) {
      seqmode=SEQ_PUZZLE;
    } else if(streq(argv[i],"--phyml")) {
      seqmode=SEQ_PHYML;
    } else if(streq(argv[i],"-d")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%d",&debugmode) != 1)
	byebye();
      i+=1;
    } else byebye();
  }

  /* random seed */
  smrand(seed);

  /* reading parameters */
  if(fname_pa!=NULL) {
    fp=openfp(fname_pa,fext_pa,"r",&cbuf);
    printf("\n# reading %s",cbuf);
    kk=0;
    rr=fread_vec(fp,&kk); 
    bb=fread_ivec(fp,&kk); 
    fclose(fp); FREE(cbuf);
  } else {
    if(sw_fastrep) {kk=kk00; rr=rr00; bb=bb00;}
    else {kk=kk0; rr=rr0; bb=bb0;}
  }

  for(i=0;i<kk;i++) bb[i] *= bbfact;
  printf("\n# seed:%lu (MT19937 generator)",seed);
  printf("\n# K:%d",kk);
  printf("\n# R:"); for(i=0;i<kk;i++) printf("%g ",rr[i]);
  printf("\n# B:"); for(i=0;i<kk;i++) printf("%d ",bb[i]);

  if(sw_multi) {
    if(fname_mt) fp=openfp(fname_mt,fext_svt,"r",&cbuf); else fp=STDIN;
    mf=0; infiles=fread_svec(fp,&mf);
    if(fname_mt) {fclose(fp); FREE(cbuf);}
    if(fname_rmt) {
      fp=openfp(fname_rmt,fext_svt,"r",&cbuf);
      outfiles=fread_svec(fp,&mf);
      fclose(fp); FREE(cbuf);
    } else {
      outfiles=NEW_A(mf,char*);
      for(i=0;i<mf;i++) outfiles[i]=rmvaxt(infiles[i]);
    }
    for(i=0;i<mf;i++) {
      printf("\n# %d/%d %s %s",i+1,mf,infiles[i],outfiles[i]);
      genrmt(infiles[i],outfiles[i]);
    }
  } else {
    if(fname_mt && !fname_rmt) fname_rmt=rmvaxt(fname_mt);
    genrmt(fname_mt,fname_rmt);
  }

  printf("\n# exit normally\n");
  exit(0);
}
