/*

  consel.c : assessing the confidence in selection
             using the multi-scale bootstrap

  Time-stamp: <2001-05-05 17:52:27 shimo>

  shimo@ism.ac.jp 
  Hidetoshi Shimodaira

  typical usage:
  # foo.rmt -> foo.pv, foo.ci
  consel foo
  # foo.rep -> foo.pv, foo.ci
  consel -R foo
  # foo.cnt -> foo.pv
  consel -C foo
  # foo.rmt tree.ass -> aho.pv aho.ci aho.rep aho.cnt
  consel foo aho -a tree -r -c
  # foo.rep -> aho.pv aho.ci
  consel -R foo aho
  #
  # consel [switches] infile [outfile]
  #  reads "infile.{rmt,rep,cnt}" to output "outfile.{pv,ci,rep,cnt}"
  #  (1) outfile = infile when omitted
  #  (2) infile.rmt is read if neither -R nor -C is specified
  #  (3) infile.rep is read if -R is specified
  #  (4) infile.cnt is read if -C is specified
  #
  # switches
  # -R, -C: specifies input mode
  # -a NAME: reads "NAME.ass" for association
  # -r: produces outfile.rep
  # -c: produces outfile.cnt
  # -d VAL: debug output of level VAL
  #
*/

static const char rcsid[] = "$Id$";

#include <stdio.h>
#include <math.h>
#include "rand.h"
#include "misc.h"

typedef struct {
  double *ve;
  int len;
} dvec;

/* count */
double cntdist(double *vec, int bb, double t, int modesw);

/* mc-test routines */
double *calcmaxs(double *xx, int m); 
void repminmaxs(double **repmat, double **statmat, 
		int m, int bb, int cm, 
		int **assvec, int *asslen, double *buf1, double *buf2);
int compdvec(dvec *xp, dvec *yp);
int *getcass(int *ass, int n, int m, int *cass);
double pvaldist(double *vec, int bb, double t);
double calcmcpval(double *datvec, double **repmat, int mm, int nb,
		  int *ass, int *cass, int alen, double **wt);
double calckhpval(double *datvec, double **repmat, int mm, int nb,
		  int *ass, int *cass, int alen, double **wt);
double *getseval(double *pval, int m, int nb, double *seval);

/* multiscale bootstrap routines */
int rcalpval(double *cnts, double *rr, int *bb, int kk,
	     double *pv, double *se, double *rss, double *df, 
	     double **betap, double ***vmatp);
int vcalpval(double **statps, double *rr, int *bb, int kk,
	     double threshold, double *thp,
	     double *pvp, double *sep, double *rssp, double *dfp, 
	     double **betap, double ***vmatp);
int invtpval(double **statps, double *rr, int *bb, int kk,
	     double alpha,
	     double *cival, double *seval, double *eival);

/* for main */
void putdot() {putchar('.'); fflush(STDOUT);}
void byebye() {error("error in command line");}


/* file name and their extensions */
char *fname_rmt = NULL; char *fext_rmt = ".rmt";
char *fname_cnt = NULL; char *fext_cnt = ".cnt";
char *fname_rep = NULL; char *fext_rep = ".rep";
char *fname_ass = NULL; char *fext_ass = ".ass";
char *fname_pv = NULL; char *fext_pv = ".pv";
char *fname_ci = NULL; char *fext_ci = ".ci";
char *fname_vt = NULL; char *fext_vt = ".vt";

/* for the scales */
int kk=0; /* no. of scales */
int *bb=NULL; /* kk-vector of no.'s of replicates */
double *rr=NULL; /* kk-vector of relative sample sizes */

/* for associations */
int mm=0; /* no. of items */
int cm0=0; /* no. of associations (read) */
int **assvec0=NULL; /* association vectors (read) */
int *asslen0=NULL; /* lengths of the associations (read) */
int cm=0; /* truncate small items and use only top cm itmes */
int **assvec=NULL; /* sorted assvec */
int *asslen=NULL; /* sorted asslen */
int **cassvec=NULL; /* complement of assvec used in mctest */
int *orderv=NULL; /* cm-vector of the item id's */

/* for rmt file */
double ***repmats=NULL; /* (kk,mm,bb[i])-array of replicates */
double *datvec=NULL; /* mm-vector of data */

/* for rep file */
double ***statmats=NULL; /* (kk,cm,bb[i])-array of test statistics */
double *obsvec=NULL; /* cm-vector of observed statistics */

/* for cnt file */
double **cntmat=NULL; /* (kk,cm)-matrix of counts */

/* misc */
double threshold=0.0; /* the region is defined as x <= thereshold */
double varadd=1.0; /* adding to wtmat */
#define HUGENUM 1.0e30;

/* for mc-tests: p-values and se (cm-vector) */
double *npvec, *nsvec; /* naive p-value */
double *mcpvec, *mcsvec; /* p-value of nonweighted mc-test */
double *mspvec, *mssvec; /* p-value of weighted mc-test */
double *khpvec, *khsvec; /* kh-test p-value */

/* for msboot: p-values and se (cm-vector) */
double *pvvec, *sevec ; /* multi-scale pvalue */
double *pv0vec, *se0vec; /* pvalue for the derived naive method */

/* for msboot pv */
double **betamat=NULL; /* (2,cm)-matrix of signed distance, curvature */
double *rssvec=NULL; /* cm-vec rss */
double *dfvec=NULL; /* cm-vec df */
double *pfvec=NULL; /* cm-vec pvalue of the diagnostic */
double *thvec=NULL; /* cm-vec actual values of thesholds */
double rrmin, rrmax;
double pfmin; int dfmin;

/* for msboot ci */
int nalpha=0, nalpha0=5; /* number of alpha's */
double *alphavec=NULL; /* nalpha-vector of alpha's */
double alphavec0[]={0.05,0.10,0.5,0.90,0.95}; /* default value */
double **cimat=NULL; /* (cm,nalpha)-mat of confidence limits */
double **csmat=NULL; /* se of cimat */
double **eimat=NULL; /* error index of cimat */
double cieps;

/* work */
double *buf1; /* mm-vector */
double *buf2; /* cm0-vector */

/* switches */
int sw_pvaldist;
int sw_outrep=0; /* output rep file */
int sw_outcnt=0; /* output cnt file */
int sw_inrep=0; /* input rep file */
int sw_incnt=0; /* input cnt file */
int sw_smoothcnt=0; /* smooth cnt file */
int sw_domc=1; /* dont skip mc-tests */
int sw_doau=1; /* dont skip au-tests */

/* main routines */
int do_rmtmode();
int do_repmode();
int do_cntmode();
/* sub routines */
int do_mctest();
int do_bootrep();
int do_bootcnt();
/* output routines */
int write_pv(int sw_mc, int sw_au);
int write_ci();
int write_rep();
int write_cnt();

int main(int argc, char** argv)
{
  /* working variables */
  int i,j;
  char *fname_in=NULL, *fname_out=NULL;

  printf("# %s",rcsid);

  /* args */
  for(i=j=1;i<argc;i++) {
    if(argv[i][0] != '-') {
      switch(j) {
      case 1: fname_in=fname_out=argv[i]; break;
      case 2: fname_out=argv[i]; break;
      default: byebye();
      }
      j++;
    } else if(streq(argv[i],"-R")) {
      sw_inrep=1;
    } else if(streq(argv[i],"-C")) {
      sw_incnt=1;
    } else if(streq(argv[i],"-r")) {
      sw_outrep=1;
    } else if(streq(argv[i],"-c")) {
      sw_outcnt=1;
    } else if(streq(argv[i],"-a")) {
      if(i+1>=argc) byebye();
      fname_ass=argv[i+1];
      i+=1;
    } else if(streq(argv[i],"-v")) {
      if(i+1>=argc) byebye();
      fname_vt=argv[i+1];
      i+=1;
    } else if(streq(argv[i],"-d")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%d",&debugmode) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"-t")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%d",&cm) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"--skip_mc")) {
      sw_domc=0;
    } else if(streq(argv[i],"--skip_au")) {
      sw_doau=0;
    } else if(streq(argv[i],"--smooth_cnt")) {
      sw_smoothcnt=1;
    } else if(streq(argv[i],"--quick_mctest")) {
      sw_pvaldist=0; /* avoid sorting the replicates */
    } else if(streq(argv[i],"--threshold")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%lf",&threshold) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"--cieps")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%lf",&cieps) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"--rmin")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%lf",&rrmin) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"--rmax")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%lf",&rrmax) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"--pmin")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%lf",&pfmin) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"--dmin")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%d",&dfmin) != 1)
	byebye();
      i+=1;
    } else byebye();
  }

  if(sw_incnt+sw_inrep>1) error("only one of -C and -R can be specified");
  if(sw_incnt+sw_inrep==0) { /* rmt mode */
    fname_rmt=fname_in;
    fname_pv=fname_ci=fname_cnt=fname_rep=fname_out;
    do_rmtmode();
  } else if(sw_inrep) { /* rep mode */
    fname_rep=fname_in;
    fname_pv=fname_ci=fname_cnt=fname_out;
    do_repmode();
  } else if(sw_incnt) { /* cnt mode */
    fname_cnt=fname_in;
    fname_pv=fname_out;
    do_cntmode();
  }

  printf("\n# exit normally\n");
  exit(0);
}

/*
  (1) input rmt file
  (2) read association file and sort the items
  (3) call mctest and bootrep
*/
int do_rmtmode()
{
  int i,j;
  char *cbuf;
  FILE *fp;
  dvec **dvbuf;

  /* reading rmt */
  mm=kk=0;
  if(fname_rmt){ /* binary read from file */
    fp=openfp(fname_rmt,fext_rmt,"r",&cbuf);
    printf("\n# reading %s",cbuf);
    datvec=fread_bvec(fp,&mm); 
    rr=fread_bvec(fp,&kk); 
    bb=fread_bivec(fp,&kk); 
    i=fread_bi(fp); if(i != kk) error("wrong size in rmt");
    repmats=NEW_A(kk,double**);
    for(i=0;i<kk;i++) {
      repmats[i]=fread_bmat(fp,&mm,&(bb[i])); putdot();
    }
    fclose(fp); FREE(cbuf);
  } else { /* ascii read from stdin */
    printf("\n# reading from stdin");
    datvec=read_vec(&mm);
    rr=read_vec(&kk);
    bb=read_ivec(&kk);
    i=read_i(); if(i != kk) error("wrong size in rmt");
    repmats=NEW_A(kk,double**);
    for(i=0;i<kk;i++) {
      repmats[i]=read_mat(&mm,&(bb[i])); putdot();
    }
  }
  printf("\n# K:%d",kk);
  printf("\n# R:"); for(i=0;i<kk;i++) printf("%g ",rr[i]);
  printf("\n# B:"); for(i=0;i<kk;i++) printf("%d ",bb[i]);
  printf("\n# M:%d",mm);

  /* reading association vector */
  if(fname_ass){
    fp=openfp(fname_ass,fext_ass,"r",&cbuf);
    printf("\n# reading %s",cbuf);    
    cm0=fread_i(fp);
    assvec0=NEW_A(cm0,int*); asslen0=NEW_A(cm0,int);
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
    for(i=0;i<cm0;i++) {
      asslen0[i]=1;
      assvec0[i]=NEW_A(1,int);
      assvec0[i][0]=i;
    }
  }
  printf("\n# CM:%d",cm0);
  if(cm<=0 || cm>cm0) cm=cm0;
  if(cm<cm0) printf(" -> %d",cm);

  /* allocating the memory */
  buf1=new_vec(mm); buf2=new_vec(cm0); orderv=new_ivec(cm0);
  dvbuf=NEW_A(cm0,dvec*);
  statmats=NEW_A(kk,double**);
  for(i=0;i<kk;i++) statmats[i]=new_mat(cm,bb[i]);
  cntmat=new_mat(cm,kk);
  obsvec=new_vec(cm);
  assvec=NEW_A(cm,int*); asslen=NEW_A(cm,int);
  cassvec=NEW_A(cm,int*);

  /* sorting the association by the test statistics */
  for(i=0;i<mm;i++) buf1[i]=datvec[i];
  calcmaxs(buf1,mm);
  for(i=0;i<cm0;i++) {
    dvbuf[i]=NEW(dvec); dvbuf[i]->len=asslen0[i];
    dvbuf[i]->ve=new_vec(asslen0[i]);
    for(j=0;j<asslen0[i];j++) dvbuf[i]->ve[j]=buf1[assvec0[i][j]];
    sort(dvbuf[i]->ve,NULL,dvbuf[i]->len);
  }
  psort((void**)dvbuf,orderv,cm0,(int (*)(void *, void *))&compdvec);
  dprintf(1,"\n# observed statistics for the associations");
  dprintf(1,"\n# %4s %4s ","rank","item");
  for(i=0;i<cm0;i++) {
    dprintf(1,"\n# %4d %4d ",i+1,orderv[i]+1);
    for(j=0;j<dvbuf[i]->len;j++) dprintf(1," %6.2f",dvbuf[i]->ve[j]);
  }

  /* copying and truncating the sorted associations */
  for(i=0;i<cm;i++) {
    obsvec[i]=dvbuf[i]->ve[0]; /* take the min */
    asslen[i]=asslen0[orderv[i]];
    assvec[i]=NEW_A(asslen[i],int);
    for(j=0;j<asslen[i];j++) assvec[i][j]=assvec0[orderv[i]][j];
    cassvec[i]=getcass(assvec[i],asslen[i],mm,NULL); /* for mctest */
  }

  /* freeing unnecessary associations */
  for(i=0;i<cm0;i++) {
    free_vec(dvbuf[i]->ve); FREE(dvbuf[i]);
    free_ivec(assvec0[i]);
  }
  FREE(dvbuf); FREE(assvec0); free_ivec(asslen0);

  /* reading alphas (for bootrep) */
  if(fname_vt) {
      fp=openfp(fname_vt,fext_vt,"r",&cbuf);
      printf("\n# reading %s",cbuf);
      nalpha=0; alphavec=fread_vec(fp,&nalpha);
      fclose(fp); FREE(cbuf);
  } else {
    nalpha=nalpha0; alphavec=alphavec0;
  }

  /* conventional methods */
  if(sw_domc) do_mctest();

  if(sw_doau) {
    /* calculate the statistics for the replicates */
    printf("\n# calculate replicates of the statistics");
    for(i=0;i<kk;i++) {
      repminmaxs(repmats[i],statmats[i],mm,bb[i],cm,
		 assvec,asslen,buf1,buf2);
      putdot();
    }

    /* multie-scale bootstrap */
    do_bootrep();
  }

  /* output results */
  write_pv(sw_domc,sw_doau);
  if(sw_doau) {
    write_ci();
    if(sw_outrep) write_rep();
    if(sw_outcnt) write_cnt();
  }

  return 0;
}

int do_repmode()
{
  int i;
  char *cbuf;
  FILE *fp;

  /* reading rep */
  cm=kk=0;
  if(fname_rep){ /* binary read from file */
    fp=openfp(fname_rep,fext_rep,"r",&cbuf);
    printf("\n# reading %s",cbuf);
    orderv=fread_bivec(fp,&cm); obsvec=fread_bvec(fp,&cm);
    rr=fread_bvec(fp,&kk); bb=fread_bivec(fp,&kk);
    i=fread_bi(fp);
    if(i!=kk) error("wrong number of matrices");
    statmats=NEW_A(kk,double**);
    for(i=0;i<kk;i++) {
      statmats[i]=fread_bmat(fp,&cm,bb+i); putdot();
    }
    fclose(fp); FREE(cbuf);
  } else { /* ascii read from stdin */
    printf("\n# reading from stdin");
    orderv=read_ivec(&cm); obsvec=read_vec(&cm);
    rr=read_vec(&kk); bb=read_ivec(&kk);
    i=read_i();
    if(i!=kk) error("wrong number of matrices");
    statmats=NEW_A(kk,double**);
    for(i=0;i<kk;i++) {
      statmats[i]=read_mat(&cm,bb+i); putdot();
    }
  }
  printf("\n# K:%d",kk);
  printf("\n# R:"); for(i=0;i<kk;i++) printf("%g ",rr[i]);
  printf("\n# B:"); for(i=0;i<kk;i++) printf("%d ",bb[i]);
  printf("\n# CM:%d",cm);

  /* reading alphas (for bootrep) */
  if(fname_vt) {
      fp=openfp(fname_vt,fext_vt,"r",&cbuf);
      printf("\n# reading %s",cbuf);
      nalpha=0; alphavec=fread_vec(fp,&nalpha);
      fclose(fp); FREE(cbuf);
  } else {
    nalpha=nalpha0; alphavec=alphavec0;
  }

  if(sw_doau) {
    /* multie-scale bootstrap */
    do_bootrep();

    /* output results */
    write_pv(0,1);
    write_ci();
    if(sw_outrep) write_rep();
    if(sw_outcnt) write_cnt();
  }

  return 0;
}

int do_cntmode()
{
  int i;
  char *cbuf;
  FILE *fp;

  if(fname_cnt){
    fp=openfp(fname_cnt,fext_cnt,"r",&cbuf);
    printf("\n# reading %s",cbuf);
  } else {
    printf("\n# reading from stdin"); fp=STDIN;
  }
  kk=cm=0;
  orderv=fread_ivec(fp,&cm); obsvec=fread_vec(fp,&cm);
  rr=fread_vec(fp,&kk); bb=fread_ivec(fp,&kk);
  cntmat=fread_mat(fp,&cm,&kk);
  if(fname_cnt) {fclose(fp); FREE(cbuf);}

  printf("\n# K:%d",kk);
  printf("\n# R:"); for(i=0;i<kk;i++) printf("%g ",rr[i]);
  printf("\n# B:"); for(i=0;i<kk;i++) printf("%d ",bb[i]);
  printf("\n# CM:%d",cm);


  if(sw_doau) {
    /* multie-scale bootstrap */
    do_bootcnt();
    
    /* output results */
    write_pv(0,1);
  }

  return 0;
}


/*
  calculate p-values of conventional methods

  cm-vector of p-values (se's)
  (1) npvec (nsvec) : naive bootstrap
  (2) khpvec (khsvec) : Kishino-Hasegawa test
  (3) mcpvec (mcsvec) : Multiple Comparisons test
  (4) mspvec (mssvec) : MC-test with standardization weight

 */
int do_mctest()
{
  double *buf,x,*xp,*yp;
  int i,j,ib,i0;
  int nb;
  double **repmat; /* replicates of data */
  double **statmat; /* replicates of statistics */
  double *mnvec; /* mean of the replicates */
  double **wtmat; /* weight for ms-pvec */

  printf("\n# MC-TEST STARTS");
  
  /* find i0 s.t. rr[i0] is close to 1.0 */
  buf=new_vec(kk);
  for(i=0;i<kk;i++) buf[i]=fabs(rr[i]-1.0);
  i0=argmin_vec(buf,kk); nb=bb[i0];
  printf("\n# I0:%d R0:%g B0:%d",i0+1,rr[i0],bb[i0]);

  /* alloc memory */
  repmat=new_mat(mm,nb); statmat=new_mat(cm,nb);
  npvec=new_vec(cm); khpvec=new_vec(cm); 
  mcpvec=new_vec(cm); mspvec=new_vec(cm);
  mnvec=new_vec(mm); wtmat=new_mat(mm,mm);

  /* copy repmats[i0] to repmat */
  for(i=0;i<mm;i++) for(j=0;j<nb;j++) repmat[i][j]=repmats[i0][i][j];

  /* calculate the naive p-value */
  printf("\n# calculating the naive p-values");
  repminmaxs(repmat,statmat,mm,nb,cm,assvec,asslen,buf1,buf2);
  for(j=0;j<cm;j++) npvec[j]=1.0-pvaldist(statmat[j],nb,threshold);
  nsvec=getseval(npvec,cm,nb,NULL);
  
  /* centering the replicates */
  printf("\n# centering the replicates");
  for(j=0;j<mm;j++) {
    x=0.0; for(i=0;i<nb;i++) x+=repmat[j][i];
    x=mnvec[j]=x/nb;
    for(i=0;i<nb;i++) repmat[j][i]-=x;
  }

  /* calculate kh-pvalue */
  printf("\n# calculating kh-pvalue");    
  for(i=0;i<cm;i++) {
    khpvec[i]=calckhpval(datvec,repmat,mm,nb,
			 assvec[i],cassvec[i],asslen[i],NULL);
    putdot();
  }
  khsvec=getseval(khpvec,cm,nb,NULL);
  
  /* calculate mc-pvalue */
  printf("\n# calculating mc-pvalue");  
  for(i=0;i<cm;i++) {
    mcpvec[i]=calcmcpval(datvec,repmat,mm,nb,
			 assvec[i],cassvec[i],asslen[i],NULL);
    putdot();
  }
  mcsvec=getseval(mcpvec,cm,nb,NULL);

  /* calculate the weights */
  printf("\n# calculating the variances");
  for(i=0;i<mm;i++) {
     wtmat[i][i]=0.0; xp=repmat[i];
    for(j=0;j<i;j++) {
      yp=repmat[j]; x=0.0;
      for(ib=0;ib<nb;ib++) x+=(xp[ib]-yp[ib])*(xp[ib]-yp[ib]);
      wtmat[i][j]=wtmat[j][i]=1.0/sqrt((x+varadd)/(nb-1));
    }
    putdot();
  }

  /* calculate ms-pvalue */
  printf("\n# calculating ms-pvalue");  
  for(i=0;i<cm;i++) {
    mspvec[i]=calcmcpval(datvec,repmat,mm,nb,
			 assvec[i],cassvec[i],asslen[i],wtmat);
    putdot();
  }
  mssvec=getseval(mspvec,cm,nb,NULL);

  printf("\n# MC-TEST DONE");

  free_vec(buf); free_mat(repmat); free_mat(statmat);
  free_vec(mnvec); free_mat(wtmat);
  return 0;
}

int do_bootrep()
{
  int i,j,k;
  double **statp,*beta,**vmat,x;

  printf("\n# AU-TEST STARTS");

  /* memory alloc for pv */
  pvvec=new_vec(cm); sevec=new_vec(cm);
  pv0vec=new_vec(cm); se0vec=new_vec(cm);
  betamat=new_mat(cm,2); 
  rssvec=new_vec(cm); dfvec=new_vec(cm); pfvec=new_vec(cm);
  thvec=new_vec(cm); 
  statp=NEW_A(kk,double*);

  /* memory alloc for ci */
  cimat=new_mat(cm,nalpha); csmat=new_mat(cm,nalpha);
  eimat=new_mat(cm,nalpha);

  /* sorting the replicates */
  printf("\n# sorting the replicates");    
  for(i=0;i<kk;i++) {
    for(j=0;j<cm;j++) sort_vec(statmats[i][j],bb[i]);
    putdot();
  }

  /* calculate au-pvalues */
  printf("\n# calculating approximately unbiased p-values");
  for(i=0;i<cm;i++) {
    dprintf(1,"\n# rank=%d item=%d",i+1,orderv[i]+1);
    for(k=0;k<kk;k++) statp[k]=statmats[k][i];
    j=vcalpval(statp,rr,bb,kk,threshold,thvec+i,
	       pvvec+i,sevec+i,rssvec+i,dfvec+i,&beta,&vmat);
    dprintf(1,"\n# ret=%d",j);
    if(vmat==NULL) {
      pfvec[i]=0.0;
      warning("regression degenerated: df[%d]=%g",i+1,dfvec[i]);
      betamat[i][0]=betamat[i][1]=0.0;
      pv0vec[i]=se0vec[i]=0.0;
    } else {
      pfvec[i]=pochisq(rssvec[i],(int)(dfvec[i]));
      if(pfvec[i]<0.001)
	warning("theory does not fit well: pfit[%d]=%.4f",i+1,pfvec[i]);
      betamat[i][0]=beta[0]; betamat[i][1]=beta[1];
      pv0vec[i]=1.0-poz(beta[0]+beta[1]);
      x=dnorm(beta[0]+beta[1]);
      se0vec[i]=sqrt(x*x*(vmat[0][0]+vmat[1][1]+vmat[0][1]+vmat[1][0]));
    }
    /* make cntmat for output */
    if(sw_outcnt){
      for(k=0;k<kk;k++)
	cntmat[i][k]=cntdist(statmats[k][i],bb[k],threshold,
			     sw_smoothcnt?3:2);
    }
    putdot();
  }

  /* ci */
  printf("\n# ALPHA:");
  for(i=0;i<nalpha;i++) printf("%g ",alphavec[i]);
  printf("\n# calculating confidence intervals");
  for(j=0;j<cm;j++) {
    for(i=0;i<kk;i++) statp[i]=statmats[i][j];
    for(i=0;i<nalpha;i++) {
      dprintf(1,"\n# invtpval item=%d alpha=%g ",j+1,alphavec[i]);
      invtpval(statp,rr,bb,kk,alphavec[i],
	       &(cimat[j][i]),&(csmat[j][i]),&(eimat[j][i]));
    }
    putdot();
  }


  FREE(statp);
  printf("\n# AU-TEST DONE");

  return 0;
}

int do_bootcnt()
{
  int i,j;
  double *beta,**vmat,x;

  printf("\n# AU-TEST STARTS");

  /* memory alloc for pv */
  pvvec=new_vec(cm); sevec=new_vec(cm);
  pv0vec=new_vec(cm); se0vec=new_vec(cm);
  betamat=new_mat(cm,2); 
  rssvec=new_vec(cm); dfvec=new_vec(cm); pfvec=new_vec(cm);
  thvec=new_vec(cm); 

  for(i=0;i<cm;i++) {
    dprintf(1,"\n# rank=%d item=%d",i+1,orderv[i]+1);
    j=rcalpval(cntmat[i],rr,bb,kk,
	       pvvec+i,sevec+i,rssvec+i,dfvec+i,&beta,&vmat);
    thvec[i]=0.0; /* unused in cntmode */
    dprintf(1,"\n# ret=%d",j);
    if(vmat==NULL) {
      pfvec[i]=0.0;
      warning("regression degenerated: df[%d]=%g",i+1,dfvec[i]);
      betamat[i][0]=betamat[i][1]=0.0;
      pv0vec[i]=se0vec[i]=0.0;
    } else {
      pfvec[i]=pochisq(rssvec[i],(int)(dfvec[i]));
      if(pfvec[i]<0.001)
	warning("theory does not fit well: pfit[%d]=%.4f",i+1,pfvec[i]);
      betamat[i][0]=beta[0]; betamat[i][1]=beta[1];
      pv0vec[i]=1.0-poz(beta[0]+beta[1]);
      x=dnorm(beta[0]+beta[1]);
      se0vec[i]=sqrt(x*x*(vmat[0][0]+vmat[1][1]+vmat[0][1]+vmat[1][0]));
    }
  }

  printf("\n# AU-TEST DONE");  
  return 0;
}

#define MCPVNUM 4
#define AUPVNUM 2
#define MCAUXNUM 0
#define AUAUXNUM 6
int write_pv(int sw_mc, int sw_au)
{
  FILE *fp;
  char *cbuf;
  double **pvmat,**semat,**auxmat;
  int pvnum,auxnum,i,j;

  if(fname_pv) {
    fp=openfp(fname_pv,fext_pv,"w",&cbuf);
    printf("\n# writing %s",cbuf);
  } else {
    fp=STDOUT;
    printf("\n# writing to stdout\n");
  }

  fprintf(fp,"\n# ITEM:\n"); fwrite_ivec(fp,orderv,cm);
  fprintf(fp,"\n# STAT:\n"); fwrite_vec(fp,obsvec,cm);

  pvnum=(sw_mc?MCPVNUM:0)+(sw_au?AUPVNUM:0);
  pvmat=new_mat(cm,pvnum); semat=new_mat(cm,pvnum);
  for(i=0;i<cm;i++) {
    j=0;
    if(sw_mc) {
      pvmat[i][j]=npvec[i]; semat[i][j]=nsvec[i]; j++;
      pvmat[i][j]=khpvec[i]; semat[i][j]=khsvec[i]; j++;
      pvmat[i][j]=mcpvec[i]; semat[i][j]=mcsvec[i]; j++;
      pvmat[i][j]=mspvec[i]; semat[i][j]=mssvec[i]; j++;
    }
    if(sw_au){
      pvmat[i][j]=pvvec[i]; semat[i][j]=sevec[i]; j++;
      pvmat[i][j]=pv0vec[i]; semat[i][j]=se0vec[i]; j++;
    }
  }
  fprintf(fp,"\n# PV:\n"); fwrite_mat(fp,pvmat,cm,pvnum);  
  fprintf(fp,"\n# SE:\n"); fwrite_mat(fp,semat,cm,pvnum);  

  auxnum=(sw_mc?MCAUXNUM:0)+(sw_au?AUAUXNUM:0);
  auxmat=new_mat(cm,auxnum);
  for(i=0;i<cm;i++) {
    j=0;
    if(sw_mc) {
    }
    if(sw_au){
      auxmat[i][j]=pfvec[i]; j++;  /* pvalue of diagnostic */
      auxmat[i][j]=rssvec[i]; j++; /* rss */
      auxmat[i][j]=dfvec[i]; j++;  /* df */
      auxmat[i][j]=betamat[i][0]; j++; /* signed distance */
      auxmat[i][j]=betamat[i][1]; j++; /* curvature */
      auxmat[i][j]=thvec[i]; j++; /* threshold */
    }
  }
  fprintf(fp,"\n# AX:\n"); fwrite_mat(fp,auxmat,cm,auxnum);  

  if(fname_pv) {fclose(fp); FREE(cbuf);}
  return 0;
}

int write_ci()
{
  FILE *fp;
  char *cbuf;

  if(fname_ci) {
    fp=openfp(fname_ci,fext_ci,"w",&cbuf);
    printf("\n# writing %s",cbuf);
  } else {
    fp=STDOUT;
    printf("\n# writing to stdout\n");
  }
  fprintf(fp,"\n# ITEM:\n"); fwrite_ivec(fp,orderv,cm);
  fprintf(fp,"\n# STAT:\n"); fwrite_vec(fp,obsvec,cm);
  fprintf(fp,"\n# ALPHA:\n"); fwrite_vec(fp,alphavec,nalpha);

  fprintf(fp,"\n# CI:\n"); fwrite_mat(fp,cimat,cm,nalpha);
  fprintf(fp,"\n# SE:\n"); fwrite_mat(fp,csmat,cm,nalpha);
  fprintf(fp,"\n# EI:\n"); fwrite_mat(fp,eimat,cm,nalpha);

  if(fname_ci) {fclose(fp); FREE(cbuf);}
  return 0;
}

int write_rep()
{
  FILE *fp;
  char *cbuf;
  int i;

  if(fname_rep) { /* binary write to file */
    fp=openfp(fname_rep,fext_rep,"w",&cbuf);
    printf("\n# writing %s",cbuf);
    fwrite_bivec(fp,orderv,cm); fwrite_bvec(fp,obsvec,cm);
    fwrite_bvec(fp,rr,kk); fwrite_bivec(fp,bb,kk);
    fwrite_bi(fp,kk);
    for(i=0;i<kk;i++) {
      fwrite_bmat(fp,statmats[i],cm,bb[i]); putchar('.');
    }
    fclose(fp); FREE(cbuf);
  } else { /* ascii write to stdout */
    printf("\n# writing to stdout\n");
    write_ivec(orderv,cm); write_vec(obsvec,cm);
    write_vec(rr,kk); write_ivec(bb,kk);
    write_i(kk);
    for(i=0;i<kk;i++) {
      write_mat(statmats[i],cm,bb[i]);
    }
  }

  return 0;
}

int write_cnt()
{
  FILE *fp;
  char *cbuf;
  
  if(fname_cnt) {
    fp=openfp(fname_cnt,fext_cnt,"w",&cbuf);
    printf("\n# writing %s",cbuf);
  } else {
    fp=STDOUT;
    printf("\n# writing to stdout\n");
  }
  fprintf(fp,"\n# ITEM:\n"); fwrite_ivec(fp,orderv,cm);
  fprintf(fp,"\n# STAT:\n"); fwrite_vec(fp,obsvec,cm);
  fprintf(fp,"\n# R:\n"); fwrite_vec(fp,rr,kk);
  fprintf(fp,"\n# B:\n"); fwrite_ivec(fp,bb,kk);
  fprintf(fp,"\n# CNT:\n"); fwrite_mat(fp,cntmat,cm,kk);

  if(fname_cnt) {fclose(fp); FREE(cbuf);}
  return 0;
}


/*
 *  COUNTING ROUTINES
 */

/* returns #{vec[i] <= t} */
int cntdist1(double *vec, int bb, double t)
{
  int i,x;
  for(x=i=0;i<bb;i++) if(vec[i]<=t) x++;
  return x;
}

/* binary search for a sorted vector 
   find k s.t. vec[k-1] <= t < vec[k]
 */
int cntdist2(double *vec, int bb, double t)
{
  int i,i0,i1;

  i0=0; i1=bb-1;
  if(t < vec[0]) return (double)0;
  else if(vec[bb-1] <= t) return (double)bb;

  while(i1-i0>1) {
    i=(i0+i1)/2;
    if(vec[i] <= t) i0=i;
    else i1=i;
  }

  return i1;
}
/*
  smoothing the counting for a sorted vector
  the piecewise linear function connecting
  F(v[i]) =  1/(2n) + i/n, for i=0,...,n-1
  F(1.5v[0]-0.5v[1]) = 0
  F(1.5v[n-1]-0.5v[n-2]) = 1.

  1. F(x)=0 for x<=1.5v[0]-0.5v[1] 

  2. F(x)=1/(2n) + (1/n)*(x-v[0])/(v[1]-v[0])
  for 1.5v[0]-0.5v[1] < x <= v[0]

  3. F(x)=1/(2n) + i/n + (1/n)*(x-v[i])/(v[i]-v[i+1])
  for v[i] < x <= v[i+1], i=0,...,

  4. F(x)=1-(1/2n) + (1/n)*(x-v[n-1])/(v[n-1]-v[n-2])
  for v[n-1] < x <= 1.5v[n-1]-0.5v[n-2]

  5. F(x)=1 for x > 1.5v[n-1]-0.5v[n-2]
 */
double cntdist3(double *vec, int bb, double t)
{
  double p,n;
  int i;
  i=cntdist2(vec,bb,t)-1; /* to find vec[i] <= t < vec[i+1] */
  n=(double)bb;
  if(i<0) {
    if(vec[1]>vec[0]) p=0.5+(t-vec[0])/(vec[1]-vec[0]);
    else p=0.0;
  } else if(i<bb-1) {
    if(vec[i+1]>vec[i]) p=0.5+(double)i+(t-vec[i])/(vec[i+1]-vec[i]);
    else p=0.5+(double)i; /* <- should never happen */
  } else {
    if(vec[bb-1]-vec[bb-2]>0) p=n-0.5+(t-vec[bb-1])/(vec[bb-1]-vec[bb-2]);
    else p=n;
  }
  if(p>n) p=n; else if(p<0.0) p=0.0;
  return p;
}

double cntdist(double *vec, int bb, double t, int modesw)
{
  /* vec must be sorted for modesw=2,3 */
  switch(modesw) {
  case 1: return (double)cntdist1(vec,bb,t);
  case 2: return (double)cntdist2(vec,bb,t);
  case 3: return cntdist3(vec,bb,t);
  default: 
    sort_vec(vec,bb);
    return cntdist3(vec,bb,t);
  }
}

int sw_pvaldist=1; /* smoothing in pvaldist */
double pvaldist(double *vec, int bb, double t)
{
  double x;

  if(sw_pvaldist) x=cntdist(vec,bb,t,4); /* smoothing */
  else x=cntdist(vec,bb,t,1); /* nonsmoothing */

  return (1.0-x/bb);
}

/*
 *  MC-TEST ROUTINES
 */

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

void repminmaxs(double **repmat,  /* m x bb replicate mat */
		double **statmat, /* cm x bb output */
		int m,   /* number of itmes for input */
		int bb,  /* number of replicates */
		int cm,  /* number of itmes for output */
		int **assvec, /* ass vectors of size cm */
		int *asslen, /* length of each assvec[i] */
		double *buf1, /* buf of size m */
		double *buf2  /* buf of size cm */
		)
{
  int i,j;

  for(i=0;i<bb;i++) {
    for(j=0;j<m;j++) buf1[j]=repmat[j][i];
    calcmaxs(buf1,m);
    calcassmins(buf1,assvec,asslen,cm,buf2);
    for(j=0;j<cm;j++)
      statmat[j][i]=buf2[j];
  }
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
  datvec = mm-vector of data
  repmat = mm x nb matrix of centered replicates
  ass = alen-vector of association
  cass = (mm-alen)-vector of complements
  wt = mm x mm matrix of weights
 */
double calcmcpval(double *datvec, double **repmat, int mm, int nb,
		  int *ass, int *cass, int alen, double **wt)
{
  int i,j;
  double x;
  static int alen0=0, nb0=0, mm0=0;
  static double *pvec=NULL, *buf=NULL, **pmat=NULL, *stat=NULL;

  if(nb>nb0) {
    for(i=0;i<alen0;i++) pmat[i]=renew_vec(pmat[i],nb);
    nb0=nb;
  }
  if(alen>alen0) {
    pvec=renew_vec(pvec,alen); 
    stat=renew_vec(stat,alen); 
    pmat=RENEW_A(pmat,alen,double*);
    for(i=alen0;i<alen;i++) pmat[i]=new_vec(nb0);
    alen0=alen;
  }
  if(mm>mm0) {
    buf=renew_vec(buf,mm); 
    mm0=mm;
  }

  for(i=0;i<nb;i++) {
    for(j=0;j<mm;j++) buf[j]=repmat[j][i];
    calcmaxass(buf,wt,mm,ass,cass,alen,stat);
    for(j=0;j<alen;j++) pmat[j][i]=stat[j];
  }

  calcmaxass(datvec,wt,mm,ass,cass,alen,stat);
  for(i=0;i<alen;i++) pvec[i]=pvaldist(pmat[i],nb,stat[i]);

  x=0.0;
  for(i=0;i<alen;i++) if(pvec[i]>x) x=pvec[i];

  return x;
}

/*
  datvec = mm-vector of data
  repmat = mm x nb matrix of centered replicates
  ass = alen-vector of association
  cass = (mm-alen)-vector of complements
  wt = mm x mm matrix of weights
 */
double calckhpval(double *datvec, double **repmat, int mm, int nb,
		  int *ass, int *cass, int alen, double **wt)
{
  int i,j,ia,ja,i0,j0,jj;
  double x,y,x0;
  static int alen0=0, nb0=0;
  static double *stat=NULL, *repvec=NULL;

  if(nb>nb0) {repvec=renew_vec(repvec,nb); nb0=nb;}
  if(alen>alen0) {stat=renew_vec(stat,alen); alen0=alen;}

  /* find the "contrast" that is effective (min-max) */
  if(wt) {
    x0=HUGENUM;
    for(i=0;i<alen;i++) {
      ia=ass[i]; x=-HUGENUM;
      for(j=0;j<mm-alen;j++) {
	ja=cass[j]; y=wt[ia][ja]*(datvec[ja]-datvec[ia]);
	if(y>x) {x=y; jj=j;} /*max*/
      }
      if(x<x0) {x0=x; i0=i; j0=jj;} /*min*/
    }
  } else {
    x=-HUGENUM;
    for(j=0;j<mm-alen;j++) {
      y=datvec[cass[j]]; if(y>x) {x=y; j0=j;} /*max*/
    }
    x=HUGENUM;
    for(i=0;i<alen;i++) {
      y=-datvec[ass[i]]; if(y<x) {x=y; i0=i;} /*min*/
    }
  }

  ja=cass[j0]; ia=ass[i0]; x0=datvec[ja]-datvec[ia];
  for(i=0;i<nb;i++) repvec[i]=repmat[ja][i]-repmat[ia][i];
  x=pvaldist(repvec,nb,x0);
  return x;
}

double *getseval(double *pval, int m, int nb, double *seval)
{
  int i;
  if(seval==NULL) seval=new_vec(m);
  for(i=0;i<m;i++) seval[i]=sqrt(pval[i]*(1.0-pval[i])/nb);

  return seval;
}

/*
 *  MS-BOOT ROUTINES
 */

double zvaleps=1.0e-6;
double zvalbig=6.0;
int calcpval(double *cnts, double *rr, int *bb, int kk,
	     double *pv, double *se, double *rss, double *df, 
	     double **betap, double ***vmatp, /* reference only */
	     double rmin, double rmax)
{
  int i;
  double x,s,**vmat;
  static double beta[2];
  static int kk0=0;
  static double *npv=NULL,*zval=NULL,*wt=NULL; /* kk-vector */
  static int *unused=NULL;
  static double **X=NULL; /* 2 x kk matrix */
  int m; /* # of used items */

  if(kk>kk0){
    X=renew_mat(X,2,kk);
    npv=renew_vec(npv,kk);
    zval=renew_vec(zval,kk);
    wt=renew_vec(wt,kk);    
    unused=renew_ivec(unused,kk);
    kk0=kk;
  }

  m=0;
  for(i=0;i<kk;i++) {
    s=sqrt(1.0/rr[i]);
    X[0][i]=1.0/s; X[1][i]=s;
    npv[i]=cnts[i]/bb[i];
    if((rr[i]<rmin) || (rr[i]>rmax)) {
      unused[i]=1; zval[i]=0.0; wt[i]=0.0;}
    else if(npv[i]<zvaleps) {
      unused[i]=1; zval[i]=zvalbig; wt[i]=0.0;}
    else if(npv[i]>1.0-zvaleps) {
      unused[i]=1; zval[i]=-zvalbig; wt[i]=0.0;}
    else {
      unused[i]=0; zval[i]=critz(1.0-npv[i]); 
      x=dnorm(zval[i]);
      wt[i]=bb[i]*x*x/((1.0-npv[i])*npv[i]);
      m++;
    }
  }
  *df=(double)(m-2);

  if(m<2) {
    for(x=0.0,i=0;i<kk;i++) x+=zval[i];
    if(x>=0.0) *pv=0.0; else *pv=1.0;
    *se=0.0; *rss=0.0;
    beta[0]=beta[1]=0.0; vmat=NULL;
    i=1; /* degenerate */
  } else {
    lsfit(X,zval,wt,2,kk,beta,rss,&vmat);
    *pv=1.0-poz(beta[0]-beta[1]);
    x=dnorm(beta[0]-beta[1]);
    *se=sqrt(x*x*(vmat[0][0]+vmat[1][1]-vmat[0][1]-vmat[1][0]));
    i=0; /* non-degenerate */
  }

  if(betap) *betap=beta;  if(vmatp) *vmatp=vmat; 
  return i;
}


double pfmin=0.001; /* minimum p-value of diagnostics */
double rrmin=0.0; /* minimum rr */
double rrmax=100.0; /* maximum rr */
int rcalpval(double *cnts, double *rr, int *bb, int kk,
	     double *pv, double *se, double *rss, double *df, 
	     double **betap, double ***vmatp) /* reference only */
{
  static int kk0=0;
  static double *rr0=NULL;
  int i,i1,i2,j;
  double pfit;

  if(kk>kk0){
    rr0=renew_vec(rr0,kk);
    kk0=kk;
  }
  for(i=0;i<kk;i++) rr0[i]=rr[i];
  sort_vec(rr0,kk);
  /* find the minimum rr[i] which is not less than rmin */
  for(i=0;i<kk;i++) if(rr0[i]>=rrmin) break;
  i1=i;
  /* find the maximum rr[i] which is not less than 1.0 */
  for(;i<kk;i++) if(rr0[i]>=1.0) break;
  i2=i;
  /* call calcpval until the theory fits well */
  for(i=i1;i<i2;i++) {
    j=calcpval(cnts,rr,bb,kk,pv,se,rss,df,betap,vmatp,rr0[i],rrmax);
    if(j) break;
    pfit=pochisq(*rss,(int)(*df));
    if(pfit>=pfmin) break;
  }

  return j;
}


#define VCLOOPMAX 30
#define VCSCALE 0.5     
#define VCEPS1 0.01
#define VCEPS2 0.1
int dfmin=1;
int vcalpval(double **statps, double *rr, int *bb, int kk,
	     double threshold, double *thp,
	     double *pvp, double *sep, double *rssp, double *dfp, 
	     double **betap, double ***vmatp) /* reference only */
{
  int i,j,idf,idf0;
  double x,xn,pv,se,rss,df,z0,z,ze0,ze,th;
  double *beta, **vmat;
  static double *cntvec=NULL, *buf=NULL;
  static int kk0=0;

  if(kk>kk0){
    cntvec=renew_vec(cntvec,kk);
    buf=renew_vec(buf,kk);
    kk0=kk;
  }

  /* get approximate qunatile for the initial values */
  for(i=0;i<kk;i++) buf[i]=fabs(rr[i]-1.0);
  i=argmin_vec(buf,kk); xn=statps[i][bb[i]/2];

  *pvp=*sep=*rssp=*dfp=*thp=0.0;*betap=NULL;*vmatp=NULL;
  z=ze=0.0; th=threshold;
  /* iteration */
  for(j=0;j<VCLOOPMAX;j++) {
    x=xn;
    /* get cnts */
    for(i=0;i<kk;i++) cntvec[i]=cntdist(statps[i],bb[i],x,3);
    /* get pvalue */
    i=calcpval(cntvec,rr,bb,kk,&pv,&se,&rss,&df,&beta,&vmat,rrmin,rrmax);
    idf=df;
    if(vmat) {
      z=beta[0]-beta[1];
      ze=sqrt(vmat[0][0]+vmat[1][1]-vmat[0][1]-vmat[1][0]);      
    }
    dprintf(1,"\n# %3d %6.3f %6.3f %9.6f %9.6f %3.0f %6.3f %6.3f %6.3f %6.3f",
	    j,th,x,pv,se,df,beta[0],beta[1],z,ze);

    if(idf < dfmin && idf0 < dfmin) return 1; /* degenerated */
    if((idf < dfmin) || 
       (idf0 >= dfmin && (z-z0)*(x-*thp) > 0.0 && fabs(z-z0)>VCEPS2*ze0) ) {
      dprintf(1,"\n# non-monotone");
      th=x; xn=VCSCALE*x+(1.0-VCSCALE)*(*thp);
      continue;
    }
    if(idf0 >= dfmin && (fabs(z-z0)<VCEPS1*ze0)) {
      if(th==threshold) xn=th; else th=x;
    } else xn=VCSCALE*th+(1.0-VCSCALE)*x;
    *pvp=pv;*sep=se;*rssp=rss;*dfp=df;*betap=beta;*vmatp=vmat;*thp=x;
    z0=z;ze0=ze; idf0=idf;
    if(x==th) return 0;
  }
  
  dprintf(1,"\n# no-convergence");
  return 2;
}


/*
  statps[i] must be sorted
 */
#define CILOOPMAX 30
double cieps=0.1;
int invtpval(double **statps, double *rr, int *bb, int kk,
	     double alpha, 
	     double *cival, double *seval, double *eival)
{
  int i,j;
  double x1,x2,x,x0,e,pv,pv0,se,se0,rss,df,*beta;
  static double *cntvec=NULL, *buf=NULL;
  static int kk0=0;

  if(alpha<0.0) alpha=0.0; else if(alpha>1.0) alpha=1.0;

  if(kk>kk0){
    cntvec=renew_vec(cntvec,kk);
    buf=renew_vec(buf,kk);
    kk0=kk;
  }

  /* get approximate qunatile for the initial values */
  for(i=0;i<kk;i++) buf[i]=fabs(rr[i]-1.0);
  i=argmin_vec(buf,kk);
  x1=statps[i][0]; x2=statps[i][bb[i]-1];

  /* iteration */
  x=pv=0.0; se=1.0;
  for(j=0;j<CILOOPMAX;j++) {
    x0=x; x=0.5*(x1+x2);
    dprintf(2,"\n## %d %g ",j,x);
    /* get cnts */
    for(i=0;i<kk;i++) cntvec[i]=cntdist(statps[i],bb[i],x,3);
    /* get pvalue */
    se0=se; pv0=pv;
    i=calcpval(cntvec,rr,bb,kk, &pv,&se,&rss,&df,&beta,NULL,rrmin,rrmax);
    e=pv-alpha;
    dprintf(2,"%g %g %g %g %g %g %g %g ",
	    pv,se,rss,df,beta[0],beta[1],e,e/se);
    /* check */
    if((pv0-alpha) * (pv-pv0) > 0.0)
      warning("non-monotone in pval (%g %g) -> (%g %g)",x0,pv0,x,pv);
    if(i) {
      warning("degenerecy in pval df=%g",df);
      *cival=x0; *seval=0.0; *eival=0.0;
      return 1;
    }
    if((j>0) && (fabs(e) < se*cieps)) break;
    if(e>0.0) x2=x; else x1=x; /* update */
  }
  if(j==CILOOPMAX) warning("no convergence: ei=%g",e/se);

  *cival=x;
  *seval=se*(x-x0)/(pv-pv0);
  *eival=e/se;
  return 0;
}
