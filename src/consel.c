/*

  consel.c : assessing the confidence in selection
             using the multi-scale bootstrap

  Time-stamp: <2011-05-12 16:01:09 shimo>

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

static const char rcsid[] = "$Id: consel.c,v 1.20 2011/05/12 07:23:08 shimo Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rand.h"
#include "misc.h"
#include "opt.h"

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
void repminmaxsrs(double **repmat, double **statmat, 
		  int m, int bb, int cm, 
		  int **assvec, int *asslen, double *buf1, double *buf2,
		  double *repmean, double rs);
int compdvec(dvec *xp, dvec *yp);
int *getcass(int *ass, int n, int m, int *cass);
double *calcmaxass(double *xx, double **wt, int m, 
		   int *ass, int *cass, int n, double *out);
double pvaldist(double *vec, int bb, double t);
double calcmcpval(double *datvec, double **repmat, int mm, int nb,
		  int *ass, int *cass, int alen, double **wt);
double calckhpval(double *datvec, double **repmat, int mm, int nb,
		  int *ass, int *cass, int alen, double **wt);
double *getseval(double *pval, int m, double nb, double *seval);

/* multiscale bootstrap routines */
int rcalpval(double *cnts, double *rr, double *bb, int kk,
	     double *pv, double *se, double *pv0, double *se0, 
	     double *rss, double *df, 
	     double **betap, double ***vmatp, double kappa);
int vcalpval(double **statps, double *rr, int *bb, int kk,
	     double threshold, double *thp,
	     double *pvp, double *sep, double *pv0p, double *se0p, 
	     double *rssp, double *dfp, 
	     double **betap, double ***vmatp, double kappa);
int invtpval(double **statps, double *rr, int *bb, int kk,
	     double alpha, 
	     double *cival, double *seval, double *eival, double kappa);

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
char *fname_pa = NULL; char *fext_pa = ".pa";
char *fname_svt = NULL; char *fext_svt = ".svt";

/* for the scales */
int kk=0; /* no. of scales */
int *bb=NULL; /* kk-vector of no.'s of replicates */
double *rr=NULL; /* kk-vector of relative sample sizes */

/* for associations */
int mm=0; /* no. of items, i.e., no. of trees */
int cm0=0; /* no. of associations, edges typically (read) */
int **assvec0=NULL; /* association vectors (read) */
int *asslen0=NULL; /* lengths of the associations (read) */
int **cassvec0=NULL; /* complement of assvec0 */
int cm=0; /* truncate small items and use only top cm itmes */
int **assvec=NULL; /* sorted assvec */
int *asslen=NULL; /* sorted asslen */
int **cassvec=NULL; /* complement of assvec used in mctest */
int *orderv=NULL; /* cm-vector of the item id's */


/* chidim */
double chidimarg=0.0; /* default dimension */
int chinumgrid=10; /* number of grids */
int sw_chidimopt=1; /* optimize dimension */
int sw_chidimgs=1; /* use golden section search */
double chieps=1e-50;
int chiloopmax=50;
int chibackmax=10;
double chidimmin=2.0;

/* for multiple rmt */
int mf=1;

/* for rmt file */
double ***repmats=NULL; /* (kk,mm,bb[i])-array of replicates */
double *datvec=NULL; /* mm-vector of data */
double *datvec00=NULL; /* save original */

/* for rep file */
double ***statmats=NULL; /* (kk,cm,bb[i])-array of test statistics */
double *obsvec=NULL; /* cm-vector of observed statistics */

/* for cnt file */
double **rrmat=NULL;
double **bbmat=NULL;
double **cntmat=NULL; /* (cm,kk)-matrix of counts */

/* misc */
double threshold=0.0; /* the region is defined as x <= thereshold */
double varadd=1.0; /* adding to wtmat */
#define VARADDWARN 100.0 /* warning if var < varadd*VARADDWARN */
#define HUGENUM 1.0e30;
double kappa=1.0; /* weight for the curvature */
double vceps2;
double mleeps;
double congbpadd;

/* mspar */
int kk0=10;
double rr0[]={0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4};
int kk1=0;
double *rr1=NULL;
int *bb1=NULL;

/* for mc-tests: p-values and se (cm-vector) */
double *npvec, *nsvec; /* naive p-value */
double *mcpvec, *mcsvec; /* p-value of nonweighted mc-test */
double *mcwpvec, *mcwsvec; /* p-value of weighted mc-test */
double *khpvec, *khsvec; /* kh-test p-value */
double *khwpvec, *khwsvec; /* weighted kh-test p-value */

/* for msboot: p-values and se (cm-vector) */
double *pvvec, *sevec ; /* multi-scale pvalue */
double *pv0vec, *se0vec; /* pvalue for the derived naive method */

/* bayes posterior prob */
double *bapvec=NULL;
double bapcoef=1.0; // factor used for logL
double bapn=1.0; // sample size n used as log(n)=log(bapn)

/* for msboot pv */
double **betamat=NULL; /* (3,cm)-matrix of signed distance, curvature, dim */
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
double cieps;
/* au */
double **cimat=NULL; /* (cm,nalpha)-mat of confidence limits */
double **csmat=NULL; /* se of cimat */
double **eimat=NULL; /* error index of cimat */
/* bp */
double **ci0mat=NULL; /* (cm,nalpha)-mat of confidence limits */
double **cs0mat=NULL; /* se of cimat */
double **ei0mat=NULL; /* error index of cimat */

/* work */
double *buf1; /* mm-vector */
double *buf2; /* cm0-vector */
double *buf3; /* mm-vector */

/* switches */
int sw_pvaldist;
int sw_outrep=0; /* output rep file */
int sw_outcnt=0; /* output cnt file */
int sw_inrep=0; /* input rep file */
int sw_incnt=0; /* input cnt file */
int sw_smoothcnt=0; /* smooth cnt file */
int sw_doba=1; /* dont skip ba-test */
int sw_domc=1; /* dont skip mc-tests */
int sw_doau=1; /* dont skip au-tests */
int sw_dobp=1; /* dont skip bp-tests */
int sw_nosort=0; /* dont sort the items */
#define FITMODE_WLS 0
#define FITMODE_MLE 1
#define FITMODE_CHI 2
int sw_fitmode=FITMODE_MLE; /* use mle DEFAULT*/
int sw_fastrep=0; /* rescaling from r=1 */
int sw_multi=0; /* multiple input files */
int sw_cong=1; /* overall congruence */
int find_i0(double *rr, int kk);

/* AIC correction added on 2010/01/26 
pdm = the number of parameters for evolution model
*/
char *fname_pdm = NULL; char *fext_pdm = ".vt";
double *pdmvec=NULL; /* vector of pdm for each "tree" */
int npdm=0;
int aictype=2; /* 0=none, 1=ordinary, 2=new aic, 3=rell correction, 4=yet another aic */

/* main routines */
int do_rmtmode();
int do_repmode();
int do_cntmode();
int do_rmtmultimode();
/* sub routines */
int do_bptest();
int do_mctest();
int do_bootrep();
int do_bootcnt();
int do_batest();
/* output routines */
int write_pv(int sw_bp, int sw_ba, int sw_mc, int sw_au);
int write_ci();
int write_rep();
int write_cnt1();
int write_cnt2();

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
      case 1: fname_in=argv[i]; break;
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
    } else if(streq(argv[i],"--pdm")) {
      if(i+1>=argc) byebye();
      fname_pdm=argv[i+1];
      i+=1;
    } else if(streq(argv[i],"-p")) {
      if(i+1>=argc) byebye();
      fname_pa=argv[i+1];
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
    } else if(streq(argv[i],"--no_sort")) {
      sw_nosort=1;
    } else if(streq(argv[i],"--no_bp")) {
      sw_dobp=0;
    } else if(streq(argv[i],"--no_pp")) {
      sw_doba=0;
    } else if(streq(argv[i],"--no_sh")) {
      sw_domc=0;
    } else if(streq(argv[i],"--no_au")) {
      sw_doau=0;
    } else if(streq(argv[i],"--smooth_cnt")) {
      sw_smoothcnt=1;
    } else if(streq(argv[i],"--quick_mctest")) {
      sw_pvaldist=0; /* avoid sorting the replicates */
    } else if(streq(argv[i],"--th")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%lf",&threshold) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"--kappa")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%lf",&kappa) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"--chi_dim")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%lf",&chidimarg) != 1)
	byebye();
      sw_fitmode=FITMODE_CHI;
      i+=1;
    } else if(streq(argv[i],"--chi_dimmin")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%lf",&chidimmin) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"--chi_grid")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%d",&chinumgrid) != 1)
	byebye();
      sw_fitmode=FITMODE_CHI;
      i+=1;
    } else if(streq(argv[i],"--chi_opt")) {
      sw_chidimopt=1; sw_chidimgs=1;
      sw_fitmode=FITMODE_CHI;
    } else if(streq(argv[i],"--chi_noopt")) {
      sw_chidimopt=0; sw_chidimgs=0;
      sw_fitmode=FITMODE_CHI;
    } else if(streq(argv[i],"--chi_nogs")) {
      sw_chidimopt=1; sw_chidimgs=0;
      sw_fitmode=FITMODE_CHI;
    } else if(streq(argv[i],"--chieps")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%lf",&chieps) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"--ppcoef")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%lf",&bapcoef) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"--ppn")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%lf",&bapn) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"--aictype")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%d",&aictype) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"--cieps")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%lf",&cieps) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"--vceps2")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%lf",&vceps2) != 1)
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
    } else if(streq(argv[i],"--mle")) {
      sw_fitmode=FITMODE_MLE;
    } else if(streq(argv[i],"--wls")) {
      sw_fitmode=FITMODE_WLS;
    } else if(streq(argv[i],"--chi")) {
      sw_fitmode=FITMODE_CHI;
    } else if(streq(argv[i],"-f")) {
      sw_fastrep=1;
    } else if(streq(argv[i],"-g")) {
      sw_multi=1;
    } else if(streq(argv[i],"--mleeps")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%lf",&mleeps) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"--congbp")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%lf",&congbpadd) != 1)
	byebye();
      i+=1;
    } else byebye();
  }

  if(fname_in && !fname_out) fname_out=fname_in;
  if(fname_out) fname_out=rmvaxt(fname_out);
  if(sw_incnt+sw_inrep>1) error("only one of -C and -R can be specified");
  if(sw_incnt+sw_inrep==0) { /* rmt mode */
    fname_rmt=fname_in;
    fname_pv=fname_ci=fname_cnt=fname_rep=fname_out;
    if(sw_multi) do_rmtmultimode();
    else do_rmtmode();
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
  input: datvec, mm
*/
int read_assoc()
{
  int i,j;
  char *cbuf;
  FILE *fp;
  dvec **dvbuf;

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
  obsvec=new_vec(cm); orderv=new_ivec(cm0);
  assvec=NEW_A(cm,int*); asslen=NEW_A(cm,int); cassvec=NEW_A(cm,int*);

  /* sorting the association by the test statistics */
  dvbuf=NEW_A(cm0,dvec*);
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

  return 0;
}


int read_alpha()
{
  char *cbuf;
  FILE *fp;

  /* reading alphas (for bootrep) */
  if(fname_vt) {
      fp=openfp(fname_vt,fext_vt,"r",&cbuf);
      printf("\n# reading %s",cbuf);
      nalpha=0; alphavec=fread_vec(fp,&nalpha);
      fclose(fp); FREE(cbuf);
  } else {
    nalpha=nalpha0; alphavec=alphavec0;
  }

  return 0;
}

int read_scale()
{
  int i;
  char *cbuf;
  FILE *fp;

  /* reading parameters */
  if(fname_pa!=NULL) {
    fp=openfp(fname_pa,fext_pa,"r",&cbuf);
    printf("\n# reading %s",cbuf);
    kk1=0; rr1=fread_vec(fp,&kk1);
    fclose(fp); FREE(cbuf);
  } else {
    kk1=kk0; rr1=rr0;
  }
  printf("\n# RESCALING");
  printf("\n# K:%d",kk1);
  printf("\n# R:"); for(i=0;i<kk1;i++) printf("%g ",rr1[i]);

  return 0;
}


/* AIC correction 2011/01/26 

   for the theory see my note on 2011/01/25

   n = sample size for original data
   m = sample size for bootstrap data (R = m/n)
   pdm = number of model parameters
   datvec contains
       sum_{t=1}^n log p(x_t|\hat theta) 
   repmats contains
      (n/m) * sum_{t=1}^m log p(x*_t|\hat theta)

   We want to compute
   AIC(X) = AIC for the dataset
   AIC(X*) = AIC for a bootstrap dataset
   AIC'(X) and AIC'(X*) are those with a new correction term

   Note:
   Minimizing AIC value is equivalent to maximizing
         sum_{t=1}^n log p(x_t|\hat theta) - pdm.
   In the definitions of AIC below, an arbiterary constant factor 
   such as (1/n) is  multiplied to this value, which may be different
   from the standard definition of AIC such as
     -2 * { sum_{t=1}^n log p(x_t|\hat theta) - pdm }.

   --------------------------------------------------------
   aictype==0: do no aic correciton (comparing log-likelihoods)
   RECOMMENDED ONLY FOR TESTING TREES
   
   This method (aictype==0) is standard for testing trees.
   For comparing models with different number of parameters,
   AIC should be used instead of log-likleihdoods.

   --------------------------------------------------------
   aictype==1: a plug-in of multiscale bootstrap
   NOT RECOMMENDED FOR ANY SITUATIONS
   
   AIC value for sample size n is
     -2 * { sum_{t=1}^n log p(x_t|\hat theta) - pdm },
   and AIC for sample size m is
     -2 * { sum_{t=1}^m log p(x_t|\hat theta) - pdm }.

   Except for the arbiterary constant factor, AIC value for
   the original data and that for the bootstrap data are

   AIC(X) = -(1/n) *{
        sum_{t=1}^n log p(x_t|\hat theta) - 1.0*pdm  }

   AIC(X*) = -(1/n)*{
       (n/m)*sum_{t=1}^m log p(x*_t|\hat theta) - (n/m)*1.0*pdm }

   --------------------------------------------------------
   aictype==2: use a new aic correction
   RECOMMENDED FOR GENERAL SITUATIONS (when number of parameters
    are different for models).
   
   AIC'(X) = -(1/n) *{
        sum_{t=1}^n log p(x_t|\hat theta) - 0.5*pdm  }

   AIC'(X*) = -(1/n)*{
       (n/m)*sum_{t=1}^m log p(x*_t|\hat theta) - 0.5*pdm }
  
   --------------------------------------------------------
   Warnings: do not use aictype==3.
   aictype==3: use a new aic correction with rell taken into account
   NOT RECOMMNDED. IMPLEMENTED ONLY FOR EXPERIMENTAL PURPOSE
   
   AIC'(X) = -(1/n) *{
        sum_{t=1}^n log p(x_t|\hat theta) - 0.5*pdm  }

   AIC'(X*) = -(1/n)*{
       (n/m)*sum_{t=1}^m log p(x*_t|\hat theta) - 0.5*pdm* (1 - (n/m)) }
  
   note: AIC'(X*) above approximates
   AIC'(X*) = -(1/n)*{
       (n/m)*sum_{t=1}^m log p(x*_t|\hat theta*) - 0.5*pdm }
     
   the difference is between \hat theta and \hat theta*

   --------------------------------------------------------
   aictype==4: similar to aictype==2, but using AIC instead of RISK
   
   AIC'(X) = -(1/n) *{
        sum_{t=1}^n log p(x_t|\hat theta) - pdm  }

   AIC'(X*) = -(1/n)*{
       (n/m)*sum_{t=1}^m log p(x*_t|\hat theta) - pdm }
  
   --------------------------------------------------------
*/
int read_pdm()
{
  char *cbuf;
  FILE *fp;
  int i;

  /* reading pdm, i.e., number of parameters for each item */
  npdm=0;

  if(fname_pdm) {
      fp=openfp(fname_pdm,fext_pdm,"r",&cbuf);
      printf("\n# reading %s",cbuf);
      pdmvec=fread_vec(fp,&npdm);
      fclose(fp); FREE(cbuf);

      printf("\n# pdm:");
      for(i=0;i<npdm;i++) printf("%g ",pdmvec[i]);

      printf("\n# aictype:%d",aictype);
  }

  return 0;
}

int  do_pdm()
{
  int i,j,k;
  double x;

  if(npdm<=0) return 0;
  if(npdm != mm) error("size of pdm should be M");

  if(aictype==0) return 0;

  // update datvec
  for(i=0;i<mm;i++) {
    if(aictype==1) {
      x= -1.0*pdmvec[i];
    } else if(aictype==2) {
      x= -0.5*pdmvec[i];
    } else if(aictype==3) {
      x= -0.5*pdmvec[i];
    } else if(aictype==4) {
      x= -1.0*pdmvec[i];
    }
    datvec[i]+=x;
  }

  // update repmats
  for(k=0;k<kk;k++) {
    for(i=0;i<mm;i++) {
      if(aictype==1) {
	x= -1.0*pdmvec[i]/rr[k];
      }	else if(aictype==2) {
	x= -0.5*pdmvec[i];
      }	else if(aictype==3) {
	x= -0.5*pdmvec[i]*(1.0 - 1.0/rr[k]);
      }	else if(aictype==4) {
	x= -1.0*pdmvec[i];
      }
      for(j=0;j<bb[k];j++) repmats[k][i][j]+=x;
    }
  }

  return 0;
}

/*
  input1: cm,mm,assvec,asslen, rr,bb,kk,repmats,
  input2: sw_fastrep,rr1,kk1
  output1: statmats
  output2: rr,bb,kk
  work: buf1,buf2,buf3

  freeing repmats!!!
 */
int genrep()
{
  int i,j,i0;
  int nb;
  double x;

  if(sw_fastrep) { /* rescaling approximation */
    i0=find_i0(rr,kk); nb=bb[i0];
    bb1=new_ivec(kk1); for(i=0;i<kk1;i++) bb1[i]=nb;
    statmats=NEW_A(kk1,double**);
    for(i=0;i<kk1;i++) statmats[i]=new_lmat(cm,nb);
    for(j=0;j<mm;j++) {
      x=0.0; for(i=0;i<nb;i++) x+=repmats[i0][j][i];
      buf3[j]=x/nb;
    }
    for(i=0;i<kk1;i++) {
      repminmaxsrs(repmats[i0],statmats[i],mm,nb,cm,
		   assvec,asslen,buf1,buf2,buf3,
		   sqrt(rr[i0]/rr1[i]));
      putdot();
    }
    kk=kk1; rr=rr1; bb=bb1;
  } else { /* without rescaling approximation */
    // this is default
    statmats=NEW_A(kk,double**);
    for(i=0;i<kk;i++) {
      statmats[i]=new_lmat(cm,bb[i]);
      /* calculate the statistics for the replicates */
      repminmaxs(repmats[i],statmats[i],mm,bb[i],cm,
		 assvec,asslen,buf1,buf2);
      free_lmat(repmats[i],mm);
      putdot();
    }
  }

  return 0;
}


/*
  (1) input rmt file
  (2) read association file and sort the items
  (3) call mctest and bootrep

  THIS IS TYPICALLY USED FOR PHYLOGENETIC ANALYSIS

  pdm is added on 2010/01/26

*/
int do_rmtmode()
{
  int i;
  char *cbuf;
  FILE *fp;

  /* reading rmt */
  mm=kk=0;
  if(fname_rmt){ /* binary read from file */
    fp=openfp(fname_rmt,fext_rmt,"rb",&cbuf);
    printf("\n# reading %s",cbuf);
    datvec=fread_bvec(fp,&mm); 
    rr=fread_bvec(fp,&kk); bb=fread_bivec(fp,&kk); 
    i=fread_bi(fp); if(i != kk) error("wrong size in rmt");
    repmats=NEW_A(kk,double**);
    for(i=0;i<kk;i++) {
      repmats[i]=fread_blmat(fp,&mm,&(bb[i])); putdot();
    }
    fclose(fp); FREE(cbuf);
  } else { /* ascii read from stdin */
    printf("\n# reading from stdin");
    datvec=read_vec(&mm);
    rr=read_vec(&kk); bb=read_ivec(&kk);
    i=read_i(); if(i != kk) error("wrong size in rmt");
    repmats=NEW_A(kk,double**);
    for(i=0;i<kk;i++) {
      repmats[i]=read_lmat(&mm,&(bb[i])); putdot();
    }
  }
  datvec00=new_vec(mm); for(i=0;i<mm;i++) datvec00[i]=datvec[i]; 

  printf("\n# K:%d",kk); // number of scales
  printf("\n# R:"); for(i=0;i<kk;i++) printf("%g ",rr[i]); // scales
  printf("\n# B:"); for(i=0;i<kk;i++) printf("%d ",bb[i]); // no of reps
  printf("\n# M:%d",mm); // number of trees
  buf1=new_vec(mm); buf3=new_vec(mm); 

  if(kk<1) return 0;
  if(mm<2) {
    warning("M should be at least two");
    return 0;
  }

  /*  pdm correction for AIC 2011/01/26 */
  read_pdm();
  do_pdm();

  /* reading association vector */
  read_assoc();
  buf2=new_vec(cm0);

  /* reading other parameters */
    read_alpha();
  if(sw_fastrep) read_scale();

  
  /////// COMPUTATION STARTS HERE /////////

  /* bayes */
  if(sw_doba) do_batest();

  /* conventional methods */
  if(sw_domc) do_mctest();

  if(sw_doau||sw_dobp||sw_outrep||sw_outcnt) {
    printf("\n# calculate replicates of the statistics");
    genrep();
  }
  if(kk<2) sw_doau=0;

  /* naive bootstrap */
  if(sw_dobp) do_bptest();
  /* multi-scale bootstrap */
  if(sw_doau) do_bootrep();

  /* output results */
  if(sw_dobp||sw_domc||sw_doau||sw_doba)
    write_pv(sw_dobp,sw_doba,sw_domc,sw_doau);
  if(sw_doau) write_ci();
  if(sw_outrep) write_rep();
  if(sw_outcnt) write_cnt1();

  return 0;
}

/*
  input1: cm,mm,assvec,asslen, "file"
  input2: sw_fastrep,rr1,kk1
  output1: cntmat
  output2: rr,bb,kk
  work: buf1,buf2,buf3
 */
int gencnt(char *file)
{
  int i,j,i0;
  double x,**repmat,**statmat;
  char *cbuf;
  FILE *fp;
  int kk0, *bb0;
  double *rr0;

  /* reading rmt index */
  fp=openfp(file,fext_rmt,"rb",&cbuf);
  printf("\n# reading %s",cbuf);
  datvec=fread_bvec(fp,&mm); 
  kk=0; rr=fread_bvec(fp,&kk); bb=fread_bivec(fp,&kk); 
  i=fread_bi(fp); if(i != kk) error("wrong size in rmt");

  /* prepare */
  if(sw_fastrep) { /* rescaling approximation */
    i0=find_i0(rr,kk); 
    for(i=0;i<=i0;i++) {
      repmat=fread_blmat(fp,&mm,&(bb[i]));
      if(i!=i0) free_lmat(repmat,mm);
    }
    for(j=0;j<mm;j++) {
      x=0.0; for(i=0;i<bb[i0];i++) x+=repmat[j][i];
      buf3[j]=x/bb[i0];
    }
    bb1=new_ivec(kk1); for(i=0;i<kk1;i++) bb1[i]=bb[i0];
    kk0=kk; rr0=rr; bb0=bb;
    kk=kk1; rr=rr1; bb=bb1;
  }

  cntmat=new_mat(cm,kk);
  for(i=0;i<kk;i++) {
    statmat=new_lmat(cm,bb[i]);
    /* calc rep */
    if(sw_fastrep) {
      repminmaxsrs(repmat,statmat,mm,bb0[i0],cm,
		   assvec,asslen,buf1,buf2,buf3,
		   sqrt(rr0[i0]/rr[i]));
    } else {
      repmat=fread_blmat(fp,&mm,&(bb[i])); 
      repminmaxs(repmat,statmat,mm,bb[i],cm,
		 assvec,asslen,buf1,buf2);
    }
    putdot();
    /* calc cnt */
    for(j=0;j<cm;j++) {
      if(sw_smoothcnt) {
	sort_vec(statmat[j],bb[i]);
	cntmat[j][i]=cntdist(statmat[j],bb[i],threshold,3);
      } else {
	cntmat[j][i]=cntdist(statmat[j],bb[i],threshold,1);	
      }
    }
    if((!sw_fastrep) || (i==kk-1)) free_lmat(repmat,mm);
    free_lmat(statmat,cm);
  }

  fclose(fp); FREE(cbuf);
  return 0;
}

/*
  input1: cntmatm, bbm, mf
  input2: kk, cm
  output: cntmat, bbmat
 */
double congbpadd=1.0;
int congcnt(double ***cntmatm, int **bbm, int mf)
{
  int i,j,f,cm2;
  double x,y,z;
  double ***pvmatm,**pvmat,**vamat;

  /* alloc */
  cm2=cm+sw_cong;
  pvmatm=NEW_A(mf,double**);
  cntmat=new_mat(cm2,kk);
  bbmat=new_mat(cm2,kk);
  pvmat=new_mat(cm2,kk);
  vamat=new_mat(cm2,kk);

  /* calc bp with a correction */
  for(f=0;f<mf;f++){
    pvmatm[f]=new_mat(cm,kk);
    for(j=0;j<kk;j++) for(i=0;i<cm;i++) {
      x=cntmatm[f][i][j]; y=bbm[f][j];
      z=y*0.5+((y-2.0*congbpadd)/y)*(x-y*0.5);
      pvmatm[f][i][j]=z/y;
    }
  }
  for(j=0;j<kk;j++) for(i=0;i<cm;i++) {
    x=1.0; for(f=0;f<mf;f++) x*=pvmatm[f][i][j];
    pvmat[i][j]=x;
  }
  if(sw_cong) {
    for(j=0;j<kk;j++) {
      x=0.0; for(i=0;i<cm;i++) x+=pvmat[i][j];
      pvmat[cm][j]=x;
    }
  }

  /* calc var */
  for(j=0;j<kk;j++) {
    z=y=0.0; for(f=0;f<mf;f++) y+=1.0/bbm[f][j];
    for(i=0;i<cm;i++) {
      x=0.0;
      if(pvmat[i][j] > 0.0) {
	for(f=0;f<mf;f++) x+=1.0/(bbm[f][j]*pvmatm[f][i][j]);
      }
      x*=fsquare(pvmat[i][j]);
      z+=x;
      vamat[i][j]=x-y*fsquare(pvmat[i][j]);
    }
    if(sw_cong) vamat[cm][j]=z-y*fsquare(pvmat[cm][j]);
  }

  /* calc Beff */
  for(j=0;j<kk;j++) for(i=0;i<cm2;i++) {
    x=pvmat[i][j];
    bbmat[i][j]=x*(1.0-x)/vamat[i][j];
  }

  /* re-calc bp without the correction */
  for(f=0;f<mf;f++){
    for(j=0;j<kk;j++) for(i=0;i<cm;i++) 
      pvmatm[f][i][j]=cntmatm[f][i][j] / bbm[f][j];
  }
  for(j=0;j<kk;j++) for(i=0;i<cm;i++) {
    x=1.0; for(f=0;f<mf;f++) x*=pvmatm[f][i][j];
    pvmat[i][j]=x;
  }
  if(sw_cong) {
    for(j=0;j<kk;j++) {
      x=0.0; for(i=0;i<cm;i++) x+=pvmat[i][j];
      pvmat[cm][j]=x;
    }
  }

  /* calc cntmat */
  for(j=0;j<kk;j++) for(i=0;i<cm2;i++) {
    cntmat[i][j]=bbmat[i][j]*pvmat[i][j];
  }

  /* clean */
  free_mat(pvmat); free_mat(vamat);
  for(f=0;f<mf;f++) free_mat(pvmatm[f]);
  free(pvmatm);

  return 0;
}

/*
  (1) input rmt file
  (2) read association file and sort the items
  (3) call mctest and bootrep

  for multiple genes analysis
  (untouched since 2002)

  THIS IS YET EXPERIMENTAL

*/
int do_rmtmultimode()
{
  int i,j,f,*ip,cm2,i0;
  char *cbuf,**infiles;
  FILE *fp;
  double x;
  /* save */
  double *datvec0;
  /* pval */
  double **bapvecm;
  /* for multiple rmt file */
  double **datvecm=NULL; /* mf x mm-vector of data */
  int **bbm=NULL; /* mf x kk-vector of no.'s of replicates */
  double **rrm=NULL; /* mf x kk-vector of relative sample sizes */
  double ***cntmatm; /* mf x cm x kk matrix of counts */

  /* read file names */
  if(fname_rmt) fp=openfp(fname_rmt,fext_svt,"r",&cbuf); else fp=STDIN;
  mf=0; infiles=fread_svec(fp,&mf);
  if(fname_rmt) {fclose(fp); FREE(cbuf);}
  printf("\n# MF:%d",mf);

  /* pre-reading multi-rmt */
  datvecm=NEW_A(mf,double*); 
  rrm=NEW_A(mf,double*); bbm=NEW_A(mf,int*);
  mm=kk=0;
  for(f=0;f<mf;f++) {
    fp=openfp(infiles[f],fext_rmt,"rb",&cbuf);
    printf("\n# checking %s",cbuf);
    datvecm[f]=fread_bvec(fp,&mm);
    rrm[f]=fread_bvec(fp,&kk);
    bbm[f]=fread_bivec(fp,&kk); 
    printf("\n# R:"); for(i=0;i<kk;i++) printf("%6.2f ",rrm[f][i]);
    printf("\n# B:"); for(i=0;i<kk;i++) printf("%6d ",bbm[f][i]);
    for(i=0;i<kk;i++)
      if( fabs(rrm[f][i]-rrm[0][i])>rrm[0][i]*0.1) error("rr mismatch");
    i=fread_bi(fp); if(i != kk) error("wrong size in rmt");
    fclose(fp); FREE(cbuf);
  }
  printf("\n# K:%d",kk);
  printf("\n# M:%d",mm);
  if(kk<1) return 0;
  buf1=new_vec(mm); buf3=new_vec(mm); 

  /* compose total lik */
  datvec0=new_vec(mm);
  for(i=0;i<mm;i++) {
    x=0.0; for(f=0;f<mf;f++) x+= datvecm[f][i];
    datvec0[i]=x;
  }

  /* reading association vector */
  datvec=datvec0;
  read_assoc();
  buf2=new_vec(cm0);

  /* check if disjoint */
  ip=new_ivec(mm); for(i=0;i<mm;i++) ip[i]=0;
  for(i=0;i<cm;i++) for(j=0;j<asslen[i];j++) ip[assvec[i][j]]++;
  for(i=0;i<mm;i++) if(ip[i]>1) {
    printf("\n# non-disjoint");
    sw_cong=0; break;
  }
  free_ivec(ip);
  cm2=cm+sw_cong;

  /* copy obsvec and orderv */
  if(sw_cong) {
    obsvec=renew_vec(obsvec,cm2);
    orderv=renew_ivec(orderv,cm2);
    obsvec[cm]=0.0; orderv[cm]=cm;
  }

  /* reading parameters */
  read_alpha();
  if(sw_fastrep) read_scale();
  
  /* bayes */
  if(sw_doba) {
    bapvecm=NEW_A(mf,double*);
    for(f=0;f<mf;f++) {
      datvec=datvecm[f]; do_batest(); bapvecm[f]=bapvec;
    }
    bapvec=new_vec(cm2);
    for(i=0;i<cm;i++) {
      x=1.0;
      for(f=0;f<mf;f++) x*=bapvecm[f][i];
      bapvec[i]=x;
    }
    if(sw_cong) {
      x=0.0; for(i=0;i<mf;i++) x+=bapvec[i];
      bapvec[cm]=x;
    }
    free_lmat(bapvecm,mf);
  }

  /* reading multi-rmt files */
  if(sw_doau||sw_dobp||sw_outcnt) {
    cntmatm=NEW_A(mf,double**);
    for(f=0;f<mf;f++) {
      gencnt(infiles[f]);
      cntmatm[f]=cntmat; bbm[f]=bb; rrm[f]=rr;
    }
    rrmat=new_mat(cm2,kk);
    for(i=0;i<kk;i++) { /* averaging rr's */
      x=0.0; for(f=0;f<mf;f++) x+=rrm[f][i];
      x=x/mf; for(j=0;j<cm2;j++) rrmat[j][i]=x;
    }
    congcnt(cntmatm,bbm,mf);
  }


  /* change cm */
  cm=cm2;

  /* bp-test */
  if(sw_dobp) {
    npvec=new_vec(cm); nsvec=new_vec(cm);
    for(i=0;i<cm;i++) {
      i0=find_i0(rrmat[i],kk);
      npvec[i]=cntmat[i][i0]/bbmat[i][i0];
      nsvec[i]=sqrt(npvec[i]*(1.0-npvec[i])/bbmat[i][i0]);
    }
  }

  /* au-test */
  if(sw_doau) {
    do_bootcnt();
  }

  if(sw_dobp||sw_doau||sw_doba)
    write_pv(sw_dobp,sw_doba,0,sw_doau);
  if(sw_outcnt) write_cnt2();

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
    fp=openfp(fname_rep,fext_rep,"rb",&cbuf);
    printf("\n# reading %s",cbuf);
    orderv=fread_bivec(fp,&cm); obsvec=fread_bvec(fp,&cm);
    rr=fread_bvec(fp,&kk); bb=fread_bivec(fp,&kk);
    i=fread_bi(fp);
    if(i!=kk) error("wrong number of matrices");
    statmats=NEW_A(kk,double**);
    for(i=0;i<kk;i++) statmats[i]=fread_blmat(fp,&cm,bb+i);
    fclose(fp); FREE(cbuf);
  } else { /* ascii read from stdin */
    printf("\n# reading from stdin");
    orderv=read_ivec(&cm); obsvec=read_vec(&cm);
    rr=read_vec(&kk); bb=read_ivec(&kk);
    i=read_i();
    if(i!=kk) error("wrong number of matrices");
    statmats=NEW_A(kk,double**);
    for(i=0;i<kk;i++) {
      statmats[i]=read_lmat(&cm,bb+i); putdot();
    }
  }
  printf("\n# K:%d",kk);
  printf("\n# R:"); for(i=0;i<kk;i++) printf("%g ",rr[i]);
  printf("\n# B:"); for(i=0;i<kk;i++) printf("%d ",bb[i]);
  printf("\n# CM:%d",cm);

  read_alpha();

  if(sw_dobp) do_bptest();
  if(sw_doau) do_bootrep();

  write_pv(sw_dobp,0,0,sw_doau);
  if(sw_doau) {
    write_ci();
    if(sw_outrep) write_rep();
    if(sw_outcnt) write_cnt1();
  }

  return 0;
}

int do_cntmode()
{
  int j,i0;
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
  rrmat=fread_mat(fp,&cm,&kk); bbmat=fread_mat(fp,&cm,&kk);
  cntmat=fread_mat(fp,&cm,&kk);
  if(fname_cnt) {fclose(fp); FREE(cbuf);}

  printf("\n# CM:%d",cm);
  printf("\n# K:%d",kk);

  if(sw_dobp) {
    npvec=new_vec(cm); nsvec=new_vec(cm);
    for(j=0;j<cm;j++) {
      i0=find_i0(rrmat[j],kk);
      npvec[j]=cntmat[j][i0]/bbmat[j][i0];
      nsvec[j]=sqrt(npvec[j]*(1.0-npvec[j])/bbmat[j][i0]);
    }
  }
  if(sw_doau) do_bootcnt();
  write_pv(sw_dobp,0,0,sw_doau);

  return 0;
}

/*
  calculate naive bootstrap p-value

  cm-vector of p-value (se's)
  npvec (nsvec) : naive bootstrap

 */
int do_bptest()
{
  int j,i0;
  int nb;

  printf("\n# BP-TEST STARTS");
  i0=find_i0(rr,kk); nb=bb[i0];

  /* alloc memory */
  npvec=new_vec(cm); 

  /* calculate the naive p-value */
  for(j=0;j<cm;j++) npvec[j]=1.0-pvaldist(statmats[i0][j],nb,threshold);
  nsvec=getseval(npvec,cm,(double)nb,NULL);

  printf(" - DONE");
  return 0;
}

/* calculate approximate bayes posterior probability */
/*
  input: datvec,mm,assvec,asslen,cm
  output: bapvec
  work: buf1
 */
int do_batest()
{
  int i,j;
  double x,y;

  bapvec=new_vec(cm);

  if(npdm!=mm) {
    // bapcoef = 1
    // P = const*exp{bapcoef * logL } is computed
    // the following code is for avoiding overflow
    x=-HUGENUM; for(i=0;i<mm;i++) if(datvec[i]>x) x=datvec[i]; // find min
    y=0.0; for(i=0;i<mm;i++) {y+=buf1[i]=exp(bapcoef*(datvec[i]-x));}
    for(i=0;i<mm;i++) buf1[i]=buf1[i]/y; // normalization
  } else {
    // BIC correction 2011/01/26
    // P = const*exp{bapcoef*logL - 0.5*pdm*log(n) }
    printf("\n# ppn:%g",bapn);
    for(i=0;i<mm;i++) buf1[i]=bapcoef*datvec00[i] - 0.5*pdmvec[i]*log(bapn);
    x=-HUGENUM; for(i=0;i<mm;i++) if(buf1[i]>x) x=buf1[i];
    y=0.0; for(i=0;i<mm;i++) {y+=buf1[i]=exp(buf1[i]-x);}
    for(i=0;i<mm;i++) buf1[i]=buf1[i]/y; // normalization
  }

  for(i=0;i<cm;i++) {
    x=0.0; for(j=0;j<asslen[i];j++) x+=buf1[assvec[i][j]];
    bapvec[i]=x; 
  }
  return 0;
}

/*
  calculate p-values of conventional methods

  cm-vector of p-values
  (1) khpvec  : Kishino-Hasegawa test
  (2) khwpvec : Kishino-Hasegawa test with standardization weight
  (3) mcpvec  : Multiple Comparisons test
  (4) mcwpvec : MC-test with standardization weight

 */
int do_mctest()
{
  double x,*xp,*yp;
  int i,j,ib;
  int nb,i0;
  double **repmat; /* replicates of data */
  double **statmat; /* replicates of statistics */
  double *mnvec; /* mean of the replicates */
  double **wtmat; /* weight for ms-pvec */

  printf("\n# MC-TEST STARTS");
  i0=find_i0(rr,kk); nb=bb[i0];

  /* alloc memory */
  repmat=new_lmat(mm,nb);
  statmat=new_lmat(cm,nb);
  khpvec=new_vec(cm); khwpvec=new_vec(cm); 
  mcpvec=new_vec(cm); mcwpvec=new_vec(cm);
  mnvec=new_vec(mm); 
  wtmat=new_lmat(mm,mm);

  /* copy repmats[i0] to repmat */
  for(i=0;i<mm;i++) for(j=0;j<nb;j++) repmat[i][j]=repmats[i0][i][j];

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
  khsvec=getseval(khpvec,cm,(double)nb,NULL);
  
  /* calculate mc-pvalue */
  printf("\n# calculating mc-pvalue");  
  for(i=0;i<cm;i++) {
    mcpvec[i]=calcmcpval(datvec,repmat,mm,nb,
			 assvec[i],cassvec[i],asslen[i],NULL);
    putdot();
  }
  mcsvec=getseval(mcpvec,cm,(double)nb,NULL);

  /* calculate the weights */
  printf("\n# calculating the variances");
  for(i=0;i<mm;i++) {
     wtmat[i][i]=0.0; xp=repmat[i];
    for(j=0;j<i;j++) {
      yp=repmat[j]; x=0.0;
      for(ib=0;ib<nb;ib++) x+=(xp[ib]-yp[ib])*(xp[ib]-yp[ib]);
      if(x<varadd*VARADDWARN) 
	warning("small variance [%d,%d]: %g<%g*%g",i+1,j+1,
		x,varadd,VARADDWARN);
      wtmat[i][j]=wtmat[j][i]=1.0/sqrt((x+varadd)/(nb-1));
    }
    putdot();
  }

  /* calculate khw-pvalue */
  printf("\n# calculating weighted kh-pvalue");    
  for(i=0;i<cm;i++) {
    khwpvec[i]=calckhpval(datvec,repmat,mm,nb,
			  assvec[i],cassvec[i],asslen[i],wtmat);
    putdot();
  }
  khwsvec=getseval(khwpvec,cm,(double)nb,NULL);

  /* calculate mcw-pvalue */
  printf("\n# calculating weighted mc-pvalue");  
  for(i=0;i<cm;i++) {
    mcwpvec[i]=calcmcpval(datvec,repmat,mm,nb,
			  assvec[i],cassvec[i],asslen[i],wtmat);
    putdot();
  }
  mcwsvec=getseval(mcwpvec,cm,(double)nb,NULL);

  printf("\n# MC-TEST DONE");
  free_vec(mnvec); 
  free_lmat(repmat,mm); free_lmat(statmat,cm); free_lmat(wtmat,mm);
  return 0;
}

int do_bootrep()
{
  int i,j,k;
  double **statp,*beta;
  double t0,t1;

  printf("\n# AU-TEST STARTS");

  /* memory alloc for pv */
  pvvec=new_vec(cm); sevec=new_vec(cm);
  pv0vec=new_vec(cm); se0vec=new_vec(cm);
  betamat=new_mat(cm,3); 
  rssvec=new_vec(cm); dfvec=new_vec(cm); pfvec=new_vec(cm);
  thvec=new_vec(cm); 
  statp=NEW_A(kk,double*);

  /* memory alloc for ci */
  cimat=new_mat(cm,nalpha); csmat=new_mat(cm,nalpha);
  eimat=new_mat(cm,nalpha);
  ci0mat=new_mat(cm,nalpha); cs0mat=new_mat(cm,nalpha);
  ei0mat=new_mat(cm,nalpha);

  /* memory alloc for cnt */
  if(sw_outcnt) cntmat=new_mat(cm,kk);

  /* sorting the replicates */
  printf("\n# sorting the replicates");    
  for(i=0;i<kk;i++) {
    for(j=0;j<cm;j++) sort_vec(statmats[i][j],bb[i]);
    putdot();
  }

  /* calculate au-pvalues */
  printf("\n# calculating approximately unbiased p-values by ");
  switch(sw_fitmode) {
  case FITMODE_WLS:
    printf("WLS (very fast)"); break;
  case FITMODE_MLE:
    printf("MLE (fast)"); break;
  case FITMODE_CHI:
    printf("CHI (slow)"); break;
  }
  printf(" fitting"); fflush(STDOUT);
  t0=get_time();
  for(i=0;i<cm;i++) {
    mydprintf(1,"\n# rank=%d item=%d",i+1,orderv[i]+1);
    for(k=0;k<kk;k++) statp[k]=statmats[k][i];
    j=vcalpval(statp,rr,bb,kk,threshold,thvec+i,
	       pvvec+i,sevec+i,pv0vec+i,se0vec+i,
	       rssvec+i,dfvec+i,&beta,NULL,kappa);
    mydprintf(1,"\n# ret=%d",j);
    if(j) {
      pfvec[i]=0.0;
      warning("regression degenerated: df[%d]=%g",i+1,dfvec[i]);
      betamat[i][0]=betamat[i][1]=betamat[i][2]=0.0;
    } else {
      pfvec[i]=pochisq(rssvec[i],(int)(dfvec[i]));
      if(pfvec[i]<0.01)
	warning("theory does not fit well: pfit[%d]=%.4f",i+1,pfvec[i]);
      betamat[i][0]=beta[0]; betamat[i][1]=beta[1]; betamat[i][2]=beta[2];
    }
    /* make cntmat for output */
    if(sw_outcnt){
      for(k=0;k<kk;k++)
	cntmat[i][k]=cntdist(statmats[k][i],bb[k],threshold,
			     sw_smoothcnt?3:2);
    }
    putdot();
  }
  t1=get_time();
  printf("\n# time elapsed for AU test is t=%g sec",t1-t0);

  /* ci */
  printf("\n# ALPHA:");
  for(i=0;i<nalpha;i++) printf("%g ",alphavec[i]);
  printf("\n# calculating confidence intervals");
  if(sw_fitmode==FITMODE_CHI) sw_fitmode=FITMODE_MLE;
  for(j=0;j<cm;j++) {
    for(i=0;i<kk;i++) statp[i]=statmats[i][j];
    for(i=0;i<nalpha;i++) {
      mydprintf(1,"\n# invtpval item=%d alpha=%g ",j+1,alphavec[i]);
      k=invtpval(statp,rr,bb,kk,alphavec[i], 
		 &(cimat[j][i]),&(csmat[j][i]),&(eimat[j][i]),kappa);
      if(k) printf(" au item=%2d, alpha=%4.2f",
		   j+1,alphavec[i]);
      k=invtpval(statp,rr,bb,kk,alphavec[i], 
		 &(ci0mat[j][i]),&(cs0mat[j][i]),&(ei0mat[j][i]),-1.0);
      if(k) printf(" bp item=%2d, alpha=%4.2f",
		   j+1,alphavec[i]);
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
  double *beta;
  double t0,t1;

  printf("\n# AU-TEST STARTS");

  /* memory alloc for pv */
  pvvec=new_vec(cm); sevec=new_vec(cm);
  pv0vec=new_vec(cm); se0vec=new_vec(cm);
  betamat=new_mat(cm,3); 
  rssvec=new_vec(cm); dfvec=new_vec(cm); pfvec=new_vec(cm);
  thvec=new_vec(cm); 

  printf("\n# calculating approximately unbiased p-values by ");
  switch(sw_fitmode) {
  case FITMODE_WLS:
    printf("WLS (very fast)"); break;
  case FITMODE_MLE:
    printf("MLE (fast)"); break;
  case FITMODE_CHI:
    printf("CHI (slow)"); break;
  }
  printf(" fitting"); fflush(STDOUT);
  t0=get_time();
  for(i=0;i<cm;i++) {
    mydprintf(1,"\n# rank=%d item=%d",i+1,orderv[i]+1);
    j=rcalpval(cntmat[i],rrmat[i],bbmat[i],kk,
	       pvvec+i,sevec+i,pv0vec+i,se0vec+1,
	       rssvec+i,dfvec+i,&beta,NULL,kappa);
    thvec[i]=threshold; /* unused in cntmode */
    mydprintf(1,"\n# ret=%d",j);
    if(j) {
      pfvec[i]=0.0;
      warning("regression degenerated: df[%d]=%g",i+1,dfvec[i]);
      betamat[i][0]=betamat[i][1]=betamat[i][2]=0.0;
    } else {
      pfvec[i]=pochisq(rssvec[i],(int)(dfvec[i]));
      if(pfvec[i]<0.01)
	warning("theory does not fit well: pfit[%d]=%.4f",i+1,pfvec[i]);
      betamat[i][0]=beta[0]; betamat[i][1]=beta[1]; betamat[i][2]=beta[2];
    }
    putdot();
  }
  t1=get_time();
  printf("\n# time elapsed for AU test is t=%g sec",t1-t0);

  printf("\n# AU-TEST DONE");  
  return 0;
}


/* define bit sw */
#define BPPVSW 1
#define BAPVSW 2
#define MCPVSW 4
#define AUPVSW 8

/* number of entries */
#define BPPVNUM 1
#define BAPVNUM 1
#define MCPVNUM 4
#define AUPVNUM 2

/* number of aux entries */
#define BPAUXNUM 0
#define BAAUXNUM 0
#define MCAUXNUM 0
#define AUAUXNUM 7

int write_pv(int sw_bp, int sw_ba, int sw_mc, int sw_au)
{
  FILE *fp;
  char *cbuf;
  double **pvmat,**semat,**auxmat;
  int pvnum,auxnum,i,j;
  int outbit;

  if(fname_pv) {
    fp=openfp(fname_pv,fext_pv,"w",&cbuf);
    printf("\n# writing %s",cbuf);
  } else {
    fp=STDOUT;
    printf("\n# writing to stdout\n");
  }

  outbit=(sw_bp?BPPVSW:0)+(sw_ba?BAPVSW:0)
    +(sw_mc?MCPVSW:0)+(sw_au?AUPVSW:0);
  pvnum=(sw_bp?BPPVNUM:0)+(sw_ba?BAPVNUM:0)
    +(sw_mc?MCPVNUM:0)+(sw_au?AUPVNUM:0);
  auxnum=(sw_bp?BPAUXNUM:0)+(sw_ba?BAAUXNUM:0)+
    (sw_mc?MCAUXNUM:0)+(sw_au?AUAUXNUM:0);

  fprintf(fp,"\n# ID:\n%d\n",1); 

  fprintf(fp,"\n# ITEM:\n"); fwrite_ivec(fp,orderv,cm);
  fprintf(fp,"\n# STAT:\n"); fwrite_vec(fp,obsvec,cm);

  fprintf(fp,"\n# BIT:\n%d\n",outbit); 

  pvmat=new_mat(cm,pvnum); semat=new_mat(cm,pvnum);
  auxmat=new_mat(cm,auxnum);

  for(i=0;i<cm;i++) {
    j=0;
    if(sw_bp) {
      pvmat[i][j]=npvec[i]; semat[i][j]=nsvec[i]; j++;
    }
    if(sw_ba) {
      pvmat[i][j]=bapvec[i]; semat[i][j]=0.0; j++;
    }
    if(sw_mc) {
      pvmat[i][j]=khpvec[i]; semat[i][j]=khsvec[i]; j++;
      pvmat[i][j]=mcpvec[i]; semat[i][j]=mcsvec[i]; j++;
      pvmat[i][j]=khwpvec[i]; semat[i][j]=khwsvec[i]; j++;
      pvmat[i][j]=mcwpvec[i]; semat[i][j]=mcwsvec[i]; j++;
    }
    if(sw_au){
      pvmat[i][j]=pvvec[i]; semat[i][j]=sevec[i]; j++;
      pvmat[i][j]=pv0vec[i]; semat[i][j]=se0vec[i]; j++;
    }
  }

  for(i=0;i<cm;i++) {
    j=0;
    if(sw_au){
      auxmat[i][j]=pfvec[i]; j++;  /* pvalue of diagnostic */
      auxmat[i][j]=rssvec[i]; j++; /* rss */
      auxmat[i][j]=dfvec[i]; j++;  /* df */
      auxmat[i][j]=betamat[i][0]; j++; /* signed distance */
      auxmat[i][j]=betamat[i][1]; j++; /* curvature */
      auxmat[i][j]=thvec[i]; j++; /* threshold */
      auxmat[i][j]=betamat[i][2]; j++; /* dim */
    }
  }

  fprintf(fp,"\n# PV:\n"); fwrite_mat(fp,pvmat,cm,pvnum);  
  fprintf(fp,"\n# SE:\n"); fwrite_mat(fp,semat,cm,pvnum);  
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

  fprintf(fp,"\n# CI0:\n"); fwrite_mat(fp,ci0mat,cm,nalpha);
  fprintf(fp,"\n# SE0:\n"); fwrite_mat(fp,cs0mat,cm,nalpha);
  fprintf(fp,"\n# EI0:\n"); fwrite_mat(fp,ei0mat,cm,nalpha);

  if(fname_ci) {fclose(fp); FREE(cbuf);}
  return 0;
}

int write_rep()
{
  FILE *fp;
  char *cbuf;
  int i;

  if(fname_rep) { /* binary write to file */
    fp=openfp(fname_rep,fext_rep,"wb",&cbuf);
    printf("\n# writing %s",cbuf);
    fwrite_bivec(fp,orderv,cm); fwrite_bvec(fp,obsvec,cm);
    fwrite_bvec(fp,rr,kk); fwrite_bivec(fp,bb,kk);
    fwrite_bi(fp,kk);
    for(i=0;i<kk;i++) {
      fwrite_bmat(fp,statmats[i],cm,bb[i]); putdot();
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

int write_cnt1()
{
  int i,j;
  rrmat=new_mat(cm,kk); bbmat=new_mat(cm,kk);
  for(i=0;i<cm;i++) for(j=0;j<kk;j++) {
    rrmat[i][j]=rr[j]; bbmat[i][j]=(double)bb[j];
  }
  write_cnt2();
  free_mat(rrmat); free_mat(bbmat);
  return 0;
}

int write_cnt2()
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
  fprintf(fp,"\n# R:\n"); fwrite_mat(fp,rrmat,cm,kk);
  fprintf(fp,"\n# B:\n"); fwrite_mat(fp,bbmat,cm,kk);
  fprintf(fp,"\n# CNT:\n"); fwrite_mat(fp,cntmat,cm,kk);

  if(fname_cnt) {fclose(fp); FREE(cbuf);}
  return 0;
}


/* find i0 s.t. rr[i0] is close to 1.0 */
int find_i0(double *rr, int kk)
{
  double *buf;
  int i,i0;
  buf=new_vec(kk);
  for(i=0;i<kk;i++) buf[i]=fabs(rr[i]-1.0);
  i0=argmin_vec(buf,kk);
  free_vec(buf);

  return i0;
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
  if(t < vec[0]) return 0;
  else if(vec[bb-1] <= t) return bb;

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

void repminmaxsrs(double **repmat,  /* m x bb replicate mat */
		  double **statmat, /* cm x bb output */
		  int m,   /* number of itmes for input */
		  int bb,  /* number of replicates */
		  int cm,  /* number of itmes for output */
		  int **assvec, /* ass vectors of size cm */
		  int *asslen, /* length of each assvec[i] */
		  double *buf1, /* buf of size m */
		  double *buf2,  /* buf of size cm */
		  double *repmean, /* m x 1 mean of repmat */
		  double rs /* rescaling constnat */
		  )
{
  int i,j;

  for(i=0;i<bb;i++) {
    for(j=0;j<m;j++)
      buf1[j]=rs*(repmat[j][i]-repmean[j])+repmean[j];
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

  i0=j0=jj=0; /* supress warning */
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

double *getseval(double *pval, int m, double nb, double *seval)
{
  int i;
  if(seval==NULL) seval=new_vec(m);
  for(i=0;i<m;i++) seval[i]=sqrt(pval[i]*(1.0-pval[i])/nb);

  return seval;
}

/*
 *  MS-BOOT BY WLS
 */

int wlscalcpval(double *cnts, double *rr, double *bb, int kk,
		double *pv, double *se,    /* au */
		double *pv0, double *se0,  /* bp */
		double *rss, double *df, 
		double **betap, double ***vmatp, /* reference only */
		double rmin, double rmax, double kappa)
{
  int i;
  double x,s,**vmat;
  static double beta[3];
  static int kk0=0;
  static double *npv=NULL,*zval=NULL,*wt=NULL; /* kk-vector */
  static double **X=NULL; /* 2 x kk matrix */
  int m,m2; /* # of used items */

  if(kk>kk0){
    X=renew_mat(X,2,kk);
    npv=renew_vec(npv,kk);
    zval=renew_vec(zval,kk);
    wt=renew_vec(wt,kk);    
    kk0=kk;
  }

  m=m2=0;
  for(i=0;i<kk;i++) {
    s=sqrt(1.0/rr[i]);
    X[0][i]=1.0/s; X[1][i]=s;
    npv[i]=cnts[i]/bb[i];
    if((rr[i]<rmin)||(rr[i]>rmax)){
      zval[i]=0.0; wt[i]=0.0;
    } else if((npv[i]>=1.0)||(npv[i]<=0.0)){
      zval[i]=0.0; wt[i]=0.0; m2++;
    } else {
      zval[i]=-qnorm(npv[i]); 
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
    *pv0=*pv; *se0=*se;
    i=1; /* degenerate */
    mydprintf(1,"\n# error in wls");
  } else {
    lsfit(X,zval,wt,2,kk,beta,rss,&vmat);
    *pv=pnorm( -(beta[0]-kappa*beta[1]) );
    x=dnorm(beta[0]-kappa*beta[1]);
    *se=sqrt(x*x*(vmat[0][0]+kappa*kappa*vmat[1][1]
		  -kappa*vmat[0][1]-kappa*vmat[1][0]));
    *pv0=pnorm( -(beta[0]+beta[1]) );
    x=dnorm(beta[0]+beta[1]);
    *se0=sqrt(x*x*(vmat[0][0]+vmat[1][1]+vmat[0][1]+vmat[1][0]));
    i=0; /* non-degenerate */
  }
  beta[2]=0.0;
  if(betap) *betap=beta;  if(vmatp) *vmatp=vmat; 
  return i;
}


/* MS-BOOT by the MLE */
double log3(double x)
{
  double y,z1,z2,z3,z4,z5;
  if(fabs(x)>1.0e-3) {
    y=-log(1.0-x);
  } else {
    z1=x; z2=z1*x; z3=z2*x; z4=z3*x; z5=z4*x;
    y=((((z5/5.0)+z4/4.0)+z3/3.0)+z2/2.0)+z1;
  }
  return y;
}
int mleloopmax=30;
double mleeps=1e-10;
int mlecoef(double *cnts, double *rr, double *bb, int kk,
	    double *coef0, /* set initinal value (size=2) */
	    double **vmat, /* variance (2x2) */
	    double *lrt, double *df, /* LRT statistic */
	    double rmin, double rmax)
{
  int i,m,loop;
  double coef[2], update[2];
  double d1f, d2f, d11f, d12f, d22f; /* derivatives */
  double v11, v12, v22; /* inverse of -d??f */
  double a,e;
  static int kk0=0;
  static double *s=0,*r=0,*c=0,*b=0,*z=0,*p=0,*d=0,*g=0,*h=0;

  if(kk>kk0) {
    s=renew_vec(s,kk); r=renew_vec(r,kk); c=renew_vec(c,kk);
    b=renew_vec(b,kk); z=renew_vec(z,kk); p=renew_vec(p,kk); 
    d=renew_vec(d,kk); g=renew_vec(g,kk); h=renew_vec(h,kk); 
    kk0=kk;
  }

  m=0;
  for(i=0;i<kk;i++)
    if(rr[i]>=rmin && rr[i]<=rmax) {
      r[m]=rr[i]; s[m]=sqrt(rr[i]); c[m]=cnts[i]; b[m]=bb[i];
      m++;
    }
  if(m<2) return 1;

  coef[0]=coef0[0]; /* signed distance */
  coef[1]=coef0[1]; /* curvature */
  mydprintf(3,"\n### mlecoef: deg=%d, coef0=(%g,%g)",m-2,coef[0],coef[1]);

  for(loop=0;loop<mleloopmax;loop++) {
    d1f=d2f=d11f=d12f=d22f=0.0;
    for(i=0;i<m;i++) {
      z[i]=coef[0]*s[i]+coef[1]/s[i];
      p[i]=pnorm(-z[i]);
      d[i]=dnorm(z[i]);
      if(p[i]>0.0 && p[i]<1.0) {
	g[i]=d[i]*( d[i]*(-c[i]+2.0*c[i]*p[i]-b[i]*p[i]*p[i])/
		    (p[i]*p[i]*(1.0-p[i])*(1.0-p[i]))
		    + z[i]*(c[i]-b[i]*p[i])/(p[i]*(1.0-p[i])) );
	h[i]=d[i]*(c[i]-b[i]*p[i])/(p[i]*(1.0-p[i]));
      } else { g[i]=h[i]=0.0; }
      d1f+= -h[i]*s[i]; d2f+= -h[i]/s[i];
      d11f+= g[i]*r[i]; d12f+= g[i]; d22f+= g[i]/r[i];
      mydprintf(3,"\n### mlecoef: %d %g %g %g %g %g",i,z[i],p[i],d[i],g[i],h[i]);
    }

    a=d11f*d22f-d12f*d12f;
    if(a==0.0) {
      mydprintf(1,"\n# singular matrix in mle");
      return 2;
    }
    v11=-d22f/a; v12=d12f/a; v22=-d11f/a;

    /* Newton-Raphson update */
    update[0]=v11*d1f+v12*d2f; update[1]=v12*d1f+v22*d2f;
    coef[0]+=update[0]; coef[1]+=update[1];

    /* check convergence */
    e=-d11f*update[0]*update[0]-2.0*d12f*update[0]*update[1]
      -d22f*update[1]*update[1];

    mydprintf(2,"\n## mlecoef: %d %g %g %g %g %g",
	    loop,coef[0],coef[1],update[0],update[1],e);
    if(e<mleeps) break;
  }

  /* calc log-likelihood */
  *lrt=0.0; *df=0.0;
  for(i=0;i<m;i++) {
    if(p[i]>0.0 && p[i]<1.0) {
      *df+=1.0;
      if(c[i]>0.0) a=c[i]*log(c[i]/b[i]/p[i]); else a=0.0;
      if(c[i]<b[i]) a+=(b[i]-c[i])*(log3(p[i])-log3(c[i]/b[i]));
      *lrt += a;
    }
  }
  *lrt *= 2.0; *df -= 2.0;

  /* write back the results */
  coef0[0]=coef[0]; coef0[1]=coef[1];
  vmat[0][0]=v11;vmat[0][1]=vmat[1][0]=v12;vmat[1][1]=v22; 
  if(loop==mleloopmax || *df< -0.01) i=1; else i=0;
  return i;
}

int mlecalcpval(double *cnts, double *rr, double *bb, int kk,
		double *pv, double *se,    /* au */
		double *pv0, double *se0,  /* bp */
		double *rss, double *df, 
		double **betap, double ***vmatp, /* reference only */
		double rmin, double rmax, double kappa)
{
  int i;
  double x,*coef0;
  static double coef[3], **vmat=NULL;
  if(!vmat) vmat=new_mat(2,2);

  i=wlscalcpval(cnts,rr,bb,kk,pv,se,pv0,se0,
		rss,df,&coef0,vmatp,rmin,rmax,kappa);
  if(betap) *betap=coef0;
  if(i) return i; /* error */
  
  coef[0]=coef0[0], coef[1]=coef0[1];
  i=mlecoef(cnts,rr,bb,kk,coef,vmat,rss,df,rmin,rmax);
  if(i) {mydprintf(0,"\n# error in mle, use wls instead"); return 0;}

  *pv=pnorm(-(coef[0]-kappa*coef[1]) );
  x=dnorm(coef[0]-kappa*coef[1]);
  *se=sqrt(x*x*(vmat[0][0]+kappa*kappa*vmat[1][1]
		-kappa*vmat[0][1]-kappa*vmat[1][0]));
  *pv0=pnorm(-(coef[0]+coef[1]));
  x=dnorm(coef[0]+coef[1]);
  *se0=sqrt(x*x*(vmat[0][0]+vmat[1][1]+vmat[0][1]+vmat[1][0]));
  if(betap) *betap=coef;  if(vmatp) *vmatp=vmat;
  coef[2]=0.0;

  return 0;
}

/*
 *  MS-BOOT BY SPHERICAL MODEL
 */


double chilrt(double *cc, double *rr, double *bb, int kk,
	      double sid, double a, double m)
/* sid = signed distance. a is the radius of the
 spherical region. m is the dimension of the space */
{
  double aa,cv,x,xx,lrt,lik,bp,r;
  int k,aneg;

  aneg=a<0; x=a+sid;
  xx=x*x; aa=a*a;
  lrt=0.0;
  for(k=0;k<kk;k++) {
    if(xx*rr[k]<1e6) {
      if(!aneg) bp=pchisqnc(aa*rr[k],m,xx*rr[k]);
      else bp=tchisqnc(aa*rr[k],m,xx*rr[k]);
    } else {
      r=sqrt(rr[k]); 
      cv=(m-1.0)*(0.5/x+sid*0.25/xx);
      bp=pnorm(-sid*r-cv/r);
    }
    if(bp>0.0 && bp<1.0) {
      if(cc[k]>0.0) lik=cc[k]*log(cc[k]/bb[k]/bp); else lik=0.0;
      if(cc[k]<bb[k]) lik+=(bb[k]-cc[k])*(log3(bp)-log3(cc[k]/bb[k]));
      lrt += lik;
    }
  }
  lrt *= 2.0;
  return lrt;
}

double chipval(double sid, double a, double m)
/* sid = signed distance. a is the radius of the
 spherical region. m is the dimension of the space */
{
  double aa,cv,x,xx,pval;
  int aneg;

  aneg=a<0; x=a+sid;
  xx=x*x; aa=a*a;
  if(aa<1e6) {
    if(!aneg) pval=tchisqnc(xx,m,aa);
    else pval=pchisqnc(xx,m,aa);
  } else {
    cv=(m-1.0)*(0.5/a-sid*0.25/aa);
    pval=pnorm(-sid+cv);
  }
  return pval;
}

double chibp(double sid, double a, double m)
/* sid = signed distance. a is the radius of the
 spherical region. m is the dimension of the space */
{
  double aa,cv,x,xx,bp;
  int aneg;

  aneg=a<0; x=a+sid;
  xx=x*x; aa=a*a;
  if(xx<1e6) {
    if(!aneg) bp=pchisqnc(aa,m,xx);
    else bp=tchisqnc(aa,m,xx);
  } else {
    cv=(m-1.0)*(0.5/x+sid*0.25/xx);
    bp=pnorm(-sid-cv);
  }
  return bp;
}

#define CHIFUNCP chifuncp
int chifuncp=2;
double chixh1z[]={0.001,0.001,0.001};
double chixh2z[]={0.0002,0.0002,0.0002};
double chixh1[3], chixh2[3];
double chidim=0.0; /* dimension of the space */
int chikk; double *chicc,*chirr,*chibb;
#define PAR2DIM(P) (1.0+exp(P))
#define DIM2PAR(D) (log((D)-1.0))
#define PARM0 parm[0]
#define PARM2 (chifuncp==2)?chidim:PAR2DIM(parm[2])
#define PARM1 0.5*((PARM2)-1.0)/parm[1]

double chimodel(double *parm)
{
  double x;
  x=chilrt(chicc,chirr,chibb,chikk,PARM0,PARM1,PARM2);
  if(x==0.0) x=1e30;
  mydprintf(4,"\n#### P: %g %g = %g",parm[0],parm[1],x);

  return x;
}
double chiobj(double *parm)
{
  return chipval(PARM0,PARM1,PARM2);
}
double chiobj0(double *parm)
{
  return chibp(PARM0,PARM1,PARM2);
}
void chiseth(double *parm)
{
  int i; double x;
  for(i=0;i<CHIFUNCP;i++) {
    x=fabs(parm[i]); if(x<10.0) x=10.0;
    chixh1[i]=x*chixh1z[i];
    chixh2[i]=x*chixh2z[i];
  }
}
void chidmodel(double *parm, double *diff)
{
  chiseth(parm);
  dfmpridr(parm,diff,chixh1,CHIFUNCP,chimodel);
  mydprintf(3,"\n### D: %g %g",diff[0],diff[1]);
}
void chiddmodel(double *parm, double **diff)
{
  chiseth(parm);
  ddfsimple(parm,diff,chixh2,CHIFUNCP,chimodel);
}
void chidobj(double *parm, double *diff)
{
  chiseth(parm);
  dfmpridr(parm,diff,chixh1,CHIFUNCP,chiobj);
}
void chidobj0(double *parm, double *diff)
{
  chiseth(parm);
  dfmpridr(parm,diff,chixh1,CHIFUNCP,chiobj0);
}

#define NEWTON {chidim=PAR2DIM(parm[2]);\
               dfnmin(parm,2,chieps,chiloopmax,chibackmax,\
               &loop,&y,&vmat,chimodel, chidmodel, chiddmodel);}
int chicoef(double *cc, double *rr, double *bb, int kk,
	    double *parm0, /* set initinal value (size=3) */
	    double ***vmatp, /* variance ptr (3x3) */
	    double *lrt, double *df, /* LRT statistic */
	    double rmin, double rmax)
{
  int i,m,loop,is;
  double parm[3],x,y,**vmat;
  double parm2max, parm2min, parm2grid;
  double ***vmata, *ya, **parma;
  int im,ik;
  int im1,im2,im3,im4;
  double x1,x2,x3;
  static int kk0=0;
  static double *cc0=0,*rr0=0,*bb0=0;

  if(kk>kk0) {
    cc0=renew_vec(cc0,kk); rr0=renew_vec(rr0,kk);
    bb0=renew_vec(bb0,kk); kk0=kk;
  }

  m=0; x=0.0;
  for(i=0;i<kk;i++)
    if(rr[i]>=rmin && rr[i]<=rmax) {
      rr0[m]=rr[i]; cc0[m]=cc[i]; bb0[m]=bb[i];
      x+=bb[i]; m++;
    }

  chifuncp=2;
  if(m<CHIFUNCP) return 1;
  chicc=cc; chirr=rr0; chibb=bb0; chikk=m;
  *df=m-CHIFUNCP;

  for(i=0;i<3;i++) parm[i]=parm0[i];  

  NEWTON;
  mydprintf(1,"\n# opt grid=%d loop=%d sid=%g cv=%g dim=%g lrt=%g",
	  0,loop,parm[0],parm[1],PAR2DIM(parm[2]),y);

  if(sw_chidimopt && (*df)>=1 && PAR2DIM(parm0[2])>chidimmin+0.1) {
    (*df)--;
    /* save the parameters */
    ik=chinumgrid+1;
    vmata=NEW_A(ik,double**); ya=new_vec(ik); parma=new_mat(ik,3);
    ya[0]=y; for(i=0;i<3;i++) parma[0][i]=parm[i]; vmata[0]=vmat;
    im=0;
    /* find the optimal dimension by grid search */
    parm2max=parm0[2]; parm2min=DIM2PAR(chidimmin);
    parm2grid=(parm2max-parm2min)/(ik-1);
    for(is=1;is<ik;is++) {
      parm[2]=parm2max-is*parm2grid;
      NEWTON;
      ya[is]=y; for(i=0;i<3;i++) parma[is][i]=parm[i]; vmata[is]=vmat;
      if(y<ya[im]) im=is;
      mydprintf(1,"\n# opt grid=%d loop=%d sid=%g cv=%g dim=%g lrt=%g",
	      is,loop,parm[0],parm[1],PAR2DIM(parm[2]),y);
    }
    if(sw_chidimgs && im>0 && im<ik-1) {
#define GOLDRATIO 0.38196601125010515
#define GOLDLOOPMAX 10
      im1=im+1; im2=im; im3=im-1;
      /* golden section search */
      for(is=0;is<GOLDLOOPMAX;is++) {
	x1=parma[im1][2]; x2=parma[im2][2]; x3=parma[im3][2];
	if(x3-x2>=x2-x1) {
	  x=x2+GOLDRATIO*(x3-x2);
	} else {
	  x=x2-GOLDRATIO*(x2-x1);
	}
	for(i=0;i<3;i++) parm[i]=parma[im2][i];	parm[2]=x;
	NEWTON;
	mydprintf(1,"\n# opt gold=%d loop=%d sid=%g cv=%g dim=%g lrt=%g",
		is,loop,parm[0],parm[1],PAR2DIM(parm[2]),y);
	if(y<ya[im2]){
	  /* new point becomes the bracketted point */
	  if(x>x2) {
	    /* x1 will be deleted: x1, (x2, x, x3) */
	    i=im1; im1=im2; im2=i;
	  } else {
	    /* x3 will be deleted: (x1, x, x2), x3 */
	    i=im3; im3=im2; im2=i;
	  }
	  im4=im2;
	} else {
	  /* new point becomes a bracketing point */
	  if(x>x2) {
	    /* x3 will be deleted: (x1, x2, x), x3 */
	    im4=im3;
	  } else {
	    /* x1 will be deleted: x1, (x, x2, x3) */
	    im4=im1;
	  }
	}
	free_mat(vmata[im4]); vmata[im4]=vmat; ya[im4]=y;
	for(i=0;i<3;i++) parma[im4][i]=parm[i];
      }
      im=im2;
#undef GOLDRATIO      
#undef GOLDLOOPMAX
    } 
    y=ya[im]; for(i=0;i<3;i++) parm[i]=parma[im][i]; vmat=vmata[im];
    mydprintf(1,"\n# opt im=%d sid=%g cv=%g dim=%g lrt=%g",
	    im,parm[0],parm[1],PAR2DIM(parm[2]),y);
    for(i=0;i<ik;i++) {
      if(vmata[i]!=vmat) free_mat(vmata[i]);
    }
    free(vmata); free_mat(parma); free_vec(ya);
  }

  *lrt=y; *vmatp=vmat; for(i=0;i<3;i++) parm0[i]=parm[i]; 
  return 0;
}
#undef NEWTON

int chicalcpval(double *cnts, double *rr, double *bb, int kk,
		double *pv, double *se,    /* au */
		double *pv0, double *se0,  /* bp */
		double *rss, double *df, 
		double **betap, double ***vmatp, /* reference only */
		double rmin, double rmax, double kappa)
{
  int i,j;
  double x,*beta0,**vmat,rss0,df0;
  static double parm[3], diff[3], beta[3];

  i=wlscalcpval(cnts,rr,bb,kk,pv,se,pv0,se0,
		rss,df,&beta0,vmatp,rmin,rmax,kappa);
  if(betap) *betap=beta0;
  if(i) return i; /* error */

  if(chidimarg!=0.0) chidim=chidimarg;
  else chidim=(mm-1.0)*mf;
  if(chidim<chidimmin) chidim=chidimmin;
  parm[0]=beta0[0]; /* signed distance */
  parm[1]=beta0[1]; /* curvature */
  parm[2]=DIM2PAR(chidim);

  i=chicoef(cnts,rr,bb,kk,parm,&vmat,&rss0,&df0,rmin,rmax);
  if(vmatp) *vmatp=vmat;
  if(i) {mydprintf(0,"\n# error in chi, use wls instead"); return 0;}
  *rss=rss0; *df=df0;

  *pv=chiobj(parm); chidobj(parm,diff);
  x=0.0; 
  for(i=0;i<CHIFUNCP;i++) for(j=0;j<CHIFUNCP;j++)
    x+=diff[i]*diff[j]*vmat[i][j];
  if(x<0.0) { warning("negative variance of au in chicalcpval"); x=0.0; }
  *se=sqrt(x);

  chidim=PAR2DIM(parm[2]);
  *pv0=chiobj0(parm); chidobj0(parm,diff);
  x=0.0; 
  for(i=0;i<CHIFUNCP;i++) for(j=0;j<CHIFUNCP;j++)
    x+=diff[i]*diff[j]*vmat[i][j];
  if(x<0.0) { warning("negative variance of np in chicalcpval"); x=0.0; }
  *se0=sqrt(x);

  /* for compatibility */
  beta[0]=parm[0]; beta[1]=parm[1]; beta[2]=PAR2DIM(parm[2]);
  mydprintf(1,"\n# pv=%g (%g) pv0=%g (%g) sid=%g cv=%g a=%g m=%g",
	  *pv,*se,*pv0,*se0,beta[0],beta[1],PARM1,beta[2]);
  if(betap) *betap=beta;
  if(!vmatp) free_mat(vmat);
  return 0;
}

/*

  MS-BOOT INTERFACE

 */

/* switch the wls and the mle for calculating au test */
int calcpval(double *cnts, double *rr, double *bb, int kk,
	     double *pv, double *se,    /* au */
	     double *pv0, double *se0,  /* bp */
	     double *rss, double *df, 
	     double **betap, double ***vmatp, /* reference only */
	     double rmin, double rmax, double kappa)
{
  switch(sw_fitmode) {
  case FITMODE_WLS:
    return wlscalcpval(cnts,rr,bb,kk,pv,se,pv0,se0,rss,df, 
		       betap,vmatp,rmin,rmax,kappa);
  case FITMODE_MLE:
    return mlecalcpval(cnts,rr,bb,kk,pv,se,pv0,se0,rss,df, 
		       betap,vmatp,rmin,rmax,kappa);
  case FITMODE_CHI:
    return chicalcpval(cnts,rr,bb,kk,pv,se,pv0,se0,rss,df, 
		       betap,vmatp,rmin,rmax,kappa);
  default:
    return 2;
  }
}


double pfmin=0.0; /* minimum p-value of diagnostics */
double rrmin=0.0; /* minimum rr */
double rrmax=100.0; /* maximum rr */
int rcalpval(double *cnts, double *rr, double *bb, int kk,
	     double *pv, double *se, 
	     double *pv0, double *se0, 
	     double *rss, double *df, 
	     double **betap, double ***vmatp, /* reference only */
	     double kappa)
{
  static int kk0=0;
  static double *rr0=NULL;
  int i,i1,i2,j=1;
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
  /* find the minimum rr[i] which is not less than 1.0 */
  for(;i<kk;i++) if(rr0[i]>=1.0) break;
  i2=i;
  /* call calcpval until the theory fits well */
  for(i=i1;i<i2;i++) {
    j=calcpval(cnts,rr,bb,kk,pv,se,pv0,se0,rss,df,betap,vmatp,
	       rr0[i],rrmax,kappa);
    if(j) break;
    pfit=pochisq(*rss,(int)(*df));
    if(pfit>=pfmin) break;
  }

  return j;
}

int vcloopmax=30;
double vcscale=0.5;
double vceps1=0.01;
double vceps2=0.1;
int dfmin=0; /* kk >= 2 */
int vcalpval(double **statps, double *rr, int *bb, int kk,
	     double threshold, double *thp,
	     double *pvp, double *sep, double *pv0p, double *se0p, 
	     double *rssp, double *dfp, 
	     double **betap, double ***vmatp, /* reference only */
	     double kappa)
{
  int i,j,idf,idf0;
  double x,xn,pv,se,rss,df,z0,z,ze0,ze,th,pv0,se0;
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

  *pvp=*sep=*pv0p=*se0p=*rssp=*dfp=*thp=0.0;
  if(betap) *betap=NULL; if(vmatp) *vmatp=NULL;
  z=ze=0.0; th=threshold; idf0=-2;
  z0=ze0=0.0;
  /* iteration */
  for(j=0;j<vcloopmax;j++) {
    x=xn;
    /* get cnts */
    for(i=0;i<kk;i++) {
      cntvec[i]=cntdist(statps[i],bb[i],x,3);
      buf[i]=(double)bb[i];
    }
    /* get pvalue */
    i=calcpval(cntvec,rr,buf,kk,&pv,&se,&pv0,&se0,
	       &rss,&df,&beta,&vmat,rrmin,rrmax,kappa);
    if(!i) { z=-pv; ze=se; }
    idf=(int)df;
    mydprintf(2,
"\n## vcalpval: %3d %6.3f %6.3f %9.6f %9.6f %3.0f %6.3f %6.3f %6.3f %6.3f %g",
	    j,th,x,pv,se,df,beta[0],beta[1],z,ze,rss);

    if(idf < dfmin && idf0 < dfmin) return 1; /* degenerated */
    if((idf < dfmin) || 
       (idf0 >= dfmin && (z-z0)*(x-*thp) > 0.0 && fabs(z-z0)>vceps2*ze0) ) {
      mydprintf(1,"\n# non-monotone");
      th=x; xn=vcscale*x+(1.0-vcscale)*(*thp);
      continue;
    }
    if(idf0 >= dfmin && (fabs(z-z0)<vceps1*ze0)) {
      if(fabs(th-threshold)<1e-10) xn=th; else th=x;
    } else xn=vcscale*th+(1.0-vcscale)*x;
    *pvp=pv;*sep=se;*rssp=rss;*dfp=df;*thp=x; 
    *pv0p=pv0;*se0p=se0;
    if(betap) *betap=beta;  if(vmatp) *vmatp=vmat;

    z0=z;ze0=ze; idf0=idf;
    if(fabs(x-th)<1e-10) return 0;
  }
  
  mydprintf(1,"\n# no-convergence");
  return 2;
}


/*
  statps[i] must be sorted
 */
#define CILOOPMAX 30
double cieps=0.1;
int invtpval(double **statps, double *rr, int *bb, int kk,
	     double alpha, 
	     double *cival, double *seval, double *eival,
	     double kappa)
{
  int i,j;
  double x1,x2,x,x0,e,pv,pv0,se,se0,rss,df,*beta,pv00,se00;
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
    mydprintf(2,"\n## invtpval: j=%d x=%g ",j,x);
    /* get cnts */
    for(i=0;i<kk;i++) {
      cntvec[i]=cntdist(statps[i],bb[i],x,3);
      buf[i]=(double)bb[i];
    }
    /* get pvalue */
    se0=se; pv0=pv;
    debugmode--;
    i=calcpval(cntvec,rr,buf,kk,&pv,&se,&pv00,&se00,
	       &rss,&df,&beta,NULL,rrmin,rrmax,kappa);
    debugmode++;
    e=pv-alpha;
    mydprintf(2,"%g %g %g %g %g %g %g %g ",
	    pv,se,rss,df,beta[0],beta[1],e,e/se);
    /* check */
    if((pv0-alpha) * (pv-pv0) > 0.0)
      mydprintf(1,"non-monotone in ci (%g %g)->(%g %g)",x0,pv0,x,pv);
    if(i) {
      warning("degenerecy in pval df=%g",df);
      *cival=x0; *seval=0.0; *eival=0.0;
      return 1;
    }
    if((j>0) && (fabs(e) < se*cieps)) break;
    if(e>0.0) x2=x; else x1=x; /* update */
  }
  *cival=x;
  *seval=se*(x-x0)/(pv-pv0);
  *eival=e/se;

  if(j==CILOOPMAX) {
    warning("no convergence in ci: ei=%6.3f",e/se);
    return 1;
  }

  return 0;
}
