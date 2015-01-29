/* randrep.c

  Time-stamp: <2011-01-25 16:57:26 shimo>

  shimo@ism.ac.jp 
  Hidetoshi Shimodaira

  typical usage:
  # foo.vt -> foo.rep
  randrep foo

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rand.h"
#include "misc.h"

static const char rcsid[] = "$Id: randrep.c,v 1.6 2011/05/12 07:22:30 shimo Exp $";

void putdot() {putchar('.'); fflush(STDOUT);}
void byebye() {error("error in command line");}


char *fnamenum(char *basename, int num, int max)
{
  int len1,len2,len;
  char *buf;
  if(basename==NULL) return NULL;
  max--; if(max<1) max=1;
  len1=(int)floor(log10((double)max)+2.000001); len2=strlen(basename);
  len=len1+len2+10;  buf=NEW_A(len,char);
  sprintf(buf,"%s%0*d",basename,len1,num);
  return buf;
}

unsigned long seed=0;
char *fname_pa = NULL;
char *fname_rep = NULL;
char *fname_rmt = NULL;
char *fname_vt = NULL;
char *fext_pa = ".pa";
char *fext_rep = ".rep";
char *fext_rmt = ".rmt";
char *fext_vt = ".vt";
char *fname_in=NULL, *fname_out=NULL;


/* for the scales */
int kk=0; /* no. of scales */
double *rr=NULL; /* kk-vector of relative sample sizes */
int *bb=NULL; /* kk-vector of no.'s of replicates */

/* scale parameters */
int kk0=10;
double rr0[]={0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4};
int bb0[]={10000,10000,10000,10000,10000,10000,10000,10000,10000,10000};

double *rr1=0; /* for save */

/* for rep file */
double ***statmats=NULL; /* (kk,cm,bb[i])-array of test statistics */
double *obsvec=NULL; /* cm-vector of observed statistics */
int *orderv=NULL; 
int cm=0; /* number of itmes */

/* for rmt file */
double ***repmats=NULL; /* (kk,mm,bb[i])-array of replicates */
double *datvec=NULL; /* mm-vector of data */
int mm=0; /* number of items */

/* for repeat */
int sw_repeat=0; int num_repeat=1;
double flucscale=0.0;
int mode_rmt=0;
int sw_nonpara=0; 

enum repmodel {REP_CHI, REP_EXP};
enum rmtmodel {RMT_NORM};
int repmode=REP_CHI;
int rmtmode=RMT_NORM;

/* sub routines */
int do_repchi();
int do_repexp();
int do_rmtnorm();
int write_rep(char *name);
int write_rmt(char *name);

int main(int argc, char** argv)
{
  /* working variables */
  int i,j;
  char *cbuf;
  FILE *fp;

  printf("# %s",rcsid);

  /* args */
  for(i=j=1;i<argc;i++) {
    if(argv[i][0] != '-') {
      switch(j) {
      case 1:
	fname_in=argv[i];
	fname_out=rmvaxt(argv[i]);
	break;
      case 2: fname_out=argv[i]; break;
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
    } else if(streq(argv[i],"-d")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%d",&debugmode) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"-r")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%d",&num_repeat) != 1)
	byebye();
      sw_repeat=1;
      i+=1;
    } else if(streq(argv[i],"-f")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%lf",&flucscale) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"-m")) {
      mode_rmt=1;
    } else if(streq(argv[i],"--nonpara")) {
      sw_nonpara=1;
    } else if(streq(argv[i],"--model")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%d",&repmode) != 1)
	byebye();
      rmtmode=repmode;
      i+=1;
    } else byebye();
  }

  /* random seed */
  smrand(seed);

  /* read parameter file */
  if(fname_pa) {
    fp=openfp(fname_pa,fext_pa,"r",&cbuf);
    printf("\n# reading %s",cbuf);
    kk=0;
    rr=fread_vec(fp,&kk);
    bb=fread_ivec(fp,&kk);
    fclose(fp); FREE(cbuf);
  } else {
    kk=kk0; rr=rr0; bb=bb0; 
  }
  rr1=new_vec(kk);
  for(i=0;i<kk;i++) rr1[i]=rr[i]; /* save the original */

  printf("\n# seed:%lu",seed);
  printf("\n# K:%d",kk);
  printf("\n# R:"); for(i=0;i<kk;i++) printf("%g ",rr[i]);
  printf("\n# B:"); for(i=0;i<kk;i++) printf("%d ",bb[i]);

  /* random number generation */
  if(mode_rmt) {
    fname_vt=fname_in; fname_rmt=fname_out;
    do_rmtnorm();
  } else {
    fname_vt=fname_in; fname_rep=fname_out;
    switch(repmode) {
    case REP_CHI: do_repchi(); break;
    case REP_EXP: do_repexp(); break;
    default: byebye();
    }
  }

  printf("\n# exit normally\n");
  exit(0);
}

int do_repchi()
{
  int i,j,k,ib;
  FILE *fp;
  char *cbuf;
  double *degvec,*xp,*xbuf,*ncvec;
  double nc,x,y,r2,t0,t1,sg;
  int deg,irep;

  if(fname_vt){ 
    fp=openfp(fname_vt,fext_vt,"r",&cbuf);
    printf("\n# reading %s",cbuf);
  } else {
    fp=STDIN;
    printf("\n# reading from stdin\n");
  }
  cm=0; ncvec=fread_vec(fp,&cm); degvec=fread_vec(fp,&cm);
  if(fname_vt) {fclose(fp); FREE(cbuf);}

  printf("\n# CM:%d",cm);
  printf("\n# OBS:"); for(i=0;i<cm;i++) printf("%g ",ncvec[i]);
  printf("\n# DEG:"); for(i=0;i<cm;i++) printf("%g ",degvec[i]);
  deg=0; for(i=0;i<cm;i++) {
    if((degvec[i] < 0.99) || fabs(degvec[i]-floor(degvec[i]+0.5))>0.01)
      error("degrees of freedom must be positive integer");
    if(degvec[i]>deg) deg=(int)degvec[i];
  }
  xbuf=new_vec(deg+1); obsvec=new_vec(cm);
  statmats=NEW_A(kk,double**);
  for(i=0;i<kk;i++) statmats[i]=new_mat(cm,bb[i]);
  orderv=new_ivec(cm); for(i=0;i<cm;i++) orderv[i]=i;


  for(i=j=0;i<kk;i++) j+=bb[i];
  printf("\n# start generating %d repeat%s of\
 total %d replicates for %d items",
	 num_repeat,(num_repeat>1)?"s":"",j,cm); fflush(STDOUT);
  t0=get_time();
  for(irep=0;irep<num_repeat;irep++) {
    for(i=0;i<cm;i++) {
      deg=(int)degvec[i]; nc=ncvec[i]; sg=(nc>=0)?1.0:-1.0;
      for(j=0;j<deg;j++) xbuf[j]=flucscale*rnorm();
      xbuf[0]+=nc;
      x=0.0; for(j=0;j<deg;j++) x+=xbuf[j]*xbuf[j];
      obsvec[i]=sg*sqrt(x);
      for(k=0;k<kk;k++) {
	r2=1.0/sqrt(rr[k]); xp=statmats[k][i];
	for(ib=0;ib<bb[k];ib++) {
	  x=0.0; for(j=0;j<deg;j++) {y=xbuf[j]+r2*rnorm(); x+=y*y;}
	  xp[ib]=sg*sqrt(x);
	}
      }
      putdot();
    }
    if(sw_repeat) {
      cbuf=fnamenum(fname_rep,irep,num_repeat);
      write_rep(cbuf); FREE(cbuf);
    } else write_rep(fname_rep);    
  }
  t1=get_time();
  printf("\n# time elapsed for replicate generation t=%g sec",t1-t0);
  fflush(STDOUT);

  return 0;
}

int do_repexp()
{
  int i,j,k,ib;
  FILE *fp;
  char *cbuf;
  double *degvec,*xp,*xbuf,*ncvec;
  double nc,x,r2,t0,t1,sg,y;
  int deg,irep,bn;

  if(fname_vt){ 
    fp=openfp(fname_vt,fext_vt,"r",&cbuf);
    printf("\n# reading %s",cbuf);
  } else {
    fp=STDIN;
    printf("\n# reading from stdin\n");
  }
  cm=0; ncvec=fread_vec(fp,&cm); degvec=fread_vec(fp,&cm);
  if(fname_vt) {fclose(fp); FREE(cbuf);}

  printf("\n# CM:%d",cm);
  printf("\n# MEAN:"); for(i=0;i<cm;i++) printf("%g ",ncvec[i]);
  printf("\n# SAMPLE SIZE:"); for(i=0;i<cm;i++) printf("%g ",degvec[i]);
  deg=0; for(i=0;i<cm;i++) {
    if((degvec[i] < 0.99) || fabs(degvec[i]-floor(degvec[i]+0.5))>0.01)
      error("sample size must be positive integer");
    if(degvec[i]>deg) deg=(int)degvec[i];
  }
  xbuf=new_vec(deg+1); obsvec=new_vec(cm);
  statmats=NEW_A(kk,double**);
  for(i=0;i<kk;i++) statmats[i]=new_mat(cm,bb[i]);
  orderv=new_ivec(cm); for(i=0;i<cm;i++) orderv[i]=i;


  for(i=j=0;i<kk;i++) j+=bb[i];
  printf("\n# start generating %d repeat%s of\
 total %d replicates for %d items",
	 num_repeat,(num_repeat>1)?"s":"",j,cm); fflush(STDOUT);
  t0=get_time();
  for(irep=0;irep<num_repeat;irep++) {
    for(i=0;i<cm;i++) {
      deg=(int)degvec[i]; nc=ncvec[i];
      sg=ncvec[i]-flucscale; /* mean =  flucscale +  sg */
      for(j=0;j<deg;j++) {
	while((x=runif())==0.0);
	xbuf[j]=sg-flucscale*log(x); /* exponentinal(flucscale)+sg */
      }
      mydprintf(1,"\n# $d %d: ",irep,i);
      x=0.0; for(j=0;j<deg;j++) {
	x+=xbuf[j];
	mydprintf(1,"%g ",xbuf[j]);
      }
      obsvec[i]=x/deg; /* average */
      if(sw_nonpara) {
	/* non-parametric */
	for(k=0;k<kk;k++) {
	  bn=deg*rr1[k];
	  rr[k]=(double)bn / deg;
	  r2=1.0/bn; xp=statmats[k][i];
	  for(ib=0;ib<bb[k];ib++) {
	    x=0.0; for(j=0;j<bn;j++) x+=xbuf[(int)(runif()*deg)];
	    xp[ib]=r2*x;
	  }
	}
      } else {
	/* parametric */
	for(k=0;k<kk;k++) {
	  bn=deg*rr1[k];
	  rr[k]=(double)bn / deg;
	  r2=obsvec[i]/bn; xp=statmats[k][i];
	  for(ib=0;ib<bb[k];ib++) {
	    x=0.0; for(j=0;j<bn;j++) {
	      while((y=runif())==0.0); x+=-log(y);}
	    xp[ib]=r2*x;
	  }
	}
      }
      putdot();
    }
    if(sw_repeat) {
      cbuf=fnamenum(fname_rep,irep,num_repeat);
      write_rep(cbuf); FREE(cbuf);
    } else write_rep(fname_rep);    
  }
  t1=get_time();
  printf("\n# time elapsed for replicate generation t=%g sec",t1-t0);
  fflush(STDOUT);

  return 0;
}

int do_rmtnorm()
{
  int i,j,k,ib,irep;
  FILE *fp;
  char *cbuf;
  double *xp,*meanvec;
  double mn,r2,t0,t1;

  if(fname_vt){ 
    fp=openfp(fname_vt,fext_vt,"r",&cbuf);
    printf("\n# reading %s",cbuf);
  } else {
    fp=STDIN;
    printf("\n# reading from stdin\n");
  }
  mm=0; meanvec=fread_vec(fp,&mm);
  if(fname_vt) {fclose(fp); FREE(cbuf);}

  printf("\n# M:%d",mm);
  printf("\n# OBS:"); for(i=0;i<mm;i++) printf("%g ",meanvec[i]);
  datvec=new_vec(mm);
  repmats=NEW_A(kk,double**);
  for(i=0;i<kk;i++) repmats[i]=new_mat(mm,bb[i]);

  for(i=j=0;i<kk;i++) j+=bb[i];
  printf("\n# start generating %d repeat%s of\
 total %d replicates for %d items",
	 num_repeat,(num_repeat>1)?"s":"",j,mm); fflush(STDOUT);
  t0=get_time();
  for(irep=0;irep<num_repeat;irep++) {
    for(i=0;i<mm;i++) datvec[i]=meanvec[i]+flucscale*rnorm();
    for(k=0;k<kk;k++) {
      r2=1.0/sqrt(rr[k]);
      for(i=0;i<mm;i++) {
	mn=datvec[i]; xp=repmats[k][i];
	for(ib=0;ib<bb[k];ib++) xp[ib]=mn+r2*rnorm();
      }
      putdot();
    }
    if(sw_repeat) {
      cbuf=fnamenum(fname_rmt,irep,num_repeat);
      write_rmt(cbuf); FREE(cbuf);
    } else write_rmt(fname_rmt);    
  }
  t1=get_time();
  printf("\n# time elapsed for bootstrap t=%g sec",t1-t0);
  fflush(STDOUT);

  return 0;
}

int write_rep(char *fname)
{
  FILE *fp;
  char *cbuf;
  int i;

  if(fname) { /* binary write to file */
    fp=openfp(fname,fext_rep,"wb",&cbuf);
    printf("\n# writing %s",cbuf);
    fwrite_bivec(fp,orderv,cm); fwrite_bvec(fp,obsvec,cm);
    fwrite_bvec(fp,rr,kk); fwrite_bivec(fp,bb,kk);
    fwrite_bi(fp,kk);
    for(i=0;i<kk;i++) {
      fwrite_bmat(fp,statmats[i],cm,bb[i]); putdot();
    }
    fclose(fp); FREE(cbuf); putchar(';');
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

int write_rmt(char *fname)
{
  FILE *fp;
  char *cbuf;
  int i;

  if(fname) { /* binary write to file */
    fp=openfp(fname,fext_rmt,"wb",&cbuf);
    printf("\n# writing %s",cbuf);
    fwrite_bvec(fp,datvec,mm);
    fwrite_bvec(fp,rr,kk); fwrite_bivec(fp,bb,kk);
    fwrite_bi(fp,kk);
    for(i=0;i<kk;i++) {
      fwrite_bmat(fp,repmats[i],mm,bb[i]); putdot();
    }
    fclose(fp); FREE(cbuf); putchar(';');
  } else { /* ascii write to stdout */
    printf("\n# writing to stdout\n");
    printf("\n# OBS:\n"); write_vec(datvec,mm);
    printf("\n# R:\n"); write_vec(rr,kk);
    printf("\n# B:\n"); write_ivec(bb,kk);
    printf("\n# RMAT:\n");
    printf("%d\n",kk);
    for(i=0;i<kk;i++) {
      printf("\n## RMAT[%d]:\n",i); write_mat(repmats[i],mm,bb[i]);
    }
  }

  return 0;
}
