/*

  catrep.c : join and select rep files

  Time-stamp: <2001-06-23 13:24:08 shimo>

  shimo@ism.ac.jp 
  Hidetoshi Shimodaira

  typical usage:
  # foo1.rep foo2.rep -> foo.rep
  catrep foo1 foo2 foo
  # selecting replicates with scales in pa file
  catrep -p new.pa fooin fooout

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "misc.h"

static const char rcsid[] = "$Id: catrep.c,v 1.6 2001/08/10 06:03:10 shimo Exp $";

typedef struct {
  int kk; /* number of scales */
  int *bb; /* kk-vector: numnbers of the replicates */
  double *rr; /* kk-vector: scales = (rep size)/(orig size) */
  int cm; /* number of itmes (mm in rmt) */
  double ***mats; /* kk-array of {cm,mm} x bb[i] matrices */
  double *obs; /* {cm,rmt}-vector of observations */
  int *ord; /* cm-vector of id's (only in rep) */
} repstat;

int mode_rmt=0;

void printrepf(repstat *rp) {
  int i;
  printf("\n# K:%d",rp->kk);
  printf("\n# R:"); for(i=0;i<rp->kk;i++) printf("%g ",rp->rr[i]);
  printf("\n# B:"); for(i=0;i<rp->kk;i++) printf("%d ",rp->bb[i]);
  if(mode_rmt) printf("\n# MM:%d",rp->cm);
  else printf("\n# CM:%d",rp->cm);
}

void putdot() {putchar('.'); fflush(STDOUT);}
void byebye() {error("error in command line");}

int mode_nop=0;
int sw_sort=1;
char *fname_pa = NULL;
char **fnamev_in = NULL;
int nfile;
char *fname_out = NULL;
char *fext_pa = ".pa";
char *fext_rep = ".rep";
char *fext_rmt = ".rmt";
double rreps=0.5e-2;

repstat **repins=NULL;
int nrepins=0;
repstat reppar, repout, repin;
int sw_asciiout=0;
int sw_lmat=0;

int main(int argc, char** argv)
{
  /* working variables */
  int i,j,k,l,ifile;
  int *ibuf;
  FILE *fp;
  char *cbuf;
  repstat *rp;

  printf("# %s",rcsid);

  fnamev_in=NEW_A(argc-1,char*);
  nfile=0;
  /* args */
  for(i=1;i<argc;i++) {
    if(argv[i][0] != '-') {
      fnamev_in[nfile]=argv[i];
      nfile++;
    } else if(streq(argv[i],"-p")) {
      if(i+1>=argc) byebye();
      fname_pa=argv[i+1];
      i+=1;
    } else if(streq(argv[i],"-d")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%d",&debugmode) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"-e")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%lf",&rreps) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"-n")) {
      mode_nop=1;
    } else if(streq(argv[i],"-m")) {
      mode_rmt=1;
    } else if(streq(argv[i],"-a")) {
      sw_asciiout=1;
    } else if(streq(argv[i],"-L")) {
      sw_lmat=1;
    } else if(streq(argv[i],"--no_sort")) {
      sw_sort=0;
    } else byebye();
  }

  if(nfile<2) {
    error("must specify input-file and output-file");
    byebye();
  } else {
    nfile--;
    fname_out = fnamev_in[nfile];
  }

  /* reading parameter file */
  if(fname_pa!=NULL) {
    fp=openfp(fname_pa,fext_pa,"r",&cbuf);
    printf("\n# reading %s",cbuf);
    reppar.kk=0;
    reppar.rr=fread_vec(fp,&(reppar.kk));
    fclose(fp);
  }

  /* reading the replicates */
  repins=NEW_A(nfile,repstat*);

  for(ifile=0;ifile<nfile;ifile++){
    fp=openfp(fnamev_in[ifile],mode_rmt?fext_rmt:fext_rep,"rb",&cbuf);
    printf("\n# reading %s",cbuf);
    rp=repins[ifile]=NEW_A(1,repstat);
    rp->cm=rp->kk=0;
    if(!mode_rmt) rp->ord=fread_bivec(fp,&(rp->cm));
    rp->obs=fread_bvec(fp,&(rp->cm));
    rp->rr=fread_bvec(fp,&(rp->kk));
    rp->bb=fread_bivec(fp,&(rp->kk));
    i=fread_bi(fp);
    if(i!=rp->kk) error("wrong number of matrices");
    rp->mats=NEW_A(rp->kk,double**);
    if(!mode_nop){
      for(i=0;i<rp->kk;i++) {
	if(sw_lmat) rp->mats[i]=fread_blmat(fp,&(rp->cm),(rp->bb)+i);
	else rp->mats[i]=fread_bmat(fp,&(rp->cm),(rp->bb)+i); 
	putdot();
      }
    } else {
      printf(" --- skipped.");
    }
    printrepf(rp);
    FREE(cbuf); fclose(fp);
  }

  /* check if read files are compatible to each other */
  repin.cm=repins[0]->cm;
  repin.obs=new_vec(repin.cm);
  if(!mode_rmt) repin.ord=new_ivec(repin.cm);
  for(i=0;i<repin.cm;i++) {
    repin.obs[i]=repins[0]->obs[i];
    if(!mode_rmt) repin.ord[i]=repins[0]->ord[i];
  }
  for(k=ifile=0;ifile<nfile;ifile++){
    if(repins[ifile]->cm != repin.cm) error("cm mismatch");
    for(i=0;i<repin.cm;i++)
      if(fabs(repins[ifile]->obs[i] - repin.obs[i])>1e-6)
	error("obs mismatch");
    k+=repins[ifile]->kk;
  }
  /* find all the rr values */
  repin.rr=new_vec(k); repin.kk=0; 
  repin.bb=new_ivec(k);
  for(ifile=0;ifile<nfile;ifile++){
    rp=repins[ifile];
    for(i=0;i<rp->kk;i++){
      for(j=0;j<repin.kk;j++)
	if(fabs(repin.rr[j]-rp->rr[i])<rreps) break;
      if(j==repin.kk) {
	/* not found */
	repin.rr[j] = rp->rr[i];
	repin.bb[j] = rp->bb[i];
	repin.kk++;
      } else {
	/* found */
	repin.bb[j] += rp->bb[i];
      }
    }
  }
  /* memory allocation */
  repin.mats = NEW_A(repin.kk,double**);
  ibuf=new_ivec(repin.kk);
  for(i=0;i<repin.kk;i++) {
    if(sw_lmat) repin.mats[i]=new_lmat(repin.cm,repin.bb[i]);
    else repin.mats[i]=new_mat(repin.cm,repin.bb[i]);
    ibuf[i]=0;
  }
  /* now merging all the repmats */
  for(ifile=0;ifile<nfile;ifile++){
    rp=repins[ifile];
    for(i=0;i<rp->kk;i++){
      for(j=0;j<repin.kk;j++)
	if(fabs(repin.rr[j]-rp->rr[i])<rreps) break;
      if(j==repin.kk) {
	/* not found */
	error("internal error");
      } else {
	/* found */
	for(k=0;k<repin.cm;k++)
	  for(l=0;l<rp->bb[i];l++)
	    repin.mats[j][k][ibuf[j]+l]=rp->mats[i][k][l];
	ibuf[j]+=rp->bb[i];
      }
    }
  }
  /* check if everything is ok */
  if(nfile>1) {
    printf("\n# merged input");
    printrepf(&repin);
  }
  for(i=0;i<repin.kk;i++) 
    if(ibuf[i] != repin.bb[i]) error("internal error");
  FREE(ibuf);

  /* selecting rr */
  repout = repin;
  if(fname_pa!=NULL) {
    ibuf=new_ivec(reppar.kk);
    /* find the matchings */
    repout.kk=0;
    for(i=0;i<reppar.kk;i++) {
      for(j=0;j<repin.kk;j++)
	if(fabs(repin.rr[j] - reppar.rr[i])<rreps) break;
      if((ibuf[i]=j)!=repin.kk) repout.kk++; /* if found */
      else warning("scale not found for rr=%g",reppar.rr[i]);
    }
    repout.rr=new_vec(repout.kk);
    repout.bb=new_ivec(repout.kk);
    repout.mats = NEW_A(repout.kk,double**);
    /* finally copying the selected rr's */
    for(j=i=0;i<reppar.kk;i++)
      if(ibuf[i]!=repin.kk) {
	repout.rr[j]=repin.rr[ibuf[i]];
	repout.bb[j]=repin.bb[ibuf[i]];
	repout.mats[j]=repin.mats[ibuf[i]];
	j++;
      } 
    printf("\n# selected output");
    printrepf(&repout);
  }
  if(sw_sort) {
    printf("\n# sorting"); 
    for(i=0;i<repout.kk;i++) {
      for(j=0;j<repout.cm;j++)
	sort_vec(repout.mats[i][j],repout.bb[i]);
      putdot();
    }
  }

  /* output */
  fp=openfp(fname_out,mode_rmt?fext_rmt:fext_rep,
	    sw_asciiout?"w":"wb",&cbuf);
  printf("\n# writing %s",cbuf);
  if(sw_asciiout) {
    if(!mode_rmt) fwrite_ivec(fp,repout.ord,repout.cm);
    fwrite_vec(fp,repout.obs,repout.cm);
    fwrite_vec(fp,repout.rr,repout.kk);
    fwrite_ivec(fp,repout.bb,repout.kk);
    fwrite_i(fp,repout.kk);
    for(i=0;i<repout.kk;i++) {
      fwrite_mat(fp,repout.mats[i],repout.cm,repout.bb[i]); putdot();
    }
  } else {
    if(!mode_rmt) fwrite_bivec(fp,repout.ord,repout.cm);
    fwrite_bvec(fp,repout.obs,repout.cm);
    fwrite_bvec(fp,repout.rr,repout.kk);
    fwrite_bivec(fp,repout.bb,repout.kk);
    fwrite_bi(fp,repout.kk);
    for(i=0;i<repout.kk;i++) {
      fwrite_bmat(fp,repout.mats[i],repout.cm,repout.bb[i]); putdot();
    }
  }
  printf("\n# exit normally\n");
  exit(0);
}
