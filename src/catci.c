/*

  catci.c : cat ci-files

  Time-stamp: <2001-05-05 15:24:02 shimo>

  shimo@ism.ac.jp 
  Hidetoshi Shimodaira

  typical usage:
  # foo.pv -> foo.out
  catpv foo

*/

#include <stdio.h>
#include <math.h>
#include "misc.h"

static const char rcsid[] = "$Id$";

char *fext_ci = ".ci";

void byebye() {error("error in command line");}

int sw_verpose=0;

int main(int argc, char** argv)
{
  /* working variables */
  int i,j,ifile,nfile,cm;
  FILE *fp;
  char *fname;
  double *alphavec=NULL;
  double **cimat=NULL, **semat=NULL, **eimat=NULL;
  int nalpha;
  int *orderv; double *obsvec; /* auxiliary info */
  char **fnamev;

  fnamev=NEW_A(argc-1,char*);
  nfile=0;

  /* args */
  for(i=1;i<argc;i++) {
    if(argv[i][0] != '-') {
      fnamev[nfile]=argv[i];
      nfile++;
    } else if(streq(argv[i],"-d")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%d",&debugmode) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"-v")) {
      sw_verpose=1;
    } else byebye();
  }

  for(ifile=0;ifile<nfile;ifile++) {
    fname=mstrcat(fnamev[ifile],fext_ci);
    fp=fopen(fname,"r");
    if(fp==NULL) error("cant open %s",fname);
    else printf("\n# read %s",fname);
    cm=nalpha=0;
    orderv=fread_ivec(fp,&cm); obsvec=fread_vec(fp,&cm);
    alphavec=fread_vec(fp,&nalpha);
    cimat=fread_mat(fp,&cm,&nalpha); 
    semat=fread_mat(fp,&cm,&nalpha);
    eimat=fread_mat(fp,&cm,&nalpha);
    fclose(fp);

    printf("\n# %4s %4s","rank","item");
    printf(" %6s","obs");
    for(j=0;j<nalpha;j++) {
      printf(" %6.3f",alphavec[j]);
      if(sw_verpose) printf(" %4s %4s","se","ei");
    }
    for(i=0;i<cm;i++) {
      printf("\n# %4d %4d",i+1,orderv[i]+1);
      printf(" %6.1f",obsvec[i]);
      for(j=0;j<nalpha;j++) {
	printf(" %6.1f",cimat[i][j]);
	if(sw_verpose) printf(" %4.1f %4.1f",semat[i][j],eimat[i][j]);
      }
    }
    printf("\n");
    free_vec(alphavec);
    free_mat(cimat); free_mat(semat); free_mat(eimat);
    free_ivec(orderv); free_vec(obsvec);
  }
  return 0;
}
