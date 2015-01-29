/*

  catci.c : cat ci-files

  Time-stamp: <2001-06-27 09:52:53 shimo>

  shimo@ism.ac.jp 
  Hidetoshi Shimodaira

  typical usage:
  # display foo.ci
  catci foo

*/

#include <stdio.h>
#include <math.h>
#include "misc.h"

static const char rcsid[] = "$Id: catci.c,v 1.3 2001/08/10 06:11:44 shimo Exp $";

char *fext_ci = ".ci";

void byebye() {error("error in command line");}
void repchar(char c, int n) { int i; for(i=0;i<n;i++) putchar(c);}
  
int sw_verpose=0;
int sw_au=1;
int sw_bp=1;

int main(int argc, char** argv)
{
  /* working variables */
  int i,j,ifile,nfile,cm,w,w2;
  FILE *fp;
  char *cbuf;
  double *alphavec=NULL;
  double **cimat=NULL, **semat=NULL, **eimat=NULL;
  double **ci0mat=NULL, **se0mat=NULL, **ei0mat=NULL;
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
    } else if(streq(argv[i],"--no_au")) {
      sw_au=0;
    } else if(streq(argv[i],"--no_np")) {
      sw_bp=0;
    } else byebye();
  }

  for(ifile=0;ifile<nfile;ifile++) {
    fp=openfp(fnamev[ifile],fext_ci,"r",&cbuf);
    printf("\n# reading %s",cbuf);
    cm=nalpha=0;
    orderv=fread_ivec(fp,&cm); obsvec=fread_vec(fp,&cm);
    alphavec=fread_vec(fp,&nalpha);
    cimat=fread_mat(fp,&cm,&nalpha); 
    semat=fread_mat(fp,&cm,&nalpha);
    eimat=fread_mat(fp,&cm,&nalpha);
    ci0mat=fread_mat(fp,&cm,&nalpha); 
    se0mat=fread_mat(fp,&cm,&nalpha);
    ei0mat=fread_mat(fp,&cm,&nalpha);
    fclose(fp);

    printf("\n#"); repchar(' ',17);
    w=nalpha*(sw_verpose?17:7);
    w2=w/2-2;
    if(sw_au) {
      repchar('-',w2); printf(" au "); repchar('-',w-w2-4);
    }
    printf(" |");
    if(sw_bp) {
      repchar('-',w2); printf(" np "); repchar('-',w-w2-4);
    }
    printf("\n# %4s %4s","rank","item");
    printf(" %6s","obs");
    if(sw_au) {
      for(j=0;j<nalpha;j++) {
	printf(" %6.3f",alphavec[j]);
	if(sw_verpose) printf(" %4s %4s","se","ei");
      }
    }
    printf(" |");
    if(sw_bp) {
      for(j=0;j<nalpha;j++) {
	printf(" %6.3f",alphavec[j]);
	if(sw_verpose) printf(" %4s %4s","se","ei");
      }
    }
    for(i=0;i<cm;i++) {
      printf("\n# %4d %4d",i+1,orderv[i]+1);
      printf(" %6.1f",obsvec[i]);
      if(sw_au){
	for(j=0;j<nalpha;j++) {
	  printf(" %6.1f",cimat[i][j]);
	  if(sw_verpose) printf(" %4.1f %4.1f",semat[i][j],eimat[i][j]);
	}
      }
      printf(" |");
      if(sw_bp){
	for(j=0;j<nalpha;j++) {
	  printf(" %6.1f",ci0mat[i][j]);
	  if(sw_verpose) printf(" %4.1f %4.1f",se0mat[i][j],ei0mat[i][j]);
	}
      }
    }
    printf("\n");
    free_vec(alphavec);
    free_mat(cimat); free_mat(semat); free_mat(eimat);
    free_mat(ci0mat); free_mat(se0mat); free_mat(ei0mat);
    free_ivec(orderv); free_vec(obsvec);
  }
  return 0;
}
