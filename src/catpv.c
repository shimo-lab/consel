/*

  catpv.c : cat pv-files

  Time-stamp: <2001-04-16 15:59:54 shimo>

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

int main(int argc, char** argv)
{
  /* working variables */
  int i,j,ifile,nfile,cm,mx1size;
  FILE *fp;
  char *fname;
  double *pvvec=NULL; /* pvalue */
  double *sevec=NULL; /* stderr of pv */
  int *ix1; double *ax1; double **mx1; /* auxiliary info */

  nfile=argc;

  for(ifile=1;ifile<nfile;ifile++) {
    fname=argv[ifile];
    fp=fopen(fname,"r");
    if(fp==NULL) error("cant open %s",fname);
    else printf("\n# read %s",fname);
    cm=mx1size=0;
    pvvec=fread_vec(fp,&cm);
    sevec=fread_vec(fp,&cm);
    ix1=fread_ivec(fp,&cm);
    ax1=fread_vec(fp,&cm);
    mx1=fread_mat(fp,&cm,&mx1size);

    printf("\n# %4s %4s %6s %6s %6s %6s %6s %6s",
	   "rank","item","obs","pval","pvse","perr","x","c");
    for(i=0;i<cm;i++) {
      printf("\n# %4d %4d %6.1f %6.3f %6.3f %6.3f %6.3f %6.3f",
	     i+1,ix1[i]+1,ax1[0]-ax1[i],pvvec[i],sevec[i],mx1[i][0],
	     mx1[i][3],mx1[i][4]);
    }
    fclose(fp);
    printf("\n");
  }

}
