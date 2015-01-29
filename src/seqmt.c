/*  seqmt.c : sequence to mt file converter  */

static const char rcsid[] = "$Id: seqmt.c,v 1.6 2010/01/29 16:46:04 shimo Exp $";

#include <stdio.h>
#include <stdlib.h>
#include "misc.h"
#include "freadmat.h"

void byebye() {error("error in command line");}

int main(int argc, char** argv)
{
  /* working variables */
  int i,j;
  char *fname_in=NULL, *fname_out=NULL;
  char *cbuf,*fext;
  FILE *fp;
  double **mat;
  int n,m;

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
  if(fname_in) {
    fp=openfp(fname_in,fext,"r",&cbuf);
    printf("\n# reading %s",cbuf);
  } else {
    fp=STDIN;
    printf("\n# reading from stdin");
  }

  /* read file */
  n=m=0;
  switch(seqmode) {
  case SEQ_MOLPHY: 
    mat = fread_mat_lls(fp, &m, &n); break;
  case SEQ_PAML: 
    mat = fread_mat_lfh(fp, &m, &n); break;
  case SEQ_PAUP: 
    mat = fread_mat_paup(fp, &m, &n); break;
  case SEQ_PUZZLE: 
    mat = fread_mat_puzzle(fp, &m, &n); break;
  case SEQ_PHYML: 
    mat = fread_mat_phyml(fp, &m, &n); break;
  case SEQ_MT: 
  default: 
    mat = fread_mat(fp, &m, &n); break;  
  }
  if(fname_in) {fclose(fp);  FREE(cbuf);}
  printf("\n# M:%d N:%d",m,n);

  /* write file */
  if(fname_in && !fname_out) fname_out=rmvaxt(fname_in);  
  if(fname_out) {
    fp=openfp(fname_out,fext_mt,"w",&cbuf);
    printf("\n# writing %s",cbuf);
  } else {
    fp=STDOUT;
    printf("\n# writing to stdout\n");
  }
  fwrite_mat(fp, mat, m, n);
  if(fname_out) {fclose(fp);  FREE(cbuf);}
  putchar('\n');

  exit(0);
}

