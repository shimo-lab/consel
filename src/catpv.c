/*

  catpv.c : cat pv-files

  Time-stamp: <2001-05-05 17:23:32 shimo>

  shimo@ism.ac.jp 
  Hidetoshi Shimodaira

  typical usage:
  # foo.pv -> foo.out
  catpv foo

*/

#include <stdio.h>
#include <math.h>
#include "misc.h"

static const char rcsid[] = "$Id: catpv.c,v 1.1 2001/04/16 06:59:57 shimo Exp shimo $";

char *fext_pv = ".pv";

void putdot() {putchar('.'); fflush(STDOUT);}
void byebye() {error("error in command line");}

int sw_printse=0;

void print_pvname(char *name)
{
  printf(" %6s",name);
  if(sw_printse)   printf(" %7s","(se)");
}

void print_pval(double pv, double se)
{
  printf(" %6.3f",pv);
  if(sw_printse)   printf(" (%5.3f)",se);
}

int sw_verpose=0;

#define MCPVNUM 4
#define AUPVNUM 2
#define MCAUXNUM 0
#define AUAUXNUM 6
int main(int argc, char** argv)
{
  /* working variables */
  int i,j,ifile,nfile,cm,pvnum,auxnum,jau,jmc;
  FILE *fp;
  char **fnamev,*cbuf;
  int *orderv; double *obsvec;
  double **pvmat,**semat,**auxmat;
  int sw_mc, sw_au;

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
    } else if(streq(argv[i],"-e")) {
      sw_printse=1;
    } else if(streq(argv[i],"-v")) {
      sw_verpose=1;
    } else byebye();
  }

  if(nfile<1) {
    error("must specify input-file");
    byebye();
  }
  
  for(ifile=0;ifile<nfile;ifile++) {
    fp=openfp(fnamev[ifile],fext_pv,"r",&cbuf);
    printf("\n# reading %s",cbuf);
    cm=pvnum=auxnum=0;
    orderv=fread_ivec(fp,&cm); obsvec=fread_vec(fp,&cm);
    pvmat=fread_mat(fp,&cm,&pvnum);
    semat=fread_mat(fp,&cm,&pvnum);
    auxmat=fread_mat(fp,&cm,&auxnum);

    sw_au=sw_mc=j=0;
    if((pvnum==MCPVNUM || pvnum==(MCPVNUM+AUPVNUM)) &&
	(auxnum==MCAUXNUM || auxnum==(MCAUXNUM+AUAUXNUM))) {
      sw_mc=1; jmc=j; j+=MCPVNUM;
    }
    if((pvnum==AUPVNUM || pvnum==(MCPVNUM+AUPVNUM)) &&
	      (auxnum==AUAUXNUM || auxnum==(MCAUXNUM+AUAUXNUM))) {
      sw_au=1; jau=j; j+=AUPVNUM;
    }
    if(sw_au + sw_mc == 0) {
      warning("file type unknown (pv:%d, aux:%d)",pvnum,auxnum);
    }

    printf("\n# %4s %4s %6s", "rank","item","obs");
    if(sw_au) {
      print_pvname("au");
      print_pvname("bp");
    }
    if(sw_mc) {
      print_pvname("bp0");
      print_pvname("kh");
      print_pvname("mc");
      print_pvname("ms");
    }
    if(sw_verpose && sw_au) {
      printf(" %6s %6s %2s %6s %6s %6s","pf","rss","df","x","c","th");
    }

    for(i=0;i<cm;i++) {
      printf("\n# %4d %4d %6.1f",i+1,orderv[i]+1,obsvec[i]);
      j=0;
      if(sw_au) {
	j=jau;
	print_pval(pvmat[i][j+0],semat[i][j+0]);
	print_pval(pvmat[i][j+1],semat[i][j+1]);
      }
      if(sw_mc) {
	j=jmc;
	print_pval(pvmat[i][j+0],semat[i][j+0]);
	print_pval(pvmat[i][j+1],semat[i][j+1]);
	print_pval(pvmat[i][j+2],semat[i][j+2]);
	print_pval(pvmat[i][j+3],semat[i][j+3]);
      }
      if(sw_verpose && sw_au) {
	printf(" %6.3f %6.3f %2.0f %6.3f %6.3f %6.3f",
	       auxmat[i][0],auxmat[i][1],auxmat[i][2],
	       auxmat[i][3],auxmat[i][4],auxmat[i][5]);
      } 
    }

    free_mat(pvmat); free_mat(semat); free_mat(auxmat);
    free_ivec(orderv); free_vec(obsvec);
    fclose(fp);  FREE(cbuf);
    printf("\n");
  }

  return 0;
}
