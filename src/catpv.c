/*

  catpv.c : cat pv-files

  Time-stamp: <2001-05-30 10:08:07 shimo>

  shimo@ism.ac.jp 
  Hidetoshi Shimodaira

  typical usage:
  # foo.pv -> foo.out
  catpv foo

*/

#include <stdio.h>
#include <math.h>
#include "misc.h"

static const char rcsid[] = "$Id: catpv.c,v 1.4 2001/05/29 06:35:57 shimo Exp shimo $";

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

int *permrev(int *order, int *rev, int n)
{
  int i;
  if(!rev) rev=new_ivec(n);
  for(i=0;i<n;i++) rev[i]=i;
  isort(order,rev,n);
  for(i=0;i<n;i++) if(order[i]!=i) error("invalid perm vector");
  return rev;
}

int sw_verpose=0;
int sw_pau=1;
int sw_pmc=1;
int sw_cat=0;
int sw_prt=1;

char *fname_cat=NULL;
char *fext_cat=".out";

double ***pvmats, ***semats, ***auxmats;
double **obsvecs;
int **revords;

int trc=0; /* number of items for aggregate */
int trc1=0; /* start for aggregate */
int trc2=0; /* end for aggregate */


#define MCPVNUM 5
#define AUPVNUM 2
#define MCAUXNUM 0
#define AUAUXNUM 6
int main(int argc, char** argv)
{
  /* working variables */
  int i,j,k,ir,ifile,nfile,cm,pvnum,auxnum,jau,jmc,cm0=0;
  FILE *fp;
  char **fnamev,*cbuf;
  int *orderv; double *obsvec;
  double **pvmat,**semat,**auxmat;
  int sw_mc, sw_au;
  double **outmat;

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
    } else if(streq(argv[i],"-t")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%d",&trc) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"-i")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%d",&trc1) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"-o")) {
      sw_cat=1;
      if(i+1<argc) {fname_cat=argv[i+1]; i++;}
    } else if(streq(argv[i],"-e")) {
      sw_printse=1;
    } else if(streq(argv[i],"-v")) {
      sw_verpose=1;
    } else if(streq(argv[i],"--no_au")) {
      sw_pau=0;
    } else if(streq(argv[i],"--no_sh")) {
      sw_pmc=0;
    } else if(streq(argv[i],"--no_print")) {
      sw_prt=0;
    } else byebye();
  }

  if(nfile<1) {
    error("must specify input-file");
    byebye();
  }

  if(sw_cat) {
    pvmats=NEW_A(nfile,double**);
    semats=NEW_A(nfile,double**);
    auxmats=NEW_A(nfile,double**);
    obsvecs=NEW_A(nfile,double*);
    revords=NEW_A(nfile,int*);
  }

  for(ifile=0;ifile<nfile;ifile++) {
    fp=openfp(fnamev[ifile],fext_pv,"r",&cbuf);
    if(sw_prt) printf("\n# reading %s",cbuf);
    cm=pvnum=auxnum=0;
    orderv=fread_ivec(fp,&cm); obsvec=fread_vec(fp,&cm);
    pvmat=fread_mat(fp,&cm,&pvnum);
    semat=fread_mat(fp,&cm,&pvnum);
    auxmat=fread_mat(fp,&cm,&auxnum);
    if(!cm0 || cm<cm0) cm0=cm; /* take the min */

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

    sw_au = sw_au && sw_pau && sw_prt;
    sw_mc = sw_mc && sw_pmc && sw_prt;

    if(sw_prt) printf("\n# %4s %4s %6s", "rank","item","obs");
    if(sw_au) {
      print_pvname("au");
      print_pvname("bp");
    }
    if(sw_prt) printf(" |");
    if(sw_mc) {
      print_pvname("nbp");
      print_pvname("kh");
      print_pvname("sh");
      print_pvname("wkh");
      print_pvname("wsh");
    }
    if(sw_prt) printf(" |");
    if(sw_verpose && sw_au) {
      printf(" %6s %6s %2s %6s %6s %6s","pf","rss","df","d","c","th");
    }

    for(i=0;i<cm;i++) {
      if(sw_prt) printf("\n# %4d %4d %6.1f",i+1,orderv[i]+1,obsvec[i]);
      j=0;
      if(sw_au) {
	j=jau;
	print_pval(pvmat[i][j+0],semat[i][j+0]);
	print_pval(pvmat[i][j+1],semat[i][j+1]);
      }
      if(sw_prt) printf(" |");
      if(sw_mc) {
	j=jmc;
	print_pval(pvmat[i][j+0],semat[i][j+0]);
	print_pval(pvmat[i][j+1],semat[i][j+1]);
	print_pval(pvmat[i][j+2],semat[i][j+2]);
	print_pval(pvmat[i][j+3],semat[i][j+3]);
	print_pval(pvmat[i][j+4],semat[i][j+4]);
      }
      if(sw_prt) printf(" |");
      if(sw_verpose && sw_au) {
	printf(" %6.3f %6.3f %2.0f %6.3f %6.3f %6.3f",
	       auxmat[i][0],auxmat[i][1],auxmat[i][2],
	       auxmat[i][3],auxmat[i][4],auxmat[i][5]);
      } 
    }

    if(sw_cat){ /* saving the file */
      pvmats[ifile]=pvmat;  semats[ifile]=semat;
      auxmats[ifile]=auxmat; obsvecs[ifile]=obsvec;
      revords[ifile]=permrev(orderv,NULL,cm); free_ivec(orderv);
    } else { /* discard them */
      free_mat(pvmat); free_mat(semat); free_mat(auxmat);
      free_ivec(orderv); free_vec(obsvec);
    }
    fclose(fp);  FREE(cbuf);
    if(sw_prt) printf("\n");
  }

  if(sw_cat) { /* aggregating the results */
    if(trc1>0) trc1--; if(trc1>=cm0) trc1=cm0-1;
    if(!trc) trc=cm0;
    trc2=trc1+trc-1; if(trc2>=cm0) trc2=cm0-1;
    printf("\n# aggregating from item %d to item %d",trc1+1,trc2+1);

    if(fname_cat) {
      fp=openfp(fname_cat,fext_cat,"w",&cbuf);
      printf("\n# writing %s",cbuf);
    } else {
      fp=STDOUT;
      printf("\n# writing to stdout\n");
    }
    outmat=new_mat(pvnum*2+auxnum+1,nfile);
    for(i=trc1;i<=trc2;i++) {
      /* item (i+1) */
      for(ifile=0;ifile<nfile;ifile++) {
	ir=revords[ifile][i]; k=0;
	for(j=0;j<pvnum;j++) outmat[k++][ifile]=pvmats[ifile][ir][j];
	for(j=0;j<pvnum;j++) outmat[k++][ifile]=semats[ifile][ir][j];
	for(j=0;j<auxnum;j++) outmat[k++][ifile]=auxmats[ifile][ir][j];
	outmat[k++][ifile]=obsvecs[ifile][ir];
      }
      fwrite_mat(fp,outmat,pvnum*2+auxnum+1,nfile); putdot();
    }
    if(fname_cat) {fclose(fp); FREE(cbuf);}
  }
  printf("\n");

  return 0;
}
