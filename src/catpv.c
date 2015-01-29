/*

  catpv.c : cat pv-files

  Time-stamp: <2011-01-25 16:47:31 shimo>

  shimo@ism.ac.jp 
  Hidetoshi Shimodaira

  typical usage:
  # foo.pv -> foo.out
  catpv foo

*/

#include <stdio.h>
#include <math.h>
#include "misc.h"

static const char rcsid[] = "$Id: catpv.c,v 1.15 2011/05/12 07:19:26 shimo Exp $";

char *fext_pv = ".pv";
int sw_help=0;
int sw_verpose=0;
int sw_pba=1;
int sw_pau=1;
int sw_pbp=1;
int sw_pmc=1;
int sw_cat=0;
int sw_prt=1;
int sw_sort=0;
int sw_cong=0;
int labeladd=1;
int sw_prank=0;
int sw_printse=0;

void putdot() {putchar('.'); fflush(STDOUT);}
void byebye() {error("error in command line");}

void print_real(double x)
{
  if(x<100.0 && x>-10.0 && fabs(x)>0.0005) printf(" %6.3f",x);
  else if(x==0.0) printf(" %6s","0");
  else if(fabs(x)<0.0005 || x>999999 || x<-99999) printf(" %6.0e",x);
  else  printf(" %6.0f",x);
}

void print_pvname(char *name)
{
  printf(" %6s",name);
  if(sw_printse)   printf(" %7s","(se)");
}

char indicate_sort(int col)
{
  if(col==sw_sort) return '+';
  else if(col==-sw_sort) return '-';
  else return ' ';
}

void print_pvcol(int col)
{
  printf(" %5d%c",col,indicate_sort(col));
  if(sw_printse)   printf(" %7s","");
}


void print_pval(double pv, double se)
{
  print_real(pv);
  if(sw_printse)   printf(" (%5.3f)",se);
}


char *fname_cat=NULL; char *fext_cat=".out";
char *fname_cong=NULL; char *fext_cong=".pv";

double ***pvmats, ***semats, ***auxmats;
double **obsvecs;
int **revords;

int trc=0; /* number of items for aggregate */
int trc1=0; /* start for aggregate */
int trc2=0; /* end for aggregate */


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

int main(int argc, char** argv)
{
  /* working variables */
  int i,j,k,ir,ifile,nfile,cm,pvnum,auxnum,jau,jmc,jbp,jba,cm0=0,cm1=0;
  FILE *fp;
  char **fnamev,*cbuf;
  int *orderv; double *obsvec; int *revordv=NULL;
  double **pvmat,**semat,**auxmat;
  double **pvmat2;
  int sw_bp,sw_ba,sw_mc,sw_au,fileid,outbit;
  double **outmat;
  double x,y,z,*xp;

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
    } else if(streq(argv[i],"-l")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%d",&labeladd) != 1)
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
    } else if(streq(argv[i],"-c")) {
      sw_cong=1;
      if(i+1<argc) {fname_cong=argv[i+1]; i++;}
    } else if(streq(argv[i],"-s")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%d",&sw_sort) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"-r")) {
      sw_prank=1;
    } else if(streq(argv[i],"-e")) {
      sw_printse=1;
    } else if(streq(argv[i],"-v")) {
      sw_verpose=1;
    } else if(streq(argv[i],"-h")) {
      sw_help=1;
    } else if(streq(argv[i],"--no_au")) {
      sw_pau=0;
    } else if(streq(argv[i],"--no_bp")) {
      sw_pbp=0;
    } else if(streq(argv[i],"--no_pp")) {
      sw_pba=0;
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

  if(sw_cat||sw_cong) {
    pvmats=NEW_A(nfile,double**);
    semats=NEW_A(nfile,double**);
    auxmats=NEW_A(nfile,double**);
    obsvecs=NEW_A(nfile,double*);
    revords=NEW_A(nfile,int*);
  }

  for(ifile=0;ifile<nfile;ifile++) {
    fp=openfp(fnamev[ifile],fext_pv,"r",&cbuf);
    if(sw_prt) printf("\n# reading %s",cbuf);
    fileid=fread_i(fp); if(fileid!=1) error("wrong file type");
    cm=pvnum=auxnum=0;
    orderv=fread_ivec(fp,&cm); obsvec=fread_vec(fp,&cm);
    outbit=fread_i(fp);
    pvmat=fread_mat(fp,&cm,&pvnum);
    semat=fread_mat(fp,&cm,&pvnum);
    auxmat=fread_mat(fp,&cm,&auxnum);
    if(!cm0 || cm<cm0) cm0=cm; /* take the min */
    if(cm>cm1) cm1=cm; /* take the max */

    revordv=new_ivec(cm);
    for(i=0;i<cm;i++) revordv[i]=i;
    if(sw_sort) {
      xp=new_vec(cm);
      /* sw_sort = 1:item, 2:lik, 3..(pvnum+2):pvs... */
      j=sw_sort; if(j<0) {j=-j;};
      if(j==1) { 
	for(i=0;i<cm;i++) xp[i]=(double)orderv[i];
      } else if(j==2) {
	for(i=0;i<cm;i++) xp[i]=obsvec[i];
      } else if(j<pvnum+3) {
	for(i=0;i<cm;i++) xp[i]=-pvmat[i][j-3];
      } else {
	error("-s # <= %d\n",pvnum+2);
      }
      if(sw_sort<0) for(i=0;i<cm;i++) xp[i]=-xp[i];
      sort(xp,revordv,cm);
      free_vec(xp);
    }

    sw_bp=sw_ba=sw_au=sw_mc=j=0;
    if(outbit & BPPVSW) {sw_bp=1; jbp=j; j+=BPPVNUM;}
    if(outbit & BAPVSW) {sw_ba=1; jba=j; j+=BAPVNUM;}
    if(outbit & MCPVSW) {sw_mc=1; jmc=j; j+=MCPVNUM;}
    if(outbit & AUPVSW) {sw_au=1; jau=j; j+=AUPVNUM;}
      
    sw_bp = sw_bp && sw_pbp && sw_prt;
    sw_ba = sw_ba && sw_pba && sw_prt;
    sw_au = sw_au && sw_pau && sw_prt;
    sw_mc = sw_mc && sw_pmc && sw_prt;

    if(sw_prt && (sw_sort||sw_verpose)) {
      printf("\n# %3d%c %3d%c %5d%c", 
	     0,indicate_sort(0),1,indicate_sort(1),2,indicate_sort(2));
      if(sw_au) 
	for(i=jau;i<jau+AUPVNUM;i++) print_pvcol(3+i);
      printf(" |");
      if(sw_bp) 
	for(i=jbp;i<jbp+BPPVNUM;i++) print_pvcol(3+i);
      if(sw_ba) 
	for(i=jba;i<jba+BAPVNUM;i++) print_pvcol(3+i);
      if(sw_mc) 
	for(i=jmc;i<jmc+MCPVNUM;i++) print_pvcol(3+i);
    }
    if(sw_prt) printf("\n# %4s %4s %6s", "rank","item","obs");
    if(sw_au) {
      print_pvname("au");
      print_pvname("np");
    }
    if(sw_prt) printf(" |");
    if(sw_bp) {
      print_pvname("bp");
    }
    if(sw_ba) {
      print_pvname("pp");
    }
    if(sw_mc) {
      print_pvname("kh");
      print_pvname("sh");
      print_pvname("wkh");
      print_pvname("wsh");
    }
    if(sw_prt) printf(" |");
    if(sw_verpose && sw_au) {
      printf(" %6s %6s %2s %6s %6s %6s %6s %6s","pf","rss","df","x","c",
	     "th","a","dim");
    }

    for(i=0;i<cm;i++) {
      ir=revordv[i];
      if(sw_prt) printf("\n# %4d %4d %6.1f",i+1,orderv[ir]+labeladd,obsvec[ir]); 
      j=0;
      if(sw_au) {
	j=jau;
	print_pval(pvmat[ir][j+0],semat[ir][j+0]);
	print_pval(pvmat[ir][j+1],semat[ir][j+1]);
      }
      if(sw_prt) printf(" |");
      if(sw_bp) {
	j=jbp;
	print_pval(pvmat[ir][j+0],semat[ir][j+0]);
      }
      if(sw_ba) {
	j=jba;
	print_pval(pvmat[ir][j+0],semat[ir][j+0]);
      }
      if(sw_mc) {
	j=jmc;
	print_pval(pvmat[ir][j+0],semat[ir][j+0]);
	print_pval(pvmat[ir][j+1],semat[ir][j+1]);
	print_pval(pvmat[ir][j+2],semat[ir][j+2]);
	print_pval(pvmat[ir][j+3],semat[ir][j+3]);
      }
      if(sw_prt) printf(" |");
      if(sw_verpose && sw_au) {
	print_real(auxmat[ir][0]);     /* pf */
	print_real(auxmat[ir][1]);     /* rss */
	printf(" %2.0f",auxmat[ir][2]); /* df */
	print_real(auxmat[ir][3]);     /* sid */
	print_real(auxmat[ir][4]);     /* cv */
	print_real(auxmat[ir][5]);     /* th */
	if(auxnum>6) {
	  if(auxmat[ir][6] != 0.0 && auxmat[ir][4] != 0.0)
	    x=0.5*(auxmat[ir][6]-1.0)/auxmat[ir][4];
	  else x=0.0;
	  print_real(x);                 /* a */
	  print_real(auxmat[ir][6]);     /* dim */
	}
      } 
    }

    if(sw_prank) {
      printf("\n\n# %6s %6s %6s","rank","ordr","item");
      for(i=0;i<cm;i++) {
      ir=revordv[i];
      printf("\n  %6d %6d %6d",i+1,ir+1,orderv[ir]+labeladd);
      }
    }

    if(sw_cat||sw_cong){ /* saving the file */
      pvmats[ifile]=pvmat;  semats[ifile]=semat;
      auxmats[ifile]=auxmat; obsvecs[ifile]=obsvec;
      revords[ifile]=revordv; free_ivec(orderv);
    } else { /* discard them */
      free_mat(pvmat); free_mat(semat); free_mat(auxmat);
      free_ivec(orderv); free_vec(obsvec);
    }
    fclose(fp);  FREE(cbuf);
    if(sw_prt) printf("\n");
  }

  if(sw_help && sw_prt) {
    printf("\n# ABBREVIATIONS");
    printf("\n# rank: rank of the item");
    printf("\n# ordr: the input order of the item");
    printf("\n# item: the label of the item");
    printf("\n# obs:  observed statistic value");
    printf("\n# --- p-values using the multiscale bootstrap ---");
    printf("\n# au:  the approximately unbiased test");
    printf("\n# np:  the naive p-value");
    printf("\n# --- p-values using the unscaled bootstrap ---");
    printf("\n# bp:  the bootstrap probability");
    printf("\n# pp:  the Bayesian posterior probability");
    printf("\n# kh:  the Kishino-Hasegawa test");
    printf("\n# sh:  the Shimodaira-Hasegawa test");
    printf("\n# wkh: the weighted KH-test");
    printf("\n# wsh: the weighted SH-test");
    printf("\n# --- details of the au test ---");
    printf("\n# pf:  diagnostic p-value");
    printf("\n# rss: fitting error");
    printf("\n# df:  fitting degrees of freedom");
    printf("\n# x:   signed distance");
    printf("\n# c:   curvature");
    printf("\n# th:  threshold for the region");
    printf("\n#\n# OPTIONS");
    printf("\n# -s: sort by the item");
    printf("\n# -v: show the details");
    printf("\n# -e: show the standard errors");
    printf("\n# --no_au: suppress au test");
    printf("\n# --no_bp: suppress bp test");
    printf("\n# --no_pp: suppress pp test");
    printf("\n# --no_sh: suppress sh test");
    printf("\n# --no_print: suppress printing");
    printf("\n# -o file: aggregating output");
    printf("\n# -c file: testing congruence");
    printf("\n");
  }

  if(trc1>0) trc1--; if(trc1>=cm0) trc1=cm0-1;
  if(!trc) trc=cm0;
  trc2=trc1+trc-1; if(trc2>=cm0) trc2=cm0-1;

  if(sw_cat) { /* aggregating the results */
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
      fwrite_mat(fp,outmat,pvnum*2+auxnum+1,nfile); 
      if(fname_cat) putdot();
    }
    if(fname_cat) {fclose(fp); FREE(cbuf);}
  }

  if(sw_cong) { /* congruence */
    cm=trc2-trc1+2;
    pvmat=new_mat(cm,pvnum);
    pvmat2=new_mat(cm,pvnum);
    semat=new_mat(cm,pvnum);
    auxmat=new_mat(cm,auxnum);
    orderv=new_ivec(cm);
    obsvec=new_vec(cm);
    
    /* calculating p-value */
    for(j=0;j<pvnum;j++) {
      z=0.0;
      for(i=0;i<cm-1;i++) {
	y=1.00;
	for(ifile=0;ifile<nfile;ifile++) {
	  ir=revords[ifile][i+trc1]; 
	  x=pvmats[ifile][ir][j];
	  if(x<y) y=x; /* y=min(x) */
	}
	pvmat[i][j]=y;
	if(y>z) z=y; /* z=max(y) */
      }
      pvmat[i][j]=z;
    }

    /* calculating log-likelihood */
    for(i=0;i<cm-1;i++) {
      y=0.0;
      for(ifile=0;ifile<nfile;ifile++) {
	ir=revords[ifile][i+trc1];
	x=obsvecs[ifile][ir];
	if(x<0.0) x=0.0;
	y+=x;
      }
      obsvec[i]=y;
      orderv[i]=i;
    }
    orderv[i]=i;

    /* sorting items */
    sort(obsvec,orderv,cm-1);
    x=obsvec[0];
    for(i=0;i<cm-1;i++) {
      obsvec[i]=obsvec[i]-x;
    }
    obsvec[0]=-obsvec[1];
    for(i=0;i<cm;i++) {
      for(j=0;j<pvnum;j++) 
	pvmat2[i][j]=pvmat[orderv[i]][j];
    }

    if(fname_cong) {
      fp=openfp(fname_cong,fext_cong,"w",&cbuf);
      printf("\n# writing %s",cbuf);
    } else {
      fp=STDOUT;
      printf("\n# writing to stdout\n");
    }

    fprintf(fp,"\n# ID:\n%d\n",1); 
    fprintf(fp,"\n# ITEM:\n"); fwrite_ivec(fp,orderv,cm);
    fprintf(fp,"\n# STAT:\n"); fwrite_vec(fp,obsvec,cm);
    fprintf(fp,"\n# BIT:\n%d\n",outbit); 
    fprintf(fp,"\n# PV:\n"); fwrite_mat(fp,pvmat2,cm,pvnum);  
    fprintf(fp,"\n# SE:\n"); fwrite_mat(fp,semat,cm,pvnum);  
    fprintf(fp,"\n# AX:\n"); fwrite_mat(fp,auxmat,cm,auxnum);  

    if(fname_cong) {fclose(fp); FREE(cbuf);}
  }

  printf("\n");
  return 0;
}
