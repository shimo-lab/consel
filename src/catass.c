/*
  catass.c : manipulate association files

  Time-stamp: <2001-06-22 14:05:28 shimo>

  shimo@ism.ac.jp 
  Hidetoshi Shimodaira

  typical usage:
  # foo1.ass foo2.ass foo3.ass -> goo.ass; merge
  catass foo1 foo2 foo3 goo
  # foo.ass, hoo.vt -> goo.ass; extract
  catass -x hoo foo goo
  # foo.ass -> goo.ass; intersection
  catass -i foo goo
  # foo.ass -> goo.ass; union
  catass -u foo goo
  # foo.ass -> goo.ass; not
  catass -n foo goo

*/

#include <stdio.h>
#include <stdlib.h>
#include "misc.h"

static const char rcsid[] = "$Id: catass.c,v 1.1 2001/06/22 05:05:55 shimo Exp $";

void putdot() {putchar('.'); fflush(STDOUT);}
void byebye() {error("error in command line");}

char *fname_vt = NULL;
char *fext_vt = ".vt";
char *fext_ass = ".ass";
char **fnamev_in, *fname_out;
int nfile;

/* input associations */
int *cmv;
int ***assvecv;
int **asslenv;

/* output association */
int cm;
int **assvec;
int *asslen;

/* extract */
int cm1;
int **assvec1;
int *asslen1;
int *extvec;
int sw_extid=0;

/* switches */
int sw_intersection=0;
int sw_union=0;
int sw_not=0;
int sw_maxm=0;

/* routines */
void union_ivec(int **vecs, int *lens, int m, int **outvec, int *outl);
void intersection_ivec(int **vecs, int *lens, int m, int **outvec, int *outl);
void subtract_ivec(int *vec1, int len1, int *vec2, int len2, 
		   int **veco, int *leno);

int main(int argc, char** argv)
{
  int i,j,k,ifile;
  FILE *fp;
  char *cbuf;
  int *vec,len;
  
  printf("# %s",rcsid);

  fnamev_in=NEW_A(argc-1,char*);
  nfile=0;

  /* args */
  for(i=1;i<argc;i++) {
    if(argv[i][0] != '-') {
      fnamev_in[nfile]=argv[i];
      nfile++;
    } else if(streq(argv[i],"-i")) {
      sw_intersection=1;
    } else if(streq(argv[i],"-u")) {
      sw_union=1;
    } else if(streq(argv[i],"-n")) {
      sw_not=1;
    } else if(streq(argv[i],"-m")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%d",&sw_maxm) != 1)
	byebye();
      i+=1;
    } else if(streq(argv[i],"-x")) {
      if(i+1>=argc) byebye();
      fname_vt=argv[i+1];
      i+=1;
    } else if(streq(argv[i],"-X")) {
      if(i+1>=argc) byebye();
      fname_vt=argv[i+1]; sw_extid=1;
      i+=1;
    } else byebye();
  }

  if(nfile<1) {
    error("must specify output-file");
    byebye();
  } else {
    nfile--;
    fname_out=fnamev_in[nfile];
  }

  if(nfile==0) {
    /* identity association */
    printf("\n# identity association");
    if(sw_maxm) cm=sw_maxm; else cm=1;
    assvec=NEW_A(cm,int*); asslen=NEW_A(cm,int);
    for(i=0;i<cm;i++) {
      asslen[i]=1;
      assvec[i]=NEW_A(1,int);
      assvec[i][0]=i;
    }
  } else {
    cmv=new_ivec(nfile); 
    assvecv=NEW_A(nfile,int**); asslenv=NEW_A(nfile,int*);

    /* reading the associations */
    for(ifile=0;ifile<nfile;ifile++) {
      fp=openfp(fnamev_in[ifile],fext_ass,"r",&cbuf);
      printf("\n# %3d %s",ifile+1,cbuf);
      cm=cmv[ifile]=fread_i(fp);
      assvecv[ifile]=NEW_A(cm,int*); asslenv[ifile]=NEW_A(cm,int);
      for(i=0;i<cm;i++) {
	asslenv[ifile][i]=0;
	assvecv[ifile][i]=fread_ivec(fp,&(asslenv[ifile][i]));
      }
      fclose(fp); FREE(cbuf);
    }
    /* join associations */
    for(cm=ifile=0;ifile<nfile;ifile++) cm+=cmv[ifile];
    assvec=NEW_A(cm,int*); asslen=NEW_A(cm,int);
    for(j=ifile=0;ifile<nfile;ifile++) {
      printf("\n# %3d %3d items =%4d -%4d",
	     ifile+1,cmv[ifile],j+1,j+cmv[ifile]);
      for(i=0;i<cmv[ifile];i++,j++) {
	asslen[j]=asslenv[ifile][i];
	assvec[j]=assvecv[ifile][i];
      }
    }
  }
  printf("\n# CM:%d",cm);

  k=-1;
  for(i=0;i<cm;i++) for(j=0;j<asslen[i];j++)
    if(assvec[i][j]>k) k=assvec[i][j];
  printf("\n# MAX:%d",k+1);
  if(!sw_maxm) sw_maxm=k+1;

  if(fname_vt) { /* extract */
    fp=openfp(fname_vt,fext_vt,"r",&cbuf);
    printf("\n# reading %s",cbuf);
    cm1=0; extvec=fread_ivec(fp,&cm1);
    fclose(fp); FREE(cbuf);
    if(sw_extid) for(i=0;i<cm1;i++) extvec[i]--;
    assvec1=NEW_A(cm,int*); asslen1=NEW_A(cm,int);
    printf("\n# extract %d items", cm1);
    for(i=0;i<cm1;i++) {
      assvec1[i]=assvec[extvec[i]]; asslen1[i]=asslen[extvec[i]];
      printf("\n# %3d <- %3d",i+1,extvec[i]+1);
    }
    printf("\n# CM:%d",cm1);
    cm=cm1; assvec=assvec1; asslen=asslen1;
  }

  if(sw_union) {
    /* union */
    union_ivec(assvec,asslen,cm,&vec,&len);
    assvec[0]=vec; asslen[0]=len; cm=1; 
    printf("\n# union");
  }

  if(sw_intersection) {
    /* intersection */
    intersection_ivec(assvec,asslen,cm,&vec,&len);
    assvec[0]=vec; asslen[0]=len; cm=1; 
    printf("\n# intersection");
  }

  if(sw_not) {
    /* not */
    printf("\n# not");
    vec=new_ivec(sw_maxm);
    for(i=0;i<cm;i++) {
      for(j=0;j<sw_maxm;j++) vec[j]=j;
      subtract_ivec(vec,sw_maxm,assvec[i],asslen[i],
		    &(assvec[i]),&(asslen[i]));
    }
  }
  
  if(fname_out) {
    fp=openfp(fname_out,fext_ass,"w",&cbuf);
    printf("\n# writing %s",cbuf);
  } else {
    fp=STDOUT;
    printf("\n# writing to stdout\n");
  }
  printf("\n# CM:%d",cm);
  fprintf(fp,"\n# ASSOCIATIONS\n%d\n",cm);
  for(i=0;i<cm;i++) {
    printf("\n# %3d (%d)",i+1,asslen[i]);
    fprintf(fp,"\n# item: %d\n",i);
    fwrite_ivec(fp,assvec[i],asslen[i]);
  }

  if(fname_out) {fclose(fp); free(cbuf);}

  printf("\n# exit normally\n");
  exit(0);
}

void union_ivec(int **vecs, int *lens, int m, int **outvec, int *outl)
{
  int i,j,k,len,*vec;
 
  for(len=i=0;i<m;i++) len+=lens[i];
  vec=new_ivec(len);
  for(j=i=0;i<m;i++) for(k=0;k<lens[i];k++,j++) vec[j]=vecs[i][k];
  isort(vec,NULL,len); /* then the output will naturally be sorted */
  for(i=0;i<len;i++) { /* "uniq" operation */
    for(j=i+1;j<len;j++) if(vec[i]!=vec[j]) break;
    if(j-i-1>0) for(k=j;k<len;k++) vec[k-(j-i-1)]=vec[k];
    len-=j-i-1;
  }

  *outvec=vec; *outl=len;
}

void intersection_ivec(int **vecs, int *lens, int m, int **outvec, int *outl)
{
  int i,j,k,l,*vec,len;

  len=lens[0]; vec=new_ivec(len);
  for(i=0;i<len;i++) vec[i]=vecs[0][i];

  for(i=1;i<m;i++) {
    for(j=0;j<len;j++) {
      for(k=0;k<lens[i];k++)
	if(vec[j]==vecs[i][k]) break;
      if(k==lens[i]) {
	for(l=j;l<len-1;l++) vec[l]=vec[l+1];
	len--; j--;
      }
    }
  }

  union_ivec(&vec,&len,1,outvec,outl);
  free_ivec(vec);
}

void subtract_ivec(int *vec1, int len1, int *vec2, int len2, 
		   int **veco, int *leno)
{
  int i,j,k, *vec, len;

  len=len1; vec=new_ivec(len);
  for(i=0;i<len;i++) vec[i]=vec1[i];

  for(i=0;i<len;i++) {
    for(j=0;j<len2;j++)
      if(vec[i]==vec2[j]) break;
    if(j<len2) {
      for(k=i;k<len-1;k++) vec[k]=vec[k+1];
      len--; i--;
    }
  }
  
  *veco=vec; *leno=len;
}
