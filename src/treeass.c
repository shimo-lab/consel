/*
  treeass.c : find the associations of the trees

  Time-stamp: <2001-04-25 15:59:23 shimo>

  shimo@ism.ac.jp 
  Hidetoshi Shimodaira

  typical usage:
  # foo.tpl -> foo.ass
  treeass foo
  # foo.tpl -> aho.ass
  treeass foo aho
*/

#include <stdio.h>
#include "misc.h"
#include "tree.h"

static const char rcsid[] = "$Id$";

void putdot() {putchar('.'); fflush(STDOUT);}
void byebye() {error("error in command line");}

char *fname_tpl = NULL;
char *fname_ass = NULL;
char *fext_tpl = ".tpl";
char *fext_ass = ".ass";

int ntree,nleaf,nsplit; /* numbers of trees, leaves, and splits */
Tnode **treevec; /* list of trees */
ivec **treesplit; /* split decompositions for the trees */
Wnode wlist; /* word list for labels */
Snode slist; /* split list */
Snode **splitvec; /* pointers to slist */

ivec *splitrev; /* indicate the equialent revserse split */
ivec *splitcom; /* splits common to all the trees */
ivec *splitbase; /* base splits */

ivec **splittree; /* the table from split to tree */

int sw_nleaf=0;

int main(int argc, char** argv)
{
  /* working variables */
  int i,j,k,len;
  FILE *fp,*fpl,*fpr;
  Snode *sp;
  Wnode *wp;
  ivec *iv;
  int *ibuf;

  printf("# %s",rcsid);

  /* args */
  for(i=j=1;i<argc;i++) {
    if(argv[i][0] != '-') {
      switch(j) {
      case 1: fname_tpl=fname_ass=argv[i]; break;
      case 2: fname_ass=argv[i]; break;
      default: byebye();
      }
      j++;
    } else if(streq(argv[i],"-l")) {
      sw_nleaf=1;
    } else if(streq(argv[i],"-d")) {
      if(i+1>=argc ||
	 sscanf(argv[i+1],"%d",&debugmode) != 1)
	byebye();
      i+=1;
    } else byebye();
  }

  /* open log-file */
  fpl=STDOUT;

  /* read tree topologies */
  if(fname_tpl){
    fname_tpl=mstrcat(fname_tpl,fext_tpl);
    if((fp=fopen(fname_tpl,"r"))==NULL) error("cant open %s",fname_tpl);
    fprintf(fpl,"\n# reading %s",fname_tpl);
  } else {
    fprintf(fpl,"\n# reading from stdin"); fp=STDIN;
  }
  /* tpl file starts with the number of trees */
  if(sw_nleaf) nleaf=fread_i(fp);
  ntree=fread_i(fp);
  /* init vectors and lists */
  treevec=NEW_A(ntree,Tnode*);
  treesplit=NEW_A(ntree,ivec*);
  wlist.next=NULL; slist.next=NULL;
  /* read trees */
  for(i=0;i<ntree;i++) {
    treevec[i]=fread_rtree(fp,&wlist); /* read */
    sp=treetosplits(treevec[i]); /* get the split decomposition */
    treesplit[i]=splitstoids(sp,&slist); /* convert it to ivec */
  }
  fclose(fp);

  /* count leaves */
  for(i=0,wp=wlist.next; wp!=NULL; i++,wp=wp->next);
  if(sw_nleaf)
    if(nleaf!=i) warning("no. of leaf=%d is different from %d\n",
			 i,nleaf);
  nleaf=i;

  fprintf(fpl,"\n# %d trees, %d leaves",ntree,nleaf);

  /* find the equivalent reverse split for unrooted trees */
  splitrev=unrootsplits(&slist);
  nsplit=splitrev->len;
  splitvec=NEW_A(nsplit,Snode*); /* this is easier to access than slist */
  for(i=0,sp=slist.next; i<nsplit; i++,sp=sp->next)
    splitvec[i]=sp;
  fprintf(fpl,"\n# splits total: %d",nsplit);
  for(i=j=0;i<nsplit;i++) if(splitrev->ve[i]!=i) j++;
  fprintf(fpl,"\n# equivalent reverse splits: %d",j);

  /* split decompositions are converted to unrooted version */
  for(i=0;i<ntree;i++) {  
    for(j=0;j<treesplit[i]->len;j++)
      treesplit[i]->ve[j]=splitrev->ve[treesplit[i]->ve[j]];
    isort(treesplit[i]->ve,NULL,treesplit[i]->len); /* looks neat */
  }

  /* find common splits */
  splitcom=isec_ivec(treesplit,ntree);

  /* subtract */
  for(i=0;i<ntree;i++) {
    iv=subt_ivec(treesplit[i],splitcom);
    FREE(treesplit[i]->ve); FREE(treesplit[i]);
    treesplit[i]=iv;
  }

  /* union */
  splitbase=unin_ivec(treesplit,ntree);

  /* renumbering of the base-splits */
  for(i=0;i<ntree;i++)
    for(j=0;j<treesplit[i]->len;j++)
      treesplit[i]->ve[j]=bsearch_ivec(splitbase,treesplit[i]->ve[j])-1;

  /* make reverse table from base-split to tree */  
  splittree=NEW_A(splitbase->len,ivec*);
  ibuf=NEW_A(ntree,int); 
  for(i=0;i<splitbase->len;i++){
    len=0;
    for(j=0;j<ntree;j++) {
      for(k=0;k<treesplit[j]->len;k++)
	if(treesplit[j]->ve[k]==i) break;
      if(k<treesplit[j]->len) {
	ibuf[len++]=j;
      }
    }
    splittree[i]=NEW(ivec);
    splittree[i]->len=len;
    if(len>0) {
      splittree[i]->ve=NEW_A(len,int);
      for(j=0;j<len;j++) splittree[i]->ve[j]=ibuf[j];
    } else splittree[i]->ve=NULL;
  }
  FREE(ibuf);


  /****************************************/

  fprintf(fpl,"\n# %d base-splits, %d common-splits, %d root-split",
	  splitbase->len,splitcom->len-1,1);

  /* print split->tree association */
  fprintf(fpl,"\n\n# split->tree\n%d",splitbase->len);
  for(i=0;i<splitbase->len;i++) {
    fprintf(fpl,"\n%3d %3d ",i+1,splittree[i]->len);
    for(j=0;j<splittree[i]->len;j++) fprintf(fpl," %d",splittree[i]->ve[j]+1);
  }

  /* print tree->split association */
  fprintf(fpl,"\n\n# tree->split\n%d",ntree);
  for(i=0;i<ntree;i++) {
    fprintf(fpl,"\n%3d %3d  ",i+1,treesplit[i]->len);
    for(j=0;j<treesplit[i]->len;j++)
      fprintf(fpl,"%d ",treesplit[i]->ve[j]+1);
  }


  /* print the splits */
  fprintf(fpl,"\n\n# base splits\n%d %d",splitbase->len,nleaf);
  fprintf(fpl,"\n    ");
  for(i=1;i<=nleaf;i++) fprintf(fpl,"%d",i%10);
  for(i=0;i<splitbase->len;i++) {
    sp=splitvec[splitbase->ve[i]];
    fprintf(fpl,"\n%3d ",i+1);
    for(j=1;j<sp->len;j++) 
      if(sp->set[j]!=0) fprintf(fpl,"%d",sp->set[j]);
      else fprintf(fpl,"-");
    for(;j<=nleaf;j++) fprintf(fpl,"-");
  }
  

  /* print the common splits */
  fprintf(fpl,"\n\n# common splits\n%d %d",splitcom->len-1,nleaf);
  fprintf(fpl,"\n    ");
  for(i=1;i<=nleaf;i++) fprintf(fpl,"%d",i%10);
  for(i=1;i<splitcom->len;i++) { /* discard the first split, i.e. root */
    sp=splitvec[splitcom->ve[i]];
    fprintf(fpl,"\n%3d ",i);
    for(j=1;j<sp->len;j++) 
      if(sp->set[j]!=0) fprintf(fpl,"%d",sp->set[j]);
      else fprintf(fpl,"-");
    for(;j<=nleaf;j++) fprintf(fpl,"-");
  }

  /* print leaves */
  fprintf(fpl,"\n\n# leaves: %d",nleaf);
  fprintf(fpl,"\n%d\n",nleaf);
  for(i=0,wp=wlist.next; i<nleaf; i++,wp=wp->next)
    fprintf(fpl,"%3d %s\n",i+1,wp->label);

  /* print trees */
  fprintf(fpl,"\n\n# trees: %d",ntree);
  fprintf(fpl,"\n%d\n",ntree);
  for(i=0;i<ntree;i++) {
    fwrite_rtree(fpl,treevec[i],NULL); /* write to log */
    fprintf(fpl," %d\n",i+1); 
  }


  /* OUTPUT ASSOCIATION */
  if(fname_ass) {
    fname_ass = mstrcat(fname_ass,fext_ass);
    if((fpr=fopen(fname_ass,"w"))==NULL) error("cant open %s",fname_ass);
    fprintf(fpl,"\n# writing %s",fname_ass);
  } else {
    fpr=STDOUT;
    fprintf(fpl,"\n# writing to stdout\n");
  }
  /* print split->tree association */
  fprintf(fpr,"\n# SPLIT->TREE\n%d\n",splitbase->len);
  for(i=0;i<splitbase->len;i++) {
    fwrite_ivec(fpr,splittree[i]->ve,splittree[i]->len);
  }
  /* print tree->split association */
  fprintf(fpr,"\n# TREE->SPLIT\n%d\n",ntree);
  for(i=0;i<ntree;i++) {
    fwrite_ivec(fpr,treesplit[i]->ve,treesplit[i]->len);
  }
  if(fname_ass) fclose(fpr);

  printf("\n# exit normally\n");
  exit(0);
}

