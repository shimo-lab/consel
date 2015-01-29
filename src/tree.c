/*

  tree.c: tree handling routines

  Time-stamp: <2002-04-16 13:54:45 shimo>

  shimo@ism.ac.jp 
  Hidetoshi Shimodaira

*/

static const char rcsid[] = "$Id: tree.c,v 1.3 2002/04/16 15:58:50 shimo Exp $";

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "misc.h"
#include "tree.h"


int fskipgetc(FILE *fp)
{
    int c;
    fskipjunk(fp);
    c=getc(fp);
    if(c == EOF) error("unexpected EOF");
    return c;
}

char *fread_word(FILE *fp)
{
  static char buf[BUFSIZ];
  char *word;
  int i,c;
    
  fskipjunk(fp);
  for(i=0; i<BUFSIZ; i++){
    c=getc(fp);
    if(!isalnum(c) && c!='_') break;
    buf[i]=c;
  }
  if(i<BUFSIZ) ungetc(c,fp);
  else error("too long word");

  word=NEW_A(i+1,char);
  word[i]='\0';
  while(i-- > 0) word[i]=buf[i];

  return word;
}

int fread_branch(FILE *fp, double *bl)
{
  int c;
    
  fskipjunk(fp);
  c=getc(fp);
  if(c != ':') {
    ungetc(c,fp);
    return 0;
  }

  return fscanf(fp,"%lf",bl);
}


int wordtoid(char *word, Wnode *list)
{
  int id; /* id starts from 1 */
  Wnode *wp;

  /* the first node in the list is a dummy */
  for(id=1,wp=list; wp->next != NULL; id++,wp=wp->next) {
    if(strcmp(word,wp->next->label)==0) break;
  }
  if(wp->next == NULL) { /* not found */
    wp->next = NEW(Wnode);
    wp->next->next=NULL;
    wp->next->label=NEW_A(strlen(word)+1,char);
    strcpy(wp->next->label,word);
  }

  return id;
}

char *idtoword(int id, Wnode *list)
{
  int i;
  Wnode *wp;

  for(i=1,wp=list;(i<id) && (wp->next!=NULL);i++,wp=wp->next);
  if(wp->next==NULL) return NULL;
  return wp->next->label;
}

Tnode *fread_tree(FILE *fp, Wnode *list)
{
  int c;
  Tnode *tp;
  char *buf;
  double x;

  tp=NEW(Tnode);
  tp->parent=NULL; tp->child=NULL; tp->id=0;
  tp->nchild=0;

  c=fskipgetc(fp);
  if(c != '(') { /* leaves */
    ungetc(c,fp);
    buf=fread_word(fp);
    tp->id=wordtoid(buf,list);
    FREE(buf);
    fread_branch(fp,&x);
  } else { /* internal nodes */
    do {
      tp->child=RENEW_A(tp->child,tp->nchild+1,Tnode*);
      tp->child[tp->nchild]=fread_tree(fp,list);
      (tp->child[tp->nchild])->parent=tp;
      tp->nchild++;
      c=fskipgetc(fp);
    } while(c==',');
    if(c != ')') error(" ')' expected");
    fread_branch(fp,&x);
  }

  return tp;
}

int fwrite_tree(FILE *fp, Tnode *tree, Wnode *list)
{
  int i;
  char *cp;

  if(tree->nchild == 0) { /* leave */
    if(list!=NULL) {
      if((cp=idtoword(tree->id,list))) fputs(cp,fp);
      else error("non-id at leaves");
    } else {
      printf("%d",tree->id);
    }
  } else { /* internal node */
    putc('(',fp); i=0;
    do {
      fwrite_tree(fp,tree->child[i],list);
    } while ((++i < tree->nchild) && putc(',',fp));
    putc(')',fp);    
  }
  return 0;
}

Tnode *fread_rtree(FILE *fp, Wnode *list)
{
  int c;
  Tnode *tp;

  tp=fread_tree(fp,list);
  c=fskipgetc(fp);
  if(c != ';') error("\";\" is expected");
  while ( (c=getc(fp)) != '\n' ); /* skip until eol */

  return tp;
}

int fwrite_rtree(FILE *fp, Tnode *tree, Wnode *list)
{
  int i;
  i=fwrite_tree(fp,tree,list);
  putc(';',fp);

  return i;
}

Snode *treetosplits(Tnode *tree)
{
  int i,j,lenmax;
  Snode **spvec,*sq,*sp;

  sq=NEW(Snode);
  sq->next=NULL;
  if(tree->nchild==0) { /* leave */
    /* create a new split node with the leave */
    sq->len=tree->id+1;
    sq->set=NEW_A(sq->len,char);
    for(i=0;i<sq->len;i++) sq->set[i]=0;
    sq->set[tree->id]=1;
  } else { /* internal node */
    /* collect the split nodes from the children */
    spvec=NEW_A(tree->nchild, Snode*);
    lenmax=tree->id+1;
    for(i=0;i<tree->nchild;i++) {
      spvec[i]=treetosplits(tree->child[i]);
      if(spvec[i]->len>lenmax) lenmax=spvec[i]->len;
    }
    /* create a new split node from the children */
    sq->len=lenmax;
    sq->set=NEW_A(sq->len,char);
    for(i=0;i<sq->len;i++) sq->set[i]=0;
    sq->set[tree->id]=1;
    for(i=0;i<tree->nchild;i++) 
      for(j=0;j<spvec[i]->len;j++) sq->set[j]+=spvec[i]->set[j];
    /* adjoin sq, and spvec's */
    sq->next=spvec[0];
    for(i=0;i<tree->nchild-1;i++) {
      /* find the tail of spvec[i] */
      for(sp=spvec[i]; sp->next != NULL; sp=sp->next);
      sp->next = spvec[i+1];
    }
    FREE(spvec);
  }
  return sq;
}

int splittoid(Snode *split, Snode *list)
{
  int id; /* id starts from zero (has been changed from 1) */
  int i;
  Snode *sp,*sq;

  /* the first node in the list is a dummy */
  for(id=0,sq=list; (sp=sq->next)!=NULL; id++,sq=sp) {
    if(sp->len != split->len) continue;
    for(i=0;i<split->len;i++)
      if(sp->set[i] != split->set[i]) break;
    if(i==split->len) break; /* found */
  }
  if(sp==NULL) {/* not found */
    sp = sq->next = NEW(Snode);
    sp->next=NULL;
    sp->len=split->len;
    sp->set=NEW_A(sp->len,char);
    for(i=0;i<split->len;i++) sp->set[i]=split->set[i];
  }

  return id;
}

ivec *splitstoids(Snode *splits, Snode *list)
{
  int i,*iv;
  Snode *sp,*sq;
  ivec *ids;

  /* we IGNORE the number of unlabelled nodes */
  for(i=0,sp=splits; sp!=NULL; i++,sp=sp->next)
    if(sp->len>0) sp->set[0]=0; /* set to zero! */
  /* create a vector of the same size as splits */
  iv=NEW_A(i,int);
  ids=NEW(ivec); ids->len=i; ids->ve=iv;
  for(i=0,sp=splits; sp!=NULL; i++,sp=sq) {
    iv[i]=splittoid(sp,list);
    sq=sp->next;
    FREE(sp->set); FREE(sp); /* freeing up the inputs */
  }

  return ids;
}

ivec *unrootsplits(Snode *list)
{
  int i,j,k,len;
  char *set;
  Snode *sp,*sq,*sref;
  ivec *iv;

  /* use id=0 as the root of the tree (has been changed from 1) */
  sref=list->next;
  if(sref==NULL) error("empty list");
  set=NEW_A(sref->len,char);
  for(i=0,sp=list->next; sp!=NULL; i++,sp=sp->next);
  /* make a ivec */
  iv=NEW(ivec);  
  iv->len=i; iv->ve=NEW_A(iv->len,int);
  /* find the reverse split */
  for(i=0,sp=list->next; sp!=NULL; i++,sp=sp->next) {
    for(k=0;k<sref->len;k++) set[k]=sref->set[k];
    for(k=0;k<sp->len;k++) set[k]-=sp->set[k];
    for(len=sref->len; (len>0) && (set[len-1]==0);len--);
    for(j=0,sq=list->next; j<i; j++,sq=sq->next) {
      if(sq->len != len) continue;
      /* see if reversed split equlas sq */
      for(k=0;k<len;k++) if(sq->set[k]!=set[k]) break;
      if(k==len) break;
    }
    iv->ve[i]=j; /* no reversed split if i=j */
  }

  return iv;
}


ivec* isec_ivec(ivec **vecs, int m)
{
  int i,j,k,l,*vec,len;
  ivec *isec,iv;

  len=vecs[0]->len; vec=NEW_A(len,int);
  for(i=0;i<len;i++) vec[i]=vecs[0]->ve[i];

  for(i=1;i<m;i++) {
    for(j=0;j<len;j++) {
      for(k=0;k<vecs[i]->len;k++)
	if(vec[j]==vecs[i]->ve[k]) break;
      if(k==vecs[i]->len) {
	for(l=j;l<len-1;l++) vec[l]=vec[l+1];
	len--; j--;
      }
    }
  }

  isec=&iv; isec->len=len; isec->ve=vec;
  isec=unin_ivec(&isec,1); /* for "uniq" operation */
  FREE(vec);
  return isec;
}

ivec* subt_ivec(ivec *x, ivec *y)
{
  int i,j,k, *vec, len;
  ivec *sub;

  len=x->len; vec=NEW_A(len,int);
  for(i=0;i<len;i++) vec[i]=x->ve[i];

  for(i=0;i<len;i++) {
    for(j=0;j<y->len;j++)
      if(vec[i]==y->ve[j]) break;
    if(j<y->len) {
      for(k=i;k<len-1;k++) vec[k]=vec[k+1];
      len--; i--;
    }
  }

  sub=NEW(ivec);
  if((sub->len=len)>0) sub->ve=NEW_A(sub->len,int);
  else sub->ve=NULL;
  for(i=0;i<len;i++) sub->ve[i]=vec[i];
  FREE(vec);
  return sub;
}

ivec* unin_ivec(ivec **vecs, int m)
{
  int i,j,k,*vec,len;
  ivec *unin;

  len=0; for(i=0;i<m;i++) len+=vecs[i]->len;
  vec=NEW_A(len,int);
  j=0; for(i=0;i<m;i++) for(k=0;k<vecs[i]->len;k++,j++)
    vec[j]=vecs[i]->ve[k];
  isort(vec,NULL,len); /* then the output will naturally be sorted */
  for(i=0;i<len;i++) { /* "uniq" operation */
    for(j=i+1;j<len;j++) if(vec[i]!=vec[j]) break;
    if(j-i-1>0) for(k=j;k<len;k++) vec[k-(j-i-1)]=vec[k];
    len-=j-i-1;
  }
  unin=NEW(ivec);
  unin->len=len; unin->ve=NEW_A(unin->len,int);
  for(i=0;i<len;i++) unin->ve[i]=vec[i];
  FREE(vec);
  return unin;
}

/* binary search for a sorted ivec
   find k s.t. vec[k-1] <= t < vec[k]
 */
int bsearch_ivec(ivec *iv, int t)
{
  int i,i0,i1;
  int *vec,bb;

  vec=iv->ve; bb=iv->len;
  i0=0; i1=bb-1;
  if(t < vec[0]) return 0;
  else if(vec[bb-1] <= t) return bb;

  while(i1-i0>1) {
    i=(i0+i1)/2;
    if(vec[i] <= t) i0=i;
    else i1=i;
  }

  return i1;
}
