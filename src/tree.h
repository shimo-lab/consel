/*

  tree.h: header for tree.c

  Time-stamp: <2001-04-29 18:52:12 shimo>
  $Id: tree.h,v 1.1 2001/05/05 09:06:43 shimo Exp $

  shimo@ism.ac.jp 
  Hidetoshi Shimodaira

*/

typedef struct tnode {
  struct tnode *parent;
  struct tnode **child;
  int nchild;
  int id;
} Tnode;

typedef struct wnode {
  struct wnode *next;
  char *label;
} Wnode;

typedef struct snode {
  struct snode *next;
  char *set;
  int len;
} Snode;

typedef struct {
  int *ve;
  int len;
} ivec;

int fskipgetc(FILE *fp);
char *fread_word(FILE *fp);
int wordtoid(char *word, Wnode *list);
char *idtoword(int id, Wnode *list);
Tnode *fread_tree(FILE *fp, Wnode *list);
int fwrite_tree(FILE *fp, Tnode *tree, Wnode *list);
Tnode *fread_rtree(FILE *fp, Wnode *list);
int fwrite_rtree(FILE *fp, Tnode *tree, Wnode *list);
Snode *treetosplits(Tnode *tree);
ivec *splitstoids(Snode *splits, Snode *list);
ivec *unrootsplits(Snode *list);
ivec* isec_ivec(ivec **vecs, int m);
ivec* subt_ivec(ivec *x, ivec *y);
ivec* unin_ivec(ivec **vecs, int m);
int bsearch_ivec(ivec *iv, int t);

