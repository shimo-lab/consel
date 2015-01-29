/*

  misc.c: miscellaneous functions

  Time-stamp: <2011-01-25 16:57:00 shimo>

  shimo@ism.ac.jp 
  Hidetoshi Shimodaira

  Meschach library:  v_sort, skipjunk (David E. Stewart)

*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <time.h>
#include "misc.h"

static const char rcsid[] = "$Id: misc.c,v 1.12 2011/05/12 07:21:18 shimo Exp $";

/*
  error message handling
*/

int debugmode=0;
void mydprintf(int n, char *fmt, ...)
{
  va_list args;
  va_start(args,fmt);
  if(debugmode>=n) {
    vprintf(fmt, args);  
  }
  va_end(args);
}

void error(char *fmt, ...)
{
  va_list args;
  va_start(args,fmt);
  printf("\n#! Terminated: ");
  vprintf(fmt, args);
  printf(".\n");
  va_end(args);
  exit(1);
}

void warning(char *fmt, ...)
{
  va_list args;
  va_start(args,fmt);
  printf("\n#! Warning: ");
  vprintf(fmt, args);
  printf(". ");
  va_end(args);
}

/*
  get time from the system clock
*/

double get_time(void)
{
  clock_t t;
  double x;

  t = clock();
  x = (double) t / CLOCKS_PER_SEC;
  return x;
}

/*
  date
 */

long get_date(void)
{
  time_t t;

  t = time(0);
  return (long)t;
}

/*
  memory allocation: free, alloc, realloc
*/
void myfree(void *ptr)
{
  free(ptr);
}

void *myalloc(size_t size)
{
    void *ptr;
    if(size==0) warning("zero size in malloc");
    ptr=calloc(size,1); /* use calloc to initialize for safty */
    if(ptr == NULL) error("cant allocate memory size=%d",size);
    return ptr;
}

void *myrealloc(void *old, size_t size)
{
    void *ptr;
    ptr=realloc(old,size);
    if(ptr == NULL) error("cant reallocate memory size=%d",size);
    return ptr;
}

/*
  matrix
*/
double **new_mat(int m, int n)
{
  double **base,*xp;
  int i;

  if(m*n==0) return NULL;
  xp = new_vec(m*n);
  base = (double **)MALLOC(sizeof(double*) * m);
  for(i=0;i<m;i++) base[i] = xp + i*n;

  return base;
}

double **renew_mat(double **old, int m, int n)
{
  double *buf,**base;
  int i;

  if(old) buf=old[0]; else buf=NULL;
  buf=(double *)REALLOC(buf,sizeof(double)*(m*n));
  base=(double **)REALLOC(old,sizeof(double*) * m);
  for(i=0;i<m;i++) base[i] = buf + i*n;

  return base;
}

void free_mat(double **buf)
{
  if(buf && buf[0]) FREE(buf[0]);
  if(buf) FREE(buf);
}

/* large matrix (useful when large linear memory is not allocated */
double **new_lmat(int m, int n)
{
  double **base;
  int i;

  if(m*n==0) return NULL;
  base = (double **)MALLOC(sizeof(double*) * m);
  for(i=0;i<m;i++) base[i] = new_vec(n);

  return base;
}

void free_lmat(double **buf, int m)
{
  int i;
  if(!buf) return;
  for(i=0;i<m;i++) FREE(buf[i]);
  FREE(buf);
}



/* skipjunk -- skips white spaces and strings of the form #....\n
   Here .... is a comment string 
   borrowed from from "Meschach" library
*/
int fskipjunk(FILE *fp)
{
  int c;
     
  for ( ; ; ) {       /* forever do... */
    /* skip blanks */
    do c = getc(fp);
    while ( isspace(c) );

    /* skip comments (if any) */
    if ( c == '#' )
      /* yes it is a comment (line) */
      while ( (c=getc(fp)) != 0xa ); /* UNIX '\n'->0xA but DOS '\n'->0xD 0xA */
    else {
      ungetc(c,fp);
      break;
    }
  }
  return 0;
}

/*
  string comparison returns 1 if they are identical otherwise 0
*/
int streq(char const *s1,char const *s2)
{
    int i=0;
    while(s1[i] == s2[i]) {
        if(!s1[i]) return 1;
        i++;
    }
    return 0;
}
/*
  catenate string with memory allocation
*/
char *mstrcat(char *str1, char *str2)
{
  int len1,len2;
  char *str;
  len1 = strlen(str1); len2 = strlen(str2);
  str = (char *)MALLOC(len1+len2+1);
  strcpy(str,str1);
  strcpy(str+len1,str2);
  return str;
}
/*
  check if same extension
*/
int chkext(char *name, char *ext)
{
  int i,len1,len2;
  len1=strlen(name); len2=strlen(ext);
  if(len1<len2) return 0; /* differ */
  for(i=0;i<len2;i++) if(name[len1-len2+i] != ext[i]) break;
  if(i != len2) return 0; /* differ */
  return 1; /* same */
}
/*
  check if any extension
*/
int chkaxt(char *name)
{
  int i,len;
  len=strlen(name);
  for(i=len-1;i>=0;i--) if(!isalnum(name[i])) break;
  if(i>0 && name[i]=='.' && 
     (name[i-1]!='/' && name[i-1]!='.')) return 1; /* some extension */
  return 0; /* no extension */
}
/*
  get basename by removing extension if any
*/
char *rmvaxt(char *name)
{
  int i,len;
  char *out;
  len=strlen(name);
  for(i=len-1;i>=0;i--) if(!isalnum(name[i])) break;
  if(i>0 && name[i]=='.' && 
     (name[i-1]!='/' && name[i-1]!='.')) 
    len=i; /* some extension */
  out=NEW_A(len+1,char);
  for(i=0;i<len;i++) out[i]=name[i];
  out[i]=0;
  return out;
}
/*
  open file
*/
FILE *openfp(char *name, char *ext, char *mode, char **fnamep)
{
  char *fname;
  FILE *fp;

  /*  if(chkext(name,ext)) fname = mstrcat(name,"");  */
  if(chkaxt(name)) fname = mstrcat(name,"");
  else fname = mstrcat(name,ext); 
  fp=fopen(fname,mode);
  if(fp==NULL) error("cant open %s",fname);
  if(fnamep != NULL) *fnamep = fname; else FREE(fname);
  return fp;
}

/*
  ascii read/write
*/

int fread_i(FILE *fp)
{
  int x;
  fskipjunk(fp);
  if(fscanf(fp,"%d",&x)!=1) error("cant read int");
  return x;
}

double fread_d(FILE *fp)
{
  double x;
  fskipjunk(fp);
  if(fscanf(fp,"%lf",&x)!=1) error("cant read double");
  return x;
}

/* handling of gap is contributed by Jeff Craig 20080827 */
double fread_d_paup2(FILE *fp)
{
  int c;
  fskipjunk(fp);
  c = getc(fp);
  if (c == '-') {
    return -1.0;
  } else {
    ungetc(c, fp);
    return fread_d(fp);
  }
}

char *fread_s(FILE *fp)
{
  static char buf[BUFSIZ];
  char *word;
  int i,c;
    
  fskipjunk(fp);
  for(i=0; i<BUFSIZ; i++){
    c=getc(fp);
    if(isspace(c)) break;
    buf[i]=c;
  }
  if(i<BUFSIZ) ungetc(c,fp);
  else error("too long word");

  word=NEW_A(i+1,char);
  word[i]='\0';
  while(i-- > 0) word[i]=buf[i];

  return word;
}

double *fread_vec(FILE *fp, int *mp)
{
  int i,m;
  double *A;

  m=fread_i(fp); /* number of items */
  if(*mp>0 && *mp != m) error("size mismatch in vec");

  if(m>0) {
    A = new_vec(m);
    for(i=0;i<m;i++) A[i]=fread_d(fp);
  } else A = NULL;

  *mp = m;
  return A;
}

int *fread_ivec(FILE *fp, int *mp)
{
  int i,m;
  int *A;

  m=fread_i(fp); /* number of items */
  if(*mp>0 && *mp != m) error("size mismatch in ivec");

  if(m>0) {
    A = new_ivec(m);
    for(i=0;i<m;i++) A[i]=fread_i(fp);
  } else A = NULL;

  *mp = m;
  return A;
}

char **fread_svec(FILE *fp, int *mp)
{
  int i,m;
  char **A;

  m=fread_i(fp); /* number of items */
  if(*mp>0 && *mp != m) error("size mismatch in svec");

  if(m>0) {
    A = NEW_A(m,char*);
    for(i=0;i<m;i++) A[i]=fread_s(fp);
  } else A = NULL;

  *mp = m;
  return A;
}

double **fread_mat(FILE *fp, int *mp, int *np)
{
  int i,j,m,n;
  double **A;

  m=fread_i(fp); /* number of items (rows) */
  n=fread_i(fp); /* number of samples (columns) */
  if(*mp>0 && *mp != m) error("size of rows mismatch in mat");
  if(*np>0 && *np != n) error("size of columns mismatch in mat");

  if(m*n==0) return NULL;

  A = new_mat(m,n);

  for(i=0;i<m;i++) for(j=0;j<n;j++) A[i][j]=fread_d(fp);

  *mp = m;
  *np = n;
  return A;
}

double **freread_mat(FILE *fp, int *mp, int *np, double **old)
{
  int i,j,m,n;
  double **A;

  m=fread_i(fp); /* number of items (rows) */
  n=fread_i(fp); /* number of samples (columns) */
  if(*mp>0 && *mp != m) error("size of rows mismatch in mat");
  if(*np>0 && *np != n) error("size of columns mismatch in mat");

  A = renew_mat(old,m,n);

  for(i=0;i<m;i++) for(j=0;j<n;j++) A[i][j]=fread_d(fp);

  *mp = m;
  *np = n;
  return A;
}

double **fread_lmat(FILE *fp, int *mp, int *np)
{
  int i,j,m,n;
  double **A;

  m=fread_i(fp); /* number of items (rows) */
  n=fread_i(fp); /* number of samples (columns) */
  if(*mp>0 && *mp != m) error("size of rows mismatch in mat");
  if(*np>0 && *np != n) error("size of columns mismatch in mat");

  if(m*n==0) return NULL;

  A = new_lmat(m,n);

  for(i=0;i<m;i++) for(j=0;j<n;j++) A[i][j]=fread_d(fp);

  *mp = m;
  *np = n;
  return A;
}


static int mcol=5;
static char    *format = "%14.9g ";

int fwrite_i(FILE *fp, int x)
{
  fprintf(fp,"%d\n",x);
  return 0;
}

int fwrite_d(FILE *fp, double x)
{
  fprintf(fp,"%g\n",x);
  return 0;
}

int fwrite_vec(FILE *fp, double *A, int m)
{
  int i,k;

  fprintf(fp,"#!VEC:\n%d\n",m);
  k=0;
  for(i=0;i<m;i++) {
    if(++k > mcol) {
      putc('\n',fp);
      k=1;
    }
    fprintf(fp,format, A[i]);
  }
  putc('\n',fp);
  return 0;
}

int fwrite_ivec(FILE *fp, int *A, int m)
{
  int i,k;

  fprintf(fp,"#!VEC:\n%d\n",m);
  k=0;
  for(i=0;i<m;i++) {
    if(++k > mcol) {
      putc('\n',fp);
      k=1;
    }
    fprintf(fp,"%10d ", A[i]);
  }
  putc('\n',fp);
  return 0;
}

int fwrite_mat(FILE *fp, double **A, int m, int n)
{
  int i,j,k;
  double *xp;

  fprintf(fp,"#!MAT:\n%d %d\n",m,n);

  if(A==NULL) return 0;

  for(i=0;i<m;i++) {
    fprintf(fp,"\n# row: %d\n",i);
    k=0;
    xp=A[i];
    for(j=0;j<n;j++) {
      if(++k > mcol) {
	putc('\n',fp);
	k=1;
      }
      fprintf(fp,format,xp[j]);
    }
    putc('\n',fp);
  }
  return 0;
}

/*
  binary read/write
*/

int fread_bi(FILE *fp)
{
  int x;
  if(fread(&x,sizeof(int),1,fp)!=1) error("cant read binary int");
  return x;
}

double fread_bd(FILE *fp)
{
  double x;
  if(fread(&x,sizeof(double),1,fp)!=1) error("cant read binary double");
  return x;
}

int *fread_bivec(FILE *fp, int *mp)
{
  int m;
  int *A;

  m=fread_bi(fp); /* number of itmes */
  if(*mp>0 && *mp != m) error("size mismatch in bivec");
  A = new_ivec(m);
  if(fread(A,sizeof(int),m,fp)!=m) error("cant read bivec");
  *mp = m;
  return A;
}

double *fread_bvec(FILE *fp, int *mp)
{
  int m;
  double *A;

  m=fread_bi(fp); /* number of itmes */
  if(*mp>0 && *mp != m) error("size mismatch in bvec");
  A = new_vec(m);
  if(fread(A,sizeof(double),m,fp)!=m) error("cant read bvec");
  *mp = m;
  return A;
}

double **fread_bmat(FILE *fp, int *mp, int *np)
{
  int m,n;
  double **A;

  m=fread_bi(fp); /* number of items (rows) */
  n=fread_bi(fp); /* number of samples (columns) */
  if(*mp>0 && *mp != m) error("size of rows mismatch in bmat");
  if(*np>0 && *np != n) error("size of columns mismatch in bmat");
  A = new_mat(m,n);
  if(fread(A[0],sizeof(double),m*n,fp)!=(m*n)) error("cant read bmat");
  *mp = m; *np = n;
  return A;
}

double **fread_blmat(FILE *fp, int *mp, int *np)
{
  int m,n,i;
  double **A;

  m=fread_bi(fp); /* number of items (rows) */
  n=fread_bi(fp); /* number of samples (columns) */
  if(*mp>0 && *mp != m) error("size of rows mismatch in bmat");
  if(*np>0 && *np != n) error("size of columns mismatch in bmat");
  A = new_lmat(m,n);
  for(i=0;i<m;i++)
    if(fread(A[i],sizeof(double),n,fp)!=n) error("cant read bmat");
  *mp = m; *np = n;
  return A;
}

double **freread_bmat(FILE *fp, int *mp, int *np, double **old)
{
  int m,n;
  double **A;

  m=fread_bi(fp); /* number of items (rows) */
  n=fread_bi(fp); /* number of samples (columns) */
  if(*mp>0 && *mp != m) error("size of rows mismatch in bmat");
  if(*np>0 && *np != n) error("size of columns mismatch in bmat");
  A = renew_mat(old,m,n);
  if(fread(A[0],sizeof(double),m*n,fp)!=(m*n)) error("cant read bmat");
  *mp = m; *np = n;
  return A;
}

int fwrite_bi(FILE *fp, int x)
{
  if(fwrite(&x,sizeof(int),1,fp)!=1) error("cant write binary int");
  return 0;
}

int fwrite_bd(FILE *fp, double x)
{
  if(fwrite(&x,sizeof(double),1,fp)!=1) error("cant write binary double");
  return 0;
}

int fwrite_bivec(FILE *fp, int *A, int m)
{
  fwrite_bi(fp,m);
  if(fwrite(A,sizeof(int),m,fp)!=m) error("cant write bivec");
  return 0;
}

int fwrite_bvec(FILE *fp, double *A, int m)
{
  fwrite_bi(fp,m);
  if(fwrite(A,sizeof(double),m,fp)!=m) error("cant write bvec");
  return 0;
}

int fwrite_bmat(FILE *fp, double **A, int m, int n)
{
  int i;
  fwrite_bi(fp,m);
  fwrite_bi(fp,n);
  for(i=0;i<m;i++)
    if(fwrite(A[i],sizeof(double),n,fp)!=n) error("cant write bmat");
  return 0;
}


/* borrowed from mesch and then corrected a bug (almost rewritten by shimo) */

#define	MAX_STACK	60

/* v_sort -- sorts vector x, and generates permutation that gives the order
	of the components; x = [1.3, 3.7, 0.5] -> [0.5, 1.3, 3.7] and
	the permutation is order = [2, 0, 1].
	-- if order is NULL on entry then it is ignored
	-- the sorted vector x is returned */
void sort(double *xve, int *order, int dim)
{
  double tmp, v;
  int i, j, l, r, tmp_i,k;
  int stack[MAX_STACK], sp;


  if(order != NULL) for(i=0;i<dim;i++) order[i]=i;
  if ( dim <= 1 ) return;
  sp = 0; l = 0; r = dim-1;
  for ( ; ; ) {
    while ( r > l ) {
      /* sort xve[l],xve[l+1],...,xve[r] */
      k=(r+l)/2; /* often a good divider */
      v = xve[k]; i = l-1; j = r+1;
      for ( ; ; ) {
	/* make 
	   xve[l],...,xve[j] <= v;
	   xve[j+1],...,xve[i-1] == v;
	   xve[i],...,xve[r] >= v;
	*/
	while (i<r && xve[++i]<=v);
	while (j>l && xve[--j]>=v);
	if ( i > j ) break;
	if(i==r) j=k; else if(j==l) i=k;
	tmp = xve[i]; xve[i] = xve[j]; xve[j] = tmp;
	if ( order != NULL ) {
	  tmp_i = order[i]; order[i] = order[j]; order[j] = tmp_i;}
	/* note that l <= k < r, so checking i==r must be first */
	if(i==r) {j=r-1; break;} else if(j==l) {i=l+1; break;}
      }

      if ( j-l > r-i ) {
	  stack[sp++] = l; stack[sp++] = j; l = i;
	} else {
	  stack[sp++] = i; stack[sp++] = r; r = j; 
	}
    }

    /* recursion elimination */
    if ( sp == 0 ) break;
    r = stack[--sp]; l = stack[--sp];
  }
}

void isort(int *xve, int *order, int dim)
{
  int tmp, v;
  int i, j, l, r, tmp_i,k;
  int stack[MAX_STACK], sp;


  if(order != NULL) for(i=0;i<dim;i++) order[i]=i;
  if ( dim <= 1 ) return;
  sp = 0; l = 0; r = dim-1;
  for ( ; ; ) {
    while ( r > l ) {
      /* sort xve[l],xve[l+1],...,xve[r] */
      k=(r+l)/2; /* often a good divider */
      v = xve[k]; i = l-1; j = r+1;
      for ( ; ; ) {
	/* make 
	   xve[l],...,xve[j] <= v;
	   xve[j+1],...,xve[i-1] == v;
	   xve[i],...,xve[r] >= v;
	*/
	while (i<r && xve[++i]<=v);
	while (j>l && xve[--j]>=v);
	if ( i > j ) break;
	if(i==r) j=k; else if(j==l) i=k;
	tmp = xve[i]; xve[i] = xve[j]; xve[j] = tmp;
	if ( order != NULL ) {
	  tmp_i = order[i]; order[i] = order[j]; order[j] = tmp_i;}
	/* note that l <= k < r, so checking i==r must be first */
	if(i==r) {j=r-1; break;} else if(j==l) {i=l+1; break;}
      }

      if ( j-l > r-i ) {
	  stack[sp++] = l; stack[sp++] = j; l = i;
	} else {
	  stack[sp++] = i; stack[sp++] = r; r = j; 
	}
    }

    /* recursion elimination */
    if ( sp == 0 ) break;
    r = stack[--sp]; l = stack[--sp];
  }
}

void mypsort(void **xve, int *order, int dim,int (*compar)(void *, void *))
{
  void *tmp, *v;
  int i, j, l, r, tmp_i,k;
  int stack[MAX_STACK], sp;


  if(order != NULL) for(i=0;i<dim;i++) order[i]=i;
  if ( dim <= 1 ) return;
  sp = 0; l = 0; r = dim-1;
  for ( ; ; ) {
    while ( r > l ) {
      /* sort xve[l],xve[l+1],...,xve[r] */
      k=(r+l)/2; /* often a good divider */
      v = xve[k]; i = l-1; j = r+1;
      for ( ; ; ) {
	/* make 
	   xve[l],...,xve[j] <= v;
	   xve[j+1],...,xve[i-1] == v;
	   xve[i],...,xve[r] >= v;
	*/
	while (i<r && (*compar)(xve[++i],v)<=0);
	while (j>l && (*compar)(xve[--j],v)>=0);
	if ( i > j ) break;
	if(i==r) j=k; else if(j==l) i=k;
	tmp = xve[i]; xve[i] = xve[j]; xve[j] = tmp;
	if ( order != NULL ) {
	  tmp_i = order[i]; order[i] = order[j]; order[j] = tmp_i;}
	/* note that l <= k < r, so checking i==r must be first */
	if(i==r) {j=r-1; break;} else if(j==l) {i=l+1; break;}
      }

      if ( j-l > r-i ) {
	  stack[sp++] = l; stack[sp++] = j; l = i;
	} else {
	  stack[sp++] = i; stack[sp++] = r; r = j; 
	}
    }

    /* recursion elimination */
    if ( sp == 0 ) break;
    r = stack[--sp]; l = stack[--sp];
  }
}

int *perm_ivec(int *px, int *iv, int n) 
{
  /*
    iv[i] := iv[px[i]]
   */

  int s,i,j;
  int x0;

  for(s=0;s<n;s++) {
    if(px[s]>=n) continue;

    x0=iv[s];
    i=s;   

    while( j = px[i], px[i] += n, j != s ) {
      iv[i] = iv[j];
      i = j;
    }

    iv[i] = x0;
  }
  for(s=0;s<n;s++) px[s] -= n;

  return iv;
}

int sort_vec(double *v, int n)
{
  int i;
  for(i=1;i<n;i++) if(v[i] < v[i-1]) break;
  if(i==n) return 1; /* already sorted */
  sort(v,NULL,n);
  return 0;
}


/*
  OTHERS
*/	


int argmin_vec(double *vec, int n)
{
  int i,k;
  double x;
  k=0; x=vec[0];
  for(i=1;i<n;i++) if(vec[i]<x) {x=vec[i]; k=i;}
  return k;
}

int argmax_vec(double *vec, int n)
{
  int i,k;
  double x;
  k=0; x=vec[0];
  for(i=1;i<n;i++) if(vec[i]>x) {x=vec[i]; k=i;}
  return k;
}


double fsquare(double x)
{
  return x*x;
}
