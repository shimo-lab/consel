/*

  misc.h: header for misc.c

  Time-stamp: <2001-05-03 16:20:03 shimo>
  $Id: misc.h,v 1.3 2001/04/16 07:00:45 shimo Exp shimo $

  shimo@ism.ac.jp 
  Hidetoshi Shimodaira

*/

#define STDERR stdout
#define STDIN stdin
#define STDOUT stdout

/* system */
extern int debugmode;
void dprintf(int n, char *fmt, ...);
void error(char *fmt, ...);
void warning(char *fmt, ...);
double get_time(void);

/* memory */
void myfree(void *ptr);
void *myalloc(size_t size);
void *myrealloc(void *old, size_t size);
#define FREE myfree
#define MALLOC myalloc
#define REALLOC myrealloc
#define NEW_A(m,t) (t*)MALLOC((m)*sizeof(t))
#define NEW(t) NEW_A(1,t)
#define RENEW_A(A,m,t) (t*)REALLOC(A,(m)*sizeof(t))
#define FREE_A(A,t) FREE(A);
#define new_vec(m) NEW_A(m,double)
#define new_ivec(m) NEW_A(m,int)
#define renew_vec(A,m) RENEW_A(A,m,double)
#define renew_ivec(A,m) RENEW_A(A,m,int)
#define free_vec(A) FREE_A(A,double)
#define free_ivec(A) FREE_A(A,int)

double **renew_mat(double **old, int m, int n);
double **new_mat(int m, int n);
void free_mat(double **buf);

/* ascii read/write */
int fskipjunk(FILE *fp);
#define skipjunk() fskipjunk(STDIN)
int streq(char const *s1,char const *s2);
char *mstrcat(char *str1, char *str2);
FILE *openfp(char *name, char *ext, char *mode, char **fnamep);
int fread_i(FILE *fp);
#define read_i() fread_i(STDIN)
double fread_d(FILE *fp);
#define read_d() fread_d(STDIN)
double *fread_vec(FILE *fp, int *mp);
#define read_vec(mp) fread_vec(STDIN,mp)
int *fread_ivec(FILE *fp, int *mp);
#define read_ivec(mp) fread_ivec(STDIN,mp)
double **fread_mat(FILE *fp, int *mp, int *np);
#define read_mat(mp,np) fread_mat(STDIN,mp,np)
double **freread_mat(FILE *fp, int *mp, int *np, double **buf);
#define reread_mat(mp,np,buf) fread_mat(STDIN,mp,np,buf)
int fwrite_i(FILE *fp, int x);
#define write_i(x) fwrite_i(STDOUT,x)
int fwrite_d(FILE *fp, double x);
#define write_d(x) fwrite_d(STDOUT,x)
int fwrite_vec(FILE *fp, double *A, int m);
#define write_vec(A, m) fwrite_vec(STDOUT,A,m)
int fwrite_ivec(FILE *fp, int *A, int m);
#define write_ivec(A,m) fwrite_ivec(STDOUT,A,m)
int fwrite_mat(FILE *fp, double **A, int m, int n);
#define write_mat(A,m,n) fwrite_mat(STDOUT,A,m,n)

/* binary read/write */
int fread_bi(FILE *fp);
int fread_bd(FILE *fp);
int *fread_bivec(FILE *fp, int *mp);
double *fread_bvec(FILE *fp, int *mp);
double **fread_bmat(FILE *fp, int *mp, int *np);
double **freread_bmat(FILE *fp, int *mp, int *np, double **old);
int fwrite_bi(FILE *fp, int x);
int fwrite_bd(FILE *fp, double x);
int fwrite_bivec(FILE *fp, int *A, int m);
int fwrite_bvec(FILE *fp, double *A, int m);
int fwrite_bmat(FILE *fp, double **A, int m, int n);

/* misc */
void luinverse(double **omtrx, double **imtrx, int size); /* From Molphy */
void sort(double *xve, int *order, int dim);
void isort(int *xve, int *order, int dim);
void psort(void **xve, int *order, int dim,
	   int (*compar)(void *, void *));
int *perm_ivec(int *px, int *iv, int n);

double sym_mat(double **mat, int m);
double *lsfit(double **X, double *Y, double *W,
	      int m, int n,
	      double *beta, double *rss, double ***acmat);
int sort_vec(double *v, int n);
int argmin_vec(double *vec, int n);
int argmax_vec(double *vec, int n);

