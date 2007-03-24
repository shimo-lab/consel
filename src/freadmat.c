/* freadmat.c May 29 2001 H.Mine */
/* modified by shimo May 29 */
/* $Id: freadmat.c,v 1.6 2005/09/20 07:57:16 shimo Exp shimo $ */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "misc.h"

int fskipline(FILE *fp)
{
  int c;

  while( 1 ) {
    c = getc(fp);
    if( ferror(fp) )
      return 1;
    if( feof(fp) )
      return 0;
    if( c == 0xa )
      return 0;
  }
}

int fskipword(FILE *fp)
{
  int c;

  while( 1 ) {
    c = getc(fp);
    if( ferror(fp) )
      return 1;
    if( feof(fp) )
      return 0;
    if(isspace(c))
      return 0;
  }
}

double **fread_mat_lls(FILE *fp, int *mp, int *np)
{
  int i,j,m,n;
  double **A;

  m=fread_i(fp); /* number of items (rows) */
  n=fread_i(fp); /* number of samples (columns) */
  fskipline(fp);
  if(*mp>0 && *mp != m) error("size of rows mismatch in mat");
  if(*np>0 && *np != n) error("size of columns mismatch in mat");

  if(m*n==0) return NULL;

  A = new_mat(m,n);

  for(i=0;i<m;i++) for(j=0;j<n;j++) A[i][j]=fread_d(fp);

  *mp = m;
  *np = n;
  return A;
}

double **fread_mat_puzzle(FILE *fp, int *mp, int *np)
{
  int i,j,m,n;
  double **A;

  m=fread_i(fp); /* number of items (rows) */
  n=fread_i(fp); /* number of samples (columns) */
  dprintf(2,"\nm=%d n=%d",m,n);
  fskipline(fp);
  if(*mp>0 && *mp != m) error("size of rows mismatch in mat");
  if(*np>0 && *np != n) error("size of columns mismatch in mat");

  if(m*n==0) return NULL;

  A = new_mat(m,n);

  for(i=0;i<m;i++) {
    fskipjunk(fp);
    fskipword(fp);
    dprintf(2,"\ni=%d",i+1);
    for(j=0;j<n;j++) {
      A[i][j]=fread_d(fp);
      dprintf(2,"\n%d %lg",j+1, A[i][j]);
    }
  }

  *mp = m;
  *np = n;
  return A;
}


double **fread_mat_lfh(FILE *fp, int *mp, int *np)
     /* modified by shimo */
{
  int i,j,m,n,t,p,k,jj;
  double **A,x;

  m=fread_i(fp); /* number of items (rows) */
  n=fread_i(fp); /* number of sites (columns for output) */
  p=fread_i(fp); /* number of patterns (columns for input) */
  dprintf(2,"\nm=%d n=%d p=%d",m,n,p);
  fskipline(fp);
  if(*mp>0 && *mp != m) error("size of rows mismatch in mat");
  if(*np>0 && *np != n) error("size of columns mismatch in mat");

  if(m*n==0) return NULL;

  A = new_mat(m,n);

  for(i=0;i<m;i++) {
    t = fread_i(fp);
    dprintf(2,"\ni=%d",t);
    if( t != i + 1) error("wrong row index");
    for(jj=j=0;j<p;j++) {
      t = fread_i(fp);
      dprintf(2,"\n%d",t);
      if( t != j + 1) error("wrong column index");
      t = fread_i(fp); /* number of repeats */
      x = fread_d(fp); /* site-wise log-likelihood */
      dprintf(2," %d %lg",t,x);
      fskipline(fp);
      for(k=0;k<t;k++) {
	if(jj>=n) error("too many repeats");
	A[i][jj++]=x;
      }
    }
    if(jj != n) error("too few repeats");
  }

  *mp = m;
  *np = n;
  return A;
}

int fread_line(FILE *fp, char *buff, int length)
{
  int i, c;

  for( i = 0; i < length; i++ ) {
    c = getc(fp);
    if( ferror(fp) || feof(fp) || c == 0xa )
      break;
    buff[i] = c;
  }
  return i;
}

#define ID_STRING "Single-site"
int skip_id(FILE *fp)
{
  char buff[sizeof ID_STRING];
  int r;

  while( !ferror( fp ) && !feof( fp ) ) {
    fskipjunk(fp);
    r = fread_line( fp, buff, sizeof ID_STRING );
    fskipline(fp);
    if( r == sizeof ID_STRING 
	&& strncmp( ID_STRING, buff, sizeof ID_STRING - 1 ) == 0 )
      return 1;
  }
  return 0;
}

int fread_i_noerror(FILE *fp)
{
  int x;
  fskipjunk(fp);
  x = 0;
  fscanf(fp,"%d",&x);
  return x;
}

int fread_d_noerror(FILE *fp)
{
  double x;
  fskipjunk(fp);
  x = 0.0;
  fscanf(fp,"%lf",&x);
  return x;
}

#define INIT_VEC_SIZE 500
double **fread_mat_paup1(FILE *fp, int *mp, int *np)
{
  int m,n,t,len;
  double **A;
  double *V;

  m = 0;
  len = INIT_VEC_SIZE;
  A = NULL; V = NULL;
  while( skip_id(fp) ) {
    for( n = 0; !ferror(fp) && !feof(fp); fskipjunk(fp) ) {
      t = fread_i_noerror(fp);
      if( t != n + 1 )
	break;
      if( V == NULL || t > len ) {
	len *= 2;
	V = (double *)renew_vec(V, len);
      }
      V[n++] = fread_d(fp);
    }
    if(*np>0 && *np != n) error("size of columns mismatch in mat");
    *np = n;
    A = (double **)renew_mat( A, m + 1, n );
    memcpy( A[m++], V, n * sizeof(double) );
  }
  if(*mp>0 && *mp != m) error("size of rows mismatch in mat");
  *mp = m;
  free( V );
  return A;
}

double **fread_mat_paup2(FILE *fp, int *mp, int *np)
/* by shimo (20020820) */
{
  int m,n,len,c;
  double **A;
  double *V;

  n = m = 0;
  len = INIT_VEC_SIZE;
  A = NULL; V = NULL;
  
  while(!ferror(fp) && !feof(fp)) {
    c=getc(fp);
    if(c=='#') {fskipline(fp); continue; }
    if(c=='\t') {
      fread_i_noerror(fp);
      if( V == NULL || n >= len ) {
	len *= 2;
	V = (double *)renew_vec(V, len);
      }
      V[n++] = -fread_d(fp);
      dprintf(2,"\n%d %lg",n,V[n-1]);
      fskipline(fp);
    } else {
      fskipline(fp);
      if(n>0) {
	if(*np>0 && *np != n) error("size of columns mismatch in mat");
	*np = n;
	dprintf(2,"\nn=%d",n);
	A = (double **)renew_mat( A, m + 1, n );
	memcpy( A[m++], V, n * sizeof(double) );
	n=0;
      }
    }
  }
  if(*mp>0 && *mp != m) error("size of rows mismatch in mat");
  *mp = m;
  dprintf(2,"\nm=%d",m);
  free( V );
  return A;
}


#define PAUP2HEAD "Tree\t-lnL\tSite"
double **fread_mat_paup(FILE *fp, int *mp, int *np)
{
  char buff[sizeof PAUP2HEAD];
  int r;
  
  r = fread_line( fp, buff, sizeof PAUP2HEAD );
  fskipline(fp);
  if( r == sizeof PAUP2HEAD
      && strncmp( PAUP2HEAD, buff, sizeof PAUP2HEAD - 1 ) == 0 )
    return fread_mat_paup2(fp, mp, np);
  else
    return fread_mat_paup1(fp, mp, np);
}
