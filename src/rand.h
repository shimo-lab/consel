/* rand.h 
$Id$
*/
double mrand(void);      /* returns random # in [0,1) */
void smrand(int seed),   /* seeds mrand() */
  mrandlist(double *x, int len);  /* generates len random numbers */

#define runif() mrand() /* defines uniform dist on [0,1] */
double rnorm(void);  /* normal dist */
double rchisq(double df); /* chi-square */

/* perlman */
double poz(double z);
double critz(double p);
double pochisq(double x, int df);
double critchi(double p, int df);
double pof(double F, int df1, int df2);
double critf(double p, int df1, int df2);