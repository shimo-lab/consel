/* rand.h 
$Id: rand.h,v 1.4 2002/02/20 08:53:47 shimo Exp shimo $
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

/* misc */
double dnorm(double x);
double qnorm(double x);

/* recipe */
double gammln(double xx);
double gammp(double a, double x);
double gammq(double a, double x);
double erfn(double x);
double erfcn(double x);
double erfcn2(double x);
double pnorm(double x);
double pgammadist(double x, double nu);
double tgammadist(double x, double nu);
double pchisq(double x, double df);
double tchisq(double x, double df);
double pchisqnc(double x, double df, double nc);
double tchisqnc(double x, double df, double nc);
