/* rand.c 


  (I) basic random number generation functions (Meschach library)
  double mrand(void);      generate a random number 
  void smrand(int seed)    set seed 
  mrandlist(double *x, int len);   vector

  (II) density types
  runif() : unit random number U[0,1]
  rnorm() : standard normal N(0,1)
  rchisq(df) : chi-square of df degrees of freedom

 $Id: rand.c,v 1.1 1998/11/30 09:03:13 shimo Exp shimo $

 */

#include <math.h>
#include <limits.h>
#include "rand.h"


/* Pseudo random number generator data structures */
/* Knuth's lagged Fibonacci-based generator: See "Seminumerical Algorithms:
   The Art of Computer Programming" sections 3.2-3.3 */

#ifdef LONG_MAX
#define MODULUS LONG_MAX
#else
#define MODULUS 1000000000L     /* assuming long's at least 32 bits long */
#endif

static long mrand_list[56];
static int  started = 0;
static int  inext = 0, inextp = 31;

double mrand(void)
{
    long        lval;
    static double  factor = 1.0/((double)MODULUS);
    
    if ( ! started )
        smrand(3127);
    
    inext = (inext >= 54) ? 0 : inext+1;
    inextp = (inextp >= 54) ? 0 : inextp+1;

    lval = mrand_list[inext]-mrand_list[inextp];
    if ( lval < 0L )
        lval += MODULUS;
    mrand_list[inext] = lval;
    
    return (double)lval*factor;
}


/* mrandlist -- fills the array a[] with len random numbers */
void    mrandlist(a, len)
double    a[];
int     len;
{
    int         i;
    long        lval;
    static double  factor = 1.0/((double)MODULUS);
    
    if ( ! started )
        smrand(3127);
    
    for ( i = 0; i < len; i++ )
    {
        inext = (inext >= 54) ? 0 : inext+1;
        inextp = (inextp >= 54) ? 0 : inextp+1;
        
        lval = mrand_list[inext]-mrand_list[inextp];
        if ( lval < 0L )
            lval += MODULUS;
        mrand_list[inext] = lval;
        
        a[i] = (double)lval*factor;
    }
}


/* smrand -- set seed for mrand() */
void smrand(seed)
int     seed;
{
    int         i;

    mrand_list[0] = (123413*seed) % MODULUS;
    for ( i = 1; i < 55; i++ )
        mrand_list[i] = (123413*mrand_list[i-1]) % MODULUS;

    started = 1;

    /* run mrand() through the list sufficient times to
       thoroughly randomise the array */
    for ( i = 0; i < 55*55; i++ )
        mrand();
}


/* rnorm -- generating normal random variable N(0,1) using
            the "RATIO-OF-UNIFORMS" method propsed by 
            Kinderman and Monahan (1977)
 REF: "Non-Uniform Random Variate Generation" Luc Devroye (1986, p.199)
 PROGRAMMER: H. Shimodaira <shimo@ism.ac.jp> 1997-01
 COPYRIGHT: NONE
*/
double rnorm(void)
{
    static int init=0;
    static double a2;
    double u,v,x,x2,u2;

    if(!init) {
        a2=sqrt(8.0/exp(1.0)); /* 2 a */
        init=1;
    }

    while(1) {
        u=runif(); if(u==0.0) continue;
        v=(runif()-0.5)*a2;
        x=v/u;
        x2=x*x;
        u2=u*u;
        if(x2 <= 6.0-8.0*u+2.0*u2) break;
        if(x2 >= 2.0*(1.0-u2)/u) continue;
        if(u2 <= exp(-x2/2.0)) break;
    }

    return x;
}

static double rgammaXG(double a)
/*
VALID FOR a (shape parameter) >= 1
Luc Devroye, Uon-Uniform Random Variate Generation 1986, p.413 
NOTE: the algorithm in p.413 is wrong, though, the method described
      in the previous pages is correct.
*/
{
  static double a0=-1,b,c,d;
  static double con1=1.38629436111989061883; /* log(4) */
  static double con2=2.50407739677627407337; /* 1+log(9/2) */
  double u,v,x,y,z,r;

  if(a != a0) {
    a0 = a;
    b = a - con1;
    u = sqrt(2.0*a - 1.0);
    c = a + u;
    d = 1.0/u;
  }

  while(1) {
    u = runif();  v = runif();
    y = d*log(v/(1.0-v));
    x = a*exp(y); 
    z = u*v*v;
    r = b + c*y - x;

    if(r >= 4.5*z - con2) break;
    if(r >= log(z)) break;
  }

  return x;
}

double rchisq(double df)
{
  double r;
  if(df < 2) {
    /* use the rgamma for shape <= 1, or */
    /* assume df = 1 (since integer) */
    r=rnorm();
    return r*r;
  } else {
    return rgammaXG(0.5*df) * 2.0;
  }
}


/*

PERLMAN STARTS HERE

 */

/*HEADER
	Module:       z.c
	Purpose:      compute approximations to normal z distribution probabilities
	Programmer:   Gary Perlman
	Organization: Wang Institute, Tyngsboro, MA 01879
	Tester:       compile with -DZTEST to include main program
	Copyright:    none
	Tabstops:     4
*/

#define	Z_EPSILON      0.000001       /* accuracy of critz approximation */
#define	Z_MAX          6.0            /* maximum meaningful z value */

/*FUNCTION poz: probability of normal z value */
/*ALGORITHM
	Adapted from a polynomial approximation in:
		Ibbetson D, Algorithm 209
		Collected Algorithms of the CACM 1963 p. 616
	Note:
		This routine has six digit accuracy, so it is only useful for absolute
		z values < 6.  For z values >= to 6.0, poz() returns 0.0.
*/
double            /*VAR returns cumulative probability from -oo to z */
poz (z)
double	z;        /*VAR normal z value */
	{
	double	y, x, w;
	
	if (z == 0.0)
		x = 0.0;
	else
		{
		y = 0.5 * fabs (z);
		if (y >= (Z_MAX * 0.5))
			x = 1.0;
		else if (y < 1.0)
			{
			w = y*y;
			x = ((((((((0.000124818987 * w
				-0.001075204047) * w +0.005198775019) * w
				-0.019198292004) * w +0.059054035642) * w
				-0.151968751364) * w +0.319152932694) * w
				-0.531923007300) * w +0.797884560593) * y * 2.0;
			}
		else
			{
			y -= 2.0;
			x = (((((((((((((-0.000045255659 * y
				+0.000152529290) * y -0.000019538132) * y
				-0.000676904986) * y +0.001390604284) * y
				-0.000794620820) * y -0.002034254874) * y
				+0.006549791214) * y -0.010557625006) * y
				+0.011630447319) * y -0.009279453341) * y
				+0.005353579108) * y -0.002141268741) * y
				+0.000535310849) * y +0.999936657524;
			}
		}
	return (z > 0.0 ? ((x + 1.0) * 0.5) : ((1.0 - x) * 0.5));
	}

/*FUNCTION critz: compute critical z value to produce given probability */
/*ALGORITHM
	Begin with upper and lower limits for z values (maxz and minz)
	set to extremes.  Choose a z value (zval) between the extremes.
	Compute the probability of the z value.  Set minz or maxz, based
	on whether the probability is less than or greater than the
	desired p.  Continue adjusting the extremes until they are
	within Z_EPSILON of each other.
*/
double        /*VAR returns z such that fabs (poz(p) - z) <= .000001 */
critz (p)
double	p;    /*VAR critical probability level */
	{
	double	minz = -Z_MAX;    /* minimum of range of z */
	double	maxz = Z_MAX;     /* maximum of range of z */
	double	zval = 0.0;       /* computed/returned z value */
	double	poz (), pval;     /* prob (z) function, pval := poz (zval) */
	
	if (p <= 0.0 || p >= 1.0)
		return (0.0);
	
	while (maxz - minz > Z_EPSILON)
		{
		pval = poz (zval);
		if (pval > p)
			maxz = zval;
		else
			minz = zval;
		zval = (maxz + minz) * 0.5;
		}
	return (zval);
	}


/*
	Module:       chisq.c
	Purpose:      compute approximations to chisquare distribution probabilities
	Contents:     pochisq(), critchi()
	Uses:         poz() in z.c (Algorithm 209)
	Programmer:   Gary Perlman
	Organization: Wang Institute, Tyngsboro, MA 01879
	Tester:       compile with -DCHISQTEST to include main program
	Copyright:    none
	Tabstops:     4
*/

/*LINTLIBRARY*/
#define	CHI_EPSILON     0.000001    /* accuracy of critchi approximation */
#define	CHI_MAX     99999.0         /* maximum chi square value */

#define	LOG_SQRT_PI     0.5723649429247000870717135 /* log (sqrt (pi)) */
#define	I_SQRT_PI       0.5641895835477562869480795 /* 1 / sqrt (pi) */
#define	BIGX           20.0         /* max value to represent exp (x) */
#define	ex(x)             (((x) < -BIGX) ? 0.0 : exp (x))

/*FUNCTION pochisq: probability of chi sqaure value */
/*ALGORITHM Compute probability of chi square value.
	Adapted from:
		Hill, I. D. and Pike, M. C.  Algorithm 299
		Collected Algorithms for the CACM 1967 p. 243
	Updated for rounding errors based on remark in
		ACM TOMS June 1985, page 185
*/
double
pochisq (x, df)
double	x;       /* obtained chi-square value */
int 	df;      /* degrees of freedom */
	{
	double	a, y, s;
	double	e, c, z;
	double	poz ();   /* computes probability of normal z score */
	int 	even;     /* true if df is an even number */
	
	if (x <= 0.0 || df < 1)
		return (1.0);
	
	a = 0.5 * x;
	even = (2*(df/2)) == df;
	if (df > 1)
		y = ex (-a);
	s = (even ? y : (2.0 * poz (-sqrt (x))));
	if (df > 2)
		{
		x = 0.5 * (df - 1.0);
		z = (even ? 1.0 : 0.5);
		if (a > BIGX)
			{
			e = (even ? 0.0 : LOG_SQRT_PI);
			c = log (a);
			while (z <= x)
				{
				e = log (z) + e;
				s += ex (c*z-a-e);
				z += 1.0;
				}
			return (s);
			}
		else
			{
			e = (even ? 1.0 : (I_SQRT_PI / sqrt (a)));
			c = 0.0;
			while (z <= x)
				{
				e = e * (a / z);
				c = c + e;
				z += 1.0;
				}
			return (c * y + s);
			}
		}
	else
		return (s);
	}

/*FUNCTION critchi: compute critical chi square value to produce given p */
double
critchi (p, df)
double	p;
int 	df;
	{
	double	minchisq = 0.0;
	double	maxchisq = CHI_MAX;
	double	chisqval;
	
	if (p <= 0.0)
		return (maxchisq);
	else if (p >= 1.0)
		return (0.0);
	
	chisqval = df / sqrt (p);    /* fair first value */
	while (maxchisq - minchisq > CHI_EPSILON)
		{
		if (pochisq (chisqval, df) < p)
			maxchisq = chisqval;
		else
			minchisq = chisqval;
		chisqval = (maxchisq + minchisq) * 0.5;
		}
	return (chisqval);
	}


/*
	Module:       f.c
	Purpose:      compute approximations to F distribution probabilities
	Contents:     pof(), critf()
	Programmer:   Gary Perlman
	Organization: Wang Institute, Tyngsboro, MA 01879
	Tester:       compile with -DFTEST to include main program
	Copyright:    none
	Tabstops:     4
*/

/*LINTLIBRARY*/

#ifndef	I_PI        /* 1 / pi */
#define	I_PI        0.3183098861837906715377675
#endif
#define	F_EPSILON     0.000001       /* accuracy of critf approximation */
#define	F_MAX      9999.0            /* maximum F ratio */

/*FUNCTION pof: probability of F */
/*ALGORITHM Compute probability of F ratio.
	Adapted from Collected Algorithms of the CACM
	Algorithm 322
	Egon Dorrer
*/
double
pof (F, df1, df2)
double	F;
int 	df1, df2;
	{
	int	i, j;
	int	a, b;
	double	w, y, z, d, p;
	
	if (F < F_EPSILON || df1 < 1 || df2 < 1)
		return (1.0);
	a = df1%2 ? 1 : 2;
	b = df2%2 ? 1 : 2;
	w = (F * df1) / df2;
	z = 1.0 / (1.0 + w);
	if (a == 1)
		if (b == 1)
			{
			p = sqrt (w);
			y = I_PI; /* 1 / 3.14159 */
			d = y * z / p;
			p = 2.0 * y * atan (p);
			}
		else
			{
			p = sqrt (w * z);
			d = 0.5 * p * z / w;
			}
	else if (b == 1)
		{
		p = sqrt (z);
		d = 0.5 * z * p;
		p = 1.0 - p;
		}
	else
		{
		d = z * z;
		p = w * z;
		}
	y = 2.0 * w / z;

	for (j = b + 2; j <= df2; j += 2)
		{
		d *= (1.0 + a / (j - 2.0)) * z;
		p = (a == 1 ? p + d * y / (j - 1.0) : (p + w) * z);
		}

	y = w * z;
	z = 2.0 / z;
	b = df2 - 2;
	for (i = a + 2; i <= df1; i += 2)
		{
		j = i + b;
		d *= y * j / (i - 2.0);
		p -= z * d / j;
		}
	/* correction for approximation errors suggested in certification */
	if (p < 0.0)
		p = 0.0;
	else if (p > 1.0)
		p = 1.0;
	return (1.0-p);
	}

/*FUNCTION critf: compute critical F value t produce given probability */
/*ALGORITHM
	Begin with upper and lower limits for F values (maxf and minf)
	set to extremes.  Choose an f value (fval) between the extremes.
	Compute the probability of the f value.  Set minf or maxf, based
	on whether the probability is less than or greater than the
	desired p.  Continue adjusting the extremes until they are
	within F_EPSILON of each other.
*/
double
critf (p, df1, df2)
double	p;
int 	df1;
int 	df2;
	{
	double	fval;
	double	fabs ();          /* floating absolute value */
	double	maxf = F_MAX;     /* maximum F ratio */
	double	minf = 0.0;       /* minimum F ratio */
	
	if (p <= 0.0 || p >= 1.0)
		return (0.0);
	
	fval = 1.0 / p;             /* the smaller the p, the larger the F */
	
	while (fabs (maxf - minf) > F_EPSILON)
		{
		if (pof (fval, df1, df2) < p)     /* F too large */
			maxf = fval;
		else                              /* F too small */
			minf = fval;
		fval = (maxf + minf) * 0.5;
		}
	
	return (fval);
	}


/**** additional functions ****/
#define I_SQRT_2_PI 0.398942280401432677939946059934 /* 1/sqrt(2pi) */
double dnorm(double x) {
  return  I_SQRT_2_PI*exp(-0.5*x*x);
}
