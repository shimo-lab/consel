/*
  opt.h
  Time-stamp: <2002-02-28 16:47:56 shimo>

  $Id$
*/



void luinverse(double **omtrx, double **imtrx, int size); /* From Molphy */

double sym_mat(double **mat, int m);

double *lsfit(double **X, double *Y, double *W,
	      int m, int n,
	      double *beta, double *rss, double ***acmat);

void dfsimple(double parm[], double diff[], double xh[],
	      int n, double(*func)(double []));

void dfmpridr(double parm[], double diff[], double xh[],
	      int n, double(*func)(double []));

void ddfsimple(double parm[], double **diff, double xh[],
	       int n, double(*func)(double []));

void dfnmin(double p[], int n, double gtol, int *iter, double *fret,
	    double ***hesinv,
	    double(*func)(double []), void (*dfunc)(double [], double []),
	    void (*ddfunc)(double [], double **));
