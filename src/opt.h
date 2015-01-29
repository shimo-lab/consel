/*
  opt.h
  Time-stamp: <2002-03-01 13:32:40 shimo>

  $Id: opt.h,v 1.2 2002/03/01 09:38:22 shimo Exp $
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

void dfnmin(double p[], int n, double gtol, int itmax, int maxback,
	    int *iter, double *fret,
	    double ***hesinv,
	    double(*func)(double []), void (*dfunc)(double [], double []),
	    void (*ddfunc)(double [], double **));
