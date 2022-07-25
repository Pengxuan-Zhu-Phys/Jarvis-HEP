#ifndef __POLINT__
#define __POLINT__

extern double polint3(double x, int dim, const  double *xa, const double *ya);
extern double polint1(double x, int n,  double *xa, double *ya);
extern double polint1Exp(double x, int n,  double *xa, double *ya);
extern double  polintN(double x, int n, const  double *xa, const double *ya);
extern double polintDiff(int n, const  double *xa, const double *ya, double * dxdy);
#endif
