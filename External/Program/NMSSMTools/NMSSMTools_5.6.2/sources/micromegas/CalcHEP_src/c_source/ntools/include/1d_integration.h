#ifndef  __SIMPSON__
#define  __SIMPSON__
extern double gauss( double (*func)(double),double a,double b, int n);
extern double gauss345(double (*func)(double), double a, double b, double eps,int * err_code);
extern double gauss345_arg( double (*func)(double,void*),void*par,double a,double b, double eps,int * err_code);
extern double gauss_arg( double (*func)(double,void*par),void*par,double a,double b, int n);

extern int simpson_err;
extern double simpson( double (*func)(double),double a,double b, double  eps, int *err);
extern double simpson_arg( double (*func)(double,void*par), void*par,double a,double b,double eps, int*err);

extern double peterson21(double (*func)(double), double a, double b, double *aerr);
extern double peterson21_arg(double (*F)(double,void*),void*par, double a, double b, double *aerr);
#endif
