#ifndef __NUM_IN_
#define __NUM_IN_

#include<stdlib.h>
#include<string.h> 
#include<math.h>
#include "nType.h"

extern double alpha_2(double);
#ifndef __cplusplus
typedef REAL (DNN)(double, REAL *,REAL*,double*, int *);
typedef REAL (FNN)(double, REAL*,REAL*,REAL*,COMPLEX*,REAL*,REAL*,int gsw,int gtw);
#endif

extern  REAL Helicity[2];
extern  REAL N_pol_11_,N_pol_12_,N_pol_21_,N_pol_22_;

extern  int    CalcConst;
extern  int    indx_(int k,int l);
extern  void   sprod_(int ntot, REAL * momenta, REAL*DP);
#ifndef __cplusplus
extern  int    prepDen(int nden, int nin, double BWrange2,   
                       REAL * dmass,REAL * dwidth, char **q,REAL * mom, 
                       REAL *Q0,COMPLEX*Q1,REAL*Q2);
#endif
#endif
