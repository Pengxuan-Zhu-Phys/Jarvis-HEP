#ifndef __NUM_OUT_int
#define __NUM_OUT_int

#include<stdlib.h>
#include<string.h> 
#include<math.h>

#include"nType.h"
#include"num_out.h"

extern  int nin_int;
extern  int nout_int;
extern  int nprc_int;
extern  int nvar_int;
extern  int nfunc_int;

extern char * (*pinf_int)(int nsub, int nprtcl,  REAL* pmass, int*pnum);
extern int (*pinfAux_int)(int nsub, int nprtcl, int* spin2, int*color,int*neutral,int*ndf);
extern char ** polarized_int;
extern char ** varName_int;

extern double (*sqme_int)(int nsub,double GG,REAL*momenta,REAL*cb_coeff,int*err);
extern int (*calcFunc_int)(void);
extern int *twidth_int, *gswidth_int, *gtwidth_int;
extern double *BWrange_int;
extern REAL *va_int;
extern char*hiddenf;

extern colorBasis * cb_int;
extern REAL * cb_coeff_int;


#define  DENOMINATOR_ERROR   2
#endif

