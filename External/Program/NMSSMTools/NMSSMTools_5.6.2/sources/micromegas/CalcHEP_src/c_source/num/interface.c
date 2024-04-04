#include"interface.h"
#include"subproc.h"
#include"rw_sess.h"
#include"kininpt.h"
#include"alphas2.h"
#include"runVegas.h"
#include"dynamic_cs.h"

#include "nType.h"

int nin_int;
int nout_int;
int nprc_int;
int nvar_int;
int nfunc_int;

char * (*pinf_int)(int nsub, int nprtcl, REAL * pmass,int * num);
int (*pinfAux_int)(int nsub, int nprtcl, int*spin2,int*color,int*neutral,int*ndf);
char ** polarized_int;
char ** varName_int;

double (*sqme_int)(int nsub,double GG,REAL*momenta,REAL*cb_coeff,int * err);
int (*calcFunc_int)(void);
double *BWrange_int;
int *twidth_int, *gtwidth_int, *gswidth_int;
REAL *va_int;

colorBasis * cb_int;

#include"../../include/num_out.h"

char * hiddenf=NULL;

void link_process( CalcHEP_interface * interface)
{ int i;
  nin_int =       interface->nin;
  nout_int=       interface->nout;
  nprc_int=       interface->nprc;
  nvar_int=       interface->nvar;
  nfunc_int=      interface->nfunc;

  pinf_int=       interface->pinf;
  pinfAux_int=    interface->pinfAux;
  polarized_int=  interface->polarized;
  varName_int=    interface->varName;

  sqme_int=       interface->sqme;
  calcFunc_int=   interface->calcFunc;
  BWrange_int=    interface->BWrange;
  gtwidth_int=    interface->gtwidth;
  twidth_int=     interface->twidth;
  gswidth_int=    interface->gswidth;
  va_int=         interface->va;

  cb_int=   interface->cb;
  
  *(interface->aWidth)=&aWidth;
  hiddenf=realloc(hiddenf,nfunc_int);
  for(i=0; i<nfunc_int; i++) if(va_int[1+i+nvar_int]>0.5)hiddenf[i]=1; else hiddenf[i]=0;

  Nsub=1;
  wrtprc_();

  stdkin_();
  i_alphaQCD();
  i_Scales();
  nSess = 1;

  integral.n_it=0;

  *BWrange_int=2.7;
  *gswidth_int=0;
  *gtwidth_int=0;
  *twidth_int=0;
}
