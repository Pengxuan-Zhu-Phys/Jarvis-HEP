#include<stdlib.h>

#ifndef __Variables_and_Particles__
#define __Variables_and_Particles__

#include "nType.h"

extern int nModelParticles;
 
typedef struct
{ 
   char* name; char* aname;  int selfC; int NPDG; char* mass; char* width; int spin2; int cdim;  int g;  int q3;
}  ModelPrtclsStr;

extern ModelPrtclsStr*ModelPrtcls;
extern int nModelVars;
extern int nModelFunc;
extern int*currentVarPtr;
extern char**varNames;
extern REAL *varValues;
extern int calcMainFunc(void);

/*  VandPgate */

extern double usrFF_(int n_in, int n_out,double * pvect,char**pnames,int*pdg);
extern double usrfun_(char * name,int n_in, int n_out, double * pvect,char**pnames,int*pdg);
extern double usrFF(int n_in, int n_out,double * pvect,char**pnames,int*pdg);
extern double usrfun(char * name,int n_in, int n_out, double * pvect,char**pnames,int*pdg);

extern double alpha_lha(double q );
extern void   getdatapath(char* dirpath, int len);
extern void   getlhapdfversion(char* s, size_t len);
extern void   initpdfsetbynamem(int *P,char *name, int len);
extern void   evolvePDFm(int P,double x,double Q,double *f);
extern int    has_photon(void);
extern void   initpdfm(int* P,int * nSet,double*xMin,double*xMax,double*qMin,double*qMax );
extern void   numberpdfm(int* P,int * nMax);
extern double alphaspdfm(int*S,double*Q);
#endif
