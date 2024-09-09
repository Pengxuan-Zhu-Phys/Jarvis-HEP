#ifndef __SLHA_PLUS_
#define __SLHA_PLUS_

#include<stdio.h>
#include<math.h>
#include<stdarg.h>
#include<stdlib.h>
#include<sys/wait.h>
#include<unistd.h>
#include<string.h>
#include<ctype.h>


#ifdef __cplusplus

#include <cmath>
#include <complex>
using std::complex;

typedef complex<double> Complex;

extern "C" {
#else 

#include<complex.h>
#define Complex double complex

#endif 

#include "aList.h"

extern int FError;


extern int slhaBasicReader( int mode, int (*readLn)(int, char*),int * anydate,char * end );

extern int slhaRead(char *fname,int mode);
int slhaReadStream(FILE*f, int mode, char * end );

extern double  slhaVal(char * Block, double Q, int nKey, ...);
extern Complex cslhaVal(char * Block, double Q, int nKey, ...);
extern int slhaValExists(char * Block, int nKey, ...);
extern double slhaWidth(int pNum);
extern int slhaWrite(char *fname);
extern int slhaWarnings(FILE*f);
extern int slhaDecayExists(int pNum);
extern double slhaBranch(int pNum,int N, int * nCh);

extern char* slhaComment;
extern int findQnumbers(int pdg, int *eQ3,int * spinDim,int*cDim,int *neutral);
extern int allQnumbers(int i, int *pdg,int*eQ3,int*spinDim,int*cDim,int*neutral);
extern int allBlocks(int i,int j,char*name,int*Len,int*key, Complex * val);
extern int allDecays(int i,int j,int* pdg, int*Len,int*decay,double*width,double*br);

extern int rJacobi(double*  a, int n, double *d, double * v);
extern int rJacobiA(double*  a, int n, double *d, double* u,double * v);
extern int cJacobiH(Complex* a, int n, double *d, Complex* v);
extern int cJacobiS(Complex* a, int n, double *d, Complex* v);
extern int cJacobiA(Complex* a, int n, double *d, Complex* u,Complex * v);

extern int initDiagonal(void);
extern double rDiagonal(double nDim,...);
extern double rDiagonalA(double nDim,...);
extern double MassArray(int id,  int i);
extern double MixMatrix(int id, int i,int j);
extern double MixMatrixU(int id, int i,int j);
extern int cDiagonalH(int Dim,...);
extern int cDiagonalA(int Dim,...);
extern int cDiagonalS(int Dim,...);
extern Complex cMixMatrix(int id,int i,int j);
extern Complex cMixMatrixU(int id,int i,int j);


extern int rDiagonal2(aList3(double));
extern int rDiagonal4(aList10(double));
extern int rDiagonal5(aList15(double));

extern int rDiagonalA2(aList4(double));
extern int rDiagonalA4(aList16(double));
extern int rDiagonalA5(aList25(double));

extern int  cDiagonalH2(aList3(Complex)); 
extern int  cDiagonalH4(aList10(Complex));
extern int  cDiagonalH5(aList15(Complex));

extern int  cDiagonalS2(aList3(Complex)); 
extern int  cDiagonalS4(aList10(Complex));
extern int  cDiagonalS5(aList15(Complex));

extern int  cDiagonalA2(aList4(Complex)); 
extern int  cDiagonalA4(aList16(Complex));
extern int  cDiagonalA5(aList25(Complex));

extern int System(char * format, ...);
extern int openAppend(char * fileName);
extern int aPrintF(char * format,...);

extern int System1(char * format);
extern int System2(char * format,char*path);

extern int aPrintF0(char * format);
extern int aPrintF1(char*format,double x1);
extern int aPrintF2(char*format,double x1,double x2);
extern int aPrintF3(char * format, double x1,double x2,double x3);
extern int aPrintF4(char * format, double x1,double x2,double x3,double x4);
extern int aPrintF5(char * format, double x1,double x2,double x3,double x4,double x5);

extern double initQCD(double MZalphaS,double McMc,double MbP,double MtP);
extern double alphaQCD(double Q);
extern double MbRun(double Q);
extern double MbEff(double Q);
extern double MtRun(double Q);
extern double MtEff(double Q);
extern double McRun(double Q);
extern double McEff(double Q);
extern double MbPole;

#include "delList.h"

#ifdef __cplusplus
}
#endif 


#endif
