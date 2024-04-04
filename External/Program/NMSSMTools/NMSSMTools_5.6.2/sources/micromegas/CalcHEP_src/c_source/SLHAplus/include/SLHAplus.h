#ifndef __SLHA_PLUS_
#define __SLHA_PLUS_

#include"../../../include/nType.h"

#include<stdio.h>
#include<stdarg.h>
#include<stdlib.h>
#include<sys/wait.h>
#include<unistd.h>
#include<string.h>
#include<ctype.h>
#include<sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif 

#include "aList.h"

extern int FError;
extern int slhaBasicReader( int mode, int (*readLn)(int, char*),int * anydate,char * end );
extern int slhaReadStream(FILE*f,int mode, char*end);
extern void cleanSLHAdata(void);
extern int slhaRead(char *fname,int mode);

extern double  slhaVal(char * Block, double Q, int nKey, ...);
extern double  slhaValFormat(char * Block, double Q, char * format);
extern int   slhaSTRFormat(char * Block, char * format, char *txt);
#ifndef __cplusplus
extern double complex cslhaVal(char * Block, double Q, int nKey, ...);
#endif
extern int slhaValExists(char * Block, int nKey, ...);
extern double slhaWidth(int pNum);
extern int slhaWrite(char *fname);
extern int slhaWarnings(FILE*f);
extern char * Warnings;
extern int slhaDecayExists(int pNum);
extern double slhaBranch(int pNum,int N, int * nCh);
extern double slhaBr(int pNum, int len, ...);
extern char* slhaComment;
extern int findQnumbers(int pdg, int *eQ3,int * spinDim,int*cDim,int *neutral);
extern int allQnumbers(int i, int *pdg,int*eQ3,int*spinDim,int*cDim,int*neutral);
#ifndef __cplusplus
extern int allBlocks(int i,int j,char*name,int*Len,int*key, double complex * val);
#endif
extern int allDecays(int i,int j,int* pdg, int*Len,int*decay,double*width,double*br);

extern int rJacobi(REAL* a, int n, REAL *d, REAL * v);
extern int rJacobiA(REAL*  a, int n, REAL *d, REAL* u,REAL * v);
#ifndef __cplusplus
extern int cJacobiH(COMPLEX* a, int n, REAL *d, COMPLEX* v);
extern int cJacobiS(COMPLEX* a, int n, REAL *d, COMPLEX* v);
extern int cJacobiA(COMPLEX* a, int n, REAL *d, COMPLEX* u,COMPLEX * v);
#endif
extern int initDiagonal(void);
extern int rDiagonal(int nDim,...);
extern int rDiagonalA(int nDim,...);
extern REAL MassArray(int id,  int i);
extern REAL MixMatrix(int id, int i,int j);
extern REAL MixMatrixU(int id, int i,int j);
extern int cDiagonalH(int Dim,...);
extern int cDiagonalA(int Dim,...);
extern int cDiagonalS(int Dim,...);
#ifndef __cplusplus
extern COMPLEX cMixMatrix(int id,int i,int j);
extern COMPLEX cMixMatrixU(int id,int i,int j);
#endif

extern int rDiagonal2(aList3(double));
extern int rDiagonal3(aList6(double));
extern int rDiagonal4(aList10(double));
extern int rDiagonal5(aList15(double));

extern int rDiagonalA2(aList4(double));
extern int rDiagonalA3(aList9(double));
extern int rDiagonalA4(aList16(double));
extern int rDiagonalA5(aList25(double));

#ifndef __cplusplus
extern int  cDiagonalH2(aList3(double complex));
extern int  cDiagonalH3(aList6(double complex)); 
extern int  cDiagonalH4(aList10(double complex));
extern int  cDiagonalH5(aList15(double complex));

extern int  cDiagonalS2(aList3(double complex));
extern int  cDiagonalS3(aList6(double complex)); 
extern int  cDiagonalS4(aList10(double complex));
extern int  cDiagonalS5(aList15(double complex));

extern int  cDiagonalA2(aList4(double complex)); 
extern int  cDiagonalA3(aList9(double complex));
extern int  cDiagonalA4(aList16(double complex));
extern int  cDiagonalA5(aList25(double complex));

#endif

extern unsigned sysTimeLim;
extern unsigned sysTimeQuant;

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
extern int aPrintF6(char * format, double x1,double x2,double x3,double x4,double x5,double x6);
extern int aPrintF7(char * format, double x1,double x2,double x3,double x4,double x5,double x6,double x7);
extern int aPrintF8(char * format, double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8);
extern int aPrintF9(char * format, double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8,double x9);


extern double initQCD(double MZalphaS,double McMc,double MbP,double MtP);
extern double initQCD5(double MZalphaS,double McMc,double MbMb,double MtP);

extern double alphaQCD(double Q);
extern double MbRun(double Q);
extern double MbEff(double Q);
extern double MtRun(double Q);
extern double MtEff(double Q);
extern double McRun(double Q);
extern double McEff(double Q);
extern double MbPole;
extern double poleQmass(double M_M_, double alpha, int nf);

extern double bPoleMass(void);
extern double tPoleMass(void);
    
extern double MqRun(double mass2GeV, double Q);
extern double MqEff(double mass2GeV, double Q);
extern double nfQCD(double Q);

#ifndef __cplusplus
extern double complex HggF(double z);
extern double complex HggS(double z);
extern double complex HggV(double z);
extern double complex HggA(double z);

extern double complex Hgam1F(double z);
extern double complex Hgam1S(double z);
extern double complex Hgam1A(double z);

extern double complex hGGeven(double MH, double alphaMH, int Nitems, ...);
extern double complex hAAeven(double MH, double alphaMH, int Nitems, ...);
extern double complex  hGGodd(double MH, double alphaMH, int Nitems, ...);
extern double complex  hAAodd(double MH, double alphaMH, int Nitems, ...);
#endif

extern int blockExists( char * block);
extern int writeBlock(char*block,  FILE*f);

#include "delList.h"

#ifdef __cplusplus
}
#endif 


#endif
