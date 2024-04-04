#include <stdio.h>
#include <stdio.h>
#include <math.h> 

#include"nType.h"
        
extern int slhaReadW(char *fname,int mode);
extern int slhaRead(char *fname,int mode);
extern double slhaVal(char * Block, double Q, int nKey, ...);
extern int slhaValExists(char * Block, int nKey, ...);
extern double slhaWidth(int pNum);
extern int slhaWrite(char *fname);
extern int slhaWarnings(FILE*f);
extern int slhaDecayExists(int pNum);
extern double slhaBranch(int pNum,int N, int * nCh);
extern int initDiagonal(void);
extern int rDiagonal(int nDim,...); 
extern int rDiagonalA(int nDim,...);
extern REAL MassArray(int id,  int i);
extern REAL MixMatrix(int id, int i,int j); 
extern REAL MixMatrixU(int id, int i,int j);
extern int cDiagonalH(int Dim,...);
extern int cDiagonalA(int Dim,...);
extern int cDiagonalS(int Dim,...);

extern COMPLEX cMixMatrix(int id,int i,int j);
extern COMPLEX cMixMatrixU(int id,int i,int j);

extern int System(char * format, ...);
extern int openAppend(char * fileName);
extern int aPrintF(char * format,...);

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


extern double initQCD(double,double,double,double);
extern double initQCD5(double,double,double,double);
extern double MbEff(double);
extern double MtEff(double);
extern double McEff(double);
extern double MqEff(double mass2GeV, double Q);
extern double alphaQCD(double);
extern double nfQCD(double Q);
extern double poleQmass(double M_M_, double alpha, int nf);

extern double complex HggF(double z);
extern double complex HggS(double z);
extern double complex HggV(double z);
extern double complex HggA(double z);

extern double complex Hgam1F(double z);
extern double complex Hgam1S(double z);
extern double complex Hgam1V(double z);
extern double complex Hgam1A(double z);
extern double Mbp(void);


extern double dbl(REAL a);
extern double slhaVal0(char *block, double scale);
extern double slhaVal1(char *block, double scale,int i1);
extern double slhaVal2(char *block, double scale,int i1,int i2);
extern double slhaVal3(char *block, double scale,int i1,int i2,int i3);
extern double slhaVal4(char *block, double scale,int i1,int i2,int i3,int i4);

extern int  slhaValExists0(char *block);
extern int  slhaValExists1(char *block,int i1);
extern int  slhaValExists2(char *block,int i1,int i2);
extern int  slhaValExists3(char *block,int i1,int i2,int i3);
extern int  slhaValExists4(char *block,int i1,int i2,int i3,int i4);



extern double complex  hGGeven(double MH, double alphaMH, int Nitems, ...);
extern double complex  hAAeven(double MH, double alphaMH, int Nitems, ...);
extern double complex  hGGodd(double MH,  double alphaMH, int Nitems, ...);
extern double complex  hAAodd(double MH,  double alphaMH, int Nitems, ...);

extern double complex lAAhiggs(double Q, char*hName);
extern double complex lGGhiggs(double Q, char*hName);
extern double complex lAA5higgs(double Q,char*hName);
extern double complex lGG5higgs(double Q,char*hName);


/* To avoid avto-prototyping  

extern double  sqrt(double);
extern double  sin(double);
extern double  cos(double);
extern double  tan(double);
extern double  asin(double);
extern double  acos(double);
extern double  atan(double);
extern double  exp(double);
extern double  log(double); 
extern double  pow(double,double);
extern double  fabs(double);
extern double  atan2(double,double);
extern double  log10(double);
extern double  sinh(double);
extern double  cosh(double);
extern double  tanh(double);
extern double  asinh(double);
extern double  acosh(double);
extern double  atanh(double);

extern double  creal(double complex);
extern double  cimag(double complex);
extern double  carg(double complex);
extern double  cabs(double complex);
extern double complex conj(double complex);
extern double complex cacos(double complex);
extern double complex casin(double complex);
extern double complex catan(double complex);
extern double complex ccos(double complex);
extern double complex csin(double complex);
extern double complex ctan(double complex);
extern double complex cacosh(double complex);
extern double complex casinh(double complex);
extern double complex catanh(double complex);
extern double complex ccosh(double complex);
extern double complex csinh(double complex);
extern double complex ctanh(double complex);
extern double complex cexp(double complex);
extern double complex clog(double complex);
extern double complex clog10(double complex);
extern double complex cpow(double complex,double complex);
extern double complex csqrt(double complex);
extern double complex cproj(double complexl);

extern int printf(char*, ...); 

extern double  if(double,double,double);

*/
