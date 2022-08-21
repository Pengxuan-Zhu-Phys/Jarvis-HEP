#ifndef __MICRO_AUX__
#define __MICRO_AUX__

#define MPlanck      1.22091E19   /* Planck mass in [GeV] */
#define EntropyNow  2.8912E9     /* Present day entropy   [1/m^3] */
#define RhoCrit100  10.537       /* A critical density of Universe
                                    for Hubble rate  100 km/s/Mpc [GeV/m^3] */
                                    
#include"micromegas.h"
#include"../CalcHEP_src/include/VandP.h"



#ifdef __cplusplus
extern "C" {
#endif 

#include "../CalcHEP_src/c_source/dynamicME/include/dynamic_cs.h"
#include "../CalcHEP_src/c_source/dynamicME/include/vp.h"
#include "../CalcHEP_src/c_source/ntools/include/polint.h"
#include "../CalcHEP_src/c_source/ntools/include/vegas.h"
#include "../CalcHEP_src/c_source/ntools/include/1d_integration.h"

extern  char * trim(char *);

//==================== displayPlot =================
#include "../CalcHEP_src/c_source/plot/include/plot.h"

extern void displayPlotN(char * title,  char*xName, double xMin, double xMax, int Scale, int N,char**Y, int*Dim,double**f,double**ff);
extern int displayPlot(char*title,char*xName,double xMin,double xMax,  int Scale, int N, ...);

#define displayFunc(title,Func,xMin,xMax,xName,Scale)  displayPlot(title,xName,xMin,xMax,Scale, 1, title, 0, Func,NULL)

// integration 
extern double simpson(double (*func)(double), double a, double b, double eps,int*err);
extern double simpson_arg( double (*func)(double,void*par), void*par,double a,double b,double  eps,int*err);
extern int simpson_err;
extern double gauss(double (*func)(double), double a, double b, int n);
extern double gauss_arg( double (*func)(double,void*par),void*par,double a,double b, int n);
extern double gauss345_arg( double (*func)(double,void*),void*par,double a,double b, double eps,int * err_code);
extern double peterson21(    double (*F)(double),                double a, double b, double *aerr);
extern double peterson21_arg(double (*F)(double,void*),void*par, double a, double b, double *aerr);



extern int solveLinEq( int N,double *A, double *c);
extern double detAl( int N, double  *A);
/*======================
   Lock/UnLock service
========================*/
extern int checkLockFile(int *delay);
extern void removeLock(int fd);
/*
extern char *  prepareWorkPlace(void);
extern int cleanworkPlace(void);
*/
extern char*WORK;

/*=======================
      Odd particles
========================*/
extern int Nodd;
extern  ModelPrtclsStr * OddPrtcls;
extern int createTableOddPrtcls(void);
extern char * OddParticles(int mode);
extern int  * oddPpos;
extern char * EvenParticles(void); 
extern int maxPlistLen;
    
/*=============================================
C->F and F->C    string convertation     
===============================================*/

extern void  cName2f(char*c_name,char*f_name,int len);
extern void  fName2c(char*f_name,char*c_name,int len);

/*=============================================
Fortran output
===============================================*/
extern void fortreread_(int* N, char * fname, int len);

/*====================================
   Tools for integration 
=====================================*/

extern double vegas_chain(int ndim, double (*Integrand)(double*, double),
                          int N0, double Nfact, double eps,double * dI,void (*clean)(void));


extern int  odeint(double*ystart,int nvar,double x1,double x2, double eps,
                   double h1, void(*derivs)(double,double*,double *));

extern int stiff( int first,double xstart, double xend, int n, double*y, double *yscal, double eps, double*htry,
    void (*derivs)(double, double*, double*, double,double*,double*));
extern int stifbs(int first,double xstart, double xend, int nv, double*y, double *yscal, double eps, double*htry,
    void (*derivs)(double,double*,double*,double,double*,double*));
   
    
/*==== Tool  for interpolation  ====*/
extern int buildInterpolation( double (*Fun)(double), double x1,double x2, 
            double eps, double delt, int * N, double ** xa, double **ya);
            
extern void spline(double*x,double*y,int n,double*y2);            
extern void splint(double*xa, double*ya, double*y2a, int n, double x, double *y);
            
/*======= special functions ========*/
extern double bessI0(double x);
extern double bessK0(double x);
extern double bessK1(double x);
extern double bessK2(double x);
extern double K2pol(double x); /*bessk1(1/x)*exp(1/x)*sqrt(2/M_PI/x);*/
extern double K1pol(double x); /*bessk2(1/x)*exp(1/x)*sqrt(2/M_PI/x);*/

extern double K1to2(double m1,double m2,double m3, double s1,double s2,double s3);
extern double Stat2(double P, double M, double m1,double m2,double eta1, double eta2);

/*======== read Table =============*/

extern int readTable(char * fileName, int *Ncolumn, double **tab);


/* Hidden interface with CalcHEP */
extern int FError;
extern int OnlyTEQ0;

extern int pname2lib(char*pname, char * libname);
extern numout*newProcess_(int twidth, int model,int UG,char*Process,
          char * excludeVirtual,char*excludeOut,char*lib,int usr);


//=====  2->2 processes ===========
extern void  all22procList(void);
extern double (*sqme22)(int nsub, double GG, REAL *pvect, REAL*cb_coeff, int * err_code); 

extern int     kin22(double PcmIn,REAL * pmass);
extern double  kinematic_23(double Pcm,int i3, double M45, double cs1, double cs2,  double fi,REAL*pmass, REAL*pvect);
extern double  kinematic_24(double Pcm,int i3, int i4, double M1, double M2, double xcos,double xcos1, double xcos2, double fi1, double fi2,
                            REAL*pmass, REAL * P);
                            
extern double  dSigma_dCos(double  cos_f);
extern int  nsub22;
extern double  vcs22(numout * cc,int nsub,int * err); 

typedef  struct{ numout*cc; 
                 int nsub,s13,s14;
                 int err;
                 int sqme_err; 
                 double GG,M13,M14,e13,e14,sqrtSmin,totFactor,E;
                 int spin2[4],pdg[4],ndf[4];
                 double eta[4]; 
                 REAL PcmIn, PcmOut, pmass[4], pvect[16];
                 double T,ch,sh,n[3];
               } par_22;
#define ErrE   4
#define ErrC   1
#define ErrD   2
#define ErrN   8
#define ErrIA 16
#define ErrIS 32
#define ErrIT 64
#define ErrP 128
#define Errt 256

extern int init22_par(par_22*arg, numout*cc, int nsub);
extern void mass22_par(par_22*arg,double T);
extern int  kin22_par(par_22*arg, REAL sqrtS,double GG);
extern double sqmeInt(par_22*arg,double eps);

extern void mass22_parDel(par_22*arg,double T);
extern double sqmeIntDel(par_22*arg,double eps);

#define NTOF(X) extern forCalchep1 X; double X(double ok){return findValW(#X);}

typedef  double (forCalchep1)(double);
typedef  double (forCalchep2)(double,double);

/*  Loop integrals I1 ... I5   DreesNojiri */
extern double   LintIk(int II,double MSQ,double MQ,double MNE);

extern int readVarSpecial(char *fname, int nVar, char ** names);

extern double parton_x( int pNum, double  Q);

extern double (*parton_alpha)(double q);
extern double (*parton_distr)(int pNum, double x, double q);

extern double convStrFun2(double x, double q, int pc1, int pc2, int ppFlag ); 
                           /* result of convolution of structure functions 
                              of pc1 and pc2 particles  */

extern int  wimpPos(void);


extern int  vPolar( char**N,  double*lng);
extern double plazmaWidth(char *process,double T);

extern double cs23MM(numout*cc, int nsub, double Pcm, int fast,int ii3, double M45min, double M45max);

extern double cs23(numout*cc, int nsub, double Pcm, int ii3,int*err);

extern double amoeba(double *p, double *dp, int ndim, double (*f)(double *), double eps, int *nCalls);
               
extern REAL *Qaddress;

extern double lGGhSM(double Mh, double aQCDh, double Mcp,double Mbp,double Mtp,double vev);
extern double lAAhSM(double Mh, double aQCDh, double Mcp,double Mbp,double Mtp,double vev);
extern int hbBlocksMO(char*fname,int *nHch);
extern int  LilithMO(char * fname);

extern int makePdtConv(void);
extern int initPDFconv(void);

// Statistics
extern double FeldmanCousins(int n0, double b, double cl);
extern double ch2pval(int nexp, double ch2obs);

extern double plr2pval(double l);
extern double pval2plr(double p); 

extern double uConversion(int u1,int u2);

#define  _GeV_  0
#define  _g_    1
#define  _K_    2   
#define  _cm_1_ 3
#define  _s_1_  4
#define  _MP_   5
#define  _erg_  6

extern double FSRdNdE(double E, double p,double m, double q, int spin2);

//  python
void  pythonversion_(int *n1,int *n2);

//  LHAPDF

extern double tWidth21(char *name, double T,int show);

extern float UpperLim(float CL,int If, int N, float* FC, float muB,float*FB,int *Iflag);

extern void  addErrorMess( char** All, char * one);
extern void delInterval(double x1,double x2,double **intervals, int*n);

#include"../CalcHEP_src/include/num_in.h"

#include"microPath.h"


#ifdef __cplusplus
}
#endif 


#endif
