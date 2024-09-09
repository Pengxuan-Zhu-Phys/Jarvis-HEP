#include<math.h>
#include<stdio.h>
#include<unistd.h>
#include<sys/stat.h>
#include<sys/types.h>

#include <dlfcn.h>
#include <sys/wait.h> 
             

#include"nType.h"
#include"num_in.h"
#include"num_out.h"
#include"VandP.h"
#include"dynamic_cs.h"
#include"rootDir.h"
#include"vegas.h"

static REAL m[4];   // particle masses
static int pdg[4];  // particle codes
static double Pcm;  // collider beam momentum
static double PtCut;//  Pt cut
static numout*cc;   // address of code of matrix element

// LHAPDF  functions   
extern double evolvepdf_(double *x,double *Q,double*pdf);
extern double alphaspdf_(double *Q); 
extern void initpdfsetbyname_(const char* setname, int len);

static double integrand(double *x, double w)
{
    double cs=2*x[0]-1;  // cosine of scattering angle
    double J=2;          // Jacobian 
    REAL sn=sqrt(1-cs*cs);
    REAL pin=Pcm*sqrt(x[1]*x[2]);
    REAL e1=sqrt(m[0]*m[0]+pin*pin),e2=sqrt(m[1]*m[1]+pin*pin);
    REAL s_hat= pow(e1+e2,2);
    REAL  ms=m[2]+m[3],md=m[2]-m[3];
    REAL pout=sqrt((s_hat - ms*ms) * (s_hat - md*md))/(2*sqrt(s_hat));
    REAL pvec[16];
    for(int i=0;i<16;i++) pvec[i]=0;   
// p1
    pvec[0]= e1;
    pvec[3]= pin;
//  p2    
    pvec[4]= e2;
    pvec[7]= -pin;
//  p3    
    pvec[8]=  sqrt(pout*pout + m[2]*m[2]);
    pvec[10]= pout*sn;
    pvec[11]=  pout*cs;
    if(pvec[10]<PtCut) return 0;
// p4    
    pvec[12]=  sqrt(pout*pout + m[3]*m[3]);
    pvec[14]= -pvec[10];
    pvec[15]=  -pvec[11];

    double Scale=sqrt(s_hat);
    Scale=100;
    double GG=sqrt(alphaspdf_(&Scale)*4*M_PI);
    double totcoef=3.8937966E8 * pin /(32.0 * M_PI * pout * s_hat);
    int err=0; 
    double dSigmadCos=totcoef*cc->interface->sqme(1,GG,pvec,NULL,&err);
    double f1=0,f2=0,ff[14];
    
    evolvepdf_(&(x[1]),&Scale,ff);
    switch(pdg[0])
    {
        case 2 :           f1=ff[8];  break;
        case 1 :           f1=ff[7];  break;
        case 3 : case -3 : f1=ff[9];  break;
        case 4 : case -4 : f1=ff[10]; break;
        case 5 : case -5 : f1=ff[11]; break; 
        case 21: case -21: f1=ff[6];  break;
        case -1:           f1=ff[5];  break;
        case -2:           f1=ff[4];  break;                             
    } 

    evolvepdf_(&(x[2]),&Scale,ff);
    switch(pdg[1])
    {
        case 2 :           f2=ff[8];  break;
        case 1 :           f2=ff[7];  break;
        case 3 : case -3 : f2=ff[9];  break;
        case 4 : case -4 : f2=ff[10]; break;
        case 5 : case -5 : f2=ff[11]; break; 
        case 21: case -21: f2=ff[6];  break;
        case -1:           f2=ff[5];  break;
        case -2:           f2=ff[4];  break;                             
    } 

    return f1/x[1]*f2/x[2]*dSigmadCos*J; 
}    


int main(void)
{ int err;

// Specify model for work: Model directory and model number.
  setModel("models" ,2 ); 

// SQME code generation 
   cc= newProcess("u,U->m,M");

// Calculation of public constraints  
  err=calcMainFunc();
  if(err) { printf("Can not calculate constrained parameter %s\n",varNames[err]);return err;}
  
  err=passParameters(cc);
  if(err) { printf("Can not calculate constrained parameter %s\n",cc->interface->varName[err]); return err;}
 
/* find masses and PDG codes */
  for(int i=0;i<4;i++) cc->interface->pinf(1,i+1,m+i,pdg+i); 
  
  Pcm=6500; // Energy of LHC Run 2
  PtCut=200; // Cut for  transverse momentum

  char * pdfName="cteq6l1";
  initpdfsetbyname_(pdfName, strlen(pdfName)); 
  {  double x=0.1,q=100,ff[14];  evolvepdf_( &x,&q,ff);alphaspdf_(&q);  } // to avoid PDF initialisation  in multi-core calculation
// integration by Vegas
  vegasGrid *vegPtr=vegas_init(3/*dimension of integration*/,integrand,50); 
// First Vegas call
  double ti,dti;
  vegas_int(vegPtr, 100000 // number of calls
                  , 1.5   // grid improving parameter
                  ,  4    // number of processors
                  , &ti   // evaluated integral
                  , &dti  // statistical uncertainty
                   );
printf("Cross section  %.2E +/- %.1E\n",ti,dti);
// Vegas call with improved grid
  vegas_int(vegPtr, 100000, 1.5,4, &ti , &dti);
printf("Cross section  %.2E +/- %.1E\n",ti,dti);                                                                                           


  return 0;
}
