#include<string.h>
#include <pthread.h>
#include"VandP.h"


double usrFF_(int n_in, int n_out,double * pvect,char**pnames,int*pdg)
{ return  usrFF(n_in,n_out,pvect,pnames,pdg);}
double usrfun_(char * name,int n_in, int n_out,double * pvect,char**pnames,int*pdg) 
{ return  usrfun(name,n_in,n_out,pvect,pnames,pdg);}


extern void   getdatapath_(char* dirpath, int len);
extern void  initpdfsetbynamem_(int *P,char *name, int len);
extern void  getlhapdfversion_(char* s, size_t len);
extern void  evolvepdfphotonm_(int* P,double *x,double *Q,double *f,double *fph);
extern void  initpdfm_(int* P,int * nSet );
extern void  numberpdfm_(int* P,int * nMax);
extern double   alphaspdfm_(int*,double*);
extern void getxmaxm_(int*,int*,double *);
extern void getxminm_(int*,int*,double *);
extern void getq2maxm_(int*,int*,double *);                     
extern void getq2minm_(int*,int*,double *);
extern int has_photon_(void);
extern void  evolvepdfm_(int*P,double* x, double*q, double*fxq);

double alphaspdfm(int*S,double*Q){ return alphaspdfm_(S,Q);}

static int photon[5];
void   getdatapath(char* dirpath, int len){ getdatapath_(dirpath, len);}
void   initpdfsetbynamem(int *P,char *name, int len)
       { initpdfsetbynamem_(P,name,len);
         photon[*P]=has_photon_();
       
       }
void   getlhapdfversion(char* s, size_t len){ getlhapdfversion_(s,len);}
void   evolvePDFm(int P,double x,double Q,double *f)
       {  if(photon[P]) { double fph; evolvepdfphotonm_(&P,&x,&Q,f,&fph); f[13]=fph;}
          else { evolvepdfm_(&P,&x,&Q,f); f[13]=0; }  
       }

void   initpdfm(int* P,int * nSet,double*xMin,double*xMax,double*qMin,double*qMax)
{
  initpdfm_(P,nSet ) ;
  getxmaxm_(P,nSet,xMax); 
  getxminm_(P,nSet,xMin);
  getq2maxm_(P,nSet,qMax);  *qMax=sqrt(fabs(*qMax));
  getq2minm_(P,nSet,qMin);  *qMin=sqrt(fabs(*qMin));
  double x=0.1, q=100, f[15],fph;
  evolvepdfphotonm_(P,&x,&q,f,&fph);
  alphaspdfm_(P,&q);
//  printf("evolvepdfphotonm_ ok\n");
//  printf("limits: %E %E %E %E\n", *xMin,*xMax,*qMin,*qMax);
}
int has_photon(void) { return has_photon_();}
void  numberpdfm(int* P,int * nMax){ numberpdfm_(P,nMax);}

extern int findval(char *name,double *val);
extern int qnumbers(char*pname, int *spin2, int * charge3, int * cdim);

int findval(char *name,double *val)
{ int i;
  for(i=0;i<nModelVars+nModelFunc;i++) 
  if(!strcmp(name,varNames[i])){ *val= varValues[i]; return 0;}  
  return 1;
}


int qnumbers(char*pname, int *spin2, int * charge3, int * cdim)
{
   int n,sign;
   for(n=0;n<nModelParticles;n++)
   { 
     if(!strcmp(pname,ModelPrtcls[n].name )) {sign=1; break;} 
     if(!strcmp(pname,ModelPrtcls[n].aname)) {sign=-1;break;}
   }
   if(n==nModelParticles) return 0;

   if(spin2)   *spin2  =ModelPrtcls[n].spin2;
   if(charge3) *charge3=sign*ModelPrtcls[n].q3;
   if(cdim)    
   {  *cdim   =ModelPrtcls[n].cdim; 
      if(sign==-1 &&(*cdim==3 || *cdim==-3)) (*cdim)*=-1;
   }
   return sign*ModelPrtcls[n-1].NPDG;
}
