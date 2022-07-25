#include"micromegas.h"
#include"micromegas_aux.h"
#include"micromegas_f.h"


double hprofilezhao_(double *r) {return hProfileZhao(*r);}
double noclamps_(double *r) {return 1;}


void setprofilezhao_(double *alpha, double*beta, double*gamma, double*rc)
{ setProfileZhao(*alpha, * beta, *gamma, *rc );} 

double hprofileeinasto(double * r){ return hProfileEinasto(*r);}
void  setprofileeinasto_(double *alpha){ setProfileEinasto(*alpha);}

static double (*dFdF1)(double *);
static double dFdC1(double x) {return (*dFdF1)(&x);}

static double (*dFdF2)(double *);
static double dFdC2(double x) {return (*dFdF2)(&x);}

void sethaloprofile_(double(*hProfile)(double*))
{
  dFdF1=hProfile;
  setHaloProfile(dFdC1); 
}

void setrhoclumps_(double(*cProfile)(double*))
{
  dFdF2=cProfile;
  setRhoClumps(dFdC2); 
}

void setclumpconst_(double*f,double*rho) {  setClumpConst(*f,*rho); }


double calcspectrum_(int*key, double *Sg, double *Se, double *Sp, double *Sne, double *Snm,double *Snl , int *err,int len)
{  
  double *Sg_=NULL, *Se_=NULL, *Sp_=NULL, *Sne_=NULL, *Snm_=NULL, *Snl_=NULL;
  
  if(Sg !=mocommon_.par) Sg_=Sg;
  if(Se !=mocommon_.par) Se_=Se;
  if(Sp !=mocommon_.par) Sp_=Sp;
  if(Sne!=mocommon_.par) Sne_=Sne;
  if(Snm!=mocommon_.par) Snm_=Snm;
  if(Snl!=mocommon_.par) Snl_=Snl;

  return calcSpectrum(*key,Sg_, Se_, Sp_, Sne_, Snm_, Snl_ ,err); 
}

double  spectrinfo_(double *Xmin,double*tab, double*Etot)
{ return spectrInfo(*Xmin,tab, Etot);}


double zinterp_(double*x, double*tab) {  return zInterp(*x, tab);}

double spectdnde_(double *E, double *tab){ return SpectdNdE(*E, tab); }   
 

double halofactor_(double *fi,double * dfi){ return HaloFactor(*fi,*dfi); }

void gammafluxtab_(double *fi, double * dfi, double *sigmaV, double *Sp, double *Sobs)
{  gammaFluxTab(*fi, *dfi, *sigmaV,Sp, Sobs);}

double gammaflux_(double *fi, double *dfi,  double *dSigmadE)
{ return gammaFlux(*fi, *dfi, *dSigmadE);}

void posifluxtab_(double* Emin, double*sigmav, double *tab ,double *tabOut)
{  posiFluxTab(*Emin, *sigmav, tab, tabOut); }   

void pbarfluxtab_(double* Emin, double*sigmav, double *tab,double *tabOut)
{    pbarFluxTab(*Emin, *sigmav, tab,tabOut); }   


extern void  solarmodulation_(double *PHI, double *mass, double * inTab, double * outTab)
{   solarModulation(*PHI, *mass, inTab, outTab);}

int basicspectra_(double *Mass, int *pdgN, int *outN, double * tab)
{ return basicSpectra(*Mass,*pdgN, *outN, tab);}

int basicnuspectra_(int*forSun, double *Mass, int *pdgN, int*pol, double * nu, double * nuB)
{ return basicNuSpectra(*forSun,*Mass, *pdgN,*pol, nu, nuB);}


double  spectrint(double *Emin,double *Emax, double * tab)
{
  return spectrInt(*Emin,*Emax, tab);
}

static   double(*fFunc)(double *x);
static   double cFunc(double x){ return fFunc(&x);}


void spectrmult_( double *spect, double(*func)(double*))
{ 
  fFunc=func; 
  spectrMult(spect, cFunc);
}
  