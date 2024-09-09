/*  Direct Detection */

#include "micromegas.h"
#include "micromegas_aux.h"
#include "micromegas_f.h"


void calcscalarff_(double *muDmd,double *msDmd,double *sigmaPiN,double *sigma0)
{ calcScalarFF(*muDmd,*msDmd,*sigmaPiN,*sigma0);}

void calcscalarquarkff_(double *muDmd,double *msDmd,double *sigmaPiN,double *sigmaS)
{ calcScalarQuarkFF(*muDmd,*msDmd,*sigmaPiN,*sigmaS);}


static double(*_XYloop)(double*,double*,double*,double*);

double XYloop_(double sing, double mq, double msq, double mne)
{
   return (*_XYloop)(&sing, &mq, &msq, &mne);
}

double noloop_(double *sing, double *mq, double *msq, double *mne)
{
   return 0;
}


int nucleonamplitudes_(char * cmdF,double*pA0,double*pA5,double*nA0,double*nA5,int len)
{
  char cdmC[20];
  fName2c(cmdF,cdmC,len);
  return nucleonAmplitudes(cdmC,pA0,pA5,nA0,nA5);
}

double fermiff_(int *A, double * Qfermi)
{ return  FermiFF(*A, *Qfermi);}
 

double maxwell_(double *v){ return Maxwell(*v); }



static double(*_fDv)(double*);
static void(*_Sxx)(double*,double*,double*,double*);

static double fDv_(double v){ return (*_fDv)(&v);}
static void Sxx_(double p,double*S00,double*S01,double*S11){ (*_Sxx)(&p,S00,S01,S11);}


double nucleusrecoil_(
     double(*fDv)(double*),int*A, int*Z, double*J, 
     void(*Sxx)(double*,double*,double*,double*), double * dNdE )
{  
  _Sxx=Sxx;
  
  if(fDv  == maxwell_)  return  nucleusRecoil(Maxwell, *A,*Z,*J,Sxx_,dNdE);
  else  
  { _fDv=fDv;
    return  nucleusRecoil(fDv_, *A,*Z,*J,Sxx_,dNdE);
  }    
}

double nucleusrecoilcs_(
     double(*fDv)(double*),int*A, int*Z, double*J, 
     void(*Sxx)(double*,double*,double*,double*),
     double *LmbdP, double*XiP, double *LmbdN, double*XiN,  double * dNdE )
{  
   double (*c_fDv)(double);
  _Sxx=Sxx;
  
  if(fDv  == maxwell_) c_fDv=Maxwell;
  else { _fDv=fDv; c_fDv=fDv_;}
   
    return  nucleusRecoilCS(c_fDv, *A,*Z,*J,Sxx_,*LmbdP,*XiP,*LmbdN,*XiN, dNdE);
}



double nucleusrecoil0_( double (*fDv)(double*),
 int*A,int*Z,double*J,double*Sp,double*Sn,
 double (*LF)(double*,double*,double*,double*),double*dNdE)
{

  if(fDv  == maxwell_) return nucleusRecoil0(Maxwell,*A,*Z,*J,*Sp,*Sn,dNdE);
  else 
  { _fDv=fDv;  
    return nucleusRecoil0( fDv_,*A,*Z,*J,*Sp,*Sn,dNdE);
  }
}

double nucleusrecoil0cs_( double (*fDv)(double*),
 int*A,int*Z,double*J,double*Sp,double*Sn,
 double *LmbdP, double*XiP, double *LmbdN, double*XiN, double*dNdE)
{
  double (*c_fDv)(double);

  if(fDv  == maxwell_) c_fDv=Maxwell;
  else { _fDv=fDv; c_fDv=fDv_;}

  return nucleusRecoil0CS(c_fDv,*A,*Z,*J,*Sp,*Sn,*LmbdP,*XiP,*LmbdN,*XiN,dNdE);
}


double cutrecoilresult_(double *tab, double *E1, double *E2)
{
   return  cutRecoilResult(tab, *E1, *E2);
}

double dnderecoil_( double *E,double *tab)
{  
    return  dNdERecoil(*E,tab);
} 

 void sxxf19_   (double *p,double*S00,double*S01,double*S11){ SxxF19(*p,S00,S01,S11);}
 void sxxna23_  (double *p,double*S00,double*S01,double*S11){ SxxNa23(*p,S00,S01,S11);}
 void sxxal27_  (double *p,double*S00,double*S01,double*S11){ SxxAl27(*p,S00,S01,S11);}
 void sxxsi29_  (double *p,double*S00,double*S01,double*S11){ SxxSi29(*p,S00,S01,S11);}
 void sxxK39_   (double *p,double*S00,double*S01,double*S11){ SxxK39(*p,S00,S01,S11);}
 void sxxge73_  (double *p,double*S00,double*S01,double*S11){ SxxGe73(*p,S00,S01,S11);}
 void sxxnb93_  (double *p,double*S00,double*S01,double*S11){ SxxNb93(*p,S00,S01,S11);}
 void sxxte125_ (double *p,double*S00,double*S01,double*S11){ SxxTe125(*p,S00,S01,S11);}
 void sxxi127_  (double *p,double*S00,double*S01,double*S11){ SxxI127(*p,S00,S01,S11);}
 void sxxxe129_ (double *p,double*S00,double*S01,double*S11){ SxxXe129(*p,S00,S01,S11);}
 void sxxxe131_ (double *p,double*S00,double*S01,double*S11){ SxxXe131(*p,S00,S01,S11);}
 void sxxpb207_ (double *p,double*S00,double*S01,double*S11){ SxxPb207(*p,S00,S01,S11);}
 void sxxna23a_ (double *p,double*S00,double*S01,double*S11){ SxxNa23A(*p,S00,S01,S11);}
 void sxxsi29a_ (double *p,double*S00,double*S01,double*S11){ SxxSi29A(*p,S00,S01,S11);}
 void sxxte125a_(double *p,double*S00,double*S01,double*S11){ SxxTe125A(*p,S00,S01,S11);}
 void sxxi127a_ (double *p,double*S00,double*S01,double*S11){ SxxI127A(*p,S00,S01,S11);}
 void sxxxe129a_(double *p,double*S00,double*S01,double*S11){ SxxXe129A(*p,S00,S01,S11);}
 void sxxxe131a_(double *p,double*S00,double*S01,double*S11){ SxxXe131A(*p,S00,S01,S11);}
 void sxxge73a_ (double *p,double*S00,double*S01,double*S11){ SxxGe73A(*p,S00,S01,S11);}
 void sxxxe131b_(double *p,double*S00,double*S01,double*S11){ SxxXe131B(*p,S00,S01,S11);}
