/*====== Modules ===============
   Keys to switch on 
   various modules of micrOMEGAs  
================================*/

#define XENON1T
//#define DS_50
//#define PICO
//#define CRESST_III
//#define LowMassMediator  
//#define milliCharge
//#define VELOCITY
//#define SD_LIMITS
//#define Xenon1tFit

/*===== End of DEFINE  settings ===== */

#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"

typedef double ( * funcArg)(double x, void * arf);


static double CRESST_SI_NEW(double M)
{ 
  double X[14]= {134.0, 135.5, 139.5, 145.5, 152.5, 161.5, 169.5, 180.5, 199.0, 220.5, 235.0, 261.0, 279.0, 293.0 };  
  double Y[14]= { 91.0,  92.5, 101.5, 115.5, 128.0, 141.0, 149.0, 159.5, 169.0, 176.5, 183.0, 196.0, 203.5, 208.0 };

  double  x=293+log10(M)*(293-93);
  if(x<X[0] || x>X[13]) return sqrt(-1);
  double  y=polint3(x,14,X,Y);
  return  1E-36*pow(10, -10*(y-156.5)/(399-156.5));
}

static double CRESST_SD_NEW(double M)
{ 
  double X[13]= { 115.5, 121.5, 128.5, 136.0, 142.0, 150.5, 161.0, 177.5, 194.0, 220.0, 243.5, 260.5, 274.5};  
  double Y[13]= {  88.5, 103.0, 120.0, 139.0, 148.5, 157.5, 171.0, 179.5, 187.5, 200.5, 210.5, 219.5, 224.5};

  double  x=274+log10(M)*(274-74);
  if(x<X[0] || x>X[12]) return sqrt(-1);
  double  y=polint3(x,13,X,Y);
  return  1E-36*1E4*pow(10, 6*(y-228.5)/(65-228.5));
}

static double PICO_SI_1(double M)
{
  double X[10]={274, 284,306,360,425,489,560,627,697};
  double Y[10]={43,  140,248,355,403,424,427,424,414};
  
  double x=119+log10(M)*(704-119)/2;
  if(x<X[0] || x>X[8]) return NAN;
  double y=polint3(x,9,X,Y);
  return 1E-44*pow(10,6*(y-493)/(43-493));
}

#if defined(Xenon1tFit)
#define FIT
#endif

#ifdef FIT

#define  exposure1t (900*279)

static double (*SigmaExp)(double M)=XENON1T_90;
static double Exposure=exposure1t; 

static double (*dNdE_fact)(double E)=NULL;

double kappaReg=1;
int chi2Flag=0;

static double sigma90_1t(double M)
{ 
// 1805.12562  
/*
  double   data_1t_M[25]=  { 5.988E+00, 6.589E+00, 7.433E+00, 8.565E+00, 9.942E+00, 1.158E+01, 1.344E+01, 1.617E+01, 1.931E+01, 2.257E+01, 2.716E+01, 3.130E+01, 3.712E+01, 4.450E+01, 5.667E+01, 7.140E+01, 9.761E+01, 1.349E+02, 1.774E+02, 2.307E+02, 3.199E+02, 4.437E+02, 5.834E+02, 6.993E+02, 9.837E+02};
  double   data_1t_S[25]=  { 2.830E-44, 1.028E-44, 3.667E-45, 1.307E-45, 5.570E-46, 2.757E-46, 1.585E-46, 8.946E-47, 6.033E-47, 4.861E-47, 4.262E-47, 4.143E-47, 4.340E-47, 4.632E-47, 5.637E-47, 6.860E-47, 8.997E-47, 1.214E-46, 1.577E-46, 1.974E-46, 2.738E-46, 3.799E-46, 4.936E-46, 6.008E-46, 8.334E-46};
*/
  double   data_1t_lnM[25]={ 1.790E+00, 1.885E+00, 2.006E+00, 2.148E+00, 2.297E+00, 2.449E+00, 2.598E+00, 2.783E+00, 2.961E+00, 3.117E+00, 3.302E+00, 3.444E+00, 3.614E+00, 3.796E+00, 4.037E+00, 4.268E+00, 4.581E+00, 4.904E+00, 5.178E+00, 5.441E+00, 5.768E+00, 6.095E+00, 6.369E+00, 6.550E+00, 6.891E+00};
  double   data_1t_lnS[25]={ 1.040E+00, 2.794E-02,-1.003E+00,-2.035E+00,-2.888E+00,-3.591E+00,-4.144E+00,-4.717E+00,-5.111E+00,-5.326E+00,-5.458E+00,-5.486E+00,-5.440E+00,-5.375E+00,-5.178E+00,-4.982E+00,-4.711E+00,-4.412E+00,-4.150E+00,-3.925E+00,-3.598E+00,-3.270E+00,-3.009E+00,-2.812E+00,-2.485E+00};

  return 1E-44*exp( polint3(log(M),25,data_1t_lnM,  data_1t_lnS)); 
} 


static double recoilSignal(double E,double *dNdE) // integrate recoil energy  with different distributions.     
{
   double r= dNdERecoil(E,dNdE);
   if(dNdE_fact) r*=dNdE_fact(E);
   return r;
}

static int Ntab=12;
static double Etab[100]={3.000E+00,4.000E+00,5.000E+00,6.000E+00,7.000E+00,8.000E+00,9.000E+00,1.000E+01,1.100E+01,1.200E+01,1.300E+01,1.400E+01};
static double Ptab[100]={1.188E-01,2.088E-01,3.331E-01,4.562E-01,5.176E-01,5.134E-01,4.735E-01,3.896E-01,2.863E-01,1.851E-01,8.795E-02,5.774E-04};
static double Ptab_[100];

static void printTabs(void)
{  printf("  int Ntab=%d;\n",Ntab);
   printf("  double Etab[%d]={%.3E",Ntab,Etab[0]);
   for(int i=1;i<Ntab;i++) printf(",%.3E",Etab[i]); 
   printf("};\n");
   printf("  double Ptab[%d]={%.3E",Ntab,Ptab[0]);
   for(int i=1;i<Ntab;i++) printf(",%.3E",Ptab[i]); 
   printf("};\n");
}

static double p0Z_fit(double E)
{ 
  if(E<Etab[0]) return 0;
  if(E>Etab[Ntab-1]) return 0;
  double r= polint3(E,Ntab,Etab,Ptab);
  if(r<=0) return 0; 
  return r;
} 

static double p0Z_fit_(double E)
{ 
  if(E<Etab[0]) return 0;
  if(E>Etab[Ntab-1]) return 0;
  double r= polint3(E,Ntab,Etab,Ptab_);
  if(r<=0) return 0; 
  return r;
} 


 
static double maxErecoil(double M) { return 1E6*2*131*0.94* pow( M*(vEsc+vEarth)/299792./(131*0.94+M),2);}


static double maxDev(int prn) 
{
  double dNdE[RE_DIM];  
  double chi2=0;
  double McdmMem=Mcdm;
  dNdE_fact=p0Z_fit;
  for(Mcdm=6; Mcdm<1000; Mcdm*=1.2)  
  { double cs=SigmaExp(Mcdm)*1E36;

    nucleusRecoilCS(Maxwell,131,Z_Xe,0*J_Xe131,SxxXe131,cs,cs,0,0,dNdE);
    double Emax=maxErecoil(Mcdm);
    if(Emax>Etab[Ntab-1]) Emax=Etab[Ntab-1];
    double nEvents=Exposure*simpson_arg(recoilSignal,dNdE,Etab[0]+0.001, Emax, 1E-3,NULL);
  
    if(prn) printf("Mcdm=%.3E nEvents=%.3E  \n",Mcdm,nEvents );
    double dchi2=fabs(nEvents  -log(10));
    if(chi2Flag) chi2+=dchi2*dchi2;
    else {if(dchi2>chi2) chi2=dchi2;}
  }
  if(prn) printf("maxDev=%.2e\n",chi2); 
  Mcdm=McdmMem;  
  if(chi2Flag)  return sqrt(chi2)/Ntab; else return chi2;
}

static double fitFunc(double *x) 
{ 
   Etab[0]=fabs(x[0]);
   for(int i=1;i<Ntab-1;i++) Ptab[i]=fabs(x[i]);
   double r1= maxDev(0),r2=0;
   for(int i=1;i<=Ntab-2;i++)  r2+= 0.5*(Etab[i+1]-Etab[i-1])
   *pow(
          Ptab[i-1]/(Etab[i-1]-Etab[i]  )/(Etab[i-1]-Etab[i+1])
        + Ptab[i]  /(Etab[i]  -Etab[i-1])/(Etab[i]  -Etab[i+1])  
        + Ptab[i+1]/(Etab[i+1]-Etab[i])  /(Etab[i+1]-Etab[i-1]) 
       ,2);
   return r1+kappaReg*r2;
}



#endif 


static double Xenon1TE_lm(double M)
{ 
// 1907.11485 
double Y[22]= {6.5 , 29,   50.8, 76.2,  100.0, 121.5, 141.8, 160.0, 174.8, 185.8, 190.2, 191.2, 191.2, 191.2, 199.5, 209.8, 210.5, 210.8, 211.0, 211.0, 207.5, 202.0};
double X[22]= {44.2, 49.0, 56.0, 65.5,  78.2,  91.8,  110.5, 132.2, 160.8, 194.5, 227.0, 246.2, 250.2, 253.5, 255.0, 255.2, 258.2, 262.2, 269.8, 299.0, 336.2, 388.2};

double  x=49+ log10(M/3)*(369-49);
if(x<X[0]|| x>X[21]) return sqrt(-1);
double  y=polint3(x,22,X,Y);
return 1E-36*pow(10, -3*(y-35)/(197-35));
}




static double medMass=100; // mediator mass. To test case of light mediator.

static double myDDfactor( double E, int A) {  return pow( 1/(1+2*(A*0.94)*(E*1E-6)/medMass/medMass),2); }
#define m_ph 1
static double myPhotonFactor(double E,  int A) { return pow( m_ph*m_ph/(2*(A*0.94)*(E*1E-6)),2); } 

static double pXe1T_D(double E)
{

// 1805.12562 : Fig.1 S1 Detection
  double X[]={ 91.5,  97.0, 102.0, 106.0, 111.0, 116.5, 122.0, 127.5, 136.5, 147.5, 159.0, 180.5, 206.5, 242.5, 292.0, 363.0, 520.0, 756.0 };
  double Y[]={461.5, 459.5, 451.5, 436.5, 417.5, 384.0, 349.5, 309.5, 254.0, 190.0, 142.5,  95.5,  71.5,  53.5,  48.0,  43.3,  43.5,  45.5};
 
  double x=81+ (580-81)*E/50.;
  
  if(x<X[0]) return 0;
  if(x>X[17])return 0; 
//printf("x=%E\n",x);
  double y=polint3(x,18,X,Y);
  return  (462.5-y)/(462.5-30); 

}


static double pXe1T(double E)
{

// 1805.12562 : Fig.1 S1 detection and selection
  double X[]={ 96.5, 107,   115.5, 121.5,  127.5, 133.0, 140.5, 150.5, 161.0, 173.5, 193.5, 258, 254.5, 304.5, 346.5, 393.5, 429.5, 462, 490, 514,529,548,  571, 598, 616.5, 644.5, 672,  688 };
  double Y[]={462.5, 448.5, 421.5, 389.5,  353.5, 319.5, 274.0, 220.5, 183  , 151.0, 125.5,  95,  95.5,  91.5,  90,   91.5 , 104.5, 134, 182, 243,282,332.5,381, 420, 437,   451.5, 458.5,459 };
 
  double x=83+ (779-83)*E/70.;
  if(x<X[0]) return 0;
  if(x>X[27])return 0; 

  double y=polint3(x,28,X,Y);
  return  (462.5-y)/(462.5-31); 

}

static double XenonEff1(double E)
{
  double X[]={19.975  ,18.05   ,16.925  ,15.625  ,14.075  ,12.6    ,11.05   ,9.975   ,8.95    ,8.275   ,7.625   ,7.2     ,6.8     ,6.4     ,6.025   ,5.725   ,5.525   ,5.125   ,4.575   ,4.125   ,3.675   ,3.225   ,2.8     ,2.35    ,1.95};
  double Y[]={ 0.8480726, 0.84126985, 0.8344671, 0.82539684, 0.81179136, 0.7913832, 0.7619048, 0.7324263, 0.69161, 0.6507937, 0.6099773, 0.5782313, 0.53968257, 0.49433106, 0.4376417, 0.39909297, 0.3537415, 0.292517, 0.2244898, 0.1632653, 0.10204082, 0.06575964, 0.03628118, 0.015873017, 0.};

  if(E>X[0]) return NAN;
  if(E<X[24]) return NAN;
  return polint3(E,25,X,Y);
}

static double XenonEff2(double E)
{

  double X[]={1.15 ,1.7    ,2.2   ,2.7    ,3.15   ,3.675  ,4.0    ,4.375  ,4.75   ,5.15   ,5.725  ,6.175  ,6.65   ,7.25   ,7.8    ,8.25   ,8.875  ,9.725  ,10.725 ,11.825 ,12.975 ,13.9   ,15.1   ,16.275 ,17.425 ,19.45,20.475}; 
  double Y[]={ 0.004535147, 0.013605442, 0.034013607, 0.07256236, 0.11337868, 0.18367347, 0.23356009, 0.2834467, 0.35147393, 0.414966, 0.49659863, 0.5555556, 0.6122449, 0.6553288, 0.6938776, 0.7210884, 0.7505669, 0.7777778, 0.8027211, 0.8208617, 0.83900225, 0.8503401, 0.8639456, 0.86848074, 0.8730159, 0.8798186, 0.8798186};
  if(E>X[26]) return NAN;
  if(E<X[0]) return NAN;
  return polint3(E,27,X,Y);
}




static double dXdE_milli(double E)
     {
//  Earth          O        Na         Mg         Al       Si          P           S         Ca         Cr         Fe         Ni           

/*               O  Na  Mg   Al  Si  P   S   Ca  Cr, Fe  Ni  */
       double  A[11] = {16, 23 ,24, 27, 28, 30, 32, 40, 52, 56, 58};
       double  Z[11] = { 8, 11 ,12, 13, 14, 15, 16, 20, 24, 26, 28};
       double ro[11]=  { 4.40E-01,  2.70E-03,  2.28E-01,  2.35E-02,  2.10E-01,  9.00E-05,  2.50E-04,  2.53E-02,  2.60E-03,  6.26E-02,  1.96E-03 };

       double sum=0;
       if(E<=0) return 0;
       double v=sqrt(2*E/Mcdm);
//printf("E=%E v=%E\n",E,v);       
       for(int i=0;i<11;i++)
       {  double MA=A[i]*0.94;
          double RA= 0.8853/(pow(Z[i],1./3.)*0.0005/137);
          double mu=Mcdm*MA/(Mcdm+MA);
          double k=pow(2*v*mu*RA,2);
//printf("k=%E\n",k);          
          if(k>0.01) sum+= Z[i]*Z[i]/(MA*MA)*(log(1+k) -k/(1+k))*ro[i];    
          else       sum+= Z[i]*Z[i]/(MA*MA)*(k*k/2.)*ro[i];
       }
       sum*=2*M_PI* pow( 4*M_PI/137/v,2)*2.7;
       sum*=3.8937966E-28/1.782E-24*1E2; 
//       sum-=Mcdm*1E-11;//  1E-13;
//printf("v=%E sum=%E\n",v,sum);       
       return 1/sum; 
     }


int main(int argc,char** argv)
{  int err;

   dNdEfact=NULL;
   double dNdEarr[6][RE_DIM];
   char mess[6][30];
   double sigma[6][100];


   double mMin=3,mMax=3000;    // DM mass range
   double csTest=5E-40*1E36;   // impact cross section. Precise value does not matter
   double dNdE[RE_DIM];
   double nEvents,f;

   double csTabX[100], csTab0[100],csTab1[100],csTab2[100],csTab3[100],csTab4[100];
   char*expName="";

   sortOddParticles(NULL);

#ifdef XENON1T
{

displayPlot("XENON", "E", 1,20,0,4,"Detection",0, pXe1T_D,NULL,"Selection",0, pXe1T,NULL,"down",0,XenonEff1,NULL,"up",0, XenonEff2,NULL); 

displayPlot("XENON", "E", 1,20,0,5, "up",0,XenonEff2,NULL,"down",0, XenonEff1,NULL,
                         "pEff0",0, Xe1TpEff0,NULL, "pEff1",0, Xe1TpEff1,NULL,"pEff2",0, Xe1TpEff2,NULL);  



Mcdm=9.891E+02;

csTest=XENON1T_90(Mcdm)*1E36;


// the best-fit point
   Mcdm=200;
   csTest=4.7E-11; //pb Table I
   nucleusRecoilCS(Maxwell, 130,Z_Xe,0,NULL,csTest,csTest,0,0,dNdEarr[0]);
   for(int i=0;i<RE_DIM;i++) { double E=RE_START*pow(RE_STEP,i);  dNdEarr[0][i]*=1300*279*pXe1T(E);}
   nEvents=simpson_arg((funcArg)dNdERecoil,dNdEarr[0], 0, 60, 1E-3,NULL);
   printf("nEvents for the best-fit point in 1.3t volume(before cut) = %.2E ( 3.56 Table I)\n", nEvents);
   
   printf("\n");
   printf("Number of events after cuts for 90%% excluded cross section %.2E\n",1.16/*Table I*/ *XENON1T_90(Mcdm)*1E36/csTest);
   printf("Number of events for 90%% exclusion according to   Feldman-Cousins %.2E\n",FeldmanCousins(2,1.62,0.9));

// low mass
   Mcdm=7; 
   csTest=XENON1T_90(Mcdm)*1E36;
   nucleusRecoilCS(Maxwell, 130,Z_Xe,0,NULL,csTest,csTest,0,0,dNdEarr[1]);
   for(int i=0;i<RE_DIM;i++) { double E=RE_START*pow(RE_STEP,i);  dNdEarr[1][i]*=900*279*pXe1T(E);}
   nEvents=simpson_arg((funcArg)dNdERecoil,dNdEarr[1], 1,10 , 1E-3,NULL);
   printf("\n(nEvents for  Mcdm=%.2E and   90%% excluded  cross section %.2E[pb] in 0.9t volume) = %.2E : not enough for 90%% exclusion\n"
                      ,Mcdm,csTest,nEvents);

// M35-M200 recoil energy spectra

   Mcdm=35;
   csTest=XENON1T_90(Mcdm)*1E36;
   nucleusRecoilCS(Maxwell, 130,Z_Xe,0,NULL,csTest,csTest,0,0,dNdEarr[2]);
   sprintf(mess[2],"Mcdm=%.2E cs=%.2E[pb]",Mcdm,csTest);
   Mcdm=200;
   csTest=XENON1T_90(Mcdm)*1E36;
   nucleusRecoilCS(Maxwell, 130,Z_Xe,0,NULL,csTest,csTest,0,0,dNdEarr[3]);
   sprintf(mess[3],"Mcdm=%.2E cs=%.2E[pb]",Mcdm,csTest);
    
   for(int i=0;i<RE_DIM;i++)
   {  
     double E=RE_START*pow(RE_STEP,i);
     dNdEarr[2][i]*=900*279*pXe1T(E);
     dNdEarr[3][i]*=900*279*pXe1T(E);   
   }

   displayPlot("Xenon signals*pXe","E[keV]",1,33,0,2  // Fig.1 Right  
     ,mess[2], 0,(funcArg)dNdERecoil,dNdEarr[2]
     ,mess[3], 0,(funcArg)dNdERecoil,dNdEarr[3] 
              );                     
//  improved acceptance 
   displayPlot("Efficiency","E[keV]",1,25,0,4,"pXe",0,pXe1T,NULL,"pEff0",0, Xe1TpEff0,NULL, "pEff1",0, Xe1TpEff1,NULL,"pEff2",0, Xe1TpEff2,NULL);  //  Fig.2 Left

   mMin=6; mMax=1000;
   csTest=XENON1T_90(35)*1E36;
   
   for(int i=0;i<100;i++)
   { Mcdm=mMin*pow(mMax/mMin, (i+0.5)/100.);

     Xe1TnEvents=0;
      f= DD_factorCS(XENON1T_2018,0.1,Maxwell, csTest, csTest,0,0,&expName);
      csTab1[i]=1E-36*f*csTest;
     Xe1TnEvents=1;
      f= DD_factorCS(XENON1T_2018,0.1,Maxwell, csTest, csTest,0,0,&expName);
      csTab2[i]=1E-36*f*csTest;
     Xe1TnEvents=2; 
      f= DD_factorCS(XENON1T_2018,0.1,Maxwell, csTest, csTest,0,0,&expName);
      csTab3[i]=1E-36*f*csTest; 
   }

   displayPlot("Xenon1T SI. 90\% and 3 sigma ","M[GeV]", mMin,mMax,1,4 
                                                 ,"Xenon 1T(exp.)  90\%",0,XENON1T_90,NULL          
                                                 ,"micrO(0 event 90\%)",100,csTab1,NULL
                                                 ,"micrO(1 event 90\%)",100,csTab2,NULL
                                                 ,"micrO(2 event 90\%)",100,csTab3,NULL
                                             ); 

}
#endif

#ifdef DS_50
{ 

  extern int DS_neStart;         //   =4;    first experimental bin.  
  extern double DS_E1ion;        //    0.0195 KeV DS_50 ;  15.76 ev - WWW. More likely DS takes into account probability of recombination.  
  extern int neFluct;            //  1 - binominal fluctuations, 2 - improved Poisson.

  dNdEfact=NULL;  
  
  csTest=1E-5;

  mMin=1.8; mMax=6;
  for(int i=0;i<100;i++)
  { Mcdm=mMin+ (i+0.5)/100.*(mMax-mMin);

DS_neStart=7;
    f= DD_factorCS(DarkSide_2018,0.1,Maxwell, csTest, csTest,0,0,&expName);
    csTab1[i]=1E-36*f*csTest;                

DS_neStart=4;
     double f= DD_factorCS(DarkSide_2018,0.1,Maxwell, csTest, csTest,0,0,&expName);
    csTab0[i]=1E-36*f*csTest;
            
  } 

  displayPlot("DS_50 90\% exclusion (Bin.Fluct.)","M[GeV]", mMin,mMax,0,4
                                                        ,"DS_50_Bin",0,DS50_90,NULL
                                                        ,"DS_50_noBin",0,DS50_90_noB,NULL            
                                                        ,"micrOMEGAs ne>=4",100,csTab0,NULL
                                                        ,"micrOMEGAs ne>=7",100,csTab1,NULL
             );
                                                                                     
  mMin=0.4,mMax=2;
  for(int i=0;i<100;i++)
  { Mcdm=mMin+ (i+0.5)/100.*(mMax-mMin);
    neFluct=1;
    double f= DD_factorCS(DarkSide_2018,0.1,Maxwell, csTest, csTest,0,0,&expName);
    csTab1[i]=1E-36*f*csTest;
    neFluct=2;
    f= DD_factorCS(DarkSide_2018,0.1,Maxwell, csTest, csTest,0,0,&expName);
    csTab2[i]=1E-36*f*csTest; 
  } 
                    
  displayPlot("DS_50 90\% exclusion (No Bin.Fluct.)","M[GeV]", mMin,mMax,0,4
                                                        ,"DS_50",0,DS50_90,NULL
                                                        ,"micrO Binominal",100,csTab1,NULL
                                                        ,"micrO Poisson",100,csTab2,NULL
                                                        ,"CRESST_III",0,CRESST_III_90,NULL
                                                      ); 
   dNdEfact=NULL;
}
#endif 


#ifdef PICO
{ 
   
//  displayPlot("PICO_1","M", 3,100,1,1,"exp",0,PICO_SI_1,NULL);
//  extern double PICO_Eff_F0(double E);
//  extern double PICO_Eff_C0(double E);
//  displayPlot("Eff0","E",0,10,0,2,"EffF0",0, PICO_Eff_F0,NULL, "EffC0",0, PICO_Eff_C0,NULL);

  mMin=3;mMax=1000;
  csTest=PICO60_90(10)*1E36;
  for(int i=0;i<100;i++)
  {  Mcdm=mMin*pow(mMax/mMin, (i+0.5)/100.);
     PICO60Flag=0;
     f= DD_factorCS(PICO_2019,0.1,Maxwell, csTest, csTest,0,0,&expName); 
     csTab1[i]=1E-36*f*csTest; 
     PICO60Flag=1;
     f= DD_factorCS(PICO_2019,0.1,Maxwell, csTest, csTest,0,0,&expName);
      csTab2[i]=1E-36*f*csTest;    
  } 

                      
  displayPlot("PICO-60-Spin Independent 90%% exclusion","M[GeV]", mMin,mMax,1,3  //Fig.8
                 ,"PICO60_SI_new",0, PICO60_90,NULL
                 ,"micrO(FC)",     100, csTab1 ,NULL
                 ,"micrO(N)", 100, csTab2,NULL 
                                                );
// spinDepFlag=EFT;
#define SD 
#ifdef SD                     
  setSpinDepFF(PICO_2019,EFT);
  csTest= PICO60_SDp_90(10) * 1E36;
  for(int i=0;i<100;i++)
  {  Mcdm=mMin*pow(mMax/mMin, (i+0.5)/100.);
     PICO60Flag=0;
     f= DD_factorCS(PICO_2019,0.1,Maxwell, 0,0, csTest,0,&expName);  // Feldman-Cousins    
     csTab2[i]=1E-36*f*csTest;    
  } 

  displayPlot("PICO-60-Spin Dependent 90%% exclusion","M[GeV]", mMin,mMax,1,2
                 ,"PICO60_SDp_90",0, PICO60_SDp_90,NULL
                 ,"micrO", 100, csTab2,NULL 
//                 ,"micrOMEGAs(LL)", 100, csTab3,NULL

                                                );
#endif
}
#endif 

#ifdef CRESST_III
{
   mMin=0.140;mMax=30;
   

  for(int i=0;i<100;i++)
  {  Mcdm=mMin*pow(mMax/mMin, (i+0.5)/100.);
     csTest=9.E-7;
     f= DD_factorCS(CRESST_2019,0.1,Maxwell, csTest, csTest,0,0, &expName);     
     csTab1[i]=1E-36*f*csTest;  
  }

extern double CRESST_III_90_(double);
   
  displayPlot("CRESST-Spin Independent 90%% exclusion","M[GeV]", mMin,mMax,1,2
                 ,"CRESST_III_90",0, CRESST_III_90,NULL
//                 ,"CRESST_SI_OLD",0, CRESST_III_90_,NULL
                 ,"micrOMEGAs", 100, csTab1,NULL 
             );

  for(int i=0;i<100;i++)
  {  Mcdm=mMin*pow(mMax/mMin, (i+0.5)/100.);
     csTest=6.7E3;
     f= DD_factorCS(CRESST_2019,0.1,Maxwell, 0,0,0, csTest, &expName);     
     csTab2[i]=1E-36*f*csTest;     
  } 
  
extern double  CRESST_III_SDn_90_(double);

  displayPlot("CRESST-Spin Dependent(N) 90%% exclusion","M[GeV]", mMin,mMax,1,2
                 ,"CRESST_III_SDn_90",0, CRESST_III_SDn_90,NULL
//                 ,"CRESST_SD_OLD",0, CRESST_III_SDn_90_,NULL
                 ,"micrOMEGAs", 100, csTab2,NULL
             );

/*             
   Mcdm=0.19;
   double dNdE[RE_DIM];
   csTest=CRESST_III_90(Mcdm)*1E36;
   nucleusRecoilCS(Maxwell,16,Z_O,0,NULL,csTest,csTest, 0,0,dNdE);

   displayPlot("dNdE","E[kev]", 0.015,0.035,0,1,"O16 Mcdm=0.19",0,(funcArg)dNdERecoil,dNdE);             
*/
double X[20]={0.16   , 0.208     ,0.2704    ,0.35152   ,0.456976  ,0.594069  ,0.772289  ,1         ,1.5       ,2.25      ,3.375     ,5.0625    ,7.59375   ,11.3906   ,17.0859   ,25.6289   ,38.4434   ,57.665    ,86.4976   ,129.746};

for(int i=0;i<20;i++) { 
Mcdm=X[i];
  csTest=6.7E3;
  f= DD_factorCS(CRESST_2019,0.1,Maxwell, csTest, csTest,0,0, &expName);     
  double resSI=1E-36*f*csTest;
  f= DD_factorCS(CRESST_2019,0.1,Maxwell, 0, 0,0, csTest, &expName);
  double resSD=1E-36*f*csTest;  
//  printf("M=%.2E   diffSI=%6.3f   diffSD=%6.3f \n", X[i],  resSI/CRESST_III_90(X[i])-1, resSD/CRESST_III_SDn_90(X[i])-1); 
}

}
#endif


#ifdef LowMassMediator
{  
   dNdEfact=myDDfactor;
//   double dNdEarr[10][RE_DIM];   
   
   double MA=130*0.94;

   int N;

   csTest=XENON1T_90(30);
   for(int k=0;k<3;k++)
   {  Mcdm=10*pow(3,k);
      sprintf(mess[2*k]," Mcdm=%.2E, no event",Mcdm);
      sprintf(mess[2*k+1],"1 event");
      N=0;
      for(medMass=0.001;medMass<1;medMass*=pow(10,0.1))      
      { 
        Xe1TnEvents=0; 
        f= DD_factorCS(XENON1T_2018,0.1,Maxwell, csTest, csTest,0,0,&expName);
        sigma[2*k][N]=csTest*f*1E-36;

        Xe1TnEvents=1;
        f= DD_factorCS(XENON1T_2018,0.1,Maxwell, csTest, csTest,0,0,&expName);
        sigma[2*k+1][N]=csTest*f*1E-36;

        N++;      
      }
   }

   displayPlot("Xenon1T:  90%% exclusion for light mediator","M[MeV]",pow(10,-0.05) ,1000*medMass*pow(10,-0.05), 1,6 
                ,mess[0],  N, sigma[0], NULL
                ,mess[1],  N, sigma[1], NULL
                ,mess[2],  N, sigma[2], NULL
                ,mess[3],  N, sigma[3], NULL
                ,mess[4],  N, sigma[4], NULL
                ,mess[5],  N, sigma[5], NULL
              );      


    Mcdm=1000;
    N=0;
    for(medMass=0.001;medMass<1;medMass*=pow(10,0.1))      
    {  
        Xe1TnEvents=0;
        f= DD_factorCS(XENON1T_2018,0.1,Maxwell, csTest, csTest,0,0,&expName);
        csTest*=f;  
        sigma[0][N]=csTest*1E-36;

        Xe1TnEvents=1;
        f= DD_factorCS(XENON1T_2018,0.1,Maxwell, csTest, csTest,0,0,&expName);
        csTest*=f;  
        sigma[1][N]=csTest*1E-36;

        Xe1TnEvents=2;
        f= DD_factorCS(XENON1T_2018,0.1,Maxwell, csTest, csTest,0,0,&expName);
        csTest*=f;  
        sigma[2][N]=csTest*1E-36;
        
        N++;      
    }
                                                    
   displayPlot("Xenon1T:  90%% exclusion for light mediator. Mcdm=1000","M[MeV]",pow(10,-0.05) ,1000*medMass*pow(10,-0.05), 1,3 
                ,"0 events", N, sigma[0], NULL
                ,"1 event ", N, sigma[1], NULL
                ,"2 events", N, sigma[2], NULL
              );      


   mMin=0.2; mMax=30;
    medMass=0.001;

       
  for(int i=0;i<100;i++)
  { Mcdm=mMin*pow(mMax/mMin, (i+0.5)/100.);
    f= DD_factorCS(DarkSide_2018,0.1,Maxwell, csTest, csTest,0,0,&expName);
    csTab0[i]=1E-36*f*csTest;
    f= DD_factorCS(CRESST_2019,0.1,Maxwell, csTest, csTest,0,0,&expName);
    csTab1[i]=1E-36*f*csTest;     
    f= DD_factorCS(PICO_2019,0.1,Maxwell, csTest, csTest,0,0,&expName);
    csTab2[i]=1E-36*f*csTest;
    f= DD_factorCS(XENON1T_2018,0.1,Maxwell, csTest, csTest,0,0,&expName);
    csTab3[i]=1E-36*f*csTest;         
  } 
   
 displayPlot("mu^4 cs MeV^4*cm^2","M[GeV]", mMin,mMax,1,5
                                    ,"Xenon1T_El",0, Xenon1TE_lm,NULL   
                                    ,"DarkSide",100, csTab0,NULL
                                   ,"CRESST",  100, csTab1,NULL
                                    ,"PICO",    100, csTab2,NULL
                                    ,"XENON1T_N", 100, csTab3,NULL
            ); 


   medMass=0.01;
   
   Mcdm=10;
   csTest=XENON1T_90(Mcdm);
   nucleusRecoilCS(Maxwell,130,Z_Xe,0,NULL,csTest,csTest,0,0,dNdEarr[0]);
   Mcdm=15;
   csTest=XENON1T_90(Mcdm);
   nucleusRecoilCS(Maxwell,130,Z_Xe,0,NULL,csTest,csTest,0,0,dNdEarr[2]);
   for(int i=0;i<RE_DIM;i++) 
   {  double E=RE_START*pow(RE_STEP,i); 
      f=Xe1TpEff0(E);
      double g=myDDfactor(E,130);
//      printf("E=%e g=%E\n",E,g);
      dNdEarr[0][i]*=f;
      dNdEarr[1][i]= dNdEarr[0][i]*g;
      dNdEarr[2][i]*=f;
      dNdEarr[3][i]= dNdEarr[2][i]*g;
   }
   
   for(int n=0;n<4;n++)
   {  nEvents=simpson_arg((funcArg)dNdERecoil,dNdEarr[n],1.15,14,1E-3,&expName);
      for(int i=0;i<RE_DIM;i++)dNdEarr[n][i]*=log(10)/nEvents;   // normalization on 90% exclusion
   }  
   
   displayPlot("Xenon1T: normalized  recoil sinal",  "E[keV]",1.5, 10, 0, 3
                       ,"M_{MD}=10  Large mu", 0,  (funcArg)dNdERecoil, dNdEarr[0]
                       ,"M_{DM}=15  Large mu", 0,  (funcArg)dNdERecoil, dNdEarr[2]
                       ,"M_{DM}=15  mu=10MeV",  0, (funcArg) dNdERecoil, dNdEarr[3]
                       ); 
                    
   dNdEfact=NULL;

}
#endif 


#ifdef milliCharge
{   
   dNdEfact=myPhotonFactor;
   double qXep[100],qXem[2][100], qDSp[100],qDSm[100] ;
   double Mp=0.9383, alpha=1./137., qtest=1;
   
   int N=0;
   for(Mcdm=6;Mcdm<1000;Mcdm*=1.1,N++)
   {  
      double cs=3.8937966E8*  16*M_PI* pow( alpha*qtest* Mcdm*Mp/(Mcdm+Mp)/m_ph/m_ph ,2);
      Xe1TnEvents=0; 
      f=DD_factorCS(XENON1T_2018,0.1,Maxwell, cs,0,0,0,NULL);      
      qXem[0][N]= qtest*sqrt(f);
      Xe1TnEvents=1; 
      f=DD_factorCS(XENON1T_2018,0.1,Maxwell, cs,0,0,0,NULL);      
      qXem[1][N]= qtest*sqrt(f);          
   }
//   displayPlot("Xenon millicharge ", "Mcdm",  6/sqrt(1.1), Mcdm/sqrt(1.1), 1, 1,  "q",N, qXem,NULL); 

   N=0;
   for(Mcdm=1.8;Mcdm<6;Mcdm+=0.1,N++)
   {  
      double cs=3.8937966E8*  16*M_PI* pow( alpha*qtest* Mcdm*Mp/(Mcdm+Mp)/m_ph/m_ph ,2);
      f=DD_factorCS(DarkSide_2018,0.1,Maxwell, cs,0,0,0,NULL);      
      qDSm[N]= qtest*sqrt(f);  
   }
//   displayPlot("DarkSide millicharge", "Mcdm",  1.8-0.1/2, Mcdm-0.1/2, 0, 1,  "q",N, qDSm,NULL); 
          
   dNdEfact=NULL; 

//  Upper limit 
   // Xenon
   double MD=130;
   double Etr=1.6E-6;
       
   N=0;
   for(Mcdm=6;Mcdm<1000;Mcdm*=1.1,N++)
   {  
     double muD=Mcdm*MD/(Mcdm+MD);
     double Emin=Etr*MD*Mcdm/(4*muD*muD);
     double Emax=pow( (vEsc+vEarth)/299792.458,2)*Mcdm/2;

     double zz=simpson(dXdE_milli,Emin,Emax,1E-3,NULL);   
     qXep[N]= sqrt(zz/1400);  
   }
   displayPlot("Xenon millicharge ", "Mcdm",  6/sqrt(1.1), Mcdm/sqrt(1.1), 1, 3
              ,  "zXep",N, qXep,NULL
              ,  "zXem0",N, qXem[0],NULL
              ,  "zXem1",N, qXem[1],NULL); 

   N=0;
   MD=40;
   Etr=1E-7;
   for(Mcdm=1.8;Mcdm<6;Mcdm+=0.1,N++)
   { 
     double muD=Mcdm*MD/(Mcdm+MD); 
     double Emin=Etr*MD*Mcdm/(4*muD*muD);
     double Emax=pow( (vEsc+vEarth)/299792.458,2)*Mcdm/2;
     double zz=simpson(dXdE_milli,Emin,Emax,1E-3,NULL);   
     qDSp[N]= sqrt(zz/1400);  
   }
   displayPlot("DarkSide millicharge", "Mcdm",  1.8-0.1/2, Mcdm-0.1/2, 0, 2
            ,  "zDSp",N, qDSp,NULL
            ,  "zDSp",N, qDSm,NULL);                 
}
#endif 


#ifdef milliCharge2
{

     double z_milli=2.113E-5 ;
z_milli*=1.0001;     
     Mcdm=1000;     

     
double MD=130; //Xenon    
double Etr=1.6E-6;
double muD=Mcdm*MD/(Mcdm+MD);
double E=1.6E-6*MD*Mcdm/(4*muD*muD);


double Emin=Etr*MD*Mcdm/(4*muD*muD);
double Emax=pow( (vEsc+vEarth)/299792.458,2)*Mcdm/2;

double zz=simpson(dXdE_milli,Emin,Emax,1E-3,NULL);

printf("z-milli= %E\n", sqrt(zz/1400));
exit(0);     
//     double E= pow( (vEsc+vEarth)/299792.458,2)*Mcdm/2;
// E=1E-8;     
     double lnE=log(E);

 int n;
 double l0=0.1;
 double vv[1000],ll[1000];

 double dEdX(double E)
     {
//  Earth          O        Na         Mg         Al       Si          P           S         Ca         Cr         Fe         Ni           

/*               O  Na  Mg   Al  Si  P   S   Ca  Cr, Fe  Ni  */
       double  A[11] = {16, 23 ,24, 27, 28, 30, 32, 40, 52, 56, 58};
       double  Z[11] = { 8, 11 ,12, 13, 14, 15, 16, 20, 24, 26, 28};
       double ro[11]=  { 4.40E-01,  2.70E-03,  2.28E-01,  2.35E-02,  2.10E-01,  9.00E-05,  2.50E-04,  2.53E-02,  2.60E-03,  6.26E-02,  1.96E-03 };

       double sum=0;
       if(E<=0) return 0;
       double v=sqrt(2*E/Mcdm);
//printf("E=%E v=%E\n",E,v);       
       for(int i=0;i<11;i++)
       {  double MA=A[i]*0.94;
          double RA=  pow(Z[i],1./3.)/(0.0005/137); //0.8853
          double mu=Mcdm*MA/(Mcdm+MA);
          double k=pow(2*v*mu*RA,2);
//printf("k=%E\n",k);          
          if(k>0.01) sum+= Z[i]*Z[i]/(MA*MA)*(log(1+k) -k/(1+k))*ro[i];    
          else       sum+= Z[i]*Z[i]/(MA*MA)*(k*k/2.)*ro[i];
       }
       sum*=2*M_PI* pow( 4*M_PI/137*z_milli/v,2)*2.7;
       sum*=3.8937966E-28/1.782E-24*1E2; 
//       sum-=Mcdm*1E-11;//  1E-13;
//printf("v=%E sum=%E\n",v,sum);       
       return  sum; 
     }



//#define ODEINT
#ifdef ODEINT     
     void deriv(double x,double *lnE, double *dE)
     {  double E=exp(*lnE);
        *dE=dEdX(E)/(E);
     }
     odeint(&lnE,1 ,0 , 1.4 , 1.E-4, 0.1, deriv);
#else        
     void stiff_derivs(double x, double*lnE, double*dlnE, double h,double*dFdx,double*dFdY)
     {  double E=exp(*lnE);
        *dlnE= dEdX(E)/E;
        if(dFdx)  *dFdx=0;
        if(dFdY)  *dFdY= ( dEdX(E*1.1)/E/1.1 -dEdX(E/1.1)/E*1.1)/log(1.1*1.1)  ;
//if(dFdY)         printf("x=%E  E=%E   *dlnE=%E FdY=%E \n",x,E,*dlnE,*dFdY ) ;
     }  
 double h=0.1;
 double Yscal=1;
 
 printf("vEarth+vEsc=%E  =%E\n",vEarth+vEsc,(vEarth+vEsc)/300000);
 
 vv[0]=sqrt(2*exp(lnE)/Mcdm);
 printf("vv[0]=%E\n",vv[0]);
 ll[0]=log(l0);
 n=0;
 for(n=1; exp(lnE)< pow((vEarth+vEsc)/299792.458,2)*Mcdm/2 ;n++) 
 {
    double l=l0*pow(1.1,n-1);
    stiff(1,l, l*1.1, 1, &lnE, &Yscal, 1E-4, &h, stiff_derivs); 
    vv[n]= sqrt( 2*exp(lnE)/Mcdm);
    ll[n]=log(l*1.1);
printf("n=%d ll=%e (%e) vv=%E\n",n,ll[n],l*1.1,vv[n]);    
 } 
#endif    
  displayPlot("v(l)","l[m]",l0/sqrt(1.1),l0*pow(1.1,n-0.5),1,1,"v(l)",n,vv,NULL);
printf("n=%d\n",n);   

double vSurface(double vD, double l)
{
//   if( vD<vv[0] || vD>  (vEarth+vEsc)/299792.458  ) return 0;
   
   double l1=exp(polint3(vD,n,vv,ll));
   return polint3(log(l1+l),n,ll,vv);

// printf("vD=%E, l1=%E log(l1)=%E  vD'=%E vS=%e\n",vD,l1,log(l1), polint3(log(l1),n,ll,vv) ,    vS);

}

 double fDarr[1000];
 
 
 double vDstat;
 
 double cosIntegrand(double cs)
 {
  double vSp= vSurface(vDstat*(1.01), 1400/cs);
  if(vSp>(vEarth+vEsc)/299792.458) return 0; 
  double vSm= vSurface(vDstat*(0.99), 1400/cs);
  double vS = vSurface(vDstat       , 1400/cs);
  double r=  299792.458*Maxwell(vS*299792.458)*(vSp-vSm)/(0.02*vDstat);
//  printf("cs=%E  vSp=%E vSm=%E\n", cs,vSp, vSm);
  return r;
  
 }

 for(int i=0;i<n;i++){ vDstat=vv[i]; fDarr[i]=0.5*simpson(cosIntegrand, 0.1,1,1E-3,NULL);
           printf("f(%E)=%E\n", vv[i], fDarr[i]);
 }
 
 double MaxwellD(double v) 
 {  double v1=v/299792.458;
 printf("v1=%e vv[0]=%E vSurface(vv[0], 1400)=%e\n",v1,vv[0],vSurface(vv[0], 1400));
 
    if(v1<vv[0]) return 0;
//    if(v1>vSurface(vv[0], 1400)) return 0;
    double r=  polint3(v1,n,vv,fDarr)/299792.458 ;
    if(r<0) r=0;
    return r;
 }    
 
 double ff=MaxwellD(100);
 printf("ff(100)=%e\n",ff);
 
 displayPlot("f(v)","v", vv[0]*299792.458, vEarth+vEsc,0,2,"Naxwell",0, Maxwell,NULL,"MaxwellD", 0,MaxwellD,NULL);

 double Mp=0.9383, alpha=1./137.; 
 dNdEfact=myPhotonFactor;
   double cs=3.8937966E8*  16*M_PI* pow( alpha*z_milli* Mcdm*Mp/(Mcdm+Mp)/m_ph/m_ph ,2);
   double pval=DD_pvalCS(XENON1T_2018,MaxwellD, cs,0,0,0,NULL);
   
   printf("pval=%E\n",pval);
            
//    E=exp(lnE);
 
//printf("Eend=%E Edm=%E\n", E, pow( (vEsc+vEarth)/299792.458,2)*Mcdm/2  );   


//     printf("vStart=%E  Estart=%E    E_end=%E, vEnd=%E  vtr=%E\n",  (vEsc+vEarth)/299792.458,  pow( (vEsc+vEarth)/299792.458,2)*Mcdm/2,   E, sqrt(2*E/Mcdm),sqrt(1.6E-6*MD/2)/muD);        
  
}
#endif 

#ifdef VELOCITY
{  
   double cs_min[100],cs_max[100],cs_std[100],cs_pp[100];
//   double E0=1.644270;
   rhoDM=0.3;
   int N=0; 
   int dim=5;
   csTest=1.E-9;
   dNdEfact=NULL;
   rhoDM=0.3; vRot=220; vEsc=544;  vEarth=232;  //  default
      
// DarkSide   

   for(Mcdm=1.8;Mcdm<6;Mcdm+=0.1,N++)
   {
      double Emax=maxRecoil(Mcdm);
      if(Emax>14) Emax=14;
      double nEvents,cs_;
      rhoDM=0.3; 
      for(int n1=0;n1<dim;n1++) for(int n2=0;n2<dim;n2++) for(int n3=0;n3<dim;n3++)
      {     
        vRot  =220*(1-n1/(dim-1.)) +248*n1/(dim-1.);  // vRot= 220;
        vEarth=232*(1-n2/(dim-1.)) +252*n2/(dim-1.);  // vEarth=232;
        vEsc  =498*(1-n3/(dim-1.)) +608*n3/(dim-1.);  // vEsc=544;
 
        f=DD_factorCS(DarkSide_2018,0.1,Maxwell, csTest,csTest,0,0,NULL); 
        double cs=csTest*f*1E-36;
        if(n1==0 && n2==0 && n3==0) { cs_min[N]=cs; cs_max[N]=cs;}
        else 
        { if(cs>cs_max[N]) cs_max[N]=cs;
          if(cs<cs_min[N]) cs_min[N]=cs;
        } 
      }
       
      rhoDM=0.55;  vRot=233; vEsc=580;            // SHMpp parameters
      f=DD_factorCS(DarkSide_2018,0.1,SHMpp, csTest,csTest,0,0,&expName);
      cs_pp[N]=1E-36*csTest*f;
       
      rhoDM=0.3; vRot=220; vEsc=544;  vEarth=232;  // return to default
      f=DD_factorCS(DarkSide_2018,0.1,Maxwell, csTest,csTest,0,0,&expName);  
      cs_std[N]=1E-36*csTest*f;  
   }

   displayPlot("velocity distribution uncertainty for DarkSide", "Mcdm",1.8-0.1/2, Mcdm-0.1/2,0,4
                                                                                         ,"std",   N,cs_std,NULL
                                                                                          ,"min",  N,cs_min,NULL
                                                                                          ,"max",  N,cs_max,NULL                                                                                          ,"SHM++",N,cs_pp,NULL
                                                                                          ); 
//  XENON                                                                                          
   N=0;
   for(Mcdm=6;Mcdm<1000;Mcdm*=1.2,N++)
   {
      double Emax=maxRecoil(Mcdm);
      if(Emax>14) Emax=14;
      double nEvents,cs_;
      rhoDM=0.3; 
      for(int n1=0;n1<dim;n1++) for(int n2=0;n2<dim;n2++) for(int n3=0;n3<dim;n3++)
      { vRot  =220*(1-n1/(dim-1.)) +248*n1/(dim-1.);  // vRot= 220;
        vEarth=232*(1-n2/(dim-1.)) +252*n2/(dim-1.);  // vEarth=234;
        vEsc  =498*(1-n3/(dim-1.)) +608*n3/(dim-1.);  // vEsc=544;
 
        f=DD_factorCS(XENON1T_2018,0.1,Maxwell, csTest,csTest,0,0,NULL);
        double cs=csTest*f*1E-36;
        if(n1==0 && n2==0 && n3==0) { cs_min[N]=cs; cs_max[N]=cs;}
        else 
        { if(cs>cs_max[N]) cs_max[N]=cs;
          if(cs<cs_min[N]) cs_min[N]=cs;
        } 
      }
       
      rhoDM=0.55;  vRot=233; vEsc=580;            // SHMpp parameters
      f=DD_factorCS(XENON1T_2018,0.1,SHMpp, csTest,csTest,0,0,NULL);
      cs_pp[N]=1E-36*csTest*f;
      
      rhoDM=0.3; vRot=220; vEsc=544;  vEarth=232;  // return to default
      f=DD_factorCS(XENON1T_2018,0.1,Maxwell, csTest,csTest,0,0,NULL);  
      cs_std[N]=1E-36*csTest*f;    
   }
      
   displayPlot("velocity distribution uncertainty for Xenon", "Mcdm",6/sqrt(1.2), Mcdm/sqrt(1.2),1,4
                                                                                         ,"std",   N,cs_std,NULL
                                                                                          ,"min",  N,cs_min,NULL
                                                                                          ,"max",  N,cs_max,NULL
                                                                                          ,"SHM++",N,cs_pp,NULL
                                                                                          );
}
#endif


#ifdef SD_LIMITS
{

   Mcdm=30;
   
  mMin=3; mMax=1000;
  
  for(int i=0;i<100;i++)
  {  Mcdm=mMin*pow(mMax/mMin, (i+0.5)/100.);
     csTest=PICO60_SDp_90(Mcdm);
//    spinDepFlag=EFT;
     setSpinDepFF(PICO_2019,EFT);     
     f= DD_factorCS(PICO_2019,0.1,Maxwell, 0,0, csTest*1E36,0,&expName); 
     csTab1[i]=f*csTest;  
//    spinDepFlag=SHELL;
     setSpinDepFF(PICO_2019,SHELL);
     f= DD_factorCS(PICO_2019,0.1,Maxwell, 0,0, csTest*1E36,0,&expName);
     csTab2[i]=csTest*f; 
//    spinDepFlag=SHELLm;
     setSpinDepFF(PICO_2019,SHELLm);  
     f= DD_factorCS(PICO_2019,0.1,Maxwell, 0,0, csTest*1E36,0,&expName);
     csTab3[i]=csTest*f; 
  } 

  displayPlot("PICO-60-Spin Independent 90%% exclusion","M[GeV]", mMin,mMax,1,4
                 ,"PICO60_SDp_90",0, PICO60_SDp_90,NULL
                 ,"MO(FF=Fitzpatrick et.al.)", 100, csTab1,NULL 
                 ,"MO(FF=Klos et.al.)", 100, csTab2,NULL 
                 ,"MO(FF=Klos^{min} et.al.)", 100, csTab3,NULL
                                                           );                       
  mMin=6; mMax=1000;
  for(int i=0;i<100;i++)
  {  Mcdm=mMin*pow(mMax/mMin, (i+0.5)/100.);
     csTest=XENON1T_SDp_90(Mcdm);
//     spinDepFlag=SHELL;
       setSpinDepFF(XENON1T_2018,SHELL);
       f= DD_factorCS(XENON1T_2018,0.1,Maxwell, 0,0, csTest*1E36,0,&expName);    
       csTab1[i]=f*csTest;
//     spinDepFlag=SHELLm;
      setSpinDepFF(XENON1T_2018,SHELLm);
       f= DD_factorCS(XENON1T_2018,0.1,Maxwell, 0,0, csTest*1E36,0,&expName);    
       csTab2[i]=f*csTest; 
//     spinDepFlag=EFT;
       setSpinDepFF(XENON1T_2018,EFT);
       f= DD_factorCS(XENON1T_2018,0.1,Maxwell, 0,0, csTest*1E36,0,&expName);    
       csTab3[i]=f*csTest;     
  } 
 
  displayPlot("Xenon1T-Spin Dependent #sigma_p  90%% exclusion","M[GeV]", mMin,mMax,1,4
                 ,"experiment",0, XENON1T_SDp_90,NULL
                 ,"MO (FF=Klos et.al.)", 100, csTab1,NULL
                 ,"MO (FF=Klos^{min} et.al.)", 100, csTab2,NULL
                 ,"MO (FF=Fitzpatrick et.al.)", 100, csTab3,NULL 
                                                );  

  mMin=6; mMax=1000;
  for(int i=0;i<100;i++)
  {  Mcdm=mMin*pow(mMax/mMin, (i+0.5)/100.);
     csTest=XENON1T_SDn_90(Mcdm);
//     spinDepFlag=SHELL;
     setSpinDepFF(XENON1T_2018,SHELL); 
       f= DD_factorCS(XENON1T_2018,0.1,Maxwell, 0,0,0, csTest*1E36,&expName);    
       csTab1[i]=f*csTest; 
//     spinDepFlag=SHELLm;
     setSpinDepFF(XENON1T_2018,SHELLm);
       f= DD_factorCS(XENON1T_2018,0.1,Maxwell, 0,0,0, csTest*1E36,&expName);    
       csTab2[i]=f*csTest; 
//     spinDepFlag=EFT;
     setSpinDepFF(XENON1T_2018,EFT);  
       f= DD_factorCS(XENON1T_2018,0.1,Maxwell, 0,0,0, csTest*1E36,&expName);    
       csTab3[i]=f*csTest;    
  } 
 
  displayPlot("Xenon1T-Spin Dependent #sigma_n  90%% exclusion","M[GeV]", mMin,mMax,1,4
                 ,"experiment",0, XENON1T_SDn_90,NULL
                 ,"MO (FF=Klos et.al.)", 100, csTab1,NULL
                 ,"MO (FF=Klos^{min} et.al.)", 100, csTab2,NULL
                 ,"MO (FF=Fitzpatrick et.al.)", 100, csTab3,NULL 
                                                );  
}
#endif 


#ifdef Xenon1tFit
{  // auxillary code for fitting

   double X[100],dX[100],chi2;

   dNdE_fact=p0Z_fit;
//   SigmaExp=XENON1T_90;
   SigmaExp=sigma90_1t;
   Exposure=exposure1t;
    
   double    step=0.45;
   Ntab=35;
   Etab[0]=1.4;
   for(int i=1;i<Ntab;i++) Etab[i]= Etab[0] + step*i; 
chi2Flag=1;
   Ptab[0]=0;
   for(int i=1;i<Ntab-1;i++)  Ptab[i]= /* 0; //(i<Ntab/2)? 0.01*i : 0.01*(Ntab-1-i) ;  // */  pXe1T(Etab[i]);
   Ptab[Ntab-1]=0;

//   displayPlot("Detector acceptance 1","E",0,Etab[Ntab-1],0,3,"Xenon acceptance",0, pXe1T ,NULL,"p0Z'",0,p0Z_fit ,NULL,"p0Z",0,p0Z1T_best,NULL); 
   
   X[0]=Etab[0];
   for(int i=1;i<Ntab-2;i++) X[i]=Ptab[i];

   maxDev(1);
   kappaReg=0.1;

for(int k=0;k<30;k++)
{ 
//chi2Flag=1;
//   for(int i=0;i<Ntab-1; i++) dX[i]=0.1; dX[0]/=10;
//int nCall=1000;   
//   chi2= amoeba(X, dX, Ntab-1,fitFunc,0.001, &nCall);
   for(int i=0;i<Ntab-1; i++) dX[i]=0.1; dX[0]/=10;
chi2Flag=0;
int     nCall=1000;
 chi2= amoeba(X, dX, Ntab-1,fitFunc,0.001, NULL);    
   displayPlot("Detector acceptance 1","E",0,Etab[Ntab-1],0,2,"Xenon acceptance",0, pXe1T ,NULL,"p0Z",0,p0Z_fit,NULL);  
   maxDev(1);
} 


//   chi2= amoeba(X, dX, Ntab-1,fitFunc,0.001, NULL);
//   double difO=maxDev(1);
//   char messN[30],messO[30];
//   sprintf(messN,"kappa=0.1, maxDiff=%.2E",difN);
//   sprintf(messO,"kappa=0.005, maxDiff=%.2E",difO);
   
//   displayPlot("Detector acceptance","E",0,Etab[Ntab-1],0,3,"Xenon1T",0, pXe1T ,NULL,
//   messN,0,p0Z_fit_ ,NULL,messO,0,p0Z_fit,NULL);

     
    
   printTabs();

}
#endif 

 

  killPlots();

  return 0;
}
  

        