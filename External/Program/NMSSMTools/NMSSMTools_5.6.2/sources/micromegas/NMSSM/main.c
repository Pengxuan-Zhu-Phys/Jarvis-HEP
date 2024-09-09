/*=========   NMSSM scenario  ==========
  One can define SUGRA  for GUT scale scenario
  or EWSB  for low scale input.
  Otherwise the program should read SLHA file.
=======================================*/ 

//#define SUGRA 
#define EWSB

/*====== Modules ===============
   Keys to switch on 
   various modules of micrOMEGAs
================================*/

#define MASSES_INFO      
      /* Display information about SUSY and Higgs masses 
      */
#define CONSTRAINTS     
      /* Display  deltarho, B_>sgamma, Bs->mumu, gmuon and
         check LEP mass limits 
      */ 
//#define HIGGSBOUNDS
//#define HIGGSSIGNALS
//#define SUPERISO 
//#define LILITH
//#define SMODELS
 
#define OMEGA            
      /* Calculate relic density and display contribution of
         individual channels 
      */
#define INDIRECT_DETECTION  
      /* Compute spectra of gamma/positron/neutrinos
         for DM annihilation; calculate <sigma*v> and
         integrate gamma signal over DM galactic squared
         density for given line of sight.  
      */

//#define LoopGAMMA
      
/*#define RESET_FORMFACTORS  */
      /* Modify default nucleus form factors, 
         DM velocity distribution,
         A-dependence of Fermi-dencity
      */     
//#define CDM_NUCLEON     
      /* Calculate amplitudes and cross-sections for 
         CDM-mucleon collisions 
      */  
/*#define CDM_NUCLEUS */ 
      /* Calculate number of events for 1kg*day 
         and recoil energy distibution for various nuclei
      */
#define NEUTRINO
 /*  Neutrino signal of DM annihilation in Sun and Earth */
             
#define DECAYS
    /* Calculate decay widths and branchings  */

#define CROSS_SECTIONS */
      /* Calculate cross sections and widths for 
         reactions specified by the user
      */
/*===== end of Modules  ======*/

#define CLEAN

/*===== Options ========*/
/* #define SHOWPLOTS */
     /* Display  graphical plots on the screen */ 

/*===== End of DEFINE  settings ===== */


#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"
#include"lib/pmodel.h"


int main(int argc,char** argv)
{
  int err,nw;
   char cdmName[10];
   int spin2, charge3,cdim;
   double laMax;   

  ForceUG=0;  /* to Force Unitary Gauge assign 1 */
//   VWdecay=0; VZdecay=0;
   
#ifdef SUGRA
{
  double m0,mhf,a0,tb;
  double Lambda, aLambda,aKappa,sgn;
  double mXiF=0,mXiS=0,muP=0,msP=0,m3h=0;

  if(argc<7) 
  { 

    printf(" This program needs 6 parameters:\n"
           "   m0      common scalar mass at GUT scale\n"
           "   mhf     common gaugino mass at GUT scale\n"
           "   a0      trilinear soft breaking parameter at GUT scale\n"
           "   tb      tan(beta) \n"
           "   Lambda   Lambda parameter at SUSY\n"
           "   aKappa  aKappa parameter at GUT\n"
           );
    printf(" Auxiliary parameters are:\n"
           "   sgn     +/-1,  sign of Higgsino mass term (default 1)\n" 
           "   aLambda at GUT (default aLambda=a0)\n"    
           "   Mtp     top quark pole mass\n"
           "   MbMb    Mb(Mb) scale independent b-quark mass\n"
           "   alfSMZ  strong coupling at MZ\n");
    printf("Example:  ./main 135 600 -1300 2 0.5 -1400\n");
      exit(1); 
  } else  
  {  double Mtp,MbMb,alfSMZ;
     sscanf(argv[1],"%lf",&m0);
     sscanf(argv[2],"%lf",&mhf);
     sscanf(argv[3],"%lf",&a0);
     sscanf(argv[4],"%lf",&tb);
     sscanf(argv[5],"%lf",&Lambda);
     sscanf(argv[6],"%lf",&aKappa);
     if(argc>7)  sscanf(argv[7],"%lf",&sgn); else sgn=1;
     if(argc>8)  sscanf(argv[8],"%lf",&aLambda); else aLambda=a0;

     if(argc>9){ sscanf(argv[9],"%lf",&Mtp);    assignValW("Mtp",Mtp);      }
     if(argc>10){ sscanf(argv[10],"%lf",&MbMb);   assignValW("MbMb",MbMb);    }
     if(argc>11){ sscanf(argv[11],"%lf",&alfSMZ); assignValW("alfSMZ",alfSMZ);}
  }

  err=nmssmSUGRA( m0,mhf, a0,tb, sgn,  Lambda, aLambda, aKappa,mXiF,mXiS,muP,msP,m3h);
}
#elif defined(EWSB)
{
  if(argc!=2)
  { 
      printf(" Correct usage:  ./main <file with NMSSM parameters> \n");
      printf(" Example      :  ./main  data1.par \n");
      exit(1);
  }

  err=readVar(argv[1]);
  
  if(err==-1)     {printf("Can not open the file\n"); exit(1);}
  else if(err>0)  { printf("Wrong file contents at line %d\n",err);exit(1);}          
  err=nmssmEWSB();  
}
#else
{
   int mode =0;
            
   
   printf("\n========= SLHA file input =========\n");
   

   if(argc <2) 
   {  printf("The program needs one argument:the name of SLHA input file.\n"
            "Example: ./main spectr.dat \n");
      exit(1);
   }  
   
   printf("Initial file  \"%s\"\n",argv[1]);
   

/*
  mode= 1*m1+2*m2+4*m4+8*m8+16*m16
  
  m1  0 overwrite all;  1 keep old data
  m2  0 ignore mistake  1: stop in case of mistake in input file.
  m4  0 read DECAY      1: don't read   Decay 
  m8  0 read BLOCK      1: don't read   Blocks
  m16 0 read QNUMBERS   1: don't read   QNUMBERS 
*/
     mode=4;
   err=readSLHA(argv[1],mode);
     
   if(err) exit(2);
}

#endif
    slhaWarnings(stdout);
  if(err) exit(1);

// to load decay tables produced by  NMSSMTools uncomment next line
//   slhaRead("decay",1);
  err=sortOddParticles(cdmName);
  if(err) { printf("Can't calculate %s\n",cdmName); return 1;}

  qNumbers(cdmName,&spin2, &charge3, &cdim);
  printf("\nDark matter candidate is '%s' with spin=%d/2\n",
  cdmName,       spin2); 
  if(charge3) { printf("Dark Matter has electric charge %d/3\n",charge3); exit(1);}
  if(cdim!=1) { printf("Dark Matter is a color particle\n"); exit(1);}
  if(strcmp(cdmName,"~o1")) printf(" ~o1 is not CDM\n"); 
                    else o1Contents(stdout);

/*  printVar(stdout);  */

 

#ifdef MASSES_INFO
{
  printf("\n=== MASSES OF HIGGS AND SUSY PARTICLES: ===\n");
  printHiggs(stdout);
  printMasses(stdout,1);  
}
#endif

#ifdef CONSTRAINTS
{ double constr0,constrM, constrP,csLim;

  printf("\n\n==== Physical Constraints: =====\n");

  constr0=bsgnlo(&constrM,&constrP);
  printf("B->s,gamma = %.2E (%.2E ,  %.2E  ) \n",constr0,constrM, constrP );

  constr0= bsmumu(&constrM,&constrP);
  printf("Bs->mu,mu  = %.2E (%.2E ,  %.2E  ) \n",constr0,constrM, constrP );
  
  constr0=btaunu(&constrM,&constrP);
  printf("B+->tau+,nu= %.2E (%.2E ,  %.2E  ) \n",constr0, constrM, constrP );
  
  constr0=deltaMd(&constrM,&constrP);
  printf("deltaMd    = %.2E (%.2E ,  %.2E  ) \n",constr0,constrM, constrP );

  constr0=deltaMs(&constrM,&constrP);
  printf("deltaMs    = %.2E (%.2E ,  %.2E  ) \n",constr0,constrM, constrP );

  constr0=gmuon(&constrM,&constrP);
  printf("(g-2)/BSM = %.2E (%.2E ,  %.2E  ) \n",constr0,constrM, constrP );

  if(Zinvisible()) printf("Excluded by Z->invisible\n");
  if(LspNlsp_LEP(&csLim)) printf("LEP excluded by e+,e- -> DM q qbar cross section=%.2E pb \n",csLim);
     
}
#endif


#if defined(HIGGSBOUNDS) || defined(HIGGSSIGNALS)
{  int NH0=5,NHch=1;

   int HB_id[3],HB_result[3];
   double  HB_obsratio[3],HS_observ,HS_chi2, HS_pval;
   char HB_chan[3][100]={""}, HB_version[50], HS_version[50]; 

   system("cat spectr decay > HB.in");
//   NH0=hbBlocksMO("HB.in",&NHch);
   system("echo 'BLOCK DMASS\n 25  2\n 35  2\n 45 2\n 36 2\n 46 2\n'>> HB.in");   
#include "../include/hBandS.inc"
#ifdef HIGGSBOUNDS
   printf("  HB(%s)\n", HB_version);
   for(int i=0;i<3;i++) printf("  id= %d  result = %d  obsratio=%.2E  channel= %s \n", HB_id[i],HB_result[i],HB_obsratio[i],HB_chan[i]);
#endif
#ifdef HIGGSSIGNALS
   printf("  HS(%s)\n",HS_version); 
   printf(" Nobservables=%.0f chi^2 = %.2E pval= %.2E\n",HS_observ,HS_chi2, HS_pval);
#endif
}   
#endif

#ifdef SUPERISO
{
  int err= callSuperIsoSLHA();
  if(err==0)
  { printf("\nSuperIso Flavour MSSM and ( SM)  observables :\n");
    printf("  BR(b->s gamma)                     %.3E  (%.3E)\n", slhaValFormat("FOBS",0.,"    5    1  %lf    0     2     3    22        "), slhaValFormat("FOBSSM",0.,"    5    1  %lf    0     2     3    22        ")); 
    printf("  Delta0(B->K* gamma)                %.3E  (%.3E)\n", slhaValFormat("FOBS",0.,"  521    4  %lf    0     2   313    22        "), slhaValFormat("FOBSSM",0.,"  521    4  %lf    0     2   313    22        ")); 
    printf("  BR(B_s->mu+ mu-)                   %.3E  (%.3E)\n", slhaValFormat("FOBS",0.,"  531    1  %lf    0     2    13   -13        "), slhaValFormat("FOBSSM",0.,"  531    1  %lf    0     2    13   -13        ")); 
    printf("  BR(B_u->tau nu)                    %.3E  (%.3E)\n", slhaValFormat("FOBS",0.,"  521    1  %lf    0     2   -15    16        "), slhaValFormat("FOBSSM",0.,"  521    1  %lf    0     2   -15    16        ")); 
    printf("  R(B_u->tau nu)                     %.3E  (%.3E)\n", slhaValFormat("FOBS",0.,"  521    2  %lf    0     2   -15    16        "), slhaValFormat("FOBSSM",0.,"  521    2  %lf    0     2   -15    16        ")); 
    printf("  BR(D_s->tau nu)                    %.3E  (%.3E)\n", slhaValFormat("FOBS",0.,"  431    1  %lf    0     2   -15    16        "), slhaValFormat("FOBSSM",0.,"  431    1  %lf    0     2   -15    16        ")); 
    printf("  BR(D_s->mu nu)                     %.3E  (%.3E)\n", slhaValFormat("FOBS",0.,"  431    1  %lf    0     2   -13    14        "), slhaValFormat("FOBSSM",0.,"  431    1  %lf    0     2   -13    14        ")); 
    printf("  BR(B+->D0 tau nu)                  %.3E  (%.3E)\n", slhaValFormat("FOBS",0.,"  521    1  %lf    0     3   421   -15    16  "), slhaValFormat("FOBSSM",0.,"  521    1  %lf    0     3   421   -15    16  ")); 
    printf("  BR(B+->D0 tau nu)/BR(B+-> D0 e nu) %.3E  (%.3E)\n", slhaValFormat("FOBS",0.,"  521   11  %lf    0     3   421   -15    16  "), slhaValFormat("FOBSSM",0.,"  521   11  %lf    0     3   421   -15    16  ")); 
    printf("  BR(K->mu nu)/BR(pi->mu nu)         %.3E  (%.3E)\n", slhaValFormat("FOBS",0.,"  321   11  %lf    0     2   -13    14        "), slhaValFormat("FOBSSM",0.,"  321   11  %lf    0     2   -13    14        ")); 
    printf("  R_mu23                             %.3E  (%.3E)\n", slhaValFormat("FOBS",0.,"  321   12  %lf    0     2   -13    14        "), slhaValFormat("FOBSSM",0.,"  321   12  %lf    0     2   -13    14        ")); 
  }
}
#endif

#ifdef LILITH
{  double m2logL, m2logL_reference=0,pvalue;
   int exp_ndf,n_par=0,ndf;
   char call_lilith[100], Lilith_version[20];
//LilithMO("Lilith_in.xml");
   if(LilithMDL("Lilith_in.xml"))
   {        
#include "../include/Lilith.inc"
      if(ndf)
      {
        printf("LILITH(DB%s):  -2*log(L): %.2f; -2*log(L_reference): %.2f; ndf: %d; p-value: %.2E \n", 
        Lilith_version,m2logL,m2logL_reference,ndf,pvalue);
      }  
   } else printf("LILITH: there is no Higgs candidate\n");
}     
#endif


#ifdef SMODELS
{  int result=0;
   double Rvalue=0;
   char analysis[30]={},topology[30]={};
   int LHCrun=LHC8|LHC13;  //  LHC8  - 8TeV; LHC13  - 13TeV; 
#include "../include/SMODELS.inc" 
}   
#endif 



#ifdef OMEGA
{ int fast=1;
  double Beps=1.E-4, cut=0.01;
  double Omega,Xf;   
  printf("\n==== Calculation of relic density =====\n");
// to exclude processes with virtual W/Z in DM   annihilation
  VWdecay=1; VZdecay=0; cleanDecayTable();  
  Omega=darkOmega(&Xf,fast,Beps,&err);
// Omega=darkOmega2(fast,Beps);
  
  printf("Xf=%.2e Omega=%.2e\n",Xf,Omega);
  if(Omega>0)printChannels(Xf,cut,Beps,1,stdout);
  
// to restore default switches  
  VZdecay=1; VWdecay=1; cleanDecayTable();  
}
#endif

#ifdef INDIRECT_DETECTION
{ 
  int err,i;
  double Emin=1,/* Energy cut  in GeV   */  sigmaV;
  double vcs_gz,vcs_gg;
  char txt[100];
  double SpA[NZ],SpE[NZ],SpP[NZ];
  double FluxA[NZ],FluxE[NZ],FluxP[NZ];
//  double * SpNe=NULL,*SpNm=NULL,*SpNl=NULL;
double SpNe[NZ],SpNm[NZ],SpNl[NZ];
  double Etest=Mcdm/2;
  
printf("\n==== Indirect detection =======\n");  

  sigmaV=calcSpectrum(1+2+4,SpA,SpE,SpP,SpNe,SpNm,SpNl ,&err);
    /* Returns sigma*v in cm^3/sec.     SpX - calculated spectra of annihilation.
       Use SpectdNdE(E, SpX) to calculate energy distribution in  1/GeV units.
       
       First parameter 1-includes W/Z polarization
                       2-includes gammas for 2->2+gamma
                       4-print cross sections             
    */
  printf("sigmav=%.2E[cm^3/s]\n",sigmaV);  

  { 
     double fi=0.1,dfi=0.05; /* angle of sight and 1/2 of cone angle in [rad] */ 

     gammaFluxTab(fi,dfi, sigmaV, SpA,  FluxA);
     printf("Photon flux  for angle of sight f=%.2f[rad]\n"
     "and spherical region described by cone with angle %.2f[rad]\n",fi,2*dfi);
     
#ifdef SHOWPLOTS
     sprintf(txt,"Photon flux[cm^2 s GeV]^{1} at f=%.2f[rad], cone angle %.2f[rad]",fi,2*dfi);
     displayPlot(txt,"E[GeV]",Emin,Mcdm,0,1,"flux",0,SpectdNdE,FluxA);
#endif
     printf("Photon flux = %.2E[cm^2 s GeV]^{-1} for E=%.1f[GeV]\n",SpectdNdE(Etest,FluxA), Etest);       
  }

  { 
    posiFluxTab(Emin, sigmaV, SpE,  FluxE);
#ifdef SHOWPLOTS     
    displayPlot("positron flux [cm^2 s sr GeV]^{-1}","E[GeV]",Emin,Mcdm,0,1,"flux",0,SpectdNdE,FluxE);
#endif
    printf("Positron flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxE),  Etest);           
  }
  
  { 
    pbarFluxTab(Emin, sigmaV, SpP, FluxP  ); 
#ifdef SHOWPLOTS    
     displayPlot("antiproton flux [cm^2 s sr GeV]^{-1}" ,"E[GeV]",Emin,Mcdm,0,1,"flux",0,SpectdNdE,FluxP);
#endif
    printf("Antiproton flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxP),  Etest);             
  }
}  
#endif

#ifdef LoopGAMMA  
  { double vcs_gg,vcs_gz;
    if(loopGamma(&vcs_gg, &vcs_gz )==0)
    {
      printf("Gamma  ray lines:\n");
      printf("E=%.2E[GeV]  vcs(Z,A)= %.2E[cm^3/s]\n",Mcdm-91.19*91.19/4/Mcdm,vcs_gz);  
      printf("E=%.2E[GeV]  vcs(A,A)= %.2E[cm^3/s]\n",Mcdm,vcs_gg);
    }
  }
#endif       



#ifdef RESET_FORMFACTORS
{
/* 
   The user has approach to form factors  which specifies quark contents 
   of  proton and nucleon via global parametes like
      <Type>FF<Nucleon><q>
   where <Type> can be "Scalar", "pVector", and "Sigma"; 
         <Nucleon>     "P" or "N" for proton and neutron
         <q>            "d", "u","s"

   calcScalarQuarkFF( Mu/Md, Ms/Md, sigmaPiN[MeV], sigmaS[MeV])  
   calculates and rewrites Scalar form factors
*/

  printf("protonFF (default) d %E, u %E, s %E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);                               
  printf("neutronFF(default) d %E, u %E, s %E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);

//                    To restore default form factors of  version 2  call 
     calcScalarQuarkFF(0.553,18.9,55.,243.5);

  printf("protonFF (new)     d %E, u %E, s %E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);                               
  printf("neutronFF(new)     d %E, u %E, s %E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);

//                    To restore default form factors  current version  call 
//  calcScalarQuarkFF(0.56,20.2,34,42);




/* Option to change parameters of DM velocity  distribution  */   
   SetfMaxwell(220.,600.);
/* 
    dN  ~  exp(-v^2/arg1^2)*Theta(v-arg2)  d^3v     
    Earth velocity with respect to Galaxy defined by 'Vearth' parameter.
    All parameters are  in [km/s] units.       
*/
}
#endif

#ifdef CDM_NUCLEON
{ double pA0[2],pA5[2],nA0[2],nA5[2];
  double Nmass=0.939; /*nucleon mass*/
  double SCcoeff;        

printf("\n==== Calculation of CDM-nucleons amplitudes  =====\n");   

    nucleonAmplitudes(CDM1, pA0,pA5,nA0,nA5);
    printf("CDM-nucleon micrOMEGAs amplitudes:\n");
    printf("proton:  SI  %.3E  SD  %.3E\n",pA0[0],pA5[0]);
    printf("neutron: SI  %.3E  SD  %.3E\n",nA0[0],nA5[0]); 

  SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm/(Nmass+ Mcdm),2.);
    printf("CDM-nucleon cross sections[pb]:\n");
    printf(" proton  SI %.3E  SD %.3E\n",SCcoeff*pA0[0]*pA0[0],3*SCcoeff*pA5[0]*pA5[0]);
    printf(" neutron SI %.3E  SD %.3E\n",SCcoeff*nA0[0]*nA0[0],3*SCcoeff*nA5[0]*nA5[0]);

}
#endif
  
#ifdef CDM_NUCLEUS
{ double dNdE[300];
  double nEvents;

printf("\n======== Direct Detection ========\n");    

  nEvents=nucleusRecoil(Maxwell,73,Z_Ge,J_Ge73,SxxGe73,dNdE);

  printf("73Ge: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));
                                                                                                         
#ifdef SHOWPLOTS
    displayPlot("Distribution of recoil energy of 73Ge","E[keV]",1,50,0,1,"dNdE",0,dNdERecoil,dNdE);
#endif

  nEvents=nucleusRecoil(Maxwell,131,Z_Xe,J_Xe131,SxxXe131,dNdE);

  printf("131Xe: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   
#ifdef SHOWPLOTS
    displayPlot("Distribution of recoil energy of 131Xe","E[keV]",1,50,0,1,"dNdE",0,dNdERecoil,dNdE);
#endif

  nEvents=nucleusRecoil(Maxwell,23,Z_Na,J_Na23,SxxNa23,dNdE);

  printf("23Na: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   
#ifdef SHOWPLOTS
    displayPlot("Distribution of recoil energy of 23Na","E[keV]",1,50,0,1,"dNdE",0,dNdERecoil,dNdE);
#endif

  nEvents=nucleusRecoil(Maxwell,127,Z_I,J_I127,SxxI127,dNdE);

  printf("I127: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   
#ifdef SHOWPLOTS
    displayPlot("Distribution of recoil energy of 127I","E[keV]",1,50,0,1,"dNdE",0,dNdERecoil,dNdE);
#endif
  
}
#endif 

#ifdef NEUTRINO
{ double nu[NZ], nu_bar[NZ],mu[NZ];
  int forSun=1;
  double Emin=1;
  
 printf("\n===============Neutrino Telescope=======  for  "); 
 if(forSun) printf("Sun\n"); else printf("Earth\n");  

  err=neutrinoFlux(Maxwell,forSun, nu,nu_bar);
  if(err==0)
  {
#ifdef SHOWPLOTS
    displayPlot("neutrino fluxes [1/Year/km^2/GeV]","E[GeV]",Emin,Mcdm,0,2,"nu",0,SpectdNdE,nu, "nu_bar",0,SpectdNdE,nu_bar);
#endif

    printf(" E>%.1E GeV neutrino/anti-neutrino fluxes   %.2E/%.2E [1/Year/km^2]\n",Emin,
          spectrInfo(Emin,nu,NULL), spectrInfo(Emin,nu_bar,NULL));
 
//ICE CUBE  
    if(forSun) printf("IceCube22 exclusion confidence level = %.2E%%\n", 100*exLevIC22(nu,nu_bar,NULL));
  
/* Upward events */
 
    Emin=1;  
    muonUpward(nu,nu_bar, mu);
#ifdef SHOWPLOTS  
    displayPlot("Upward muons[1/Year/km^2/GeV]","E[GeV]",Emin,Mcdm/2,0,1,"flux",0,SpectdNdE,mu);
#endif
    printf(" E>%.1E GeV Upward muon flux    %.2E [1/Year/km^2]\n",Emin,spectrInfo(Emin,mu,NULL));
   
/* Contained events */
    muonContained(nu,nu_bar,1., mu);
#ifdef SHOWPLOTS  
    displayPlot("Contained  muons[1/Year/km^3/GeV]","E[GeV]",Emin,Mcdm,0,1,"flux",0,SpectdNdE,mu); 
#endif
    printf(" E>%.1E GeV Contained muon flux %.2E [1/Year/km^3]\n",Emin,spectrInfo(Emin,mu,NULL));
  }
}        
#endif 


#ifdef DECAYS
{  
  txtList L;
   double width,br;
   char * pname;
   
   printf("\nParticle decays\n"); 
   pname = "h1";
    width=pWidth(pname,&L);
    printf("%s->  :   total width=%E \n and Branchings:\n",pname,width);
    printTxtList(L,stdout);

   pname = "h2";
    width=pWidth(pname,&L);
    printf("%s-> :   total width=%E \n and Branchings:\n",pname,width);
    printTxtList(L,stdout);

   pname = "~o2";
    width=pWidth(pname,&L);
    printf("%s-> :   total width=%E \n and Branchings:\n",pname,width);
    printTxtList(L,stdout);    
}
#endif

#ifdef CROSS_SECTIONS
{
  char* next,next_;
  double nextM;
    
  next=nextOdd(1,&nextM); 
  if(next && nextM<1000)  
  { 
     double cs, Pcm=6500, Qren, Qfact, pTmin=0;
     int nf=3;
     char*next_=antiParticle(next);
     Qren=Qfact=nextM; 
 
     printf("\npp > nextOdd  at sqrt(s)=%.2E GeV\n",2*Pcm);  
  
     Qren=Qfact;
     cs=hCollider(Pcm,1,nf,Qren, Qfact, next,next_,pTmin,1);
     printf("Production of 'next' odd particle: cs(pp-> %s,%s)=%.2E[pb]\n",next,next_, cs);
  }  
}
#endif

#ifdef CLEAN
//  system("rm -f inp decay spectr nngg.out omega");
  system("rm -f HB.* HS.* hb.* hs.*  debug_channels.txt debug_predratio.txt  Key.dat");
  system("rm -f Lilith_*   particles.py*");
  system("rm -f   smodels.in  smodels.log  smodels.out  summary.*");  

#endif

  return 0;
}
