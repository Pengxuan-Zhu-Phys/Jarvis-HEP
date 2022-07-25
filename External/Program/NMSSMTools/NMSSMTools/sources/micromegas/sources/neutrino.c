#include<stdio.h>
#include<stdlib.h>
#include"micromegas.h"
#include"micromegas_aux.h"
//#define PRINT

#define Gconst (0.7426E-30) /* m/g */
#define Vlight (299792458.) /* m/s */  
#define mp_g   (1.673E-24)  /* proton mass in grams */
#define mp_gev (0.939)   
#define GeVfm  (5.067731162)


/*====================  CAPTURE RATES ===============================*/


static double Rcm, *rTab,*tTab,*rhoTab, *phiTab, *aFraction;
static int nTab;


static double RS=6.96E10;   /* cm */

static double (*fvStat)(double)=NULL;
static double v_stat,MA, muX, FFalpha;
static char*CDM;
/*
double Tcapture(double T, double w, double v)
{ double mu=M_cdm/MA;
  double muP= (mu+1)/2;
  double muM= (mu-1)/2;
  double aP=sqrt(MA/(2*T))*(muP*v + muM*w);
  double aM=sqrt(MA/(2*T))*(muP*v - muM*w);  
  double bP=sqrt(MA/(2*T))*(muM*v + muP*w);
  double bM=sqrt(MA/(2*T))*(muM*v - muP*w);  
   

 return - T/(MA*mu*mu)*(  mu*(-aP*exp(-aM*aM) - aM*exp(-aP*aP))/sqrt(M_PI) 
                         +(mu/2 -mu*aP*aM - muP*muM)*( erf(aP)-erf(-aM)) + muP*muP*( erf(bP)-erf(-bM))*exp(-M_cdm/(2*T)*(v*v-w*w)));  

}
*/

static numout* cc23; 


static double r3Integrand(double r3)
{
  double  r,vescQ,rhoGeVcm3,mfrac,w2rate;
  r=pow(r3,0.333333333333333);
  vescQ=2*polint1(r/Rcm,nTab,rTab,phiTab);
  w2rate= vescQ-v_stat*v_stat*muX;
  if(w2rate<=0) return 0;
  rhoGeVcm3=polint1(r/Rcm,nTab,rTab,rhoTab)/mp_g;   // in proton mass units
  mfrac=polint1(r/Rcm,nTab,rTab,aFraction);
  if(FFalpha==0) return rhoGeVcm3*mfrac*w2rate; 
  else  return  rhoGeVcm3*mfrac*(1+muX)*( exp(-v_stat*v_stat*FFalpha) - exp(-(v_stat*v_stat+vescQ)*FFalpha/(1+muX)))/FFalpha;
}

static double vIntegrand(double v)
{ double u=fvStat(v); 
  double phiv,Rmax;
  if(!u)return 0;

  v_stat=v*1.E3/Vlight;  // km -> c units 
  phiv=muX*v_stat*v_stat/2;

  if(phiv> phiTab[0]) return 0;
  if(phiv <= phiTab[nTab-1]) Rmax=Rcm; else   Rmax=Rcm*polint1(phiv,nTab,phiTab,rTab); 
  return u*(4*M_PI/3.)*simpson(r3Integrand,0,Rmax*Rmax*Rmax,1.E-4,NULL)*Vlight*Vlight*1.E-1; // in cm
}



static double sIntegrand(double r3)
{
  double  r,rhoGeVcm3,mfrac;
  r=pow(r3,0.333333333333333);
  rhoGeVcm3=polint1(r/Rcm,nTab,rTab,rhoTab)/mp_g;   // in proton mass units
  mfrac=polint1(r/Rcm,nTab,rTab,aFraction);
  return rhoGeVcm3*mfrac;   
}

static void derivePhi(double r, double *x, double *dx)
{  double r2=r*r;
   dx[0]=4*M_PI*r2*polint1(r/Rcm,nTab,rTab,rhoTab);
   dx[1]=Gconst*100*x[0]/r2;
//printf(" r=%E dx[0]=%e dx[1]=%E\n", r,dx[0],dx[1]);    
}


static void fillGravPotential(void)
{ double x[2],phi0,r1,r0,m1;
  int i;
  r0=rTab[0]*Rcm;
  if(r0==0)
  {  phiTab[0]=0;
     r1=rTab[1]*Rcm;
     x[0]=4./3.*M_PI*rhoTab[0]*r1*r1*r1;
     phiTab[1]=x[1]=2./3.*M_PI*rhoTab[0]*r1*r1*Gconst*100;
     i=2;
  } else 
  { 
    x[0]=4./3.*M_PI*rhoTab[0]*r0*r0*r0;        
    phiTab[0]=x[1]=2./3.*M_PI*rhoTab[0]*r0*r0*Gconst*100;
//printf("phiTab[0]=%E\n",phiTab[0]);    
    i=1;
  }  
  for( ;i<nTab;i++)
  {
    odeint(x,2,Rcm*rTab[i-1], Rcm*rTab[i],1.E-3,Rcm*(rTab[i]-rTab[i-1])/3,derivePhi);
    phiTab[i]=x[1];       
  }    
  phi0=x[0]*Gconst*100/Rcm;
//printf("VescSurface=%E\n", sqrt(2*phi0)*Vlight); 
  for(i=0;i<nTab;i++) phiTab[i]= phiTab[nTab-1]-phiTab[i]+phi0;
//printf("VescCenter=%E\n", sqrt(2*phiTab[0])*Vlight);
}



#define SUNPOINTS 1268 
static double  rSun[SUNPOINTS], tSun[SUNPOINTS], rhoSun[SUNPOINTS], Hsun[SUNPOINTS],He4sun[SUNPOINTS],
He3sun[SUNPOINTS],C12sun[SUNPOINTS],N14sun[SUNPOINTS],O16sun[SUNPOINTS],phiSun[SUNPOINTS];


#define EARTHPOINTS 101

static double OEarth[EARTHPOINTS],NaEarth[EARTHPOINTS],MgEarth[EARTHPOINTS],AlEarth[EARTHPOINTS],
SiEarth[EARTHPOINTS],PEarth[EARTHPOINTS],SEarth[EARTHPOINTS],CaEarth[EARTHPOINTS],CrEarth[EARTHPOINTS],
FeEarth[EARTHPOINTS],NiEarth[EARTHPOINTS];

static double rEarth[EARTHPOINTS],tEarth[EARTHPOINTS],rhoEarth[EARTHPOINTS],phiEarth[EARTHPOINTS];


static int readSunData(void)
{
  static int rdOK=0;
  int i;
  if(!rdOK) 
  { FILE*f;
    char fname[300];
    sprintf(fname,"%s/sources/data_nu/SSM.dat",micrO);
    
    f=fopen(fname,"r"); 
    fscanf(f,"%*[^\n]"); 
    for(i=0;i<SUNPOINTS;i++)
    {
      fscanf(f,"%*lf %lf %lf %lf %*lf %*lf %lf %lf %lf %lf %lf %lf ",
      rSun+i, tSun+i, rhoSun+i, Hsun+i,He4sun+i, He3sun+i, C12sun+i,N14sun+i,O16sun+i);
//    printf("%E %E %E %E %E %E %E %E %E\n", rsun[i], tsun[i], rhosun[i], Hsun[i],He4sun[i], He3sun[i], C12sun[i],N14sun[i],O16sun[i]);
//      printf("frac= %E\n",(1-( Hsun[i]+He4sun[i]+ He3sun[i]+ C12sun[i]+N14sun[i]+O16sun[i]))/O16sun[i] ); 
    }
    fclose(f);
  }    
  rTab=rSun;
  tTab=tSun;
  rhoTab=rhoSun;
  phiTab=phiSun;
  nTab=SUNPOINTS;
  Rcm=RS; 

  if(!rdOK) {fillGravPotential(); rdOK=1;}
  return 0;
}

static double captureSun(double(*vfv)(double), double M_cdm, double pA0, double nA0, double pA5, double nA5)
{ 
  /*                 H    He      C     N     O    Ne     Mg    Si    S     Fe    Na    Al   Cl    Ar     Ca   Cr   Ni */
  int      A[17] ={  1,   4   , 12   ,14   ,   16,   20,   24,   28,   32,   56,   23,   27,   35,   40,   40,  52,  59  };
  int      Z[17] ={  1,   2   ,  6   , 7   ,    8,   10,   12,   14,   16,   26,   11,   13,   17,   18,   20,  24,  28  };
  double P10[17] ={ 12,  10.93,  8.39, 7.78, 8.66, 7.84, 7.53, 7.51, 7.14, 7.45, 6.17, 6.37, 5.50, 6.18, 6.31, 5.64, 6.23};  
//  double P10[17] ={ 12,  10.93,  8.39, 7.78, 8.66, 8.08, 7.58, 7.55, 7.33, 7.50 , 6.33, 6.47, 5.50, 6.40, 6.36, 5.64, 6.25};  
  double  sumI,csSDp,mu;
  double * sunFractions[5]={Hsun,He4sun, C12sun,N14sun,O16sun};
  int i;

  readSunData();
  fvStat=vfv;
  
  mu=M_cdm*mp_gev/(M_cdm+mp_gev); 
  csSDp=12/M_PI*(pA5*mu)*(pA5*mu)*3.8937966E8*1E-36;  // cm^2
  
  for(sumI=0,i=0;i<17;i++) 
  { double si,FF,fr, cI;
    double vmaxC,vmaxQ;
    MA=A[i]*mp_gev;
    mu=M_cdm*MA/(M_cdm+MA);
    si= (Z[i]*pA0+(A[i]-Z[i])*nA0)*mu;
    si=4/M_PI*si*si*3.8937966E8*1E-36;  // cm^2
    FF=1;
    if(i==0) si+=csSDp;
    if(i<4) {aFraction=sunFractions[i];fr=1;}  else { aFraction=O16sun; fr=A[i]/((double)A[4])*pow(10., P10[i]-P10[4]);}
    muX=(M_cdm-MA)*(M_cdm-MA)/(4*M_cdm*MA);
    vmaxC=(vEsc+vRot)/(1.E-3*Vlight);
    vmaxQ=vmaxC*vmaxC;
    if(vmaxQ*muX> 2*phiSun[0] ) vmaxQ=2*phiTab[0]/muX;

    if(i==0) FFalpha=0; else FFalpha=M_cdm*MA*pow((0.91*pow(MA,0.3333333) +0.3)*GeVfm,2)/3;
    cI=fr*FF/A[i]*si*simpson(vIntegrand,0,1.E-3*Vlight*sqrt(vmaxQ),1.E-3,NULL);     
    sumI+=cI;
//printf("C%d = %E  \n",A[i],rhoDM/M_cdm*cI);
  }
//printf("Sun capture =%E\n",rhoDM*sumI/M_cdm);  
  return  rhoDM*sumI/M_cdm;
}



static double sigmaSun(double M_cdm,double pA0, double nA0, double pA5, double nA5, double R3)
{ 
  double si,fr, cI,sumI,csSDp,mu;
  /*                 H    He      C     N     O    Ne     Mg    Si    S     Fe    Na    Al   Cl    Ar     Ca   Cr   Ni */
  int      A[17] ={  1,   4   , 12   ,14   ,   16,   20,   24,   28,   32,   56,   23,   27,   35,   40,   40,  52,  59  };
  int      Z[17] ={  1,   2   ,  6   , 7   ,    8,   10,   12,   14,   16,   26,   11,   13,   17,   18,   20,  24,  28  };
  double P10[17] ={ 12,  10.93,  8.39, 7.78, 8.66, 7.84, 7.53, 7.51, 7.14, 7.45, 6.17, 6.37, 5.50, 6.18, 6.31, 5.64, 6.23};  
  double * sunFractions[5]={Hsun,He4sun, C12sun,N14sun,O16sun};
  int i;
  
  readSunData();
     
  mu=M_cdm*mp_gev/(M_cdm+mp_gev); 
  csSDp=12/M_PI*(pA5*mu)*(pA5*mu)*3.8937966E8*1E-36;  // cm^2
  double rSunInt;
  for(sumI=0,i=0;i<5;i++) 
  {
    MA=A[i]*mp_gev;
    mu=M_cdm*MA/(M_cdm+MA);
    si= (Z[i]*pA0+(A[i]-Z[i])*nA0)*mu;
    si=4/M_PI*si*si*3.8937966E8*1E-36;  // cm^2
    if(i==0) si+=csSDp;
/*
    if(i<4) {aFraction=sunFractions[i];fr=1;}  else { aFraction=O16sun; fr=A[i]/((double)A[4])*pow(10., P10[i]-P10[4]);}
    cI=si/A[i]*(4./3*M_PI)*simpson(sIntegrand,0,R3,1.E-3,NULL);   
    sumI+=cI;
*/ 
    aFraction=sunFractions[i];
    rSunInt=simpson(sIntegrand,0,R3,1.E-3,NULL);
    sumI+=si/A[i]*(4./3*M_PI)*rSunInt;    
  }
  
//  cI/=si/A[i]*pow(10,P10[i]);
  
  for(i=5;i<17;i++) 
  { MA=A[i]*mp_gev;
    mu=M_cdm*MA/(M_cdm+MA);
    si= (Z[i]*pA0+(A[i]-Z[i])*nA0)*mu;
    si=4/M_PI*si*si*3.8937966E8*1E-36;  // cm^2
//    sumI+=cI*si/A[i]*pow(10,P10[i]);
    sumI+= si/A[i]*(4./3*M_PI)*rSunInt*pow(10,P10[i]-P10[4]);
  } 
  return  sumI;
}

#define RE 6378E5   /* cm */

static double rhoFun(double r) { return polint1(r,nTab,rTab,rhoTab);}

static int readEarthData(void)
{
  static int rdOK=0;
  int i;
  if(!rdOK) 
  { FILE*f;
    char fname[300];
    sprintf(fname,"%s/sources/data_nu/EarthModel.dat",micrO);
    f=fopen(fname,"r"); 
    fscanf(f,"%*[^\n]"); 
    for(i=0;i<EARTHPOINTS;i++)
    {
      fscanf(f,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
      rEarth+i, tEarth+i, rhoEarth+i, OEarth+i,  NaEarth+i,  MgEarth+i,   AlEarth+i,  SiEarth+i,  PEarth+i, 
      SEarth+i,   CaEarth+i,  CrEarth+i, FeEarth+i, NiEarth+i);
//    printf("%E %E %E %E %E %E %E %E %E\n", rsun[i], tsun[i], rhosun[i], Hsun[i],He4sun[i], He3sun[i], C12sun[i],N14sun[i],O16sun[i]);
//      printf("frac= %E\n",(1-( Hsun[i]+He4sun[i]+ He3sun[i]+ C12sun[i]+N14sun[i]+O16sun[i]))/O16sun[i] ); 
    }
    fclose(f);
  }   
  rTab= rEarth;
  tTab=tEarth;
  rhoTab=rhoEarth;
  phiTab=phiEarth;
  nTab=EARTHPOINTS;
  Rcm=RE; 
  if(!rdOK) {fillGravPotential(); rdOK=1;} 
  return 0;
}

//double fracFun(double r){ return polint1(r,nTab, rTab,aFraction);}


static double captureEarth(double(*vfv)(double),double M_cdm,  double pA0, double nA0, double pA5, double nA5)
{ 
/*            O  Na  Mg   Al  Si  P   S   Ca  Cr, Fe  Ni  Cu  */
  int A[12] ={16, 23 ,24, 27, 28, 30, 32, 40, 52, 56, 58, 64};
  int Z[12] ={ 8, 11 ,12, 13, 14, 15, 16, 20, 24, 26, 28, 29};

  double  sumI,mu;
  double * earthFractions[11]={OEarth,NaEarth,MgEarth,AlEarth,SiEarth,PEarth,SEarth,CaEarth,CrEarth,FeEarth,NiEarth };
  int i;
  
  readEarthData();
  fvStat=vfv;
  
  for(sumI=0,i=0;i<11;i++) 
  { double si,cI;
    double vmaxC,vmaxQ;
    MA=A[i]*mp_gev;
    aFraction=earthFractions[i]; 
    mu=M_cdm*MA/(M_cdm+MA);
    si= (Z[i]*pA0+(A[i]-Z[i])*nA0)*mu;
    si=4/M_PI*si*si*3.8937966E8*1E-36;  // cm^2
    muX=(M_cdm-MA)*(M_cdm-MA)/(4*M_cdm*MA);
    vmaxC=(vEsc+vRot)/(1.E-3*Vlight);
    vmaxQ=vmaxC*vmaxC;
    if(vmaxQ*muX> 2*phiTab[0] ) vmaxQ=2*phiTab[0]/muX;
    FFalpha=M_cdm*MA*pow((0.91*pow(MA,0.3333333) +0.3)*GeVfm,2)/3; 
    cI=si/A[i]*simpson(vIntegrand,0,1.E-3*Vlight*sqrt(vmaxQ),1.E-3,NULL);     
    sumI+=cI;
//printf("C%d = %E  \n",A[i],rhoDM/M_cdm*cI);
  }
  return  rhoDM*sumI/M_cdm;
}


static double sigmaEarth(double M_cdm,double pA0, double nA0, double pA5, double nA5, double R3)
{ 
  double si,cI,sumI,mu;
/*               O  Na  Mg   Al  Si  P   S   Ca  Cr, Fe  Ni  Cu  */
  double A[12] ={16, 23 ,24, 27, 28, 30, 32, 40, 52, 56, 58, 64};
  double Z[12] ={ 8, 11 ,12, 13, 14, 15, 16, 20, 24, 26, 28, 29};

  double*earthFractions[11]={OEarth,NaEarth,MgEarth,AlEarth,SiEarth,PEarth,SEarth,CaEarth,CrEarth,FeEarth,NiEarth };
  int i;
  
  readEarthData();
   
  for(sumI=0,i=0;i<11;i++) 
  {
    MA=A[i]*mp_gev;
    mu=M_cdm*MA/(M_cdm+MA);
    si= (Z[i]*pA0+(A[i]-Z[i])*nA0)*mu;
    si=4/M_PI*si*si*3.8937966E8*1E-36;  // cm^2
    aFraction=earthFractions[i];
    cI=si/A[i]*(4./3*M_PI)*simpson(sIntegrand,0,R3,1.E-3,NULL);   
    sumI+=cI;
  }  
  return  sumI;
}

double captureCS(double(*vfv)(double),int forSun, double M_cdm, double csIp,double csIn,double csDp, double csDn)
{ 
  double pA0, nA0, pA5, nA5;
  double MN=0.939;
  double Mr=MN*M_cdm/(MN+M_cdm);
  double sCoeff= 2*sqrt(3.8937966E8/M_PI)*Mr;
           
  pA0=sqrt(fabs(csIp))/sCoeff;   if(csIp<0) pA0*=-1;
  nA0=sqrt(fabs(csIn))/sCoeff;   if(csIn<0) nA0*=-1;
  pA5=sqrt(fabs(csDp/3))/sCoeff; if(csDp<0) pA5*=-1;
  nA5=sqrt(fabs(csDn/3))/sCoeff; if(csDn<0) nA5*=-1;
  
  if(forSun) return captureSun(vfv,M_cdm,pA0,nA0,pA5,nA5);
       else  return captureEarth(vfv,M_cdm,pA0,nA0,pA5,nA5);               
}


/* ================ Basic Spectra  =====  */


typedef struct
{ int nX,nM,nCh;
  double *x; 
  double *m;
  int *ch;
  float *data;
// data[ iCh+nCh*(iX+nX*iM)]
  double*buff;        
}  nuTableStr;
  
static   void cleanNuTab(nuTableStr*t)
  { if(t) { if(t->x) free(t->x);
            if(t->m) free(t->m);
            if(t->ch)free(t->ch);
            if(t->buff)free(t->buff);
            if(t->data) free(t->data);
            free(t);
           }              
  }          


static nuTableStr*readNuTable( FILE*F)
{  
  
  nuTableStr*T=malloc(sizeof(nuTableStr));
  T->x=T->m=T->buff=NULL; T->ch=NULL; T->data=NULL;
  T->nX=T->nM=T->nCh=0;
  
  fscanf(F,"%*[^\n]\n");
  int posLn=ftell(F);
  fseek(F,0,SEEK_SET);
  fscanf(F,"%*s %*s");

  while(ftell(F)<posLn) 
  { T->nCh++;
    T->ch=realloc(T->ch,T->nCh*(sizeof(int)));
    if(1!=fscanf(F," %d ",T->ch+T->nCh-1))
    { printf("Format Error: file position %d\n",ftell(F));  cleanNuTab(T); return NULL;} 
  }
  

  int id=0,ix=0;
  for(int ix=0;;)
  {  
    double m,x;

    if(2!=fscanf(F," %lf %lf ", &m,&x)) break;
//printf("m=%e x=%e\n",m,x);    
    if(!T->m) 
    {  T->m=malloc(sizeof(double)); T->nM=1; T->m[0]=m;
//printf("T->m[0]=%E\n", T->m[0]);    
       T->x=malloc(sizeof(double)); T->nX=1; T->x[0]=x;
       ix=0;
    }
    else if(m!=T->m[T->nM-1])
    { 
      T->nM++;
      T->m=realloc(T->m,T->nM*sizeof(double));
      T->m[T->nM-1]=m;
      if(ix!=T->nX-1) { printf("Format Error: file position %d\n",ftell(F));cleanNuTab(T); return NULL;}
      ix=0;
    } else
    {
      if(T->nM==1)
      { T->nX++;
        T->x=realloc(T->x,T->nX*(sizeof(double)));
        T->x[T->nX-1]=x;
      }
      ix++;    
    }
    if(x!=T->x[ix]) { printf("Format Error: file position %d\n", ftell(F)); cleanNuTab(T); return NULL;} 
    
    T->data=realloc(T->data,(id+T->nCh)*sizeof(float));
    for(int i=0;i<T->nCh;i++) fscanf(F," %f ", T->data+(id++));
  }
  T->buff=malloc(T->nX*sizeof(double));
/*  
  { int i,j,k;
    id=0;
    for(i=0;i<T->nM;i++) for(j=0;j<T->nX;j++)
    { printf(" %.2E %.2E | ", T->m[i],T->x[j]);
      for(k=0;k<T->nCh;k++) printf(" %.2E(%.2E)",T->data[id++],T->data[k+T->nCh*(j+T->nX*i)] );
      printf("\n");
    }    
  }   
*/  
  return T;  
}
  

static int  nuTabInterpolation(nuTableStr*T,double Mass,int pdgN,int pol, int smooth, double*tab)
{  
  int i;
  
  tab[0]=Mass;
  for(i=1;i<NZ;i++) tab[i]=0;
  if(!T) return -1;
  
  int iCh,iM;
  for(iCh=0;iCh<T->nCh;iCh++) if(T->ch[iCh]==abs(pdgN)) break;
  
  if(iCh==T->nCh) { return -2;}

  int err=0;

  for(iM=0;iM<T->nM && Mass>T->m[iM];iM++);
  if(iM>0)iM--;
  if(iM>T->nM-2)iM=T->nM-2;
{  
  double a[4],m[4], aL[3],aR[4];
  int j1=0,j2=4;   
  
  double mm=log(Mass);
  if(iM>0 && iM<T->nM-2)
  { 
    double alpha; 
    j1=0;j2=4;
  
    m[0]=log(T->m[iM-1]);
    m[1]=log(T->m[iM]);
    m[2]=log(T->m[iM+1]);
    m[3]=log(T->m[iM+2]);
    alpha=(m[2]-mm)/(m[2]-m[1]);
    
    aL[0]=           (mm-m[1])*(mm-m[2])/            (m[0]-m[1])/(m[0]-m[2]);
    aL[1]= (mm-m[0])*          (mm-m[2])/(m[1]-m[0])/            (m[1]-m[2]);
    aL[2]= (mm-m[0])*(mm-m[1])          /(m[2]-m[0])/(m[2]-m[1])           ;

    aR[1]=           (mm-m[2])*(mm-m[3])/            (m[1]-m[2])/(m[1]-m[3]);
    aR[2]= (mm-m[1])*          (mm-m[3])/(m[2]-m[1])/            (m[2]-m[3]);
    aR[3]= (mm-m[1])*(mm-m[2])          /(m[3]-m[1])/(m[3]-m[2]); 
    
    a[0]=alpha*aL[0];
    a[1]=alpha*aL[1]+(1-alpha)*aR[1];
    a[2]=alpha*aL[2]+(1-alpha)*aR[2];
    a[3]=(1-alpha)*aR[3]; 

/*  
  a[0]=           (mm-m[1])*(mm-m[2])*(mm-m[3])/            (m[0]-m[1])/(m[0]-m[2])/(m[0]-m[3]); 
  a[1]= (mm-m[0])*          (mm-m[2])*(mm-m[3])/(m[1]-m[0])/            (m[1]-m[2])/(m[1]-m[3]);
  a[2]= (mm-m[0])*(mm-m[1])*          (mm-m[3])/(m[2]-m[0])/(m[2]-m[1])/            (m[2]-m[3]);
  a[3]= (mm-m[0])*(mm-m[1])*(mm-m[2])          /(m[3]-m[0])/(m[3]-m[1])/(m[3]-m[2]); 
*/     
  } else if(iM==0)
  { j1=1;j2=4;
    m[1]=log(T->m[iM]);  
    m[2]=log(T->m[iM+1]);
    m[3]=log(T->m[iM+2]);
         
    a[1]=           (mm-m[2])*(mm-m[3])/            (m[1]-m[2])/(m[1]-m[3]);
    a[2]= (mm-m[1])*          (mm-m[3])/(m[2]-m[1])/            (m[2]-m[3]);
    a[3]= (mm-m[1])*(mm-m[2])          /(m[3]-m[1])/(m[3]-m[2]); 
  } else             
  { 
    j1=0;j2=3;
    m[0]=log(T->m[iM-1]);
    m[1]=log(T->m[iM]);  
    m[2]=log(T->m[iM+1]);
    a[0]=           (mm-m[1])*(mm-m[2])/            (m[0]-m[1])/(m[0]-m[2]);
    a[1]= (mm-m[0])*          (mm-m[2])/(m[1]-m[0])/            (m[1]-m[2]);
    a[2]= (mm-m[0])*(mm-m[1])          /(m[2]-m[0])/(m[2]-m[1]); 
  } 
  
  for(i=0;i<T->nX;i++) T->buff[i]=0;

  int poleData= iCh<T->nCh-1 && T->ch[iCh]==T->ch[iCh+1]; 
  
  if(abs(pol)>1 || (pol!=0 && !poleData)){ err=2; pol=0;}
 
  if(pol==0 && poleData)  for(i=0;i<T->nX;i++)
  { 
    for(int j=j1;j<j2;j++)  T->buff[i]+= a[j]*(T->data[iCh+T->nCh*(i+T->nX*(iM+j-1))] + T->data[1+iCh+T->nCh*(i+T->nX*(iM+j-1))])/2;
  }     
  else
  { 
    if(pol==1) iCh++;
    for(i=0;i<T->nX;i++) for(int j=j1;j<j2;j++)  T->buff[i]+= a[j]*T->data[iCh+T->nCh*(i+T->nX*(iM+j-1))];
  }
}                                    

//printf("alpha=%E  T->ch[%d]=%d\n", alpha,iCh,T->ch[iCh]  );  
//  for(i=0;i<T->nX;i++) printf("x=%e  Y=%e\n", T->x[i], T->buff[i]);
  for(i=1;i<NZ;i++)
  { double x=exp(Zi(i));
    tab[i]=smooth?  x*polint3(x,T->nX, T->x ,T->buff): x*polint1(x,T->nX, T->x ,T->buff);
  }
  
  
  if(Mass<T->m[0] || Mass>T->m[T->nM-1]) err++;             
                                    
  return err;
}



int basicNuSpectra(int forSun, double Mass, int pdgN, int pol, double*nu, double*nuB)
{
  static nuTableStr *allTab[3][2][2]={{{NULL,NULL},{NULL,NULL}},  {{NULL,NULL},{NULL,NULL}}, {{NULL,NULL},{NULL,NULL}}};
  nuTableStr *currentTab;
  if(nu)  {nu[0]=Mass;  for(int i=1;i<NZ;i++) nu[i]=0; }
  if(nuB) {nuB[0]=Mass; for(int i=1;i<NZ;i++) nuB[i]=0;}
  if(!pdgN) return -2;
  pdgN=abs(pdgN);       
  if(forSun)  forSun=1;
  int outN,err;
  int pkgN;
  switch(WIMPSIM)
  { case 1: pkgN=2;
     if(pdgN<4) pdgN=1;
     break;
    case 0:
     if(forSun)
     {
       pkgN=1; 
       if(pdgN<4) pdgN=1;
       break;
     } else printf("PPPC4DM$nu does not contain spectra of neutrinos produced in the center of Earth.\n"
                   "DM nu stectra are used instead\n");   
    default:pkgN=0;
    if( pdgN==11 || pdgN==13) return 0;
    if(pdgN==21 || pdgN<4) pdgN=1;
    else if(pdgN==12 || pdgN==14 || pdgN==16)
    { if(forSun) pdgN=14; 
      else switch(pdgN)
      { case 12: case 16: return 0;
        case 14: if(nu)  nu[1] =2/(Zi(1)-Zi(2));
                 if(nuB) nuB[1]=2/(Zi(1)-Zi(2));  
                 return 0;
      }
    }       
  }
  for(outN=0;outN<2;outN++)
  if(!allTab[pkgN][forSun][outN])
  {
     char * packages[3]={"EVOL","MC","WIMPSIM"};
     char * source[2]={"Earth","Sun"};
     char * yeald[2] ={"numub","numu"};
     char * fname=malloc(50+strlen(micrO));    
     sprintf(fname,"%s/sources/data_nu/%s_%s_%s.txt",micrO,source[forSun],yeald[outN], packages[pkgN]);
     FILE *F=fopen(fname,"r");
     if(!F)
     { printf("Can not open file %s\n",fname); free(fname); return -1;}      
     allTab[pkgN][forSun][outN]=readNuTable(F);
     fclose(F);
     
     if(!allTab[pkgN][forSun][outN]) 
     { printf("Wrong file format: %s\n",fname);
       free(fname);
       return -2;
     }
     free(fname);  
  }
  { 
    nuTableStr*T=allTab[pkgN][forSun][0];
    if(Mass<T->m[0] || Mass>T->m[T->nM-1]) return 1;
  }  
  err=0;
  if(nu)  err=  nuTabInterpolation(allTab[pkgN][forSun][1],Mass,pdgN, pol , pkgN!=2, nu);
  if(nuB) err=  nuTabInterpolation(allTab[pkgN][forSun][0],Mass,pdgN, pol , pkgN!=2, nuB);  
  return err;
}


/*  ===================  Spectra ========================== */

static void getSpectrum(int forSun, double M, double m1,double m2,char*n1,char*n2, int N1, int N2,int outP, double *tab)
{
  int i,k;
  char* nn[2];
  int pdg[2];
  double mm[2],E[2],p2;

  tab[0]=M/2;  
  for(i=1;i<NZ;i++) tab[i]=0; 

  if(abs(N1)==abs(N2)) switch(abs(N1))
  { case 22: case 11: case 13: return;}

  if(N1+N2==0  || (N1==21 && N2==21) || (N1==22 && N2==22) ||  (N1==23 && N2==23)|| (N1==25 && N2==25) )  
  { 
    if(outP>0) {  if(basicNuSpectra(forSun, M/2, N1,0, tab,NULL)==0)  return; }
     else      {  if(basicNuSpectra(forSun, M/2, N1,0, NULL,tab)==0)  return; }
  } 

  nn[0]=n1;  nn[1]=n2;
  mm[0]=m1;  mm[1]=m2;  
  pdg[0]=N1; pdg[1]=N2;          
  
  
  if(M>m1+m2) p2=sqrt((M*M-(m1+m2)*(m1+m2))*(M*M-(m1-m2)*(m1-m2)))/(2*M);
  else 
  { p2=0; 
    if(abs(N1)==abs(N2)) { mm[0]=M/2; mm[1]=M/2;} else
    {
         if(N1==23 || abs(N1)==24)   mm[0]=M-m2; else mm[1]=M-m1;
    }   
  }
    
  E[0]=sqrt(mm[0]*mm[0]+p2*p2);
  E[1]=sqrt(mm[1]*mm[1]+p2*p2);

  for(k=0;k<2;k++)
  {   double dY;
      double tabAux[NZ];
      int err;
      if(abs(pdg[k])==11 || abs(pdg[k])==13 || abs(pdg[k])==22) continue;
      if(outP>0) err=basicNuSpectra(forSun,E[k],pdg[k],0,tabAux,NULL);
       else      err=basicNuSpectra(forSun,E[k],pdg[k],0,NULL,tabAux);
      if(err==0)
      { double kf; 
        dY=log(M/E[k]/2);
        switch(abs(pdg[k]))
        { case 12: case 14: case 16: 
          if(pdg[k]*outP>1) kf=1; else kf=0; break;
          default: kf=0.5;
        }  
        for(i=1;i<NZ;i++)tabAux[i]*=kf; 
        addSpectrum(tab,tabAux);
      }
      else
      { 
        double w=0,Qstat;
        numout * d2Proc;
        int l; 
        char* n[4];
        REAL m[4];
        double tab_p[NZ];
        char process[40],plib[40];
        int ntot;

        if(mm[k]==0) { fprintf(stderr,"Can not hadronize BSM zero mass %s\n",nn[k]); continue;}
                                         
        strcpy(plib,"2width_");
        sprintf(process,"%s->2*x",nn[k]);
        pname2lib(nn[k],plib+7);
        tabAux[0]=E[k];
        for(i=1;i<NZ;i++) tabAux[i]=0;
                     
        d2Proc=getMEcode(0,ForceUG,process,NULL,NULL,plib);
        if(!d2Proc) { fprintf(stderr,"Can not find decay modes for  mass %s\n",nn[k]); continue; }
        procInfo1(d2Proc,&ntot,NULL,NULL);
        
        if(Qaddress) { Qstat=*Qaddress; setQforParticle(Qaddress,nn[k]); }  
        for(l=1;l<=ntot ;l++)
        {    
            double wP=pWidth2(d2Proc,l);
            
            if(wP>0)
            { int N2=d2Proc->interface->pinfAux(l,2,NULL,NULL,NULL,NULL);
              int N3=d2Proc->interface->pinfAux(l,3,NULL,NULL,NULL,NULL); 
              procInfo2(d2Proc,l,n,m); 
              getSpectrum(forSun,E[k],m[1],m[2],n[1],n[2],N2,N3,outP, tab_p);
              for(i=1;i<NZ;i++) tabAux[i]+=wP*tab_p[i];
              w+=wP;
            }
        }
        if(Qaddress){ *Qaddress=Qstat; calcMainFunc();}
          
        if(w==0) { if(abs(pdg[k])!= abs(pNum(CDM)))   fprintf(stderr,"Can't find decays for  %s\n",nn[k]);
                   continue;
                 }
                 
       for(i=1;i<NZ;i++)tabAux[i]/=w;          
       addSpectrum(tab,tabAux);
     }           
  } 
}

static double calcSpectrum0(char *name1, char*name2, int forSun, double *Spectranu, double *SpectraNu, double *alpha)
{
  int i,k;
  double vcsSum=0,vcsSum1=0; 
  int ntot,err;
  double * v_cs;
  double M1=pMass(name1),M2=pMass(name2);
  
  char name1L[10],name2L[10], lib[20];
  char *process=malloc(maxPlistLen+20);
  numout * libPtr;
  
  Spectranu[0]=SpectraNu[0]=0.5*(M1+M2);  
  for(i=1;i<NZ;i++) Spectranu[i]=SpectraNu[i]=0;  

  pname2lib(name1,name1L);
  pname2lib(name2,name2L);
  sprintf(lib,"omg_%s%s",name1L,name2L);
  sprintf(process,"%s,%s->AllEven,1*x{%s",name1,name2,EvenParticles());
// Warning!!   in should be done in the same manner as annihilation libraries for Omega
  libPtr=getMEcode(0,ForceUG,process,NULL,NULL,lib);
  free(process);
  
  if(!libPtr) return 0;
  if(Qaddress && *Qaddress!=M1+M2) 
  { *Qaddress=M1+M2;
     calcMainFunc();
  }   
  passParameters(libPtr);
  procInfo1(libPtr,&ntot,NULL,NULL); 
  
  v_cs=malloc(sizeof(double)*ntot);
  (*libPtr->interface->twidth)=0;
  
  for(k=0;k<ntot;k++)
  { REAL m[4];
    char *N[4];
    int pdg[4];
    int l,l_;
    double br,wV;

    for(i=0;i<4;i++) N[i]=libPtr->interface->pinf(k+1,i+1,m+i,pdg+i);
    cc23=NULL;
    v_cs[k]=0;
    if(VZdecay||VWdecay)
    {  int nVV;
       int vd[4]={0,0,0,0};
       for(l=2;l<4;l++) if((pdg[l]==23&&VZdecay) || (abs(pdg[l])==24&&VWdecay)) vd[l]=1;
            
       for(l=2;l<4;l++) if(vd[l]) break;
       if(l<3)
       {  l_=5-l; 
          if(vd[l_])
          { nVV=2;
            if(m[l_]>m[l]) { l=l_; l_=5-l;}
          } else nVV=1;
          
          if(m[0]+m[1] >  m[l_] +20  && m[0] + m[1] <  m[2]+m[3] + 4*nVV)
           cc23=xVtoxll(2,2,N,pdg,l,&wV,&br);                
       }
    }

    if(cc23)
    { int i3W;  
      double  r,m1,v0=0.001;
      for(i3W=2;i3W<5;i3W++) if(strcmp(cc23->interface->pinf(1,i3W+1,NULL,NULL),N[l_])==0) break;
      r=v0*cs23(cc23,1,v0*M1*M2/(M1+M2),i3W,NULL)/br;
       
      if(pdg[l_]==23 || abs(pdg[l_])==24)
      { double wV2;
        
        wV2=pWidth(N[l_],NULL);
        r*=decayPcmW(M1+M2,m[l],m[l_],wV,wV2,0)/decayPcmW(M1+M2,m[l],m[l_],wV,0,0);
        if(pdg[l]==pdg[l_]) r/=2;
      }
      v_cs[k]=r;
      vcsSum+=r;                     
    }
    else if((m[2]+m[3])<m[0]+m[1])
    { 
#ifdef V0    
      v_cs[k]=V0*cs22(libPtr,k+1,V0*m[0]/2,-1.,1.,&err);
#else 
      v_cs[k]= vcs22(libPtr,k+1,&err);
#endif 
      if(v_cs[k]<0) v_cs[k]=0; 
      vcsSum+=v_cs[k];
    } else v_cs[k]=-1;
  }
  
  for(k=0;k<ntot ;k++) if(v_cs[k]>0)
  { char * N[4];
    REAL m[4];
    int l,pdg[2];
    int PlusAok=0;
    double tab2[NZ];

    procInfo2(libPtr,k+1,N,m);
    for(l=0;l<2;l++)  pdg[l]=qNumbers(N[2+l],NULL,NULL,NULL);
    if(N[2][0]=='~' || N[3][0]=='~') vcsSum1+=v_cs[k]; 
#ifdef PRINT
       { char txt[100];
         sprintf(txt,"%s,%s -> %s %s", N[0],N[1],N[2],N[3]);
         printf("  %-20.20s  %.2E\n",txt,v_cs[k]*2.9979E-26);
       }
#endif   
    getSpectrum(forSun, pMass(name1)+pMass(name2), m[2], m[3],N[2],N[3],pdg[0],pdg[1], 1,tab2); 
    for(i=1;i<NZ;i++) Spectranu[i]+=tab2[i]*v_cs[k]/vcsSum;
    getSpectrum(forSun, pMass(name1)+pMass(name2), m[2], m[3],N[2],N[3],pdg[0],pdg[1],-1,tab2);
    for(i=1;i<NZ;i++) SpectraNu[i]+=tab2[i]*v_cs[k]/vcsSum;  
  } 
  free(v_cs);
  if(alpha) { if(vcsSum>0) *alpha=vcsSum1/vcsSum; else *alpha=0;}
  return  vcsSum*2.9979E-26;
}

/* ========  Sun[Earth] neutrino fluxes ============= */

#define Gconst (0.7426E-30) /* m/g */
#define KelvinEv (8.61734E-05) /* ev */
#define Etime  (1.5E17)     /* time of existence of Sun and Earth in seconds */ 

static double C[2],An[2],Ev[2],Alpha;

static void deriv1(double t,double*n,double*dn) { dn[0]=C[0]-An[0]*n[0]*n[0]-Ev[0]*n[0];

//printf(" C=%E An=%E Ev=%E\n", C[0], An[0]*n[0]*n[0], Ev[0]*n[0]);
} 
  
static void deriv2(double t,double*n,double*dn) 
{ dn[0]=C[0]+ 0.5*Alpha*An[1]*n[1]*n[1] - An[1]*n[0]*n[0] - An[0]*n[0]*n[1] -Ev[0]*n[0]; 
  dn[1]=C[1]+ 0.5*Alpha*An[1]*n[0]*n[0] - An[1]*n[1]*n[1] - An[0]*n[0]*n[1] -Ev[1]*n[1];
}

int neutrinoFlux(double (* fvf)(double), int forSun, double* nu, double * Nu)
{
  int i,n,err;
  double vcs0,vcs1;

  double nu_[NZ],Nu_[NZ];
  char *name, *aname;
  double pA0[2],pA5[2],nA0[2],nA5[2];
  double R,Prop,Cr0,Cr1,Dv; 
  double Veff;
  double rho,T;

  for(i=1;i<NZ;i++) nu[i]=Nu[i]=0;
  
  if(CDM1&&CDM2) { printf(" The 'neutrinoFlux' code is still not upgrated for 2DM case\n");  return 1;}
  if(CDM1) CDM=CDM1; else CDM=CDM2;

  double M_cdm=pMass(CDM);
     
  if(nu)nu[0]=M_cdm;
  if(Nu)Nu[0]=M_cdm;

  err=basicNuSpectra(forSun,M_cdm,5,0,NULL,NULL);
  if(err) return err;
  
  err=nucleonAmplitudes(CDM,pA0,pA5,nA0,nA5);

  if(err) return err;
  
  Cr0=forSun? captureSun(fvf,M_cdm,pA0[0],nA0[0],pA5[0],nA5[0]):captureEarth(fvf,M_cdm, pA0[0],nA0[0],pA5[0],nA5[0]); 

  if(pA0[0]==pA0[1] && nA0[0]==nA0[1] &&pA5[0]==pA5[1]&&nA5[0]==nA5[1]) Cr1=Cr0; else
  Cr1=forSun? captureSun(fvf,M_cdm,pA0[1],nA0[1],pA5[1],nA5[1]):captureEarth(fvf,M_cdm,pA0[1],nA0[1],pA5[1],nA5[1]);  
                         
  name=CDM;
  aname=pdg2name(-pNum(CDM));
  if(!aname) aname=name;
  if(forSun) R=150E6;  else R=6378.1; /* Distance to Sun/Earth in [km] */  
  Prop=31556925.2/(4*M_PI*R*R);       /* for Year*km^2 */
    
  vcs0= calcSpectrum0(name,aname,forSun, nu,Nu,NULL);
  { 
     double r_,v_,r095,S,ph,Veff1;
     for(i=0,r_=0;i<10;i++) 
     { T=polint1(r_/Rcm,nTab,rTab,tTab)*KelvinEv*1E-9;
       rho= polint1(r_/Rcm,nTab,rTab,rhoTab);
       r_= sqrt(6*T/(Gconst*100*rho*M_cdm))/M_PI;
       if(r_>Rcm) r_=Rcm;
     }

     rho=polint1(r_/Rcm,nTab,rTab,rhoTab);
     T=polint1(r_/Rcm,nTab,rTab,tTab);
     r095=polint1(T*0.95,nTab,tTab,rTab);
     if(r095>1) r095=1;
     T*=KelvinEv*1E-9;
     r095*=Rcm;
     Veff1=pow(r_*M_PI,3)/8;
     if(forSun) S=sigmaSun(M_cdm,pA0[0],nA0[0],pA5[0],nA5[0],r095*r095*r095);
     else       S=sigmaEarth(M_cdm,pA0[0],nA0[0],pA5[0],nA5[0],r095*r095*r095);
     v_=sqrt(8*T/(M_PI*M_cdm));
     ph=phiTab[0]*M_cdm/T;
     ph=ph*exp(-ph);
     
     Ev[0]=v_*S/Veff1*ph*Vlight*100;
     if(strcmp(name,aname))
     { if(forSun)S=sigmaSun(M_cdm,pA0[1],nA0[1],pA5[1],nA5[1],r095);
       else S=sigmaEarth(M_cdm,pA0[1],nA0[1],pA5[1],nA5[1],r095);
       Ev[1]=v_*S/Veff1*ph*Vlight*100; 
     } 
     Veff=pow(Gconst*100*rho*M_cdm/T/3,-1.5);
  } 

  if(strcmp(name,aname))
  { double G01,G00,G11;
    double N[2]={0,0};
    int err;
    vcs1=calcSpectrum0(name,name,forSun, nu_,Nu_,&Alpha);
    
    C[0]=Cr0*(1+dmAsymm)/2;                             
    C[1]=Cr1*(1-dmAsymm)/2;       
       
    An[0]=vcs0/Veff;
    An[1]=vcs1/Veff;
    err=odeint(N,2, 0 ,Etime , 1.E-3,  Etime/10, deriv2);
//printf("err=%d N={%E,%E}\n",err, N[0],N[1]);
    G00=0.5*An[1]*N[0]*N[0];
    G11=0.5*An[1]*N[1]*N[1];
    G01=    An[0]*N[0]*N[1];    
    { /* symbolic solution for large Etime */
      double alpha,beta,x,G01_,G00_,G11_;        
      alpha=vcs1/vcs0;
      beta= Cr0/Cr1;
      x=(beta-1 + sqrt((beta-1)*(beta-1) + 4*beta*alpha*alpha))/(2*alpha);
      /* x= rho_particle/rho_antiparticle  */ 
      G01_=Cr0/(1+alpha*x);
      G00_=0.5*G01_*alpha*x;
      G11_=0.5*G01_*alpha/x;
//      printf("x=%E G00 = %E/%E  G01 = %E/%E G11 = %E/%E\n",x, G00,G00_,G01,G01_,G11,G11_);
    }
    for(i=1;i<NZ;i++) 
    { nu[i]=Prop*(G01*nu[i]+G00*nu_[i]+G11*Nu_[i]); 
      Nu[i]=Prop*(G01*Nu[i]+G00*Nu_[i]+G11*nu_[i]);
    }  
    vcs1=calcSpectrum0(aname,aname,forSun, nu_,Nu_,NULL);     
    for(i=1;i<NZ;i++)
    { nu[i]+=Prop*(G11*nu_[i]);
      Nu[i]+=Prop*(G11*Nu_[i]);
    }              

  } else   
  {   
    int err;  
    double N=0;
    double rf;
    C[0]=Cr0;
    An[0]=vcs0/Veff;
    err=odeint(&N,1, 0 ,Etime , 1.E-3,  Etime/10, deriv1);
    rf=0.5*An[0]*N*N*Prop;

//printf("Rate Factor = %E(num.sol), =%E(formula)\n",rf,   0.5*Cr0*pow(tanh(Etime*sqrt(Cr0*vcs0/Veff)),2)*Prop);
    for(i=1;i<NZ;i++) { nu[i]*=rf;  Nu[i]*=rf;}
  }
  return 0;
}

static double cosFi_stat=0.1;//   0.34;

// ========================== Earth Attenuation =====================


static double inte_grand(double x){ return polint1(sqrt(1+x*x-2*x*cosFi_stat) ,EARTHPOINTS,rEarth,rhoEarth); }

void makeTabl(void)
{ int i;

  readEarthData();
  for(i=0;i<=20;i++)
  { cosFi_stat=i/20.;
    printf("cos=%E rho=%E\n", cosFi_stat,2*RE*simpson(inte_grand,0,cosFi_stat,1.E-3,NULL));  
  }
}


double nuAttenuation(int nu, double cs,double E)  // 1/pb
{ double NA=6.002141E23;
  double GF=1.16637E-5;  /* GeV^(−2) */
  double GeVcm=0.50677E14;
  double C=2*GF*GF*mp_gev/M_PI/GeVcm/GeVcm;    

  double cos[21]={0.000000E+00,5.000000E-02,1.000000E-01,1.500000E-01,2.000000E-01,2.500000E-01,3.000000E-01,3.500000E-01,4.000000E-01,4.500000E-01,5.000000E-01,5.500000E-01,6.000000E-01,6.500000E-01,7.000000E-01,7.500000E-01,8.000000E-01,8.500000E-01,9.000000E-01,9.500000E-01,1.000000E+00};
  double rho[21]={0.000000E+00,1.775227E+08,3.755063E+08,6.110139E+08,8.372967E+08,1.056885E+09,1.293592E+09,1.538506E+09,1.841714E+09,2.167974E+09,2.553811E+09,2.913635E+09,3.278042E+09,3.656525E+09,4.047799E+09,4.461401E+09,4.899167E+09,6.159285E+09,7.876339E+09,9.341818E+09,1.096080E+10};
   
  double  tag=NA*polint3(fabs(cs),  21,cos,rho);  // particles/cm^2
  double sigma;
    
  if(nu>0)  sigma=C*0.5*E*((0.15+0.25) +(0.04+0.06)/3);  // sigma[cm^2]
  else      sigma=C*0.5*E*((0.04+0.06) +(0.15+0.25)/3);
  return exp(-tag*sigma);
}

#ifdef QQ
double sigmaNu(int nu,double E)  //  pb
{ 
   double GF=1.16637E-5;  /* GeV^(−2) */
   double GeVcm=0.50677E14;
   double C=2*GF*GF*mp_gev/M_PI/GeVcm/GeVcm;    
       
    if(nu>0) return C*0.5*E*((0.15+0.25) +(0.04+0.06)/3);  
    else     return C*0.5*E*((0.04+0.06) +(0.15+0.25)/3);     
}
#endif


/* ==================  Muon flux =========  0906.4364  =============  */

static double *SpN_stat=NULL;
static double tabNuSpectrum(double E) {return SpectdNdE(E,SpN_stat);}
static double (*nuSpectrum)(double)=tabNuSpectrum;
static double Enu_stat,Emu_stat,alpha_stat,beta_stat;


static int inPr=1; /* proton;    0 for neutron */
static int inNu=1; /* neutrino; -1 for anti-neutrino */


static double dSigmadE_nu2mu(int inNu,  double Enu,double Emu)   /*result in   1/cm^2/GeV */ 
{  double GF=1.16637E-5;  /* GeV^(-2) */
   double GeVcm=0.50677E14;
   double C=GF*GF*mp_gev/M_PI/GeVcm/GeVcm;
   double ap=0,an=0,bp=0,bn=0;
   switch (inNu)
   { case  1: ap=0.15,bp=0.04; an=0.25,bn=0.06; break; 
     case -1: ap=0.06,bp=0.25; an=0.04,bn=0.15; break; 
   }
             
   return C*((ap+an)+(bp+bn)*Emu*Emu/(Enu*Enu));       
} 

static double integrandX(double x)
{  double g=alpha_stat/beta_stat;  
   double ex=exp(x*beta_stat);
//   double q= (0.10566 /* muon mass*/)/65865.4 /* (muon_life_time)*c in cm */)/alpha_stat;
   double q= (0.10566 /* muon mass*/)/(6665865.4 /* (muon_life_time)*c in cm */)/alpha_stat;
   double Emu_prim=(Emu_stat+g)*ex-g;
   double Psurv=pow(Emu_stat*(Emu_prim+g)/(Emu_prim*(Emu_stat+g)),q);
   return  dSigmadE_nu2mu(inNu,  Enu_stat, Emu_prim )*ex*Psurv; 
} 

static double integrandEnuUpward(double lnE)
{ double E;

  E=exp(lnE);
  Enu_stat=E;
  double g=alpha_stat/beta_stat;
  double xMax=log((Enu_stat+g)/(Emu_stat+g))/beta_stat ;  
   
  return  E*nuSpectrum(E)*simpson(integrandX,0,xMax,1.E-4,NULL);
}


static double integrandEnuUpward_sigma(double e)
{ double E;
  if(e==0) return 0;
  E=1/e;
  Enu_stat=E;
  double g=alpha_stat/beta_stat;
  double xMax=log((Enu_stat+g)/(Emu_stat+g))/beta_stat ;  
  double sigma=0.0122/pow(1E-3*E,0.7);
  sigma*=sigma; 
  return  sigma*E*E*nuSpectrum(E)*simpson(integrandX,0,xMax,1.E-4,NULL);
}


void muonUpward(double*nu, double*Nu, double*mu)
{
   int i,k,l;
   double C;
   double  NA=6.002141E23 /* mol-1*/;
   double*Sp[2]={nu,Nu};
   double rho,M;
   double *sigma=NULL;
   if(nu && Nu) M=nu[0]>Nu[0]?nu[0]:Nu[0];
   else if(nu) M=nu[0]; else if(Nu) M=Nu[0]; else return;
          
   
   if(forRocks)   
   {   rho=2.6;
      alpha_stat=0.002*rho;
      beta_stat=3.0E-6*rho;
   }   
   else      
   {
      rho=1;  
      alpha_stat=0.00262*rho;
      beta_stat=3.5E-6*rho;
   }
   nuSpectrum=tabNuSpectrum;  

   mu[0]=M;
   mu[1]=0;
   if(sigma){ sigma[0]=nu[0]; sigma[1]=0; for(i=1;i<NZ;i++)sigma[i]=0;}   
   for(i=2;i<NZ;i++) 
   {  double sigma_mem;
      Emu_stat=M*exp(Zi(i));
      mu[i]=0;
      if(sigma) sigma[i]=0;
      if(Emu_stat>0.01)
      {  for(k=0;k<2;k++) if(Sp[k])
         {
           inNu=1-2*k;
           SpN_stat=Sp[k];
           mu[i]+=simpson(integrandEnuUpward,log(Emu_stat),log(M),1.E-3,NULL);
           if(sigma) sigma[i]+=simpson(integrandEnuUpward_sigma,1/M,1/Emu_stat,1.E-4,NULL); 
         }    
         if(sigma) {  sigma_mem=sigma[i]= sqrt(sigma[i]/mu[i]); sigma[i]*=Emu_stat;}
         mu[i]*=rho*NA*Emu_stat;        
      } else { if(sigma) sigma[i]=sigma_mem*Emu_stat;}  
   }
}

//#define QQ
#ifdef QQ

static double integrandX_I(double Emu_prim)
{  double g=alpha_stat/beta_stat;  
   double q= (0.10566 /* muon mass*/)/(6665865.4 /* (muon_life_time)*c in cm */)/alpha_stat;
   double Psurv=pow(Emu_stat*(Emu_prim+g)/(Emu_prim*(Emu_stat+g)),q);

   return  dSigmadE_nu2mu(inNu,  Enu_stat, Emu_prim )*Psurv; 
} 


static double integrandEnuUpward_I(double E)
{ 

  Enu_stat=E;   
  
  return  nuSpectrum(E)*simpson(integrandX_I,Emu_stat,E,1.E-4,NULL);
}



void muonUpward_I(double*nu, double*Nu, double*mu)
{
   int i,k,l;
   double C;
   double  NA=6.002141E23 /* mol-1*/;
   double*Sp[2]={nu,Nu};
   double rho,M;
   double *sigma=NULL;
   if(nu && Nu) M=nu[0]>Nu[0]?nu[0]:Nu[0];
   else if(nu) M=nu[0]; else if(Nu) M=Nu[0]; else return;
          
   
   if(forRocks)   
   {   rho=2.6;
      alpha_stat=0.002*rho;
      beta_stat=3.0E-6*rho;
   }   
   else      
   {
      rho=1;  
      alpha_stat=0.00262*rho;
      beta_stat=3.5E-6*rho;
   }
   nuSpectrum=tabNuSpectrum;  

   mu[0]=M;
   mu[1]=0;
   for(i=2;i<NZ;i++) 
   {  double sigma_mem;
      Emu_stat=M*exp(Zi(i));
      mu[i]=0;
      if(Emu_stat>0.01)
      {  for(k=0;k<2;k++) if(Sp[k])
         {
           inNu=1-2*k;
           SpN_stat=Sp[k];
           mu[i]+=simpson(integrandEnuUpward_I,Emu_stat,M,1.E-4,NULL)/(alpha_stat+beta_stat*Emu_stat); 
         }    
         mu[i]*=rho*NA*Emu_stat;        
      } 
   }
}


#endif 



static double integrandEnuContained(double e) 
{  double E; 
   if(e==0) return 0; 
   E=1/e; return  E*E*nuSpectrum(E)*dSigmadE_nu2mu(inNu,E,Emu_stat) ; 
}

void muonContained(double*nu,double*Nu,double rho, double*mu)
{
  int i,k,l;
  double C;
  double  NA=6.002141E23 /* mol-1*/;
  double*Sp[2]={nu,Nu};
  double M_cdm;
  
  if(nu && Nu) M_cdm=nu[0]>Nu[0]?nu[0]:Nu[0];
  else if(nu)  M_cdm=nu[0]; else if(Nu) M_cdm=Nu[0]; else return;
       
  
  nuSpectrum=tabNuSpectrum; 
  mu[0]=nu[0]; 
  mu[1]=0;
  for(i=2;i<NZ;i++)
  {  Emu_stat=M_cdm*exp(Zi(i));
     mu[i]=0;
     if(Emu_stat>0.01) for(k=0;k<2;k++)
     { inNu=1-2*k;
       if(Sp[k])
       {
         SpN_stat=Sp[k];
         mu[i]+=simpson(integrandEnuContained,1/M_cdm, 1/Emu_stat,1.E-4,NULL);
       }  
     }  
     mu[i]*=rho*NA*Emu_stat*1E5;
  }                  
}

// Background 


double  ATMdNudE(double E) // for NuBar *1.35/1.95 
{
  int i;
  const double gamma=1.74, a=0.018, b=0.024,c=0.0069, e=0.00139;
  return 1.95E17*pow(E,-gamma-1)*(a/(1+b*E*cosFi_stat)+c/(1+e*E*cosFi_stat));
} 

double  atmNuFlux_(int nu,double cs, double E)
{ double tmp=cosFi_stat;
  double r;
  cosFi_stat=cs;
  r=ATMdNudE(E);
  cosFi_stat=tmp;
  if(nu<0) r*=1.35/1.95;
  return r;
}    


double  ATMmuonUpward(double cosFi, double E)
{ 
   int k,l;
   double  NA=6.002141E23 /* mol-1*/;
   double Nrate[2]={0.5,0.5}; /* proton , neutron */
   double rho;
   double mu=0;
   if(forRocks)    
   {  rho=2.6;
      alpha_stat=0.002*rho;
      beta_stat=3.0E-6*rho;
   } else 
   {    
      rho=1;
      alpha_stat=0.00262*rho;
      beta_stat=3.5E-6*rho;
   }
   Emu_stat=E;
   cosFi_stat=cosFi;
   nuSpectrum=ATMdNudE;
   
   for(k=0;k<2;k++) for(l=0;l<2;l++)
   {
      inNu= 1-2*k;
      inPr=-1+2*l;
      mu+= (k==0?1: 1.35/1.95)*simpson(integrandEnuUpward,log(E),log(1E6), 1.E-3,NULL)*Nrate[l];
   }   
   return  mu*rho*NA;
}

double  ATMmuonContained(double cosFi, double E,double rho)
{ 
   int k,l;
   double  NA=6.002141E23 /* mol-1*/;
   double Nrate[2]={0.5,0.5}; /* proton , neutron */
   double mu=0;
   
   alpha_stat=0.002*rho;
   beta_stat=3.0E-6*rho;
   
   
   Emu_stat=E;
   cosFi_stat=cosFi;
   nuSpectrum=ATMdNudE;
   
   for(k=0;k<2;k++) for(l=0;l<2;l++)
   {
      inNu= 1-2*k;
      inPr=-1+2*l;
      mu+= (k==0?1: 1.35/1.95)*simpson(integrandEnuContained,0,1/E,1.E-4,NULL)*Nrate[l];
   }   
   return  mu*rho*NA*1.E5;      
}

typedef  struct 
{ double E[31];
  double data[31][10];
}  atmDataStr;

static atmDataStr *atmNu=NULL, *atmNuBar=NULL;

static int readATMnuTab(void)
{ FILE*f;
  char fname[300];
  int i;

  if(atmNu && atmNuBar) return 0;  
    
  sprintf(fname,"%s/sources/data_nu/atm_nu_tab.dat",micrO);
  f=fopen(fname,"r");
  if(!f) return 1;
  atmNu=malloc(sizeof(atmDataStr)); 
  for(i=0;i<5;i++) fscanf(f," %*[^\n]");
  for(i=0;i<31;i++)
  {  int j;
     double norm;
     fscanf(f," %lf", atmNu->E+i);
     for(j=0;j<10;j++) fscanf(f," %lf",  atmNu->data[i]+j);
     fscanf(f," %lf",&norm);
     for(j=0;j<10;j++)   atmNu->data[i][j]*=norm;
  }
  fscanf(f," %*[^\n]");
  atmNuBar=malloc(sizeof(atmDataStr));
  for(i=0;i<31;i++)
  {  int j;
     double norm;
     fscanf(f," %lf", atmNuBar->E+i);
     for(j=0;j<10;j++) fscanf(f," %lf",  atmNuBar->data[i]+j);
     fscanf(f," %lf",&norm);    
     for(j=0;j<10;j++)   atmNuBar->data[i][j]*=norm;
  }
  fclose(f);
  return 0;
} 

double atmNuFlux(int nu,double cs, double E)
{
  double alpha;
  double cstab[10]={0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95};
  atmDataStr *data;
  int i;

  if(E<10 || E>1E4)  printf("ATM neutrino flux: Energy %E out of table range\n",E);

  if((!atmNu || !atmNuBar) && readATMnuTab() ) { printf(" Can not read data file for atmospheric neutrino\n"); return -1;}  
  if(nu>=0) data=atmNu; else data=atmNuBar; 
  
  i=10*log10(E/10);

  if(i<0) i=0;else if(i>29) i=29; 
  
  alpha= (log(data->E[i+1])-log(E) )/(log(data->E[i+1])-log(data->E[i]));

  return   1.E6*365*24*60*60*pow(polint1(cs, 10,cstab,data->data[i]),alpha)*pow(polint1(cs, 10,cstab,data->data[i+1]),1-alpha);    
}  

 