#include<stdio.h>
#include<stdlib.h>
#include"micromegas.h"
#include"micromegas_aux.h"
#include"ic22.h"

//#define NOCHAN

//#define DS

double IC22nuAr_(double E)
{ if(E<=36) return 0;
  else if(E<72.7108)  return  4.403e-06 *1E-6;
  else if(E<105.737)  return  9.525e-05 *1E-6;
  else if(E<153.765)  return  0.0007292 *1E-6;
  else if(E<223.607)  return  0.002949  *1E-6;
  else if(E<325.172 ) return  0.009057  *1E-6;
  else if(E<472.871 ) return  0.02201   *1E-6;
  else if(E<687.656 ) return  0.04820   *1E-6;
  else /*if(E<1000 ) */return 0.1002    *1E-6;  
}

double IC22nuAr(double E)
{ 
  double lnE_[8]={  4.099, 4.474, 4.848, 5.223, 5.597, 5.972, 6.346, 6.720};
  double lnAr[8]={-12.333,-9.259,-7.224,-5.826,-4.704,-3.816,-3.032,-2.300};

  if(E<36) return 0;
  if(E>1E4)E=1.E4;
  return exp(polint3(log(E),8,lnE_,lnAr))*1E-6; 
}



double IC22sigma(double E)
{  double  lnE_[8]= {  4.099, 4.474, 4.848, 5.223, 5.597, 5.972, 6.346, 6.720};
   double  fiTab[8]={  3.69,  3.20,  2.96 ,  2.67, 2.47,  2.32,  2.24,  2.05}; 
//    double  fiTab[8]={  4.5,  3.9,  3.5 ,  3.2, 3.0, 2.8,  2.7,  2.5};
  return  (M_PI/180)*polint3(log(E),8,lnE_,fiTab);
}
  
double IC22nuBarAr_(double E)
{ if(E<50) return 0;
  else if(E<72.7108)  return  3.154E-06 *1E-6;
  else if(E<105.737)  return  8.150E-05 *1E-6;
  else if(E<153.765)  return  0.0005786 *1E-6;
  else if(E<223.607)  return  0.002301 *1E-6;
  else if(E<325.172 ) return  0.006607 *1E-6;
  else if(E<472.871 ) return  0.01649 *1E-6;
  else if(E<687.656 ) return  0.03537 *1E-6;
  else /*if(E<1000 ) */return 0.06991*1E-6; 
}




double IC22nuBarAr(double E)
{ 
   double lnE_[8]={  4.099, 4.474, 4.848, 5.223, 5.597, 5.972, 6.346, 6.720};
   double lnAr[8]={-12.667,-9.415,-7.455,-6.074,-5.020,-4.105,-3.342,-2.660};
   if(E<50) return 0;
   if(E>1E4)E=1.E4;
   return exp(polint3(log(E),8,lnE_,lnAr))*1E-6;      
}



double IC22BGdCos(double cs)
{
  double fi_bg[25]={2.034E+00,6.017E+00,9.895E+00,1.370E+01,1.725E+01,2.086E+01,2.431E+01,2.698E+01,2.937E+01,3.160E+01,
                    3.386E+01,3.595E+01,3.821E+01,4.022E+01,4.188E+01,4.368E+01,4.547E+01,4.727E+01,4.902E+01,5.080E+01,
                    5.227E+01,5.371E+01,5.399E+01,5.184E+01,4.927E+01};


// fi_bg angular distribution of backgraund simulation. From file  DarkSUSY/IC_data/BG_distributions_IC22.dat
                     
  int i;
  for(i=0;i<25;i++) fi_bg[i]/=cos(i*M_PI/180)-cos((i+1)*M_PI/180.);
  double cs_arr[25];
  for(i=0;i<25;i++) cs_arr[i]=cos((i+0.5)*M_PI/180.);
  return  polint1(cs,25,cs_arr,fi_bg);
}

static double *nuStat,*nuBarStat;
static double cs_stat;

static double E_integ(double E)
{  double res=0,s;
   if(nuStat)   res+=SpectdNdE(E,nuStat)*IC22nuAr(E);
   if(nuBarStat)res+=SpectdNdE(E,nuBarStat)*IC22nuBarAr(E);
   s=IC22sigma(E);
      
   s*=s;   
   return res*exp((cs_stat-1)/s)/s;
}

static double cs_tab[40];
static double dNsigTab[40];

void  getdNSdCos( double *nu,double *nu_)
{
   double M=0;
   int i;
   if(nu) M=nu[0];
   if(nu_ && nu_[0]>M) M=nu_[0];

   nuStat=nu;
   nuBarStat=nu_;


   for(i=0;i<40;i++) 
   { cs_stat=1-i*(1-cos(10./180.*M_PI))/40.;
     cs_tab[i]=cs_stat;
     dNsigTab[i]=simpson(E_integ,1,M,1.E-3,NULL); 
   }  
}

double dNSdCos(double cs) {return 104./365.*polint3(cs,40,cs_tab,dNsigTab);} 


IC22chanStr*IC22chan=NULL;

int IC22histRead(void)
{  FILE*F;
   char fname[300];
   int i;
   double E1,E2;

   if(IC22chan) return 0;
   
   sprintf(fname,"%s/sources/data_nu/ic22hist.dat",micrO);
        
   F=fopen(fname,"r");
   if(!F) return 1;
   

   IC22chan=malloc(21*sizeof(IC22chanStr));
         
   for(i=0;2==fscanf(F," %*s %lf %lf",&E1,&E2); i++)
   {  int j; 
      double s;
      IC22chan[i].E1=pow(10,E1);  IC22chan[i].E2=pow(10,E2);
      for(j=0,s=0;j<17;j++) {fscanf(F," %lf", IC22chan[i].prob+j); s+=IC22chan[i].prob[j];}
//      printf("sum(prob)=%E\n",s);      
   }   
   fclose(F);                                                    
}

static   double pp(double cs,  int n, double a)
#ifdef NOCHAN
   { return IC22BGdCos(cs) +a*dNSdCos(cs)/1.2;}

#else
   { double s= IC22BGdCos(cs)*IC22chan[0].prob[n];
     int i;
     for(i=1;i<=20;i++) if(IC22chan[i].n>0)
     { double s2=IC22chan[i].s2;
       s+= a*IC22chan[i].n*exp((cs-1)/s2)/s2*IC22chan[i].prob[n];
     } 
     return s;  
   }
#endif                 
static   double L[21];

static    double P(double x)
   { double xx[21]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};
     return exp(polint3(x,21,xx,L));
   }
                  
//double fiMax=8;   

   double exLevIC22(double * nu, double*NU, double *B)
   {  int i;
      double cs,cs_;
      typedef struct{double cs; int n;} eventStr;
      static  eventStr*  events=NULL;
      static int Nev;
      double fiMax=8; 

#ifdef NOCHAN 
      getdNSdCos(nu,NU);
#endif   

#ifdef DS
      getdNSdCos(nu,NU);
#endif
      

      if(fiMax>10) { printf("Maximum angel reset 10 degrees\n"); fiMax=10;}
      if(!events)
      {  char fname[300];
         double dfi;
         int nch;
         events=malloc(180*sizeof(eventStr));
         sprintf(fname, "%s/sources/data_nu/IC22_events_25.dat",micrO); 
         FILE*F=fopen(fname,"r");
         cs_=cos(10./180.*M_PI);    
         for(Nev=0; fscanf(F," %lf %lf %d",&cs,&dfi, &nch)==3; ) if(cs>=cs_)
         { if(Nev==180) break;
           events[Nev].cs=cs;
           events[Nev].n=nch;
           Nev++; 
         }
         fclose(F);
         IC22histRead();  
      }

      cs_=cos(fiMax/180*M_PI);

      double  nu_[NZ], NU_[NZ];
      for(i=0;i<NZ;i++) {nu_[i]=nu[i];NU_[i]=NU[i];} 
      spectrMult(nu_,IC22nuAr);
      spectrMult(NU_,IC22nuBarAr); 
      addSpectrum(nu_,NU_);
      for(i=0;i<NZ;i++) NU_[i]=nu_[i];
      spectrMult(NU_,IC22sigma); spectrMult(NU_,IC22sigma);
      
      double M=nu_[0];
      for(i=1; i<=20;i++)
      {  
         double E1=IC22chan[i].E1,E2=IC22chan[i].E2;
         if(E1<M)
         {
           if(E2>M) E2=M;
           IC22chan[i].n=spectrInt(E1,E2,nu_);
           if(IC22chan[i].n<=0) IC22chan[i].s2=0; else
           IC22chan[i].s2=spectrInt(E1,E2,NU_)/IC22chan[i].n;
         } else {IC22chan[i].n=0; IC22chan[i].s2=0;} 
         IC22chan[i].n*=104/365./1.2; 
      }
      double Nbg=simpson(IC22BGdCos,cs_,1,1E-3,NULL);
      double Ns=0;
#ifdef NOCHAN
      Ns=simpson(dNSdCos,cs_,1,1.E-3,NULL)/1.2;
//printf("Ns(noChan=%E ",Ns);      
#else   
      Ns=0;           
      for(i=1;i<20;i++) if(IC22chan[i].n>0) Ns+=IC22chan[i].n*(1-exp((cs_-1)/IC22chan[i].s2));
//      printf("Ns(with Ch)=%E\n",Ns);
#endif

#ifdef DS 
   {  double rate=0,Ns,Nbg,cs_=-1,c,s;
      int nData=0;      
      for(cs=1;cs>0.99;cs-=0.0001)       
      { 
         Ns= simpson(dNSdCos,cs,1,1.E-3,NULL)/1.2;
         Nbg=simpson(IC22BGdCos,cs,1,1E-3,NULL);
        if(rate<Ns/sqrt(Nbg)) {rate=Ns/sqrt(Nbg); cs_=cs;}   
      } 
      Ns= simpson(dNSdCos,cs_,1,1.E-3,NULL)/1.2;
      Nbg=simpson(IC22BGdCos,cs_,1,1E-3,NULL);   
      for(i=0;i<Nev;i++) if(events[i].cs>cs_) nData++;

printf("nData=%d Ns=%E Nbg=%E fi=%E\n",nData,Ns,Nbg,180/M_PI*acos(cs_));      
      
      return 1-exp(-Ns)*pow(1+Ns/Nbg,nData);
      c=exp(-Ns-Nbg);
      for(i=0,s=0;i<nData;i++)
      { 
        s+=c;
        c*=(Ns+Nbg)/(i+1);
      }
printf("s=%E\n",s);   
//exit(0);   
      return 1-s;
   }      
#endif      
      double nData=0;
      for(i=0;i<Nev;i++) if(events[i].cs>cs_) nData++;
      int k;
      for(k=0;k<21;k++)
      { double a= 0.1*k;
        L[k]=nData*log(Nbg+a*Ns)-Nbg-a*Ns;
        for(i=0;i<Nev;i++)if(events[i].cs>cs_) L[k]+=log(pp(events[i].cs,events[i].n-10,a)/(Nbg+a*Ns)); 
      }
      
      for(k=1;k<21;k++) L[k]-=L[0];
      L[0]=0;
      if(L[20]>L[15]) return 0;     
      double dI=exp(L[20])/(2*(L[15]-L[20]));   
      double int1=simpson(P,0.,1,1E-3,NULL);
      double int2=simpson(P,1.,2,1E-2,NULL);
      if(B) *B=P(1)/P(0);
      return (int1)/(int1+int2+dI);
   }
static int prn=0;   
      typedef struct{double cs; int n;} eventStr;

static  eventStr*  events=NULL;
static int Nev;

//#define TEST_ANOMALY   
#ifdef TEST_ANOMALY   
   static double exLevIC22_random(double * nu, double*NU,int rand)
   {  int i;
      double cs,cs_;
      double Nev_;
      cs_=cos(8./180.*M_PI);    
         

      if(!events)
      {  char fname[300];
         double dfi;
         int nch;
          
         events=malloc(1000*sizeof(eventStr));
         sprintf(fname, "%s/sources/data_nu/IC22_events_25.dat",micrO); 
         FILE*F=fopen(fname,"r");
         for(Nev=0; fscanf(F," %lf %lf %d",&cs,&dfi, &nch)==3; ) //    if(cs>=cs_) // to reproduce  experimental result
         { if(Nev==1000) break;
           events[Nev].cs=cs;
           events[Nev].n=nch;
           Nev++; 
         }
         fclose(F);
         IC22histRead(); 
      }

      Nev_=simpson(IC22BGdCos,cs_,1,1E-3,NULL);    
if(rand)
{
// generate random number of events in Nev_ \pm 2*sqrt(Nev_) interval
   for(;;)
   { double  c=1;
     Nev= Nev_ + 2*sqrt(Nev_)*(2*drand48()-1); 

     for(i=Nev_;i<=Nev;i++) c*=Nev_/i;
     for(i=Nev_;i>Nev;i--)  c/=Nev_/i;
     if(c>drand48()) break; 
   }

   for(i=0;i<Nev;i++)
   { // printf("IC22BGdCos(0.999)/IC22BGdCos(1)=%E\n", IC22BGdCos(0.999)/IC22BGdCos(1));
     double cs;
     for(;;)
     { cs=cs_+drand48()*(1-cs_);
       if(IC22BGdCos(cs)/IC22BGdCos(1)>drand48())  break; 
     }
     events[i].cs=cs;
   }
/*   
   double xx[20],ex[20];
   for(i=0;i<10;i++)  xx[i]=0;
   for(i=0;i<Nev;i++)
   { int k= 10*(1-events[i].cs)/(1 - cs_);
      xx[k]+=1;
   }
   for(i=0;i<Nev;i++) ex[i]=sqrt(xx[i]);
   displayPlot("", cs_, 1, "cs",10,0,1,xx,ex,"dN/dcos");
*/                               
} 

      double  nu_[NZ], NU_[NZ];
      for(i=0;i<NZ;i++) {nu_[i]=nu[i];NU_[i]=NU[i];} 
      spectrMult(nu_,IC22nuAr);
      spectrMult(NU_,IC22nuBarAr); 
      addSpectrum(nu_,NU_);
      for(i=0;i<NZ;i++) NU_[i]=nu_[i];
      spectrMult(NU_,IC22sigma); spectrMult(NU_,IC22sigma);
      
      double M=nu_[0];
      for(i=1; i<=20;i++)
      {  
         double E1=IC22chan[i].E1,E2=IC22chan[i].E2;
         if(E1<M)
         {
           if(E2>M) E2=M;
           IC22chan[i].n=spectrInt(E1,E2,nu_);
           if(IC22chan[i].n<=0) IC22chan[i].s2=0; else
           IC22chan[i].s2=spectrInt(E1,E2,NU_)/IC22chan[i].n;
         } else {IC22chan[i].n=0; IC22chan[i].s2=0;} 
         IC22chan[i].n*=104/365./1.2; 
      }
      
      double Nbg=simpson(IC22BGdCos,cs_,1,1E-3,NULL);
      double Ns=0;
      for(i=1;i<20;i++) if(IC22chan[i].n>0) Ns+=IC22chan[i].n*(1-exp((cs_-1)/IC22chan[i].s2));      
      double nData=0;
      for(i=0;i<Nev;i++) if(events[i].cs>cs_) nData++;

      int k;
      for(k=0;k<21;k++)
      { double a= 0.1*k;
        L[k]=nData*log(Nbg+a*Ns)-Nbg-a*Ns;
        for(i=0;i<Nev;i++)if(events[i].cs>cs_) L[k]+=log(pp(events[i].cs,events[i].n-10,a)/(Nbg+a*Ns)); 
      }
      
      for(k=1;k<21;k++) L[k]-=L[0];
      L[0]=0;     
      if(L[20]>=L[15]) return 0;

      double dI=exp(L[20])/(2*(L[15]-L[20]));  
      double int1=simpson(P,0.,1.,1E-3,NULL);
      double int2=simpson(P,1.,2.,1E-2,NULL);
      double exLev=int1/(int1+int2+dI);
//      if(prn)
//      {  displayFunc(P,0,2,"PPP");
//         exit(0);
//      }      
      return int1/(int1+int2+dI);
   }
#endif

int  IC22events(double *nu, double * nuB, double phi_cut, double *Nsig,double *Nbg, int*Nobs)
{ double dfi,cs,cs_=cos(phi_cut*M_PI/180);                           
  if(phi_cut>25) { printf(" Too large angle cut\n");  return 1;}   
  if(Nbg) *Nbg=simpson(IC22BGdCos,cs_,1,1E-3,NULL);                        
  if(Nobs)                                                            
  {  FILE*F;                                                          
     char fname[300];                                                 
     int i;                                                           
     sprintf(fname,"%s/sources/data_nu/IC22_events_25.dat",micrO);    
     F=fopen(fname,"r");                                              
     for(*Nobs=0;fscanf(F," %lf %lf %*d",&cs,&dfi)==2; ) if(cs>=cs_)(*Nobs)++;
     fclose(F);                                                       
  }                                                                   
  if(Nsig)                                                            
  { double Mdm,E,s,n[NZ];                                             
    int i;                                                            
    *Nsig=0;                                                          
    if(nu)                                                            
    { Mdm= nu[0];                                                      
      n[0]=Mdm;                                                       
      for(i=1;i<NZ;i++)                                               
      { E= Mdm*exp(Zi(i));                                            
        s=IC22sigma(E); 
//        printf("i=%d E=%e sigma=%e (1-exp((cs_-1)/s))=%e IC22nuAr(E)=%E n[i]=%e\n",
//        i,E,s,1-exp((cs_-1)/s), IC22nuAr(E),n[i]);                                                
        n[i]=nu[i]*IC22nuAr(E)*(1-exp((cs_-1)/(s*s))) ;
      }
      *Nsig+=spectrInt(50,Mdm,n);                                     
    }                                                                 
    if(nuB)                                                           
    { Mdm= nuB[0];                                                     
      n[0]=Mdm;                                                       
      for(i=1;i<NZ;i++)                                               
      { E= Mdm*exp(Zi(i));                                            
        s=IC22sigma(E);                                                
        n[i]=nuB[i]*IC22nuBarAr(E)*(1-exp((cs_-1)/(s*s)));             
      }                                                               
      *Nsig+=spectrInt(50,Mdm,n);                                     
    }
     (*Nsig)*=104./365.; //  Explosure time                                                                 
  }                                                                   
  return 0;                                                           
}                                                                     


double  fluxFactorIC22(double cl,double *NU,double*NUbar)
{ 
  double nu[NZ],nu_[NZ];
  double f=1,cl0; 
  double fmin,fmax,exmin,exmax;
  
  nu[0]=NU[0]; nu_[0]=NUbar[0];  

  cl0=exLevIC22(NU,NUbar, NULL);
  fmax=1;exmax=cl0;
  fmin=1;exmin=cl0; 
  if(cl0>cl) while( exmin>cl) 
  { 
       fmax=fmin;
       exmax=exmin;
       fmin/=2;
       for(int i=1;i<NZ;i++) { nu[i]=NU[i]*fmin; nu_[i]=NUbar[i]*fmin;}
       exmin=exLevIC22(nu,nu_, NULL);
  } else 
  while( exmax<cl)
  { 
       fmin=fmax;
       exmin=exmax;
       fmax*=2;
       for(int i=1;i<NZ;i++) { nu[i]=NU[i]*fmax; nu_[i]=NUbar[i]*fmax;}
       exmax=exLevIC22(nu,nu_, NULL);
  }  
     
  while(exmax-exmin> 0.005*(exmax+exmin) )
  {
     f=0.5*(fmin+fmax);
     for(int i=1;i<NZ;i++) { nu[i]=NU[i]*f; nu_[i]=NUbar[i]*f;}
     cl0=exLevIC22(nu,nu_, NULL);
     if(cl0> cl)
     { fmax=f;exmax=cl0;} else { fmin=f;exmin=cl0;}
  }   
  return 0.5*(fmin+fmax);

}

#ifdef TEST_ANOMALY 
double  fluxFactorIC22_random(double cl,double *NU,double*NUbar)
{ 
  double nu[NZ],nu_[NZ];
  double f=1,cl0; 
  double fmin,fmax,exmin,exmax;

  nu[0]=NU[0]; nu_[0]=NUbar[0];  

  cl0=exLevIC22_random(NU,NUbar, 1);
  fmax=1;exmax=cl0;
  fmin=1;exmin=cl0; 
  if(cl0>cl) while( exmin>cl) 
  { 
       fmax=fmin;
       exmax=exmin;
       fmin/=2;
       for(int i=1;i<NZ;i++) { nu[i]=NU[i]*fmin; nu_[i]=NUbar[i]*fmin;}
       exmin=exLevIC22_random(nu,nu_,0);
  } else 
  while( exmax<cl)
  { 
       fmin=fmax;
       exmin=exmax;
       fmax*=2;
       for(int i=1;i<NZ;i++) { nu[i]=NU[i]*fmax; nu_[i]=NUbar[i]*fmax;}
       exmax=exLevIC22_random(nu,nu_,0);
  }  
     
  while(exmax-exmin> 0.005*(exmax+exmin) )
  {
     f=0.5*(fmin+fmax);
     for(int i=1;i<NZ;i++) { nu[i]=NU[i]*f; nu_[i]=NUbar[i]*f;}
     cl0=exLevIC22_random(nu,nu_,0);
     if(cl0> cl)
     { fmax=f;exmax=cl0;} else { fmin=f;exmin=cl0;}
  }   
  return 0.5*(fmin+fmax);
}

#endif

double fluxFactIC22_FC(double cl,double * nu, double*NU)
{  int i;
   double cs,cs_;
   typedef struct{double cs; int n;} eventStr;
   static  eventStr*  events=NULL;
   static int Nev;
   double fiMax=8; 

   getdNSdCos(nu,NU);

   if(fiMax>10) { printf("Maximum angel reset 10 degrees\n"); fiMax=10;}
   if(!events)
   {  char fname[300];
      double dfi;
      int nch;
      events=malloc(180*sizeof(eventStr));
      sprintf(fname, "%s/sources/data_nu/IC22_events_25.dat",micrO); 
      FILE*F=fopen(fname,"r");
      cs_=cos(10./180.*M_PI);    
      for(Nev=0; fscanf(F," %lf %lf %d",&cs,&dfi, &nch)==3; ) if(cs>=cs_)
      { if(Nev==180) break;
         events[Nev].cs=cs;
         events[Nev].n=nch;
         Nev++; 
      }
      fclose(F);
      IC22histRead();  
   }

   getdNSdCos(nu,NU);
   {  double rate=0,Ns,Nbg,cs_=-1;
      int nData=0;      
      for(cs=1;cs>0.99;cs-=0.0001)       
      { 
         Ns= simpson(dNSdCos,cs,1,1.E-3,NULL)/1.2;
         Nbg=simpson(IC22BGdCos,cs,1,1E-3,NULL);
        if(rate<Ns/sqrt(Nbg)) {rate=Ns/sqrt(Nbg); cs_=cs;}   
      } 
      Ns= simpson(dNSdCos,cs_,1,1.E-3,NULL)/1.2;
      Nbg=simpson(IC22BGdCos,cs_,1,1E-3,NULL);   
      for(i=0;i<Nev;i++) if(events[i].cs>cs_) nData++;
//      printf("cs_=%E nData=%d, Nbg=%E,cl=%E\n",cs_,nData,Nbg, cl);
      return  FeldmanCousins(nData, Nbg, cl)/Ns;
   }      
}
      




double ic22nuar_(double* E)  { return IC22nuAr(*E);  }
double ic22nubarar_(double*E){ return IC22nuBarAr(*E);}

double ic22bgdcos_(double*cs){ return IC22BGdCos(*cs); }
double ic22sigma_(double*E)  { return IC22sigma(*E); }

int  ic22events_(double *nu, double * nuB, double * phi_cut, double *Nsig,double *Nbg, int*Nobs) 
{ return IC22events(nu,nuB, *phi_cut,Nsig,Nbg,Nobs);}
                                     

double exlevic22_(double*nu,double*NU,double*L)
                             { return  exLevIC22(nu,NU,L); }

double fluxfactoric22_( double* pval,double *NU,double*NUbar)
                             { return fluxFactorIC22(*pval,NU,NUbar);} 
