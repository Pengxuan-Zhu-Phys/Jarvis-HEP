#include "micromegas.h"
#include "micromegas_aux.h"

//====================   g1eff  & h1eff ================

typedef   struct{ double mu,eta;}    geff_struct;


static double g1eff_int(double x,void *par_)
{  if(x==0 || x==1) return 0;
   geff_struct*par=par_; 
   double p=2*log(x);
   double e=sqrt(p*p+par->mu*par->mu);
   double den;
   
   if(par->eta==-1 && e<1E-8) den=e*(1+e/2); else den=exp(e)+par->eta;
   
   return   e*p*p/den/x/M_PI/M_PI;
}


static double h1eff_int(double x,void *par_)
{  if(x==0 || x==1) return 0;
   geff_struct*par=par_; 
   double p=2*log(x);
   double e=sqrt(p*p+par->mu*par->mu);
   double den;
   
   if(par->eta==-1 && e<1E-8) den=e*(1+e/2); else den=exp(e)+par->eta;
   return  (par->mu*par->mu+4./3.*p*p)*p*p/den/x/M_PI/M_PI/e;
}


double g1eff(double mu,int eta)
{ 
  geff_struct par;
  par.eta=eta;
  par.mu=mu;
  return simpson_arg(g1eff_int,&par, 0, 1,1E-3,NULL)/(M_PI*M_PI/30);
}   


double h1eff(double mu,int eta)
{ 
  geff_struct par;
  par.eta=eta;
  par.mu=mu;
  return simpson_arg(h1eff_int,&par, 0, 1,1E-3,NULL)/(2*M_PI*M_PI/45);
}   


double Hubble(double T) { return sqrt(8*M_PI/3.*M_PI*M_PI/30.*gEff(T))*T*T/MPlanck;}


static int Tdim=0;
static double *t_=NULL, *heff_=NULL, *geff_=NULL, *s3_=NULL,
 *heff_ln_diff_=NULL;

int loadHeffGeff(char*fname)
{  double *tabs[3];
   int nRec,nCol,i;

   if(fname==NULL)
   {
     char*fName=malloc(strlen(micrO)+40);
     sprintf(fName,"%s/sources/data/%s",micrO,"std_thg.tab");
     int err=loadHeffGeff(fName);
     free(fName);
     return err;     
   }
   nRec=readTable(fname,&nCol,tabs);
   if(nRec<=0) return nRec;
   if(nCol!=3)
   { for(i=0;i<nCol;i++) free(tabs[i]);
     return 0;
   }
   if(t_)free(t_);  if(heff_)free(heff_);   if(geff_)free(geff_);  if(s3_)free(s3_);
   if(heff_ln_diff_) free(heff_ln_diff_);
   
   t_=tabs[0];
   heff_=tabs[1];
   geff_=tabs[2];
   s3_=malloc(nRec*sizeof(double));
   for(i=0;i<nRec;i++) s3_[i]=t_[i]*pow(2*M_PI*M_PI/45.*heff_[i],0.3333333333333333);  
   heff_ln_diff_=malloc(nRec*sizeof(double));
   polintDiff(nRec,t_, heff_,  heff_ln_diff_);
   for(i=0;i<nRec;i++) heff_ln_diff_[i]*=t_[i]/heff_[i]; 
   Tdim=nRec;
   return nRec;
}

double gEff(double T) 
{
  if(Tdim==0) loadHeffGeff(NULL);    
  if(T< t_[0]) T=t_[0];
  if(T> t_[Tdim-1]) T=t_[Tdim-1];
  return polint1(T,Tdim,t_,geff_);
}

double hEffLnDiff(double T)
{ 
  if(Tdim==0) loadHeffGeff(NULL);
  if(T<t_[0]) return heff_ln_diff_[0];
  if(T>t_[Tdim-1]) return 0;
  return polint3(T,Tdim,t_,heff_ln_diff_);
}
         

double hEff(double T) 
{
  if(Tdim==0) loadHeffGeff(NULL); 
  if(T< t_[0]) T=t_[0];
  if(T> t_[Tdim-1]) T=t_[Tdim-1];
  return polint1(T,Tdim,t_,heff_);
} 

double T_s3(double s3) {if(s3>s3_[Tdim-1]) return t_[Tdim-1]*(s3/s3_[Tdim-1]); else  return polint3(s3,Tdim,s3_,t_);} 
double s3_T(double T)  { if(T> t_[Tdim-1]) return s3_[Tdim-1]*T/t_[Tdim-1];    else  return polint3(T,Tdim,t_,s3_);}

