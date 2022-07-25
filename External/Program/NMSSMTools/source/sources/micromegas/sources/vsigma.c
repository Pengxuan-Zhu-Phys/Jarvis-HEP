#include <sys/utsname.h>
#include "micromegas.h"
#include "micromegas_aux.h"
#include "micromegas_f.h"

//======================  vSigmaCC =======================

static double T_;
static double M1,M2,sqrtSmin,sqrtSmax;
static CalcHEP_interface * CI;
static int i3=2,i4=3,i5=4,i6=5;
static REAL pmass[6];
static int pdg[6];

static double sing2(char *s, int nout, double m, double w)
{  int i;
   int nin=2;
   
   double sum0=0,sum1=0;
   if(s[0]<=nin) return 0;
   if(strlen(s)<2) return 0;

   for(i=0;s[i];i++) sum1+=pmass[s[i]-1];
   if(m<=sum1) return 0;
   for(i=nin;i<nin+nout;i++) sum0+=pmass[i];      
   if(m>= sqrtSmax -sum0 +sum1) return 0;   
   return 1/(m*w);     
}

static double s_integrandT_(double  sqrtS )
{  double sv_tot,t,bess, x1,x2,y;
   REAL ms,md,PcmIn;
   double Rm;
   
   ms = M1 + M2; 
   if(ms>=sqrtS)  return 0;
   x1=M1/T_; x2=M2/T_; y=sqrtS/T_;

   if(y-x1-x2>50) return 0;   
      
   md = M1 - M2;
   PcmIn = Sqrt((sqrtS-ms)*(sqrtS+ms)*(sqrtS-md)*(sqrtS+md))/(2*sqrtS);
   kin22(PcmIn,pmass);
   
   sv_tot=simpson(dSigma_dCos,-1.,1.,1E-3,NULL); 
   improveCrossSection(pdg[0],pdg[1],pdg[2],pdg[3],PcmIn,&sv_tot);
   bess=    sqrt(2*x1*x2/y/M_PI)*exp(-(y-x1-x2))*K1pol(1/y)/(K2pol(1/x1)*K2pol(1/x2));
   Rm=PcmIn*sqrtS/M1/M2;   
   return  bess*Rm*Rm*sv_tot/T_;    
}   

/*
bessK2(x) = exp(-x)*sqrt(M_PI/2/x)*K2pol(1/x)
bessK1(x) = exp(-x)*sqrt(M_PI/2/x)*K1pol(1/x) 
*/

static double u_integrand_( double u)
{  double z,y,sv_tot,w;
   double Xf_1;
   double ms,md,sqrtS,PcmIn,res0;
   
   if(u==0. || u==1.) return 0.;
   z=1-u*u;
   sqrtS=M1+M2-3*T_*log(z);
   if(sqrtS<=M1+M2 || sqrtS<=pmass[2]+pmass[3]) return 0;
   return s_integrandT_(sqrtS )*6*T_*u/z;
}

static double vsigma23integrandT(double *x, double w)
{
   double pcmIn,sqrtS,M45;
   double M45_min=pmass[i4]+pmass[i5],M45_max=sqrtSmax-pmass[i3];
   int err;
   double r, x1,x2,y,z, bess,Rm;
   double GG;
   REAL pvect[20];

   z=1-x[0]*x[0];
   sqrtS=M1+M2-3*T_*log(z);
   
   if(sqrtS<=sqrtSmin || sqrtS>=sqrtSmax) return 0;    
   pcmIn=decayPcm(sqrtS,pmass[0], pmass[1]);

   M45=M45_min+x[1]*(M45_max-M45_min);
   
   x1=M1/T_; x2=M2/T_; y=sqrtS/T_;
   if(y-x1-x2>50) return 0;

   r=kinematic_23(pcmIn,i3, M45, 2*x[2]-1 ,2*x[3]-1,M_PI*x[4],pmass,pvect)*8*M_PI*(M45_max-M45_min); //  /pcmIn
   if(r==0) return 0;
   GG=sqrt(4*M_PI*alphaQCD(sqrtS));
   r*= CI->sqme(1,GG, pvect,NULL,&err);
   bess= sqrt(2*x1*x2/y/M_PI)*exp(-(y-x1-x2))*K1pol(1/y)/(K2pol(1/x1)*K2pol(1/x2));
//   bess=bessK1(sqrtS/T_)/bessK2(M1/T_)/bessK2(M2/T_);

   Rm=sqrtS/M1/M2;   
   
   r*= pcmIn*bess*Rm*Rm*6*x[0]/z;
   return r; 
}

static double vsigma23integrand0(double *x, double w)
{
   double r,sqrtS=M1+M2, M45_min=pmass[i4]+pmass[i5],M45_max=M1+M2-pmass[i3];
   int err;
   double GG;
   REAL pvect[20];
   
   r=kinematic_23(0.,i3,M45_min+x[0]*(M45_max-M45_min),0.5,2*x[1]-1,M_PI*x[2],pmass,pvect)*8*M_PI*(M45_max-M45_min);
   if(r==0) return 0;
   GG=sqrt(4*M_PI*alphaQCD(sqrtS));
   r*= CI->sqme(1,GG, pvect,NULL,&err);
   r*= (M1+M2)/M1/M2;
//printf("r=%e\n",r);   
   return r; 
}


static double vsigma24integrandT(double *x, double w)
{
   double pcmIn,sqrtS,M34,M56;
   double M34_min=pmass[i3]+pmass[i4],M34_max=sqrtSmax-pmass[i5]-pmass[i6];
   double M56_min=pmass[i5]+pmass[i6],M56_max=sqrtSmax-pmass[i3]-pmass[i4];
   int err;
   double r, x1,x2,y,z, bess,Rm;
   double GG;
   REAL pvect[24];

   z=1-x[0]*x[0];
   sqrtS=M1+M2-3*T_*log(z);
   
   if(sqrtS<=sqrtSmin || sqrtS>=sqrtSmax) return 0;    
   pcmIn=decayPcm(sqrtS,pmass[0], pmass[1]);

   M34=M34_min+x[1]*(M34_max-M34_min);
   M56=M56_min+x[2]*(M56_max-M56_min);
   
   x1=M1/T_; x2=M2/T_; y=sqrtS/T_;
   if(y-x1-x2>50) return 0;

   r=kinematic_24(pcmIn,i3,i4, M34, M56,  2*x[3]-1 ,2*x[4]-1,2*x[5]-1, 2*M_PI*x[6], 2*M_PI*x[7],pmass,pvect)
                  *(M34_max-M34_min)*(M56_max-M56_min)*2*4*M_PI*4*M_PI; //  /pcmIn
   GG=sqrt(4*M_PI*alphaQCD(sqrtS));
   r*= CI->sqme(1,GG, pvect,NULL,&err);
   bess=    sqrt(2*x1*x2/y/M_PI)*exp(-(y-x1-x2))*K1pol(1/y)/(K2pol(1/x1)*K2pol(1/x2));
//   bess=bessK1(sqrtS/T_)/bessK2(M1/T_)/bessK2(M2/T_);

   Rm=sqrtS/M1/M2;   
   
   r*= pcmIn*bess*Rm*Rm*6*x[0]/z;
   return r; 
}


static double vsigma24integrand0(double *x, double w)
{
   double r,sqrtS=M1+M2,M34,M56;
   double  M34_min=pmass[i3]+pmass[i4],M34_max=M1+M2-pmass[i5]-pmass[i6],
           M56_min=pmass[i5]+pmass[i6],M56_max=M1+M2-pmass[i3]-pmass[i4];
   int err;
   double GG;
   REAL pvect[24];
   
   M34=M34_min+x[0]*(M34_max-M34_min);
   M56=M56_min+x[1]*(M56_max-M56_min);
   
   r=kinematic_24(0., i3,i4,M34, M56, 0.5 ,2*x[2]-1,2*x[3]-1, 2*M_PI*x[4], 2*M_PI*x[5],pmass,pvect)
                     *(M34_max-M34_min)*(M56_max-M56_min)*2*4*M_PI*4*M_PI; 
   
   GG=sqrt(4*M_PI*alphaQCD(sqrtS));
   r*= CI->sqme(1,GG, pvect,NULL,&err);
   r*= (M1+M2)/M1/M2;
//printf("r=%e\n",r);   
   return r; 
}


double vSigmaCC(double T,numout* cc, int mode)
{
  int i,err,n,n0,m,w;
  char*s, *pname[6];
  double msum;
  double a=0,factor,dMax=0;
  int spin2,cdim,neutral1,neutral2;
  double oldQ;
  
  double bEps=1.E-4;
  double dI;
  int dmOut;
      
  CI=cc->interface; 
  T_=T;
  


  if(passParameters(cc)) return -1;
  if(Qaddress && CI->nout==2) 
  {  oldQ=*Qaddress;
     for(i=0;i<2;i++) pname[i]=CI->pinf(1,i+1,pmass+i,pdg+i);
     *Qaddress=pmass[0]+pmass[1];
     calcMainFunc();
     if(passParameters(cc)) return -1;
  }
  
  for(i=0;i<2+CI->nout;i++) pname[i]=CI->pinf(1,i+1,pmass+i,pdg+i);  

  M1=pmass[0];
  M2=pmass[1];
  

  if(mode) 
  { if(pname[0][0]!='~' || pname[1][0]!='~') return 0;
    if(T==0 && (M1< Mcdm || M2<Mcdm)) return 0;
    dmOut=0; 
    for(i=2;i<2+CI->nout;i++) if(pname[i][0]=='~') dmOut++;
    if(dmOut==2) return 0; 
  }    
 
  for(i=2,msum=0;i<CI->nout;i++) msum+=pmass[i];
  
  sqrtSmin=M1+M2;
  
  if(msum > sqrtSmin)
  { if(T==0) return 0; else sqrtSmin=msum; }
  sqrtSmax=sqrtSmin-T*log(bEps); 

  n0=0; 
  if(CI->nout>2) for(n=1;(s=CI->den_info(1,n,&m,&w,NULL));n++)
  { double mm=0, ww=0;
    if(m) mm=fabs(CI->va[m]); if(w) ww=CI->va[w];
    double d=sing2(s,CI->nout,mm,ww); 
    if(!isfinite(d)) { printf("non-integrable pole\n"); return 0;}
    if(d>dMax){ dMax=d; n0=n;} 
  }

  switch(CI->nout)
  { 
     case 2:
       if(T==0) a=vcs22(cc,1,&err); else
       {  double eps=1.E-3;
          sqme22=CI->sqme;
          nsub22=1;
          a=simpson(u_integrand_,0.,1.,eps,NULL)*3.8937966E8;
       }   
       break;
     case 3:
     {  
        if(n0)
        {  s=CI->den_info(1,n0,&m,&w,NULL);
           for(i3=2;i3<5;i3++) if(i3!=s[0]-1 && i3!=s[1]-1) break;
           for(i4=2;i4<4;i4++) if(i4!=i3) break;   
           for(i5=i4+1;i5<=4;i5++) if(i5!=i3) break;    
        } else {i3=2;i4=3;i5=4;}
   
        
        if(T==0) a=vegas_chain(3, vsigma23integrand0 ,2000,1., 0.03,&dI,NULL);
        else     a=vegas_chain(5, vsigma23integrandT ,2000,1., 0.03,&dI,NULL);
        break;
     }
     case 4:
         if(n0) 
         {  s=CI->den_info(1,n0,&m,&w,NULL);
            i3=s[0]-1;
            i4=s[1]-1;
            for(i5=2;i5<5;i5++)  if(i5!=i3 && i5!=i4) break;
            for(i6=i5+1;i6<=5;i6++) if(i6!=i3 && i6!=i4) break;
         }else { i3=2;i4=3;i5=4;i6=5;} 
//         printf("i3,i4,i5,i6= %d %d %d %d\n", i3,i4,i5,i6); 
         if(T==0) a=vegas_chain(6, vsigma24integrand0 ,4000,1., 0.03,&dI,NULL);
         else     a=vegas_chain(8, vsigma24integrandT ,20000,1., 0.001,&dI,NULL);
                                
     break;
     default:
        printf("Too many outgoing particles\n");
        a=0;  
   }  

//   WIDTH_FOR_OMEGA=0;
  
   if(mode)
   {  a*= 1-0.5*dmOut;
      char *p0=pname[0],*p1=pname[1],*c0=NULL,*c1=NULL;
      double s=0,k=M1*M1*M2*M2/sqrt(M1*M2)*K2pol(T/M1)*K2pol(T/M2);            
      for(i=0;i<nModelParticles;i++) if(ModelPrtcls[i].name[0]=='~')
      {  double m=pMass(ModelPrtcls[i].name);
         int dim=ModelPrtcls[i].cdim*(ModelPrtcls[i].spin2+1);
         
              if(strcmp(p0,ModelPrtcls[i].name) ==0){ c0=ModelPrtcls[i].aname;k*=dim;} 
         else if(strcmp(p0,ModelPrtcls[i].aname)==0){ c0=ModelPrtcls[i].name; k*=dim;}
              if(strcmp(p1,ModelPrtcls[i].name) ==0){ c1=ModelPrtcls[i].aname;k*=dim;} 
         else if(strcmp(p1,ModelPrtcls[i].aname)==0){ c1=ModelPrtcls[i].name; k*=dim;}
         if(strcmp(ModelPrtcls[i].name,ModelPrtcls[i].aname)!=0) dim*=2;
         if(0.5*(M1+M2)-m >30*T){ k=0;s=1;break;}
         s+=dim*m*m/sqrt(m)*K2pol(T/m)*exp(-(m-0.5*(M1+M2))/T); 
      }
      if(k)      
      {  if(strcmp(p0,p1)) k*=2;
         if(! (  (strcmp(p0,c0)==0 && strcmp(p1,c1)==0)
               ||(strcmp(p0,c1)==0 && strcmp(p1,c0)==0)
              ) ) k*=2;   
      } 
      a*=k/s/s;
   }  

   if(Qaddress && CI->nout==2) { *Qaddress=oldQ; calcMainFunc();}
   
   return a;
}

double vsigmacc_(double *T, int *ccf,int*mode)
{
  numout*cc;
  memcpy(&cc,ccf,sizeof(cc));
  return vSigmaCC(*T,cc,*mode); 
}
