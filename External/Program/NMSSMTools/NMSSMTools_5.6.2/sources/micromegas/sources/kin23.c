#include"micromegas.h"
#include"micromegas_aux.h"

//#include"../CalcHEP_src/c_source/ntools/include/vegas.h"
//#include"../CalcHEP_src/c_source/ntools/include/simpson.h"

//#include<math.h>

static REAL pmass[6]={10,20,30,40,50};
static int i3=2,i4=3,i5=4,i6=5;
static numout* cc23_=NULL,*cc24_=NULL;
static double eps=1.E-2;
static int nsub_;

static double MZ=91.1884 ,wZ=2.4944,K=1;
static double GG=1.23;

double  kinematic_23(double Pcm,int i3, double M45, double cs1, double cs2,  double fi,REAL*pmass, REAL*pvect)
{
  int i;
  double pcm1,pcm2,p0,p2,p3,chY,shY,sn1,sn2,r;
  int err,i4,i5;
  
  switch(i3)
  { 
    case 3: i4=2;i5=4; break;
    case 4: i4=2;i5=3; break;
    default: i3=2;i4=3;i5=4;
  }
  
  sn1=Sqrt(1-cs1*cs1);
  sn2=Sqrt(1-cs2*cs2);
  for(i=0;i<20;i++)pvect[i]=0;

  pvect[0]= Sqrt(Pcm*Pcm+pmass[0]*pmass[0]);
  pvect[3]= Pcm;
  pvect[4]= Sqrt(Pcm*Pcm+pmass[1]*pmass[1]);
  pvect[7]=-Pcm;
  
  if(pvect[0]+pvect[4]<=M45+pmass[i3]) return 0;
  pcm1=decayPcm(pvect[0]+pvect[4],pmass[i3], M45);
//printf("pcm1=%E\n",pcm1);  
  if(!pcm1) return 0;
  pvect[i3*4]=Sqrt(pmass[i3]*pmass[i3]+pcm1*pcm1);
  pvect[i3*4+3]=pcm1*cs1;
  pvect[i3*4+2]=pcm1*sn1;
  
  pcm2=decayPcm(M45,pmass[i4],pmass[i5]);
//printf("pcm2=%E  M45=%E  pmass[%d]=%E pmass[%d]=%E \n",pcm2,M45,i4,pmass[i4],i5,pmass[i5]);

  if(!pcm2) return 0;
  
  pvect[i4*4]=Sqrt(pcm2*pcm2+pmass[i4]*pmass[i4]);
  pvect[i4*4+3]=pcm2*cs2;
  pvect[i4*4+2]=pcm2*sn2*cos(fi);
  pvect[i4*4+1]=pcm2*sn2*sin(fi);

//printf("fi=%E pcm2=%E sn2=%E\n",fi, pcm2,sn2);
  
  pvect[i5*4]=Sqrt(pcm2*pcm2+pmass[i5]*pmass[i5]);
  pvect[i5*4+3]=-pvect[i4*4+3];
  pvect[i5*4+2]=-pvect[i4*4+2];    
  pvect[i5*4+1]=-pvect[i4*4+1]; 
 
  shY=pcm1/M45;
  chY=Sqrt(1+shY*shY);
    
  p0=pvect[i4*4]; p3=pvect[i4*4+3];
  pvect[i4*4  ]=p0*chY+p3*shY;
  pvect[i4*4+3]=p0*shY+p3*chY;
   
  p0=pvect[i5*4]; p3=pvect[i5*4+3];
  pvect[i5*4  ]=p0*chY+p3*shY;
  pvect[i5*4+3]=p0*shY+p3*chY;
  
  p2=pvect[i4*4+2];p3=pvect[i4*4+3];
  pvect[i4*4+2]=-(p2*cs1+p3*sn1);
  pvect[i4*4+3]=-(-p2*sn1+p3*cs1);

  p2=pvect[i5*4+2];p3=pvect[i5*4+3];
  pvect[i5*4+2]= -(p2*cs1+p3*sn1);
  pvect[i5*4+3]= -(-p2*sn1+p3*cs1);

  return 3.8937966E8*pcm1*pcm2/(512*M_PI*M_PI*M_PI*M_PI*(pvect[0]+pvect[4])*(pvect[0]+pvect[4]));
  return r; 
}  

int NCALL;

static double Pcm_, Mmax, Mmin, cs1_, cs2_, fi_,M_;
static int Np=5;

static double Idfi(double cs2) 
{ int i,err;
  double s,factor;
  REAL pvect[20];
  for(s=0,i=0;i<Np;i++)
  {  double fi=M_PI*(0.5+i)/Np,ds;
     factor=kinematic_23(Pcm_,i3, M_, cs1_, cs2,fi,pmass,pvect)/Pcm_;
     if(factor==0.) return 0;   
     ds=factor*(cc23_->interface->sqme)(nsub_,GG, pvect,NULL,&err);
     s+=ds;
  } 
  return 2*M_PI*s/Np;  
}

static double dIdcs2(double cs1) { cs1_=cs1; return gauss( Idfi ,-1.,1.,Np); }
static double dIdM(double M) { NCALL++; M_=M; return gauss( dIdcs2 ,-1.,1.,Np); }

static double q=1./2.;
static double  m0,w0;
static double yMin,yMax;

static double q2=1./2.;
static double zMin,zMax,uMax;



static double dIdM_mu(double mu)
{  if(mu==0) return 0;
   return 1./q *dIdM(Mmax - pow(mu,1/q))*pow(mu,(1-q)/q); 
}

static double dIdM_w(double y)
{ double m=sqrt(m0*(m0+w0*tan(y)));
  return dIdM(m)* m0*w0*(1+tan(y)*tan(y))/2/m;
}


static double dIdM_w_mu(double mu)
{  return 1./q *dIdM_w(yMax - pow(mu,1/q))*pow(mu,(1-q)/q); }

static double dIdM_w_mu_mu(double mu)
{  if(mu<=0 || mu>=uMax) return 0;
   return 1./q2 *dIdM_w_mu(zMax - pow(mu,1/q2))*pow(mu,(1-q2)/q2); 
}

static double dIdM_mu_mu(double mu)
{  return 1./q2 *dIdM_w_mu(zMax - pow(mu,1/q2))*pow(mu,(1-q2)/q2); }


double cs23(numout*cc, int nsub, double Pcm, int ii3,int*err) 
{ int i,n;
  char*s,s0[3];
  int m,w;
  double muMax;

  if(err) *err=0;

  passParameters(cc);
  GG=sqrt(4*M_PI*alphaQCD(GGscale));
  Pcm_=Pcm;
  nsub_=nsub;
  cc23_=cc;

  switch(ii3)
  {  
    case 3: i3=3; i4=2;i5=4; break;
    case 4: i3=4; i4=2;i5=3; break;
    default:i3=2; i4=3;i5=4; 
  } 
                       
  for(i=0;i<5;i++) cc->interface->pinf(nsub,1+i,pmass+i,NULL);

  s0[0]=i4+1;
  s0[1]=i5+1;
  s0[2]=0;

  for(n=1;;n++)
  {  s=cc->interface->den_info(nsub,n,&m,&w,NULL);
     if(s==NULL) break;
     if(strcmp(s,s0)==0) 
     { 
//       printf("%E %E \n", cc->interface->va[m], cc->interface->va[w]); 
       
      if(m) m0=cc->interface->va[m]; else m0=0;
      if(w) w0=cc->interface->va[w]; else w0=0;
       break; 
     } 
  }
     
  Mmax=Sqrt(Pcm*Pcm+pmass[0]*pmass[0])+Sqrt(Pcm*Pcm+pmass[1]*pmass[1])-pmass[i3];
  Mmin=pmass[i4]+pmass[i5];

  if(Mmin>=Mmax) { return 0;}
  muMax=pow(Mmax-Mmin,q);
  
  if(s==NULL)
  { int err_;
    double r= simpson(dIdM_mu, 0., muMax*0.99999, 0.01,&err_);
    if(err_) { if(err) *err=err_;  printf("error in simpson kin23.c line 189\n");}
    return  r;
  }

  yMin=atan((Mmin*Mmin-m0*m0)/(m0*w0));
  yMax=atan((Mmax*Mmax-m0*m0)/(m0*w0)); 
  
//  displayFunc(dIdM_w,yMin+0.001, yMax-0.001,"dIdM_mu");

  zMax=pow(yMax-yMin,q);
  
  uMax=pow(zMax,q2);
  
//  displayFunc(dIdM_w_mu_mu, 0., uMax,"dIdM_w_mu_mu");
  { int err_; 
    double r=simpson(dIdM_w_mu_mu,0, uMax, 0.01,&err_);
    if(err_){ if(err) *err=err_; printf("error in simpson cs23.c line 206\n");}
    return r;
  }
}


static double func_23(double *x, double wgt)
{ int err; 
  REAL pvect[20];
  double  factor;
  factor=kinematic_23(Pcm_, i3, Mmin+x[0]*(Mmax-Mmin),2*(x[1]-0.5) ,2*(x[2]-0.5), M_PI*x[3],pmass,pvect)/Pcm_;
  if(factor==0) return 0; 
  return  factor*2*(Mmax-Mmin)*4*M_PI*(cc23_->interface->sqme)(nsub_,GG, pvect,NULL,&err);
}

double cs23Vegas(numout * cc, int nsub, double Pcm, int ii3, 
    int nSess1, int nCall1,  int nSess2, int nCall2, 
    double*dcs, double *chi2, int * err) 
{
  int i,k, nCall[2]={nCall1,nCall2}, nSess[2]={nSess1,nSess2};
  double rVal=0; 
  vegasGrid *vegPtr=NULL;
  GG=sqrt(4*M_PI*alphaQCD(GGscale));  
  
  if(passParameters(cc)) return 0;  

  Pcm_=Pcm;
  cc23_=cc;  
  if(cc->interface->nin!=2 || cc->interface->nout!=3) {*err=-1; return 0;}
  i3=ii3;
  switch(ii3)
  { case 2: i4=3;i5=4; break; 
    case 3: i4=2;i5=4; break;
    case 4: i4=2;i5=3; break;
    default: i3=2;i4=3;i5=4; 
  } 
                    
  for(i=0;i<5;i++) cc->interface->pinf(nsub,1+i,pmass+i,NULL);  
  nsub_=nsub;
  
  Mmax=Sqrt(Pcm*Pcm+pmass[0]*pmass[0])+Sqrt(Pcm*Pcm+pmass[1]*pmass[1])-pmass[i3];
  Mmin=pmass[i4]+pmass[i5];

  vegPtr=vegas_init(4,func_23,50);
  for(k=0;k<2;k++)
  if(nSess[k] && nCall[k]) 
  { double ti,dti,s0=0,s1=0,s2=0;
     
    for(i=0;i<nSess[k];i++)
    {  *err=vegas_int(vegPtr, nCall[k], 1.5,  nPROCSS, &ti, &dti); 
        printf("ti=%E, dti=%E\n", ti,dti);    
        if(*err) { vegas_finish(vegPtr);return 0;}
        dti=1/(dti*dti);                                                            
        s0+=dti;                                                         
        s1+=ti*dti;
        s2+=ti*ti*dti;
    }
    rVal=s1/s0;
    *dcs=1/sqrt(s0);
    if(nSess[k]<=1) *chi2=0; else *chi2=(s2-s1*s1/s0)/(nSess[k]-1);
  }
  vegas_finish(vegPtr);
  return rVal;
}

double kinematic_24(double Pcm, int i3,int i4, double M1, double M2, double xcos,double xcos1, double xcos2, double fi1, double fi2,
                      REAL*pmass, REAL * P)
{ 
  double factor,pcm,p1cm,p2cm,chY,shY,xsin,sqrtS;
  double p0,p3;
  int i,j,i5,i6;  

  for(i=0;i<24;i++)P[i]=0;
  for(i5=2;i5<6;i5++)  if(i5!=i3 && i5!=i4) break;
  for(i6=i5+1;i6<6;i6++)  if(i6!=i3 && i6!=i4) break;  
  
  P[0]= Sqrt(Pcm*Pcm+pmass[0]*pmass[0]);
  P[3]= Pcm;
  P[4]= Sqrt(Pcm*Pcm+pmass[1]*pmass[1]);
  P[7]=-Pcm;
  sqrtS=P[0]+P[4];
  
  if(sqrtS<=M1+M2) return 0;
  if(sqrtS <= pmass[2]+pmass[3]+pmass[4]+pmass[5]) return 0;
  
  pcm=decayPcm(sqrtS,M1,M2);             
  
//printf("decayPcm(%E %E %E) ? %E\n", sqrtS,M1,M2, 1- sqrt(M1*M1+pcm*pcm)/sqrtS- sqrt(M2*M2+pcm*pcm)/sqrtS);  
  
  p1cm=decayPcm(M1,pmass[i3],pmass[i4]);     
  p2cm=decayPcm(M2,pmass[i5],pmass[i6]);
    
  P[i3*4+0]=Sqrt(pmass[i3]*pmass[i3]+p1cm*p1cm);    P[i4*4+0]=Sqrt(pmass[i4]*pmass[i4]+p1cm*p1cm);
  P[i3*4+2]=p1cm*Sqrt(1-xcos1*xcos1);
  P[i3*4+1]=sin(fi1)*P[i3*4+2];                     P[i4*4+1]=-P[i3*4+1];                                      
  P[i3*4+2]*=cos(fi1);                              P[i4*4+2]=-P[i3*4+2];
  P[i3*4+3]=p1cm*xcos1;                             P[i4*4+3]=-P[i3*4+3];

  P[i5*4+0]=Sqrt(pmass[i5]*pmass[i5]+p2cm*p2cm);    P[i6*4+0]=Sqrt(pmass[i6]*pmass[i6]+p2cm*p2cm);
  P[i5*4+2]=p2cm*Sqrt(1-xcos2*xcos2);
  P[i5*4+1]=sin(fi2)*P[i5*4+2];                     P[i6*4+1]=-P[i5*4+1];                                      
  P[i5*4+2]*=cos(fi2);                              P[i6*4+2]=-P[i5*4+2];
  P[i5*4+3]=-p2cm*xcos2;                            P[i6*4+3]=-P[i5*4+3];

//printf("dM1=%E %E %E %E \n", 1-(P[i3*4+0]+P[i4*4+0])/M1, P[i3*4+1]+ P[i4*4+1], P[i3*4+2]+ P[i4*4+2],P[i3*4+3]+ P[i4*4+3]);
//printf("dM2=%E %E %E %E\n",  1-(P[i5*4+0]+P[i6*4+0])/M2, P[i5*4+1]+ P[i6*4+1], P[i5*4+2]+ P[i6*4+2],P[i5*4+3]+ P[i6*4+3]);
      
  shY=pcm/M1;
  chY=Sqrt(1+shY*shY);
  for(i=2;i<6;i++) if(i==i3||i==i4)
  { double  p0=P[4*i];
    double  p3=P[4*i+3];
    P[4*i]=chY*p0+shY*p3;
    P[4*i+3]=shY*p0+chY*p3;
  }

  shY=-pcm/M2;            
  chY=Sqrt(1+shY*shY);  
  for(i=2;i<6;i++) if(i==i5||i==i6)
  { double  p0=P[4*i];
    double  p3=P[4*i+3];
    P[4*i]=chY*p0+shY*p3;
    P[4*i+3]=shY*p0+chY*p3;
  }

//printf("Energy conservation\n");
for(i=0;i<4;i++)
{ double sum=P[0+i]+P[4+i]-P[8+i]-P[12+i]-P[16+i]-P[20+i];
  if(fabs(sum/P[0]) > 1.E-4)
  { printf(" %E %E -> %E %E %E %E\n",P[0+i],P[4+i],P[8+i],P[12+i],P[16+i],P[20+i]);
    printf("No Energy conservation %E i=%d  \n",sum/P[0],i);    
    exit(22);
  }  
}    

  for(i=0;i<6;i++)
  { double m;
    m=Sqrt(fabs(P[4*i]*P[4*i]-P[4*i+1]*P[4*i+1]-P[4*i+2]*P[4*i+2]-P[4*i+3]*P[4*i+3]));
    if(fabs(m-pmass[i])>pmass[0]*1.E-5) { printf("wrong mass %d (%E != %E) \n",i,m,pmass[i]); exit(33);}
  }



  xsin=Sqrt(1-xcos*xcos);
  for(i=2;i<6;i++) 
  { 
    double  p2=P[4*i+2];     
    double  p3=P[4*i+3]; 
    P[4*i+2] = xcos*p2+xsin*p3;
    P[4*i+3] =-xsin*p2+xcos*p3;
  }  
    


  return pcm*p1cm*p2cm*3.8937966E8/pow(2*M_PI,12-5)/pow(2,4+2)/sqrtS/sqrtS ;
}

static double func_24(double *x, double wgt)
{ int err; 
  REAL pvect[24];
  double  factor,sqrtS,M1,M2;
  sqrtS=Sqrt(Pcm_*Pcm_ +pmass[0]*pmass[0])+Sqrt(Pcm_*Pcm_ +pmass[1]*pmass[1]);
  M1= (1-x[0])*(pmass[i3]+pmass[i4])+ x[0]*(sqrtS-pmass[i5]-pmass[i6]);
  factor=sqrtS-pmass[i3]-pmass[i4] -pmass[i5]-pmass[i6];
  M2= (1-x[1])*(pmass[i5]+pmass[i6])+ x[1]*(sqrtS-M1);
  factor*=sqrtS-M1-pmass[i5]-pmass[i6];
  
  factor*=2*2*2*(2*M_PI)*(2*M_PI)*kinematic_24(Pcm_,i3,i4,M1,M2,1-2*x[2],1-2*x[3],1-2*x[4],2*M_PI*x[5],2*M_PI*x[6],pmass,pvect);
                             
  if(factor==0) return 0; 
  return  factor*(cc24_->interface->sqme)(nsub_,GG, pvect,NULL,&err);
}


double cs24Vegas(numout * cc, int nsub, double Pcm, int ii3, int ii4, 
    int nSess1, int nCall1,  int nSess2, int nCall2, 
    double*dcs, double *chi2, int * err) 
{
  int i,k, nCall[2]={nCall1,nCall2}, nSess[2]={nSess1,nSess2}, nOut=4;
  double rVal=0,C,sqrtS;
   
  vegasGrid *vegPtr=NULL;
  GG=sqrt(4*M_PI*alphaQCD(GGscale));  
//  link_process(cc->interface);
  Pcm_=Pcm;
  cc24_=cc;  
  if(cc->interface->nin!=2 || cc->interface->nout!=4) {*err=-1; return 0;}
  
  i3=ii3;
  i4=ii4;
  for(i5=2;i5==i3||i5==i4;i5++) continue;
  for(i6=2;i6==i3||i6==i4||i6==i5;i6++) continue;
// printf("i3=%d i4=%d i5=%d i6=%d\n",i3,i4,i5,i6);
                    
  for(i=0;i<6;i++) cc->interface->pinf(nsub,1+i,pmass+i,NULL);  
  nsub_=nsub;
  sqrtS=Sqrt(Pcm*Pcm+pmass[0]*pmass[0])+Sqrt(Pcm*Pcm+pmass[1]*pmass[1]);
  *err=0;

  vegPtr=vegas_init(7,func_24,50);
  for(k=0;k<2;k++)
  if(nSess[k] && nCall[k]) 
  { double ti,dti,s0=0,s1=0,s2=0;
     
    for(i=0;i<nSess[k];i++)
    {  *err=vegas_int(vegPtr, nCall[k], 1.5,  nPROCSS, &ti, &dti); 
        if(*err) { vegas_finish(vegPtr);return 0;}
        dti=1/(dti*dti);                                                            
        s0+=dti;                                                         
        s1+=ti*dti;
        s2+=ti*ti*dti;
    }
    rVal=s1/s0;
    *dcs=1/sqrt(s0);
    if(nSess[k]<=1) *chi2=0; else *chi2=(s2-s1*s1/s0)/(nSess[k]-1);
  }
  vegas_finish(vegPtr);
  
  C=3.8937966E8* pow(2*M_PI,4)
  /(4*Pcm*sqrtS*sqrtS)
  /pow(2*pow(2*M_PI,3),nOut) * (2*M_PI);   
  (*dcs)*=C;
  rVal*=C;
  return rVal;
}
