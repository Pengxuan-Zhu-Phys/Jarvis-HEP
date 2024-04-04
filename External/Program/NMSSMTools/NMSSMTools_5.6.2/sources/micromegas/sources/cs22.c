#include "micromegas.h"
#include "micromegas_aux.h"
#include "micromegas_f.h"
#include"../CalcHEP_src/c_source/ntools/include/vegas.h"

double (*sqme22)(int nsub, double GG, REAL *pvect, REAL*cb_coeff, int * err_code)=NULL; 
int  nsub22=0;

/*===========================================================*/
static double Q_ren,Q_fact;
static double totcoef;
static REAL PcmOut;
static REAL pvect[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
static  int PC[5];  
static  int chan=0; 


static double eps=0.0001;

double GGscale=91.187;
/*
double  decayPcm(double am0,  double  am1,  double  am2)
{
  double  summ, diffm, pout;
  summ = am1 + am2;
  diffm = am1 - am2;
  if(am0<summ) return 0;
  return sqrt((am0-summ)*(am0+ summ)*(am0-diffm)*(am0+diffm))/(am0*2);
}
*/          

int  kin22(double PcmIn,REAL * pmass)
{  
   REAL sqrtS;
   int i;
   for(i=0;i<16;i++) pvect[i]=0;
   sqrtS=Sqrt(pmass[0]*pmass[0]+PcmIn*PcmIn)+Sqrt(pmass[1]*pmass[1]+PcmIn*PcmIn);
   PcmOut = decayPcm(sqrtS,pmass[2],pmass[3]);
//printf(" PcmOut =%E (%E %E %E) \n", PcmOut,sqrtS,pmass[2],pmass[3] );   
   if(PcmOut<sqrtS*1.E-4) return 1;
   totcoef = PcmOut /(32.0*M_PI*PcmIn*sqrtS*sqrtS);
   pvect[3] = PcmIn;
   pvect[7] =-PcmIn;
   pvect[0] = Sqrt(PcmIn*PcmIn   + pmass[0]*pmass[0]);
   pvect[4] = Sqrt(PcmIn*PcmIn   + pmass[1]*pmass[1]);
   pvect[8] = Sqrt(PcmOut*PcmOut + pmass[2]*pmass[2]);
   pvect[12]= Sqrt(PcmOut*PcmOut + pmass[3]*pmass[3]);

   return 0;
}

double  dSigma_dCos(double  cos_f)
{
   double  r;
   REAL sin_f=Sqrt(Fabs((1-cos_f)*(1+cos_f)));
   int err_code=0;
   
   pvect[11]=PcmOut*cos_f;
   pvect[15]=-pvect[11];
   pvect[10]=PcmOut*sin_f;
   pvect[14]=-pvect[10];
      
   r = (*sqme22)(nsub22,sqrt(4*M_PI*alphaQCD(GGscale)),pvect,NULL,&err_code);
    
   err_code=0;
   return r * totcoef;
}

static double kinematic_22(double PcmIn, double cs, REAL*pmass, REAL*pvect)
{  int i;
   for(i=0;i<16;i++) pvect[i]=0;
   REAL sqrtS=Sqrt(pmass[0]*pmass[0]+PcmIn*PcmIn)+Sqrt(pmass[1]*pmass[1]+PcmIn*PcmIn);
   REAL PcmOut = decayPcm(sqrtS,pmass[2],pmass[3]);
//printf(" PcmOut =%E (%E %E %E) \n", PcmOut,sqrtS,pmass[2],pmass[3] );   
   totcoef =  PcmOut /(32*M_PI*PcmIn*sqrtS*sqrtS);
   pvect[3] = PcmIn;
   pvect[7] =-PcmIn;
   pvect[0] = Sqrt(PcmIn*PcmIn   + pmass[0]*pmass[0]);
   pvect[4] = Sqrt(PcmIn*PcmIn   + pmass[1]*pmass[1]);
   pvect[8] = Sqrt(PcmOut*PcmOut + pmass[2]*pmass[2]);
   pvect[12]= Sqrt(PcmOut*PcmOut + pmass[3]*pmass[3]);

   pvect[11]=PcmOut*cs;
   pvect[15]=-pvect[11];  
   pvect[10]=Sqrt((PcmOut+pvect[11])*(PcmOut-pvect[11]));
   pvect[14]=-pvect[10];
   return totcoef*3.8937966E8;                  
}


double cs22(numout * cc, int nsub, double P, double cos1, double cos2 , int * err) 
{
  int i,k;
  REAL pmass[4];

  passParameters(cc);
  
  *(cc->interface->gtwidth)=0;
  *(cc->interface->twidth)=0;
  *(cc->interface->gswidth)=0;
  
  for(i=0;i<4;i++) cc->interface->pinf(nsub,1+i,pmass+i,NULL);  
  if(err) *err=0;

  double (*sqme22_mem)(int nsub, double GG, REAL *pvect, REAL*cb_coeff, int * err_code)=sqme22; 
  int  nsub22_mem=nsub22;
   
  sqme22=cc->interface->sqme;
  nsub22=nsub; 
  
  double res=0;
  if(!kin22(P,pmass)) res= 3.8937966E8*simpson(dSigma_dCos,cos1,cos2,0.3*eps,err);
  else { if(err) *err=-1; else printf("cs22: incoming momentum P=%.2E is too small\n",P);}  
 
  sqme22=sqme22_mem;
  nsub22=nsub22_mem;
  
  return res;
  
}

/*===================  Collider production ==========*/

static double sMin,sMax,pcmOut;
static REAL pmass[5];
static int pc1_,pc2_;
static int ppFlag;
static int i3,i4,i5;
static double pTmin_,METmin_;

static numout * colliderProduction(char * name1,char *name2, int nf, int J)
{ 
  char libname[100], process[100], lName1[20], lName2[20];
  numout *cc;
  int i,first;

  
  if(name1==NULL && name2==NULL) return NULL;
  if(name1) pname2lib(name1,lName1);
  if(name2) pname2lib(name2,lName2);
  if(name1==NULL)
  {
    sprintf(libname,"PP_%s",lName2);
    sprintf(process,"proton,proton->%s",name2);
  } else if(name2==NULL)
  { 
    sprintf(libname,"PP_%s",lName1);
    sprintf(process,"proton,proton->%s",name1); 
  } else
  {               
     if(strcmp(lName1,lName2)>0)sprintf(libname,"PP_%s%s",lName1,lName2);
     else                       sprintf(libname,"PP_%s%s",lName2,lName1); 
     sprintf(process,"proton,proton->%s,%s",name1,name2);
  }
     
  sprintf(libname+strlen(libname),"_nf%d",nf);
 
  if(J)
  {  sprintf(process+strlen(process),",proton");
     strcat(libname,"_J");
  }
  
  sprintf(process+strlen(process),"{");
  
  for(i=0,first=1;i<nModelParticles;i++) 
  { int pdg=abs(ModelPrtcls[i].NPDG);
    if(pdg && (pdg==21 ||(nf>=2 &&  pdg<=nf)|| (nf==1 && pdg==2) )) 
    { 
       if(!first) strcat(process,","); else  first=0;
       sprintf(process+strlen(process),"%s",ModelPrtcls[i].name);
       if(strcmp(ModelPrtcls[i].name,ModelPrtcls[i].aname))
       sprintf(process+strlen(process),",%s",ModelPrtcls[i].aname);
    }                              
  }
  
  cc=getMEcode(0,ForceUG,process,NULL,"",libname);
  
  return cc;
}


static double  cos_integrand(double xcos)
{ int err=0;
  REAL xsin=Sqrt((1-xcos)*(1+xcos));
  double q;
  pvect[9]=pcmOut*xcos;
  pvect[10]=pcmOut*xsin;
  pvect[13]=-pvect[9];
  pvect[14]=-pvect[10];
  q=Q_ren>0? Q_ren: pvect[10];  
  return  sqme22(nsub22,sqrt(4*M_PI*parton_alpha(q)),pvect,NULL,&err);  
}


static double  s_integrand(double y)
{  double r,q;
   REAL pcmIn;
   int  err;
   double pp=1.5;
   
   REAL  s=sMin*Pow(1+ y*( Pow(sMax/sMin,1-pp) -1) ,1/(1-pp));
   double x0=s/sMax;
   
   pcmIn=decayPcm(Sqrt(s),pmass[0], pmass[1]);
   if(pcmIn==0) return 0;
   pvect[0]=Sqrt(pmass[0]*pmass[0]+pcmIn*pcmIn);
   pvect[1]=pcmIn; pvect[2]=0; pvect[3]=0;
   pvect[4]=Sqrt(pmass[1]*pmass[1]+pcmIn*pcmIn);
   pvect[5]=-pcmIn; pvect[6]=0; pvect[7]=0;
   pcmOut=decayPcm(Sqrt(s),pmass[2], pmass[3]);
   pvect[8]=Sqrt(pmass[2]*pmass[2]+pcmOut*pcmOut);
   pvect[11]=0;
   pvect[12]=Sqrt(pmass[3]*pmass[3]+pcmOut*pcmOut);
   pvect[15]=0;

   if(pcmOut<=pTmin_) return 0;
   double sn=pTmin_/pcmOut;
   double cs=Sqrt((1-sn)*(1+sn)); 
   r=  3.8937966E8*pcmOut/(32*M_PI*pcmIn*s)*simpson(cos_integrand,-cs,cs,1.E-3,&err);
   if(err)  printf("error in simpson cs22 line 210\n");
   
   q=Q_fact>0? Q_fact:Sqrt(s);
     r*=convStrFun2(x0,q,pc1_,pc2_,ppFlag);
   r*=pow(s/sMax,pp)*(1- pow(sMin/sMax,1-pp))/(1-pp);
   return r; 
}

#define pt2etRange 1.0

static double M45_min,M45_max,S34_min,S34_max,S35_min,S35_max;
static int npole34=0,npole35=0,npole45,npole12;
static double * pole34=NULL,*pole35=NULL,*pole45=NULL,*pole12=NULL;
static int MET=0;

#define METDIM 7
static int fillArr;
static double metGrig[METDIM]={250,300,350,400,450,500,550}; 
static double metArr[METDIM],dmetArr[METDIM],metArr_[METDIM],dmetArr_[METDIM];

//static double PTarr[3]={log(200),log(400),log(600)};
static double PTarr[3]={5.298317367,5.991464547,6.396929655};

//static double MDarr[5]={50,100,500,1000,3000};
static double MDarr[5]={3.912023005,4.605170186,6.214608098,6.907755279,8.006367568};
static double  map_pt2et[8][3][5][4]=

#include "data/et_tab.inc"


static double ParamInterpolation(int iPar, int ch,double PT, double MD)
{   double X[3],Y[3][5],res; 
    int i,j;

if(MD<40) MD=40;
if(MD> 3500) MD=3500;
if(PT<170) PT=170;
if(PT>800) PT=800;
   
    for(i=0;i<3;i++){ 
                      for(j=0;j<5;j++) Y[i][j]= map_pt2et[ch-1][i][j][iPar];
                      X[i]=polintN(log(MD), 5,MDarr,Y[i]);   
                    }
    res=  polintN(log(PT),3,PTarr,X);
    
    if(res< 0 && iPar==2)
    {  printf("ch+1=%d  iPar=%d res=%e  PT=%E MD=%E \n", ch+1,iPar, res,PT,MD);
       for(i=0;i<3;i++)
       {  printf( "%E ",X[i]);
          for(j=0;j<5; j++) printf(" %E",Y[i][j]);
           printf( "\n");
       }      
    }
    return res;            
}


static void getMET(double *r, double PT,double M12,double M45, int ch, double x, double w)
{ 
     static  pthread_mutex_t lockFillKey=PTHREAD_MUTEX_INITIALIZER;
     double C,delta,w0,p;
     int i;

       
     C=    ParamInterpolation(0, ch, PT, M45);
     delta=ParamInterpolation(1, ch, PT, M45); 
     w0=   ParamInterpolation(2, ch, PT, M45);  
     p=    ParamInterpolation(3, ch, PT, M45);

/*
C=0.5;
delta=0;
w0=40;
p=2;
*/     
//printf("C=%E delta=%e w0=%e p=%e\n", C,delta,w0,p);     
     
     double det=2*PT*pt2etRange*(x-0.3);
     double et= PT+delta +det;
     if(et<METmin_) {*r=0;return;}
     double ww=2*PT*pt2etRange*C*pow(w0,2*p-1)/(2.6*pow(w0*w0 +det*det,p));
if(!isfinite(ww)) {printf("ww=%E x=%E pt2etRange=%E  C=%e  pow(w0,2*p-1)=%E p=%e w0=%E pow(w0*w0 +det*det,p)=%e   \n",
                           ww,x,      pt2etRange,    C,    pow(w0,2*p-1 ), p, w0,pow(w0*w0 +det*det,p));}
if(!isfinite(*r)) {printf("*r=%E\n",*r); }                
     (*r)*= ww; 
     if(w && *r)
     {  if(nPROCSS>1)pthread_mutex_lock(&lockFillKey); 
        for(i=0;i< METDIM ;i++) if(et>metGrig[i]) metArr_[i]+=(*r)*w;
        if(nPROCSS>1)pthread_mutex_unlock(&lockFillKey);
     }   
}

static double veg2_intergrand(double *x, double w)
{
   double r;
   REAL M12,pcmIn;
   int err;
   REAL Pout,sn_,cs_;

   REAL pvect[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};        
   
   M12=Sqrt(sMin)+x[0]*(Sqrt(sMax)-Sqrt(sMin)); 
   pcmIn=decayPcm(M12,pmass[0], pmass[1]);
   if(pcmIn==0) return 0;
   
   if(M12<=pmass[2]+pmass[3]) return 0;
   Pout=decayPcm(M12,pmass[2], pmass[3]); 
   if(Pout<=pTmin_) return 0; 
   
   sn_=pTmin_/Pout;

   cs_=Sqrt(1-sn_*sn_);
   double cs=cs_*(2*x[1]-1),sn=Sqrt(1-cs*cs),PT=Pout*sn;
   
   r=  kinematic_22(pcmIn,cs,pmass, pvect);

   {   double q= Q_ren>0? Q_ren : PT;
       double x0=M12*M12/sMax;
      
       r*= sqme22(nsub22,Sqrt(4*M_PI*parton_alpha(q)),pvect,NULL,&err);

       q= Q_fact>0? Q_fact: PT; 
       r*= fabs(convStrFun2(x0,q,pc1_,pc2_,ppFlag));       

   }
   r*= 2*cs_*(Sqrt(sMax) - Sqrt(sMin))*2*M12/sMax;

   if(MET) 
   { 
   
      getMET(&r, PT,M12,pmass[i4],chan,x[2],w*fillArr);
   
   } 
   return r; 
}



static double veg3_intergrand(double *x, double w)
{
   double r;
//   double pp=1.5;
   REAL M12,pcmIn;
   int err=0;
   REAL P,M45,m3q,sn_,cs_,J45,cs45,S34,S35;
   REAL M34_min=Sqrt(S34_min),M34_max=Sqrt(S34_max),
          M35_min=Sqrt(S35_min),M35_max=Sqrt(S35_max);
   REAL M34,M35;
  
   REAL pvect[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};        
   
   M12=Sqrt(sMin)+x[0]*(Sqrt(sMax)-Sqrt(sMin)); 
   pcmIn=decayPcm(M12,pmass[0], pmass[1]);
   if(pcmIn==0) return 0;
 
   M45=M45_min+x[1]*(M45_max-M45_min);
   
   if(M12<=pmass[i3]+M45) return 0;
   P=decayPcm(M12,pmass[i3], M45); 
   if(P<=pTmin_) return 0; 
   
   sn_=pTmin_/P;

   cs_=Sqrt(1-sn_*sn_);
   double cs=cs_*(2*x[2]-1),sn=Sqrt(1-cs*cs),PT=P*sn;
   
   if(npole34==0 && npole35==0)  { cs45=(2*x[3]-1); J45=2;}  else 
   {  
      double E3= 0.5*(M12*M12-M45*M45-pmass[i3]*pmass[i3])/M45; 
      double p= Sqrt(E3*E3-pmass[i3]*pmass[i3]);
      double pcm2=decayPcm(M45,pmass[i4],pmass[i5]);
      double E4=Sqrt(pcm2*pcm2+pmass[i4]*pmass[i4]);
      double E5=Sqrt(pcm2*pcm2+pmass[i4]*pmass[i5]);

      if(npole35==0)
      { 
         M34=M34_min+x[3]*(M34_max-M34_min);
         cs45=(M34*M34-pmass[i4]*pmass[i4]-pmass[i3]*pmass[i3]-2*E3*E4)/(2*p*pcm2);
         J45=M34*(M34_max-M34_min)/(p*pcm2);
      } else if(npole34==0)
      { 
         M35=M35_min+x[3]*(M35_max-M35_min);
         cs45=-(M35*M35 -pmass[i5]*pmass[i5]-pmass[i3]*pmass[i3]-2*E3*E4)/(2*p*pcm2);
         J45=M35*(M34_max-M34_min)/(p*pcm2);
      } else 
      {  int i;
         double pr34=0,pr35=0;
         if(x[3]<0.5)
         {   M34=M34_min+2*x[3]*(M34_max-M34_min);
             S34=M34*M34;  
             cs45=(S34 -pmass[i4]*pmass[i4]-pmass[i3]*pmass[i3]-2*E3*E4)/(2*p*pcm2);
             S35=pmass[i5]*pmass[i5]+pmass[i3]*pmass[i3]+2*E3*E5 -2*p*pcm2*cs45;
             J45=2*M34*(M34_max-M34_min)/(p*pcm2);
         } else 
         {   M35=M35_min+(2*x[3]-1)*(M35_max-M35_min);
             S35=M35*M35;  
             cs45=-(S35 -pmass[i5]*pmass[i5]-pmass[i3]*pmass[i3]-2*E3*E5)/(2*p*pcm2);
             S34=pmass[i4]*pmass[i4]+pmass[i3]*pmass[i3]+2*E3*E4 +2*p*pcm2*cs45; 
             J45=2*M35*(M35_max-M35_min)/(p*pcm2);
         } 
         for(i=0;i<npole34;i++)
         { double m=pole34[2*i],w=pole34[2*i+1];
           pr34+=  M35*m*w/( (m*m-S34)*(m*m-S34) +m*m*w*w);
         }
         pr34+=1/(M34_max-M34_min); 
           
         for(i=0;i<npole35;i++)
         { double m=pole35[2*i],w=pole35[2*i+1];
           pr35+=  M35*m*w/( (m*m-S35)*(m*m-S35) +m*m*w*w);
         } 
         pr35+=1/(M35_max-M35_min); 
         if(x[3]<0.5) J45*=pr34/(pr34+pr35); else  J45*=pr35/(pr34+pr35);  
      }    
   } 
   if(fabs(cs45)>1) return 0;
   r=  kinematic_23(pcmIn,i3,M45, cs_*(2*x[2]-1) ,cs45,M_PI*x[4],pmass, pvect)*4*M_PI*(M45_max-M45_min)*J45*cs_/pcmIn;
   double qr,qf,q0;
   q0=0.5*(PT+sqrt(PT*PT+M45*M45));
   if(MET)
   { qr=q0;
     qf=q0;
   } else 
   { qr=Q_ren>0?  Q_ren: q0;   
     qf=Q_fact>0? Q_fact:q0;
   }
      
   r*= sqme22(nsub22,sqrt(4*M_PI*parton_alpha(qr)),pvect,NULL,&err);
   r*= convStrFun2(M12*M12/sMax  ,qf,pc1_,pc2_,ppFlag);       
   r*= 2*M12/sMax*(sqrt(sMax) - sqrt(sMin));

   if(MET) getMET(&r, PT,M12,M45,chan,x[5],w*fillArr);
   
   return r; 
}

   
static  void impGrid(int N, double*x, double *y)
{ int i;
  double s,rc;
  double alph=1.5;
  double X[100],r[100];

  for(s=0, i=0; i<N; i++) s+=y[i];
  for(rc=0,i=0; i<N; i++)
  {
     r[i] = 0;
     if( y[i]  > 0)
     {  double xoln = log(s/y[i]);
           if(xoln<0.01)         r[i]=1;
           else if (xoln <= 70)  r[i] = pow( (1 - exp(-xoln))/xoln, alph);
           else                  r[i] = pow(  1/xoln,               alph);
     }
     rc += r[i];  
  }
  rc /= N;
  if(rc)
  {  double dr=0,xo=0,xn=0;
     int k=0;
     for(i=1;i<N;) 
     {  for(;dr<rc;k++) dr+=r[k];
        xo=x[k-1];
        xn=x[k];
        for(;dr>=rc; i++) {dr-=rc; X[i] = xn-(xn-xo)*(x[N]-x[0])*dr/r[k-1];}
     }
     for(i=1;i<N;i++) x[i]=X[i];
  }
}


static void setGrid(int N, double*x, double x0, double x1,double M0,double M1,int Npole,double *poles)
{ 
  int i,k,nStep=10; 
  double *y=malloc(N*sizeof(double));
  for(i=0;i<=N;i++) x[i]=x0+i*(x1-x0)/N;
  for(k=0;k<nStep;k++)
  {  int l;
     for(i=0;i<N;i++)
     {  double m1=M0+(x[i]  -x0)*(M1-M0)/(x1-x0);
        double m2=M0+(x[i+1]-x0)*(M1-M0)/(x1-x0); 
        y[i]= y[i]= (m2-m1)/(M1-M0); // M0*M1/(M1-M0)*(1/m1-1/m2); ?????     
        for(l=0;l<Npole;l++) 
        {
           double m=poles[2*l],w=poles[2*l+1]; 
           double x1=(m1*m1/m -m)/w, x2=(m2*m2/m -m)/w;
           y[i]+=atan(x2)-atan(x1);
        }       
     }
     impGrid(N, x, y);           
   }
   free(y);
}

static void getPoles(numout*cc, int nsub, char * s0,double mMin, double mMax, int*npole,double**poles)
{ int i,n,m,w;
  char*s;
  double mass;
  CalcHEP_interface * CI=cc->interface;
  
//  printf("Get Pole (%d %d %d)\n",s0[0],s0[1],s0[2]);
  
  *npole=0;
  if(s0[0]==0 || s0[1]==0) return;
      
  for(i=0; s0[i+1]!=0;)
  {  if(s0[i]>s0[i+1]) 
     { int a=s0[i]; s0[i]=s0[i+1];s0[i+1]=a; 
       if(i) i--;else i++;
     } else i++;
  }
  
  for(n=1;(s=CI->den_info(nsub,n,&m,&w,NULL));n++)
  { double mm=0,ww=0;
    if(m) mm=fabs(CI->va[m]);
    if(w) ww=fabs(CI->va[w]);
            
    if(strcmp(s0,s)==0 && mm>mMin && mm < mMax)
    { 
//printf(" %s %E %E  x= %E \n",CI->varName[m],CI->va[m], CI->va[w],(CI->va[m]-mMin)/(mMax-mMin) );

     (*npole)++;

//printf("npole=%d size=%d\n", *npole, sizeof(double)*2*(*npole));     
     *poles=realloc(*poles,sizeof(double)*2*(*npole));
     (*poles)[2*((*npole)-1)]= mm; 
     (*poles)[2*(*npole-1)+1]= ww;  
    }
  }
//    else  printf("-- %s %E %E  x= %E \n",CI->varName[m],CI->va[m], CI->va[w],(CI->va[m]-mMin)/(mMax-mMin) );

}


static double vegas_cycle(vegasGrid *vegPtr, double eps, double aeps,int maxStep, int*NN, double fact,double *dI)
{ int k,l;
  double *ti=malloc(maxStep*sizeof(double));
  double *dti=malloc(maxStep*sizeof(double));
    double ii,dii,chi2;
  
  for(k=0;k<maxStep;k++)
  { 
    double s0=0,s1=0,s2=0;    

    vegas_int(vegPtr, *NN , 1.5, nPROCSS  , ti+k, dti+k);
    
//    printf("ti=%E dti=%E  NN=%d \n",ti[k], dti[k],NN);
    if(dti[k]==0){ dii=0; ii=ti[k];  break;}
    for(l=k;l>=k/2;l--)
    { s0+=1/(dti[l]*dti[l]);
      s1+=ti[l]/(dti[l]*dti[l]);
      s2+=ti[l]*ti[l]/(dti[l]*dti[l]);
      if(l!=k)
      { 
        ii=s1/s0;
        dii=1/sqrt(s0);
        chi2=(s2-s1*s1/s0)/(k-l+1);
        if(chi2> 1 )dii*=sqrt(chi2);

        if(dii<eps*fabs(ii)) break;
        if(dii<aeps) break;
      }  
    }
    if(k && (dii<eps*fabs(ii) || dii<aeps )) break;
    (*NN)*=fact;    
  }  
  free(ti); free(dti);
  *dI=dii;
  return ii;
}


static double hColliderStat(double Pcm, int pp, int nf, double Qren,double Qfact, char * name1,char *name2,double pTmin,int met, int wrt)
{ 
  double  sigma_tot=0;
  int i;
  numout *cc;
  int n1=0,n2=0;
  int nout;
  double dI,m1=0,m2=0;
  
  if(met) MET=1; else MET=0;
  if(met) for(i=0;i<METDIM;i++) {metArr[i]=0; dmetArr[i]=0;}
  ppFlag=pp;   
  Q_fact=Qfact;
  Q_ren=Qren;
    
  if(nf>5)nf=5; 
  if(nf<0)nf=0;
    
  if(name1) 
  { n1=pTabPos(name1);  
    if(n1==0) { printf("%s - no such particle\n",name1); return 0;}
    m1=pMass(name1);
  }

  if(name2) 
  {  n2=pTabPos(name2);  
     if(n2==0) { printf("%s - no such particle\n",name2); return 0;}
     m2=pMass(name2);   
  }

  if(!(n1 || n2)) return 0;
 
  
  sMax=4*Pcm*Pcm; 
  if(pTmin<0) pTmin_=0; else  pTmin_=pTmin;
  sMin=sqrt((m1+m2)*(m1+m2)+pTmin_*pTmin_)+pTmin_;  
  sMin*=sMin;
  
  cc=colliderProduction( name1,name2, nf, pTmin>0);
  if(!cc) return 0; 
  if(passParameters(cc)) return 0;
      
  sigma_tot=0;
  
  nout=cc->interface->nout;
  sqme22=cc->interface->sqme; 
  
  for(nsub22=1;nsub22<=cc->interface->nprc; nsub22++) 
  { 
    char*n[5];
    char buff[40];
    double tmp=0,dI;
    for(i=0;i<2+nout;i++) n[i]=cc->interface->pinf(nsub22,i+1,pmass+i,PC+i); 
    
    
    if(PC[0]<=PC[1])
    { pc1_=PC[0];
      pc2_=PC[1];

      chan=0;
      if(MET)
      {      
         switch(pc1_)
         { case -3: if(pc2_== 3) chan=8; else if(pc2_==21) chan=5; break;
           case -2: if(pc2_== 2) chan=7; else if(pc2_==21) chan=4; break;
           case -1: if(pc2_== 1) chan=6; else if(pc2_==21) chan=2; break;
           case  1: if(pc2_==-1) chan=6; else if(pc2_==21) chan=2; break;
           case  2: if(pc2_== 2) chan=7; else if(pc2_==21) chan=3; break;
           case  3: if(pc2_==-3) chan=8; else if(pc2_==21) chan=5; break;
         }      
         if(chan==0) continue;
      } 

      if(wrt) { buff[0]=0; for(i=0;i<2+nout;i++) {sprintf(buff+strlen(buff),"%s ",n[i]); if(i==1) sprintf(buff+strlen(buff)," -> ");}  printf("%-30.30s",buff); }
      if(nout>1)
      {
         for(i3=2;i3<5;i3++) if( (PC[i3]<6 && PC[i3]>-6) || PC[i3]==21) break;
         for(i4=2;i4<5;i4++) if(i4!=i3) break;
         for(i5=2;i5<5;i5++) if(i5!=i3 && i5!=i4) break; 
      }

       
      switch(nout)
      { case 1: 
        { 
          double pcmIn=decayPcm(pmass[2],pmass[0], pmass[1]);
          double q;
          int err=0;
          
          if(pcmIn==0) return 0;
           pvect[0]=Sqrt(pmass[0]*pmass[0]+pcmIn*pcmIn);
           pvect[1]=pcmIn; pvect[2]=0; pvect[3]=0;
           pvect[4]=Sqrt(pmass[1]*pmass[1]+pcmIn*pcmIn);
           pvect[5]=-pcmIn; pvect[6]=0; pvect[7]=0;
           pvect[8]=pmass[2];pvect[9]=pvect[10]=pvect[11]=0;
           q=Q_fact>0? Q_fact: pmass[2];
           tmp=convStrFun2(pmass[2]*pmass[2]/sMax,q,pc1_,pc2_,ppFlag);
           q=Q_ren>0? Q_ren: pmass[2];
           tmp*=cc->interface->sqme(nsub22,sqrt(4*M_PI*parton_alpha(q)),pvect,NULL,&err);
           tmp*=389379660.0*M_PI/(2*pcmIn*pmass[2]*sMax);
           if(wrt)printf("cs=%E \n",tmp); 
          break;
        }
        case 2: if(met)
                { vegasGrid *vegPtr=vegas_init(3,veg2_intergrand,50);
                  convStrFun2(0.1,100,pc1_,pc2_,ppFlag); //testing call
                  double eps=0.001,aEps=1E-6;
                  int NN=10000;
                  if(fabs(sigma_tot)*eps>aEps) aEps=sigma_tot*eps;
                  tmp=vegas_cycle(vegPtr,eps, aEps, 100,&NN,1.1,&dI);
                  for(i=0;i<METDIM;i++) {metArr_[i]=0; dmetArr_[i]=0;}
                  fillArr=1;  
//                  printf("tmp=%E dI=%e\n", tmp,dI);
                  vegas_int(vegPtr, 2*NN , 1.5, nPROCSS  , &tmp, &dI);
//                  printf("      tmp=%E dI=%e\n", tmp,dI); 
                  fillArr=0;
                  for(i=0;i<METDIM;i++) {metArr_[i]/=2*vegPtr->intCubes; metArr[i]+=metArr_[i]; }
                  for(i=0;i<METDIM;i++)  printf(" %.2E ",metArr_[i]);
                  printf("( %.2f%%)\n",100*dI/metArr_[0]);
                  vegas_finish(vegPtr);
                }
                else{int err;  tmp=simpson(s_integrand,0.,1.,1.E-2,&err); if(err) printf("error in simpson cs22 line 696\n"); }
                 
                if(!MET && wrt)printf("cs=%E \n", tmp);  
                break;
        case 3:
        {  double m3q;
           double vdim=5;
           if(met)vdim++;
           vegasGrid *vegPtr=vegas_init(vdim,veg3_intergrand,50);
           char s0[3]={0,0,0};
           double eps,aEps;
         
           M45_min=pmass[i4]+pmass[i5];
           S35_min=pmass[i3]+pmass[i5]; S35_min*=S35_min;
           S34_min=pmass[i3]+pmass[i4]; S34_min*=S34_min;
         
           m3q=pmass[i3]*pmass[i3];
           M45_max= sMax-2*sqrt((pTmin_*pTmin_+m3q)*sMax)+m3q;
           if(M45_max<=M45_min*M45_min) return 0;
           M45_max=sqrt(M45_max);
           S34_max=sqrt(sMax-pmass[i5]); S34_max*=S34_max;
           S35_max=sqrt(sMax-pmass[i4]); S35_max*=S35_max;
         
           s0[0]=i3+1;s0[1]=i4+1; getPoles(cc,nsub22,s0,sqrt(S34_min),sqrt(S34_max),&npole34,&pole34);
           s0[0]=i3+1;s0[1]=i5+1; getPoles(cc,nsub22,s0,sqrt(S35_min),sqrt(S35_max),&npole35,&pole35);        
           s0[0]=i4+1;s0[1]=i5+1; getPoles(cc,nsub22,s0,M45_min,      M45_max,      &npole45,&pole45); 
           s0[0]=1;   s0[1]=2;    getPoles(cc,nsub22,s0,sqrt(sMin),   sqrt(sMax),   &npole12,&pole12);                 

           setGrid(50,  vegPtr->x_grid[0]   , 0, 1,sqrt(sMin),sqrt(sMax),npole12,pole12);
           setGrid(50,  vegPtr->x_grid[1]   , 0, 1,M45_min,M45_max,npole45,pole45);
         
           if(npole34 && npole35)
           { setGrid(25,  vegPtr->x_grid[3], 0, 0.5, sqrt(S34_min),sqrt(S34_max),npole34,pole34);
             setGrid(25,  vegPtr->x_grid[3]+25, 0.5, 1, sqrt(S35_min),sqrt(S35_max),npole35,pole35);
           }else if(npole34)
             setGrid(50,  vegPtr->x_grid[3], 0, 1, sqrt(S34_min),sqrt(S34_max),npole34,pole34);
           else if (npole35)
             setGrid(50,  vegPtr->x_grid[3], 0, 1, sqrt(S35_min),sqrt(S35_max),npole35,pole35);

           convStrFun2(0.1,100,pc1_,pc2_,ppFlag); //testing call
           eps=0.001;
           aEps=1E-6;
           int NN=10000;
           if(fabs(sigma_tot)*eps>aEps) aEps=sigma_tot*eps;
           tmp=vegas_cycle(vegPtr,eps, aEps, 200,&NN, 1.1,&dI); 
           if(met)
           {       for(i=0;i<METDIM;i++) {metArr_[i]=0; dmetArr_[i]=0;}
                  fillArr=1;  
//                  printf("tmp=%E dI=%e\n", tmp,dI);
                  vegas_int(vegPtr, 2*NN , 1.5, nPROCSS  , &tmp, &dI);
//                  printf("      tmp=%E dI=%e\n", tmp,dI); 
                  fillArr=0;
                  for(i=0;i<METDIM;i++) {metArr_[i]/=2*vegPtr->intCubes; metArr[i]+=metArr_[i]; }
                  if(wrt) 
                  {  for(i=0;i<METDIM;i++)  printf(" %.2E ",metArr_[i]);
                     printf("( %.2f%%)\n",100*dI/metArr_[0]);
                  }   
           }
           else if(wrt)printf("cs=%E +/-%E \n", tmp,dI);
           vegas_finish(vegPtr);
           free(pole12); free(pole34); free(pole35); free(pole45);
           pole12=pole34=pole35=pole45=NULL;
           break;
        }  
     }
     sigma_tot+=tmp;
    }
  }  

  return sigma_tot;
}

double hCollider(double Pcm,int pp,int nf,double Qren,double Qfact,char*name1,char*name2,double pTmin, int wrt)
{ 
  return hColliderStat(Pcm, pp, nf, Qren, Qfact, name1, name2, pTmin, 0, wrt);
} 

#define DELPHES_FACTOR (1.15)  
#define NF 3  


double monoJet(void)
{ 
  int i,i0;
  double ret;
  double METmin=250;
  double bg[METDIM]= {51800,19600,8190,3930,2050,1040,509};
  double dBg[METDIM]={ 2000,  830, 400, 230, 150, 100, 66};
  double exp[METDIM]={52200,19800,8320,3830,1830, 934,519};
  double sum[METDIM];
  double CL=0,CL0=0;
  for(i=0;i<METDIM;i++) {metArr[i]=0; dmetArr[i]=0;}
  char oldPDF[100]={""};
  if(strcmp(pdfName,"NNPDF23_lo_as_0130_qed"))
  {  strcpy(oldPDF,pdfName);
     setPDT("NNPDF23_lo_as_0130_qed");
  } else strcpy(oldPDF,"NNPDF23_lo_as_0130_qed");
  

  METmin_=METmin;
  for(i=0;i<METDIM;i++) sum[i]=0;
  
  if(CDM1)
  {   
     ret=hColliderStat(4000, 1, NF, 0, 0, CDM1,aCDM1,METmin/(1+pt2etRange), 1, 0);
     for(i=0;i<METDIM;i++) sum[i]+=metArr[i];
     if(strcmp(CDM1,aCDM1))
     { ret=hColliderStat(4000, 1, NF, 0, 0, CDM1,CDM1,METmin/(1+pt2etRange), 1, 0); 
       for(i=0;i<METDIM;i++) sum[i]+=2*metArr[i];
     }
  }
  if(CDM2)
  {   
     ret=hColliderStat(4000, 1, NF, 0, 0, CDM2,aCDM2,METmin/(1+pt2etRange), 1, 0);
     for(i=0;i<METDIM;i++) sum[i]+=metArr[i];
     if(strcmp(CDM2,aCDM2))
     { ret=hColliderStat(4000, 1, NF, 0, 0, CDM2,CDM2,METmin/(1+pt2etRange), 1, 0); 
       for(i=0;i<METDIM;i++) sum[i]+=2*metArr[i];
     }
  }
        
  for(i=0;i<METDIM;i++) sum[i]*=DELPHES_FACTOR;
  printf("\n Analysis: CMS monojet 8 TeV arXiv:1408.3583\n");
  printf("MET [GeV]       >250      >300      >350      >400      >450      >500      >550\n"); 
//  printf("%30.30s","csSum[pb]");  for(i=0;i<METDIM;i++) printf(" %.2E ",metArr[i]);
  printf("%s","Signal events ");   for(i=0;i<METDIM;i++) printf(" %.2E ",sum[i]*19.7*1000);
//  if( abs(pNum(name1)==12)) { printf("\n%30.30s","signal*3");   for(i=0;i<METDIM;i++) printf(" %.2E ",3*metArr[i]*19.7*1000);} 
  printf("\n1-CLs expected");
  for(i=0;i<METDIM;i++)
  { double s=sum[i]*19.7E3, b=bg[i],db=dBg[i];
    double  CLs0=(1-erf((s)/sqrt(2*(s+db*db))));;
    printf(" %.2E ",1-CLs0);
    if(1-CLs0>CL0) { CL0=1-CLs0; i0=i;}
  }
  printf("\n The region most  sensitive to signal is MET>%.0f",250.+i0*50.);
  printf("\n1-CLs         ");
  CL0=0;
  for(i=0;i<METDIM;i++)
  { double s=sum[i]*19.7E3, b=bg[i],db=dBg[i], n=exp[i];
    double CLsb=0.5*(1-erf((s+b-n)/sqrt(2*(s+db*db))));
    double CLb=0.5*(1-erf((b-n)/sqrt(2*(db*db))));
    double CLs=CLsb/CLb;  
    printf(" %.2E ",1-CLs);
       if(1-CLs>CL0) { CL0=1-CLs; i0=i;}
     if(i==i0) CL=1-CLs;
  }
  printf("\n");                              
  if(strcmp(pdfName,oldPDF)) restorePDF(oldPDF);
  return CL0;
}


#ifdef plazmaWidth
static numout* plazmaWidth_cc;
static double plazmaWidth_T;
static double plazmaWidth_m[4];
static double plazmaWidth_integrand(double Pcm)
{ int err;
  double E1,E2,sqrt_s; 
  if(Pcm==0) return 0;  
  E1=sqrt(Pcm*Pcm+plazmaWidth_m[0]*plazmaWidth_m[0]);
  E2=sqrt(Pcm*Pcm+plazmaWidth_m[1]*plazmaWidth_m[1]);
  sqrt_s=E1+E2;
  if(sqrt_s<=plazmaWidth_m[2]+plazmaWidth_m[3]) return 0;
  
  return 4*bessk1(sqrt_s/plazmaWidth_T)*cs22(plazmaWidth_cc,1,Pcm, -1., 1. , &err)*pow(sqrt_s*Pcm,3.)/E1/E2; 
}

double plazmaWidth(char *process,double T)
{  char libName[40];
   plazmaWidth_T=T;
   process2Lib(process,libName);
   process2Mass(process,plazmaWidth_m);
   plazmaWidth_cc=getMEcode(0,ForceUG,process,NULL,NULL,libName);
   { int err;
     double r=simpson(plazmaWidth_integrand,0., 5*T,1.E-3,NULL)/(4*M_PI*M_PI*plazmaWidth_m[1]*plazmaWidth_m[1]*bessk2(plazmaWidth_m[1]/T))
        /3.8937966E8;
     if(err) printf("error in simpson cs22 line 873\n");
     return r;
   }      
}
#endif
/*============ Fortran ==========*/

double cs22_(int*ccf,int*nsub,double*P,double*cos1,double*cos2,int*err)
{ numout*cc;
  memcpy(&cc,ccf,sizeof(cc));
  return cs22(cc,*nsub,*P,*cos1,*cos2 ,err);
} 

void  sethelicities_(double *h1,double *h2) { Helicity[0]=*h1; Helicity[1]=*h2;}


double hcollider_(double*Pcm, int*pp, int* nf, double*Qren,double*Qfact, char * name1,char *name2,double*pTmin,int*wrt,int len1,int len2)
{ 
  char cname1[20], cname2[20];
  char *cname1_,*cname2_;
  fName2c(name1,cname1,len1);
  fName2c(name2,cname2,len2);
  if(strlen(cname1)==0) cname1_=NULL; else  cname1_=cname1;
  if(strlen(cname2)==0) cname2_=NULL; else  cname2_=cname2;
  return  hCollider(*Pcm, *pp, *nf, *Qren,*Qfact, cname1_,cname2_,*pTmin,*wrt);  
}

double monojet_(void) { return monoJet();}
