#include <sys/utsname.h>
#include <unistd.h>
#include "micromegas.h"
#include "micromegas_aux.h"
#include "micromegas_f.h"
#include "../CalcHEP_src/include/rootDir.h" 

#define ZeroCS 1E-40


char* CDM1=NULL, *CDM2=NULL,*aCDM1=NULL,*aCDM2=NULL;
aChannel* omegaCh=NULL;
aChannel* vSigmaTCh=NULL;
REAL *Qaddress=NULL;
int maxPlistLen=0; 


static int do_err=0;

static int NT;
static double *XX,*YY;
static void loadINTER(int N, double *x, double *y)
{ NT=N;XX=x;YY=y;}
double INTER(double x) { return polint3(x,NT,XX,YY);} 

static REAL Zero=0;

typedef  struct{ int dim; double *x; double *y;}   polintArgStr;
double polint_arg(double x, void* arg) { polintArgStr*arg_=arg; return polint3(x,arg_->dim,arg_->x,arg_->y);} 


typedef  struct{ int virt,i3;double br,w[2];numout*cc23; int nTab; double pcm23max; double *pcmTab; double *csTab;} processAuxRec;

typedef  processAuxRec* processAux;
static   processAux* code22Aux0, *code22Aux1,*code22Aux2;

static   processAux AUX;

double sWidth=0.01;

extern int  WIDTH_FOR_OMEGA;

static int neg_cs_flag;

static int NC=0;

static char ** inP;
static int  *  inAP;
static int  *  inG;
static int  *  inNum;
static int  *  feebleParticles=NULL;
static numout ** code22_0;
static numout ** code22_1;
static numout ** code22_2;

static int *inC0;    /* combinatoric coefficients  NCxNC*/
static int *inC1; 
static int *inC2;
static numout*cc23=NULL;

static REAL **inMassAddress;
static double *inMass;  /* masses */
static int *sort;
static int gaussInt=1;


static int LSP;

static double M1=0,M2=0;

static REAL pmass[6];
static int pdg[4];

#define XSTEP 1.1
static double eps=0.001; /* precision of integration */

static double MassCut=MPlanck;

static double s3f_;   /* to pass the Xf argument   */

static double T_;


static int Z4ch( char *name)
{  if(name[0]!='~') return 0;
   if(name[1]!='~') return 1;
   return 2;
}



#define IMPROVE


static double sigma23(double PcmIn)
{  int l,l_;
   double r;
   double brV1,MV1, MV2,wV1,wV2;
   int err;
   brV1=AUX[nsub22].br; 
   r=cs23(cc23,1,PcmIn,AUX[nsub22].i3,&err)/brV1/3.8937966E8;
   if(err) do_err=do_err|err;
   l=AUX[nsub22].virt;
   l_=5-l;  
   if(AUX[nsub22].w[l_-2])
   {  double m1,m2,sqrtS;
      m1=pmass[0];
      m2=pmass[1];
      MV1=pmass[l];
      MV2=pmass[l_];
      wV1=AUX[nsub22].w[l-2];
      wV2=AUX[nsub22].w[l_-2];
      sqrtS=sqrt(PcmIn*PcmIn+m1*m1) + sqrt(PcmIn*PcmIn+m2*m2); 
      if(wV1*wV2>0) r*=decayPcmW(sqrtS,MV1,MV2,wV1,wV2,0)/decayPcmW(sqrtS,MV1,MV2,wV1,0,0);    
   }   

   if(r<1.E-200) r=1.E-200;
   return log(r*PcmIn);
}   


static  double sigma(double PcmIn)
{ double r; 

//  if(AUX[nsub22].nTab>0 && PcmIn<=AUX[nsub22].pcmTab[AUX[nsub22].nTab-1])
  if(AUX[nsub22].cc23  && PcmIn<=AUX[nsub22].pcm23max)
  {  
     if(AUX[nsub22].nTab==0) r=exp(sigma23(PcmIn))/PcmIn;else
   //  {
   //  if(PcmIn<AUX[nsub22].pcmTab[0]) r= 0; else 
   //  { 
         r=exp(polint3(PcmIn,AUX[nsub22].nTab,AUX[nsub22].pcmTab,AUX[nsub22].csTab))/PcmIn;
        
   //  }
   //  }  
  } 
  else
  {  if(kin22(PcmIn,pmass)) return 0.; 
     if(gaussInt) r=gauss(dSigma_dCos,-1.,1.,5); else 
     { int err;
       r=simpson(dSigma_dCos,-1.,1.,0.3*eps,&err);
       if(err) {do_err= do_err|err;  printf("error in simpson omega.c line 130\n");}
     }   
     if(r<0) { neg_cs_flag=1;r=0;}

     if((VZdecay||VWdecay) && (AUX[nsub22].w[0] || AUX[nsub22].w[1] )) 
     { double f; 
       double sqrtS=sqrt(PcmIn*PcmIn+pmass[0]*pmass[0]) + sqrt(PcmIn*PcmIn+pmass[1]*pmass[1]); 
       if( sqrtS-pmass[2]-pmass[3] < 15*(AUX[nsub22].w[0]+AUX[nsub22].w[1]))   
         f=decayPcmW(sqrtS,pmass[2],pmass[3],AUX[nsub22].w[0],AUX[nsub22].w[1],5)/decayPcm(sqrtS,pmass[2],pmass[3]);     
       else   
       {   f=1; 
           if(AUX[nsub22].w[0]) f-= AUX[nsub22].w[0]/pmass[2]/M_PI;
           if(AUX[nsub22].w[1]) f-= AUX[nsub22].w[1]/pmass[3]/M_PI;
       }
       r*=f;
     }
  }  
#ifdef IMPROVE
  improveCrossSection(pdg[0],pdg[1],pdg[2],pdg[3],PcmIn,&r);
#endif  
  return r;
}  



static double geffDM(double T)
{ double sum=0; int l;

  for(l=0;l<NC;l++) if(!feebleParticles[oddPpos[sort[l]]]) 
  { int k=sort[l];
    double A=Mcdm/T*(inMass[k]-Mcdm)/Mcdm;
    if(A>15 || Mcdm +inMass[k] > MassCut) return sum;
    sum+=inG[k]*pow(inMass[k]/Mcdm,1.5)*exp(-A)*K2pol(1/(Mcdm/T+A));
  }
  return sum;
}

static double y_pass;


static double weight_integrand(double s3)
{  double x,gf;
   double T,heff,geff;

   if(s3==0.) return 0;
   T=T_s3(s3);
   heff=hEff(T);
   geff=gEff(T);
   gf=geffDM(T);
   return K1pol(T/(Mcdm*y_pass))*exp((1/T-1/T_)*(2-y_pass)*Mcdm)*sqrt(Mcdm/T)*heff/sqrt(geff)/(gf*gf*s3);
}


static double weightBuff_x[1000];
static double weightBuff_y[1000];
static int inBuff=0;


static double weight(double y)
{ int i;
  double w;
  for(i=0;i<inBuff;i++) if(y==weightBuff_x[i]) { return weightBuff_y[i]; }
  y_pass=y;
  { int err;
    w=  simpson(weight_integrand,0.,s3f_,0.3*eps,&err);
    if(err) {do_err=do_err|err; printf("error in simpson omega.c line 198\n");}
  }  
  if(inBuff<1000){weightBuff_x[inBuff]=y; weightBuff_y[inBuff++]=w;}
  return w;
}

static int exi;

static double s_integrand( double u)
{  double y,sv_tot,w;
   double Xf_1;
   double ms,md,sqrtS,PcmIn,res0;
   
   if(u==0. || u==1.) return 0.;
   
   long double u_=u,z=u_*(2-u_);
      
//   long double u_=u,z=(1-u_)*(1+u_);
   sqrtS=M1+M2-3*T_*logl(z);
   y=sqrtS/Mcdm;
   ms = M1 + M2;  if(ms>=sqrtS)  return 0;
   md = M1 - M2;
   PcmIn = sqrt((sqrtS-ms)*(sqrtS+ms)*(sqrtS-md)*(sqrtS+md))/(2*sqrtS);
   sv_tot=sigma(PcmIn);         
   res0=sqrt(2*y/M_PI)*y*(PcmIn*PcmIn/(Mcdm*Mcdm))*sv_tot*6*(1-u)*z*z; 
   if(exi) { return res0*weight(sqrtS/Mcdm); } else return  res0*K1pol(T_/sqrtS)*sqrt(Mcdm/T_);
}


static  double m2v(double m) { long double  e=expl(((M1+M2 -m)/T_)/3); if(e>=1) return e; else return e/(1+sqrtl(1-e));} 


typedef struct vGridStr
{  int n;
   double v[100];
   int pow[100];  
}  vGridStr;

static double v_max=1,v_min=0;

static vGridStr   makeVGrid(double mp,double wp)
{
  vGridStr grd;

  int pow_[6]={7,  3,  4, 4, 3,  5};
  double c[5]={ -10,-3, 0, 3, 10};

  int n,j,jmax=4;
  
  grd.v[0]=v_min;
  for(j=jmax,n=1 ;j>=0;j--)
  { double v=m2v(mp+c[j]*wp);
    if(isfinite(v) && v>v_min && v < v_max) 
    {  
      grd.v[n]=v;
      grd.pow[n-1]=pow_[j+1];
      grd.pow[n  ]=pow_[j];
      n++;
    }  
  }
  grd.v[n]=v_max;
  if(n==1) grd.pow[0]=5;
  grd.n=n;
  return grd;
}

static vGridStr   makeVGrid2(double mp,double wp)
{
  vGridStr grd;

  int pow_[6]={2, 4, 2};
  double c[5]={-3, 3};

  int n,j,jmax=1;
  
  grd.v[0]=v_min;
  for(j=jmax,n=1 ;j>=0;j--)
  { double v=m2v(mp+c[j]*wp);
    if(isfinite(v) && v>v_min && v < v_max) 
    { 
      grd.v[n]=v;
      grd.pow[n-1]=pow_[j+1];
      grd.pow[n  ]=pow_[j];
      n++;
    }  
  }
  grd.v[n]=v_max;
  if(n==1) grd.pow[0]=5;
  grd.n=n;
  return grd;
}



static vGridStr  crossVGrids(vGridStr * grid1, vGridStr * grid2)
{ vGridStr grid;
  int n=0,i1=1,i2=1,i;
  grid.v[0]=v_min;
  n=1;
  while(i1<=grid1->n && i2<=grid2->n)
  { double d1= grid1->pow[i1-1]/(grid1->v[i1]-grid1->v[i1-1]);
    double d2= grid2->pow[i2-1]/(grid2->v[i2]-grid2->v[i2-1]);
    double d = ( d1>d2? d1:d2);
    int m=(grid1->pow[i1-1] > grid2->pow[i2-1]? grid1->pow[i1-1]:grid2->pow[i2-1]);

    if(grid1->v[i1] < grid2->v[i2]) { grid.v[n]=grid1->v[i1++];}
    else                            { grid.v[n]=grid2->v[i2];
                                      if(grid1->v[i1]==grid2->v[i2])i1++; 
                                      i2++;
                                    }  
                                        
    grid.pow[n-1]=0.999+d*(grid.v[n]-grid.v[n-1]);

    if(grid.pow[n-1]>m)   grid.pow[n-1]=m;
    if(grid.pow[n-1]<2)   grid.pow[n-1]=2;

    n++;
  }
  grid.n=n-1;
  for(i=0;i<grid.n;i++) if(grid.v[i+1]-grid.v[i]>0.4 && grid.pow[i]<4)  grid.pow[i]=4;
  return grid;
}

static void printVGrid(vGridStr gr)
{
   printf("gr.n=%d\n", gr.n);
   int i;
   for(i=0;i<gr.n;i++) printf("      %d      ", gr.pow[i]);
   printf("\n");
   for(i=0;i<=gr.n;i++) printf(" %e", gr.v[i]);
   printf("\n");
}




static int testSubprocesses(void)
{
  static int first=1;
  int err,k1,k2,i,j;
  CDM1=CDM2=NULL;
  Mcdm=Mcdm1=Mcdm2=0;
  
  err=calcMainFunc();
  if(err>0) return err;
   
  if(first)
  {
    first=0;
    if(createTableOddPrtcls())
    { printf("The model contains uncoupled odd patricles\n"); exit(10);}
    
    for(i=0,NC=0;i<Nodd;i++,NC++) 
        if(strcmp(OddPrtcls[i].name,OddPrtcls[i].aname))NC++;
    oddPpos=(int*)malloc(NC*sizeof(int));        
    inP=(char**)malloc(NC*sizeof(char*));
    inAP=(int*)malloc(NC*sizeof(int));
    inG=(int*)malloc(NC*sizeof(int));
    inMassAddress=(REAL**)malloc(NC*sizeof(REAL*));
    inMass=(double*)malloc(NC*sizeof(double));
    inNum= (int*)malloc(NC*sizeof(int));
    if(!feebleParticles) 
    {  feebleParticles=(int*)malloc(sizeof(int)*nModelParticles);
       for(i=0;i<nModelParticles;i++) feebleParticles[i]=0;
    }
    sort=(int*)malloc(NC*sizeof(int));

    code22_0 = (numout**)malloc(NC*NC*sizeof(numout*));
    code22_1 = (numout**)malloc(NC*NC*sizeof(numout*));
    code22_2 = (numout**)malloc(NC*NC*sizeof(numout*));
            
    inC0=(int*)malloc(NC*NC*sizeof(int)); 
    inC1=(int*)malloc(NC*NC*sizeof(int));
    inC2=(int*)malloc(NC*NC*sizeof(int));
            
    code22Aux0=(processAux*) malloc(NC*NC*sizeof(processAux));
    code22Aux1=(processAux*) malloc(NC*NC*sizeof(processAux));
    code22Aux2=(processAux*) malloc(NC*NC*sizeof(processAux));
    
    for(i=0,j=0;i<Nodd;i++)
    {  
       inP[j]=OddPrtcls[i].name;
       for(int k=0;k<nModelParticles;k++) if(strcmp(inP[j],ModelPrtcls[k].name)==0 || strcmp(inP[j],ModelPrtcls[k].aname)==0) 
       { oddPpos[j]=k; break;}
       inNum[j]=OddPrtcls[i].NPDG;
       inG[j]=(OddPrtcls[i].spin2+1)*OddPrtcls[i].cdim;
       if(strcmp(OddPrtcls[i].name,OddPrtcls[i].aname))
       {
         inAP[j]=j+1;
         j++;
         oddPpos[j]=oddPpos[j-1];       
         inP[j]=OddPrtcls[i].aname;
         inG[j]=inG[j-1];
         inAP[j]=j-1;
         inNum[j]=-OddPrtcls[i].NPDG;
       } else inAP[j]=j;
       j++;
     }

    for(i=0;i<NC;i++) sort[i]=i;
    for(k1=0;k1<NC;k1++) for(k2=0;k2<NC;k2++) inC0[k1*NC+k2]=-1;
    for(k1=0;k1<NC;k1++) for(k2=0;k2<NC;k2++) if(inC0[k1*NC+k2]==-1)
    {  int kk1=inAP[k1];
       int kk2=inAP[k2];
       inC0[k1*NC+k2]=1;
       if(inC0[k2*NC+k1]==-1)   {inC0[k2*NC+k1]=0;   inC0[k1*NC+k2]++;}
       if(inC0[kk1*NC+kk2]==-1) {inC0[kk1*NC+kk2]=0; inC0[k1*NC+k2]++;}
       if(inC0[kk2*NC+kk1]==-1) {inC0[kk2*NC+kk1]=0; inC0[k1*NC+k2]++;}
    }

    for(k1=0;k1<NC;k1++) for(k2=0;k2<NC;k2++) 
    { int L=k1*NC+k2;
       code22_0[L]=NULL;
       code22_1[L]=NULL;

       code22_2[L]=NULL;
       code22Aux0[L]=NULL;
       code22Aux1[L]=NULL;
       code22Aux2[L]=NULL;
       inC1[L]=inC0[L];
       inC2[L]=inC0[L];
    }    

    for(i=0,j=0;i<Nodd;i++)
    {
       inMassAddress[j]=varAddress(OddPrtcls[i].mass);
       if(!inMassAddress[j]) 
       { if(strcmp(OddPrtcls[i].mass ,"0")==0)
         { printf("Error: odd particle '%s' has zero mass.\n",OddPrtcls[i].name);
//           exit(5);
          inMassAddress[j]=&Zero;
         }  
         printf(" Model is not self-consistent:\n "
                " Mass identifier '%s' for particle '%s' is absent  among parameetrs\n",OddPrtcls[i].mass, OddPrtcls[i].name);
//         exit(5);
       }

       if(strcmp(OddPrtcls[i].name,OddPrtcls[i].aname))
       {
         j++;
         inMassAddress[j]=inMassAddress[j-1];
       }
       j++;
    }
    for(Qaddress=NULL,i=0;i<nModelVars;i++) if(strcmp(varNames[i],"Q")==0) Qaddress=varValues+i;
  }

  cleanDecayTable(); 
  if(Qaddress) *Qaddress=100;
  
   
  if(Nodd==0) { printf("No odd particles in the model\n"); return -1; }

  Mcdm=fabs(*(inMassAddress[0]));
  for(i=0;i<NC;i++) 
  { inMass[i]=fabs(*(inMassAddress[i]));
    if(Mcdm>inMass[i]) Mcdm=inMass[i];
  }
  
  if(Qaddress) 
  { *Qaddress=2*Mcdm;
     err=calcMainFunc();
     if(err>0) return err;
  }
            
  GGscale=2*Mcdm/3;

  for(i=0; i<NC-1;)
   {  int i1=i+1;
      if(inMass[sort[i]] > inMass[sort[i1]])
      { int c=sort[i]; sort[i]=sort[i1]; sort[i1]=c;
        if(i) i--; else i++;
      } else i++;
   }

 
  for(i=0;i<NC;i++) 
  {
    if(Z4ch(inP[i])==1) { if(!CDM1) { Mcdm1=inMass[i]; CDM1=inP[i];} else if(Mcdm1>inMass[i]) { Mcdm1=inMass[i];CDM1=inP[i];} }
    if(Z4ch(inP[i])==2) { if(!CDM2) { Mcdm2=inMass[i]; CDM2=inP[i];} else if(Mcdm2>inMass[i]) { Mcdm2=inMass[i];CDM2=inP[i];} }
  }
  
  if(CDM1 && CDM2) for(i=0;i<Nodd;i++) if(Z4ch(OddPrtcls[i].name) != Z4ch(OddPrtcls[i].aname))
  { if(Mcdm1>Mcdm2) { Mcdm1=Mcdm2; CDM1=CDM2;} 
    CDM2=NULL;
    Mcdm2=0;
  }  
  if(CDM1){for(i=0;i<Nodd;i++) if(OddPrtcls[i].name==CDM1) { aCDM1=OddPrtcls[i].aname; break;}}
  else aCDM1=NULL;
  if(CDM2){for(i=0;i<Nodd;i++) if(OddPrtcls[i].name==CDM2) { aCDM2=OddPrtcls[i].aname; break;}}
  else aCDM2=NULL;
  
            
  if(CDM1){strcpy(CDM1_,CDM1); i=strlen(CDM1); } else i=0;
  for(;i<20;i++) CDM1_[i]=' ';
  if(CDM2){strcpy(CDM2_,CDM2); i=strlen(CDM2); } else i=0;
  for(;i<20;i++) CDM2_[i]=' ';

  LSP=sort[0];
  Mcdm=inMass[LSP];

  for(k1=0;k1<NC;k1++)  for(k2=0;k2<NC;k2++) 
  {  if(code22_0[k1*NC+k2]) code22_0[k1*NC+k2]->init=0;
     if(code22_1[k1*NC+k2]) code22_1[k1*NC+k2]->init=0;
     if(code22_2[k1*NC+k2]) code22_2[k1*NC+k2]->init=0;
  }             
  cleanDecayTable();

  for(k1=0;k1<NC;k1++) for(k2=0;k2<NC;k2++) if( code22_0[k1*NC+k2])
  { int  nprc=code22_0[k1*NC+k2]->interface->nprc;
    processAux prc=code22Aux0[k1*NC+k2];
    int n;
    for(n=0;n<=nprc;n++)
    { 
        prc[n].w[0]=0;
        prc[n].w[1]=0;
        prc[n].virt=0;
        prc[n].i3=0;  
        prc[n].br=0;  
        prc[n].cc23=NULL;
        prc[n].virt=0; 
        
        prc[n].nTab=0;   
        if(prc[n].pcmTab) { free(prc[n].pcmTab); prc[n].pcmTab=NULL;}
        if(prc[n].csTab)  { free(prc[n].csTab);  prc[n].csTab=NULL; }  
    }     
  }

  for(k1=0;k1<NC;k1++) for(k2=0;k2<NC;k2++) if( code22_1[k1*NC+k2])
  { int  nprc=code22_1[k1*NC+k2]->interface->nprc;
    processAux prc=code22Aux1[k1*NC+k2];
    int n;
    for(n=0;n<=nprc;n++)
    { 
        prc[n].w[0]=0;
        prc[n].w[1]=0;
        prc[n].virt=0;
        prc[n].i3=0;  
        prc[n].br=0;  
        prc[n].cc23=NULL;
        prc[n].virt=0; 
        
        prc[n].nTab=0;   
        if(prc[n].pcmTab) { free(prc[n].pcmTab); prc[n].pcmTab=NULL;}
        if(prc[n].csTab)  { free(prc[n].csTab);  prc[n].csTab=NULL; }  
    }     
  }

  
  return 0;
}


int toFeebleList(char*name)
{ int i;
  if(!feebleParticles) 
  {  feebleParticles=(int*)malloc(sizeof(int)*nModelParticles);
     for(i=0;i<nModelParticles;i++) feebleParticles[i]=0;
  }
  if(name==NULL) { for(i=0;i<nModelParticles;i++) feebleParticles[i]=0;  return 0;}
  char name_[20];
  strcpy(name_,name);
  trim(name_);
  for(i=0;i<nModelParticles;i++) 
  if(strcmp(name_,ModelPrtcls[i].name)==0 || strcmp(name_,ModelPrtcls[i].aname)==0 ) { feebleParticles[i]=1;  return 0;}
  printf(" particle \"%s\" is out of particle list\n", name_); return 1;
}


int isFeeble(char*name)
{ 
  int p=pTabPos(name); 
  if(p==0) return 0; 
  return feebleParticles[abs(p)-1];  
}
/*
double aWidth(char * pName)
{  txtList LL;  
   return pWidth(pName,&LL);
}
*/
 
int sortOddParticles(char * lsp)
{ int i,err,nprc;

  nprc=sysconf(_SC_NPROCESSORS_ONLN);
  if(nPROCSS>nprc) nPROCSS=nprc;
  err=loadHeffGeff(NULL);

  if(!modelNum)
  {
    int i,k,L;
    struct utsname buff;

    L=strlen(WORK);
    modelDir=malloc(L+15);  sprintf(modelDir,"%s/models",WORK);
    modelNum=1;
  
    calchepDir=malloc(strlen(rootDir)+1);strcpy(calchepDir,rootDir);
    uname(&buff);
    compDir=malloc(strlen(WORK)+strlen(buff.nodename)+25);  
    strcpy(compDir,WORK);
    sprintf(compDir+strlen(compDir),"/_%s_%d_",buff.nodename,getpid());
      
    libDir=malloc(L+15); sprintf(libDir,"%s/so_generated",WORK);
    maxPlistLen=0;
    for(i=0;i<nModelParticles;i++)
    {
       maxPlistLen+=strlen(ModelPrtcls[i].name)+1;
       if(strcmp(ModelPrtcls[i].name,ModelPrtcls[i].aname))
       maxPlistLen+=strlen(ModelPrtcls[i].aname)+1;
    }    
  }
   
  if(omegaCh)   {free(omegaCh);   omegaCh=NULL;}
  if(vSigmaTCh) {free(vSigmaTCh); vSigmaTCh=NULL;}
  if(vSigmaCh)  {free(vSigmaCh);  vSigmaCh=NULL; } 
  err=testSubprocesses();
  if(err)
  { 
    if(lsp)
    {
      if(err>0) {strcpy(lsp,varNames[err]); printf("can not calculate parameter %s\n",varNames[err]);}
       else strcpy(lsp,"Nodd=0");
       printf("sortOddparticles err=%d\n",err);
    }    
    return err;
  }

  if(lsp) strcpy(lsp,inP[LSP]);
  return 0;
}


char * OddParticles(int mode)
{ 
  static char * out[3]={NULL,NULL,NULL};
  int i,len;
  if(mode<0||mode>2) return NULL;
  
  if(out[mode]) return out[mode];

  for(i=0,len=0;i<Nodd;i++) 
  { if(mode==0 || Z4ch(OddPrtcls[i].name)==mode) 
    len+=strlen(OddPrtcls[i].name)+1;
    if(strcmp(OddPrtcls[i].name,OddPrtcls[i].aname))
    len+=strlen(OddPrtcls[i].aname)+1;
  }
  out[mode]=malloc(len); out[mode][0]=0;

  for(i=0;i<Nodd;i++) if(mode==0 || Z4ch(OddPrtcls[i].name)==mode)
  { 
    if(out[mode][0]) strcat(out[mode],","); 
    strcat(out[mode],OddPrtcls[i].name);
    if(strcmp(OddPrtcls[i].name,OddPrtcls[i].aname))
    { strcat(out[mode],",");
      strcat(out[mode],OddPrtcls[i].aname);
    } 
  }
  return out[mode];
}


char * EvenParticles(void)
{ 
  static char * out=NULL;
  int i,len;
  if(out) return out;
  for(i=0,len=0;i<nModelParticles;i++) if(ModelPrtcls[i].name[0]!='~')
  {
    len+=strlen(ModelPrtcls[i].name)+1;
    if(strcmp(ModelPrtcls[i].name,ModelPrtcls[i].aname))
    len+=strlen(ModelPrtcls[i].aname)+1;
  }
  out=malloc(len); out[0]=0;
  for(i=0;i<nModelParticles;i++)if(ModelPrtcls[i].name[0]!='~') 
  { 
    if(out[0]) strcat(out,",");
    strcat(out,ModelPrtcls[i].name);
    if(strcmp(ModelPrtcls[i].name,ModelPrtcls[i].aname))
    { strcat(out,",");
      strcat(out,ModelPrtcls[i].aname);
    } 
  }
  return out;
}


static int  new_code(int k1,int k2, int ch)
{
   char lib[40];
   char process[4000];
   char lib1[12],lib2[12];
   numout*cc;
  
   pname2lib(inP[k1],lib1);
   pname2lib(inP[k2],lib2); 
   sprintf(lib,"omg_%s%s",lib1,lib2);
   switch(ch)
   { case  0: sprintf(process,"%s,%s->AllEven,1*x{%s",inP[k1],inP[k2],EvenParticles()); break;
     case  1: sprintf(process,"%s,%s->AllOdd1,AllOdd1{%s",inP[k1],inP[k2],OddParticles(1));
             strcat(lib,"_1"); break;
     case -1: sprintf(process,"%s,%s->AllOdd2,AllOdd2{%s",inP[k1],inP[k2],OddParticles(2));
             strcat(lib,"_2"); break;
     case  2: sprintf(process,"%s,%s->AllOdd1,AllOdd2{%s{%s",inP[k1],inP[k2],OddParticles(1),OddParticles(2));
              strcat(lib,"_12"); break;                        
   } 
   cc=getMEcode(0,ForceUG,process,NULL,NULL,lib);
   if(cc && cc->interface->nprc) 
   {  int nprc,n;
      processAux prc;
      *(cc->interface->twidth)=1;
      switch(ch)
      { case  0:  code22_0[k1*NC+k2]=cc; break;
        case  1: 
        case -1:  code22_1[k1*NC+k2]=cc; break;
        case  2:  code22_2[k1*NC+k2]=cc; break;
      }  
      nprc=cc->interface->nprc;
      prc=(processAux)malloc((nprc+1)*sizeof(processAuxRec));
      
      switch(ch)
      { case  0: code22Aux0[k1*NC+k2]=prc; break;
        case  1: 
        case -1: code22Aux1[k1*NC+k2]=prc; break;
        case  2: code22Aux2[k1*NC+k2]=prc; break;
      }  
      for(n=0;n<=nprc;n++)
      { 
        prc[n].w[0]=0;
        prc[n].w[1]=0;
        prc[n].virt=0;
        prc[n].i3=0; 
        prc[n].br=0; 
        prc[n].cc23=NULL;
        prc[n].nTab=0;   
        prc[n].pcmTab=NULL;
        prc[n].csTab=NULL;                          
      }   
   } else 
   {  switch(ch)
      { case  0:  inC0[k1*NC+k2]=0; break;
        case  1:  
        case -1: /*inC1[k1*NC+k2]=0;*/ break;
        case  2:  inC2[k1*NC+k2]=0; break;
      }
   }   
   return 0;
}



   
static void gaussC2(double * c, double * x, double * f)
{
  double  A[2][2];
  int i,j;
  double det;
  double B[2];
    
  for(i=0;i<2;i++)
  { int l=1; for(j=0;j<2;j++) {A[i][j]=l*c[i+j]; l=-l;}  
     B[i]=l*c[2+i];
  }
  
  det=A[0][0]*A[1][1] - A[0][1]*A[1][0];
  
  f[0]= ( B[0]*A[1][1]-B[1]*A[0][1])/det;
  f[1]= (-B[0]*A[1][0]+B[1]*A[0][0])/det;
  
  det=sqrt(f[0]+f[1]*f[1]/4.);
   
  x[0]= -f[1]/2.-det;
  x[1]= -f[1]/2.+det;
 
  for(i=0;i<2;i++) { B[i]=c[i]; A[0][i]=1; }
  for(j=0;j<2;j++)   A[1][j]=A[0][j]*x[j];

  det= A[0][0]*A[1][1] - A[0][1]*A[1][0];
  
  f[0]= ( B[0]*A[1][1]-B[1]*A[0][1])/det;
  f[1]= (-B[0]*A[1][0]+B[1]*A[0][0])/det;
} 


static double Logvcs22(double Pcm, void*no) 
{ numout*cc=no;
  int err;
  REAL m[4];
  for(int i=0;i<4;i++) cc->interface->pinf(nsub22,i+1,m+i,NULL); 
  double ss=sqrt( m[0]*m[0]+Pcm*Pcm)+ sqrt( m[1]*m[1]+Pcm*Pcm);
  if(ss<=Fabs(m[2])+Fabs(m[3])) return NAN;
  double cs=cs22(cc,nsub22,Pcm,-1,1,&err)/3.8937966E8;
  return log(cs*Pcm);
}  


static double aRate(double X, int average,int Fast, double * alpha, aChannel ** wPrc,int *NPrc)
{
  double Sum=0.;
  double Sum1=0;
  int i,l1,l2;
  int nPrc=0;
  char* pname[5];
  vGridStr vgrid,vgrid1;  
  double MassCutOut=MassCut+Mcdm*log(1000.)/X;
  double Msmall,Mlarge;

  int nPrcTot=0;
  if(MassCutOut<Mcdm*(2+10/X)) MassCutOut=Mcdm*(2+10/X); 
  WIDTH_FOR_OMEGA=1;

  T_=Mcdm/X;
  s3f_ = s3_T(T_);
  exi=average;

  if(wPrc) *wPrc=NULL;

  for(l1=0;l1<NC;l1++)  if(!feebleParticles[oddPpos[sort[l1]]] && Mcdm+inMass[sort[l1]]<MassCut )
  for(l2=0;l2<NC;l2++)  if(!feebleParticles[oddPpos[sort[l2]]] && inMass[sort[l1]]+inMass[sort[l2]]<MassCut) 
  {    
    double Sumkk=0.;
    double Sum1kk=0;
    double x[2],f[2];
    double factor;
    int kk,k1=sort[l1],k2=sort[l2];
    CalcHEP_interface * CI;

    if(inC0[k1*NC+k2]<=0) continue;
    if(code22_0[k1*NC+k2]==NULL) new_code(k1,k2,0);
    if(inC0[k1*NC+k2]<=0) continue;
    if(!code22_0[k1*NC+k2]->init)
    {
      if(Qaddress && *Qaddress!=inMass[k1]+inMass[k2]) 
      { *Qaddress=inMass[k1]+inMass[k2];
         calcMainFunc();
      }       
      if(passParameters(code22_0[k1*NC+k2])) {FError=1; WIDTH_FOR_OMEGA=0;  return -1;}
      code22_0[k1*NC+k2]->init=1;
    }
    if(wPrc)
    {  nPrcTot+=code22_0[k1*NC+k2]->interface->nprc;
       *wPrc=(aChannel*)realloc(*wPrc,sizeof(aChannel)*(nPrcTot+1));
    }

    sqme22=code22_0[k1*NC+k2]->interface->sqme;
    inBuff=0;

    M1=inMass[k1];
    M2=inMass[k2];


    Msmall=M1>M2? M1-Mcdm*(1-sWidth): M2-Mcdm*(1-sWidth);
    Mlarge=M1>M2? M2+Mcdm*(1-sWidth): M1+Mcdm*(1-sWidth);

    v_min=m2v(MassCutOut);
    if(v_min<1E-200) v_min=1E-200; 

    factor=inC0[k1*NC+k2]*inG[k1]*inG[k2]*exp(-(M1+M2 -2*Mcdm)/T_);
    CI=code22_0[k1*NC+k2]->interface;
    AUX=code22Aux0[k1*NC+k2];
    for(nsub22=1; nsub22<= CI->nprc;nsub22++,nPrc++)
    { double smin;
      double a=0;
      double K=0;
      for(i=0;i<4;i++) pname[i]=CI->pinf(nsub22,i+1,pmass+i,pdg+i);
      if(isFeeble(pname[2]) || isFeeble(pname[3])) continue;
      if(wPrc) 
      { (*wPrc)[nPrc].weight=0;
        for(i=0;i<4;i++) (*wPrc)[nPrc].prtcl[i]=pname[i];
        (*wPrc)[nPrc].prtcl[4]=NULL;
      }

      if(pmass[0]<Mcdm/2 || pmass[1]<Mcdm/2) continue; 
      
      smin=pmass[2]+pmass[3];
      cc23=NULL;
      
      if(VZdecay||VWdecay)
      {  int l,l_,nVV;        

         if(!AUX[nsub22].virt )  for(l=2;l<4;l++) if(pdg[l]==21 ||pdg[l]==22) { AUX[nsub22].virt=-1; break;}
         
         if(!AUX[nsub22].virt)
         {  int vd[4]={0,0,0,0};
            int c_a =  (pmass[0]>Mcdm) || (pmass[1]>Mcdm);

            if(c_a){ for(l=2;l<4;l++) if((pdg[l]==23 && VZdecay>1)   || (abs(pdg[l])==24 && VWdecay>1)) vd[l]=1;} 
            else    for(l=2;l<4;l++)
            { 
            
              if((pdg[l]==23 && VZdecay)     || (abs(pdg[l])==24 && VWdecay)) vd[l]=1;
            } 

            for(l=2;l<4;l++) if(vd[l]) break; 
            if(l<4)
            {  l_=5-l; 
               if(vd[l_])
               { nVV=2;
                 if(pmass[l_]>pmass[l]) { l=l_; l_=5-l;}
               } else nVV=1; 
               AUX[nsub22].virt=l;  
               AUX[nsub22].w[l-2]=pWidth(pname[l],NULL);
               if(abs(pdg[l_])>16 && pmass[l_]> 2) AUX[nsub22].w[l_-2]=pWidth(pname[l_],NULL);
               if(AUX[nsub22].w[l_-2] < 0.1) AUX[nsub22].w[l_-2]=0;
            } else  AUX[nsub22].virt=-1;
         }

         if(AUX[nsub22].virt>0)
         {  l=AUX[nsub22].virt;
            l_=5-l; 
            if( (pmass[0]+pmass[1] < smin+ 4*AUX[nsub22].w[l-2])   && (pmass[l_]< MassCutOut))
            {
              if(AUX[nsub22].cc23) cc23=AUX[nsub22].cc23; else
              {  double  brV1,wV1;
                 int i3W;
                 AUX[nsub22].cc23=xVtoxll(2,2,pname,pdg, l, &wV1, &brV1);
                 if(pdg[l]==pdg[l_]) brV1*=2;
                 AUX[nsub22].br=brV1;
                 cc23=AUX[nsub22].cc23; 
                 if(cc23)
                 {   double Pcm0,PcmMax;
                    *(cc23->interface->BWrange)=10000; 
                    *(cc23->interface->gswidth)=1;  
                    for(i3W=2;i3W<5;i3W++) if(strcmp(cc23->interface->pinf(1,i3W+1,NULL,NULL),pname[l_])==0) break; 
                    AUX[nsub22].i3=i3W;   
                    PcmMax=decayPcm(pmass[2]+pmass[3]+10*AUX[nsub22].w[l-2], pmass[0],pmass[1]);
                    AUX[nsub22].pcm23max=PcmMax;
                    if(Fast>=0)
                    {
                        if( pmass[0]+pmass[1]>=pmass[l_]) Pcm0=(pmass[0]+pmass[1])*0.001; else 
                        Pcm0=decayPcm(pmass[l_]+1E-3, pmass[0],pmass[1]);
                        for(double csTest=sigma23(Pcm0); !isfinite(csTest); Pcm0*=1.01); 
                        buildInterpolation(sigma23, Pcm0,PcmMax, 1E-2 /* -5E-3*/, 0.0001, &(AUX[nsub22].nTab), &(AUX[nsub22].pcmTab), &(AUX[nsub22].csTab));
//#define TEST23
#ifdef TEST23
    char proc[100]; 
    sprintf(proc,"Log(sigma23) for %s,%s -> %s,%s; TabDim=%d", pname[0], pname[1],pname[2],pname[3],AUX[nsub22].nTab);  
    polintArgStr Arg;
    Arg.dim=AUX[nsub22].nTab; Arg.x=AUX[nsub22].pcmTab; Arg.y=AUX[nsub22].csTab; 
    displayPlot(proc,"Pcm[GeV]", Pcm0,PcmMax,0,3,"orig.",0,sigma23,NULL,"interpol.",0,polint_arg,&Arg,"22",0, Logvcs22, code22_0[k1*NC+k2]);    
#endif  
                   }   
                 }
              }    
              if(cc23){ smin=pmass[l_]+0.1;}
            } 
         }
      }
      
      
//if(cc23)  printf("23  %s %s -> %s %s\n", pname[0],pname[1],pname[2],pname[3]);
    
     
//if(abs(pdg[2])!=24 && abs(pdg[3])!=24) continue; 

      if(cc23==NULL) 
      {                               
         if( (pmass[2]>Mlarge && pmass[3]<Msmall)
           ||(pmass[3]>Mlarge && pmass[2]<Msmall))
            { *(CI->twidth)=1; *(CI->gtwidth)=1;} else { *(CI->twidth)=0; *(CI->gtwidth)=0;}
      } 
      *(CI->gswidth)=0;

//*(CI->twidth)=0; *(CI->gtwidth)=0;

      if(smin > pmass[0]+pmass[1])
      { 
        if((pmass[0]!=M1 || pmass[1]!=M2)&&(pmass[0]!=M2 || pmass[1]!=M1))
        { double ms=pmass[0]+pmass[1];
          double md=pmass[0]-pmass[1];
          double Pcm=sqrt((smin-ms)*(smin+ms)*(smin-md)*(smin+md))/(2*smin);
          smin=sqrt(M1*M1+Pcm*Pcm)+sqrt(M2*M2+Pcm*Pcm);
        }      
      }

      if(pmass[0]+pmass[1]> smin) smin=pmass[0]+pmass[1];
      v_max=m2v(smin); 
//printf("v_min=%E v_max=%E  smin=%e\n", v_min,v_max,smin);      
      if(v_max<=v_min) continue; 
            
repeat:
      neg_cs_flag=0;
       
      if(Fast==0)
      { int err; 
        a=simpson(s_integrand,v_min,v_max,eps,&err);
        if(err) { do_err=do_err|err; printf("error in simpson omega.c line 1004\n");}
      } else
      {
          int isPole=0;
          char * s;
          int m,w,n;
          double mass,width;

          for(n=1;(s=code22_0[k1*NC+k2]->interface->den_info(nsub22,n,&m,&w,NULL));n++)
          if(m && w && strcmp(s,"\1\2")==0 )
          { mass=fabs(code22_0[k1*NC+k2]->interface->va[m]);
            width=code22_0[k1*NC+k2]->interface->va[w];
            if(mass<MassCutOut && mass+8*width > pmass[0]+pmass[1]
                            && mass+8*width > smin)
            { if((pmass[0]!=M1 || pmass[1]!=M2)&&(pmass[0]!=M2 || pmass[1]!=M1))
              { double ms=pmass[0]+pmass[1];
                double md=pmass[0]-pmass[1];
                double Pcm=sqrt((mass-ms)*(mass+ms)*(mass-md)*(mass+md))/(2*mass);
                mass=sqrt(M1*M1+Pcm*Pcm)+sqrt(M2*M2+Pcm*Pcm);
              }
              vgrid1=makeVGrid(mass,width);
              if(isPole) vgrid=crossVGrids(&vgrid,&vgrid1); else vgrid=vgrid1;
              isPole++;
            }
          }
                    
          if(cc23)
          {  double mass,width;
             mass=pmass[2]+pmass[3];
             width= AUX[nsub22].w[0]+ AUX[nsub22].w[1];
             if(mass-width>M1+M2)
             { if(Fast_==1)  vgrid1=makeVGrid(mass,width); else  vgrid1=makeVGrid2(mass,width);
               if(isPole) vgrid=crossVGrids(&vgrid,&vgrid1); else vgrid=vgrid1;
               isPole++; 
             }
          } 
          if(isPole==0)
          {  vgrid.n=1;
             vgrid.v[0]=v_min;
             vgrid.v[1]=v_max;
             vgrid.pow[0]=5;
          }
//for(i=0;i<grid.n;i++) printf(" (%E %E) ",grid.ul[i],grid.ur[i]); printf("\n");
/*          if(grid.n==1 && pmass[0]+pmass[1]> 1.1*(smin))
                a=f[0]*sigma(x[0])+f[1]*sigma(x[1]);
          else
*/          
           for(i=0;i<vgrid.n;i++)
          {  
             double da;
             if(Fast>0)   da=gauss(s_integrand,vgrid.v[i],vgrid.v[i+1],vgrid.pow[i]);
             else         { int err; da=simpson(s_integrand,vgrid.v[i],vgrid.v[i+1],eps,&err);
                            if(err) {do_err=do_err||err; printf("error in simpson omega.c line 1066\n");}
                          }   
             a+=da;             
          }
      } 

      if(neg_cs_flag && *(CI->gswidth)==0)
      { *(CI->gswidth)=1;
         goto  repeat;
      }   

// printf("X=%.2E (%d) %.3E %s %s %s %s %E %E %E %E\n",X,average, a, pname[0],pname[1],pname[2],pname[3], pMass(pname[0]),M1,pMass(pname[1]),M2);


      for(kk=2;kk<4;kk++) if(pmass[kk]>2*Mcdm && pname[kk][0]!='~')
      {  txtList LL;
         double BrSm=1;
         pWidth(pname[kk],&LL);
         for(;LL;LL=LL->next)
         { double b;
           char proc[40];
           sscanf(LL->txt,"%lf %[^\n]",&b,proc);
           if( strchr(proc,'~'))BrSm-=b;
         }
         a*=BrSm;
      }
      
      if(pname[2][0]=='~' || pname[3][0]=='~' ) { a/=2; Sum1kk+=a;}  
      Sumkk+=a;
      if(wPrc) (*wPrc)[nPrc].weight = a*factor;
    }
    Sum+=factor*Sumkk;
    Sum1+=factor*Sum1kk;
    
/*
printf("Sum=%E\n",Sum);
*/
  }
   
  if(wPrc) 
  { 
//    for(i=0; i<nPrc;i++) printf("wPrc  %s %s -> %s %s %E\n",(*wPrc)[i].prtcl[0],(*wPrc)[i].prtcl[1],(*wPrc)[i].prtcl[2],
//    (*wPrc)[i].prtcl[3],(*wPrc)[i].weight);

    for(i=0; i<nPrc;i++)  (*wPrc)[i].weight/=Sum;    
    for(i=0;i<nPrc-1;)
    {  if((*wPrc)[i].weight >= (*wPrc)[i+1].weight) i++; 
       else
       {  aChannel buff;
          buff=(*wPrc)[i+1];(*wPrc)[i+1]=(*wPrc)[i];(*wPrc)[i]=buff;
          if(i)i--;else i++;
       }
    }          
    if(NPrc) *NPrc=nPrc; 
//    if(nPrc==0) *wPrc=(aChannel*)realloc(*wPrc,sizeof(aChannel));
       
    if(nPrc){ (*wPrc)[nPrc].weight=0; for(i=0;i<5;i++) (*wPrc)[nPrc].prtcl[i]=NULL;} 
  }  
  if(!average) { double gf=geffDM(Mcdm/X);  Sum/=gf*gf; Sum1/=gf*gf;   }
/*
exit(1);
*/
  WIDTH_FOR_OMEGA=0;
  if(alpha) {  *alpha=Sum1/Sum;  /*  printf("ALPHA=%E\n",*alpha);*/    }   
  return Sum;
}



double vSigma(double T,double Beps ,int Fast)
{
    double X=Mcdm/T;
    double res;
    if(assignVal("Q",2*Mcdm+T)==0) calcMainFunc();
    GGscale=(2*Mcdm+T)/3; 
    MassCut=Mcdm*(2-log(Beps)/X);       
    res= 3.8937966E8*aRate(X, 0 ,Fast,NULL,&vSigmaTCh,NULL);
    return res;
}


double Yeq(double T)
{  double heff;
   double X=Mcdm/T;
   heff=hEff(T); 
   return (45/(4*M_PI*M_PI*M_PI*M_PI))*X*X*geffDM(T)*sqrt(M_PI/(2*X))*exp(-X)/heff;
//   double res= (45/(4*M_PI*M_PI*M_PI*M_PI))*X*X*geffDM(T)*sqrt(M_PI/(2*X))*exp(-X)/heff;

// printf("Yeq(%E) X=%E heff=%E MassCut=%E  geffDM=%E res=%e \n",T,X, heff, MassCut, geffDM(T),res);   

//   return res;
}
                          

struct {double*data; double *alpha; double xtop; int pow,size;} vSigmaGrid={NULL,NULL,0,0,0}; 

static void checkSgridUpdate(void)
{
  if(vSigmaGrid.pow==vSigmaGrid.size)
  { vSigmaGrid.size+=20;
    vSigmaGrid.data=(double*)realloc(vSigmaGrid.data,sizeof(double)*vSigmaGrid.size);
    vSigmaGrid.alpha=(double*)realloc(vSigmaGrid.alpha,sizeof(double)*vSigmaGrid.size);
  }       
}

static double vSigmaI(double T, double Beps, int fast,double * alpha_)
{ double XX,alpha;
  int i,n;
  double X=Mcdm/T;
  if(vSigmaGrid.pow==0)
  { checkSgridUpdate();
    vSigmaGrid.pow=1;
    vSigmaGrid.xtop=X;
    MassCut=Mcdm*(2-log(Beps)/X);
    vSigmaGrid.data[0]= aRate(X,0,fast,&alpha,NULL,NULL)+ZeroCS;
    vSigmaGrid.alpha[0]=alpha;
    if(alpha_) *alpha_=alpha;     
    return vSigmaGrid.data[0];
  }
  
  while(X<vSigmaGrid.xtop*XSTEP)
  { XX=vSigmaGrid.xtop/XSTEP;
    checkSgridUpdate();
    for(i=vSigmaGrid.pow;i;i--)
    { vSigmaGrid.data[i]=vSigmaGrid.data[i-1];
      vSigmaGrid.alpha[i]=vSigmaGrid.alpha[i-1];
    }
    vSigmaGrid.xtop=XX;
    MassCut=Mcdm*(2-log(Beps)/XX);
    vSigmaGrid.data[0]=aRate(XX,0,fast,&alpha,NULL,NULL)+ZeroCS;
    vSigmaGrid.alpha[0]=alpha; 
    vSigmaGrid.pow++;
  }
  
  n=log(X/vSigmaGrid.xtop)/log(XSTEP); 

  while(n+2>vSigmaGrid.pow-1)
  { 
    XX=vSigmaGrid.xtop* pow(XSTEP,vSigmaGrid.pow)  ;
    checkSgridUpdate();
    MassCut=Mcdm*(2-log(Beps)/XX);
    vSigmaGrid.data[vSigmaGrid.pow]=aRate(XX,0,fast,&alpha,NULL,NULL)+ZeroCS;
    vSigmaGrid.alpha[vSigmaGrid.pow]=alpha;
    vSigmaGrid.pow++;
  }

  { double X0,X1,X2,X3,sigmav0,sigmav1,sigmav2,sigmav3,alpha0,alpha1,alpha2,alpha3;
    i=log(X/vSigmaGrid.xtop)/log(XSTEP);
    if(i<0)i=0; 
    if(i>vSigmaGrid.pow-2) i=vSigmaGrid.pow-2;
    X0=vSigmaGrid.xtop*pow(XSTEP,n-1); X1=X0*XSTEP;  X2=X1*XSTEP; X3=X2*XSTEP; 

    sigmav0=log(vSigmaGrid.data[n-1]); alpha0=vSigmaGrid.alpha[n-1]; 
    sigmav1=log(vSigmaGrid.data[n]);   alpha1=vSigmaGrid.alpha[n];
    sigmav2=log(vSigmaGrid.data[n+1]); alpha2=vSigmaGrid.alpha[n+1];
    sigmav3=log(vSigmaGrid.data[n+2]); alpha3=vSigmaGrid.alpha[n+2];
    X=log(X);X0=log(X0); X1=log(X1); X2=log(X2); X3=log(X3);
    
    
    if(alpha_)
    { if(alpha1==0) *alpha_=0; 
      else 
      { *alpha_=  alpha0*       (X-X1)*(X-X2)*(X-X3)/        (X0-X1)/(X0-X2)/(X0-X3)
                 +alpha1*(X-X0)*       (X-X2)*(X-X3)/(X1-X0)/        (X1-X2)/(X1-X3)
                 +alpha2*(X-X0)*(X-X1)*       (X-X3)/(X2-X0)/(X2-X1)/        (X2-X3)
                 +alpha3*(X-X0)*(X-X1)*(X-X2)       /(X3-X0)/(X3-X1)/(X3-X2)        ;
        if(*alpha_ <0) *alpha_=0;
      }                  
    }
    double  res=exp( 
    sigmav0*       (X-X1)*(X-X2)*(X-X3)/        (X0-X1)/(X0-X2)/(X0-X3)
   +sigmav1*(X-X0)*       (X-X2)*(X-X3)/(X1-X0)/        (X1-X2)/(X1-X3) 
   +sigmav2*(X-X0)*(X-X1)*       (X-X3)/(X2-X0)/(X2-X1)/        (X2-X3) 
   +sigmav3*(X-X0)*(X-X1)*(X-X2)       /(X3-X0)/(X3-X1)/(X3-X2)        
                    ) -ZeroCS;
    if(res<0) res=0;
    return res;                             
  }
}


static double dY(double s3, double Beps,double fast)  
{ double d, dlnYds3,Yeq0X, sqrt_gStar, vSig,res;;
  double epsY,alpha;
  double T,heff,geff;
  T=T_s3(s3);   
  heff=hEff(T);
  geff=gEff(T);
  MassCut=2*Mcdm-T*log(Beps);
  d=0.001*s3;  dlnYds3=( log(Yeq(T_s3(s3+d)))- log(Yeq(T_s3(s3-d))) )/(2*d);  // ???

  epsY=deltaY/Yeq(T);

//  sqrt_gStar=polint1(Mcdm/X,Tdim,t_,sqrt_gstar_);

  vSig=vSigmaI(T,Beps, fast,&alpha);
  if(vSig <=0) return 10;
  if(vSig==0){ FError=1; return 0;}
  res= dlnYds3/(pow(2*M_PI*M_PI/45.*heff,0.66666666)/sqrt(8*M_PI/3.*M_PI*M_PI/30.*geff)*vSig*MPlanck
  *(1-alpha/2)*sqrt(1+epsY*epsY))/Yeq(T);
  res=fabs(res);
  if(res>10) return 10;
  return res;
} 


static double darkOmega1(double * Xf,double Z1,double dZ1,int Fast,double Beps)
{
  double X = *Xf;
  double CCX=(Z1-1)*(Z1+1);
  double dCCX=(Z1-1+dZ1)*(Z1+1+dZ1)-CCX;
  double ddY;
  double dCC1,dCC2,X1,X2;
    
  gaussInt= Fast? 1 : 0;

  if(Beps>=1) Beps=0.999;
  vSigmaGrid.pow=0;
  
  ddY=dY(s3_T(Mcdm/X) ,Beps,Fast); 
  if(FError || ddY==0)  return -1;
  if(fabs(CCX-ddY)<dCCX) 
  { *Xf=X; MassCut=Mcdm*(2-log(Beps)/X); 
    return Yeq(Mcdm/X)*sqrt(1+ddY);
  } 
   
  dCC1=dCC2=ddY-CCX; ;X1=X2=X; 
  while(dCC2>0) 
  {  
     X1=X2;
     dCC1=dCC2;
     X2=X2/XSTEP;
     X=X2;
     dCC2=-CCX+dY(s3_T(Mcdm/X),Beps,Fast);
     if(X<2)  return -1;   
     if(Mcdm/X>1.E5) return -1;
  }
             
  while (dCC1<0)
  {  
     X2=X1;
     dCC2=dCC1;
     X1=X1*XSTEP;
     X=X1;
     dCC1=-CCX+dY(s3_T(Mcdm/X),Beps,Fast); 
  }
  for(;;)
  { double dCC;
    if(fabs(dCC1)<dCCX) 
      {*Xf=X1; MassCut=Mcdm*(2-log(Beps)/X1); return Yeq(Mcdm/X1)*sqrt(1+CCX+dCC1);}
    if(fabs(dCC2)<dCCX || fabs(X1-X2)<0.0001*X1) 
      {*Xf=X2; MassCut=Mcdm*(2-log(Beps)/X2); return Yeq(Mcdm/X2)*sqrt(1+CCX+dCC2);}
    X=0.5*(X1+X2); 
    dCC=-CCX+dY(s3_T(Mcdm/X),Beps,Fast);
    if(dCC>0) {dCC1=dCC;X1=X;}  else {dCC2=dCC;X2=X;} 
  }
}

static double Beps_;
int Fast_=1;

static void XderivLn(double s3, double *Y, double *dYdx)
{
  double y=Y[0];
  double yeq, sqrt_gStar;
  double T,heff,geff;
  
//  s3=polint1(T,Tdim,t_,s3_);  
  T=T_s3(s3);  
  heff=hEff(T);
  geff=gEff(T);
//  sqrt_gStar=polint1(T,Tdim,t_,sqrt_gstar_);
  
  MassCut=2*Mcdm -T*log(Beps_); yeq=Yeq(T);
//  if(y<yeq) *dYdx=0; else 
  { double vSig,alpha,epsY;
  
    if(deltaY) epsY=deltaY/y; else  epsY=0; 
    vSig=vSigmaI(T,Beps_,Fast_,&alpha);
//printf("T=%E alpha=%E\n", Mcdm/x, alpha);     
    *dYdx=MPlanck
    *pow(2*M_PI*M_PI/45.*heff,0.666666666666)/sqrt(8*M_PI/3.*M_PI*M_PI/30.*geff)
//    *sqrt_gStar*sqrt(M_PI/45)
    *vSig*(y*y-(1-alpha)*yeq*yeq-alpha*y*yeq)*sqrt(1+epsY*epsY);
//printf(" T=%E  y=%E   yeq=%E  epsY=%E  alpha=%E \n",T, y,  yeq, epsY, alpha);   
  }
}

static struct 
{ 
  char *CDM1, *CDM2,*aCDM1,*aCDM2;
  double  mCDM1, mCDM2,mCDM;
} mem;


static void fillCDMmem(void)
{ 

  mem.CDM1=CDM1;   mem.CDM2=CDM2;   mem.aCDM1=aCDM1; mem.aCDM2=aCDM2; 
  mem.mCDM1=Mcdm1; mem.mCDM2=Mcdm2; mem.mCDM=Mcdm;
  
  CDM1=aCDM1=CDM2=aCDM2=NULL;
  Mcdm1=Mcdm2=0;
  
  for(int i=0;i<NC;i++) if(Z4ch(inP[i])==1 && ( !feebleParticles[oddPpos[i]] && (Mcdm1==0 ||  Mcdm1 > inMass[i]))) 
  {  Mcdm1=inMass[i];
     if(Mcdm1>0)
     {  CDM1=inP[i];
       aCDM1=ModelPrtcls[oddPpos[i]].aname;
     }  
  }   
  
  for(int i=0;i<NC;i++) if(Z4ch(inP[i])==2 && ( !feebleParticles[oddPpos[i]] && (Mcdm2==0 ||  Mcdm2 > inMass[i]))) 
  {  Mcdm2=inMass[i]; 
     if(Mcdm2>0)
     {  CDM2=inP[i];
       aCDM2=ModelPrtcls[oddPpos[i]].aname;
     }
  }                        

  if(CDM1) 
  { Mcdm=Mcdm1;
    if(CDM2 && Mcdm2<Mcdm) Mcdm=Mcdm2;   
  } else Mcdm=Mcdm2;
   
}


static void restoreCDMmem(void)
{
   CDM1=mem.CDM1;   CDM2=mem.CDM2;  aCDM1=mem.aCDM1; aCDM2=mem.aCDM2; 
   Mcdm1=mem.mCDM1; Mcdm2=mem.mCDM2; Mcdm=mem.mCDM;
   MassCut=MPlanck;
}


double darkOmegaFO(double * Xf_, int Fast, double Beps)
{
  double Yf;
  double Z1=2.5;
  double dZ1=0.05;
  double Xf=25;

  if(CDM1==NULL) fracCDM2=1; else
  if(CDM2==NULL) fracCDM2=0; else 
  if(Mcdm1<Mcdm2) fracCDM2=0; else fracCDM2=1;

  fillCDMmem();
  
  if(omegaCh) {free(omegaCh); omegaCh=NULL;}
    
  if(Xf_) *Xf_=Xf; 
  if(assignVal("Q",2*Mcdm)==0) calcMainFunc();
  GGscale=2*Mcdm/3;   
  if(Beps>=1) Beps=0.999;
  
  Yf=  darkOmega1(&Xf, Z1, dZ1,Fast, Beps);
  if(FError||Xf<1||Yf<=0) {  return -1;}
  
  double iColl=( (Mcdm/Xf)*sqrt(M_PI/45)*MPlanck*aRate(Xf, 1,Fast,NULL, NULL,NULL) );

  if(Xf_) *Xf_=Xf; 

  restoreCDMmem();
  if(FError) return -1;
  return  2.742E8*Mcdm/(1/Yf +  iColl); /* 2.828-old 2.755-new 2.742 next-new */
}


static double *lYtab=NULL;
static double *Ttab=NULL;
static int Ntab=0;

double YF(double T){if(Ntab<=0 || T<Ttab[Ntab-1] || T> Ttab[0]) return NAN; return exp(polint3(T,Ntab,Ttab, lYtab)) ;}


double darkOmega(double * Xf, int Fast, double Beps,int *err)
{
  double Yt,Xt=27;
  double Tend_;
  double Z1=1.1,Z2=10,Zf=2.5; 
  int i;
  int Nt=25;
  if(err) *err=0;

  if(CDM1==NULL)  fracCDM2=1; else
  if(CDM2==NULL)  fracCDM2=0; else 
  if(Mcdm1<Mcdm2) fracCDM2=0; else fracCDM2=1;

  fillCDMmem();
  if(Mcdm<=0) { printf(" There are no  Dark Matter particles\n"); restoreCDMmem(); return 0;} 
  
  lYtab=realloc(lYtab,sizeof(double)*Nt);
  Ttab=realloc(Ttab,sizeof(double)*Nt);
  Ntab=0;

  if(assignVal("Q",2*Mcdm)==0) calcMainFunc() ;
  GGscale=2*Mcdm/3;
  if(Beps>=1) Beps=0.999;
  Beps_=Beps; Fast_=Fast;
  
  if(Z1<=1) Z1=1.1;
  
  Yt=  darkOmega1(&Xt, Z1, (Z1-1)/5,Fast, Beps);

  if(Yt<0||FError) 
  {  restoreCDMmem(); 
     if(err) *err=1; else printf("Temperature of thermal  equilibrium is too large\n");
     if(Xf) *Xf=0;   
     return -1;
  }
  
  Tstart=Mcdm/Xt;
  
  if(Yt<fabs(deltaY)*1.E-15)
  {  
     if(deltaY>0) dmAsymm=1;  else dmAsymm=-1;
     if(Xf) *Xf=Xt;
     restoreCDMmem();
     return 2.742E8*Mcdm*deltaY;
  }   
  
  Ntab=1;
  Ttab[0]=Tstart;
  lYtab[0]=log(Yt); 
  Tend_=Tstart;
  
  for(i=0; ;i++)
  { double X2=vSigmaGrid.xtop*pow(XSTEP,i+1);
    double yeq,alpha;
    double s3_t,s3_2;

    if(Xt>X2*0.999999) continue; 

    yeq=Yeq(Mcdm/Xt);
    alpha=vSigmaGrid.alpha[i];    


    if(Yt*Yt>=Z2*Z2*( alpha*Yt*yeq+(1-alpha)*yeq*yeq) || Yt<fabs(deltaY*1E-15))  break;
    

    s3_t=s3_T(Mcdm/Xt);
    s3_2=s3_T(Mcdm/X2); 
//    if(odeint(&y,1 ,Mcdm/Xt , Mcdm/X2 , 1.E-3, (Mcdm/Xt-Mcdm/X2 )/2, &XderivLn)){ printf("problem in solving diff.equation\n"); return -1;}   
    if(odeint(&Yt,1 ,s3_t , s3_2 , 1.E-3, (s3_2-s3_t)/2, &XderivLn)){ printf("problem in solving diff.equation\n"); return -1;}
    if(Ntab>=Nt)
    { Nt+=20;
      lYtab=realloc(lYtab,sizeof(double)*Nt);
      Ttab=realloc(Ttab,sizeof(double)*Nt);
    }      
    
    Tend_=Mcdm/X2;                                 
    Xt=X2;   
    lYtab[Ntab]=log(Yt);
    Ttab[Ntab]=Tend_;
    Ntab++;
  }
  
  if(Xf) 
  {  double T1,T2,Y1,Y2,dY2,dY1;
     T1=Ttab[0];
     Y1=exp(lYtab[0]);
     dY1=Zf*Yeq(T1)-Y1;
     *Xf=Mcdm/T1;
     for(i=1;i<Ntab;i++)            
     { T2=Ttab[i];
       Y2=exp(lYtab[i]); 
       dY2=Zf*Yeq(T2)-Y2;
       if(dY2<0)
       { 
         for(;;)
         {  double al,Tx,Yx,dYx,Xx;
            al=dY2/(dY2-dY1);
            Tx=al*T1+(1-al)*T2, /*Yx=al*Y1+(1-al)*Y2,*/ Yx=exp(polint3(Tx,Ntab,Ttab,lYtab)),    dYx=Zf*Yeq(Tx)-Yx;
            if(fabs(dYx)<0.01*Yx) 
            { *Xf=Mcdm/Tx;
              break;
            } else  { if(dYx>0) {T1=Tx,Y1=Yx;}  else {T2=Tx,Y2=Yx;} }
         }
         break; 
      }  
      else {dY1=dY2; T1=T2; Y1=Y2; *Xf=Mcdm/T2;}        
    }
  }
     
  if(Yt<fabs(deltaY*1E-15))  
  {  
      if(deltaY>0) dmAsymm=1; else dmAsymm=-1;
      restoreCDMmem();
      return 2.742E8*Mcdm*deltaY;
  }  


  double iColl=( (Mcdm/Xt)*sqrt(M_PI/45)*MPlanck*aRate(Xt,1,Fast,NULL,NULL,NULL));
  restoreCDMmem();
  if(FError) { if(err) *err=8;  if(Xf) *Xf=0;     return -1;}
  
  if(deltaY==0)
  { dmAsymm=0;
    return  2.742E8*Mcdm/(1/Yt  + iColl); /* 2.828-old 2.755-new,2.742 -newnew */
  } else
  {  double a,f,z0,Y0;
     a=fabs(deltaY);
     if(Yt<a*1.E-5)  f=Yt*Yt/4/a; else f=(sqrt(Yt*Yt+a*a)-a)/(sqrt(Yt*Yt+a*a)+a);   
     f*= exp(-2*a*iColl);
     z0=sqrt(f)*2*a/(1-f);
     Y0=sqrt(z0*z0+a*a);
     dmAsymm=deltaY/Y0;
     return 2.742E8*Mcdm*Y0;
  }    
}

double darkOmegaTR(double Tr, double Yr,  int Fast, double Beps,int *err)
{
  double Yt;
  int i;
  if(err) *err=0;
  if(CDM1==NULL)  fracCDM2=1; else
  if(CDM2==NULL)  fracCDM2=0; else 
  if(Mcdm1<Mcdm2) fracCDM2=0; else fracCDM2=1;

  fillCDMmem();
  if(Mcdm<=0) { printf(" There are no  Dark Matter particles\n"); restoreCDMmem(); return 0;} 
  
  if(assignVal("Q",2*Mcdm)==0) calcMainFunc() ;
  GGscale=2*Mcdm/3;
  if(Beps>=1) Beps=0.999;
  Beps_=Beps; Fast_=Fast;
  
  Tstart=Tr;
  vSigmaGrid.pow=0;
  double alpha;  vSigmaI(Tr, Beps, 1,&alpha);
  
  Ntab=100;
  lYtab=realloc(lYtab,sizeof(double)*Ntab);
  Ttab=realloc(Ttab,sizeof(double)*Ntab);
  
  Tstart=Tr;


  for(int i=0;i<Ntab;i++) Ttab[i]=Tstart-i*(Tstart-Tend)/(Ntab-1);
  lYtab[0]=log(Yr);
  Yt=Yr;
  for(i=1;i<Ntab-1 ;i++)
  { 
    double yeq,alpha;
    double s3_t,s3_2;
    s3_t=s3_T(Ttab[i-1]);
    s3_2=s3_T(Ttab[i]);

//printf("i=%d  %E => %E \n",i, Ttab[i-1],Ttab[i] );

//    alpha=vSigmaGrid.alpha[i];    

    if(odeint(&Yt,1 ,s3_t , s3_2 , 1.E-3, (s3_2-s3_t)/2, &XderivLn)){ printf("problem in solving diff.equation\n"); return -1;}
//    printf("Y=%E\n",Yt);
    lYtab[i]=log(Yt);
  }
  double Xt=Mcdm/Ttab[Ntab-2];
  double iColl=( (Mcdm/Xt)*sqrt(M_PI/45)*MPlanck*aRate(Xt,1,Fast,NULL,NULL,NULL));
  Yt=1/( 1/Yt+iColl);
  lYtab[Ntab-1]=log(Yt);
  restoreCDMmem();
  dmAsymm=0;
  if(FError) { if(err) *err=8;  return -1;}
  return  2.742E8*Mcdm*Yt;
}     




double printChannels(double Xf ,double cut, double Beps, int prcn, FILE * f)
{ int i,nPrc,nform=log10(1/cut)-2;
  double Sum,s;

  if(!Xf) return -1;
  double Mcdm_mem=Mcdm;
  Mcdm=-1;
  for(i=0;i<NC;i++) if(!feebleParticles[oddPpos[i]] && (Mcdm<0 ||  Mcdm > inMass[i])) Mcdm=inMass[i]; 
  
  if(omegaCh) {free(omegaCh); omegaCh=NULL;}

  MassCut=Mcdm*(2-log(Beps)/Xf);
  Sum=aRate(Xf, 1,1,NULL,&omegaCh,&nPrc)*(Mcdm/Xf)*sqrt(M_PI/45)*MPlanck/(2.742E8*Mcdm_mem); 
  if(Sum==0 || FError)     { return -1;}
  if(nform<0)nform=0;
   
  if(f)
  {  int j;
     fprintf(f,"# Channels which contribute to 1/(omega) more than %G%%.\n",100*cut );
     if(prcn) fprintf(f,"# Relative contributions in %% are displayed\n");
        else  fprintf(f,"# Absolute contributions are displayed\n");
     for(i=0,s=0;i<nPrc;i++)  if(fabs(omegaCh[i].weight)>=cut)
     {  s+=omegaCh[i].weight;
        if(prcn)
        { if(cut <0.000001) fprintf(f,"  %.1E%% ",100*omegaCh[i].weight);
          else              fprintf(f,"  %*.*f%% ",nform+3,nform,
                                        100*omegaCh[i].weight);
        } else fprintf(f,"  %.1E ",Sum*omegaCh[i].weight); 
        for(j=0;j<4;j++)
        {
           fprintf(f,"%s ",omegaCh[i].prtcl[j]);
           if(j==1) fprintf(f,"->");
           if(j==3) fprintf(f,"\n");
        }
     }
  }
  
  Mcdm=Mcdm_mem;
  
  return 1/Sum;
}

static int strcmp_(char * n1, char *n2) { if( n1[0]=='*' &&  n1[1]==0) return 0; return strcmp(n1, n2);}

double oneChannel(double Xf,double Beps,char*n1,char*n2,char*n3,char*n4)
{ int j,nPrc;
  aChannel *wPrc;
  double Sum,res;
  if(!Xf) return -1;
  MassCut=Mcdm*(2-log(Beps)/Xf);
  Sum=aRate(Xf, 1,1,NULL,&wPrc,&nPrc)*(Mcdm/Xf)*sqrt(M_PI/45)*MPlanck/(2.742E8*Mcdm);
  if(FError)     { return -1;}
  if(wPrc==NULL) { return  0;}  

  for(res=0,j=0;j<nPrc;j++) 
  if( ( (strcmp_(n1,wPrc[j].prtcl[0])==0 && strcmp_(n2,wPrc[j].prtcl[1])==0) ||
        (strcmp_(n2,wPrc[j].prtcl[0])==0 && strcmp_(n1,wPrc[j].prtcl[1])==0)
      ) &&
      ( (strcmp_(n3,wPrc[j].prtcl[2])==0 && strcmp_(n4,wPrc[j].prtcl[3])==0) ||
        (strcmp_(n4,wPrc[j].prtcl[2])==0 && strcmp_(n3,wPrc[j].prtcl[3])==0)
      )
    )  {res+=wPrc[j].weight;} 
      
  free(wPrc); 
  return res;
}
//================= Ext vSigma =======

static double vSigmaZero(double T){return 0;}
static double (*vSigmaStat0)(double T)=vSigmaZero;
static double (*vSigmaStat1)(double T)=vSigmaZero;


static void XderivLnExt(double s3, double *Y, double *dYdx)
{
  double y=Y[0];
  double yeq, sqrt_gStar;
  double T,heff,geff;
  double vSig,vSig0,vSig1,alpha,epsY;
  
  T=T_s3(s3);  
  yeq=Yeq(T);
  if(y<=yeq) {*dYdx=0; return;} else epsY=deltaY/y;
  
  
  vSig0=vSigmaStat0(T)/3.8937966E8;
  vSig1=vSigmaStat1(T)/3.8937966E8;
  vSig=(vSig0+vSig1);
    
  if(vSig==0) {*dYdx=0; return;} else alpha=vSig1/(vSig0+vSig1);
  
  heff=hEff(T);
  geff=gEff(T);

  *dYdx=MPlanck
      *pow(2*M_PI*M_PI/45.*heff,0.666666666666)/sqrt(8*M_PI/3.*M_PI*M_PI/30.*geff)
      *vSig*(y*y-(1-alpha)*yeq*yeq-alpha*y*yeq)*sqrt(1+epsY*epsY); 
}


static double dYExt(double s3)  
{ double d, dlnYds3,Yeq0X, sqrt_gStar, vSig,vSig0,vSig1,res;;
  double epsY,alpha,yeq;
  double T,heff,geff;
  T=T_s3(s3);   
  yeq=Yeq(T);
  if(yeq<=0) return 10;
  epsY=deltaY/yeq;
  vSig0=vSigmaStat0(T)/3.8937966E8;
  vSig1=vSigmaStat1(T)/3.8937966E8;
  vSig=vSig0+vSig1;
  if(vSig <=0) return 10; 
  alpha=vSig1/vSig;
  
  
  heff=hEff(T);
  geff=gEff(T);
  d=0.001*s3;  dlnYds3=( log(Yeq(T_s3(s3+d)))
                        -log(Yeq(T_s3(s3-d))) )/(2*d);

  res= dlnYds3/(pow(2*M_PI*M_PI/45.*heff,0.66666666)/sqrt(8*M_PI/3.*M_PI*M_PI/30.*geff)
      *vSig*MPlanck*(1-alpha/2)*sqrt(1+epsY*epsY))/Yeq(T);
  res=fabs(res);
  if(res>10) return 10;
  return res;
} 


static double darkOmega1Ext(double * Xf,double Z1,double dZ1)
{
  double X = *Xf;
  double CCX=(Z1-1)*(Z1+1);
  double dCCX=(Z1-1+dZ1)*(Z1+1+dZ1)-CCX;
  double ddY;
  double dCC1,dCC2,X1,X2;

  
  ddY=dYExt(s3_T(Mcdm/X)); 
  if(FError || ddY==0)  return -1;
  if(fabs(CCX-ddY)<dCCX) 
  { *Xf=X;
    return Yeq(Mcdm/X)*sqrt(1+ddY);
  } 
   
  dCC1=dCC2=ddY-CCX; ;X1=X2=X; 
  while(dCC2>0) 
  {  
     X1=X2;
     dCC1=dCC2;
     X2=X2/XSTEP;
     X=X2;
     dCC2=-CCX+dYExt(s3_T(Mcdm/X));
     if(Mcdm/X>1.E5) return -1;
  }
             
  while (dCC1<0)
  {  
     X2=X1;
     dCC2=dCC1;
     X1=X1*XSTEP;
     X=X1;
     dCC1=-CCX+dYExt(s3_T(Mcdm/X)); 
  }
  for(;;)
  { double dCC;
    if(fabs(dCC1)<dCCX) 
      {*Xf=X1;  return Yeq(Mcdm/X1)*sqrt(1+CCX+dCC1);}
    if(fabs(dCC2)<dCCX || fabs(X1-X2)<0.0001*X1) 
      {*Xf=X2;  return Yeq(Mcdm/X2)*sqrt(1+CCX+dCC2);}
    X=0.5*(X1+X2); 
    dCC=-CCX+dYExt(s3_T(Mcdm/X));
    if(dCC>0) {dCC1=dCC;X1=X;}  else {dCC2=dCC;X2=X;} 
  }
}


double darkOmegaExt(double * Xf,  double (*f0)(double),double (*f1)(double))
{
  double Yt,Xt=25;
  double Z1=1.1;
  double Z2=10,Zf=2.5; 
  int i;
  double Tend_;
  
  if(f0) vSigmaStat0=f0; else vSigmaStat0=vSigmaZero;
  if(f1) vSigmaStat1=f1; else vSigmaStat1=vSigmaZero;
  
  MassCut=4*Mcdm;

  int Nt=25;
  
  lYtab=realloc(lYtab,sizeof(double)*Nt);
  Ttab=realloc(Ttab,sizeof(double)*Nt);
  Ntab=0;

  if(CDM1==NULL) fracCDM2=1; else
  if(CDM2==NULL) fracCDM2=0; else 
  if(Mcdm1<Mcdm2)fracCDM2=0; else fracCDM2=1;

  Yt=  darkOmega1Ext(&Xt, Z1, (Z1-1)/5);

  if(Yt<0||FError) { return -1;}
  
  if(Yt<fabs(deltaY)*1.E-15)
  {  
     if(deltaY>0) dmAsymm=1;  else dmAsymm=-1;
     if(Xf) *Xf=Xt;   
     return 2.742E8*Mcdm*deltaY;
  }     
  
  Tstart=Mcdm/Xt;
  Ntab=1;
  Ttab[0]=Tstart;
  lYtab[0]=log(Yt);
  Tend_=Tstart;
         
  for(i=0; ;i++)
  { 
    double s3_t,s3_2,Tbeg;
    double yeq=Yeq(Mcdm/Xt);

    if(Tend_<1.E-3 || Yt<fabs(deltaY*1E-5)) break;
    Tbeg=Tend_;       
    Tend_/=1.2;                                 
    s3_t=s3_T(Tbeg);
    s3_2=s3_T(Tend_); 
    if(odeint(&Yt,1 ,s3_t , s3_2 , 1.E-3, (s3_2-s3_t)/2, &XderivLnExt)){ printf("problem in solving diff.equation\n"); return -1;}
    if(!isfinite(Yt)||FError)  return -1;
    if(Ntab>=Nt)
    { Nt+=20;   
      lYtab=realloc(lYtab,sizeof(double)*Nt);
      Ttab=realloc(Ttab,sizeof(double)*Nt);
    }      
    lYtab[Ntab]=log(Yt);
    Ttab[Ntab]=Tend_;
    Ntab++;  
    Tbeg=Tend_;                                           
  }  

  if(Xf) 
  {  double T1,T2,Y1,Y2,dY2,dY1;
     T1=Ttab[0];
     Y1=exp(lYtab[0]);
     dY1=Zf*Yeq(T1)-Y1;  
     *Xf=Mcdm/T1;
     for(i=1;i<Ntab;i++)            
     { T2=Ttab[i];
       Y2=exp(lYtab[i]); 
       dY2=Zf*Yeq(T2)-Y2;
       if(dY2<0)
       {          
         for(;;)
         {  double al,Tx,Yx,dYx,Xx;
            al=dY2/(dY2-dY1);
            Tx=al*T1+(1-al)*T2,  /*Yx=al*Y1+(1-al)*Y2,*/ Yx=exp(polint3(Tx,Ntab,Ttab,lYtab)),   dYx=Zf*Yeq(Tx)-Yx;
            if(fabs(dYx)<0.01*Yx) 
            { *Xf=Mcdm/Tx;
              break;
            } else  { if(dYx>0) {T1=Tx,Y1=Yx;}  else {T2=Tx,Y2=Yx;} }
         }
         break; 
      }  
      else {dY1=dY2; T1=T2; Y1=Y2; *Xf=Mcdm/T2;}        
    }
  }

  if(Yt<fabs(deltaY*1E-15))  
  {  
      if(deltaY>0) dmAsymm=1; else dmAsymm=-1;   
      return 2.742E8*Mcdm*deltaY;
  }  
//  Yi=1/( (Mcdm/Xt)*sqrt(M_PI/45)*MPlanck*aRate(Xt,1,Fast,NULL,NULL,NULL));
  if(deltaY==0)
  { dmAsymm=0;
    return  2.742E8*Mcdm*Yt; /* 2.828-old 2.755-new,2.742 -newnew */
  } else
  {  double a,f,z0,Y0;
     if(Yt<fabs(deltaY*1E-15))
     { if(deltaY>0) dmAsymm=1; else dmAsymm=-1;
       return 2.742E8*Mcdm*deltaY;
     }
     a=fabs(deltaY);
     if(Yt<a*1.E-5)  f=Yt*Yt/4/a; else f=(sqrt(Yt*Yt+a*a)-a)/(sqrt(Yt*Yt+a*a)+a);   
     z0=sqrt(f)*2*a/(1-f);
     Y0=sqrt(z0*z0+a*a);
     dmAsymm=deltaY/Y0;
     return 2.742E8*Mcdm*Y0;
  }   
}


/*===========  Z4  ==================*/


static double geff1_(double T)
{ 

  double massCut=Mcdm1;
  if(Beps>0)    massCut-=T*log(Beps); else  massCut+=1.E20;

   double sum=0,t; int l;
   for(l=0;l<NC;l++) if(!feebleParticles[oddPpos[sort[l]]]) 
   { int k=sort[l];
     if(Z4ch(inP[k])==1) 
     { double bsk2; 
       double M=inMass[k];
       if(M>massCut) continue;
       t=T/M;
       if(t<0.1) bsk2=K2pol(t)*exp((Mcdm1-M)/T)*sqrt(M_PI*t/2);
        else     bsk2=bessK2(1/t)*exp(Mcdm1/T);
       sum+=inG[k]*M*M*bsk2;
     }
   }      
   return sum;
}

static double geff2_(double T)
{ 
   double massCut=Mcdm2;
   if(Beps>0)  massCut-=T*log(Beps); else  massCut+=1.E20;
   double sum=0,t; int l;
   for(l=0;l<NC;l++) if(!feebleParticles[oddPpos[sort[l]]])
   { int k=sort[l];
     if(Z4ch(inP[k])==2) 
     { double bsk2; 
       double M=inMass[k];
       if(M>massCut) continue;
       t=T/M;
       if(t<0.1) bsk2=K2pol(t)*exp(-1/t+Mcdm2/T)*sqrt(M_PI*t/2);
        else     bsk2=bessK2(1/t)*exp(Mcdm2/T);
       sum+=inG[k]*M*M*bsk2;
     }
   }  
   return sum;
}


static double McdmSum;
double Beps=1.E-4;

static double s_integrandT(double  sqrtS )
{  double sv_tot,t,bsk1;
   double ms,md,PcmIn;
   double res;
   
   ms = M1 + M2; 
   if(ms>=sqrtS)  return 0;
   md = M1 - M2;
   PcmIn = sqrt((sqrtS-ms)*(sqrtS+ms)*(sqrtS-md)*(sqrtS+md))/(2*sqrtS);
   sv_tot=sigma(PcmIn);         
   t=T_/sqrtS; 
   if(t<0.1) bsk1=K1pol(t)*exp(-1/t+McdmSum/T_)*sqrt(M_PI*t/2);
   else      bsk1=bessK1(sqrtS/T_)*exp(McdmSum/T_);
      
   res= sqrtS*sqrtS*(PcmIn*PcmIn)*sv_tot*bsk1/T_;
   return res;
}

static double v_integrand( double u)
{  double y,sv_tot,w;
   double Xf_1;
   double ms,md,sqrtS,PcmIn,res0;
   
   if(u==0. || u==1.) return 0.;

   long double u_=u,z=u_*(2-u_);
   sqrtS=M1+M2-3*T_*logl(z);
   
   return s_integrandT(sqrtS )*6*T_*(1-u)/z;  
}



static  int aRate4(double T, 
  double*vs1100_,double*vs2200_,
  double*vs1110_,double*vs1120_,double*vs1210_,double*vs1220_,double*vs2210_,double*vs2220_,
  double*vs1112_,double*vs1122_,double*vs1222_,double*vs1211_,double*vs2211_,double*vs2221_)
{
  int i,l1,l2,N12;
  char* pname[5];
  double X=Mcdm/T;
 
  double Msmall,Mlarge;
  double vs1100=0, vs2200=0, vs1110=0,vs2220=0,vs1120=0,vs1210=0,vs1122=0,vs2211=0,
         vs1112=0, vs1222=0, vs1220=0,vs2210=0,vs2221=0,vs1211=0;
     
  vGridStr vgrid,vgrid1;
   
  WIDTH_FOR_OMEGA=1;
  T_=T;

for(N12=0;N12<=2;N12++)
{ 
  if( N12==1 && (!CDM1 || !CDM2) ) break;
  for(l1=0;l1<NC;l1++) if(!feebleParticles[oddPpos[sort[l1]]])
  { int k1=sort[l1];
  for(l2=0;l2<NC;l2++) if(!feebleParticles[oddPpos[sort[l2]]])
  {
    double factor;
    int k2=sort[l2];
    numout * code=NULL;
    CalcHEP_interface * CI;
    double MassCutOut;
    int*inC;
    
    switch(N12)
    { case 0: inC=inC0; break;
      case 1: inC=inC1; break;
      case 2: inC=inC2; break;
    }  

    if(inC[k1*NC+k2]<=0) continue;

    MassCut=0;    
    if(Z4ch(inP[k1])==1) MassCut+=Mcdm1; else  MassCut+=Mcdm2;
    if(Z4ch(inP[k2])==1) MassCut+=Mcdm1; else  MassCut+=Mcdm2;
    McdmSum=MassCut;
    if(Beps>0)    MassCut-=T*log(Beps); else  MassCut+=1.E20;
    
    MassCutOut=MassCut+T*log(10000.); 
    
    if(inMass[k1]+inMass[k2] > MassCut)
    {
      continue;
    }
    switch(N12)
    { case 0:  
         if(code22_0[k1*NC+k2]==NULL) new_code(k1,k2,0);
         code=code22_0[k1*NC+k2];  break;  
      case 1:                   
       if(Mcdm1 > Mcdm2)
       { if(Z4ch(inP[k1])==2 &&  Z4ch(inP[k2])==2) continue;
         if(code22_1[k1*NC+k2]==NULL) new_code(k1,k2,-1);
         code=code22_1[k1*NC+k2];
       } else 
       {  if(Z4ch(inP[k1])==1 && Z4ch(inP[k2])==1) continue;
          if(code22_1[k1*NC+k2]==NULL) new_code(k1,k2,1);
          code=code22_1[k1*NC+k2];
       } break;
      case 2:
        if(Mcdm1>Mcdm2)
        {   if(Z4ch(inP[k1])==1 && Z4ch(inP[k2])==1) { if(code22_2[k1*NC+k2]==NULL) new_code(k1,k2,2); code=code22_2[k1*NC+k2];}
            else  continue;
        } else 
        {   if(Z4ch(inP[k1])==2 && Z4ch(inP[k2])==2) { if(code22_2[k1*NC+k2]==NULL) new_code(k1,k2,2); code=code22_2[k1*NC+k2];}   
            else  continue;
        }
        break; 
    }
    
    if(!code) continue;
    if(!code->init)
    { 
      CalcHEP_interface *cdi=code->interface;
      for(i=1;i<=cdi->nvar;i++) if(code->link[i]) cdi->va[i]=*(code->link[i]);     
      if(  cdi->calcFunc()>0 ) {FError=1; WIDTH_FOR_OMEGA=0;  return -1;}
      code->init=1;
    }
    
    sqme22=code->interface->sqme;
    
    inBuff=0;

    M1=inMass[k1];
    M2=inMass[k2];

    Msmall=M1>M2? M1-Mcdm*(1-sWidth): M2-Mcdm*(1-sWidth);
    Mlarge=M1>M2? M2+Mcdm*(1-sWidth): M1+Mcdm*(1-sWidth);

    v_min=m2v(MassCutOut);
    if(v_min< 1E-200) v_min=1E-200;
    
    factor=inC[k1*NC+k2]*inG[k1]*inG[k2];

    CI=code->interface;
    switch(N12)
    { case 0:  AUX=code22Aux0[k1*NC+k2]; break;
      case 1:  AUX=code22Aux1[k1*NC+k2]; break;
      case 2:  AUX=code22Aux2[k1*NC+k2]; break; 
    }      
//printf("N12=%d   AUX=%p\n",N12,AUX);    

                    
    for(nsub22=1; nsub22<= CI->nprc;nsub22++)
    { double smin;
      double a=0;

      int z4[4];
      for(i=0;i<4;i++) pname[i]=CI->pinf(nsub22,i+1,pmass+i,pdg+i);
      if(isFeeble(pname[2]) || isFeeble(pname[3])) continue;
      if(pmass[0]==0) continue; // for the case of absence of process
      for(i=0;i<4;i++) z4[i]=Z4ch(pname[i]);
      smin=pmass[2]+pmass[3];      
      cc23=NULL;
      
      if(N12==0 && (VZdecay||VWdecay))
      {  int l,l_,nVV;        
         if(!AUX[nsub22].virt )  for(l=2;l<4;l++) if(pdg[l]==21 ||pdg[l]==22) { AUX[nsub22].virt=-1; break;}
         
         if(!AUX[nsub22].virt)
         {  int vd[4]={0,0,0,0};
            int c_a =  (pmass[0]>Mcdm) || (pmass[1]>Mcdm);

            if(c_a){ for(l=2;l<4;l++) if((pdg[l]==23 && VZdecay>1)   || (abs(pdg[l])==24 && VWdecay>1)) vd[l]=1;} 
            else    for(l=2;l<4;l++)
            { 
            
              if((pdg[l]==23 && VZdecay)     || (abs(pdg[l])==24 && VWdecay)) vd[l]=1;
            } 

            for(l=2;l<4;l++) if(vd[l]) break; 
            if(l<4)
            {  l_=5-l; 
               if(vd[l_])
               { nVV=2;
                 if(pmass[l_]>pmass[l]) { l=l_; l_=5-l;}
               } else nVV=1; 
               AUX[nsub22].virt=l;  
               AUX[nsub22].w[l-2]=pWidth(pname[l],NULL);
               if(abs(pdg[l_])>16 && pmass[l_]> 2) AUX[nsub22].w[l_-2]=pWidth(pname[l_],NULL);
               if(AUX[nsub22].w[l_-2] < 0.1) AUX[nsub22].w[l_-2]=0;
            } else  AUX[nsub22].virt=-1;
         }        
         if(AUX[nsub22].virt>0)
         {  l=AUX[nsub22].virt;
            l_=5-l; 
            if(pmass[0]+pmass[1] < smin+ 4*AUX[nsub22].w[l-2]   && pmass[l_]< MassCutOut)
            { 
              if(AUX[nsub22].cc23) cc23=AUX[nsub22].cc23; else
              {  double  brV1,wV1;
                 int i3W;
                 AUX[nsub22].cc23=xVtoxll(2,2,pname,pdg, l, &wV1, &brV1);
                 if(pdg[l]==pdg[l_]) brV1*=2;
                 AUX[nsub22].br=brV1;
                 cc23=AUX[nsub22].cc23; 
                 if(cc23)
                 {  double Pcm0,PcmMax; 
                    *(cc23->interface->BWrange)=10000; 
                    *(cc23->interface->gswidth)=0;
                    for(i3W=2;i3W<5;i3W++) if(strcmp(cc23->interface->pinf(1,i3W+1,NULL,NULL),pname[l_])==0) break; 
                    AUX[nsub22].i3=i3W;        
                    PcmMax=decayPcm(pmass[2]+pmass[3]+10*AUX[nsub22].w[l-2], pmass[0],pmass[1]);
                    if( pmass[0]+pmass[1]>=pmass[l_]) Pcm0=(pmass[0]+pmass[1])*0.001; else 
                    Pcm0=decayPcm(pmass[l_]+1.E-3, pmass[0],pmass[1]);
                    for(double csTest=sigma23(Pcm0); !isfinite(csTest) || csTest==0; Pcm0*=1.01);  
                    buildInterpolation(sigma23, Pcm0,PcmMax, 0.01,1E-4,&(AUX[nsub22].nTab), &(AUX[nsub22].pcmTab), &(AUX[nsub22].csTab));
//#define TEST23
#ifdef TEST23
    char proc[100]; 
    sprintf(proc,"Log(sigma23) for %s,%s -> %s,%s; TabDim=%d", pname[0], pname[1],pname[2],pname[3],AUX[nsub22].nTab);  
    polintArgStr Arg;
    Arg.dim=AUX[nsub22].nTab; Arg.x=AUX[nsub22].pcmTab; Arg.y=AUX[nsub22].csTab; 
    displayPlot(proc,"Pcm[GeV]", Pcm0,PcmMax,0,3,"orig.",0,sigma23,NULL,"interpol.",0,polint_arg,&Arg,"22",0, Logvcs22, code);    
#endif  


                 }
              }    
              if(cc23) smin=pmass[l_];
            } 
         }
      }

      if(smin>=MassCutOut) continue; 
      if(cc23==NULL) 
      {
        if( (pmass[2]>Mlarge && pmass[3]<Msmall)
          ||(pmass[3]>Mlarge && pmass[2]<Msmall))
             { *(CI->twidth)=1; *(CI->gtwidth)=1;}
        else { *(CI->twidth)=0; *(CI->gtwidth)=0;}
      *(CI->gswidth)=0;
      }
                                   
      if(pmass[0]+pmass[1]> smin) smin=pmass[0]+pmass[1];
      v_max=m2v(smin);
      if(v_max <= v_min) continue;
repeat:
      neg_cs_flag=0;
      T_=T; 

      if(Fast_==0)
      { int err; a=simpson(v_integrand,v_min, v_max ,eps,&err);
        if(err) { do_err=do_err|err; printf("error in simpson omega.c line 2093\n");}
      }   else
      {
          int isPole=0;
          char * s;
          int m,w,n;
          double mass,width;
          double da;
          a=0;

          for(n=1;(s=code->interface->den_info(nsub22,n,&m,&w,NULL));n++)
          if(m && w && strcmp(s,"\1\2")==0 )
          { mass=fabs(code->interface->va[m]);
            width=code->interface->va[w];
            if(mass<MassCutOut && mass+8*width > pmass[0]+pmass[1]
                            && mass+8*width > smin)
            { if((pmass[0]!=M1 || pmass[1]!=M2)&&(pmass[0]!=M2 || pmass[1]!=M1))
              { double ms=pmass[0]+pmass[1];
                double md=pmass[0]-pmass[1];
                double Pcm=sqrt((mass-ms)*(mass+ms)*(mass-md)*(mass+md))/(2*mass);
                mass=sqrt(M1*M1+Pcm*Pcm)+sqrt(M2*M2+Pcm*Pcm);
              }
              if(Fast_==1) vgrid1=makeVGrid(mass,width); else vgrid1=makeVGrid2(mass,width);
              if(isPole) vgrid=crossVGrids(&vgrid,&vgrid1); else vgrid=vgrid1;
              isPole++;
            }
          }
          if(cc23)
          {  double mass,width;
             mass=pmass[2]+pmass[3];
             width= AUX[nsub22].w[0]+ AUX[nsub22].w[1];
             if(mass-width>M1+M2)
             { if(Fast_==1)  vgrid1=makeVGrid(mass,width); else  vgrid1=makeVGrid2(mass,width);
               if(isPole) vgrid=crossVGrids(&vgrid,&vgrid1); else vgrid=vgrid1;
               isPole++;                  
             }                  
          } 
          if(isPole==0)
          {  vgrid.n=1;
             vgrid.v[0]=v_min;
             vgrid.v[1]=v_max;
             vgrid.pow[0]=5;
          }
//printVGrid(vgrid);          
          for(i=0;i<vgrid.n;i++)
          {     
             if(Fast_<=0)
             { int err; 
               da=simpson(v_integrand,vgrid.v[i] ,vgrid.v[i+1],eps,&err);
               if(err) {do_err=do_err|err; printf("error in simpson omega.c line 2142\n");}
             }
             else         da=gauss(v_integrand,vgrid.v[i] ,vgrid.v[i+1],vgrid.pow[i]); 
             a+=da;             
          }
      }

      if(neg_cs_flag && *(CI->gswidth)==0)
      { *(CI->gswidth)=1;
         goto  repeat;
      }   
      

      { double br[2][5]={0}; 
//br[i][0]-SM, br[i][1]- cdm1, br[i][2]-cmd2, bri[i][3]-2*cdm1, br[i][4]-2*cdm2      
        int l1,l2,kk;   
        for(kk=2;kk<4;kk++)
        { 
          br[kk-2][z4[kk]]=1;
          if(pmass[kk]> 2*Mcdm && z4[kk]!=1 && abs(pdg[kk])!=24 && pdg[kk]!=23 ) 
          { txtList LL;
            pWidth(pname[kk],&LL);
            for(;LL;LL=LL->next)
            { double b;
              char proc[40];
              sscanf(LL->txt,"%lf %[^\n]",&b,proc);
              if(z4[kk]==0)
              {  if(strstr(proc,"~~"))     { br[kk-2][0]-=b; br[kk-2][4]+=b;}
                 else if(strstr(proc,"~")) { br[kk-2][0]-=b; br[kk-2][3]+=b;}
              } else  if(strstr(proc,"~~")==0 && strstr(proc,"~")){ br[kk-2][2]-=b; br[kk-2][3]+=b;}

            }  
          }
        }

        for(l1=0;l1<5;l1++) for(l2=0;l2<5;l2++)
        { double b=br[0][l1]*br[1][l2];
          int NN=0;
          if(b) 
          { if(l1<3 && l2<3) NN=z4[0]+3*z4[1]+9*l1+27*l2;
            else if((l1==0 && l2==3) || (l2==0 && l1==3) ) NN=z4[0]+3*z4[1]+9*1+27*1;
            else if((l1==0 && l2==4) || (l2==0 && l1==4) ) NN=z4[0]+3*z4[1]+9*2+27*2;
            else if((l1==1 && l2==3) || (l2==1 && l1==3) ) NN=0;    //  111               
            else if((l1==1 && l2==4) || (l2==1 && l1==4) ) NN=0;    //  122
            else if((l1==2 && l2==3) || (l2==2 && l1==3) ) NN=0;    //  211 
            else if((l1==2 && l2==4) || (l2==2 && l1==4) ) NN=0;    //  222 
            else if((l1==3 && l2==3)                     ) NN=0;    //  1111
            else if((l1==3 && l2==4) || (l2==3 && l1==4) ) NN=0;    //  1122
            else if((l1==4 && l2==4)                     ) NN=0;    //  2222 
              
            a*=factor*b;
            switch(NN)
            { case 1+3:       vs1100+=a; break;
              case 2+6:       vs2200+=a; break;
              case 1+3+9: 
              case 1+3+27:    vs1110+=a; break;
              case 2+6+18:
              case 2+6+54:    vs2220+=a; break;
              case 1+3+18:
              case 1+3+54:    vs1120+=a; break;
              case 1+6+9 :
              case 1+6+27:
              case 2+3+9 :
              case 2+3+27:    vs1210+=a; break;
              case 1+3+18+54: vs1122+=a; break;
              case 2+6+9+27 : vs2211+=a; break;
//====================              
              case 1+3+ 9+54:
              case 1+3+18+27: vs1112+=a; break;
              case 1+6+18+54:
              case 2+3+18+54: vs1222+=a; break;
              case 1+6+18:
              case 2+3+18:
              case 1+6+54:
              case 2+3+54:    vs1220+=a; break;
              case 2+6+9 :
              case 2+6+27:    vs2210+=a; break;
              case 2+6+9+54:
              case 2+6+18+27: vs2221+=a; break;
              case 1+6+9+27:   
              case 2+3+9+27:  vs1211+=a; break;
              
  //          default:    printf("unexpected type of process: %d %d => %d %d\n",z4[0],z4[1],z4[2],z4[3]);
            }   
          }                  
        }
      }
    }
  }
  }
  }
  { double g1=1,g2=1;
    if(CDM1) g1=geff1_(T);
    if(CDM2) g2=geff2_(T);
    
    *vs1100_=vs1100/(g1*g1);
    *vs1120_=vs1120/(g1*g1);
    *vs1210_=vs1210/(g1*g2);
    *vs2200_=vs2200/(g2*g2);
    *vs1110_=vs1110/(g1*g1);
    *vs2220_=vs2220/(g2*g2);
    *vs2210_=vs2210/(g2*g2);
//==================    

    *vs1220_=vs1220/(g1*g2);
    *vs2210_=vs2210/(g2*g2);

    if(CDM1&&CDM2) 
    { double C=exp(-fabs(Mcdm1/T -Mcdm2/T));
      if(Mcdm1>Mcdm2)
      { 
        *vs1122_=vs1122/(g1*g1);
        *vs2211_=C*C*vs1122/(g2*g2);
        *vs1112_=vs1112/(g1*g1);
        *vs1211_=C*vs1112/(g1*g2);    
        *vs1222_=vs1222/(g1*g2);
        *vs2221_=C*vs1222/(g2*g2);
      } else 
      {
        *vs1122_=C*C*vs2211/(g1*g1);
        *vs2211_=vs2211/(g2*g2);
        *vs1112_=C*vs1211/(g1*g1);
        *vs1211_=vs1211/(g1*g2);    
        *vs1222_=C*vs2221/(g1*g2);
        *vs2221_=vs2221/(g2*g2);
      }                  
    } else *vs1122_=*vs2211_=*vs1112_=*vs1211_=*vs1222_=*vs2221_=0;
  }
  WIDTH_FOR_OMEGA=0;
  
  return 0;
}


double Yeq1(double T)
{  double heff,s;
   s=2*M_PI*M_PI*T*T*T*hEff(T)/45;
   return  (T/(2*M_PI*M_PI*s))*geff1_(T)*exp(-Mcdm1/T);
}


double Yeq2(double T)
{  double s;
   s=2*M_PI*M_PI*T*T*T*hEff(T)/45;
   return  (T/(2*M_PI*M_PI*s))*geff2_(T)*exp(-Mcdm2/T); 
}


static double Y1SQ_Y2(double T)
{ 
  double s,X1,X2,g1_,g2_,res;
  X1=Mcdm1/T;
  X2=Mcdm2/T;
  if(X2-2*X1>500) return 0;
  s=2*M_PI*M_PI*T*T*T*hEff(T)/45;
  
  g1_=geff1_(T);
  g2_=geff2_(T);
      
  res = g1_/g2_*exp(X2-2*X1)*T*g1_/(2*M_PI*M_PI*s);
//  if(!isfinite(res)) { printf("T=%E  g1_=%E g2_=%E X1=%E X2=%E\n", T,g1_,g2_,X1,X2); exit(0);}
  return res;
}                          

static double Y2SQ_Y1(double T)
{ 
  double s,X1,X2,g1_,g2_,res;
  X1=Mcdm1/T;
  X2=Mcdm2/T;
  if(X1-2*X2>500) return 0;
  s=2*M_PI*M_PI*T*T*T*hEff(T)/45;
  
  g1_=geff1_(T);
  g2_=geff2_(T);
      
//  res = g1_/g2_*exp(X1-2*X2)*T*g1_/(2*M_PI*M_PI*s);
  res = g2_/g1_*exp(X1-2*X2)*T*g2_/(2*M_PI*M_PI*s);
//  if(!isfinite(res)) { printf("T=%E  g1_=%E g2_=%E X1=%E X2=%E\n", T,g1_,g2_,X1,X2); exit(0);}
  return res;
}                          




static int DMEQ0(double T, 
  double vs1100,double vs2200,
  double vs1110,double vs1120,double vs1210,double vs1220,double vs2210,double vs2220, 
  double vs1112,double vs1122,double vs1222,double vs1211,double vs2211,double vs2221,
  double *C, double *L, double *Q)
{
  double y1,y2,y1_, y2_,coef;
  int i;
  
  if(CDM1)  y1=Yeq1(T); else y1=0;
  if(CDM2)  y2=Yeq2(T); else y2=0;
  if(CDM1&&CDM2) { y1_=Y1SQ_Y2(T); y2_=Y2SQ_Y1(T);}  else { y1_=0; y2_=0;}

//vs1120*=0.0;

/* 
  if(Mcdm1 >Mcdm2)z1122= vs1122*(y1*y1-y2*y2*y1_y2_Q(T)); else  z1122=-vs2211*(y2*y2-y1*y1*y2_y1_Q(T));
  dYdT[0]=vs1100*(y1*y1-y1_*y1_)+     vs1120*(y1*y1-y2*y1q_y2) + z1122;
  dYdT[1]=vs2200*(y2*y2-y2_*y2_)- 0.5*vs1120*(y1*y1-y2*y1q_y2) - z1122 + 0.5*vs1210*y1*(y2-y2_) ;
*/


  L[0]= y1*(2*(vs1100+vs1120+vs1122)+0.5*(vs1110+vs1112))+ 0.5*y2*(vs1220+vs1222)+0.5*y2_*vs2210; //A1_1
  if(!isfinite(L[0])) { printf("y1=%e y2=%e y2_=%e\n", y1,y2,y2_);}
  L[3]= y2*(2*(vs2200+vs2210+vs2211)+0.5*(vs2220+vs2221))+ 0.5*y1*(vs1210+vs1211)+0.5*y1_*vs1120; //A2_2

//printf("y1=%E y2=%E y1_=%E y2_=%E\n", y1,y2,y1_,y2_);

  L[1]= -0.5*y1*vs1222 - y2*2*vs2211-0.5*y1*vs1211 -y1_*vs1120 -y2*vs2210;                //A1_2
//printf(" y1*y1*vs1122/(y2*y2*vs2211)=%E  y2*2*vs2211=%E  y1_*vs1120=%E\n",y1*y1*vs1122/(y2*y2*vs2211)  ,    y2*2*vs2211,y1_*vs1120);
  L[2]= -0.5*y2*vs1211 - y1*2*vs1122-0.5*y2*vs1222 -y2_*vs2210 -y1*vs1120;                //A2_1

//printf("L[2]=%e = %e + %e  vs2210=%E   \n", L[2], - y1_*2*vs1122, -vs2210*y2q_y1, vs2210);
 
  Q[0] =  vs1100+vs1120+vs1122+0.5*(vs1110+vs1112); //Q1_11
  Q[5] =  vs2200+vs2210+vs2211+0.5*(vs2220+vs2221); //Q2_22

  Q[1] =  0.5*(vs1222+vs1220-vs1211);  //Q1_12
  Q[4] =  0.5*(vs1211+vs1210-vs1222);  //Q2_12

  Q[2] = -vs2211 -0.5*vs2221 -0.5*vs2210;          //Q1_22
  Q[3] = -vs1122 -0.5*vs1112 -0.5*vs1120;          //Q2_11

  coef=sqrt(M_PI/45)*hEff(T)*(1+hEffLnDiff(T)/3)/sqrt(gEff(T))*MPlanck;
    
//printf("T=%E Q[3]=%E coef=%E sqrt_gStar=%E \n", T,Q[3],coef,sqrt_gStar);  
  for(i=0;i<4;i++) L[i]*=coef;
  for(i=0;i<6;i++) Q[i]*=coef;
  
  
  if(CDM1)  C[0]=(Yeq1(T*1.01)-Yeq1(T/1.01))/(2*log(1.01)*T); else C[0]=0;
  if(CDM2)  C[1]=(Yeq2(T*1.01)-Yeq2(T/1.01))/(2*log(1.01)*T); else C[1]=0;

{ int i;
for(i=0;i<2;i++) if(!isfinite(C[i])) {printf("T=%E, C[%d]=%E\n",T,i,C[i]); exit(1);}
for(i=0;i<4;i++) if(!isfinite(L[i])) {printf("T=%E, L[%d]=%E\n",T,i,L[i]); exit(1);}
for(i=0;i<6;i++) if(!isfinite(Q[i])) {printf("T=%E, Q[%d]=%E\n",T,i,Q[i]); exit(1);}
}
  return 0;
}



static int Ntab2=0;
static double*Ttab2=NULL;


static double *vs1100T = NULL;
static double *vs1120T = NULL;
static double *vs1122T = NULL;
static double *vs1210T = NULL;
static double *vs2200T = NULL;
static double *vs2211T = NULL;

static double *vs1110T = NULL;
static double *vs2220T = NULL;
static double *vs1112T = NULL;
static double *vs1222T = NULL;
static double *vs1220T = NULL;
static double *vs2210T = NULL;
static double *vs2221T = NULL;
static double *vs1211T = NULL;
static int d1100, d1120, d1122, d1210, d2200, d2211, d1110, d2220, d1112, d1222, d1220, d2210, d2221, d1211;

static double *Y1T=NULL;
static double *Y2T=NULL;


static double vsInterpolation( double T,  int d, double *vsTab, int p_aux,   int d_aux, double *vsTab_aux)
{
   if(Ntab2<=0) return NAN;
   if(T<0.9*Ttab2[Ntab2-1] || T> 1.1*Ttab2[0]) return NAN;
   if(d<2 && d_aux<2) return 0;
   if(vsTab_aux)
   { 
//      printf("T=%e p_aux=%d   geff_2(T)=%E geff_1(T)=%E Mcdm1=%e Mcem2=%e C=%e   \n", T, p_aux,   geff2_(T),  geff1_(T),Mcdm1,Mcdm2, pow( geff2_(T)/geff1_(T),p_aux)*exp(-p_aux*(Mcdm2-Mcdm1)/T));
//      exit(0);
      if(p_aux>0 && Mcdm1<Mcdm2) return vsInterpolation( T,  d_aux, vsTab_aux, 0, 0, NULL)* pow( geff2_(T)/geff1_(T),p_aux)*exp(-p_aux*(Mcdm2-Mcdm1)/T);
      if(p_aux<0 && Mcdm2<Mcdm1) return vsInterpolation( T,  d_aux, vsTab_aux, 0, 0, NULL)* pow( geff2_(T)/geff1_(T),p_aux)*exp(-p_aux*(Mcdm2-Mcdm1)/T);
   }
   if(d<2) return 0;
   if(d==Ntab2 || T>Ttab2[d-2]) return polint3(T,Ntab2,Ttab2, vsTab);
   if(d<=1 || T<=Ttab2[d]) return 0;
   if(T<Ttab2[d-1]) { double alphaT=(T-Ttab2[d])/(Ttab2[d-1]-Ttab2[d]);      return vsTab[d-1]*alphaT*alphaT; }
   if(T<=Ttab2[d-2]) { double alphaT=(T-Ttab2[d-1])/(Ttab2[d-2]-Ttab2[d-1]);  return pow(vsTab[d-1],1-alphaT )*pow(vsTab[d-2],alphaT );} 
   printf("T=%E vsInterpolation out of rules d=%d Ntab2=%d\n",T,d,Ntab2);  
   return 0;   
}


char*ExcludedFor2DM=NULL;


static int DMEQ(int useTab, double T, double *C, double *L, double *Q)
{

  double  vs1100,vs2200,vs1110,vs2220,vs1120,vs1122,vs1210,vs2211,vs1112,vs1222,vs1220,vs2210,vs2221,vs1211;
  
  if(useTab)
  {
    vs1100= vsInterpolation( T,d1100, vs1100T,  0,  0, NULL);
    vs2200= vsInterpolation( T,d2200, vs2200T,  0,  0, NULL);
    vs1120= vsInterpolation( T,d1120, vs1120T,  0,  0, NULL);
    vs1210= vsInterpolation( T,d1210, vs1210T,  0,  0, NULL);
    vs1110= vsInterpolation( T,d1110, vs1110T,  0,  0, NULL);
    vs2220= vsInterpolation( T,d2220, vs2220T,  0,  0, NULL);
    vs1220= vsInterpolation( T,d1220, vs1220T,  0,  0, NULL);
    vs2210= vsInterpolation( T,d2210, vs2210T,  0,  0, NULL);
 
    vs1122= vsInterpolation( T,d1122, vs1122T,  2, d2211, vs2211T );
    vs2211= vsInterpolation( T,d2211, vs2211T, -2, d1122, vs1122T );
    vs1112= vsInterpolation( T,d1112, vs1112T,  1, d1211, vs1211T );
    vs1211= vsInterpolation( T,d1211, vs1211T, -1, d1112, vs1112T );
    vs1222= vsInterpolation( T,d1222, vs1222T,  1, d2221, vs2221T );
    vs2221= vsInterpolation( T,d2221, vs2221T, -1, d1222, vs1222T );
  } else
  {
    aRate4(T,&vs1100,&vs2200,
    &vs1110,&vs1120,&vs1210,&vs1220,&vs2210,&vs2220,
    &vs1112,&vs1122,&vs1222,&vs1211,&vs2211,&vs2221);
  }  
            
  if(ExcludedFor2DM)
  {
    if(strstr(ExcludedFor2DM,"1100")) vs1100=0;
    if(strstr(ExcludedFor2DM,"2200")) vs2200=0;
    if(strstr(ExcludedFor2DM,"1110")) vs1110=0;
    if(strstr(ExcludedFor2DM,"2220")) vs2220=0;
    if(strstr(ExcludedFor2DM,"1210")) vs1210=0;
    if(strstr(ExcludedFor2DM,"1220")) vs1220=0;
    if(strstr(ExcludedFor2DM,"2210")) vs2210=0;
    if(strstr(ExcludedFor2DM,"1120")) vs1120=0; 

    if(strstr(ExcludedFor2DM,"1122")||strstr(ExcludedFor2DM,"2211")) {vs1122=0;vs2211=0;} 
    if(strstr(ExcludedFor2DM,"1112")||strstr(ExcludedFor2DM,"1211")) {vs1112=0;vs1211=0;}
    if(strstr(ExcludedFor2DM,"1222")||strstr(ExcludedFor2DM,"2212")) {vs1222=0;vs2221=0;}
  }
                
  return  DMEQ0(T, vs1100,vs2200,
    vs1110,vs1120,vs1210,vs1220,vs2210,vs2220,
    vs1112,vs1122,vs1222,vs1211,vs2211,vs2221, 
         C,L,Q);
}

static int dYstart(double T, double * dy, double * Lmin,double *Lmax)
{
  double C[2],L[4],Q[6],D;
  
  DMEQ(0, T, C, L, Q);

  if(!CDM2) { dy[0]=C[0]/L[0]; dy[1]=0; return 0;}
  if(!CDM1) { dy[0]=0; dy[1]=C[1]/L[3]; return 0;}
    
  D=L[0]*L[3]-L[1]*L[2];
  if(dy)
  { 
    dy[0]= ( C[0]*L[3]-C[1]*L[1])/D;
    dy[1]= (-C[0]*L[2]+C[1]*L[0])/D;
//    printf("T=%E  C[0]=%e, C[1]=%e L[0-4]=%E %E %E %E    dy= %E %E\n", T, C[0],C[1],L[0],L[1],L[2],L[3], dy[0],dy[1]);
  }  
{ 
  double s,d;
  s=0.5*(L[0]+L[3]);
  d= s*s+L[1]*L[2]-L[0]*L[3];
  if(Lmin) *Lmin=s-sqrt(d);
  if(Lmax) *Lmax=s+sqrt(d);
}  
  return 0;
}


static void stiffDerives(double T, double*Y,double*f,double h,double*dfdx,double*dfdy)
{

 double C[2], L[4], Q[6],dy1,dy2;
 double dT=-0.001*T;
 int n=2;  
 DMEQ(1,T, C, L, Q);
        
 dy1=Y[0];
 dy2=Y[1];
              
  f[0]= -C[0] + L[0]*dy1 + L[1]*dy2 + Q[0]*dy1*dy1 + Q[1]*dy1*dy2 + Q[2]*dy2*dy2;
  f[1]= -C[1] + L[2]*dy1 + L[3]*dy2 + Q[3]*dy1*dy1 + Q[4]*dy1*dy2 + Q[5]*dy2*dy2;
 
  if(dfdy)
  { dfdy[0*n+0]= L[0]+2*Q[0]*dy1+Q[1]*dy2;               
    dfdy[0*n+1]= L[1]+2*Q[2]*dy2+Q[1]*dy1;
    dfdy[1*n+0]= L[2]+2*Q[3]*dy1+Q[4]*dy2;
    dfdy[1*n+1]= L[3]+2*Q[5]*dy2+Q[4]*dy1;
  }

  if(dfdx)
  { DMEQ(1,T+dT, C, L, Q);
    dfdx[0]= -C[0] + L[0]*dy1 + L[1]*dy2 + Q[0]*dy1*dy1 + Q[1]*dy1*dy2 + Q[2]*dy2*dy2;
    dfdx[1]= -C[1] + L[2]*dy1 + L[3]*dy2 + Q[3]*dy1*dy1 + Q[4]*dy1*dy2 + Q[5]*dy2*dy2;
    dfdx[0]-=f[0]; dfdx[0]/=dT;
    dfdx[1]-=f[1]; dfdx[1]/=dT;
  }
}


static void TabDmEq(double step)
{
  int i,N;
  double T;

//printf("TabDmEq:  Tstart=%E Tend=%E\n",Tstart,Tend);  
  N=log(Tstart/Tend)/log(step)+2;
  if(N!=Ntab2)
  {                 
    Ttab2    = realloc(Ttab2,   N*sizeof(double));
    vs1100T = realloc(vs1100T,N*sizeof(double));
    vs1120T = realloc(vs1120T,N*sizeof(double));
    vs1122T = realloc(vs1122T,N*sizeof(double));
    vs1210T = realloc(vs1210T,N*sizeof(double));
    vs2200T = realloc(vs2200T,N*sizeof(double));
    vs1110T = realloc(vs1110T,N*sizeof(double));
    vs2220T = realloc(vs2220T,N*sizeof(double));
    vs2211T = realloc(vs2211T,N*sizeof(double));
    vs1112T = realloc(vs1112T,N*sizeof(double)); 
    vs1222T = realloc(vs1222T,N*sizeof(double));
    vs1220T = realloc(vs1220T,N*sizeof(double));
    vs2210T = realloc(vs2210T,N*sizeof(double));
    vs2221T = realloc(vs2221T,N*sizeof(double));
    vs1211T = realloc(vs1211T,N*sizeof(double));
    Y1T     = realloc(Y1T,    N*sizeof(double));
    Y2T     = realloc(Y2T,    N*sizeof(double)); 
//    TCoeff  = realloc(TCoeff,N*sizeof(double));
    Ntab2=N;
  }  

  double * vs[14]={ vs1100T,vs1120T,vs1122T,vs1210T,vs2200T,vs1110T,vs2220T,vs2211T,vs1112T,vs1222T,vs1220T,vs2210T,vs2221T,vs1211T};
  int *     d[14]={ &d1100, &d1120, &d1122, &d1210, &d2200, &d1110, &d2220, &d2211, &d1112, &d1222, &d1220, &d2210, &d2221, &d1211};
  
  for(T=Tstart,i=0;i<N;i++)
  { 
     Ttab2[i]=T;
     Y1T[i]=NAN;
     Y2T[i]=NAN;
     aRate4(T, vs1100T+i,vs2200T+i,
        vs1110T+i,vs1120T+i,vs1210T+i,vs1220T+i,vs2210T+i,vs2220T+i,
        vs1112T+i,vs1122T+i,vs1222T+i,vs1211T+i,vs2211T+i,vs2221T+i);
     T/=step;     
  }

  for(int k=0;k<14;k++)
  {  
     for(i=N-1; i>=0 && vs[k][i]==0;i--) continue;
     *d[k]=i+1;
  }

  if(!CDM1) { d1100=0; d1120=0; d1122=0; d1210=0;          d2211=0;  d1110=0;          d1112=0; d1222=0; d1220=0; d2210=0; d2221=0; d1211=0; }
  if(!CDM2) {          d1120=0; d1122=0; d1210=0; d2200=0; d2211=0;           d2220=0; d1112=0; d1222=0; d1220=0; d2210=0; d2221=0; d1211=0; }

  if(ExcludedFor2DM)
  {
    if(strstr(ExcludedFor2DM,"1100")) {d1100=0; printf("vs1100 ignored\n");}
    if(strstr(ExcludedFor2DM,"2200")) {d2200=0; printf("vs2200 ignored\n");}
    if(strstr(ExcludedFor2DM,"1110")) {d1110=0; printf("vs1110 ignored\n");}
    if(strstr(ExcludedFor2DM,"2220")) {d2220=0; printf("vs2220 ignored\n");}
    if(strstr(ExcludedFor2DM,"1210")) {d1210=0; printf("vs1210 ignored\n");} 
    if(strstr(ExcludedFor2DM,"1220")) {d1220=0; printf("vs1220 ignored\n");}
    if(strstr(ExcludedFor2DM,"2210")) {d2210=0; printf("vs2210 ignored\n");}
    if(strstr(ExcludedFor2DM,"1120")) {d1120=0; printf("vs1120 ignored\n");}
    if(strstr(ExcludedFor2DM,"1122")||strstr(ExcludedFor2DM,"2211")) {d1122=0; d2211=0;  printf("vs1122&vs2211 ignored\n");} 
    if(strstr(ExcludedFor2DM,"1112")||strstr(ExcludedFor2DM,"1211")) {d1112=0; d1211=0;  printf("vs1112&vs1211 ignored\n");}
    if(strstr(ExcludedFor2DM,"1222")||strstr(ExcludedFor2DM,"2212")) {d1222=0; d2221=0;  printf("vs1222&vs2221 ignored\n");}
  }  
}



static void TderivZ4tab2(double T, double *Y, double *dYdT)
{
  double C[2], L[4], Q[6],dy1,dy2;

  DMEQ(1,T,C,L,Q);  
  dy2=Y[0];
  dy1=0;
   
  dYdT[0]= -C[1] + L[2]*dy1 + L[3]*dy2 + Q[3]*dy1*dy1 + Q[4]*dy1*dy2 + Q[5]*dy2*dy2;
}

static void TderivZ4tab1(double T, double *Y, double *dYdT)
{
  double C[2], L[4], Q[6],dy1,dy2;

  DMEQ(1,T,C,L,Q);
  
  dy1=Y[0];
  dy2=0;

  dYdT[0]= -C[0] + L[0]*dy1 + L[1]*dy2 + Q[0]*dy1*dy1 + Q[1]*dy1*dy2 + Q[2]*dy2*dy2;
}



double darkOmega2TR( double TR, double Y1R,double Y2R, double fast, double Beps0)
{
  Fast_=fast;
  Beps=Beps0;
  double Y[2],YY[2],T;
  double step=1.1;
  double ips=0.01,ips_=0.005;  
  int i,err,N; 

  fillCDMmem();

  if(CDM1==NULL && CDM2==NULL)
  { 
     restoreCDMmem();
     printf(" There are no  Dark Matter particles\n");
     return 0; 
  }
  
  dmAsymm=0;
  
  Tstart=TR;
  TabDmEq(step);
  Y[0]=-Yeq1(TR)+Y1R;
  Y[1]=-Yeq2(TR)+Y2R;


  Y1T[0]=Y[0];
  Y2T[0]=Y[1]; 

 
  if(!CDM1)
  {
    for(i=1,err=0,T=Tstart; T> Tend && !err;i++ )
    { double T2=T/step;
      if(T2<Tend) T2=Tend;
      err=odeint(Y+1,1 , T ,T2 , 1.E-3, (T-T2) , &TderivZ4tab2);
      if(err) { printf(" error in odeint\n");  return  -1;}
      Y1T[i]=Y1T[0];      
      Y2T[i]=Y[1];
      T=T2;
    }
    fracCDM2=1;
    restoreCDMmem();
    return (Y[1]+Yeq2(Tend))*2.742E8*Mcdm2; 
  } else if(!CDM2) 
  { 
  
    for(i=1,err=0,T=Tstart; T> Tend && !err;i++ )
    { double T2=T/step;
      if(T2<Tend) T2=Tend;
      err=odeint(Y,1 , T ,T2 , 1.E-3, (T-T2) , &TderivZ4tab1);
      if(err) { printf(" error in odeint\n");  return  -1;}
      Y1T[i]=Y[0];      
      Y2T[i]=Y2T[0];
      T=T2;
    }
    fracCDM2=0;
    restoreCDMmem();
    return (Y[0]+Yeq1(Tend))*2.742E8*Mcdm1;
  }  else  
  { double h=0.01*Tstart*(1-1/step);
    for(i=1,err=0,T=Tstart; T> Tend && !err;i++ )
    { double T2=T/step;
      double Yscal[2]={1,1};
      if(T2<Tend) T2=Tend;
      
      err=stiff(i==1,T,T2,2,Y, Yscal,1.E-3, &h, stiffDerives);
             
      if(err) { printf(" error in stiff at T=[%.2E, %.2E]\n",T2,T);  return  -1;}
      Y1T[i]=Y[0];      
      Y2T[i]=Y[1];      
      T=T2;
//printf("i=%d T=%E %E %e\n", i,T,   Y1T[i]+Yeq1(T),Y2T[i]+Yeq2(T));    

    }
    
    restoreCDMmem();    

    fracCDM2=(Y[1]+Yeq2(Tend))*Mcdm2/( (Y[0]+Yeq1(Tend))*Mcdm1 +(Y[1]+Yeq2(Tend))*Mcdm2);
    return  (Y[0]+Yeq1(Tend)) *2.742E8*Mcdm1 + (Y[1]+Yeq2(Tend))*2.742E8*Mcdm2;
  }
}




double darkOmega2( double fast, double Beps0)
{
  Fast_=fast;
  Beps=Beps0;
  double Y[2],YY[2],T;
  double Lmin,Lmax;
  double step=1.1;
  double ips=0.01,ips_=0.005;  
  int i,err,N; 
  
  fillCDMmem();

  if(CDM1==NULL && CDM2==NULL)
  { 
     restoreCDMmem();
     printf(" There are no  Dark Matter particles\n");
     return 0; 
  }
  
  dmAsymm=0;
  if(!CDM1) Tstart= Mcdm2/20; else if(!CDM2) Tstart= Mcdm1/20;
  else { if(Mcdm1>Mcdm2) Tstart=Mcdm1/20; else Tstart= Mcdm2/20;} 
  Ntab2=0;
  dYstart(Tstart,Y,&Lmin,&Lmax);
  
  if(!CDM1)
  {   
     while( !isfinite(Y[1])  ||  fabs(Y[1])>0.01 *Yeq2(Tstart)) 
     { Tstart*=1.05; 
       if(Tstart>Mcdm2) { printf(" darkOmega2 can not find a starting point with DM in thermal equilibrium with SM\n");   restoreCDMmem();  return 0;}
       dYstart(Tstart,Y,&Lmin,&Lmax);
     }
     while(  isfinite(Y[1])  &&  fabs(Y[1])<0.005*Yeq2(Tstart)) { Tstart/=1.05; dYstart(Tstart,Y,&Lmin,&Lmax);}
  } else if(!CDM2)
  {  
     while( !isfinite(Y[0])  ||  fabs(Y[0])>0.01 *Yeq1(Tstart)) 
     { Tstart*=1.05; 
       if(Tstart>Mcdm1) { printf(" darkOmega2 can not find a starting point with DM in thermal equilibrium with SM\n");   restoreCDMmem();  return 0;} 
       dYstart(Tstart,Y,&Lmin,&Lmax);
     }     
     while(  isfinite(Y[0])  &&  fabs(Y[0])<0.005*Yeq1(Tstart)) {  Tstart/=1.05;  dYstart(Tstart,Y,&Lmin,&Lmax); } 

  } else
  {  
     while( Lmin<100/Tstart  ) 
     { Tstart*=1.05; 
       if(Tstart>Mcdm1 && Tstart>Mcdm2) { printf(" darkOmega2 can not find a starting point with DM in thermal equilibrium with SM\n");   restoreCDMmem();  return 0;}
       dYstart(Tstart,Y,&Lmin,&Lmax);
     } 
     while( Lmin>200/Tstart  ) { Tstart/=1.05; dYstart(Tstart,Y,&Lmin,&Lmax);} 
  }

  restoreCDMmem();
  return darkOmega2TR(Tstart, Y[0]+Yeq1(Tstart),Y[1]+Yeq2(Tstart), fast, Beps0);
}



//============================= Thermal equilibrium  of DM components ===============

static double sqrtSmin;

static double cWidthInt(double u)
{
   if(u==1|| u==0)  return 0;
   double z;
//   z=u*(2-u);
   z= u*u*u*(4-3*u);   
   double sqrtS=sqrtSmin-3*T_*logl(z);
   
   long double  ms = pmass[0] + pmass[1];
   long double  md = pmass[0] - pmass[2];
   double  PcmIn = sqrtl((sqrtS-ms)*(sqrtS+ms)*(sqrtS-md)*(sqrtS+md))/(2*sqrtS);
            
   
//   double PcmIn=decayPcm(sqrtS,pmass[0],pmass[1]);
   kin22(PcmIn,pmass);
   int err=0;
   double cs;
   double dcs; 
//   cs=peterson21(dSigma_dCos,-1,1,&dcs);
//int err=0;
   cs=simpson(dSigma_dCos,-0.98,0.98,1E-2,NULL);
//if(err)
//{
//   displayPlot("dSigma_dCos","cos",  -1,0.99,0,1,"",0,dSigma_dCos,NULL);
//   exit(1);
//}   
   
//   if(0.3*fabs(cs)< fabs(dcs)){ printf("PcmIn=%e => 0\n", PcmIn);    return 0;}  
//   if(0.1*fabs(cs)< fabs(dcs))   return 0;  
//   cs=simpson(dSigma_dCos,-1,1,0.5*1E-2,&err);
//   if(err) return 0;
   
   if(err) { 
     static int Ntot=0;
     printf("Pcm=%E\n", PcmIn);
     displayPlot("dSigma_dCos","cos",-1,1,0,1,"",0,dSigma_dCos,NULL);
     printf("error in dCos\n"); 
     Ntot++; if(Ntot>10)exit(0);
           }

     return z*z*6*u*u*(1-u)*PcmIn*PcmIn*sqrtS*sqrt(sqrtS)*cs*K1pol(T_/sqrtS);
//   return z*z*(1-u)*PcmIn*PcmIn*sqrtS*sqrt(sqrtS)*cs*K1pol(T_/sqrtS);
}


double collisionWidth(numout*cc, double Beps,double T)
{  passParameters(cc);
   double  width=0;
   sqme22=cc->interface->sqme;
   for( nsub22=1; nsub22<=cc->interface->nprc; nsub22++)
   {  int pnum[4];
      char*pname[4];
      for (int i=0;i<4;i++) pname[i]= cc->interface->pinf(nsub22,i+1,pmass+i,pnum+i);
      int BepsOK=1;
      for(int i=1;i<4;i++) if( pname[i][0]!='~' &&  exp(-pmass[i])/T <Beps) BepsOK=0; 
      if(!BepsOK) continue; 
//      printf(" %s %s %s %s ??\n", pname[0],pname[1],pname[2],pname[3]); 
//      int ok=1;
//      for(int i=0;i<4;i++) if( pname[i][0]!='~' && !isSMP(pnum[i])) ok=0;
//      if(!ok) continue;
//      if(pnum[1]==21 || pnum[2]==21 || pnum[3]==21 ||pnum[1]==22 || pnum[2]==22 || pnum[3]==22) continue;  

//    exclude reactions which allow decay ~X->~Y,1*x          
      txtList L;     
      int inDecay=0;
      pWidth(pname[0],&L);
      for(;L;L=L->next)
      { char*ch=ch=strstr(L->txt,"->");
        ch+=2;
        char p[4][20];
//        printf("L->txt=%s\n", L->txt);
        int n=sscanf(ch," %[^,], %[^,], %[^,], %s", p[0],p[1],p[2],p[3]);
        if(n>2) continue;     
//        int ok=1;
        if( (strcmp(p[0],pname[2])==0 && strcmp(p[1],pname[3])==0) || (strcmp(p[1],pname[2])==0 && strcmp(p[0],pname[3])==0) ) inDecay=1;   
//        for(int i=0;i<2;i++) for(int j=0;j<2;j++)
//        if(pname[2+j][0]=='~' && strcmp(p[i],pname[2+j])==0   ) inDecay=1;  
//        printf("pname[2]=%s,pname[3]=%s p[0]=%s,p[1]=%s inDecay=%d\n", pname[2],pname[3],p[0],p[1], inDecay);    

      }
      
//      if(inDecay) continue;
//printf("CONTINUE\n");
      sqrtSmin=pmass[0]+pmass[1];
      if(pmass[2]+pmass[3]> sqrtSmin) sqrtSmin=pmass[2]+pmass[3];
      T_=T;
      int ndf;
      cc->interface->pinfAux(nsub22, 2,NULL, NULL,NULL,&ndf);

/*
      if(  strcmp(pname[0],"~N2")==0 && 
           strcmp(pname[1],"d1")==0 && 
           strcmp(pname[2],"~N1")==0 && 
           strcmp(pname[3],"d2")==0)
*/           
//      if(cWidthInt(0.5)!=0) {
      
//      char txt[100]; sprintf(txt, "%s %s %s %s", pname[0],pname[1],pname[2],pname[3]);
      
//       displayPlot(txt,"u",0,1,0,1,  "",0,cWidthInt,NULL);
//        exit(0);
//      }
//       else printf("continue %s %s %s %s  \n",pname[0],pname[1],pname[2],pname[3] );
//int err=0; 
//      printf(" %s %s => %s %s\n", pname[0],pname[1],pname[2],pname[3]);
      double dWidth=3*ndf*T*exp(-(sqrtSmin-pmass[0])/T)/(M_PI*M_PI*pow(pmass[0],1.5)*K2pol(T/pmass[0]))*simpson(cWidthInt,0,1, 1E-2,NULL);
//printf("simpson(cWidthInt,0,1, 1E-2,NULL)=%e\n", simpson(cWidthInt,0,1, 1E-2,NULL));

//if(err){
//   displayPlot("cWidthInt","u",0,1,0,1,"cWidthInt",0,cWidthInt,NULL);
//       exit(0);
//       } 

      width+= dWidth;
      char txt[100];
//sprintf(txt," %s %s -> %s %s ", pname[0],pname[1],pname[2],pname[3]);
//printf(" %-30.30s   %.2E \n", txt,dWidth);
   }

   return width;
} 


//#include <quadmath.h>


double checkFO(double Beps,double T)
{
   char allSM[2000]={""};
   int nP=0;
   printf("\ncheckFO: Beps=%.1E, T=%.1E Mcdm=%.1E\n", Beps,T,Mcdm);
   int * all=malloc(sizeof(int)*Nodd);
      
   for(int i=0;i<Nodd;i++) 
   {  if(isFeeble(OddPrtcls[i].name)) continue;
      if(exp((Mcdm-pMass(OddPrtcls[i].name))/T)>Beps)  all[nP++]=i; 
   }

   for(int i=0;i<nP-1;)
   { if(pMass(OddPrtcls[all[i]].name) > pMass(OddPrtcls[all[i+1]].name)) { int m=all[i]; all[i]=all[i+1];all[i+1]=m; if(i) i--;}
     else i++;   
   } 

   printf("Odd Particles:\n");
   long double*Y=malloc(sizeof(long double)*nP);
   int * ndf=malloc(sizeof(int)*nP);
   int *ch=malloc(sizeof(int)*nP);
   for(int i=0;i<nP;i++) 
   {  ch[i]=(OddPrtcls[all[i]].name!=OddPrtcls[all[i]].aname); 
      ndf[i]= OddPrtcls[all[i]].cdim*(OddPrtcls[all[i]].spin2+1)*(ch[i]+1) ;
//      double Y_=ndf[i]*exp((Mcdm-pMass(OddPrtcls[all[i]].name))/T);
      double dmMass=pMass(OddPrtcls[all[i]].name); 
      Y[i]=T/(2*M_PI*M_PI)*ndf[i]*dmMass*dmMass*bessK2(dmMass/T)/(2*M_PI*M_PI*T*T*T*hEff(T)/45)*exp(Mcdm/T);
//Y[i]=Y_;      
      printf("name=%10.10s    pMass=%.3E  ndf=%d  Yeq=%E\n", OddPrtcls[all[i]].name, pMass(OddPrtcls[all[i]].name), ndf[i],(double)Y[i]);   
   }

//   for(int i=0;i<nP;i++) printf("%s %d %d %E \n", OddPrtcls[all[i]].name,ch[i], ndf[i], (double)Y[i]);

   long double *M=malloc(sizeof(long double)*nP*nP);
   for(int i=0;i<nP*nP;i++) M[i]=0;
   
   for(int i=1;i<nP;i++)
   {  txtList L;
      double width=pWidth(OddPrtcls[all[i]].name ,&L);
      for(int j=0;j<i;j++)
      { char pattern[20];
        sprintf(pattern,"%s,*",OddPrtcls[all[j]].name);
        double br=findBr(L, pattern);
//        printf(" %s ->  %s  %E\n", OddPrtcls[all[i]].name, OddPrtcls[all[j]].name, width*br);
        M[j*nP+i]=width*br; 
        if(!ch[j]) continue; 
        sprintf(pattern,"%s,*",OddPrtcls[all[j]].aname);
        br=findBr(L, pattern);
        M[j*nP+i]+=width*br;   
      }    
   }
   for(int i=1;i<nP;i++) for(int j=0;j<i;j++) M[i*nP+j]=M[j*nP+i]*Y[i]/Y[j];
   for(int i=0;i<nP;i++) { M[i*nP+i]=0; for(int j=0;j<nP;j++) if(i!=j) M[i*nP+i]-=M[j*nP+i]; }

   printf("\nDecay Matrix:\n");
   for(int i=0;i<nP;i++) { for(int j=0;j<nP;j++) printf(" %E ",(double)M[i*nP+j]); printf("\n");}

   
   for(int i=0;i<nP;i++) { long double s=0; for(int j=0;j<nP;j++) s+=M[j*nP+i]; printf("i=%d %E\n",i,(double)s);} printf("\n");
   for(int i=0;i<nP;i++) { long double s=0; for(int j=0;j<nP;j++) s+=M[i*nP+j]*Y[j]; printf("i=%d %E\n",i,(double)s);}

   
   double omg_eff=0;
   for(int k=1;k<=(1<<(nP-1))-1;k++)
   {  int q=2*k;
      double Ys1=0,Ys2=0,omgs=0;
      for(int i=0;i<nP;i++) for(int j=0;j<nP;j++) if(((1<<i) & k) && (!((1<<j) & k))) omgs+= M[i*nP+j]*Y[j];
      for(int i=0;i<nP;i++) if((1<<i)&k) Ys1+=Y[i]; else Ys2+=Y[i];
//      printf("k=%d omgs=%E Ys1=%E Ys2=%E ", k,omgs,Ys1,Ys2);
      omgs*=1/Ys1+1/Ys2;
//      printf(" omg_eff=%E\n",  omgs);
      if(k==1) omg_eff=omgs; else if(omgs<omg_eff) omg_eff=omgs;
   }
   
   printf("omega_eff(T)=%.2E(Decays only) Hubble=%.2E\n",omg_eff,Hubble(T));
   if(omg_eff<10*Hubble(T))
   { 
      printf("Collision widths:\n");
      for(int i=1;i<nP;i++) for(int j=0;j<i;j++) 
      {  char process[4000], lib[40], lib1[12],lib2[12];
         pname2lib(OddPrtcls[all[i]].name,lib1);
         for(int k=0;k<2;k++)
         {
            if(k==0) pname2lib(OddPrtcls[all[j]].name,lib2); 
            else  if(OddPrtcls[all[i]].name!=OddPrtcls[all[i]].aname) pname2lib(OddPrtcls[all[j]].aname,lib2);
            else  break;    
           
            sprintf(process,"%s,AllEven->%s,AllEven{",OddPrtcls[all[i]].name, OddPrtcls[all[j]].name);
            int SMcode[28]={1,2,3,4,5,6,11,12,13,14,15,16,21,24,-1,-2,-3,-4,-5,-6,-11,-12,-13,-14,-15,-16,-24,23};
            for(int i=0;i<27;i++) if(pdg2name(SMcode[i])) sprintf(process+strlen(process),"%s,",pdg2name(SMcode[i]));
            process[strlen(process)-1]=0;
            sprintf(lib,"tranf_%s%s",lib1,lib2);                 
            numout*cc=getMEcode(0,ForceUG,process,NULL,NULL,lib);
            M[j*nP+i]+=collisionWidth(cc,Beps,T); 
            printf("collision width for %s is %E\n", process, collisionWidth(cc,Beps,T));
         }       
      }
      
      for(int i=1;i<nP;i++) for(int j=0;j<i;j++) M[i*nP+j]=M[j*nP+i]*Y[i]/Y[j];
      for(int i=0;i<nP;i++) { M[i*nP+i]=0; for(int j=0;j<nP;j++) if(i!=j) M[i*nP+i]-=M[j*nP+i]; }

      printf("Decay matrix improved by collisions:\n");
      for(int i=0;i<nP;i++) { for(int j=0;j<nP;j++) printf(" %E ",(double)M[i*nP+j]); printf("\n");}
      
      omg_eff=0;
      for(int k=1;k<=(1<<(nP-1))-1;k++)
      {  int q=2*k;
         double Ys1=0,Ys2=0,omgs=0;
         for(int i=0;i<nP;i++) for(int j=0;j<nP;j++) if(((1<<i) & k) && (!((1<<j) & k))) omgs+= M[i*nP+j]*Y[j];
         for(int i=0;i<nP;i++) if((1<<i)&k) Ys1+=Y[i]; else Ys2+=Y[i];
         omgs*=1/Ys1+1/Ys2;
         if(k==1) omg_eff=omgs; else if(omgs<omg_eff) omg_eff=omgs;
      }
      printf("omega_eff(2)=%.2E (with collisions)\n",omg_eff);
   }

extern double   detA( int N,long double *A);

//  for(int i=0;i<nP;i++) M[i*nP+1]+=1E-20;
//  printf("det=%E\n",  detA(nP,M));

/*
double ff(double x) 
{ 
  long double  MM[100];
  long double xx=x;
  for(int i=0;i<nP*nP;i++) MM[i]=M[i];
  for(int i=0;i<nP;i++) MM[i*nP+i]+=xx;
  int sgn=2*(nP&1)-1;
  return (sgn*detA(nP,MM)/x);
}       



   displayPlot("det(A-x)", "x", omg_eff/100,omg_eff*2,0, 1,"ff", 0,ff,NULL);
*/
//   displayPlot("det(A-x)", "x", Hubble(T)/100,Hubble(T)*10,1, 1,"ff", 0,ff,NULL);
   
//printf(" %E %E \n %E %E \n", (double)M[0], (double)M[1],(double)M[2],(double) M[3]);
         
   free(all); free(Y); free(ndf); free(ch); free(M);
   return omg_eff/Hubble(T);  
}


//  Functions for testing and visualisation 


double vs1120F(double T){ return 3.8937966E8*vsInterpolation( T,  d1120, vs1120T, 0, 0, NULL); }
double vs2200F(double T){ return 3.8937966E8*vsInterpolation( T,  d2200, vs2200T, 0, 0, NULL); }
double vs1100F(double T){ return 3.8937966E8*vsInterpolation( T,  d1100, vs1100T, 0, 0, NULL); }
double vs1210F(double T){ return 3.8937966E8*vsInterpolation( T,  d1210, vs1210T, 0, 0, NULL); }
double vs1110F(double T){ return 3.8937966E8*vsInterpolation( T,  d1110, vs1110T, 0, 0, NULL); }
double vs2220F(double T){ return 3.8937966E8*vsInterpolation( T,  d2220, vs2220T, 0, 0, NULL); }
double vs1220F(double T){ return 3.8937966E8*vsInterpolation( T,  d1220, vs1220T, 0, 0, NULL); }
double vs2210F(double T){ return 3.8937966E8*vsInterpolation( T,  d2210, vs2210T, 0, 0, NULL); }

double vs1122F(double T){ return 3.8937966E8*vsInterpolation( T,  d1122, vs1122T, 2, d2211, vs2211T ); }
double vs2211F(double T){ return 3.8937966E8*vsInterpolation( T,  d2211, vs2211T,-2, d1122, vs1122T ); }
double vs1112F(double T){ return 3.8937966E8*vsInterpolation( T,  d1112, vs1112T, 1, d1211, vs1211T ); }
double vs1211F(double T){ return 3.8937966E8*vsInterpolation( T,  d1211, vs1211T,-1, d1112, vs1112T ); }
double vs1222F(double T){ return 3.8937966E8*vsInterpolation( T,  d1222, vs1222T, 1, d2221, vs2221T ); }
double vs2221F(double T){ return 3.8937966E8*vsInterpolation( T,  d2221, vs2221T,-1, d1222, vs1222T ); } 

double Y1F(double T){  if(Ntab2<=0 || T<Ttab2[Ntab2-1] || T> Ttab2[0]) return NAN; return polint3(T,Ntab2,Ttab2, Y1T)+Yeq1(T);}
double Y2F(double T){  if(Ntab2<=0 || T<Ttab2[Ntab2-1] || T> Ttab2[0]) return NAN; return polint3(T,Ntab2,Ttab2, Y2T)+Yeq2(T);}



