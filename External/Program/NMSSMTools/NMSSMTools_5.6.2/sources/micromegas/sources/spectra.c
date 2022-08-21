#include"micromegas.h"
#include"micromegas_aux.h"
#include"micromegas_f.h"

#define NEn 28  /* Number of Energy points for interpolation */
#define Nin 16 
#define Nout 7
#define QCUT 25 //2.5
#define X1CUT (2.E-2)

#define V0 (vRot/299792.*1.5957691)  /* 2*sqrt(2/PI).  10^(-6)*Mslp < Width */

//#define V0 1E-3
 
//#define DISPLAY_SPECTRA 
/*#define ECTEST*/

static double Mcdm0;

static char* errMess=NULL;



aChannel* vSigmaCh=NULL;
static int nAnCh=0;
static double GG=1.23;

static float phidiff[Nin][Nout][NEn][NZ];


double Zi(int i) { return log(1.E-7)*pow((double)(i-1)/(double)(NZ),1.5) ;}
static  int Iz(double z) { return NZ*pow(z/log(1.E-7),1./1.5)+1; }

void  fillSpect(double (*dNdE)(double ), double Emax, double * SpectAr)
{ int i;
  SpectAr[0]=Emax;
  for(i=1;i<NZ;i++) 
  { double E=exp(Zi(i))*Emax;  
    SpectAr[i]=dNdE(E)*E;
  }
}




// copied from CalcHEP/c_source/dynamicME/kin4.c
static double kinematic_1_3(REAL *pmass, int i3, double m12, double xcos, REAL * P)
{ 
  double factor;
  REAL pout,chY,shY,xsin, E1,P12,P13,E2,P22,P23, m0,m1,m2,m3;
  int i,i1,i2;
  
  for(i=1;i<4;i++)if(i3!=i) {i1=i; break;}
  for(i++;i<4;i++)if(i3!=i) {i2=i; break;}
  
  m0=pmass[0];
  m1=pmass[i1];
  m2=pmass[i2];
  m3=pmass[i3];

  if(m12<=m1+m2) return 0;
  for(i=0;i<16;i++) P[i]=0;

  P[0]=m0; 
  factor=1/(64*M_PI*M_PI*M_PI*m0*m0);
  
  pout=decayPcm(m0,m12,m3);
  if(!pout) return 0;  
  P[i3*4]=Sqrt(pout*pout+m3*m3); P[i3*4+3]=-pout; 

  factor*=pout;  
  
  shY=pout/m12;
  chY=Sqrt(1+shY*shY);  
  pout=decayPcm(m12,m1,m2);
  if(!pout) return 0;
  factor*=pout;
  xsin=Sqrt(1-xcos*xcos);
  E1=Sqrt(m1*m1+pout*pout);    E2=Sqrt(m2*m2+pout*pout);
  P13=xcos*pout;               P23=-P13;
  P12=xsin*pout;               P22=-P12;
  
  P[4*i1]  =chY*E1 + shY*P13;  P[4*i2]  =chY*E2 + shY*P23;
  P[4*i1+3]=shY*E1 + chY*P13;  P[4*i2+3]=shY*E2 + chY*P23;
  P[4*i1+2]=P12;               P[4*i2+2]=P22;
  
  return factor;
}

typedef struct
{  
   numout*cc;
   REAL mass[4];
   int i3;
   int nsub;
   double m12;
   double GG;
   int err;
}   argFor13Int;


static double dWidth13dCos(double xcos, void*arg )
{  REAL pvect[16];
   REAL mass[4];
   argFor13Int*arg13=arg;
   
   double factor=kinematic_1_3(arg13->mass,arg13->i3,arg13->m12,xcos, pvect);
   if(factor==0) return 0;
   return   factor*(arg13->cc->interface->sqme)(arg13->nsub,arg13->GG, pvect,NULL,&(arg13->err));
}

static double dWidth13dM(double m12, void*arg )
{
   ((argFor13Int*)arg)->m12=m12;
   int err;
   return simpson_arg(dWidth13dCos, arg, -1, 1, 1E-2,&err);
}



/* v*cs22 at v=0 */  
double  vcs22(numout * cc,int nsub,int * err)
{
   int i;
   double pcm,r;
   REAL pmass[4], pvect[16];
   
//printf("Mb4a=%E\n", findValW("Mb"));    
   for(i=1;i<=cc->interface->nvar;i++) if(cc->link[i]) cc->interface->va[i]=*(cc->link[i]);
//printf("Mbb=%E Q=%E\n", findValW("Mb"),  findValW("Q")); 
   if( cc->interface->calcFunc()>0 ) {*err=4;  return 0;}
//printf("Mbb_=%E Q=%E\n", findValW("Mb"),  findValW("Q"));
   *(cc->interface->gtwidth)=0;
   *(cc->interface->twidth)=0;
   *(cc->interface->gswidth)=0;
//printf("Mb4c=%E\n", findValW("Mb"));    
   for(i=0;i<4;i++)  cc->interface->pinf(nsub,1+i,pmass+i,NULL);   
   *err=0;
   if(pmass[0]+pmass[1] <= pmass[2]+pmass[3]) return 0;
   for(i=0;i<16;i++) pvect[i]=0;
//printf("Mb4d=%E\n", findValW("Mb")); 
   pcm= decayPcm(pmass[0]+pmass[1],pmass[2],pmass[3]);
   for(i=0;i<2; i++) pvect[4*i]=pmass[i];
   for(i=2;i<4; i++) pvect[4*i]=Sqrt(pmass[i]*pmass[i] +pcm*pcm);
   pvect[8+3]=pcm;
   pvect[12+3]=-pcm;
   r=cc->interface->sqme(nsub,GG,pvect,NULL,err);
   return 3.8937966E8*r*pcm/(16*M_PI*pmass[0]*pmass[1]*(pmass[0]+pmass[1]));  
}
 
/* New 2->3 */ 
static REAL pmass[5],pvect[20];
static int code[5];
static int iA,ix,iX;
static numout* cc23;
static double X0=0.01, dSigmadE_x0, Egamma,X1=0.8 ,dSigmadE_x1,dSigmadE_x1_e;
static int PrintOn=0;
static double dSigmadCos23(double csfi)
{
  REAL pcm1,pcm2,ms,md,chY,shY,M;
  double  r;
  int i,err_code;

  for(i=0;i<20;i++) pvect[i]=0;

  pvect[0]=pmass[0];
  pvect[4]=pmass[1];

  pvect[4*iA]= Egamma;
  pvect[4*iA+3] = -Egamma;  

  ms=pmass[ix]+pmass[iX];
  md=pmass[iX]-pmass[ix];
  pcm1=Egamma;

  M=4*pmass[0]*(pmass[0]-Egamma);
  if(M<=ms*ms) return 0;
  
  M= Sqrt(M);
  
  pcm2=Sqrt((M*M-ms*ms)*(M*M-md*md))/M/2;
  
  pvect[4*ix]=Sqrt(pmass[ix]*pmass[ix]+pcm2*pcm2);
  pvect[4*ix+3]=pcm2*csfi;
  pvect[4*ix+2]=pcm2*Sqrt(1-csfi*csfi);
  
  pvect[4*iX]= Sqrt(pmass[iX]*pmass[iX]+pcm2*pcm2);
  pvect[4*iX+3]=-pvect[4*ix+3];
  pvect[4*iX+2]=-pvect[4*ix+2];  

  chY=Sqrt(1+pcm1*pcm1/M/M);
  shY= pcm1/M;  
  
  { double p0,p3;
    p0=pvect[4*ix], p3=pvect[4*ix+3];
    pvect[4*ix]=  chY*p0 + shY*p3;
    pvect[4*ix+3]=shY*p0 + chY*p3;

    p0=pvect[4*iX]; p3=pvect[4*iX+3];
    pvect[4*iX]=  chY*p0 + shY*p3;
    pvect[4*iX+3]=shY*p0 + chY*p3;
  }
  
/*
  for(i=0;i<4;i++)
  { int j; 
    double s=0;
    s=pvect[i]+pvect[4+i]-pvect[8+i]-pvect[12+i]-pvect[16+i] ;
    printf("s(%d)=%Egamma\n",i,s);
  }              
*/  
err_code=0;
  r=(cc23->interface->sqme)(1,GG,pvect,NULL,&err_code)*pcm1*pcm2/(128*M*pmass[0]*pmass[0]*M_PI*M_PI*M_PI);
if(r<0) r=0;
  if(err_code) return 0;
//printf("csfi=%E r=%e\n", csfi,r);  
  return r*3.8937966E8;
}

static double dSigmadFi23(double fi)
{if(fi==0 || fi==M_PI) return 0; return dSigmadCos23(cos(fi))*sin(fi);}

static double dSigmadE(double E)
{ double r;
   Egamma=E;
/*
   if(pmass[ix]>1.E-3*Mcdm0 && pmass[iX]>1.E-3*Mcdm0  ) 
                       r= simpson(dSigmadFi23,0,M_PI,1.E-4,NULL);
   else { printf("dSigmadCos23\n"); //displayFunc("dSigmadCos23",dSigmadCos23, -1 , 1,0,"dSigmadCos23");
                 r= simpson(dSigmadCos23,-0.9,0.9,1.E-4,NULL); printf("ok\n");  
*/
//                   r= gauss(dSigmadCos23,-1,1,7); }
   r= gauss(dSigmadCos23,-1,1,7);                 
//   r= simpson(dSigmadFi23,0,M_PI,1.E-4,NULL);
   return r;
}

static int addGamma(int pdg)
{ pdg=abs(pdg); 
  if(pdg<=6) return 0;
  if(pdg>=11 && pdg<=15) return 0;
  if(pdg==81||pdg==83) return 0;
  return 1;
}  

static double dSigmadERest(double E)
{
  double res, eps,ms= (pmass[ix]+pmass[iX]);
  int l,nn[2]={ix,iX};
  
  if(E+sqrt(E*E + ms*ms)>=1.999999*Mcdm0) return 0;
  res= dSigmadE(E); if(res==0.) return 0;
  eps=ms/2/Mcdm0;

//  if( !(addGamma(code[iX]) && addGamma(code[ix])))
  { double subtract=0,norm=0,x=E/Mcdm0;
    double csmax,pcm,pcm0;
    
     if(pmass[ix]>1.E-3*Mcdm0 && pmass[iX]>1.E-3*Mcdm0 ) csmax=1; else csmax=0.999;
     pcm=decayPcm(2*Mcdm0*Sqrt(1-x),pmass[ix],pmass[iX]);
     for(l=0;l<2;l++)  
     {  
       double kappa;
       if(!addGamma(code[nn[l]]))
       { kappa=sqrt(1/(1+pow(pmass[nn[l]]/pcm,2)));
         subtract += log((1+kappa*csmax)/(1-kappa*csmax))*(1-x+x*x/2 -eps*eps/2)-kappa*(1-x)*csmax;
       }  
       pcm0=decayPcm(2*Mcdm0,pmass[ix],pmass[iX]);
       kappa=sqrt(1/(1+pow(pmass[nn[l]]/pcm0,2)));
       norm+=        log((1+csmax*kappa)/(1-csmax*kappa))*(1 -eps*eps/2)-kappa*csmax;
     }
     res-=subtract/norm*dSigmadE_x0/x; 
  }
  if(res<0) return 0;
  
  return res;
}

static double Xe;
static double FE(double cs)
{ 
  Egamma=2*Mcdm0*(1-Xe)/(1+cs);
  if(Egamma<X0*Mcdm0) return 0;
  if(Egamma/Mcdm0 -(1-Xe) < -1.E-4) return 0;
  return dSigmadCos23(cs)*2/(1+cs);
}

static double FEfi(double fi)
{ 
  if(fi< 1.E-4  || fi > M_PI - 1.E-4) return 0;
  return FE(cos(fi))*sin(fi);
}

static double dSigmadEe(double E)
{
   double xe,x1=0.96,x2=0.98;
   xe=Xe=E/Mcdm0; 
   if(xe> x2+1.E-4)
   { 
     return   ((x2-xe)*dSigmadEe(x1*Mcdm0) + (xe-x1)*dSigmadEe(x2*Mcdm0))/(x2-x1);
   }           
   if(Xe>X0 )
   { double csMin=(2*(1-Xe) -1)*1.0001;
     double csMax=(2*(1-Xe)/X0 -1);
     double r;
/*     csMax=1;     

        if(csMin<-0.98) csMin=-0.98;
        if(csMax> 0.98) csMax=0.98;
        return simpson(FE,csMin,csMax,1.E-3,NULL);

        if(csMin<-1) csMin=-1;
        if(csMax> 1) csMax=1;
        return gauss(FE,csMin,csMax,7);     
*/
     if(csMin<-1) csMin=-1;
     if(csMax> 1) csMax=1;
     csMin+=0.001; // !!!!!!!! Problem in  case of zero mass of electron. Reason is not clear 
int err;
     r= simpson(FEfi,acos(csMax),acos(csMin),1.E-4,NULL);
//if(err) displayPlot("dSigmadEe", "fi",acos(csMax),acos(csMin),0,1,"FEfi",0,FEfi,NULL);
     return r;
   } else return 0;
}
#ifdef DISPLAY_SPECTRA
static double EdSigmadEe(double E) {return E*dSigmadEe(E);}
static double EdSigmadERest(double E) {return E*dSigmadERest(E);}
#endif

static double Spectra22A(char*name1,char*name2,double ** Spectra,txtList plusA)
{
  double vcs=0;
  int i,l;
  txtList cRec;
  int spin2Dm;
  
  for(l=0;l<6;l++) if(Spectra[l]) { Spectra[l][0]=Mcdm0; for(i=1;i<NZ;i++) Spectra[l][i]=0;}
  if(Spectra[0]==NULL && Spectra[1]==NULL) return 0;

if(Spectra[0])   
  for(l=0;l<nModelParticles;l++)
  {
     if(ModelPrtcls[l].NPDG==22) 
     { outNames[0]=ModelPrtcls[l].name;
       break;
     }
  } 
  if(l==nModelParticles) return 0;
  
  for(cRec=plusA;cRec;cRec=cRec->next)
  { 
     char lib[50]="";
     char*N[5];
     process2Lib(cRec->txt,lib);
     cc23=getMEcode(0,0,cRec->txt,NULL,NULL,lib);
     if(cc23==NULL) continue;
     if(passParameters(cc23)) continue; 
     *(cc23->interface->gtwidth)=0;
     *(cc23->interface->twidth)=0;
     *(cc23->interface->gswidth)=0;

     for(i=0;i<5;i++) N[i]=cc23->interface->pinf(1,i+1,pmass+i,code+i);
    
     cc23->interface->pinfAux(1,1,&spin2Dm,NULL,NULL,NULL);
     
     for(ix=0,i=2;i<5;i++) if(code[i]==22) iA=i; else if(!ix)ix=i;else iX=i;
     if(Spectra[0] && dSigmadE(X1*Mcdm0) > X1CUT*dSigmadE_x1 )
     {  double x2Sum; 
        dSigmadE_x0=3*X0*(dSigmadE(X0*Mcdm0)+dSigmadE(3*X0*Mcdm0)-2*dSigmadE(2*X0*Mcdm0));
#ifdef DISPLAY_SPECTRA
     {  char buff[100];
        sprintf(buff,"d(vSigma(%s))/dE(A) [pb/GeV]",cRec->txt);      
        displayFunc(dSigmadERest, Mcdm0*X0, Mcdm0,buff);
     }
#endif 
       x2Sum=0; 
       for(i=1;i<NZ;i++) 
       { double dSigmaDz=0,x=exp(Zi(i));
         if(x>X0)
         { dSigmaDz=x*Mcdm0*dSigmadERest(x*Mcdm0);
           Spectra[0][i]+=dSigmaDz;
           x2Sum+=(Zi(i)-Zi(i+1))*dSigmaDz;
         }else  Spectra[0][i]+=dSigmaDz;
       }
       { 
         vcs+=x2Sum;
         vSigmaCh=realloc(vSigmaCh, (nAnCh+2)*sizeof(aChannel));
         vSigmaCh[nAnCh].weight=x2Sum;   
         { int j; for(j=0;j<5;j++) vSigmaCh[nAnCh].prtcl[j]=N[j]; }
         nAnCh++;                     
       }
     }
     if(Spectra[1] && abs(code[iX])==11 && code[ix]+code[iX]==0  /*&& pmass[ix]==0*/ &&   // !!! was blocked for Me!=0. Reason not clear.  
       code[0]==code[1] && spin2Dm==1 && dSigmadEe(X1*Mcdm0) > X1CUT*dSigmadE_x1_e  ) 
     {
        
#ifdef DISPLAY_SPECTRA 
  displayFunc(dSigmadEe, Mcdm0*X0, Mcdm0," electron spectrum");
{ double csA,csE,xcsA,xcsE; 
  csA=simpson(dSigmadE,Mcdm0*X0,Mcdm0,1.E-3,NULL);
  csE=simpson(dSigmadEe,Mcdm0*X0,Mcdm0,1.E-3,NULL);
  xcsA=simpson(EdSigmadERest,Mcdm0*X0,Mcdm0,1.E-3,NULL);
  xcsE=simpson(EdSigmadEe,Mcdm0*X0,Mcdm0,1.E-3,NULL);
/*  
  printf("vcs(A)= %E  vcs(E)= %E\n",csA,csE);
  printf("energy  fraction lost %E\n", (2*Mcdm0- xcsA/csA-2*xcsE/csA)/Mcdm0);  
*/  
}        
#endif
       for(i=1;i<NZ;i++)
       {
           double Ee=Mcdm0*exp(Zi(i));
           if(Ee>X0*Mcdm0) Spectra[1][i]+=Ee*dSigmadEe(Ee);
       }
     }         
  }  
  return vcs;
}

/* ------------------------------------------------- */


static int readSpectra(void)
{ static int rdOk=0;
  int k,l,i,n;
  FILE *f;
  char * buff;
  char *fnames[Nin]={"gg.dat","dd.dat","uu.dat","ss.dat","cc.dat",
                     "bb.dat","tt.dat","ee.dat","mm.dat","ll.dat",
                     "zz.dat","zz_t.dat","zz_l.dat",
                     "ww.dat","ww_t.dat","ww_l.dat"};

  if(rdOk) return 0;
  buff=malloc(strlen(micrO)+100);

  for(n=0;n<Nin;n++)
  {  sprintf(buff,"%s/sources/data/%s",micrO,fnames[n]);
     f=  fopen(buff,"r");
//     printf("f=%p\n",f);       
//     fclose(f); f=NULL;
     
     if(f==NULL) { free(buff);return 1;} 
//printf("%s\n", buff);        
     for(i=0;i<6;i++)
     { fscanf(f,"%*s"); 
       for(k=0;k<NEn;k++)
       {
         for(l=0;l<NZ;l++) if(1!=fscanf(f,"%f",phidiff[n][i][k]+l)) break;
         fscanf(f,"%*s"); 
         for(;l<NZ;l++)phidiff[n][i][k][l]=0;
       }
     }
//printf("ok before fclose\n");     
     fclose(f);
  }
  
  free(buff);
  rdOk=1;   

#ifdef ECTEST
printf("Energy conservation test\n");
for(n=0;n<Nin;n++) for(k=0;k<NEn;k++)
{
    double mi[NEn]={2,5,10,25,50,80.3,85,91.2,92,95,100,110,120,125,130,140,150,176,200,250,350,500,750,1000,1500,2000,3000,5000};

  double x[6];
  for(n=0;n<Nin;n++) for(k=0;k<NEn;k++)
  { double sum=0;


  for(l=0;l<6;l++) x[l]=0;
  for(l=0;l<6;l++) if(phidiff[n][l][k])
  { 
     for(i=0;i<NZ;i++)
     { double xx= pow(1.E-7, pow((i+0.5)/NZ,1.5));
       double dx;
       if(i==0) dx=-0.5*log(1.E-7)*pow((double)(1.)/NZ,1.5);
       else dx= 0.5*log(1.E-7)*( pow((double)(i-1)/NZ,1.5) - pow((double)(i+1)/NZ,1.5) );
       if(l==2) xx+=1/mi[k];
       x[l]+=phidiff[n][l][k][i]*dx*xx;   
     }
  }  
  sum=x[0]+2*(x[1]+x[2]+x[3]+x[4]+x[5]);
   
  if(sum && (sum>2.1 || sum<1.9))   
   printf("Channel = %s Energy=%.1E  %.2e + %.2e + %.2e + %.2e + %.2e + %.2e = %.2E\n",fnames[n],mi[k], x[0],x[1],x[2],x[3],x[4],x[5],sum);
  }
}  
   
#endif  

  
  return 0;
}


double zInterp(double zz, double * tab) 
{  
   double dz,r;
   int j0;   
   if(zz>0) return 0;
   
   j0=Iz(zz); 
   if(j0<1) j0=1;
   if(j0>=NZ-1) return tab[NZ-1];
   
   dz= (zz-Zi(j0))/(Zi(j0+1)-Zi(j0));
   r=(1-dz)*tab[j0]+dz*tab[j0+1];
   if(r<0)r=0;
   return r; 
}



static void mInterp(double Nmass,  int  CHin,int  CHout, double*tab)
{  
//  float mi[NEn]={10,25,50,80.3,91.2,100,150,176,200,250,350,500,750,1000,1500,2000,3000,5000};
   double mi[NEn]={2,5,10,25,50,80.3,85,91.2,92,95,100,110,120,125,130,140,150,176,200,250,350,500,750,1000,1500,2000,3000,5000};
   int l,i0;
   double c0,c1;
   float *p0,*p1;
   for(i0=0; i0<NEn && Nmass>=mi[i0] ;i0++);
   if(i0) i0--;

   switch(CHin)
   { 
     case 14: case 15: case 16: case 17: if(i0<5)
      { i0=6; for(l=1;l<NZ;l++) tab[l]=phidiff[CHin][CHout][i0][l-1]; tab[0]=Nmass; return; }
     case 10: case 11: case 12: case 13: if(i0<7)
      { i0=9; for(l=1;l<NZ;l++) tab[l]=phidiff[CHin][CHout][i0][l-1]; tab[0]=Nmass; return; }
   } 
   p0=phidiff[CHin][CHout][i0];
   p1=phidiff[CHin][CHout][i0+1];
   if(i0==NEn-1 || Nmass<=2) for(l=1;l<NZ;l++) tab[l]= p0[l-1];
   else
   {
     c1=(Nmass*Nmass -mi[i0]*mi[i0])/(mi[i0+1]*mi[i0+1] - mi[i0]*mi[i0]);
     c0=1-c1;
     for(l=1;l<NZ;l++) tab[l]= c0*p0[l-1]+c1*p1[l-1];
   }
   tab[0]=Nmass;
}



static double zInterpE(double x, double *tab ){ return    zInterp(x,tab)*exp(x);}  


double  spectrInfo(double Emin,double*tab,double*Etot)
{
  if(Etot) *Etot=0;
  if(Emin>=tab[0]) return 0;  else 
  { int i1,i2;
    double Xmin=Emin/tab[0],zmin,zmax=0;
    if(Xmin<1.22E-7) Xmin=1.22E-7;
    zmin=log(Xmin); 
    for(i1=Iz(zmin); i1>1 && tab[i1]==0;i1--) continue;
    if(i1<NZ-1) i1++;
    if(zmin<Zi(i1)) zmin=Zi(i1)+1E-6;

    
    for(i2=1; tab[i2]==0 && i2<NZ-1 ;i2++ );
    if(i2>1) i2--;
    if(zmax>Zi(i2)) zmax=Zi(i2)-1E-6;

    if(zmax<zmin) return 0;    
  
    if(Etot)*Etot=tab[0]*simpson_arg(zInterpE,tab, zmin, zmax,1.E-4,NULL);       
    return simpson_arg(zInterp,tab, zmin, zmax,1.E-4,NULL);
  }  
}


double spectrInt(double Emin,double Emax, double * tab)
{ 
  double M=tab[0], zmin,zmax;
  int i1,i2;

  if(Emin<M*exp(Zi(NZ-1))) zmin=Zi(NZ-1); else if( Emin>=M) zmin=0; else zmin=log(Emin/M);
  if(Emax<M*exp(Zi(NZ-1))) zmax=Zi(NZ-1); else if( Emax>=M) zmax=0; else zmax=log(Emax/M);
  
  if(zmin>=zmax) return 0;
  
  for(i1=Iz(zmin) ;i1>1 && tab[i1]==0;i1--) continue;
  if(i1<NZ-1) i1++;
  if(zmin<Zi(i1)) zmin=Zi(i1)+1E-6;
 
  i2=Iz(zmax)+1; if(i2>NZ-1) i2=NZ-1;
  for( ;i2<NZ-1 && tab[i2]==0;i2--) continue;
  if(i2<2) i2--;
  if(zmax>Zi(i2)) zmax=Zi(i2)-1E-6;
  
  return simpson_arg(zInterp,tab, zmin, zmax,1.E-4,NULL);
}


void  spectrMult( double *spect, double(*func)(double))
{ 
  int i; 
  double M=spect[0];
  for(i=1;i<NZ;i++)
  { double E=M*exp(Zi(i));
    spect[i]*=func(E);
  }  
}


typedef struct
{ double m;
  double * tab;
} boostParStr;




static double FUNB(double e, void *B)
{ 
  boostParStr * par=B;
  double r;

  if(e<=0) return 0;
  r=zInterp(log(e/par->tab[0]),par->tab)/(e*sqrt(e*(e+2*par->m)));
  return r;
}



void boost(double Y, double M0, double mx, double*tab)
{ double chY=cosh(Y), shY=sinh(Y);
  int l,k;
  double tab_out[NZ];  

  boostParStr  par;
  par.m=mx;
  par.tab=tab;
  
  if(Y<0.01) return;
  
//m12=1.073733E-03 pcm=5.249327e-02 Y=4.582795e+00  tab[0]=5.368667e-04 cosh(Y)*tab[0]=2.625212E-02   M0=5.250000E-02 M=7.046556E-02

  double e=tab[0], p=sqrt(e*e-mx*mx), M=e*chY+p*shY-mx;
//  printf("Y=%e  tab[0]=%e cosh(Y)*tab[0]=%E   M0=%E M=%E\n",Y,tab[0],cosh(Y)*tab[0], M0,M);
  if(M0<M) M0=M;     
 
  double emin=tab[0]*exp(Zi(NZ-1));
  int kmin=NZ-1;
  for(; tab[kmin]==0 && kmin>0; kmin--) continue;
  if(kmin==0) { tab[0]=M0;  return;}
  if(kmin<NZ-1) kmin++;  // tab[kmin]==0, tab[kmin+1]!=0;
   
  for(l=1;l<NZ;l++)
  { double e=M0*exp(Zi(l)),e1,e2;
    if(mx==0) { e1=e*exp(-Y); e2=e*exp(Y);}  
    else 
    {
      if(e>mx) 
      {
        double YY=acosh(1+e/mx);
        e1=mx*(cosh(Y-YY)-1);
        e2=mx*(cosh(Y+YY)-1);
//printf("Y =%E YY=%E  e1=%E e2=%E\n", Y,YY,e1,e2);        
      }else 
      {
        double p=sqrt(e*(e+2*mx)); 
        e1=chY*(e+mx)-shY*p-mx;
        e2=chY*(e+mx)+shY*p-mx;
      }
    }
//printf("l=%d e=%E e1=%E e2=%E\n",l,e,e1,e2);    
    if(e1>=tab[0]) {tab_out[l]=0; continue;}
    if(e2>tab[0]) e2=tab[0];
    
    if(e1 < emin) k=NZ-1; else k=Iz(log(e1/tab[0]))+1;
    if(k>kmin) k=kmin;
    e1=tab[0]*exp(Zi(k));
        
    if(e1>=tab[0]) {tab_out[l]=0; continue;}    
    if(e2>tab[0]) e2=tab[0];
    if (e1>=e2) tab_out[l]=0; else tab_out[l]=e*simpson_arg(FUNB,&par, e1,e2,1.E-3,NULL)/2/shY;
  }
  
  for(l=1;l<NZ;l++) tab[l]=tab_out[l];
  tab[0]=M0;
//  printf("tab[0]=%E M0=%E\n",tab[0],M0);
}

static double outMass[6]= {0.,0.511E-3,0.939,0.,0.,0.};
char* outNames[6]={"gamma","e+","p-","nu_e","nu_mu","nu_tau"};

int basicSpectra(double Mass, int pdgN, int outN, double * tab)
{ 
//  printf("basicSpectra: pdg=%d outN=%d\n", pdgN, outN);
  tab[0]=Mass;
  for(int i=1;i<NZ;i++) tab[i]=0;
  int N=abs(pdgN);
      
   
  int inP=-1;

  switch(N)
  { case 21:inP=0; break;  /*glu*/ 
    case 1: case 2: case 3: case 4: case 5: case 6:  inP=N; break; /* d,u,s,c,b,t*/
    case 11:case 13: case 15: inP=(N+1)/2+1; break; /*e,m,l*/  
    case 23:    inP=10; break;  /*z*/
    case 23+'T':inP=11; break;  
    case 23+'L':inP=12; break;     
    case 24:    inP=13; break;  /*w*/
    case 24+'T':inP=14; break;
    case 24+'L':inP=15; break;
    case 12:case 14: case 16: case 22:
    {  if( (N==22 && outN==0) ||  (N-12 == 2*(outN-3)))
       {
         tab[1]=2/(-Zi(3));
         tab[2]=2/(-Zi(3))*( 1-Zi(2)/Zi(3));       
       }
       if(N==22) {tab[1]*=2;tab[2]*=2;}
       return 0;
    }
  }
  if(inP==-1) { tab[0]=Mass; for(int i=1;i<NZ;i++) tab[i]=0;    return 1;}
  readSpectra(); 
  mInterp(Mass,inP,outN,tab);
  if(pdgN==11)
  { 
     if(outN==0 && Mass<2)
     {  double me=5.11E-4;
        if(Mass<me) for(int i=1;i<NZ;i++) tab[i]=0; else 
        { 
          for(int i=1;i<NZ;i++)
          {
             double E=Mass*exp(Zi(i));
             tab[i]=2*E*FSRdNdE(E,Mass,me,1,1);    
          }  
        } 
     }
  }    
  return 0;
}

/*static*/ void getSpectrum2(int wPol, double M, char*n1,char*n2,int outP, double *tab);

void  decaySpectrum(char*pName,int outP, double*tabD)
{ 
  double M=pMass(pName); 
  double mx=outMass[outP];
  
  tabD[0]=M/2;
  for(int i=1;i<NZ;i++)tabD[i]=0;

/*          
  if( (CDM1 &&(strcmp(CDM1,pName)==0 ||strcmp(aCDM1,pName)==0))
    ||(CDM2 &&(strcmp(CDM2,pName)==0 ||strcmp(aCDM2,pName)==0))) return; 
*/    
  txtList L;   
  double w=pWidth(pName,&L);
  if(w==0) 
  { char  mess[100]; sprintf(mess," Can not decay '%s'; ",pName);
    addErrorMess(&errMess,mess);
    return;    
  }     
  
//double E,Es=0;
  double brs=0;          
  for(;L;L=L->next)
  { char p[5][20];
    double tab_p[NZ];
    double br;
    int n=sscanf(L->txt,"%lf %s -> %[^, ], %[^,], %[^,], %[^, ]",&br, p[0],p[1],p[2],p[3],p[4]);    
    n-=2;
//printf("%s %E\n",L->txt,br);    
    if(n==4) continue;
    else if(n==2) 
    { getSpectrum2(0,M,p[1],p[2],outP, tab_p);

//      spectrInfo(1E-10,tab_p,&E);
//      Es+=E*br;
//      printf("E=%e \n", E);
    }  
    else if(n==3) 
    {  char process[50];
       sprintf(process,"%s->%s,%s,%s",p[0],p[1],p[2],p[3]); 
//printf("process=%s\n",process);       
       numout*cc13=newProcess(process);
       passParameters(cc13);
       argFor13Int arg13;
       arg13.cc=cc13;
       arg13.nsub=1;
       for(int i=0;i<4;i++) arg13.cc->interface->pinf(1,i+1,arg13.mass+i,NULL); 
       arg13.GG =sqrt(4*M_PI*alphaQCD(arg13.mass[0]));
                   
       tab_p[0]=tabD[0];
       for(int i=1;i<NZ;i++) tab_p[i]=0;
       
       int i1,i2,i3;
       i3=1;
       
       for(int i=2;i<4;i++)
       {
          int pdg=abs(pNum(p[i]));
          switch(pdg)
          { case 22: if(outP==0) continue; else break;
            case 11: if(outP==1) continue; else break;
            case 12: if(outP==3) continue; else break;
            case 14: if(outP==4) continue; else break;
            case 16: if(outP==5) continue; else break; 
          }  
          if(arg13.mass[i]>arg13.mass[i3]) i3=i;
       }  
//i3=1;    
//printf("p[i3]=%s\n",p[i3]);   
       for(i1=1;;i1++) if(i1!=i3) break;
       for(i2=i1+1;;i2++) if(i2!=i3) break; 
       double mMin=0;
       for(int i=1; i<4;i++) if(i!=i3) {mMin+=arg13.mass[i]; if(abs(pNum(p[i]))<3) mMin+=0.14;} 
       if(mMin<2*mx) mMin=2*mx;
       double mMax=arg13.mass[0]-arg13.mass[i3]; if(abs(pNum(p[i3]))<3) mMax-=0.14;
       arg13.i3=i3;
       if(mMin>=mMax) continue;
       double wt=0;
       int nInt=20;       
       for(int k=0;k<nInt;k++)
       { 
         double m12=mMin+ (k+0.5)*(mMax-mMin)/nInt;
         double dw=dWidth13dM(m12, &arg13)*(mMax-mMin)/nInt;
          
         wt+=dw;
         double pcm=decayPcm(arg13.mass[0], arg13.mass[i3], m12);            
         double tab12[NZ],tabD[NZ];

         getSpectrum2(0,m12,p[i1],p[i2],outP, tab12);
   
         double Y=asinh(pcm/m12);
         boost(Y, M/2, outMass[outP], tab12);
         if(isSMP(pNum(p[i3])))
         {  double m=arg13.mass[i3];
            double e=sqrt(pcm*pcm+m*m);  
            basicSpectra(e,abs(pNum(p[i3])),outP,tabD);
            for(int i=1;i<NZ;i++) tabD[i]/=2;
         } 
         else 
         {  
            decaySpectrum(p[i3],outP, tabD);
            Y=asinh(pcm/pMass(p[i3]));
            boost(Y, M/2, outMass[outP], tabD);
         }  
           addSpectrum(tab12, tabD);
           for(int i=1;i<NZ;i++) tab12[i]*=dw;  
           addSpectrum(tab_p, tab12);    
       }
//printf("wt=%e\n", wt);       
       for(int i=1;i<NZ;i++) tab_p[i]/=wt;
    }
    for(int i=1;i<NZ;i++) tab_p[i]*=br; 
    {
       double buff[NZ];
       for(int i=0;i<NZ;i++) buff[i]=tabD[i];
       addSpectrum(tabD,tab_p);
       
//if(strcmp(p[1],"nl")==0)        displayPlot("Spectra", "E", 38, tabD[0],0,3,"tabD",0, SpectdNdE,buff,"tab_p",0,SpectdNdE,tab_p,"tabD+",0,SpectdNdE,tabD);
//       spectrInfo(1E-10,tabD,&E);
//       printf("Ess=%E\n",E);
    }    
    brs+=br;
  }      
}


void getSpectrum2(int pol, double M,char*n1,char*n2,int outP, double *tab)
{
  double m1=pMass(n1);
  double m2=pMass(n2);
  int N1=pNum(n1);
  int N2=pNum(n2);  
//printf("getSpectrum2: %s %s -> %d\n", n1,n2,outP);
  int i;
  int inP=-1;
  int N;

  if(M<=m1+m2) pol=0;
  if(pol!='T' && pol!='L') pol=0;
  tab[0]=M/2;
  for(i=1;i<NZ;i++) tab[i]=0;
   
  if( (isSMP(N1) &&  N1+N2==0 ) || (N1==21 && N2==21) || (N1==22 && N2==22) ||  (N1==23 && N2==23)   )
  {   N=abs(N1);
      if(N==23 || N==24)  N+=pol;    
      basicSpectra(M/2, N, outP, tab);
      return;
  }     
         
    char* nn[2];
    double mm[2];
    double E[2]; 
    double pcm;
    int k;

    nn[0]=n1;
    nn[1]=n2;

    mm[0]=m1;
    mm[1]=m2;
    double tabAux[NZ];
    if(M>m1+m2) pcm=decayPcm(M,m1,m2);
    else 
    { pcm=0;
      if(abs(N1)==abs(N2)) { mm[0]=M/2; mm[1]=M/2;} else
      {
        if(N1==23 || abs(N1)==24) mm[0]=M-m2; else mm[1]=M-m1;
      }  
    } 
    E[0]=sqrt(mm[0]*mm[0]+pcm*pcm);
    E[1]=sqrt(mm[1]*mm[1]+pcm*pcm);
    
    for(i=1;i<NZ;i++) tab[i]=0;
    
    for(k=0;k<2;k++)
    { N=pNum(nn[k]);
      if(isSMP(N))
      {  if(N==23||N==24) N+=pol;
         basicSpectra(E[k],N,outP,tabAux);
         for(int i=1;i<NZ-1;i++) tabAux[i]/=2;
         addSpectrum(tab, tabAux);
      }     
      else 
      {  if(mm[k]==0)
         {
           char  mess[100]; sprintf(mess, "No hadronization  for massless '%s'  ",nn[k]);
           addErrorMess(&errMess,mess);
         }  
         else
         {
          decaySpectrum(nn[k],outP, tabAux);
          double Y=asinh(pcm/mm[k]);
          boost(Y, M/2, outMass[outP], tabAux); 
          for(i=1;i<NZ;i++)tab[i]+=tabAux[i];
        }   
      } 
    }  
}


void addSpectrum(double *Spect, double * toAdd)
{
  double m1=Spect[0];
  double m2=toAdd[0];
  double buff[NZ];
  int i;
  
  if(fabs(m1-m2)<5E-3*(fabs(m1)+fabs(m2))) { for(int i=1; i<NZ; i++) Spect[i]+=toAdd[i]; return; }
  

  if(m1>m2) for(i=NZ-1;i>0;i--) 
  { 
     double E= m1*exp(Zi(i));
     if(E>m2) return; 
     Spect[i]+=  SpectdNdE(E,toAdd)*E;  
  } else  
  { for(i=0;i<NZ;i++) buff[i]=Spect[i];
    for(i=0;i<NZ;i++) Spect[i]=toAdd[i];
    for(i=NZ-1;i>0;i--) 
    {  
       double E= m2*exp(Zi(i)); 
       if(E>m1) return;
       Spect[i]+=  SpectdNdE(E,buff)*E;
    } 
  }
} 



static double calcSpectrum0(char *name1,char*name2, int key, double **Spectra, txtList*plusA)
{
  int i,k,l;
  double vcsSum=0; 
  int ntot,err;
  double * v_cs=NULL;
  char * photonName=NULL;
  char name1L[10],name2L[10], lib[20],process[4000];
  numout * libPtr;

    
  for(l=0;l<6;l++) if(Spectra[l]) { Spectra[l][0]=Mcdm0; for(i=1;i<NZ;i++) Spectra[l][i]=0;}  

  pname2lib(name1,name1L);
  pname2lib(name2,name2L);
  


int K_max= (CDM1&&CDM2)?4:1; 

for(int K=0;K<K_max;K++)
{

// Warning!!   in should be done in the same manner as annihilation libraries for Omega
  switch(K)
  {  
     case 0:
       sprintf(lib,"omg_%s%s",name1L,name2L); 
       sprintf(process,"%s,%s->AllEven,1*x{%s",name1,name2,EvenParticles());
     break;

     case 1: 
       if(pMass(name1)+pMass(name2) <= 2*Mcdm1 +1 ) continue;
       sprintf(lib,"omg_%s%s_1",name1L,name2L);
       sprintf(process,"%s,%s->AllOdd1,AllOdd1{%s",name1,name2,OddParticles(1));
     break;

     case 2:
       if(pMass(name1)+pMass(name2) <= 2*Mcdm2+1) continue;
       sprintf(lib,"omg_%s%s_2",name1L,name2L);
       sprintf(process,"%s,%s->AllOdd2,AllOdd2{%s",name1,name2,OddParticles(2));
     break;

     case 3:    
       if(pMass(name1)+pMass(name2) <= Mcdm1+Mcdm2+1) continue;      
       sprintf(lib,"omg_%s%s_12",name1L,name2L);  
       sprintf(process,"%s,%s->AllOdd1,AllOdd2{%s{%s",name1,name2,OddParticles(1),OddParticles(2));                
     break; 
  }


  libPtr=getMEcode(0,ForceUG,process,NULL,NULL,lib);

  if(!libPtr) continue;
  if(Qaddress && *Qaddress!=pMass(name1)+pMass(name2)) 
  { *Qaddress=pMass(name1)+pMass(name2);
     calcMainFunc;
  }   
  passParameters(libPtr);
  if(plusA)
  { 
    for(l=0;l<nModelParticles;l++)
    {
      if(ModelPrtcls[l].NPDG==22) 
      { photonName=ModelPrtcls[l].name;
        break;
      }
    }  
    if(!photonName) plusA=NULL;
  } 

  procInfo1(libPtr,&ntot,NULL,NULL); 

  
  v_cs=malloc(sizeof(double)*ntot);
  (*libPtr->interface->twidth)=0;
  for(k=0;k<ntot;k++)
  { REAL m[4];
    double wV,br;
    char *N[4];
    int pdg[4];
    int l,l_;
    
    for(i=0;i<4;i++) N[i]=libPtr->interface->pinf(k+1,i+1,m+i,pdg+i);
    cc23=NULL;
    v_cs[k]=0;
//printf("%s %s -> %s %s\n", N[0],N[1],N[2],N[3]);    
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
      { int err;
        r=v0*cs23(cc23,1,v0*Mcdm0/2,i3W,&err)/br;
        if(err) printf("error in simpson spectra.c line 902\n");
      }  
      if(pdg[l_]==23 || abs(pdg[l_])==24)
      { double wV2;
        
        wV2=pWidth(N[l_],NULL);
        r*=decayPcmW(2*Mcdm0,m[l],m[l_],wV,wV2,0)/decayPcmW(2*Mcdm0,m[l],m[l_],wV,0.,0);
        if(pdg[l]==pdg[l_]) r/=2;
      }
      v_cs[k]=r; 
      vcsSum+=r;                     
    }
    else  if(m[2]+m[3]< m[0]+m[1])
    {  err=0;
#ifdef V0    
      v_cs[k]=V0*cs22(libPtr,k+1,V0*m[0]*m[1]/(m[0]+m[1]),-1.,1.,&err); 
#else 
      v_cs[k]= vcs22(libPtr,k+1,&err);
#endif 
     if(v_cs[k]<0) v_cs[k]=0; 
      vcsSum+=v_cs[k];
    } else v_cs[k]=-1;
  }

   
  for(k=0;k<ntot ;k++) if(v_cs[k]>=0)
  { char * N[4];
    REAL m[4];
    int l,charge3[2],spin2[2],pdg[2];
    int PlusAok=0;


    procInfo2(libPtr,k+1,N,m);
    
    for(l=0;l<2;l++)  pdg[l]=qNumbers(N[2+l],spin2+l,charge3+l,NULL);

    if(Spectra[0] && plusA && (charge3[0] || charge3[1])&& (m[2]+m[3]< m[0]+m[1])) 
    {
       double m1=m[2], m2=m[3], Eg=X1*Mcdm0;
       double kappa=4*Mcdm0*(Mcdm0-Eg), ms=m1+m2, md=m1-m2;
       if(ms*ms<kappa)
       {  double dp=(Mcdm0-Eg/2)*sqrt((1-ms*ms/kappa)*(1-md*md/kappa));
          double p1= Eg/2*(1+ms*md/kappa)-dp, E1=sqrt(p1*p1+m1*m1);
          double p2= Eg/2*(1-ms*md/kappa)-dp, E2=sqrt(p2*p2+m2*m2);
             
          double Q1=Mcdm0*Mcdm0+m1*m1-2*Mcdm0*E1;
          double Q2=Mcdm0*Mcdm0+m2*m2-2*Mcdm0*E2;             

          double  m_min=10*Mcdm0;
          int n,m,w;
          char *s;
/*             
printf("energy conservation:0=%E=%E\n",  (E1+Eg+sqrt(pow(p2+2*dp,2)+m2*m2))/2/Mcdm0-1   ,
                                                   (E2+Eg+sqrt(pow(p1+2*dp,2)+m1*m1))/2/Mcdm0-1 );
*/                      
          for(n=1;(s=libPtr->interface->den_info(k+1,n,&m,&w,NULL));n++)
          {   double mass=0; if(m) mass=fabs( libPtr->interface->va[m]);
                if( ( strcmp(s,"\1\3")==0  || strcmp(s,"\1\4")==0) && m  && m_min> mass) m_min=mass;
          }

          if(  m_min*m_min -Q1  < QCUT*Mcdm0*Mcdm0*abs(charge3[0])/3.  
            || m_min*m_min -Q2  < QCUT*Mcdm0*Mcdm0*abs(charge3[1])/3. 
            || (abs(pdg[0])==24 && abs(pdg[1])==24  && Mcdm0 > 500 && v_cs[k]>1.E-3*vcsSum ))
          { txtList new22A=malloc(sizeof(txtListStr));
               new22A->next=*plusA; 
               new22A->txt=malloc(50);
               *plusA=new22A;
               sprintf(new22A->txt ,"%s,%s->%s,%s,%s",N[0],N[1],N[2],N[3],photonName); 
               PlusAok=1; 
          }
       } 
    }
            
    if(v_cs[k]>0) 
    {  double tab2[NZ]; 
       int N3=pdg[0], N4=pdg[1];

//       if(PrintOn )
//       { char txt[100];
//         sprintf(txt,"%s,%s -> %s %s", N[0],N[1],N[2],N[3]);
//         printf("  %-20.20s  %.2E\n",txt,v_cs[k]*2.9979E-26);
//       }
       
       vSigmaCh=realloc(vSigmaCh, (nAnCh+2)*sizeof(aChannel));
       vSigmaCh[nAnCh].weight=v_cs[k];
       { int j; 
         for(j=0;j<4;j++) vSigmaCh[nAnCh].prtcl[j]=N[j];
         vSigmaCh[nAnCh].prtcl[4]=NULL;
       }
       nAnCh++;
              
       { double lng=0;
         int noP=-1;
         if(key&1 && (abs(pdg[0])==23 || abs(pdg[0])==24|| abs(pdg[1])==23 || abs(pdg[1])==24)) noP=vPolar(N, &lng);
//if((key&1) && noP>=0 ) printf(" %s %s %s %s noP=%d lng=%E \n",  N[0],N[1],N[2],N[3],noP,lng   );
         
         for(l=0;l<6;l++) if(Spectra[l])
         {
           if(noP) getSpectrum2(0,m[0]+m[1],N[2],N[3],l,tab2);
           else
           {
             double tabT[NZ],tabL[NZ]; 
             getSpectrum2('T',m[0]+m[1],N[2],N[3],l,tabT);
             getSpectrum2('L',m[0]+m[1],N[2],N[3],l,tabL);   
             for(int i=1;i<NZ;i++) tab2[i]=(1-lng)*tabT[i]+lng*tabL[i];
             tab2[0]=tabT[0];
           }
           for(i=1;i<NZ;i++) Spectra[l][i]+=tab2[i]*v_cs[k];
         }
       }
#ifdef ADDFSR

       if(Spectra[0] && charge3[0])
       {
          for(l=0;l<2;l++) if(addGamma(pdg[l])&& m[2+l]!=0.) for(i=1;i<NZ;i++)
          {  double x=exp(Zi(i));
             if(2*Mcdm0*sqrt(1-x) > m[2]+m[3])
             {  double pcm,kappa,one_kappa,f;
                pcm=decayPcm(2*Mcdm0*sqrt(1-x), m[2],m[3]); 
                kappa=sqrt(1/(1+m[2+l]*m[2+l]/pcm/pcm)); 
                if( m[2+l]/pcm > 1.E-2) one_kappa=1-kappa; 
                                   else one_kappa=0.5*m[2+l]*m[2+l]/pcm/pcm; 
               f=(1./137.)*charge3[l]*charge3[l]/9/M_PI*log((1+kappa)/one_kappa);
          
               if(spin2[l]&1) Spectra[0][i]+=f*(1-x*(1-x*0.5))* v_cs[k];
               else           Spectra[0][i]+=f*(1-x          )* v_cs[k];
             }  
          }
       }
#endif       
    }    
  }
  free(v_cs);
}  
  return  vcsSum;
}


double calcSpectrum(int key, double *Sg,double*Se, double*Sp, double*Sne,double*Snm,double*Snl, int *errcode)
{ int n,i,j,l,err;
  double vcs;
  char  lop[20];
  double *Spectra[6],*Spectra_[6];
  double  buff[6][NZ];
  char * name, *aname;
  txtList plusA=NULL;
  txtList*plusAptr;
  
  if(errMess) { free(errMess); errMess=NULL;}
  
  
  char *WINPS[4]={CDM1,NULL,CDM2,NULL};
  double weight[4][4];
  int checkGam[4][4];
  double  chStat[3]={0,0,0};
  double NfracCDM2;
  
  err=readSpectra(); 
  if(err) { printf("calcSpectrum: Can not read data files for spectra\n");
            if(errcode) *errcode=-1; 
            return 0;
          }
          
//printf("?WINPS[2]=%s  fracCDM2=%e \n",WINPS[2],fracCDM2);
 
  if(fracCDM2==0) { WINPS[2]=NULL; NfracCDM2=0; } else
  if(fracCDM2==1) { WINPS[0]=NULL; NfracCDM2=1; } else
   NfracCDM2=  fracCDM2*Mcdm1/(fracCDM2*Mcdm1 +(1-fracCDM2)*Mcdm2); 

//printf("WINPS[2]=%s\n",WINPS[2]);
  if(WINPS[0])
  { int n=pTabPos(WINPS[0]);
    if( strcmp(ModelPrtcls[n-1].name,ModelPrtcls[n-1].aname)) WINPS[1]=ModelPrtcls[n-1].aname;
  }
  
  if(WINPS[2])
  { int n=pTabPos(WINPS[2]);
    if( strcmp(ModelPrtcls[n-1].name,ModelPrtcls[n-1].aname)) WINPS[3]=ModelPrtcls[n-1].aname;
  }
  
  double w1[4]={0.5*(1-NfracCDM2)*(1+dmAsymm)  ,0.5*(1-NfracCDM2)*(1-dmAsymm),
                0.5*NfracCDM2*(1+dmAsymm),      0.5*NfracCDM2*(1-dmAsymm)};
  if(!WINPS[1]) { w1[0]+=w1[1]; w1[1]=0;}
  if(!WINPS[3]) { w1[2]+=w1[3]; w1[3]=0;}

  for(i=0;i<4;i++) for(j=0;j<4;j++)
  {  weight[i][j]=w1[i]*w1[j];
     checkGam[i][j]=0;
     if(weight[i][j])
     { int ni=pTabPos(WINPS[i]);   
       int nj=pTabPos(WINPS[j]);
       if(ModelPrtcls[ni-1].spin2==0 && ModelPrtcls[nj-1].spin2==0) checkGam[i][j]=1;       
       if( WINPS[i]==WINPS[j]  && ModelPrtcls[nj-1].spin2==1)    checkGam[i][j]=1;
     } 
  }
  
    
  
   Mcdm0=0;
   if(weight[0][0]) Mcdm0=Mcdm1;
   if(weight[2][2] && Mcdm2>Mcdm0) Mcdm0=Mcdm2;

//for(i=0;i<4;i++) printf(" WINP=%s weight=%e\n", WINPS[i],w1[i]);  

  Spectra[0]=Sg; Spectra[1]=Se; Spectra[2]=Sp,Spectra[3]=Sne,Spectra[4]=Snm; Spectra[5]=Snl;
  for(l=0;l<6;l++) if(Spectra[l]) Spectra_[l]=buff[l];  else Spectra_[l]=NULL;
  for(l=0;l<6;l++)if(Spectra[l]) { Spectra[l][0]=Mcdm0; for(i=1;i<NZ;i++) Spectra[l][i]=0;} 
  vcs=0;


  nAnCh=0;
  vSigmaCh=realloc(vSigmaCh, 2*sizeof(aChannel));
  if(errcode) *errcode=0;
    
  if(key&4) { PrintOn=1; printf("    Channel          vcs[cm^3/s]\n");} else PrintOn=0;  

  for(i=0;i<4;i++) for(j=i;j<4;j++) if(weight[i][j])
  { double scale,c=weight[i][j];
    int nAnCh0,k;
    int key0=key&1;
    if(key&2 && checkGam[i][j]) key0+=2;
    if(i!=j) c*=2;

//printf("i=%d j=%d w=%e \n",i,j,weight[i][j]);
     
    Mcdm0=0.5*(pMass(WINPS[i])+pMass(WINPS[j]));
    err=findVal("Q",&scale);
    if(err==0 && scale!=2*Mcdm0) { assignValW("Q",2*Mcdm0); calcMainFunc();}
    
    if(key&2) plusAptr=&plusA; else plusAptr=NULL; 
    nAnCh0=nAnCh;
    { double vcsij=c*calcSpectrum0(WINPS[i],WINPS[j], key0, Spectra_,plusAptr);
      if( i<2 && j<2) chStat[0]+=vcsij;
      else if(i>1 && j>1) chStat[2]+=vcsij;
      else chStat[1]+=vcsij;
      vcs+=vcsij;      
    }  
    for(l=0;l<6;l++) if(Spectra[l])
    {  for(k=1;k<NZ;k++)Spectra_[l][k]*=c;
      addSpectrum(Spectra[l], Spectra_[l]);
       
    }
    for(k=nAnCh0;k<nAnCh;k++) vSigmaCh[k].weight*=c;     
    if(plusA)
    {  nAnCh0=nAnCh;
       if(Spectra[0]) dSigmadE_x1=zInterp(log(X1),Spectra_[0])/(X1*Mcdm0);  
       if(Spectra[1]) dSigmadE_x1_e=zInterp(log(X1),Spectra_[1])/(X1*Mcdm0);
       vcs+=c*Spectra22A(name,aname,Spectra_,plusA);

          
       for(l=0;l<6;l++) if(Spectra[l])
       { 
          for(k=1;k<NZ;k++)Spectra_[l][k]*=c;
          addSpectrum(Spectra[l], Spectra_[l]);
       }
        
       for(k=nAnCh0;k<nAnCh;k++) vSigmaCh[k].weight*=c;
       cleanTxtList(plusA);
       plusA=NULL;
    }    
  } 

  vSigmaCh[nAnCh].weight=0;
  for(j=0;j<5;j++) vSigmaCh[nAnCh].prtcl[j]=NULL;
  for(i=0;i<nAnCh-1;)
  {  if(vSigmaCh[i].weight >= vSigmaCh[i+1].weight) i++; 
     else
     {  aChannel buff;
        buff=vSigmaCh[i+1];vSigmaCh[i+1]=vSigmaCh[i];vSigmaCh[i]=buff;
        if(i)i--;else i++;
     }
  }           
 
  if(vcs)
  {  for(l=0;l<6;l++) if(Spectra[l])for(i=1;i<NZ;i++)Spectra[l][i]/=vcs;
     for(i=0;i<nAnCh;i++) vSigmaCh[i].weight/= vcs;
  } 
  
  
  vcs*=2.9979E-26;
  

  if(PrintOn )
  { int i=0;
    char txt[100];
    printf("==================================\n annihilation cross section %.2E cm^3/s\n",vcs  );
    printf(" contribution of processes\n");     
    for(i=0;vSigmaCh[i].weight>1.E-4;i++)
    {
    sprintf(txt,"%s,%s -> %s %s ", vSigmaCh[i].prtcl[0],
                                   vSigmaCh[i].prtcl[1],
                                   vSigmaCh[i].prtcl[2],
                                   vSigmaCh[i].prtcl[3]);
    if(vSigmaCh[i].prtcl[4]) strcat(txt,vSigmaCh[i].prtcl[4]);                               
    printf("  %-20.20s  %.2E\n",txt,vSigmaCh[i].weight);
    
     
    }
//    printf("chanStat: 11 -> %.2E  12-> %.2E 22->%.2E\n",    2.9979E-26*chStat[0],2.9979E-26*chStat[1],2.9979E-26*chStat[2]);
    
    
  }     

                                                                     
  if(Mcdm0 < 2) printf("WARNING! Spectra obtained at Mcdm=2GeV are used !\n");	 

  if(errMess) { printf(" %s\n", errMess); free(errMess); errMess=NULL;} 

  return vcs; 
}


double SpectdNdE(double E, double *tab){
  double z;
  if(E>tab[0]) return 0;
  z=log(E/tab[0]); 
  return 1/E*zInterp(z,tab); 
}


double FSRdNdE(double E, double p,double m, double q, int spin2)
{
    if(E>p || E<=0) return 0;
    double x=E/p;
    
    double kappa=sqrt(1/(1+m*m/p/p));
    double f=q*q/M_PI/137.*log((1+kappa)/(1-kappa))/E;
    if(spin2&1) return f*(1-x+x*x/2); else return f*(1-x);     
}
