#include<math.h>
#include<stdio.h>
#include <dirent.h>
#include <pthread.h>
#include"micromegas.h"
#include"micromegas_aux.h"

#include "../CalcHEP_src/include/rootDir.h"
#include"../CalcHEP_src/c_source/strfun/include/pdt.h"
#define dptDefault "NNPDF31_lo_as_0118"

char pdfName[50]={};

static  pdtStr *data[8]={};
static pdtList * allPDT=NULL;

static int pos(int pNum)
{ switch(pNum) 
  { case   5: case -5:  return 1;
    case   4: case -4:  return 2;
    case   3: case -3:  return 3;
    case  -1:           return 4;
    case  -2:           return 5;
    case  21: case -21: return 6;
    case   2:           return 7;
    case   1:           return 8;
    default:            return 0;
  }  
}

static double partonPdt(int pdg, double x, double q)
{ int n=pos(pdg);
  double res;
  
  static  pthread_mutex_t lockKey=PTHREAD_MUTEX_INITIALIZER;
      
  if(!pdfName[0])
  {     
    if(nPROCSS>1) pthread_mutex_lock(&lockKey);
    if(!pdfName[0])  setPDT(dptDefault);
    if(nPROCSS>1) pthread_mutex_unlock(&lockKey);
  }
  if(!n || ! data[n-1]) return 0;
  if(abs(pdg)==4 && q<1.4) return 0;
  if(abs(pdg)==5 && q<5  ) return 0;
  if(x< data[n-1]->x_min) x=data[n-1]->x_min;
  if(q< data[n-1]->q_min) q=data[n-1]->q_min;
  res=interFunc(x, q, data[n-1]);
  if(res<0) return 0; else return res;
 
}

static double alphaPdt(double q)
{ if(!pdfName[0])  setPDT(dptDefault);
  if(data[5]) return  interAlpha(q,data[5]); else return 0;
}


static void findAllPDT(void)
{ 
   char*fname;
   DIR *dirPtr;
   struct dirent * dp;
   pdtList * cpdt;
   static  pthread_mutex_t lockKey=PTHREAD_MUTEX_INITIALIZER;   
 
   if(nPROCSS>1) pthread_mutex_lock(&lockKey);
   if(allPDT) return;
   
   fname=malloc(strlen(calchepDir)+50);
   sprintf(fname,"%s/pdTables", calchepDir);
     
   dirPtr=opendir(fname);
   if(!dirPtr) { free(fname);  if(nPROCSS>1)pthread_mutex_unlock(&lockKey);  return;}
 
   while((dp=readdir(dirPtr)))
   { char *c=dp->d_name;
     int l=strlen(c);
     if(l>=4 && strcmp(c+l-4,".pdt")==0) 
     {  sprintf(fname,"%s/pdTables/%s",calchepDir,c); 
        makePdtList(fname, &allPDT);
     }
   } 
   closedir(dirPtr);
   free(fname);
   for(cpdt=allPDT;cpdt;cpdt=cpdt->next)
   { char *c=strchr(cpdt->name,'(');
     if(c) c[0]=0;
   }
   
   if(nPROCSS>1)pthread_mutex_unlock(&lockKey); 
}     

void PDTList(void)
{ 
  pdtList * cpdt;
  if(!allPDT)  findAllPDT();
  printf(" List of built-in PDT distributions\n");  
  for(cpdt=allPDT;cpdt;cpdt=cpdt->next)
  { 
     if(cpdt->beamP==2212) printf("%s\n",cpdt->name);  
  }                        
}


int setPDT(char*name)
{  pdtList * cpdt;
   int i; 
   if(!allPDT)  findAllPDT();
    
   for(cpdt=allPDT;cpdt;cpdt=cpdt->next) if(cpdt->beamP==2212 && strcmp(name,cpdt->name)==0)
   {  int pdg[8]={ 5  ,  4  ,  3  , -1  , -2  , 21  ,  2  ,  1};
      for(i=0;i<8;i++) if(data[i]){ freePdtData(data[i]); free(data[i]); data[i]=NULL;} 
      for(i=0;i<8;i++)
      { int j;
        for(j=0;cpdt->partons[j];j++) if(cpdt->partons[j]==pdg[i]) 
        { 
          data[i]=malloc(sizeof(pdtStr));
          getPdtData(cpdt->file,cpdt->items[j],data[i]); 
          break;
        }                                       
      }
      sprintf(pdfName,"PDT:%s",name);    
      parton_distr=partonPdt;
      parton_alpha=alphaPdt;
              
      return 0;    
   }
   printf("can not find PDT set  %s\n",name);
   return 1;
}

int restorePDF(char*oldPDF)
{
 
  if(strcmp(oldPDF,pdfName))
  {
    if(strstr(oldPDF,"PDT:")) setPDT(oldPDF+4);
    else if(strstr(oldPDF,"LHA:"))
    { int lhaMem;
      char lhaName[50];
      sscanf(oldPDF+4,"%[%:]:%d",lhaName,&lhaMem);
      setLHAPDF(lhaMem,lhaName);
    }
  }  
}


static double x0_,q_;
static int pc1_,pc2_;

static double conv_integrand(double y)
{ double x=exp(-y);
  return /*fabs*/(parton_distr(pc1_,x, q_)*parton_distr(pc2_,x0_/x, q_));
}

static  pdtStr *cData[36]={};

static  pdtStr *convStrFunAux(int pc1, int pc2)
{ pdtStr *D;
  int i,ix,iq,nc,i1=pos(pc1),i2=pos(pc2);

  if(i1>=i2) nc=i1*(i1-1)/2+i2-1; else nc=i2*(i2-1)/2+i1-1;
  if(cData[nc])  return cData[nc];
  
  static  pthread_mutex_t lockKey=PTHREAD_MUTEX_INITIALIZER;
  static char currentPdf[50]={"XXX"};
  
  if(nPROCSS>1) pthread_mutex_lock(&lockKey); 
  if(!pdfName[0])  setPDT(dptDefault);
  
  if(strcmp(currentPdf,pdfName)) 
  {
    for(i=0;i<36;i++) if(cData[i]) {freePdtData(cData[i]); free(cData[i]); cData[i]=NULL;}
    strcpy(currentPdf,pdfName);
  }  
  

  { int nq=20;
    double q_grid[20]={
     1.30000E+00, 1.53106E+00, 1.83098E+00, 2.22659E+00, 2.75766E+00, 3.48439E+00, 4.50000E+00, 6.12359E+00, 8.60159E+00, 1.25127E+01,
     1.89184E+01, 2.98475E+01, 4.93546E+01, 8.59491E+01, 1.58477E+02, 3.11214E+02, 6.55142E+02, 1.48904E+03, 3.68299E+03, 1.00000E+04
                      };
    int nx=96;
    double x_grid[96]={
     5.E-7      ,  1.00000E-06,  1.28121E-06,  1.64152E-06,  2.10317E-06,  2.69463E-06,  3.45242E-06,  4.42329E-06,  5.66715E-06,  7.26076E-06,
     9.30241E-06,  1.19180E-05,  1.52689E-05,  1.95617E-05,  2.50609E-05,  3.21053E-05,  4.11287E-05,  5.26863E-05,  6.74889E-05,  8.64459E-05,
     1.10720E-04,  1.41800E-04,  1.81585E-04,  2.32503E-04,  2.97652E-04,  3.80981E-04,  4.87518E-04,  6.26039E-04,  8.00452E-04,  1.02297E-03,
     1.30657E-03,  1.66759E-03,  2.12729E-03,  2.71054E-03,  3.44865E-03,  4.37927E-03,  5.54908E-03,  7.01192E-03,  8.83064E-03,  1.10763E-02,
     1.38266E-02,  1.71641E-02,  2.11717E-02,  2.59364E-02,  3.15062E-02,  3.79623E-02,  4.53425E-02,  5.36750E-02,  6.29705E-02,  7.32221E-02,
     8.44039E-02,  9.64793E-02,  1.09332E-01,  1.23067E-01,  1.37507E-01,  1.52639E-01,  1.68416E-01,  1.84794E-01,  2.01731E-01,  2.19016E-01,
     2.36948E-01,  2.55242E-01,  2.73927E-01,  2.92954E-01,  3.12340E-01,  3.32036E-01,  3.52019E-01,  3.72282E-01,  3.92772E-01,  4.13533E-01,
     4.34326E-01,  4.55495E-01,  4.76836E-01,  4.98342E-01,  5.20006E-01,  5.41818E-01,  5.63773E-01,  5.85861E-01,  6.08077E-01,  6.30459E-01,
     6.52800E-01,  6.75387E-01,  6.98063E-01,  7.20830E-01,  7.43683E-01,  7.66623E-01,  7.89636E-01,  8.12791E-01,  8.35940E-01,  8.59175E-01,
     8.82485E-01,  9.05866E-01,  9.29311E-01,  9.52817E-01,  9.76387E-01,  1.00000E+00
                      };
    D=malloc(sizeof(pdtStr));
D->index=0;
    D->nq=nq;
    D->q_grid=malloc(nq*sizeof(double));
    for(i=0;i<nq;i++)D->q_grid[i]=q_grid[i];
    
    D->lq_grid=(double*)malloc(sizeof(double)*nq);
    for(i=0;i<nq;i++)D->lq_grid[i]=log(D->q_grid[i]);
    

    D->nx=nx;
    D->x_grid=malloc(nx*sizeof(double));
    for(i=0;i<nx;i++)D->x_grid[i]=x_grid[i];
    D->lx_grid=malloc(nx*sizeof(double));
    for(i=0;i<nx;i++)D->lx_grid[i]=log(x_grid[i]);
    D->strfun=malloc(nx*nq*sizeof(double));    
    D->interpolation= int_cteq6; //biCubicLogXQ; //  int_cteq6; 
    D->mass=1;  
    D->beamP=2141;
    D->parton=nc;

    if(D->interpolation==int_cteq6)
    {
       D->x_grid_aux=malloc(nx*sizeof(double));
       for(i=0;i<D->nx;i++) D->x_grid_aux[i]=pow(D->x_grid[i],0.3);

       D->q_grid_aux=malloc(nq*sizeof(double));
       for(i=0;i<D->nq;i++) D->q_grid_aux[i]=log( D->lq_grid[i]-log(0.22) ); 
    }
    
    D->alpha=NULL; 
    D->x_min=1.E-6;
    D->q_min=1.3;
    D->q_max=10000; 
    D->pow0=D->pow1=0;       /* factor x^pow *(1-x)^pow1 applied after interpolation */ 
    D->q_threshold=1; /* threshold q for heavy quarks?           */
    D->nSmallX=D->nSmallQ=D->nLargeX=D->nLargeQ=0;
  
    D->qt0=1;        /* position of first q point above threshold*/
    D->approx=0;
  }  
  
  pc1_=pc1;
  pc2_=pc2; 
  for(iq=0;iq<D->nq;iq++)
  { q_ =D->q_grid[iq];
    for(ix=0;ix<D->nx;ix++)
    { 
      x0_=D->x_grid[ix];
      D->strfun[D->nx*iq+ix] = simpson(conv_integrand,0.,-log(x0_),1.E-3,NULL);
         
    }
  } 
  cData[nc]=D;
  if(nPROCSS>1) pthread_mutex_unlock(&lockKey);
  return D;  
}

double convStrFun2(double x, double q, int pc1, int pc2, int pp)
{ pdtStr *D1,*D2;
  int pc1c,pc2c;
  double strF;
  
  if(pp>0) pc2c=pc2; else pc2c=-pc2;
  D1=convStrFunAux(pc1,pc2c);
  strF=interFunc(x, q, D1);
  
  if(pc1!=pc2)
  { if(pp>0) pc1c=pc1; else pc1c=-pc1;
    D2=convStrFunAux(pc1c, pc2);
    if(D1==D2) strF*=2; else strF+=interFunc(x, q, D2);
  }
  return strF;
}


double convStrFun3(double x, double q, int pc1, int pc2, int pp)
{  pdtStr *D1,*D2;
   int pc1c,pc2c;
   double strF;
   q_=q; 
   x0_=x;  
   if(pp>0) pc2c=pc2; else pc2c=-pc2;
    pc1_=pc1;
    pc2_=pc2c;
    strF=simpson(conv_integrand,0.,-log(x),1.E-3,NULL);
              
    if(pc1!=pc2)
    { if(pp>0) pc1c=pc1; else pc1c=-pc1;
      pc1_=pc1;
      pc2_=pc2c;
      strF+=simpson(conv_integrand,0.,-log(x),1.E-3,NULL);
    }
    return strF;
}
                              

double (*parton_distr)(int pdg, double x,double q)=partonPdt;
double (*parton_alpha)(double q)=alphaPdt;

double alpha_2(double q) {return parton_alpha(q);}  /* it needs only  to avoid external link with functions from  num.a  */

//  For direct detection 

static double x_integrand(double x)
{  if(x==0) return 0; else return  x*parton_distr(pc1_,x,q_); }

double parton_x( int pNum, double  Q)
{
  double x1;
  if(Q>1.E4) Q=1.E4; 
  if(!pdfName[0])  setPDT(dptDefault);  
  q_=Q;

  pc1_=pNum;
  
  x1=simpson(x_integrand,1E-4,1.,1.E-4,NULL);  
  if(pNum==21) return x1;
  if(abs(pNum)>2) return 2*x1;
  pc1_=-pNum;
  return x1+simpson(x_integrand,1E-4,1.,1.E-4,NULL);
}

// FORTRAN

extern int   setpdt_(char *fname,int len);
int   setpdt_(char *fname,int len)
{ 
  char cname[50];
  fName2c(fname,cname,len);
  setPDT( cname);
}  

extern void pdtlist_(void);

void pdtlist_(void) { PDTList(); }
