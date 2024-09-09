/*
 Copyright (C) 2002, Alexander Pukhov pukhov@theory.sinp.msu.ru
*/


#include "syst.h"
#include "interface.h"
#include "crt_util.h"
#include "plot.h"
#include "kinaux.h"
#include "strfun.h"
#include "sf_pdt.h"
#include "read_func.h"
/*
#include "phys_val.h"
#include "parser.h"
#include "4_vector.h"
#include "VandP.h"
*/
#include "num_out.h"
#include "qcdScale.h"
#include "alphas2.h"


double (*sf_alpha[2])(double)={NULL,NULL};

static double alphaMZ=0.1184;
static int alphaOrder=2;
static int alphaNF=5;
static double MbMb=4.2;
static double Mtp=173;
static int alphaPDF=1;
static char Rscale_str[61]="91.187";
static char F1scale_str[61]="Qren";
static char F2scale_str[61]="Qpdf1";
static char Sscale_str[61]="Qren";

static double alpha(int nf, int odr, double lambda,  double dscale)
{
    double b0 = 11. -  (2./3.)*nf;
    double b1 = 51. - (19./3.)*nf;
    double b2 = 2857. - (5033./9.)*nf + (325./27.)*nf*nf;
    double rl = 2*log(dscale / lambda);
    double alpha0= 4*M_PI/(b0*rl);
    double d__4 = log(rl) - .5;
    double d__2 = 2*b1/(b0*b0*rl);


    if(odr==1) return alpha0;
    else if(odr==2) return  alpha0*(1 - 2*b1*log(rl)/(b0*b0*rl));
    else if(odr==3) return  alpha0*(1 - 2*b1*log(rl)/(b0*b0*rl)
         + d__2*d__2 *(d__4*d__4 + b2*b0 /(8*b1*b1) - 1.25)  );
    else { fprintf(stderr,"Can not evaluate alpha in so large oder (%d).\n",odr);
           exit(1);
         }
}

static double findLambda(int nf,int odr, double alpha0, double M)
{ double l1=0.1, l2=0.3;
  double l,a,a1,a2;

  while((a1=alpha(nf,odr,l1,M)-alpha0) > 0)  l1*=0.7;

  while((a2=alpha(nf,odr,l2,M)-alpha0) < 0) l2*=1.3;

  do{ l=(l1*a2-l2*a1)/(a2-a1);
      a=alpha(nf,odr, l,M)-alpha0;
      if(fabs(a1)>fabs(a2)){ a1=a;l1=l;} else {a2=a;l2=l;}
    } while (fabs(a) > 0.00001*alpha0);
  return l;
}

#define  MCHARM   1.3

static double L3, L4, L5, L6;
static int   nf3,nf4,nf5,nf6;

static void init_alpha(void)
{
  double lambda,al; 
  int nf=alphaNF;

  if(nf==6) nf--;
  
  lambda=findLambda(nf,alphaOrder,alphaMZ,91.187);

  if(alphaNF>6) alphaNF=6; if(alphaNF<3) alphaNF=3;
  if(alphaOrder>3)alphaOrder=3; if(alphaOrder<1)alphaOrder=1;

  nf5=nf; L5=lambda;
  if(nf5==5) 
  {  if(MbMb<1.5*L5) {L4=L5;nf4=nf5;} else 
     { al=alpha(5,alphaOrder,L5,MbMb  );
       if(alphaOrder==3) al*=(1.+(11./72.)*pow(al/M_PI,2.));
       L4=findLambda(4,alphaOrder,al,MbMb);nf4=4; ;
     }
  } else {nf4=nf;L4=lambda;}
  if(nf4==4) 
  { if(MCHARM<1.5*L4) {L3=L4;nf3=nf4;} else
    { al=alpha(4,alphaOrder,L4,MCHARM);
      if(alphaOrder==3) al*=(1.+(11./72.)*pow(al/M_PI,2.));
      L3=findLambda(3,alphaOrder,al,MCHARM);nf3=3;
    }
  } else {nf3=nf; L3=lambda;}

  if(alphaNF==6)
  {  nf6=6;
     al=alpha(5,alphaOrder,L5,Mtp);
     if(alphaOrder==3) al*=(1.-(11./72.)*pow(al/M_PI,2.));
     L6=findLambda(6,alphaOrder,al,Mtp);
  }
  else {nf6=nf5; L6=L5;}
}

static double alpha_0(double Q)
{
    if(Q<L3*1.5) Q=L3*1.5;
  
        if(Q<MCHARM)  return alpha(nf3, alphaOrder, L3, Q);
    else  if(Q<MbMb)    return alpha(nf4, alphaOrder, L4, Q);
    else  if(Q<Mtp)     return alpha(nf5, alphaOrder, L5, Q);
    else                return alpha(nf6, alphaOrder, L6, Q);
}

double alpha_2(double Q)
{
  switch(alphaPDF)
  { case 1: if(sf_alpha[0]) return (*sf_alpha[0])(Q);
            if(sf_alpha[1]) return (*sf_alpha[1])(Q);
            break;
    case 2: if(sf_alpha[1]) return (*sf_alpha[1])(Q);         
            if(sf_alpha[0]) return (*sf_alpha[0])(Q);
            break;
  }                                  
  
  return alpha_0(Q);
}


int qcdmen_(void)
{
   void * pscr=NULL;
   int mode;
   int returnCode=0;

   initStrFun(0);
L10:{  char strmen[]="\030"
      " parton dist. alpha OFF "
      " alpha(MZ)=  ZZZZ       "
      " nf =        NF         "
      " order=      NNLO       "
      " mb(mb)=     MbMb       "
      " Mtop(pole)= Mtp        "
      " Alpha(Q) plot          "
      " Qren = RRR             "
      " Qpdf1= FF1             "
      " Qpdf2= FF2             "
      " Qshow= FFS             ";

      if(alphaPDF)
      { int k=0;
//printf("alphaPDF=%d sf_alpha[0]=%p sf_alpha[1]=%p\n", alphaPDF,sf_alpha[0],sf_alpha[1]);      
        switch(alphaPDF)
        { case 1: if(sf_alpha[0]) k=1; else if(sf_alpha[1]) k=2; break;
          case 2: if(sf_alpha[1]) k=2; else if(sf_alpha[0]) k=1; break;
        }    
        if(k) improveStr(strmen,"OFF","pdf%d",k);
      }    
      improveStr(strmen,"ZZZZ","%.4f", alphaMZ);
      improveStr(strmen,"NF","%d",alphaNF);

           if(alphaOrder==1) improveStr(strmen,"NNLO","%-.4s","LO");
      else if(alphaOrder==2) improveStr(strmen,"NNLO","%-.4s","NLO");
      else alphaOrder=3;
      
      improveStr(strmen,"MbMb","%.3f", MbMb);
      improveStr(strmen,"Mtp","%.2f", Mtp);

      improveStr(strmen,"RRR","%-.16s", Rscale_str);
      
      improveStr(strmen,"FF1","%-.16s", F1scale_str);
      improveStr(strmen,"FF2","%-.16s", F2scale_str); 
      improveStr(strmen,"FFS","%-.16s", Sscale_str);
      menu1(54,8,"QCD alpha",strmen,"n_alpha",&pscr,&mode);
    }
    switch (mode)
    { case 0: if(returnCode) init_alpha(); return returnCode;
      case 1: 
         { char alphaMen[100]="\006"                                                                                    
           "  OFF "                                                                                                     
           " pdf1 "                                                                                                     
           " pdf2 ";                                                                                                    
           int k=0;                                                                                                     
           menu1(54,12,"alpha",alphaMen,"",NULL,&k);                                                                    
           if(k)alphaPDF=k-1;                                                                                           
           if(alphaPDF && !sf_alpha[alphaPDF-1]) messanykey(20,20,"WARNING! This pdf  does not define alphaQCD ");                                             
         }                                                                                                              
         break;
      case 2:
         { double alphaMZ_old=alphaMZ; 
           if(correctDouble(3,15,"Enter new value ",&alphaMZ,1)) returnCode=1;
           if(alphaMZ>0 && alphaMZ<0.3) returnCode=1; else  
           { alphaMZ=alphaMZ_old;
              messanykeyErr(5,15,"Your input is out of alphaMZ range");
           }
         }
         break;
      case 3:
         { int NF_old=alphaNF;
           if(correctInt(3,15,"Enter new value ",&alphaNF,1))
           {
              if(alphaNF<=6 && alphaNF>=3) returnCode=1;
              else { messanykeyErr(5,15,"NF out of range"); alphaNF=NF_old;}
           }   
         }
         break;
      case 4:
         {  char lomen[]="\010"
               " LO     "
               " NLO    "
               " NNLO   ";
            void *pscrlo=NULL;
            int k=0;
            menu1(52,12,"",lomen,"",&pscrlo,&k);
            if(k) { alphaOrder=k; returnCode=1; put_text(&pscrlo);}
         }
         break;
      case 5: correctDouble(3,15,"Enter new value ",&MbMb,1);  break;
      case 6: correctDouble(3,15,"Enter new value ",&Mtp,1);   break;   
      case 7:
	{ void * screen;
	  int i;

	  static double qMin=1, qMax=1000;
	  static int nPoints=100;
	  if(returnCode) init_alpha();
	  get_text(1,1,maxCol(),maxRow(),&screen);

          if(correctDouble(40 ,15 ,"Q_min=",&qMin,0)&& qMin>=0.5
          && correctDouble(40 ,16 ,"Q_max=",&qMax,0)&& qMax>qMin
          && correctInt(33,17,"number of points=" ,&nPoints,0)
          && nPoints>3&& nPoints<=150)
	  { double *f[3]={NULL,NULL,NULL};
	    double *ff[3]={NULL,NULL,NULL};
	    char buff[3][100]; 
	    char* Y[3]={buff[0],buff[1],buff[2]};
	    switch(alphaOrder)
	    { case 1: sprintf(Y[0],"MSbar LO");   break;
	      case 2: sprintf(Y[0],"MSbar NLO");  break;
	      case 3: sprintf(Y[0],"MSbar NNLO"); break;
	      default:sprintf(Y[0],"MSbar");
	    }  
	       
	    int N,k;
	    f[0]=(double*) malloc(nPoints*sizeof(double));
            int xLog= (qMin>0 && qMax/qMin >10)? 1 : 0; 

	    for(i=0;i<nPoints;i++)
	    { double z=(i+0.5)/(double)(nPoints),q;
	      if(xLog) q=pow(qMin,1-z)*pow(qMax,z); else q=qMin*(1-z)+qMax*z;
	      f[0][i]=alpha_0(q);
//	      printf("i=%d %E\n",i,f[0][i]);
	    }
	    
	    N=1;	    
	    for(k=0;k<2;k++) if(sf_alpha[k])
	    {  
	       char buff[300];  
               strFunName(k+1,Y[N]);
               char *p=strstr(Y[N],"(proton");
               if(p) p[0]=0; else
               { p=strstr(Y[N],"(anti-proton");
                 if(p) p[0]=0;
               }  
	       f[N]=(double*) malloc(nPoints*sizeof(double));
               for(i=0;i<nPoints;i++)
               { double z=(i+0.5)/(double)(nPoints),q;
                 if(xLog) q=pow(qMin,1-z)*pow(qMax,z);
                    else  q=qMin*(1-z)+qMax*z;   
                 f[N][i]=(*sf_alpha[k])(q);
	       } 
	       N++;                                                      
            }
            int dim[3];
            for(i=0;i<N;i++) dim[i]=nPoints;
            plot_Nar(NULL, "Alpha(Q)", "Q[GeV]",  qMin,qMax, xLog, N, dim,f,ff,Y);
            for(k=0;k<N;k++) free(f[k]);
	  } else  messanykey(40,18,
	          " Correct input is \n"
	          " 0.5<= Q_min <Q_max\n"
	          " number of points <=150 and >=4");
	  put_text(&screen);
	} break;	

      case 8:
         { int npos=1,rc;
           do
           { 
              char mess[200];
              goto_xy(2,12); print("Renorm. scale: ");
              if(str_redact(Rscale_str,npos,60)==KB_ENTER) returnCode=1;
              goto_xy(2,12); clr_eol();
              rc=initScales(Rscale_str,F1scale_str,F2scale_str,Sscale_str, mess);
              if(rc) messanykey(10,10,mess);
           }  while(rc);
         }
	 break;


      case 9:
         { int npos=1,rc;
           do
           { 
              char mess[200];
              goto_xy(2,12); print("Fct1.scale: ");
              if(str_redact(F1scale_str,npos,60)==KB_ENTER) returnCode=1;
              goto_xy(2,12); clr_eol();
              rc=initScales(Rscale_str,F1scale_str,F2scale_str,Sscale_str, mess);
              if(rc) messanykey(10,10,mess);
           }  while(rc);
         }
	 break;
      case 10:
         { int npos=1,rc;
           do
           { 
              char mess[200];
              goto_xy(2,12); print("Fct1.scale: ");
              if(str_redact(F2scale_str,npos,60)==KB_ENTER) returnCode=1;
              goto_xy(2,12); clr_eol();
              rc=initScales(Rscale_str,F1scale_str,F2scale_str,Sscale_str, mess);
              if(rc) messanykey(10,10,mess);
           }  while(rc);
         }
	 break;
      case 11:
         { int npos=1,rc;
           do
           { 
              char mess[200];
              goto_xy(2,12); print("Shworing scale: ");
              if(str_redact(Sscale_str,npos,60)==KB_ENTER) returnCode=1;
              goto_xy(2,12); clr_eol();
              rc=initScales(Rscale_str,F1scale_str,F2scale_str,Sscale_str,mess);
              if(rc) messanykey(10,10,mess);
           }  while(rc);
         }
	 break;
	 
    }
    goto L10;
}

int w_alphaQCD(FILE * mode)
{
  fprintf(mode,"alphaPDF=%d alpha(MZ)=%E NF=%d Order=%d MbMb=%E Mtp=%E",
                alphaPDF, alphaMZ, alphaNF, alphaOrder, MbMb, Mtp);
  return 0;
}

int w_Scales(FILE * mode)
{
  fprintf(mode,"\n Renormalization %s\n Factorization1 %s\n Factorization2 %s\n Showering %s\n",
  Rscale_str,F1scale_str,F2scale_str,Sscale_str);
  return 0;
}
                  

void i_alphaQCD(void) { init_alpha();}

void i_Scales(void)
{
  int err;  
  char mess[200];
 
  if(nin_int==1) strcpy(Rscale_str,"M1"); else strcpy(Rscale_str,"M12");
  
  strcpy(F1scale_str,"Qren");
  strcpy(F2scale_str,"Qpdf1");
  strcpy(Sscale_str,"Qren");
  err=initScales(Rscale_str,F1scale_str,F2scale_str,Sscale_str,mess);
}

int r_alphaQCD(FILE *mode)
{ int n;
  int err;
  n=fscanf(mode, "alphaPDF=%d alpha(MZ)=%lf NF=%d Order=%d MbMb=%lf Mtp=%lf",
            &alphaPDF, &alphaMZ,&alphaNF,&alphaOrder,&MbMb,&Mtp);
  if(n!=6) return 1;
  init_alpha();
  return 0;
}

int r_Scales(FILE *mode)
{ int n,err; 
  char mess[200];

  n=fscanf(mode," Renormalization %[^\n] Factorization1 %[^\n] Factorization2 %[^\n] Showering %[^\n]",Rscale_str,F1scale_str,F2scale_str,Sscale_str);
  if(n!=4) return 1; 
  trim(Rscale_str);
  trim(F1scale_str);
  trim(F2scale_str);
  err=initScales(Rscale_str,F1scale_str,F2scale_str,Sscale_str,mess);
  if(err) i_Scales();
  return 0;
}
