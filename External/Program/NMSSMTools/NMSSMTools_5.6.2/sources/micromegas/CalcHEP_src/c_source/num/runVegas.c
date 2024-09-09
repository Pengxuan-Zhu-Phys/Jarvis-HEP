/*
 Copyright (C) 2002, Alexander Pukhov
*/

#include<unistd.h>
#include<stdarg.h>
#include<math.h>
#include<stdlib.h>
#include<stdio.h>

#include "interface.h"
#include "cut.h"
#include "4_vector.h"
#include "q_kin.h"
#include "regul.h"
#include "runVegas.h"
#include "rw_sess.h"
#include "subproc.h"
#include "vegas.h"
#include "strfun.h"
#include "alphas2.h"
#include "crt_util.h"
#include "histogram.h"
#include "drandXX.h"
#include "n_calchep_.h"
#include "files.h"
#include "usrfun.h"
#include "qcdScale.h"

//========= old events.c

#include"chep_crt.h"
#include"../../include/version.h"

static FILE * events_;


static pthread_mutex_t wrt_key=PTHREAD_MUTEX_INITIALIZER;


#define buffSize 1000

static void writeEvent(double *x,  double  w)
{ 
   int i,icc;
   double GG,qF1,qF2,qR,qS,x1,x2;

   REAL pvectR[4*maxNp],cb_coeff_[buffSize],*cb_coeff;
   double pvect[4*maxNp];
   double factor_0;
   int cb_pow=cb_int[Nsub-1].pow;
   int nc=cb_int[Nsub-1].nC; 
   
   mkmom(x, &factor_0,&x1,&x2,pvectR);
   
   for(i=0;i<4*(nin_int+nout_int);i++) pvect[i]=pvectR[i];   
   Scale(Nsub,pvect,&qR,&qF1,&qF2,&qS);
   
   if(cb_pow)
   {  double sum=0;
      int err;
      if(buffSize>cb_pow) cb_coeff=cb_coeff_; else cb_coeff=malloc(sizeof(REAL)*cb_pow);
      GG=sqrt(4*M_PI*alpha_2(qR));    
      sqme_int(Nsub,GG,pvectR,cb_coeff,&err);
         
      for(i=0;i<cb_pow;i++) sum+=fabs(cb_coeff[i]);
     if(nPROCSS) pthread_mutex_lock(&drandXX_key);
      sum*=drandXX();
     if(nPROCSS) pthread_mutex_unlock(&drandXX_key);
      for(i=0;i<cb_pow;i++)
      { sum-=fabs(cb_coeff[i]);
        if(sum<=0) break;
      } 
      if(i==cb_pow) i--;
      icc=i;
      if(cb_coeff!=cb_coeff_) free(cb_coeff);
   } 

   if(nPROCSS)pthread_mutex_lock(&wrt_key);   
   fprintf(events_,"%12.3E",w);  
   if(nin_int==2) fprintf(events_," %17.10E %17.10E",pvect[3],pvect[7]);

   for(i=0;i<nout_int;i++) fprintf(events_," %17.10E %17.10E %17.10E",
   pvect[4*(i+nin_int)+1],pvect[4*(i+nin_int)+2],pvect[4*(i+nin_int)+3]);

   fprintf(events_,"| %.3E  %.3E ", qS,alpha_2(qR));
 
   if(cb_pow)
   {  int j;
      fprintf(events_,"  ");
      for(j=0;j<nc;j++)
      { int *offset=cb_int[Nsub-1].chains +4*nc*icc+4*j;
        fprintf(events_,"(%d %d %d %d)",offset[0],offset[1],offset[2],offset[3]);
      }  
   }
   fprintf(events_,"\n");
//   fflush(events_);     
   if(nPROCSS)pthread_mutex_unlock(&wrt_key);
}

//static long n_cube=1000;
#define n_cube EventGrid
int  EventGrid=10000;

int regenEvents=0; 
static int nEvents=10000;
static double max =1.2;
static double milk=0.1; 

int saveEventSettings(FILE * f)
{
  fprintf(f, "%d %d %d ",n_cube, regenEvents,  nEvents);
  return 0;
}

int  readEventSettings(FILE * f)
{
  fscanf(f, "%d %d %d ", &n_cube, &regenEvents, &nEvents);
  return 0;
}


static void write_event_cap(void) 
{
  int i,j;

  fprintf(events_,"#%s\n",VERSION_);
  fprintf(events_,"#Type %d -> %d\n",nin_int,nout_int); 
  fprintf(events_,"#Initial_state\n");
  if(nin_int==1) fprintf(events_,"  P1_3=0\n");
  else
  {  fprintf(events_,"   P1_3=%E  P2_3=%E\n" , inP1,-inP2); 
     wrt_sf__(events_);
  }
  fprintf(events_,"#PROCESS  ");
  for(i=1;i<=nin_int+nout_int; i++)
  { int pcode;
    char * pname=pinf_int(Nsub,i,NULL,&pcode);
    fprintf(events_," %d(%s)", pcode, pname);
    if(i==nin_int)  fprintf(events_," ->");
  } 
  fprintf(events_,"\n");    
  fprintf(events_,"#MASSES ");
  for(i=0;i<nin_int+nout_int;i++)
  {  REAL m;
     pinf_int(Nsub,i+1,&m,NULL); 
     fprintf(events_," %.10E", (double)m);
  }   
  fprintf(events_,"\n");
  
  if(integral.n_it) fprintf(events_,"#Cross_section(Width) %E\n",integral.s1/integral.n_it);
  else   fprintf(events_,"#Cross_section(Width) Unknown\n"); 
  fprintf(events_,"#Number_of_events %10d\n",0);
  fprintf(events_,"#Sum_of_weights %12.4E %12.4E \n",0.,0.);
  fprintf(events_,"#Events  "); 
  if(nin_int==2) fprintf(events_,"     P1_3 [Gev]        P2_3 [Gev]   ");
  for(i=1;i<=nout_int; i++) for(j=1;j<=3;j++) 
                          fprintf(events_,"     P%d_%d [Gev]   ",i+nin_int,j);
//  fprintf(events_,"  QCD SCALE    Color  !chains\n");
  fprintf(events_,"  Q_factor   alpha_QCD  Color chains\n");

}

   
static void  generateEvents( vegasGrid * vegPtr,  char *fname,  FILE * iprt)
{ long  nGenerated=0;                                                             
  double eff;                                                                     
  int nmax,mult,neg;                                                              
  char mess[200];                                                                 
  long cEvent;                                                                    
  long fileEnd;
  double sumW,sumW2;                                                                   
  event_stat stat;                                                                
                                                                                  
  events_= fopen(fname,"a");                                                      
  if(ftell(events_)==0) write_event_cap();                                        
  fileEnd=ftell(events_);                                                         
  fflush(events_);                                                                
  cEvent= vegas_events(vegPtr,nEvents,max,writeEvent,regenEvents,nPROCSS,&stat); //   &eff,&nmax,&mult,&neg);       
  eff=stat.eff; neg=stat.neg; mult=stat.lmax;                                     
  fclose(events_);                                                                
                                                                                  
  if(cEvent>0)
  {  int l;                                                                       
     sprintf(mess,"Statistic\n Events generated: %ld\n  efficiency: %.1E\n"       
               "<weight^2>/<weight>^2-1: %.1E \nNegative weight  events: %d \n", cEvent, eff,
                 cEvent*stat.sumW2/stat.sumW/stat.sumW-1., neg);
                                                                                  
     l=strlen(mess);
     strcat(mess,"---------------\n Accept events? ");                            
     if(mess_y_n(25,15,mess))                                                     
     {                                                                            
        long  nEvPos=0,nWpos=0;                                                           
        integral.old=1;                                                           
        mess[l]=0;                                                                
        events_=fopen(fname,"r+");                                                
        while(nEvPos==0 || nWpos==0)                                                          
        { char ch;                                                                
          char word[100];                                                         
          do fscanf(events_,"%c",&ch); while(ch !='#');                           
          fscanf(events_,"%s",word);                                              
          if(strcmp(word,"Number_of_events")==0) nEvPos=ftell(events_);
          if(strcmp(word,"Sum_of_weights")==0) nWpos=ftell(events_);           
        }
        fseek(events_,nEvPos,SEEK_SET);                                                                           
        fscanf(events_,"%ld",&nGenerated);
        nGenerated+=cEvent;
        fseek(events_,nEvPos,SEEK_SET);                                                                                                  
        fprintf(events_," %10ld",nGenerated); 
        
        fseek(events_,nWpos,SEEK_SET);
        fscanf(events_,"%lf %lf",&sumW,&sumW2);
        sumW+=stat.sumW;
        sumW2+=stat.sumW2;
        fseek(events_,nWpos,SEEK_SET);
        fprintf(events_," %12.4E %12.4E",sumW,sumW2); 
                                            
        fclose(events_);                                                          
        fprintf(iprt," %ld events are stored in '%s'\n",nGenerated,fname);        
        fprintf(iprt,"%s\n",mess);                                                
        fflush(iprt);                                                             
     } else  truncate(fname,fileEnd);                                             
  }                                                                               
}                                                                                 

//======== old runVegas.c
int nSess=1;

double inP1=4000, inP2=4000;

static vegasGrid * veg_Ptr=NULL;
static int hFill=0;
vegas_integral integral={{5,5},{100000,100000},0,0.,0.,0.,0.,0.,0.,0,0,0,-1}; 


char * effInfo(void)
{ static char buff[10]; 
  int k;
  double sum; 
  if( !integral.freeze || !integral.n_it || !integral.s1 || !veg_Ptr || !veg_Ptr->fMax  ) { buff[0]=0; return buff;}
  for(sum=0,k=0;k<veg_Ptr->evnCubes;k++) sum+=veg_Ptr->fMax[k];
  sprintf(buff,"%.1E", integral.s1/integral.n_it*veg_Ptr->evnCubes/sum/1.2);
  return buff; 
}

static void clearStatistics(int tp)
{
  integral.In=0;
  integral.dI=0;
  integral.khi2=0;
  integral.s0=0; 
  integral.s1=0; 
  integral.s2=0; 
  integral.n_it=0; 
  integral.nCallTot=0; 
  integral.tp=tp;
  clearHists();  
  { char fname[20];
    sprintf(fname,"distr_%d",nSess);
    unlink(fname);
  }
}                         

void clearGrid(void){ vegas_finish(veg_Ptr); veg_Ptr=NULL;}
void clearEventMax(void)
 { if(veg_Ptr && veg_Ptr->fMax) {free(veg_Ptr->fMax); veg_Ptr->fMax=NULL; }}

void newSession(void)
{
   if(integral.old)
   { char fname[20];

     messanykey(15,15,
     "Some parameters where changed.\nSo integral and statictics for\n"
     "distribushions is forgotten!\nSession number is increased.");
     integral.old=0;
     nSess++;
     clearStatistics(-1);
     sprintf(fname,"prt_%d",nSess);
     unlink(fname);
   }
}



int saveVegasGrid( FILE * f)
{
  if(veg_Ptr)
  {  int i,j;
     fprintf(f," Vegas_grid: dim=%d  size=%d\n", veg_Ptr->dim, veg_Ptr->ndmx);
     for(i=0;i<veg_Ptr->dim;i++)
     { for(j=0;j<=veg_Ptr->ndmx;j++) fprintf(f," %.15E",veg_Ptr->x_grid[i][j]);
       fprintf(f,"\n");
     }
     if(veg_Ptr->fMax)
     { long l;
       fprintf(f,"Max(%ld):\n",veg_Ptr->evnCubes);
       for(l=0;l<veg_Ptr->evnCubes;l++) fprintf(f,"%.1E\n",veg_Ptr->fMax[l]);
     } else fprintf(f,"Max(0):\n");
  }else  fprintf(f," Vegas_grid: dim=%d  size=%d\n", 0, 0);
  return 0;
}

static double func_(double *x, double wgt);

int readVegasGrid(FILE * f)
{
  int i,j,ndim,ndmx;
  int nCubes;
  
  if(veg_Ptr) {vegas_finish(veg_Ptr);veg_Ptr=NULL;}  
  fscanf(f," Vegas_grid: dim=%d  size=%d\n", &ndim, &ndmx);
  if(ndim && ndmx)
  { 
    veg_Ptr=vegas_init(ndim,func_,ndmx);
    for(i=0;i<ndim;i++)for(j=0;j<=ndmx;j++) fscanf(f," %lf",&(veg_Ptr->x_grid[i][j]) );
    fscanf(f," Max(%d):\n",&nCubes);   
    
    setEventCubes(veg_Ptr, nCubes);
    if(nCubes && veg_Ptr->evnCubes==nCubes) 
    { long l;
      veg_Ptr->fMax=malloc(sizeof(double)*veg_Ptr->evnCubes);
      for(l=0;l<veg_Ptr->evnCubes;l++) fscanf(f,"%lf",veg_Ptr->fMax+l);
    }
  }
  return 0;
}


static double badPoints;
static double negPoints;

static void printLn(FILE * iprt,int *line,char * format, ...)
{  
   va_list args;
   char dump[STRSIZ];
   va_start(args, format);
   vsprintf(dump,format,args);
   va_end(args);

   goto_xy(1,*line); print("%53s","");
   goto_xy(1,*line);
   print("%s\n",dump); 
   (*line)++;
   if (*line >= maxRow()-2 ) *line=8; else 
   {
     scrcolor(Blue, BGmain);
     print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
   }   
   if(iprt) {fprintf(iprt,"%s\n",dump); fflush(iprt);}
}

static pthread_mutex_t hist_key=PTHREAD_MUTEX_INITIALIZER;

static double func_(double *x, double wgt)
{
    double ret_val=0.;
    int err=0,nd,i;
    double factor_0;
    double x1,x2;    
    double GG,qF1,qF2,qR,qS;
    REAL pvectR[100];
    double pvect[100];
/* get momenta */
    mkmom(x, &factor_0,&x1,&x2,pvectR);

    if (!factor_0) goto exi;
    nd=4*(nin_int+nout_int);
    for(i=0;i<nd;i++) pvect[i]=pvectR[i];
    factor_0 *= calcCutFactor(pvect)*usrFF(nin_int,nout_int,pvect,p_names,p_codes); 
    if (!factor_0)   goto exi;

    Scale(Nsub,pvect,&qR,&qF1,&qF2,&qS);
/* **  structure function  multiplication */
    if (nin_int == 2) 
    {
	if(sf_num[0]) { factor_0 *= strfun_(1, x1,qF1);  if(factor_0==0.) goto exi;}
	if(sf_num[1]) { factor_0 *= strfun_(2, x2,qF2);  if(factor_0==0.) goto exi;} 
    }   
    if (!factor_0)  { goto exi;}
/* ** call for 'running strong coupling constant' */
    GG=sqrt(4*M_PI*alpha_2(qR));    
    ret_val = factor_0 * sqme_int(Nsub,GG,pvectR,NULL,&err);

    if(err)       badPoints+=  (ret_val>0 ? ret_val*wgt : - ret_val*wgt); 
    if(ret_val<0) negPoints+=ret_val*wgt;
exi:
    if(hFill)
    { if(nPROCSS) pthread_mutex_lock(&hist_key);
       fillHists(ret_val*wgt,pvect); 
      if(nPROCSS) pthread_mutex_unlock(&hist_key); 
    } 
    return ret_val;
} /* func_ */




int runVegas(void)
{
    int i;
    double sd;
    double avgi;
    char mess[25];
    FILE * iprt = NULL;
    int mode=1;
    void * pscr=NULL;
    static int n_Line=0;
    
    if(blind) vegas_control=NULL; else  vegas_control=infoLine;
    i=imkmom(inP1,inP2);
    if(veg_Ptr&&veg_Ptr->dim!=i)clearGrid();
    if(!veg_Ptr) veg_Ptr=vegas_init(i,func_,50);     

    if(nin_int == 2) strcpy(mess, "Cross section[pb]");
      else           strcpy(mess, "   Width[Gev]    ");
    
/* ** save current session parameters */
     w_sess__(NULL);
/* ** open protocol and resulting files */
       
    {  char fname[50];
       sprintf(fname,"%sprt_%d",outputDir,nSess);
       iprt=fopen(fname,"a");
    }

/* **  initkinematics */

    while(correctHistList()) editHist();

    if(integral.old) { if(!n_Line)  n_Line=9;} else
    {
       w_sess__(iprt);
       n_Line=8;
    }
                           
/* *** Main cycle */

    for(;;)
    {
        int fz;
        char strmen[]="\030"
         " nSess  = N2_1          "
         " nCalls = N1_1          "
         " Set  Distributions     "
         "*Start integration      "
         " Display Distributions  "
         " Clear statistic        "
         " Freeze grid        OFF " 
	 " Clear  grid            "
	 " Event Cubes NCUBE      "
	 " Num. of events=NE      "
	 " Generate Events        ";
	 
        if(integral.freeze)
        {  fz=1;
           improveStr(strmen,"OFF","ON");
           improveStr(strmen,"NCUBE","%d",EventGrid);
        }  
        else
        {  fz=0;
           strmen[ 030*8+2]=0;
        }

        improveStr(strmen,"N1_1","%d",integral.ncall[fz]);
        improveStr(strmen,"N2_1","%d",integral.itmx[fz]);
        improveStr(strmen,"NE","%d",nEvents);
        menu1(54,7,"",strmen,"n_veg_*",&pscr,&mode);
        switch(mode)
        {     
        case 0:
           w_sess__(NULL);
          if(iprt) fclose(iprt);
          return 0;           
        case 1: correctInt(50,12,"Enter new value ", integral.itmx+fz,1);  break;
        case 2: correctLong(50,12,"Enter new value ",integral.ncall+fz,1);break;
        case 3:  editHist(); break;
        case 4:
          if(veg_Ptr->fMax && !integral.freeze) 
          {  free(veg_Ptr->fMax); veg_Ptr->fMax=NULL; veg_Ptr->evnCubes=0; }
          if(!veg_Ptr->fMax && integral.freeze)
          {  setEventCubes(veg_Ptr, EventGrid);
             EventGrid=veg_Ptr->evnCubes;
          }

          for (i = 1; i <= integral.itmx[fz]; ++i)                                       
          { char  errtxt[100]="";
            long nCall;
            if(integral.ncall[0]==0) break;                                                                  
            negPoints=0;                                                              
            badPoints=0; 
            hFill=1;   
            nCall=vegas_int(veg_Ptr, integral.ncall[fz],1.5*(!fz),nPROCSS,&avgi, &sd);
            if(nCall<0) { messanykey(10,10,"NaN in integrand"); break;}
            if(nCall==0) break;
            
            integral.old=1;                                              
            negPoints/=nCall;                                                         
            badPoints/=nCall;                                                         
            integral.nCallTot+=nCall;                                                          
            scrcolor(FGmain,BGmain);                                                 
            printLn(iprt,&n_Line,"%4d   %12.4E %10.2E %8d %s",                     
                 ++integral.n_it, avgi,avgi? 100*sd/(double)fabs(avgi):0.,nCall,effInfo());
                                                                   
            if(negPoints<0) sprintf(errtxt+strlen(errtxt)," Negative points %.1G%%;",                
                                      -100*negPoints/(avgi-2*negPoints));             
            if(badPoints)  sprintf(errtxt+strlen(errtxt),                             
                 "Bad Precision %.1G%%;",100*badPoints/(avgi-2*negPoints));           
                                                                                      
            if(errtxt[0])                                                             
            {                                                                         
               scrcolor(Red,BGmain);                                                  
               printLn(iprt,&n_Line,"%s",errtxt);                                     
            }
                                                                                                                                                
            integral.s0+=sd*sd;                                                                  
            integral.s1+=avgi;                                                             
            integral.s2+=avgi*avgi;                                    
          } 
          
          
          integral.In=integral.s1/integral.n_it; 
          integral.dI=sqrt(integral.s0)/integral.n_it;
          if(integral.n_it<=1 || integral.s0==0 ) integral.khi2=0; else 
          integral.khi2=(integral.s2-integral.s1*integral.s1/integral.n_it)*integral.n_it/(integral.n_it-1)/fabs(integral.s0);  
          
          scrcolor(FGmain,BGmain);

          printLn(iprt,&n_Line," < >   %12.4E %10.2E %8d %7.7s %-7.1G" ,
                      integral.In, fabs(integral.In)? 100*integral.dI/(double)fabs(integral.In):0., integral.nCallTot, 
                                                              effInfo(),  integral.khi2);
          if(histTab.strings)
          { char  fname[20];
            FILE * d;
            sprintf(fname,"distr_%d",nSess);
            d=fopen(fname,"w");  
            wrt_hist2(d,Process);
            fclose(d);
          }
                    messanykey(54,11,"Integration is over");
/*          integral.freeze=0; */
          break;

        case 5: showHist(54,10,Process); break;
        case 6: clearStatistics(-1);
                messanykey(54,13,"Old results for integral\n"
                "and distributions\nare deleted.");
                break;
        case 7: 
             if(veg_Ptr->fMax && integral.freeze) 
             {  if(mess_y_n(15,15,"You have event generator prepared.\n"
                " Setting the flag \"OFF\"  will destroy it."
                " Press 'Y' to confirm.")) integral.freeze=0; 
             } else  integral.freeze=!integral.freeze; 
             break; 
        case 8: if(!integral.freeze || mess_y_n(15,15,"The information for Event Generator will be lost\n OK?"))  
                { int ndim=veg_Ptr->dim;
                  vegas_finish(veg_Ptr);
                  veg_Ptr=vegas_init(ndim,func_,50);
                  messanykey(57,11,"OK");
                }   
                break;
        case 9: 
           if(correctInt(50,12,"Enter new value ",&EventGrid,1))
           { if(veg_Ptr->fMax) {free(veg_Ptr->fMax); veg_Ptr->fMax=NULL;}
             printf("EventGrid=%d\n",EventGrid);
             setEventCubes(veg_Ptr, EventGrid);
             EventGrid=veg_Ptr->evnCubes;  
           } break;
        case 10:  correctInt(50,15,"",&nEvents,1); break;
        case 11: 
           if( !veg_Ptr || !veg_Ptr->fMax)
           { char * mess="Before event generation one has to launch  Vegas session with freezed grid\n"
                                           "to prepare generator";
                messanykeyErr(4,13,mess);
           } else  runEvents();
       }
    }    
}


int runEvents(void)
{
    FILE * iprt = NULL;
    int i;

    i=imkmom(inP1,inP2);
//    if(veg_Ptr&&veg_Ptr->ndim!=i)clearGrid();
//    if(!veg_Ptr) veg_Ptr=vegas_init(i,50);    

    w_sess__(NULL);
/* ** open protocol and resulting files */
       
    {  char fname[50];
       sprintf(fname,"%sprt_%d",outputDir,nSess);
       iprt=fopen(fname,"a");
       if(ftell(iprt)==0) 
       { fprintf(iprt,"    CalcHEP kinematics module \n The session parameters:\n");
         w_sess__(iprt);
         fprintf(iprt,"===================================\n");   
       }
    }

/* **  initkinematics */


    { char fname[50];

      sprintf(fname,"%sevents_%d.txt",outputDir,nSess);

      hFill=0;
      generateEvents(veg_Ptr,fname, iprt);
    }
    fclose(iprt);
    return 0;
}
