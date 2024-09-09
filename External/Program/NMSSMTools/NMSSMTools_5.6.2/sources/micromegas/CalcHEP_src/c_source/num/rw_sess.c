/*
 Author  Alexander Pukhov
*/
#include"interface.h"

#include"cut.h"
#include"comp.h"
#include"kinaux.h"
#include"kininpt.h"

#include"regul.h"
#include"rw_sess.h"
#include"alphas2.h"
#include"histogram.h"
#include"crt_util.h"
#include"n_calchep_.h"
#include"runVegas.h"
#include"drandXX.h"
#include"interface.h"

#include"strfun.h"
#include"files.h"
#include"subproc.h"
#include"paragraphs.h"
#include"num_in.h"
#include"VandP.h"
#include"dynamic_cs.h"

#define NITEMS 20
char *p_names[20];
int   p_codes[20];
static int rdErr;

static int w_SLHA(FILE *f) { fprintf(f," %d",useSLHAwidth);  return 0; }
static int r_SLHA(FILE *f) { fscanf(f," %d", &useSLHAwidth); return 0;  }

static int writeIntegral(FILE *f)
{
  fprintf(f," %.17E %.17E %.17E %d %d %ld %d %d", integral.s0,integral.s1,integral.s2,
  integral.n_it, integral.old,integral.nCallTot,integral.freeze,integral.tp);
  return 0;
}

static int readIntegral(FILE *f)
{
  fscanf(f," %lf %lf %lf %d %d %ld %d %d", &(integral.s0),&(integral.s1),
  &(integral.s2),&(integral.n_it),&(integral.old),&(integral.nCallTot),&(integral.freeze),&(integral.tp));
  return 0;
}

static int wnsess_(FILE *mode) { fprintf(mode,"%d",nSess); return 0; }
static int rnsess_(FILE *mode) { fscanf(mode,"%d",&nSess); return 0; }

static int w_prc__(FILE *mode)
{ return  fprintf(mode,"%d ( %s )",Nsub,Process); }

static int r_prc__(FILE *mode)
{  char proc_name[100];
   fscanf(mode,"%d ( %[^)]",&Nsub,proc_name);
   if(Nsub> nprc_int) return 2;
   trim(proc_name);
   wrtprc_();
   if (strcmp(proc_name,Process)==0) return 0; else return 2;
}

static int w_mc__(FILE *mode)
{ fprintf(mode,"%ldx%d",integral.ncall[0],integral.itmx[0]);
  fprintf(mode," %ldx%d",integral.ncall[1],integral.itmx[1]);
  return 0;
}

static int r_mc__(FILE *mode)
{ fscanf(mode,"%ldx%d",integral.ncall,integral.itmx); 
  fscanf(mode," %ldx%d",integral.ncall+1,integral.itmx+1);
  return 0;
}

static int  wrtNcores(FILE *mode)
{ 
  fprintf(mode,"  %d\n",nPROCSS);   
  return 0;
}

static int rdNcores(FILE *mode)
{ 
  fscanf(mode," %d",&nPROCSS);
  return 0;
}




static int w_widths(FILE *mode)
{ 
   fprintf(mode,"BW range    %f\n"
                "t-channel widths %d\n"
                "GI trick in s-  %d\n"
                "GI trick in t-  %d\n"
                "useSLHAwidth    %d\n",
           *BWrange_int,*twidth_int,*gswidth_int,*gtwidth_int,useSLHAwidth);
   return 0;             
}

static int r_widths(FILE *mode)
{    fscanf(mode,"BW range    %lf\n"     
                 "t-channel widths %d\n"
                 "GI trick in s-  %d\n"                  
                 "GI trick in t-  %d\n"
                 "useSLHAwidth    %d\n",   
           BWrange_int,twidth_int,gswidth_int,gtwidth_int,&useSLHAwidth);

     return 0;
}                                   


static int w_VVdecays(FILE *mode) { fprintf(mode," %d\n", VWdecay); return 0;  }

static int r_VVdecays(FILE*mode)  { fscanf(mode," %d",&VWdecay);  VZdecay=VWdecay; return 0; }

static int w_mdl__(FILE * mode)
{
  int i;

  fprintf(mode,"\n"); 
  for (i = 0; i < nModelVars; ++i) fprintf(mode,"%10s = %.15E\n",varNames[i],(double)varValues[i]);
  fprintf(mode,"----\n"); 
  return 0;
}

static int r_mdl__(FILE * mode)
{
  static double val;
  int i;
  char name1[20];    

  for(;;)
  { if(2!= fscanf(mode,"%s = %lf",name1,&val)) break;
    for(i=0;i< nModelVars;i++)if(strcmp(name1,varNames[i])==0) { varValues[i] = val; break;}
    if(i==nModelVars && strcmp(name1,"GG")!=0  )
    { char mess[100]; sprintf(mess,"session.dat contains variable '%s', which is absent in the model\n",name1);
      messanykeyErr(10,10,mess);
    }
  }  
  return 0;
} 

static double hMem[2]={0,0};

static int w_in__(FILE *mode)
{ 
  fprintf(mode," inP1=%E  inP2=%E\n",inP1,inP2);

  if(nin_int==2) {if(is_polarized(1,Nsub)) hMem[0]=Helicity[0]; if(is_polarized(2,Nsub))hMem[1]=Helicity[1];}
  
  fprintf(mode," Polarizations= { %E  %E }\n",hMem[0],hMem[1]);
  
  wrt_sf__(mode);

  return 0;
}

static int r_in__(FILE *mode)
{ 
   fscanf(mode," inP1=%lf  inP2=%lf\n",&inP1,&inP2); 

   fscanf(mode," Polarizations= { %lf  %lf }\n",hMem,hMem+1);
            
   if(nin_int==2)
   {
      if(is_polarized(1,Nsub)) Helicity[0]=hMem[0];
      if(is_polarized(2,Nsub)) Helicity[1]=hMem[1];
   }   
   rd_sf__(mode);
   return 0;
}


static int saveRandom(FILE *f)
{ fprintf(f,"%s\n",seedXX(NULL)); return 0; }
static int readRandom(FILE *f)
{ char s[20]; fscanf(f,"%s",s); seedXX(s); return 0;}


int is_polarized(int k,int  nsub)
{ char name[20];
  if(nin_int==1 || k>2 || k<1 ) return 0;
  sprintf(name,",%s,",pinf_int(nsub,k,NULL,NULL));
  if(strstr(polarized_int[k],name)) return 1; else return 0;
}

int wrtprc_(void) /* write process string */
{
  int i;
  char * s=Process;
  for(i=0; i<nin_int+nout_int; i++)p_names[i]=pinf_int(Nsub,i+1,NULL,p_codes+i);

  strcpy(s,p_names[0]);
  if(is_polarized(1,Nsub)) strcat(s,"%");
  if(nin_int==2) 
  { strcat(s,", ");
    strcat(s,p_names[1]);
    if(is_polarized(2,Nsub)) strcat(s,"%");
  }
  strcat(s," -> ");
  for (i=nin_int;i<nin_int+nout_int ;i++)
  {  if( i!=nin_int)  strcat(s,", ");
     strcat(s,p_names[i]);
  }   
  return 0;
}


int w_sess__(FILE *mode_)
{
   FILE*mode;
   rw_paragraph  wrt_array[NITEMS]=
   {
      {"Subprocess",      w_prc__},
      {"Session_number",  wnsess_},
      {"Initial_state",   w_in__},
      {"Physical_Parameters",  w_mdl__},
      {"Breit-Wigner",    w_widths}, 
      {"VVdecays",        w_VVdecays},
      {"SLHAwidth",       w_SLHA},
      {"alphaQCD",        w_alphaQCD},
      {"QCDscales",       w_Scales},      
      {"Composites",      wrtcomp_},
      {"Cuts",            wrtcut_},
      {"Parallelization", wrtNcores},
      {"Distributions",   wrt_hist},
      {"Kinematical_scheme",  wrtkin_},
      {"Regularization",  wrtreg_},
      {"Vegas_calls",     w_mc__},
      {"Vegas_integral",  writeIntegral},
      {"Events",          saveEventSettings},
      {"Random",          saveRandom},
      {"VEGAS_Grid",      saveVegasGrid}
   };

   if (mode_ == NULL)
   { char fname[100];
     sprintf( fname,"%ssession.dat",outputDir);
       mode=fopen(fname,"w");
      if(mode == NULL) return 0;
   } else mode=mode_;

   if(mode_ == NULL)
   {  char fname[100];
      writeParagraphs(mode,NITEMS,wrt_array);
      fclose(mode);
      sprintf(fname,"%saux/session.dat",outputDir);
      mode=fopen(fname,"w");
      if(mode) {writeParagraphs(mode,9,wrt_array+2); fclose(mode);}
   } else   writeParagraphs(mode,11,wrt_array);

   return 0;
}    

int r_sess__(FILE *mode_)
{
  FILE*mode;
  rdErr=0;
 rw_paragraph  rd_array[NITEMS]=
 {
   {"Subprocess",      r_prc__},
   {"Session_number",  rnsess_},
   {"Initial_state",   r_in__},
   {"Physical_Parameters",  r_mdl__},
   {"Breit-Wigner",    r_widths},
   {"VVdecays",        r_VVdecays}, 
   {"SLHAwidth",       r_SLHA},
   {"alphaQCD",        r_alphaQCD},
   {"QCDscales",       r_Scales},          
   {"Composites",      rdrcomp_},
   {"Cuts",            rdrcut_},
   {"Parallelization", rdNcores}, 
   {"Distributions",   rdr_hist},
   {"Regularization",  rdrreg_},
   {"Kinematical_scheme",  rdrkin_},
   {"Vegas_integral",  readIntegral},
   {"Vegas_calls",     r_mc__},  
   {"Events",          readEventSettings},
   {"Random",          readRandom},
   {"VEGAS_Grid",      readVegasGrid}
 };
                       
                         
 if (mode_ == NULL) 
 { 
    mode=fopen("session.dat","r"); 
    if (mode ==NULL)
    { mode=fopen("aux/session.dat","r");
      if(mode==NULL) return -1;
    }
 }else mode=mode_;
 readParagraphs(mode, NITEMS,rd_array);  
 if (mode_ == NULL) fclose(mode); 
 wrtprc_();
 return rdErr;;
}
