/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>
#include <pthread.h>

#include <dlfcn.h>

#include "files.h"
#include "interface.h"
#include "subproc.h"
#include "chep_crt.h"
#include "strfun.h"
//#include "n_calchep_.h"
#include "lha.h"
#include "alphas2.h"
#include "../../../include/VandP.h"
#include "n_proc.h"
#include "sf_lha.h"


static int VERS=0;

static double xMin[2]={0,0},xMax[2]={0,0},qMin[2]={1,1},qMax[2]={1E10,1E10};

/*static char * pdfName[2]={NULL,NULL}; */

static int nGroup[2]={3,3}, nSet[2]={41,41}, sgn[2]={1,1},lastLHA=0;

static int parton[12]={2,1,4,3,5,21,-5,-3,-4,-1,-2,22};

static int pnum[2]={0,0};


static char *param=NULL;

static double sc2=0.22*0.22;

static char* dataPath=NULL;
static void (*getpdfsetlist)(char* s, size_t len)=NULL;
static char  fname[2][30]={"",""};
static int   setNum[2]={0,0};


int mc_lha(int i) { return 2212*sgn[i-1];}

static pthread_mutex_t alpha_key=PTHREAD_MUTEX_INITIALIZER;

static double alpha_lha1(double q ) 
{ double al;
  int one=1;
   if(nPROCSS>1 && VERS<6) pthread_mutex_lock(&alpha_key);                                  
   al=alphaspdfm(&one,&q);
   if(nPROCSS>1 && VERS<6) pthread_mutex_unlock(&alpha_key);
  return al; 
}

static double alpha_lha2(double q )
{ double al;
  int two=2;
  if(nPROCSS>1 && VERS<6) pthread_mutex_lock(&alpha_key);
  al=alphaspdfm(&two,&q);    
  if(nPROCSS>1 && VERS<6) pthread_mutex_unlock(&alpha_key);
  return al;
}




static int initLHA(void)
{ 
  int i;  
  if(VERS<0) return 0;
  if(VERS>0) return 1;  

  char versionTxt[40];
  getlhapdfversion(versionTxt,39);
  if(1!=sscanf(versionTxt,"%d",&VERS)) { VERS=-1; return 0;}
  if(VERS<5) { VERS=-1; return 0;}
  
  if(VERS==5)
  {  char buff[500]; 
     getdatapath(buff,500);  
     for(i=499; i>=0 && buff[i]==' '; i--) ;
     if(i<0) { VERS=-1;  return 0;}
     buff[i+1]=0;
     dataPath=malloc(strlen(buff)+1);
     strcpy(dataPath,buff);
  } else
  { 
  
     void*h=dlopen(NULL,RTLD_NOW); 
     getpdfsetlist=dlsym(h,"lhapdf_getpdfsetlist_");
     if(!getpdfsetlist) getpdfsetlist=dlsym(h,"_lhapdf_getpdfsetlist_");
     if(!getpdfsetlist) {VERS=-1; return 0;}
  }
  return 1;
}

int p_lha(int * pNum) 
{  
  int i;
  if(!initLHA())return 0;
  for(;*pNum;pNum++) 
  { if(*pNum==22) 
    { if(VERS==5) return 0;
//      if(!has_photon()) return 0; 
    }
    for(i=0;i<12;i++) if(*pNum==parton[i]) break;
    if(i==12) return 0;
  }
  return 1;
}

void n_lha(int i, char *name) 
{  int i1=i-1;

   if(strlen(fname[i1]))  sprintf(name,"LHA:%s:%d:%d",fname[i1],setNum[i1],sgn[i1]);
   else  strcpy(name,"LHA:"); 
}



int init_lha(int i,double * be, double * mass) 
{  int k;
   int N,N1,N2;
   pinf_int( Nsub,1,NULL,&N1);
   pinf_int( Nsub,2,NULL,&N2);

   if(i==1) N=N1; else N=N2;
    
   *mass=0.9383;
   *be=1;
   if(i==1) sf_alpha[0]=alpha_lha1; else sf_alpha[1]=alpha_lha2;   
   for(k=0;k<12;k++) if(parton[k]==N) {pnum[i-1]=N; return 1;}
   pnum[i-1]=0; 
   return 0;
}
 
static int filter (const struct dirent * dp)
{  
   char *ch ;
   
/*   if(dp->d_type == DT_DIR) return 0; */

   ch=strstr(dp->d_name,".LHpdf");
   if(ch && ch[6]==0) return 1;

   ch=strstr(dp->d_name,".LHgrid");
   if(ch && ch[7]==0) return 1;
   
   return 0;
}


static char * lhaMenu(void)
{
  struct  dirent **namelist;
  int N,i,width;
  char * menutxt=NULL;

  if(getpdfsetlist)
  {  char buff[10000];
     char *ch1,*ch2;
     getpdfsetlist(buff,10000);
     for(i=9999; i>=0 && buff[i]==' '; i--) ;
     buff[i+1]=0;
     i=strlen(buff);
     width=0; N=0;
     ch1=ch2=buff;
     if(ch1[0]==0) return NULL;
     for(;ch2!=buff+i; ch1=ch2+1,N++)
     {  ch2=strchr(ch1,' ');
        if(ch2==NULL) ch2=buff+i;
        if(ch2-ch1>width) width=ch2-ch1;   
     }
     
     menutxt=malloc(N*(width+2)+2);
     N=0;
     menutxt[0]=width+2; 
     for(ch1=buff;ch1; ch1=strchr(ch1,' '),N++)
    {  char name[50];
       while(ch1[0]==' ') ch1++;
       sscanf(ch1,"%s",name);
       sprintf(menutxt+1+N*(width+2)," %-*.*s ",width,width,name);
    }
    return menutxt;  
  } 
  
  
  N=scandir(dataPath, &namelist,filter ,alphasort);   
  if(N<=0) return NULL;
                              
  for(i=0,width=0;i<N;i++) 
  {  int l= strlen(namelist[i]->d_name);
     if(width<l) width=l;
  }
  menutxt=malloc(N*(width+1)+2);
  menutxt[0]=width+1; menutxt[1]=0;
  for(i=0;i<N;i++)
  {
    sprintf(menutxt+1+(width+1)*i," %-*.*s",width,width,namelist[i]->d_name);
    free(namelist[i]);
  }
  free(namelist);
  menutxt[N*(width+1)+1]=0;
  return menutxt;
}  


int r_lha(int i, char *name)
{ int i1=i-1;
  char txt[50];
  char*men,*c;
  int max;
  if(!initLHA())return 0; 
  if(3!=sscanf(name,"LHA:%[^:]:%d:%d",txt+1,setNum+i1,sgn+i1)) return 0;  
  if(abs(sgn[i1])!=1) return 0;
  men=lhaMenu();

  if(!men) return 0;
  men[0]=' ';

  txt[0]=' ';
  txt[strlen(txt)+1]=0;
  txt[strlen(txt)]=' ';
  c=strstr(men,txt);
  trim(txt);
  free(men);
  if(!c) return 0;
  initpdfsetbynamem(&i,txt,strlen(txt));
  if(i==0) return 0;
  lastLHA=i;
  numberpdfm(&i,&max);
  if(max<setNum[i1] || setNum[i1]<0) return 0;
  initpdfm(&i,setNum+i1,xMin+i1,xMax+i1,qMin+i1,qMax+i1);
  strcpy(fname[i1],txt);
  return 1;
}

                        
int m_lha(int i,int*pString)
{ 
  int size=2+8*500;
  int size_=1;
  void *pscr=NULL;
  void *pscr0=NULL;
  static int n1=0;
  char * strmen=lhaMenu(); 
  int i1=i-1;

  if(!strmen) return 0;
  

  int n0=1,k,l;

  if(n1==0 && strlen(fname[i1]))
  { char *ch=strstr(strmen,fname[i1]);
    if(ch) n1= 1+(ch-strmen)/strmen[0];
  }
  menu1(5,10,"LHAlib menu",strmen,"",&pscr,&n1);
  if(n1)
  { char buff[50];
    sscanf(strmen+1+strmen[0]*(n1-1),"%s",buff);
    if(strcmp(buff,fname[i1]))
    { strcpy(fname[i1],buff);
      setNum[i1]=0;
      initpdfsetbynamem(&i,buff,strlen(buff));
      lastLHA=i;
      initpdfm(&i,setNum+i1,xMin+i1,xMax+i1,qMin+i1,qMax+i1);
    }
  }  
  else 
  { fname[i1][0]=0;
    setNum[i1]=0;
    sgn[i1]=1;
    return 0;
  } 


  for(;n0!=0 && n0!=3;)
  {  char buff[50];
     int nMax;
     char strmen0[]="\030"
                    " Set = 0                "   
                    " Proton                 "
                    " OK                     ";

     numberpdfm(&i,&nMax);
     if(nMax>1) improveStr(strmen0,"Set = 0","Set = %d [0,%d]",setNum[i1],nMax);
     else       improveStr(strmen0,"Set = 0","Set = 0 (only)");
     
     if(sgn[i1]<0) improveStr(strmen0,"Proton","%s","antiProton");
     
     menu1(5,10,"",strmen0,"",&pscr0,&n0);
     switch(n0) 
     { 
       case 1: if(nMax>1) 
               { correctInt(50,12,"Enter new value ",setNum+i1,1);
                 if(setNum[i1]<0) setNum[i1]=0;
                 if(setNum[i1]>nMax) setNum[i1]=nMax;
                 initpdfm(&i,setNum+i1,xMin+i1,xMax+i1,qMin+i1,qMax+i1);
               }   
               break;
       case 2: sgn[i1]=-sgn[i1]; break;
       case 3: put_text(&pscr0); break;
     }
  }
  
  put_text(&pscr); 
  free(strmen);
  return 1;
}

static pthread_mutex_t strfun_key=PTHREAD_MUTEX_INITIALIZER;

             



double c_lha(int i, double x, double q)
{
  double f[14],fph;
  int p;
  int i1=i-1;
  double z;
  
  p=pnum[i1];

  if(x<xMin[i1]) x=xMin[i1]; else if(x>xMax[i1]) x=xMax[i1];
//  if(q<qMin[i1]) q=qMin[i1]; else 
                        if(q>qMax[i1]) q=qMax[i1];
  if(nPROCSS>1 && VERS<6) pthread_mutex_lock(&strfun_key);
   evolvePDFm(i,x,q,f);
  if(nPROCSS>1 && VERS<6) pthread_mutex_unlock(&strfun_key);  

  if(sgn[i1]<0) p=-p;

  switch(p)
  {
    case 2 :           z=f[8]/x;  break;
    case 1 :           z=f[7]/x;  break;
    case 3 : case -3 : z=f[9]/x;  break;
    case 4 : case -4 : z=f[10]/x; break;
    case 5 : case -5 : z=f[11]/x; break; 
    case 21: case -21: z=f[6]/x;  break;
    case -1:           z=f[5]/x;  break;
    case -2:           z=f[4]/x;  break;
    case 22: case -22: z=f[13]/x; break;
  } 
  if(z<0) return 0;
//  if(z<=0) printf("x=%E q=%E z=%E sc2=%E   \n",x,q,z,sc2); 
  return z;  
}
