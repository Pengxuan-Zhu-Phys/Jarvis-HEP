#include <stdlib.h>

#ifdef __hpux
#include<dl.h>
#else
#include <dlfcn.h>
#endif

#include <unistd.h>
#include <ctype.h>
#include <string.h>

#include <sys/utsname.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "fcompare.h"
#include "../../include/VandP.h"
#include "vp.h"
#include "n_proc.h"
#include "dynamic_cs.h"

#include"SLHAplus.h"

extern  char * trim(char *);
  
char  * libDir=NULL;
char  * modelDir=NULL;
char  * compDir=NULL;
char  * calchepDir=NULL;
int   modelNum=0;
double BWrange=2.7;

int  prepareWorkPlace(void)
{  char * command;
   struct stat buf;
   int err,len,mknew;

   if(!compDir) return  -1;   
   mknew=stat(compDir,&buf);
   len=strlen(compDir)+500;
   if(modelDir) len+=strlen(modelDir);
   command=malloc(len);  

   if(mknew) 
   { char * dir[3]={"tmp","results","models"};
     int i;
     if(mkdir(compDir, 00755)) return -3; 
     for(i=0;i<3;i++) 
     { 
        sprintf(command,"%s/%s",compDir,dir[i]);
        mkdir(command,00755);
     } 
     if(modelDir && modelNum)
     {
       sprintf(command,
       "for FILE in vars func prtcls lgrng extlib\n do\n"
       "  cp %s/\"$FILE\"%d.mdl %s/models/\"$FILE\"1.mdl\n" 
       "done\n", modelDir,modelNum,  compDir);  
       system(command);
     } else { free(command); return -2;} 
   } else 
   {   
     sprintf(command, "dName=%s\n"
     "for FILE in $dName/tmp/* $dName/results/*\n"
     "do\n"
     " if(test ! -d $FILE) then\n"
     "   rm -f $FILE\n"
     " fi\n" 
     "done\n",compDir);     
     system(command);
   } 
   free(command);
   return mknew;
}


//#define SAVE
int cleanWorkPlace(void)
{ 
#ifndef SAVE
   char * command;
   command=malloc(strlen(compDir)+100);
   sprintf(command,"rm -fr %s",compDir);
   system(command);
   free(command);
#endif    
   return 0; 
}


int  checkWorkPlace(void)
{
  char * n1=malloc(strlen(modelDir)+50);
  char * n2=malloc(strlen(compDir)+50);
  char *fList[5]= {"vars","prtcls","extlib","func","lgrng",};
  int i;
  for(i=0;i<5;i++)
  { sprintf(n1,"%s/%s%d.mdl",modelDir,fList[i],modelNum);
    sprintf(n2,"%s/models/%s1.mdl",compDir,fList[i]);
    if(fcompare(n1,n2)) break;
  }   
  free(n1);
  free(n2);
  if(i==5) return 0;
  if(modelDir && modelNum)
  { 
     char* command=malloc(strlen(modelDir)+strlen(compDir)+200);
     sprintf(command,
     "for FILE in vars func prtcls lgrng extlib\n do\n"
     "  cp %s/\"$FILE\"%d.mdl %s/models/\"$FILE\"1.mdl\n" 
     "done\n", modelDir,modelNum,  compDir);  
     system(command);
     free(command);
     delAllLib();
  }
  return 1;
}  


int  checkMtime(char * fname)
{ int i,L;
  time_t tt;
  struct stat buff;
  char * mf[4]={"vars","func","prtcls","lgrng"};
  char *mfname;
  if(modelDir==NULL) return 0;
  stat(fname,&buff); tt=buff.st_mtime;
  L=strlen(modelDir)+20;
  mfname=malloc(strlen(modelDir)+20); 
  
  for(i=0;i<4;i++)
  { sprintf(mfname,"%s/%s%d.mdl",modelDir,mf[i],modelNum);
    stat(mfname,&buff);
    if(buff.st_mtime > tt) break;
  }
  free(mfname);
  if(i<4) {unlink(fname); return 1;}
  return 0;
}

static void* newSymbol(void*handle,char *name)
{
#ifdef __hpux
void * addr;
     if(shl_findsym((shl_t*)&handle,name,TYPE_UNDEFINED,&addr)) return NULL;
       else return addr;
#else
      return dlsym(handle, name);
#endif
}

static void * dLoad(char * libName)
{
void *q;

if(access(libName,R_OK)) return NULL;

#ifdef __hpux
   return  shl_load(libName,0,0L);
#else
   q= dlopen(libName, RTLD_NOW);
   if(!q) printf("%s\n",dlerror()); 
   return q;
#endif
}


static void dClose(void * handle)
{
#ifdef __hpux
       shl_unload(handle);
#else
       dlclose(handle);
#endif
}

extern double aWidth(char *);

//=======  SQME library generation ==================

static numout* loadLib(void* handle, char * lib)
{ numout * cc=malloc(sizeof(numout));
  char name[100];
  if(!handle) {free(cc); return NULL;}   
  cc->handle=handle;
  sprintf(name,"interface_%s",lib);
  cc->interface=newSymbol(handle, name);
  if(!cc->interface || cc->interface->nprc==0){free(cc); return NULL;}
  else
  {  int i;
     cc->init=0;
     cc->Q=NULL, cc->SC=NULL;
     cc->link=malloc(sizeof(double*)*(1+cc->interface->nvar));
     cc->link[0]=NULL;
     for(i=1;i<=cc->interface->nvar;i++) 
     { char *name=cc->interface->varName[i];
       cc->link[i]=varAddress(name);
if(cc->link==NULL) printf("No link for %s\n",name);       
       if(strcmp(name,"Q")==0) cc->Q=cc->interface->va+i; 
       else if(strcmp(name,"SC")==0) cc->SC=cc->interface->va+i;
     }  
     *(cc->interface->aWidth)=&aWidth;
  }
  return cc;
}


typedef struct  procRec 
{ struct procRec  * next;
  char * libname;
  numout * cc;
}  procRec;   

static  procRec* allProc=NULL;

#include<stdlib.h>
#include<sys/wait.h>

numout*getMEcode(int twidth,int Gauge, char*Process, char*excludeVirtual, 
                          char*excludeOut,char*lib)
{
   char *proclibf,*command;
   void * handle=NULL;
   int new=0;
   numout * cc;
   procRec*test;
   int Len;
   char * lib_;

   lib_=malloc(strlen(lib)+6);
   
   if(Gauge) sprintf(lib_,"%s_u",lib);    else  strcpy(lib_,lib); 
   if(twidth) strcat(lib_,"_t");  
   for(test=allProc;test; test=test->next) if(strcmp(lib_,test->libname)==0) 
   {
     *(test->cc->interface->gtwidth)=0;
     *(test->cc->interface->twidth) =twidth;
     *(test->cc->interface->gswidth)=0;
     *(test->cc->interface->BWrange)=BWrange; 
     free(lib_);
     return test->cc;
   }
   
   Len=strlen(compDir)+strlen(lib)+strlen(libDir)+300;
   proclibf=malloc(Len);

   if(Process) Len+=strlen(Process);
   if(excludeVirtual) Len+=strlen(excludeVirtual);
   if(excludeOut) Len+=strlen(excludeOut);
   command=malloc(Len);
 

   sprintf(proclibf,"%s/%s.so",libDir,lib_);

   if(access(proclibf,R_OK)==0 && checkMtime(proclibf)==0) handle=dLoad(proclibf);
   if(!handle)
   {  int i;
   
      for(i=0;Process[i]==' ';i++); if(Process[i]==0)
      {    
        free(command); free(proclibf); free(lib_); 
        return NULL;    
      }      
      if(!handle)
      {
        char options[20];
        char nPROCSSch[20];
        char GaugeCh[4];
        int ret;  
        int delWorkDir;
        if(twidth) strcpy(options,"5[{}");else strcpy(options,"");
        if(Gauge) strcpy(GaugeCh,"U"); else strcpy(GaugeCh,"F");
 
        delWorkDir=prepareWorkPlace();

        sprintf(command,"cd %s; %s/sbin/newProcess %s %s \"%s\" %s \"%s\"",
                       compDir, calchepDir, lib_, libDir,options,GaugeCh,Process);
  
        if(excludeVirtual) sprintf(command+strlen(command)," \"%s\"",excludeVirtual);
        else  sprintf(command+strlen(command)," \"\"");            
        if(excludeOut) sprintf(command+strlen(command)," \"%s\"",excludeOut);       

        sprintf(nPROCSSch,"%d",nPROCSS);
        setenv("nParProc",nPROCSSch,1);
        
        ret=system(command);
      
        if(ret<0 || WIFSIGNALED(ret)>0 ) exit(10);
        if(delWorkDir )cleanWorkPlace();
        
        if(ret==0) handle=dLoad(proclibf); else 
        { printf(" Can not compile %s \n", Process);
          free(command); free(proclibf); free(lib_);
          return NULL;
        } 
        if(!handle)
        { printf(" Can not load the compiled library %s \n",proclibf);
           free(command); free(proclibf); free(lib_);
          return NULL;
        }         
        new=1;   
      }
   }
   cc=loadLib(handle,lib_);
   if(!cc && new) dClose(handle);
   if(cc)
   {  test=(procRec*)malloc(sizeof(procRec));
      test->next=allProc; allProc=test;
      test->libname=(char*) malloc(strlen(lib_)+1);
      strcpy(test->libname,lib_);
      test->cc=cc;
   } else if(new) dClose(handle);  
    free(command); free(proclibf); free(lib_);
    if(cc) *(cc->interface->twidth) =twidth;
    return cc; 
}

// ================== Lagrangian vertexes library =======================

typedef struct  vertRec 
{ struct vertRec * next;
  char * libname;
  lVert * vv;
}  vertRec;   

static  vertRec* allVert=NULL;

static  lVert* loadVert(void* handle, char * lib)
{ if(!handle) return NULL;
  lVert * vv=newSymbol(handle, lib);
  
  if(!vv){free(vv); return NULL;}
  else
  {  int i;
     vv->handle=handle;   
     vv->init=0;
     vv->link=malloc(sizeof(REAL*)*(vv->nVar));
     for(i=0;i<vv->nVar;i++) 
     { char *name=vv->varNames[i];
       vv->link[i]=varAddress(name);
     }  
  }
  return vv;
}


static int vert2name(char**field ,char * lib)
{ 
  int i,n=0;
  strcpy(lib,"vert_");

  for(i=0;i<4;i++) 
  { 
    char bufflib[20];
    if(field[i]) 
    { 
      int err=pname2lib(field[i],bufflib);
      if(err) return i+1;
      strcat(lib,bufflib);
      n++;     
    }
  }
  if(n<3) return 5; else return 0;
} 


lVert* getVertex(char*n1,char*n2,char*n3,char*n4)
{
   char lib[50];
   char *proclibf,*command;
   void * handle=NULL;
   int new=0;
   lVert * vv;
   vertRec*test;
   int Len;
   char *fields[4]={n1,n2,n3,n4};

   int err=vert2name(fields,lib);
   if(err) return NULL;
   
   for(test=allVert;test; test=test->next) if(strcmp(lib,test->libname)==0)
   { if(test->vv->nTerms==0) return NULL;  else  return   test->vv;
   }  
   
   Len=strlen(compDir)+strlen(lib)+strlen(libDir)+300;

   proclibf=malloc(Len);
   command=malloc(Len);
 
   sprintf(proclibf,"%s/%s.so",libDir,lib);
   if(access(proclibf,R_OK)==0 && checkMtime(proclibf)==0) handle=dLoad(proclibf);
   if(!handle)
   {  int i;
      char vert_txt[50];
      strcpy(vert_txt,"");
      for(i=0; i<4;i++) if(fields[i]) sprintf(vert_txt+strlen(vert_txt)," %s",fields[i]);  

      if(!handle)
      {
        int ret;  
        int delWorkDir;
        
        delWorkDir=prepareWorkPlace();
        sprintf(command,"cd %s; %s/sbin/newVertex models  1  %s  %s %s",
                       compDir, calchepDir,  lib, libDir, vert_txt);
        ret=system(command);
      
//        if(ret<0 || WIFSIGNALED(ret)>0 ) exit(10);
        if(delWorkDir )cleanWorkPlace();
        
        if(ret==0) handle=dLoad(proclibf); else 
        { printf(" Can not compile codes for vertex  [%s] \n",  vert_txt);
          free(command); free(proclibf);
          return NULL;
        } 
        if(!handle)
        { printf(" Can not load the compiled library %s \n",proclibf);
           free(command); free(proclibf); 
          return NULL;
        }         
        new=1;   
      }
   }
   vv=loadVert(handle,lib);
   if(!vv && new) dClose(handle);
   if(vv)
   {  test=(vertRec*)malloc(sizeof(vertRec));
      test->next=allVert; allVert=test;
      test->libname=(char*) malloc(strlen(lib)+1);
      strcpy(test->libname,lib);
      test->vv=vv;
   } else if(new) dClose(handle);  
   free(command); free(proclibf);
   if(vv->nTerms==0) return NULL;
   else  return vv; 
}


int  getNumCoeff(lVert*vv, double * coeff)
{
   int i;   
   if(!vv) { printf("Call of  getNumCoeff with NULL argument. Program terminated\n"); exit(1);}
 
   for(i=0;i<vv->nVar;i++) if(vv->link[i] &&  ((REAL*)(vv->link[i])-varValues) < *currentVarPtr) vv->varValues[i+1]=*((REAL*)(vv->link[i]));
   else 
   { 
     printf("!Value of variable  '%s' needed for calculation of '%s' is not known yet\n",
     vv->varNames[i], varNames[*currentVarPtr]);  
     vv->varValues[i]=0; FError=1;
   }
                                     
   int err= vv->calcCoeff(coeff);
   return err;
}



//=======  generation of process list ======     

txtList  makeProcList(char ** InNames, char** OutNames, int nx)
{ 
  char lname[20],buff[200];
  char * fnameG;
  FILE *f;
  txtList List=NULL;
  char *process;
  char ** c;
  

  char* allNames[10];
  char* alias[10][10];
  int nin=0,nout=0,nnew=0;
  for(int i=0;InNames[i];i++) {allNames[i]=InNames[i];  nin++;}
  if(OutNames) for(int i=0;OutNames[i];i++){allNames[nin+i]=OutNames[i]; nout++;}
  for(int i=0;i<nin+nout;i++) 
  {  int j;
     for(j=0;  j<i;j++){ if(allNames[i]==allNames[j] || strcmp(allNames[i],allNames[j])==0) break;} 
     if(i==j){ nnew++; sprintf(alias[i],"#_%d",nnew);} else alias[i][0]=0;
  } 

  int len=0; for(int i=0;i<nin+nout;i++) if(alias[i][0]) len+=strlen(allNames[i]);
  
  process=malloc(len+50);
  strcpy(process,alias[0]);  
  for(int i=1;i<nin+nout;i++)
  {  
     if(i==nin) strcat(process,"->"); else strcat(process,",");
     for(int j=0;j<=i;j++) if(allNames[i]==allNames[j] || strcmp(allNames[i],allNames[j])==0)
     { sprintf(process+strlen(process),"%s", alias[j]); break;}
  }

  if(nx){ if(nout) strcat(process,",");  sprintf(process+strlen(process),"%d*x{",nx); } 
  for(int i=0;i<nin+nout;i++) if(alias[i][0]) sprintf(process+strlen(process),"{%s",allNames[i]);

//printf("command:\n%s\n",command);      

  char * command=malloc(strlen(compDir) + strlen(calchepDir) + strlen(process) + 100);
  int delWorkDir=prepareWorkPlace();
  sprintf(command,"cd %s; %s/bin/s_calchep -blind \"{{%s{{[[{0\" >/dev/null",
                  compDir,calchepDir,                process);   
//printf("command=%s\n", command);                  
  system(command);
  free(command); free(process);     
  fnameG=malloc(strlen(compDir)+50);
  sprintf(fnameG,"%s/results/list_prc.txt",compDir);        
  f=fopen(fnameG,"r");
  free(fnameG);

  if(!f) {if(delWorkDir) cleanWorkPlace();  return NULL;}
  for(; 1==fscanf(f,"%[^\n]\n",buff); )
  { txtList l=malloc(sizeof(txtListStr));
    l->next=List; List=l;
    l->txt=malloc(strlen(buff)+1);
    strcpy(l->txt,buff);
  }
     
  fclose(f);
       
  if(delWorkDir) cleanWorkPlace();
  
  return List;
}


txtList  makeDecayList(char * pname, int nx)
{ 
  char fnameL[50],lname[20],buff[200];
  char * fnameG;
  FILE *f;
  txtList List=NULL;

//printf("asks 1->%d for %s\n",nx,pname);
  pname2lib(pname,lname);
  sprintf(fnameL,"dList_%s_%dx",lname,nx);
  fnameG=malloc(strlen(libDir)+50);
  sprintf(fnameG,"%s/%s",libDir,fnameL);
  if(access(fnameG, R_OK)|| checkMtime(fnameG))
  {  char * command=malloc(strlen(compDir) + strlen(calchepDir)+strlen(libDir) +200);
     int delWorkDir=prepareWorkPlace();
     sprintf(command,"cd %s;"
       " %s/bin/s_calchep -blind \"{{%s->%d*x{{{[[{0\" >/dev/null ;"
       " if(test $? -eq 0) then mv results/list_prc.txt %s/%s;fi",
       compDir,calchepDir,pname,nx,libDir,fnameL);
     system(command);
     free(command);
     if(delWorkDir) cleanWorkPlace();
  }

  f=fopen(fnameG,"r");
  free(fnameG);
  if(!f) return NULL;
  for(; 1==fscanf(f,"%[^\n]\n",buff); )
  { txtList l=malloc(sizeof(txtListStr));
    l->next=List; List=l;
    l->txt=malloc(strlen(buff)+1);
    strcpy(l->txt,buff);
  }
  fclose(f);
  return List;
}


void cleanTxtList(txtList L)
{ 
   while(L) {txtList l=L; free(L->txt); L=L->next; free(l);}  
}

void printTxtList(txtList L, FILE *f)
{ for(;L;L=L->next) fprintf(f,"%s\n",L->txt);}



void delAllLib(void)
{
  procRec* curProc=allProc;
  while(curProc)
  {  procRec*tmp=curProc;
     free(curProc->libname);
     free(curProc->cc->link);
     dClose(curProc->cc->handle);     
     free(curProc->cc);
     curProc=curProc->next;
     free(tmp);
  }
 
  allProc=NULL;

  vertRec*curVert=allVert;
  
  while(curVert)
  {  vertRec*tmp=curVert;
     free(curVert->libname);
     free(curVert->vv->link);
     dClose(curVert->vv->handle);
     curVert=curVert->next;
     free(tmp);
  }
                                       
  allVert=NULL;
}

numout*newProcess(char*Process)
{  int err;
   char lib[100];
   err=process2Lib(Process,lib);
   if(err) return NULL;
   return getMEcode(1,ForceUG,Process,NULL,"",lib);
}


//============================

int HiggsLambdas(double Q, char * Higgs, double complex*lAA,double complex *lGG, double complex *lAA5, double complex *lGG5)
{               
   if(lGG)*lGG=0;  if(lGG5)*lGG5=0; if(lAA)*lAA=0; if(lAA5)*lAA5=0;
   int i;
   for(i=0;i<nModelParticles;i++) if(strcmp(ModelPrtcls[i].name,Higgs)==0) break;
   if(i==nModelParticles) { printf(" %s - no such particles\n",Higgs);   return 0;}
   if( ModelPrtcls[i].spin2!=0 ||  ModelPrtcls[i].cdim!=1 || ModelPrtcls[i].q3!=0) 
   { printf(" %s - has not quantum numbers of Higgs\n",Higgs);   return 0;}

   int pdg= ModelPrtcls[i].NPDG;

   double complex ffE=0,faE=0,ffC=0,faC=0; 
  
   txtList L = makeDecayList(Higgs,2), l=L;

   double a=alphaQCD(Q)/M_PI;
   for(l=L;l;l=l->next)
   { char Xp[10], Xm[10];
     sscanf(strstr(l->txt,"->")+2, "%[^,],%s",Xp,Xm);
     trim(Xp); trim(Xm); 

     if(strcmp(Xp,antiParticle(Xm))==0)
     {  
        int pdg,spin2,charge3,cdim;
        double mX=pMass(Xp);
        double dffE=0,dfaE=0,dffC=0,dfaC=0;
        pdg=qNumbers(Xp, &spin2, &charge3, &cdim);

        if((charge3 !=0 || cdim!=1) && mX>0)
        {  
           double coeff[10]; 
           int k;
           double mXp; // pole mass
           switch(abs(pdg))
           { case 4: mXp=1.48;         break;
             case 5: mXp=bPoleMass();  break;
             case 6: mXp=tPoleMass();  break;
             default:mXp=mX;
           }             
           double mN= (spin2&1)?  mX : mX*mX;   
           lVert *xxh=getVertex(Xm,Xp,Higgs,NULL);
           if(!xxh) continue;
           getNumCoeff(xxh,coeff);
           for(k=0;k<xxh->nTerms;k++)  
           { int addFF=0,addFA=0;
             switch(spin2)
             {
               case 0:  if(strcmp(xxh->SymbVert[k],"1")==0) addFF=1;    break;
               case 1:  if(strcmp(xxh->SymbVert[k],"1")==0) addFF=1; else
                        if(strcmp(xxh->SymbVert[k],"G5*i")==0) addFA=1; break;
               case 2:  if(strcmp(xxh->SymbVert[k],"m2.m1")==0) addFF=1;break; 
             }
               if(addFF)
               { if(lGG && cdim!=1)  *lGG+=hGGeven(Q,a,1,spin2,cdim,(REAL)mXp,(REAL)(coeff[k]/mN)); 
                 if(lAA && charge3 ) *lAA+=hAAeven(Q,a,1,spin2,cdim,(REAL)mXp,(REAL)(coeff[k]/mN))*charge3*charge3/9.;
               }
               if(addFA) 
               { if(lGG5 && cdim!=1) *lGG5+=0.5*hGGodd(Q,a,1,spin2,cdim,(REAL)mXp,(REAL)(coeff[k]/mN)); 
                 if(lAA5 && charge3) *lAA5+=0.5*hAAodd(Q,a,1,spin2,cdim,(REAL)mXp,(REAL)(coeff[k]/mN))*charge3*charge3/9.; 
               } 
             }   
          }                     
        }    
      }
      cleanTxtList(L);
      return pdg;
}

static int nHiggs=0;
typedef  struct { char hName[20]; double complex lXX[5];} lambdaRec;
static lambdaRec* allLambdas=NULL;

void cleanHiggs_AA_GG(void) { free(allLambdas); allLambdas=NULL; nHiggs=0;} 
 

double complex lAAhiggs(double Q, char * hName) { double complex lAA;  if(HiggsLambdas(Q,hName, &lAA,NULL,NULL, NULL))  return lAA; else  { FError=1;  return 0;} }
double complex lGGhiggs(double Q, char * hName) { double complex lGG;  if(HiggsLambdas(Q,hName, NULL,&lGG,NULL, NULL))  return lGG; else  { FError=1;  return 0;} } 
double complex lAA5higgs(double Q,char * hName) { double complex lAA5; if(HiggsLambdas(Q,hName, NULL,NULL,&lAA5,NULL))  return lAA5; else { FError=1;  return 0;} }
double complex lGG5higgs(double Q,char * hName) { double complex lGG5; if(HiggsLambdas(Q,hName, NULL,NULL,NULL,&lGG5))  return lGG5; else { FError=1;  return 0;} }
