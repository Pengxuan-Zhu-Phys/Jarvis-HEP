/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include <unistd.h>
#include <dirent.h>
#include <pthread.h>

#include "n_proc.h"
#include "crt_util.h"
#include "syst.h"
#include "syst2.h" 
#include "s_files.h"

 FILE * menup;
 FILE * menuq;
 FILE * diagrp;   /* file of Adiagram; */
 FILE * diagrq;   /* file of CSdiagram; */

 FILE * catalog;
 FILE * archiv;


char*mdFls[5] = {"vars","func","prtcls","lgrng","extlib"};
shortstr  pathtouser;

void wrt_menu(int menutype,int k,char*txt,int ndel,int ncalc,int nrest,long recpos)
{
   if (menutype == 1)
   {
      fseek(menup,(k - 1)*71 + 2,SEEK_SET);
      fprintf(menup,"%4d| %-44.44s|%5d|%5d|%-8ld",k,txt,ndel,nrest,recpos);
   }
   else
   {
      fseek(menuq,(k - 1)*77 + 2,SEEK_SET);
      fprintf(menuq,"%4d| %-44.44s|%5d|%5d|%5d|%-8ld",k,txt,ndel,ncalc,nrest,recpos);      
   }
}


int rd_menu(int menutype,int k,char*txt,int*ndel,int*ncalc,int*nrest,long*recpos)
{
   if (menutype == 1)
   {
      fseek(menup,(k - 1)*71 + 2,SEEK_SET);
      if(5!=fscanf(menup,"%d| %[^|]%*c%d|%d|%ld",&k,txt,ndel,nrest,recpos)) return 0;
      *ncalc=0;
   }
   else
   {
      fseek(menuq,(k - 1)*77 + 2,SEEK_SET);
      if(6!=fscanf(menuq,"%4d| %[^|]%*c%d|%d|%d|%ld",&k,txt,ndel,ncalc,nrest,recpos))
           return 0;
   }
   return 1;
}

int whichArchive(int nFile,int rw)
{ static int ArchNum=0, rw_=0;
  char archivName[40];  
  if(nFile==0)
  { if(ArchNum!=0) {ArchNum=0; fclose(archiv);}
    return 0;
  }
  
  if(ArchNum==nFile  && rw==rw_ )
  {  if(rw=='w' && ftell(archiv) >= MAXARCHIVE) 
     { fclose(archiv);
       ArchNum++;
       sprintf(archivName,ARCHIV_NAME,ArchNum);
       archiv=fopen(archivName,"w");
     } 
     return ArchNum;
  } 
  
  if(ArchNum) fclose(archiv); 
  ArchNum=nFile;  rw_=rw;
  sprintf(archivName,ARCHIV_NAME,ArchNum);
  if(rw=='w')
  { if(access(archivName,F_OK)==0) archiv=fopen(archivName,"a");
    else                           archiv=fopen(archivName,"w");
  }
  else  archiv=fopen(archivName,"r");
   
  return ArchNum;
} 

  typedef char nameRec[12];
  static nameRec * list=NULL;
  static int nFiles=0,cFile;


static pthread_mutex_t keyN=PTHREAD_MUTEX_INITIALIZER;

static int N,Ntot,nFor1,fb,fe,db,de,sn,n_prc;

static void * pCompile_cycle(int * input)
{  
  int i,ibeg,iend;
  char src[20000];
  char command[20000];
  char arname[20];
  int myPrc=input[0];
  int dN=0;      
  for(;;)
  { char *c;
    src[0]=0;
    pthread_mutex_lock(&keyN);
    if(sn)
    {  
       sn=0;
       pthread_mutex_unlock(&keyN);       
       dN=3;
       
       if(access("service.c",R_OK)==0) strcat(src," service.c");
       if(access("sqme.c",R_OK)==0)    strcat(src," sqme.c"); 
       if(access("VandP.c",R_OK)==0)   strcat(src," VandP.c");
       if(strlen(src))
       {        
//printf("sn: src=%s\n",src);
          sprintf(command,". $CALCHEP/FlagsForSh\n"
                          "SRC=\"%s\"\n" 
                          "$CC $CFLAGS -c -I$CALCHEP/include $SRC\n"
                          "rm $SRC\n",src);
          system(command);
          for(c=strstr(src,".c"); c; c=strstr(c,".c")) c[1]='o';
          sprintf(command, ". $CALCHEP/FlagsForSh\n"
                           "SRC=\"%s\"\n"
                           "ar r lib_0.a $SRC\n"
                           "$RANLIB  lib_0.a \n" 
                           "rm $SRC",src);
          system(command);
       }                     
    } else if(de>=db)
    {
      ibeg=db;
      db+=nFor1;
      if(db>de+1)db=de+1;
      iend=db;
      pthread_mutex_unlock(&keyN);
      dN=iend-ibeg;
      for(i=ibeg;i<iend;i++) sprintf(src+strlen(src)," d%d.c",i);
//printf("d: src=%s\n",src);       
      sprintf(command,". $CALCHEP/FlagsForSh\n" 
                      "SRC=\"%s\"\n $CC $CFLAGS -c -I$CALCHEP/include $SRC\n rm $SRC",src); system(command); 
      for(c=strstr(src,".c"); c; c=strstr(c,".c")) c[1]='o';
      for(;;){  sprintf(arname,"ld%d.a",myPrc); if(access(arname,R_OK)) break; else myPrc+=nPROCSS;}      
      sprintf(command,". $CALCHEP/FlagsForSh\n" 
                      "SRC=\"%s\"\n ar r %s $SRC\n rm $SRC\n ",src,arname);  system(command);  
    }  else if(fe>=fb)
    {
      ibeg=fb;
      fb+=nFor1;
      if(fb>fe+1)fb=fe+1;
      iend=fb;
      pthread_mutex_unlock(&keyN);
      dN=iend-ibeg;
      for(i=ibeg;i<iend;i++) sprintf(src+strlen(src)," f%d.c",i);
      for(;;){  sprintf(arname,"lf%d.$SO",myPrc); if(access(arname,R_OK)) break; else myPrc+=nPROCSS;}   
      sprintf(command,". $CALCHEP/FlagsForSh\n" 
                      "SRC=\"%s\"\n"
                      "if(test -z  $SONAME) then extName= ; else extName=\"$SONAME $PWD/%s\"; fi\n"
                      "$CC  $CFLAGS  -I$CALCHEP/include $SHARED  -o $PWD/%s $extName  $SRC\n rm $SRC",src,arname,arname);   system(command); 
    } 
    else { n_prc--;  pthread_mutex_unlock(&keyN); return NULL;}
    myPrc+=nPROCSS;
    pthread_mutex_lock(&keyN);
    N+=dN;
    pthread_mutex_unlock(&keyN);   
  } 
  return NULL;
}

static int comp_control(void * par_)
{  
  int err=0;
  infoLine(0.);
  for(;;)
  { double r;
    sleep(1);
    pthread_mutex_lock(&keyN);
    if(!n_prc) { pthread_mutex_unlock(&keyN); infoLine(2.);  return err; }
    r=(double)N/(double)Ntot;
    
    pthread_mutex_unlock(&keyN);
    if(infoLine(r)) 
    { 
      pthread_mutex_lock(&keyN);
       sn=0;
       db=1;de=0;
       fb=1;fe=0;  
      pthread_mutex_unlock(&keyN);
      int Y=where_y();
      goto_xy(15,Y); scrcolor(Red,White); print("Esc signal is accepted. ");
      err=1;
    }
  }   
}


int pCompile(void)
{
  DIR *dirPtr=opendir("./");
  struct dirent * dp;
  int i,k,n,err;
  char command[5000];
  char *c;
    
  if(!dirPtr) return 1;
  
  fb=0;fe=0;db=0,de=0; sn=3;
  while((dp=readdir(dirPtr)))
  { char *c=strstr(dp->d_name,".c");
    if( dp->d_name[0]!='.' && c && c[2]==0) 
    { if(sscanf(dp->d_name,"f%d",&n)==1){ if(!fb) {fb=n;fe=n;} else {if(n>fe) fe=n; else if(n<fb) fb=n;} }
      if(sscanf(dp->d_name,"d%d",&n)==1){ if(!db) {db=n;de=n;} else {if(n>de) de=n; else if(n<db) db=n;} }
    }  
  }
  if(fe==0) fb=1;
  if(de==0) db=1;
  closedir(dirPtr);
  Ntot=5+de-db+fe-fb;
  n_prc=nPROCSS;
  N=0;
//printf(" sn=%d db=%d de=%d fb=%d fe=%d\n", sn,db,de,fb,fe);   
  if(Ntot)
  { 
     for(nFor1=Ntot/nPROCSS+1;  nFor1>1500; nFor1/=2)continue;  
     pthread_t *threads = malloc(nPROCSS*sizeof(pthread_t));
     int*num=malloc(nPROCSS*sizeof(int));
     pthread_t control;  
     for (k=0;k<nPROCSS;k++) num[k]=k+1;  
     for (k=0;k<nPROCSS;k++) pthread_create(threads+k,NULL,pCompile_cycle, num+k);
     err=comp_control(NULL);     
     for (k=0;k<nPROCSS;k++)   pthread_join(threads[k],NULL);
     free(threads);
     free(num);
  }
  
  if(err==0) system(". $CALCHEP/FlagsForSh\n" 
                    ". ./EXTLIBsh\n"
                    "$CALCHEP/sbin/ld_n");
  return 0;
}
