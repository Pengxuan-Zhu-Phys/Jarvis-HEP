/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include <unistd.h>
#include <stdlib.h>
#include <dirent.h>
#include <pwd.h>   

#include"chep_crt.h"
#include "syst2.h"
#include "physics.h"
#include "screen.h"
#include "s_files.h"
#include "file_scr.h"
#include "read_mdl.h"
#include "sortarr.h"
#include "m_utils.h"


#define nhardmdl 4
#define newmodeltxt "   IMPORT MODEL       "

#define MAX_MODEL_NUM  ((STRSIZ-2)/22 -2)
  
void fillModelMenu(void)
{ int i,j;
  FILE * txt;
  char name[80];
  strcpy(modelmenu,"\026");
  maxmodel=0;
  for (i=1;i<=MAX_MODEL_NUM; i++)
  {  char fname[100];
     sprintf(fname,"./models%c%s%d.mdl",f_slash,mdFls[0],i);     
     txt=fopen(fname,"r");		
     if(txt==NULL) break;
     fgets(name,60,txt);
     trim(name);
     for (j=strlen(name); j<21;j++) name[j]=' ';
     name[21]=0;
     strcat(modelmenu," ");
     strcat(modelmenu,name);
     fclose(txt);
     maxmodel++;
  }  
  if ( maxmodel < MAX_MODEL_NUM ) strcat(modelmenu,newmodeltxt);
}



int  deletemodel(int n)
{
  int   i;
  char *from,*to;

  if( mess_y_n(56,10,"Delete model?"))
  {  from=malloc(strlen(pathtouser)+100);
     to=malloc(strlen(pathtouser)+100); 
     for(i=0;i<5;i++)
     {
        sprintf(to,"%smodels%c%s%d.mdl",pathtouser,f_slash,mdFls[i],n);
	unlink(to);
     }
     
     sprintf(to,"%sProcesses/m%d",pathtouser,n);
     if(access(to,R_OK)==0) { sprintf(to,"rm -r %sProcesses/m%d",pathtouser,n); system(to);}
     
     for(; n < maxmodel;n++)
     {
        for(i=0;i<5;i++)
        {
          sprintf(from,"%smodels%c%s%d.mdl",pathtouser,f_slash,mdFls[i],n+1);
          sprintf(to,  "%smodels%c%s%d.mdl",pathtouser,f_slash,mdFls[i],n);
          rename(from,to);
        }
        sprintf(from, "%sProcesses/m%d",pathtouser,n+1);
        if(access(from,R_OK)==0) 
        {  sprintf(to, "%sProcesses/m%d",pathtouser,n);
           rename(from,to); 
        }    
     }
     sprintf(to,"%sProcesses/models.txt",pathtouser);
     if(access(to,R_OK)==0) unlink(to);
     free(from);free(to);
     return 1;
  } 
  return 0;
}

int  makenewmodel(void)
{ char dirName[STRSIZ];
  char fname[STRSIZ];
  char buff[STRSIZ];
  int   nmdl=0;
  int  i, xmdl;
  DIR *dirPtr;
  FILE*f,*g;
  struct dirent * dp;  
  void *pscr = NULL;
  int *all=NULL;
  char*newM=NULL; 

  for(;;)
  { 
    int len=0;
    if(!findCalcHEPfile(dirName)) {free(all); free(newM);return 0;} 
    dirPtr=opendir(dirName);
    if(!dirPtr)
    {  messanykey(10,10,"Such directory is absent");
       continue;
    }
    while((dp=readdir(dirPtr)))
    { char tail[100];
      int n;
      if(sscanf(dp->d_name,"vars%d.%s",&n,tail)==2 && strcmp(tail,"mdl")==0)
      {
         for(i=1;i<4;i++)
         { sprintf(fname,"%s/%s%d.mdl",dirName,mdFls[i],n);  
           if(access(fname,R_OK)) break;
         }
      
         if(i==4) 
         {  
           len++;
           all=realloc(all,len*sizeof(int));
           all[len-1]=n;
         }
      }  
    }
    closedir(dirPtr);
    if(len==0)
    { messanykey(12,12," This directory does not contain\n"
                       " model files");
      continue;
    }

    SORTARR(all,len)
    newM=realloc(newM,(len+1)*32+2);
    newM[0]=32;
    sprintf(newM+1,"%32.32s","    Download all models         ");
    for(i=0;i<len;i++)
    { sprintf(fname,"%s/%s%d.mdl",dirName,mdFls[0],all[i]);
      f=fopen(fname,"r");
      fscanf(f,"%[^\n]",buff);
      trim(buff);
      sprintf(newM+1+32*(i+1),"%22.22s |%3d.mdl ",buff,all[i]); 
      fclose(f);               
    }
    menu1(5,12,"Choose a model",newM,NULL,&pscr,&nmdl);
    if(nmdl==0){len=0; continue;} 
    else for(xmdl=0;xmdl<len;xmdl++) if(nmdl==1 || nmdl-2==xmdl)
    { int n;
      sscanf(newM+1+newM[0]*(xmdl+1), "%22c |%d",buff,&n);
      buff[22]=0;
      trim(buff);
      if(nmdl!=1){correctStr(5,16,"Correct name ",buff,22,1);trim(buff);}
      if(1+xmdl+maxmodel > MAX_MODEL_NUM) { messanykey(10,10,"Too many models"); nmdl=-1;  break;}
      for(i=0;i<5;i++)
      {  int s; char ch[1000];
         sprintf(fname,"%s/%s%d.mdl",dirName,mdFls[i],n);
         f=fopen(fname,"r");
         sprintf(fname,"models/%s%d.mdl",mdFls[i],maxmodel+1);
         g=fopen(fname,"w"); fprintf(g,"%s\n",buff);
         if(f==NULL && i==4)
         { fprintf(g,"Libraries\n");
           fprintf(g,"External libraries  and citation                                      <|\n");
         }else 
         {
           fscanf(f,"%*[^\n]%*c");
           while((s=fread(ch,1,1000,f))) fwrite(ch,1,s,g);
           fclose(f); 
         }  
         fclose(g);  
      }
      maxmodel++;
    }
    fillModelMenu();
    if(nmdl>0)
    {
      if(nmdl==1) { messanykey(10,16,"All models are added");  put_text(&pscr);  return 0;}
      messanykey(10,16,"The model is added");
    }   
    put_text(&pscr);     
  }
  return 0;
}

int  continuetest(void)
{shortstr  txt;
 int  k, ndel, ncalc, nrest; 
 long recpos;

   menuq=fopen(MENUQ_NAME,"rb");
  
   for(k=1;k<=subproc_sq;k++)
   {
      rd_menu(2,k,txt,&ndel,&ncalc,&nrest,&recpos);
      if (nrest != 0)
      {
         fclose(menuq);
         return 0;
      }
   }
   fclose(menuq);
   return 1;
}


void  clear_tmp(void)
{  int i;
   char name[40];
   for(i=1;;i++) { sprintf(name,ARCHIV_NAME,i);    if(unlink(name)) break;}
   for(i=1;;i++) { sprintf(name,ARCHIV_NAME "2",i);if(unlink(name)) break;}  
   unlink(CATALOG_NAME);
   unlink(CATALOG_NAME "2");
}
