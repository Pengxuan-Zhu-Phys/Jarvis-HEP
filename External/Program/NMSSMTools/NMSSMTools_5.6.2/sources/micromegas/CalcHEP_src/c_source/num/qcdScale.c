#include <dlfcn.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "phys_val.h"
#include "VandP.h"
#include "parser.h"
#include "num_out.h"
#include "rootDir.h"
#include "subproc.h"
#include "interface.h"
#include "qcdScale.h"


static int NX, pos[4];
static char* sname[4]={"Qren","Qpdf1","Qpdf2","Qshow"};

static FILE*fout;
static double (*scaleCC_)(int nsub,REAL*,double (*calcPV)(char,char*,double*), double *,double*,double*,double*,double*) =NULL;
static void* scaleLib=NULL;


static void*  rd_num_(char* s)
{
   char*p;
   char key, plist[20];
   int i;
   
   if(isdigit(*s))
   { p=malloc(strlen(s)+10); 
     sprintf(p,"%s",s);
     return p;
   }
   else if(s[0]=='"') 
   {  p=malloc(strlen(s)+10);
      sprintf(p,"%s",s);
      return p;
   }      
   else 
   { p=malloc(20);
     sprintf(p,"X[%d]",NX++);

     if(checkPhysVal(s,&key,plist)) 
     { fprintf(fout,"  %s=calcPV('%c',\"",p,key);
       for(i=0;plist[i];i++) fprintf(fout,"\\%d",plist[i]);
       fprintf(fout,"\",pvect);\n");   
        return p;
     }

     for(i=0;i<nModelVars+nModelFunc;i++) if(strcmp(varNames[i],s)==0)
     { fprintf(fout,"%s = modelVal[%d];  /* %s  */ \n",p,i, varNames[i]); 
       return p;
     }
         
     for(i=0;i<3;i++) if(pos[i]>=0 && strcmp(s,sname[i])==0) 
     {  p=malloc(20);
        sprintf(p,"X[%d]",pos[i]);  
        return p; 
     }
   }  
   return NULL;
}


static void*  act_num_(char* ch,int n, void**args)
{  char*  p=NULL,*p1=NULL, *p2=NULL,*p3=NULL;
   p1= (char*)args[0];
   if(n>=2) p2=(char*)args[1];
   if(n>=3) p3=(char*)args[2];
   
   switch (ch[0])
   { 
     case '-': p=malloc(strlen(p1)+10);
               sprintf(p,"(-(%s))",p1);
               return p;
     case '+': 
     case '*': p=malloc(strlen(p1)+strlen(p2)+10);
               sprintf(p,"(%s)%c(%s)",p1,ch[0],p2);
               return p;
     case '/':p=malloc(strlen(p1)+strlen(p2)+15);
              sprintf(p,"(%s)/(double)(%s)",p1,p2);  
              return p;           
     
     case '^': p=malloc(strlen(p1)+strlen(p2)+10);
               sprintf(p,"pow(%s,%s)",p1,p2);
               return p;
     case '.': rderrcode=typemismatch; 
               return NULL; 
     default : switch(n) 
               {      
                   case 1: 
                      if(!strcmp(ch,"sqrt")||!strcmp(ch,"sin")||!strcmp(ch,"cos")||!strcmp(ch,"tan")
                       ||!strcmp(ch,"asin")||!strcmp(ch,"acos")||!strcmp(ch,"atan")||!strcmp(ch,"exp")
                       ||!strcmp(ch,"log")||!strcmp(ch,"fabs"))
                       { 
                         p=malloc(strlen(p1)+10);
                         sprintf(p,"%s(%s)",ch,p1);
                         return p;
                       }   
                       break;
                   case 2:
                      if(!strcmp(ch,"atan2"))
                      { p=malloc(strlen(p1)+strlen(p2)+10);
                        sprintf(p,"atan2(%s,%s)",p1,p2);
                        return p;   
                      }
                      break; 
                   case 3:   
                      if(!strcmp(ch,"if"))
                      { p=malloc(strlen(p1)+strlen(p2)+strlen(p3)+10);
                        sprintf(p,"(%s>0 ? %s:%s)",p1,p2,p3);
                        return p;
                      }  
               }

               if(strchr("ACDEJKMPSTUYNWZ",ch[0]) && (ch[1]==0 || ( ch[2]==0 && (ch[1]=='_' || ch[1]=='`') )
               || !strcmp(ch,"S1") || !strcmp(ch,"S1_") || !strcmp(ch,"S1`")
               || !strcmp(ch,"S2") || !strcmp(ch,"S2_") || !strcmp(ch,"S2`")      )) 
               { char buff[100];
                 int i;
                 
                 sprintf(buff,"%s(",ch);
                 for(i=0;i<n;i++) sprintf(buff+strlen(buff),"%s,",(char*)args[i]);
                 buff[strlen(buff)-1]=')';
                 int Nsub_mem=Nsub;
                 fprintf(fout,"  switch(nsub)\n  {\n");
                 for(Nsub=1;Nsub<=nprc_int;Nsub++)
                 {  fprintf(fout,"   case %d:\n   { ",Nsub);
                    char key[3];
                    physValRec *pList,*p;
                    if(!checkPhysValN(buff, key, &pList)) return NULL;
                    fprintf(fout,"char * plist[]={");
                    for(p=pList;p;p=p->next) {fprintf(fout,"\"");
                      for(i=0;i<strlen(p->pstr);i++) fprintf(fout,"\\%d",p->pstr[i]); fprintf(fout,"\",");
                     } fprintf(fout,"NULL};\n");
                    cleanPVlist(pList); 
                    fprintf(fout,"     for(i=1,ss=calcPV('%c',plist[0],pvect); plist[i];i++) ",key[0]);
                    switch (key[1]) 
                    { case 0  : fprintf(fout," ss+=calcPV('%c',plist[i],pvect);\n",key[0]);  break;
                      case '_': fprintf(fout,"{ s1=calcPV('%c',plist[i],pvect); if(s1<ss) ss=s1;} ",key[0]);  break;
                      case '`': fprintf(fout,"{ s1=calcPV('%c',plist[i],pvect); if(s1>ss) ss=s1;} ",key[0]);  break;
                    }                   
                    fprintf(fout,"   } break;\n",NX);   
                    
                 }
                 fprintf(fout,"  };\n");
                 fprintf(fout,"  X[%d]=ss;\n",NX); 
                 Nsub=Nsub_mem;             
                 
                 p=malloc(30);
                 sprintf(p,"X[%d]",NX++);
                 return p;
               }    
               
    }                    

   if(strcmp(ch,"min")==0 || strcmp(ch,"max")==0)
   { int l=0,i;
     char *q;
     for(i=0;i<n;i++) l+=strlen((char*)args[i])+8;
     p=malloc(l);
     q=malloc(l);
     strcpy(p,(char*)args[0]);
     for(i=1;i<n;i++)
     { strcpy(q,p);
       sprintf(p,"%s(%s,%s)",ch,q,(char*)args[i]); 
     }
     free(q);
     return p;
   }
   
   return NULL;
}

int initScales(char*sR, char*sF1,char*sF2, char*sS, char * mess)
{  char *ch;
   char *command;
   int k;
   char*ss[4];  
   ss[0]=sR; ss[1]=sF1; ss[2]=sF2, ss[3]=sS;
   for(k=0;k<4;k++) pos[k]=-1;
   long fpos;
       
   if(scaleLib) dlclose(scaleLib); 
   scaleCC_=NULL;
   
   fout=fopen("scale.c","w");
   if(!fout){ if(mess) sprintf(mess,"can't open file scale.c for writing");  return -1;}
   NX=0;
   fprintf(fout,"#include<math.h>\n");
   fprintf(fout,"#include<stdlib.h>\n");  
   fprintf(fout,"#include\"%s/include/nType.h\"\n",rootDir);
   fprintf(fout,"#define min(x,y) (x<y? x:y)\n");
   fprintf(fout,"#define max(x,y) (x>y? x:y)\n"); 
   fprintf(fout,"extern void ScaleCC(int nsub, REAL*, double (*calcPV)(char,char*,double*), double*,double*,double*,double*,double *);\n");
   fprintf(fout,"void ScaleCC(int nsub,REAL*modelVal, double (*calcPV)(char,char*,double*), double *pvect,double *%s, double *%s,double *%s, double*%s)\n",
            sname[0],sname[1],sname[2],sname[3]);
   fprintf(fout,"{ double ");
   fpos=ftell(fout);
   fprintf(fout,"               \n");
   fprintf(fout,"  double ss,s1; int i;\n");

   for(k=0;k<4;k++)
   {
      ch=readExpression(ss[k],rd_num_,act_num_,free);
      if(!ch) 
      { fclose(fout); 
        fout=NULL;
        if(mess)
        { 
           sprintf(mess,"Error in %s scale definition: position %d :\n%s", 
                   sname[k], rderrpos, errmesstxt(rderrcode));
        }
        return rderrpos;    
      } 
      pos[k]=NX;
      fprintf(fout,"X[%d]=%s;\n",NX++,ch);
      fprintf(fout," *%s= X[%d];\n",sname[k],pos[k]);
      free(ch);
   }   
   
   fprintf(fout,"}\n");
   fseek(fout,fpos,SEEK_SET);
   fprintf(fout,"X[%d];",NX+1);
   fclose(fout);
   
    command=malloc(strlen(rootDir)+100);
    sprintf(command,". %s/FlagsForSh; $CC $CFLAGS $SHARED -o scale.so scale.c",rootDir);
    system(command);
    free(command);
    scaleLib=dlopen("./scale.so",RTLD_NOW);
    if(!scaleLib)
    {
       if(mess) sprintf(mess,"Can't load shared library for QCD scale"); 
       return -3;
    } 
    scaleCC_=dlsym(scaleLib,"ScaleCC");
    if(!scaleCC_)
    {  
       if(mess) sprintf(mess,"Problem in library scale.so ");
       return -4;
    }
    return 0;
}


void Scale(int nsub, double*pv, double *qR, double *qF1,double *qF2,double *qS)
{ double q=91.187; 

  if(scaleCC_)  
  { 
     (*scaleCC_)(nsub,varValues, calcPhysVal, pv,qR,qF1,qF2,qS); 
     if(*qF1<0) *qF1*=-1; 
     if(*qF2<0) *qF2*=-1;
     if(*qR<0) *qR*=-1; if(*qR<1) *qR=1;
     if(qS){ if(*qS<0) *qS*=-1; if(*qS<1) *qS=1;}
  } 
  else {*qR=q; *qF1=q; *qF2=q; if(qS) *qS=q;}
}
