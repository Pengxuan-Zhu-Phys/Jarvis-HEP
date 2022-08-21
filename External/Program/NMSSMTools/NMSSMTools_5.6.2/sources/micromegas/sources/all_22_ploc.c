
#include "../include/micromegas.h"
#include "../include/micromegas_aux.h"
#define P_NAME_SIZE 11

void  all22procList(void)
{  

  char *all=NULL;
  for(int n=0;n<nModelParticles;n++)
  {   int l=0;
      if(all) l=strlen(all);
      all=realloc(all,l+2+strlen(ModelPrtcls[n].name));
      sprintf(all+l,",%s",ModelPrtcls[n].name);
      if(strcmp(ModelPrtcls[n].name,ModelPrtcls[n].aname))
      { l=strlen(all);
        all=realloc(all,l+2+strlen(ModelPrtcls[n].aname));
        sprintf(all+l,",%s",ModelPrtcls[n].aname);
      }    
  }       

//printf("all=%s\n",all); 
  char * command=malloc(strlen(all)+ strlen(compDir)+strlen(calchepDir) +strlen(libDir)   + 200);

  int delWorkDir=prepareWorkPlace();
 
  sprintf(command,"cd %s;"
           " %s/bin/s_calchep -blind \"{{allP,allP->2*x{%s{{{[[{0\" >/dev/null;"
           " if(test $? -eq 0) then mv results/list_prc.txt %s/22ProcList.txt; fi",
            compDir,calchepDir,all+1,libDir);
  system(command); 
//printf("command=%s\n", command); 

  sprintf(command,"mv %s/22ProcList.txt .",libDir);
  system(command); 
  free(command);
  free(all);
  if(delWorkDir) cleanWorkPlace();
  {
    char fname0[50];   
    sprintf(fname0,"22ProcList.txt"); 
    FILE*f1=fopen(fname0,"r");
    sprintf(fname0,"22ShortProcList.txt");
    FILE*f2=fopen(fname0,"w");
    
    if(f1)
    { char p[4][P_NAME_SIZE];
      while(4==fscanf(f1," %[^ ,] , %s -> %[^ ,] , %s",p[0],p[1],p[2],p[3]))
      {  int i,err;  
         int n[4],n_[4];
         for(i=0;i<4;i++)
         {  n[i]=pTabPos(p[i]);
            if(strcmp(ModelPrtcls[abs(n[i])-1].name, ModelPrtcls[abs(n[i])-1].aname)==0) n_[i]=n[i];else n_[i]=-n[i];
         }

         int nn;
         nn=n[0];n[0]=n_[0];n_[0]=nn;
         nn=n[1];n[1]=n_[1];n_[1]=nn;
         if(n[2]>n[3]){ nn=n[2];n[2]=n[3];n[3]=nn;}  
       
         if(n[0]<=n[1] && n[1] <= n[2]&& n[0]<=n[3] && n[1] <= n[2] && n[1]<=n[3])
         {  int i;         
            for(i=0;i<3;)
            { if(n_[i]>n_[i+1]) 
              { nn=n_[i];n_[i]=n_[i+1]; n_[i+1]=nn;
                if(i==0) i++; else i--;
              } else i++;  
            } 
                     
            for(i=0;i<4;i++) if(n[i]!=n_[i]) break;
            
            if(i>4 || n[i]>n_[i]) fprintf(f2,"%s,%s -> %s,%s\n", p[0],p[1],p[2],p[3]);
         }
      }
    }
    fclose(f1);
    fclose(f2);          
  }  
}  
  
