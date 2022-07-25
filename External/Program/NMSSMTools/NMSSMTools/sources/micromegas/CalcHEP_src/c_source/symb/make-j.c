#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <pthread.h>


static pthread_mutex_t keyN=PTHREAD_MUTEX_INITIALIZER;

typedef struct fList { struct fList* next; char fName[30];} fList;

static fList * list=NULL;


static void * pCompile_cycle(void * input)
{  
  char command[100];
  for(;;)
  { 
    char*fName=NULL;   
    pthread_mutex_lock(&keyN);
     if(list) 
     { fName=list->fName;
       list=list->next;
     }    
    pthread_mutex_unlock(&keyN);
    if(fName)
    {  sprintf(command," $CC -c $CFLAGS %s; if(test $? -eq 0) then  rm %s;fi\n",fName,fName);     
       system(command);
    } else return NULL;  
  } 
  return NULL;
}

int main(int argc, char** argv)
{  
   DIR *dirPtr=NULL;
   struct dirent * dp=NULL;
   int Ntot=0;
   static int nParProc;

   if(argc!=2) return 1;
   if(1!=sscanf(argv[1],"%d",&nParProc)) return 2;
   if(nParProc<=0) nParProc=1;
   
   dirPtr=opendir("./");   
   while((dp=readdir(dirPtr)))
   { char *c=strstr(dp->d_name,".c");
     if((dp->d_name[0]!='.') && c&&c[2]==0) { fList *tmp=malloc(sizeof(fList));
                      tmp->next=list;
                      strcpy(tmp->fName,dp->d_name);
                      list=tmp; 
                       Ntot++;
                    }   
   }  
   closedir(dirPtr);   
   if(!Ntot) return 0;
   if(nParProc>Ntot) nParProc=Ntot;
   { pthread_t *threads = malloc(nParProc*sizeof(pthread_t));
     int k;  
     for (k=0;k<nParProc;k++) pthread_create(threads+k,NULL,pCompile_cycle,NULL);  
     for (k=0;k<nParProc;k++) pthread_join(threads[k],NULL);
     free(threads);
   }
   return 0;   
}
