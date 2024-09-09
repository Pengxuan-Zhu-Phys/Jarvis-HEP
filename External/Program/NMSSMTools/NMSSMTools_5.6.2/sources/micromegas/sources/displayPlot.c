#include"micromegas.h"
#include"micromegas_aux.h"
#include"../CalcHEP_src/c_source/chep_crt/include/crt.h"
#include"../CalcHEP_src/include/rootDir.h"
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>              
#include <sys/wait.h>
#include <stdarg.h>

 
  
static int First=1;
  
static void disconnect(int N) { setsid();}
  

extern void  plot_1(double xMin, double xMax, int dim,
                       double *f, double *ff,char* upstr, 
                       char* xstr, char* ystr);
extern int blind;

  
static int newPID=0;
static int pidList[100];

extern char pathtocalchep[], pathtohelp[];
  
void displayPlotN(char * title,char*xName,  double xMin, double xMax,  int lScale, int N, char**Y, int*Dim, double**f,double**ff)
{ int pid;
  
  if(First) { First=0;   signal(SIGUSR1, disconnect);}  

  fflush(NULL);
  pid=fork();
  if(pid==0) 
  {  int err,i;
     va_list ap;   
     
     blind=0; 
     err=start1("micrOMEGAs Plot",NULL ,"calchep.ini",NULL);
     if(err) 
     { printf("Can not display plot because micromegas is compiled without X11\n");
       exit(0);
     }
     sprintf(pathtocalchep,"%s/",rootDir);
     sprintf(pathtohelp,"%s/help/",pathtocalchep);
     clearTypeAhead();     
     plot_Nar(NULL,title,xName,xMin,xMax, lScale,N, Dim,f,ff,Y);
     finish();
     exit(0);
  } else pidList[newPID++]=pid;
}

int displayPlot(char * title, char*xName, double xMin, double xMax ,  int lScale, int N, ...)
{
  int i,j;
  double **f; double**ff; char**Y;
  int *Dim;
  int * delMark;
  va_list ap;   
  if(lScale !=0 && lScale!=1)
  { printf(" 5-th parameter  is not 0/for linear scale/ or 1/for log scale/.\n");
    return 1;
  }   

  Dim=(int*)malloc(N*sizeof(int));
  f =malloc(N*sizeof(double*));
  ff=malloc(N*sizeof(double*));
  Y = malloc(N*sizeof(char*));
  delMark=(int*)malloc(N*sizeof(int));
  va_start(ap,N);
  for(i=0;i<N;i++) 
  { 
    Y[i]=va_arg(ap,char*);
    for(j=0;j<50;j++) if(Y[i][j]==0)break;
    if(j==50) 
    {  printf("Parameter %d is not a text string or too long. Curve title expected in this position.\n",7+N*4);
       free(Dim); free(f); free(ff); free(Y); free(delMark);
       return 2;
    }    
    Dim[i]=va_arg(ap,int);
    if(Dim[i]<0 || Dim[i]>300) 
    { printf("Parameter %d =%d should present number of points for curve of number of bins for histogram.\n" 
      "Values larger than 300 are forbidden. Zero input reserved for presentation of  functions\n",4+N*4,Dim[i] );
      free(Dim); free(f); free(ff); free(Y); free(delMark);
      return 3;
    } 
    if(Dim[i])
    { 
      f[i]=va_arg(ap,double*);
      ff[i]=va_arg(ap,double*);
      delMark[i]=0;
    } else 
    { void * F,*A;
      Dim[i]=150;
      f[i]=malloc(Dim[i]*sizeof(double));
      ff[i]=NULL;
      delMark[i]=1;
      F=va_arg(ap,void*);
      A=va_arg(ap,void*);
         
      for(j=0;j<Dim[i];j++)
      { double x;
        if(lScale) x=xMin*pow(xMax/xMin,(j+0.5)/Dim[i]); else x= xMin+(j+0.5)/Dim[i]*(xMax-xMin);
        if(A) { double (*func)(double,void*)=F; f[i][j]=func(x,A);} 
        else  { double (*func)(double)=F;       f[i][j]=func(x);}
      }
                 
    }     
  }   
  va_end(ap);
  
  displayPlotN(title, xName,xMin,xMax,lScale,N,Y,Dim,f,ff);
  for(i=0;i<N;i++) if(delMark[i]) free(f[i]);
  free(Dim); free(f); free(ff); free(Y); free(delMark);
}


void  killPlots(void)
{
  int  C,i,pid;

  if(!newPID) return;
  for(i=0;i<newPID;i++) if(waitpid(pidList[i],NULL,WNOHANG)==0) break; else pidList[i]=0;
  
  if(i!=newPID)
  { int id=fork();
    if(id==0)
    {  printf("Kill all plots (Y/N)? "); 
       C=getchar();
       if(C=='y'|| C=='Y')   exit(0); else exit(1); 
    } else for(;;) 
    {  sleep(1);
       int status;
       for(i=0;i<newPID;i++ ) if(pidList[i]){ if(waitpid(pidList[i],NULL,WNOHANG)==0) break; else pidList[i]=0; }
       if(i==newPID)  kill(0,SIGKILL);
       else  if( waitpid(id,&status,WNOHANG))
       {  if( !WIFEXITED(status) ||  WEXITSTATUS(status)==0)   kill(0,SIGKILL); else kill(0,SIGUSR1);
          break; 
       } 
    }
  }  
  newPID=0;
}



void  killplots_(void) { killPlots();}

