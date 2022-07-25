#include "../include/micromegas_aux.h"

int solveLinEq( int N,double *A, double *c)
{
  int *index=(int*) malloc(N*sizeof(int));
  int i,j,k;
  for(i=0;i<N;i++) index[i]=i;
  
  for(k=0;k<N;k++)
  {  double max=fabs(A[index[k]*N+k]);
     int i0=k;

     for(i=k+1;i<N;i++) 
     {  double  m=fabs(A[index[i]*N+k]);
        if(m>max) { i0=i; max=m;}
     } 
     if(max==0){ free(index); return k+1;}
     if(i0!=k) { int mem=index[k]; index[k]=index[i0];index[i0]=mem;}
     
     for(i=k+1;i<N;i++) 
     { double b=A[index[i]*N+k]/A[index[k]*N+k];
       for(j=k;j<N;j++) A[index[i]*N+j]-=b*A[index[k]*N+j];
       c[index[i]]-=b*c[index[k]];
     }
  }

  for(k=N-1;k>=0;k--)
  {  
    for(j=N-1;j>k;j--) c[index[k]]-=A[index[k]*N+j]*c[index[j]];
    c[index[k]]/=A[index[k]*N+k];
  }
                               
  for(k=0;k<N-1;k++) if(k!=index[k])
  {  
    i=index[k];
    double ck=c[i];
    c[i]=c[k];
    c[k]=ck;  
    for(j=k+1;index[j]!=k;j++) continue;
    index[j]=i;
    index[k]=k;
  }
                             
  free(index);
  return 0; 
}

//#include <math.h>
//#include <quadmath.h>


double  detA( int N,long double  *A)
{


  int *index=(int*) malloc(N*sizeof(int));
  int i,j,k;
  for(i=0;i<N;i++) index[i]=i;
/*
printf("==================\n");
  printf("det(1)=%E\n", (double) ( A[0]*A[3]-A[1]*A[2]));
  for(int i=0;i<N;i++)
  {  for(int j=0;j<N;j++) printf("%E ",(double) A[index[i]*N+j]);
     printf("\n");
  }     
printf("---------------------\n");
*/
  
  int sperm=0;  
  for(k=0;k<N;k++)
  {  long double max=fabs(A[index[k]*N+k]);
     int i0=k;

     for(i=k+1;i<N;i++) 
     {  long double  m=fabs(A[index[i]*N+k]);
        if(m>max) { i0=i; max=m;}
     }
      
     if(max==0){ free(index); return 0;}
     if(i0!=k) { sperm+=abs(index[k]-index[i0]);   int mem=index[k]; index[k]=index[i0];index[i0]=mem;   }
     
     for(i=k+1;i<N;i++) 
     { long double b=A[index[i]*N+k]/A[index[k]*N+k];
//       for(j=k;j<N;j++) A[index[i]*N+j]-=b*A[index[k]*N+j];
         for(j=0;j<N;j++) A[index[i]*N+j]-=b*A[index[k]*N+j];
     }
  }

/*
  for(int i=0;i<N;i++)
  {  for(int j=0;j<N;j++) printf("%E ",(double) A[index[i]*N+j]);
     printf("\n");
  }
*/      
  double ret=1;
  for(int i=0;i<N;i++) ret*=A[index[i]*N+i]; 
  sperm=1;
  for(int i=0;i<N-1;)
  if(index[i]>index[i+1]) 
  { int m=index[i]; index[i]=index[i+1];index[i+1]=m; sperm*=-1;
    if(i)i--;
  } else i++;      
  
  
  
   ret*=sperm;                             

//  printf("det(2)=%e\n",ret);

//printf("=================\n");

  free(index);
  return ret; 
}

