#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "polint.h"

static int  leftXN(int n,int dim, const  double * xa, const double x)
{  int k1,k2,k3;

   k1=n/2;                         
   k2=dim-(n+1)/2-1;

   if(xa[0]< xa[dim-1])
   {              
     if(x<=xa[k1]) return 0;
     if(x>=xa[k2]) return dim-n;                   
     while(k2-k1>1)                
     { k3=(k1+k2)/2;               
       if(xa[k3]>x)k2=k3; else k1=k3;
     }
   } else 
   {  
     if(x>=xa[k1]) return 0;
     if(x<=xa[k2]) return dim-n;
     while(k2-k1>1)                
     { k3=(k1+k2)/2;               
       if(xa[k3]<x)k2=k3; else k1=k3;
     }
   }
   return k1+1-n/2;
}

double  polintN(double x, int n,  const double *xa, const double *ya)
{  double z[20];
   int i,m;
   for(i=0;i<n;i++) z[i]=ya[i];
   for(m=1;m<n;m++) for(i=0;i<n-m;i++)
   z[i]=(z[i]*(xa[i+m]-x) - z[i+1]*(xa[i]-x))/(xa[i+m]-xa[i]);
   return z[0];
}


double polint3(double x, int dim, const  double *xa, const double *ya)
{ 
  if(dim<4) return  polintN(x,dim,xa, ya);
  if((xa[0]<xa[dim-1] && x<=xa[1]) || 
     (xa[0]>xa[dim-1] && x>=xa[1])       ) return  polintN(x,3,xa, ya);

  if((xa[0]<xa[dim-1] && x>=xa[dim-2])|| 
     (xa[0]>xa[dim-1] && x<=xa[dim-2])   ) return  polintN(x,3,xa+dim-3, ya+dim-3);

  int shift=leftXN(4,dim, xa, x);
  if(!isfinite(ya[shift])) return polintN(x,3,xa+shift+1,ya+shift+1);
  if(!isfinite(ya[shift+3])) return polintN(x,3,xa+shift,ya+shift);     
  double x1=xa[shift],x2=xa[shift+1],x3=xa[shift+2],x4=xa[shift+3];
  double y1=ya[shift],y2=ya[shift+1],y3=ya[shift+2],y4=ya[shift+3];
  double h=x3-x2;
  double res=  y3*(x-x2)/h - y2*(x-x3)/h;
  double d=(y3-y2)/h, 
  d2=y1*(x2-x3)/(x1-x2)/(x1-x3)+y2*(2*x2-x3-x1)/(x2-x1)/(x2-x3)+y3*(x2-x1)/(x3-x1)/(x3-x2),
  d3=y2*(x3-x4)/(x2-x3)/(x2-x4)+y3*(2*x3-x4-x2)/(x3-x2)/(x3-x4)+y4*(x3-x2)/(x4-x2)/(x4-x3);

  return res+(x-x2)*(x-x3)*((d2-d)*(x-x3)+(d3-d)*(x-x2))/h/h;
}



double polint1(double x, int n,  double *xa, double *ya)
{ int shift=leftXN(2,n, xa, x);
   return polintN(x,2,xa+shift, ya+shift);
}

double polint1Exp(double x, int n,  double *xa, double *ya)
{ int shift=leftXN(2,n, xa, x);
  double  x1=xa[shift], x2=xa[shift+1],y1=ya[shift],y2=ya[shift+1];
  double alpha= (x-x1)/(x2-x1);
  
  if(y1>0 && y2>0) return  pow(y1,1-alpha)*pow(y2,alpha);
  return  y1*(1-alpha)+y2*alpha;  
}

double polintDiff(int n, const  double *xa, const double *ya, double * dxdy)
{
  for(int i=0;i<n;i++)
  { double x0=xa[i], y0=ya[i], x1,x2,y1,y2;
    if(i==0)       { x1=xa[1];  y1=ya[1];  x2=xa[2];  y2=ya[2];}
    else if(i==n-1){ x1=xa[n-2];y1=ya[n-2];x2=xa[n-3];y2=ya[n-3];}
    else           { x1=xa[i-1];y1=ya[i-1];x2=xa[i+1];y2=ya[i+1];}
    
    dxdy[i]=y0*(2*x0-(x1+x2))/(x0-x1)/(x0-x2)
           +y1*(x0-x2)       /(x1-x0)/(x1-x2)
           +y2*(x0-x1)       /(x2-x0)/(x2-x1);
  }                   
}    
