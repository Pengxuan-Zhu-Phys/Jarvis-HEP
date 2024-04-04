#include "../include/micromegas_aux.h"
#define smallDelta 1E-10
//#define ERROR_PLOT

int init22_par(par_22*arg, numout*cc, int nsub)
{  int i;
   arg->cc=cc;
   arg->nsub=nsub;
   int ntot,nin,nout; 
   if(!cc) return 1;
   procInfo1(cc,&ntot,&nin,&nout);
   if(nsub>ntot || nin!=2 || nout!=2) return 2;    
   for(i=0;i<4;i++) arg->pdg[i]=arg->cc->interface->pinfAux(arg->nsub,1+i,arg->spin2+i,NULL,NULL,(arg->ndf)+i);
//   printf("arg->ndf = %d %d %d %d\n", arg->ndf[0], arg->ndf[1], arg->ndf[2], arg->ndf[3]);
   arg->T=0;     
   arg->err=0;
   return 0;
}


static double statFactor(par_22*arg) // for Monte Carlo
{ double c=1;
  for(int i=0;i<4;i++)if(arg->eta[i])
  { double p=0; for(int j=0;j<3;j++) p+=arg->pvect[j+1+4*i]*arg->n[j]; 
    double E=arg->pvect[4*i]*arg->ch - p*arg->sh;
    c/=1-arg->eta[i]*exp(-E/arg->T); 
  } return c;  
} 


void mass22_par(par_22*arg,double T)
{
   int i;
   numout*cc=arg->cc; 
   int nsub=arg->nsub;
   char * p[5];
passParameters(cc);   
   for(i=0;i<4;i++) p[i]=cc->interface->pinf(nsub,1+i,arg->pmass+i,NULL);
     
   
   arg->s13=arg->s14=0;
   for(int n=1;;n++)
   {  
     int m,w,pnum;
     char*s=cc->interface->den_info(nsub,n,&m,&w,&pnum);
     if(!s) break;
     if(!m || !w) continue;
     if(s[0]==1 && s[1]==3)
     { 
       double mass=cc->interface->va[m];
       if(arg->s13==0)
       { arg->s13=1- 2*((arg->spin2[0]+arg->spin2[2])&1);
         arg->M13=mass;  
       } else if(arg->M13>mass) arg->M13=mass;
       char*p=ModelPrtcls[pnum].name; 
       cc->interface->va[w]=pWidth(p,NULL)+tWidth21(p,T,0);
//printf("width(%s,%E)=%E+%E\n", p, T, pWidth(p,NULL),tWidth21(p,T,0)); 
             
       if( (   arg->pmass[0] > mass + arg->pmass[2] 
            && arg->pmass[3] > mass + arg->pmass[1]
           ) ||
           (   arg->pmass[1] > mass + arg->pmass[3]
            && arg->pmass[2] > mass + arg->pmass[0] 
           ) )  cc->interface->va[w]*=3;
       
     }
     if(s[0]==1 && s[1]==4)
     { 
       double mass=cc->interface->va[m];
       if(arg->s14==0)
       { arg->s14=1- 2*((arg->spin2[0]+arg->spin2[3])&1);
         arg->M14=mass;  
       } else if(arg->M14>mass) arg->M14=mass;
        char*p=ModelPrtcls[pnum].name;
        cc->interface->va[w]=pWidth(p,NULL)+tWidth21(p,T,0);
       if( (   arg->pmass[1] > mass + arg->pmass[2] 
            && arg->pmass[3] > mass + arg->pmass[0]
           ) ||
           (   arg->pmass[0] > mass + arg->pmass[3]
            && arg->pmass[2] > mass + arg->pmass[1] 
           ) )  cc->interface->va[w]*=3; 
     }
   }
   arg->sqrtSmin=arg->pmass[0]+arg->pmass[1];
   if(arg->sqrtSmin < arg->pmass[2]+arg->pmass[3]) arg->sqrtSmin=arg->pmass[2]+arg->pmass[3]; 
}


int  kin22_par(par_22*arg, REAL sqrtS,double GG)
{  
   REAL ms,md,lambda,Pin,Pout,m[5];
   int i;
   numout*cc=arg->cc;
   arg->err-=arg->err&ErrE;

   int nsub=arg->nsub;
//printf("kin22_par\n");   
   for(i=0;i<4;i++) m[i+1]=arg->pmass[i];
//printf("masses again %E %E %E %E\n", m[1],m[2],m[3],m[4]);
   ms= m[1]+m[2];
   md= m[1]-m[2];
   if(sqrtS<=ms) { arg->err+=ErrE; return 1;} 
   lambda = Sqrt((sqrtS-ms)*(sqrtS+ms)*(sqrtS-md)*(sqrtS+md));
   Pin=arg->PcmIn = lambda/(2*sqrtS);

   ms= m[3]+m[4];
   md= m[3]-m[4];
   if(sqrtS<=ms)  { arg->err+=ErrE; return 1;} 
   lambda = Sqrt((sqrtS-ms)*(sqrtS+ms)*(sqrtS-md)*(sqrtS+md));
   Pout=arg->PcmOut = lambda/(2*sqrtS);

   for(i=0;i<16;i++) arg->pvect[i]=0;   
   
   for(i=1;i<5;i++) m[i]*=m[i];
   
   arg->pvect[3] = Pin;
   arg->pvect[7] =-Pin;
   arg->pvect[0] = Sqrt(Pin*Pin + m[1]);
   arg->pvect[4] = Sqrt(Pin*Pin + m[2]);
      
   arg->pvect[8] = Sqrt(Pout*Pout + m[3]);
   arg->pvect[12]= Sqrt(Pout*Pout + m[4]);


   double dE,dP;
   if(arg->s13)
   { dE=(m[1]-m[3]-m[2]+m[4])/2/sqrtS; 
     dE=arg->pvect[0]-arg->pvect[8];
     dP=Pin-Pout; 
     arg->e13=(arg->M13*arg->M13+(dP-dE)*(dP+dE))/2/Pin/Pout;
   }
   if(arg->s14)
   {  
     dE=(m[1]-m[4]-m[2]+m[3])/2/sqrtS;
     dP=Pin-Pout;
     arg->e14=(arg->M14*arg->M14-dE*dE+dP*dP )/2/Pin/Pout;
   }    
   arg->totFactor=Pout /(32.0*M_PI*Pin*sqrtS*sqrtS);
   arg->GG=GG;
   return 0;
}


static double  dSqme_dCosR_arg(double  dCos,void*arg_)
{
   double  r;
   par_22*arg=arg_;
   REAL cosf=(REAL)1 -(REAL)dCos;
   REAL sinf=Sqrt(Fabs((REAL)dCos*(1+cosf)));
   int err_code=0;

//printf("arg->PcmOut=%E\n", arg->PcmOut);   
   arg->pvect[11]= arg->PcmOut*cosf;
   arg->pvect[15]=-arg->pvect[11];
   arg->pvect[10]= arg->PcmOut*sinf;
   arg->pvect[14]=-arg->pvect[10];

   r = arg->cc->interface->sqme(arg->nsub,arg->GG,arg->pvect,NULL,&err_code); 
   if(arg->T>0) r*=statFactor(arg);
   if(err_code)  arg->err= arg->err | err_code ;
   return r;
}

static double  dSqme_dCosL_arg(double  dCos,void*arg_)
{
   double  r;
   par_22*arg=arg_;
   REAL cosf=-(REAL)1 +(REAL)dCos;
   REAL sinf=Sqrt(Fabs((REAL)dCos*(1-cosf)));
   int err_code=0;

//printf("arg->PcmOut=%E\n", arg->PcmOut);   
   arg->pvect[11]= arg->PcmOut*cosf;
   arg->pvect[15]=-arg->pvect[11];
   arg->pvect[10]= arg->PcmOut*sinf;
   arg->pvect[14]=-arg->pvect[10];

   r = arg->cc->interface->sqme(arg->nsub,arg->GG,arg->pvect,NULL,&err_code);
   if(arg->T>0) r*=statFactor(arg);
   if(err_code)  arg->err= arg->err | err_code ;
   return r;
}


static double  poleSqme(double x0,double x1,double x2,double delta,double f0, double f1,double f2)
{
  if(f2-f0<=0) {  return f0*x0;}    
  double r=(f1-f0)/(f2-f0),z1=(x0+delta)/(x1+delta),z2=(x0+delta)/(x2+delta);
  double p,a,b;
  if(r>=(1/z1-1)/(1/z2-1)) p=-1;
  else if(r<=(z1*z1-1)/(z2*z2-1)) p=2;
  else 
  {  double p1=-1,p2=2;
     double a1=(pow(z1,p1)-1)/(pow(z2,p1)-1)-r , a2=(pow(z1,p2)-1)/(pow(z2,p2)-1)-r;
     if(a1>a2) { a=a1;a1=a2;a2=a; p=p1;p1=p2;p2=p;}

// printf("solve a1=%E a2=%E\n",a1,a2);
     
     while(fabs(p1-p2)>0.001) 
     { p=0.5*(p1+p2);
       a=(pow(z1,p)-1)/(pow(z2,p)-1)-r; 
       if(a<0) { a1=a;p1=p;} else { a2=a;p2=p;} 
     }   
  }   
//p=1;
//printf("p=%E\n", p);
  
  if(fabs(1-p)<0.1)
  {  p=1;
     b=(f0-f1)/(1/(x0+delta)-1/(x1+delta));
     a=f0-b/(x0+delta);
//printf("f0=%e f1=%E x0=%e x1=%e\n", f0,f1,x0,x1);     
//printf("dI=%E a=%E b=%e p=%e\n",a*x0+ b*log(x0/delta+1),a,b,p);
     return   a*x0+ b*log(x0/delta+1);
  } 
                
  b=(f0-f1)/(pow(x0+delta,-p)-pow(x1+delta,-p)); 
  a=f0-b*pow(x1+delta,-p);
//printf("dI=%E a=%E b=%e p=%e\n", a*x0+b/(1-p)*(pow(x0+delta,1-p) - pow(delta,1-p)),a,b,p);  
  return a*x0+b/(1-p)*(pow(x0+delta,1-p) - pow(delta,1-p));

}  

// cross section [pb]= sqmeInt*3.8937966E8*arg.PcmOut/(32.0*M_PI*arg.PcmIn*sqrtS*sqrtS))

double sqmeInt(par_22*arg, double eps)
{  
     
  double res1=0;
  int i;
  int err=0;
  double x0=0,f0,f1,f2,in;  
  
  int err_tmp=arg->err;
  if(arg->s13!=0 && arg->e13 < smallDelta)
  { 
     x0=smallDelta;
     for(;;)
     { arg->err=0;
       f0=dSqme_dCosR_arg(x0,arg);
       f1=dSqme_dCosR_arg(x0/2,arg);
       f2=dSqme_dCosR_arg(x0/4,arg);
//?       if(arg->err&ErrC && arg->err&ErrD) { arg->err=err_tmp|ErrP;  return 0;}
       if(x0<0.8) { in=0; break; } // lost of precision in matrix element not caused by t-channel
       if(arg->err&ErrD){x0*=2; continue;}
       if(arg->err&ErrC){x0/=2; continue;}      
       in=poleSqme(x0,x0/2,x0/4,arg->e13,f0,f1,f2);
       break;
     }     
     if(isfinite(in)) res1+=in;
  }
  
  arg->err=err_tmp;
  res1+= simpson_arg(dSqme_dCosR_arg,arg,x0,1,eps,&err); 
  if(err)  arg->err=arg->err|(2*err);
/*
  if(err)
  {
   printf("simpson error code %d. File 22_par.c line 363\n",err);
#ifdef ERROR_PLOT  
    char txt[100];
    sprintf(txt,"Process %s,%s -> %s %s", arg->cc->interface->pinf(1,1,NULL,NULL),
                                          arg->cc->interface->pinf(1,2,NULL,NULL),
                                          arg->cc->interface->pinf(1,3,NULL,NULL),
                                          arg->cc->interface->pinf(1,4,NULL,NULL));
    displayPlot("dSqme_dCos","1-cos",x0,1,0,1,"",0,dSqme_dCosR_arg,arg);
    exit(0);
#endif    
  }   
*/  
  if(arg->pdg[0]==arg->pdg[1]) return 2*res1;

  err_tmp=arg->err;
  x0=0;    
  
  if(arg->s14!=0 && arg->e14< smallDelta )
  {
    x0=smallDelta;
    for(;;)
    { arg->err=0;
      f0=dSqme_dCosL_arg(x0,arg);  
      f1=dSqme_dCosL_arg(x0/2,arg);      
      f2=dSqme_dCosL_arg(x0/4,arg);  
      
//?      if(arg->err&ErrC && arg->err&ErrD) { arg->err=err_tmp|ErrP;  return 0;}
      if(x0<0.8) { in=0; break; } // lost of precision in matrix element not caused by t-channel
      if(arg->err&ErrD){x0*=2; continue;}
      if(arg->err&ErrC){x0/=2; continue;}
                                                
      in=poleSqme(x0,x0/2,x0/4,arg->e14,f0,f1,f2);
      break; 
    }        
    if(isfinite(in))  res1+=in;
  }                                                                     
  arg->err=err_tmp;
  res1+=simpson_arg(dSqme_dCosL_arg,arg,x0,1,eps,&err); 
  if(err) arg->err=arg->err|(2*err);   
/*
  if(err)
  {
   printf("simpson error code %d. File 22_par.c line 401\n",err);
#ifdef ERROR_PLOT  
    char txt[100];
    sprintf(txt,"Process %s,%s -> %s %s", arg->cc->interface->pinf(1,1,NULL,NULL),
                                          arg->cc->interface->pinf(1,2,NULL,NULL),
                                          arg->cc->interface->pinf(1,3,NULL,NULL),
                                          arg->cc->interface->pinf(1,4,NULL,NULL));
    displayPlot("dSqme_dCos","cos+1",x0,1,0,1,"",0,dSqme_dCosL_arg,arg);
    exit(0);
#endif    
  }
*/  
  return res1;
} 


static double  dSqme_dCos_arg(double cs,void*arg_)
{
   double  r;
   par_22*arg=arg_;
   REAL cosf=cs;
   REAL sinf=Sqrt(Fabs((1-cosf)*(1+cosf)));
   int err_code=0;

//printf("arg->PcmOut=%E\n", arg->PcmOut);   
   arg->pvect[11]= arg->PcmOut*cosf;
   arg->pvect[15]=-arg->pvect[11];
   arg->pvect[10]= arg->PcmOut*sinf;
   arg->pvect[14]=-arg->pvect[10];

   r = arg->cc->interface->sqme(arg->nsub,arg->GG,arg->pvect,NULL,&err_code);
  // if(arg->T>0) r*=statFactor(arg);
   if(err_code)  arg->err= arg->err | err_code ;
   return r;
}



void mass22_parDel(par_22*arg,double T)
{
   int i;
   numout*cc=arg->cc; 
   int nsub=arg->nsub;
   char * p[5];
   for(i=0;i<4;i++) p[i]=cc->interface->pinf(nsub,1+i,arg->pmass+i,NULL);
   arg->sqrtSmin=arg->pmass[0]+arg->pmass[1];
   if(arg->sqrtSmin < arg->pmass[2]+arg->pmass[3]) arg->sqrtSmin=arg->pmass[2]+arg->pmass[3]; 
}




double sqmeIntDel(par_22*arg, double eps)
{  
  double delta=0.05;
  
  REAL *mass=arg->pmass;
  REAL E[4];
  for(int i=0;i<4;i++) E[i]=arg->pvect[4*i];
  
  double pp=2*arg->PcmIn*arg->PcmOut;
   
  double t13min =mass[0]*mass[0]+mass[2]*mass[2]-2*E[0]*E[2], t13max =t13min;
  t13min-=pp; t13max+=pp;

  double t14min =mass[0]*mass[0]+mass[3]*mass[3]-2*E[0]*E[3], t14max =t14min;;
  t14min-=pp; t14max+=pp;
  
  double T=arg->T;
  
  double cos0=1-T*T*cTcut*cTcut/pp;
  if(cos0<=0) return 0;

        
  double*intervals=malloc(2*sizeof(double));
  intervals[0]=-cos0;intervals[1]=cos0;
  int nIn=1;
  for(int n=1;;n++)
  {  
     int m,w,pnum;
     char*s=arg->cc->interface->den_info(arg->nsub,n,&m,&w,&pnum);
     if(!s) break;
     if(s[0]==1 && s[1]==3) 
     { double mt=0;
       if(m) mt=fabs(arg->cc->interface->va[m]);
double k=1E-5;
       if((mass[0]>mt+mass[2] && mass[3]>mt+mass[1])  || (mass[1]>mt+mass[3] && mass[2]>mt+mass[0])) k=0.5;
       {
         double co= -1 +2*(t13min - mt*mt )/(t13min-t13max);  
         if(delta>0 && co<1+k*delta && co>-1-k*delta) delInterval(co-k*delta,co+k*delta,&intervals, &nIn);
       }   
     } else if(s[0]==1 && s[1]==4)
     { 
       double mu=0;
       if(m) mu=fabs(arg->cc->interface->va[m]);
double k=1E-5;
       if((mass[0]>mu+mass[3] && mass[2]>mu+mass[1])  || (mass[1]>mu+mass[2] && mass[3]>mu+mass[0])) k=0.5;
       {
          double co= 1 -2*(t14min - mu*mu)/(t14min-t14max);     
          if(delta>0 && co<1+k*delta && co>-1-k*delta) delInterval(co-k*delta,co+k*delta,&intervals, &nIn); 
       }  
     }           
  }

                 
  double sum=0;
  int err;
  for(int i=0;i<nIn;i++)
  { double dI;

//    double dSum=peterson21_arg(dSqme_dCos_arg,arg, intervals[2*i], intervals[2*i+1],&dI);
        
//    if(fabs(dI)> 1E-1*fabs(dSum)) 

    double dSum=simpson_arg(dSqme_dCos_arg,arg, intervals[2*i], intervals[2*i+1],1E-2, &err );

    if(err)    
    {
    
//    displayPlot("angle integrand", "cos", intervals[2*i], intervals[2*i+1], 0,1,"",0, dSqme_dCos_arg,arg);
    
      arg->err=arg->err|8;
                                  }  
    sum+=dSum; 
  }      
  free(intervals);
  if(arg->pdg[0]==arg->pdg[1]) sum*=2;
   if(!isfinite(sum)) arg->err=arg->err|2;
  return sum;
} 
 