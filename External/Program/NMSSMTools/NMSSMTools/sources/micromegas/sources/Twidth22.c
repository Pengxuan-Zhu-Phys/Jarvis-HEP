#include<math.h>
#include<stdlib.h>
#include<stdio.h>

#include"micromegas.h"
#include"micromegas_aux.h" 
#include "../CalcHEP_src/include/VandP.h"

//#include"../CalcHEP_src/c_source/num/include/alphas2.h"

#include"../CalcHEP_src/c_source/dynamicME/include/dynamic_cs.h"

//#include "vp.h"


double tWidth21(char *name, double T, int show)
{
  REAL pvect[12];
  double width=0.;
  REAL m[3]; 
  int i,ntot,nsub,nin,nout,pIn,pOut,err_code;
  double GG;
  char proc[20];
  if(T==0) return 0;

  if(show) printf("Thermal width %s (T=%.2E GeV)\n", name);  
  double tWidth=0;

  sprintf(proc,"%s->2*x",name);
  numout *cc=newProcess(proc);
  if(cc)
  { 
     procInfo1(cc,&ntot,&nin,&nout);                                                                                                        
     if(passParameters(cc)) return -1;

     for(nsub=1;nsub<=ntot;nsub++)
     {
       REAL m[3];
       int ndf[3],spin2[3];
       char*p[3];

       for(int i=0;i<3;i++)
       {
         p[i]=cc->interface->pinf(nsub,i+1,m+i,NULL);
         cc->interface->pinfAux(nsub, i+1,spin2+i,NULL,NULL,ndf+i);
       }
       if(isFeeble(p[1]) || isFeeble(p[2])) continue;

       if(m[0]==0) continue;

            if(m[1]>m[0]+m[2]) { pIn=1; pOut=2; }
       else if(m[2]>m[0]+m[1]) { pIn=2; pOut=1; }
       else continue;

       if( m[pIn]-m[0] > 20*T) continue;

       double md=m[0]-m[pOut];
       double ms=m[0]+m[pOut];
       double pRestOut=Sqrt((m[pIn]*m[pIn] - ms*ms)*(m[pIn]*m[pIn]-md*md))/(2*m[pIn]);
       double totcoef= pRestOut/(8. * M_PI * m[pIn]*m[pIn]);

       for(i=1;i<12;i++) pvect[i]=0;

       pvect[0]=Sqrt(pRestOut*pRestOut+m[0]*m[0]);
       pvect[1]=pRestOut;

       pvect[4*pIn]=m[pIn];
       pvect[4*pOut]=Sqrt(pRestOut*pRestOut+m[pOut]*m[pOut]);
       pvect[4*pOut+1]=-pRestOut;

       GG=sqrt(4*M_PI*alphaQCD(m[pIn]));

       width = totcoef * (cc->interface->sqme)(nsub,GG,pvect,NULL,&err_code);
       if((spin2[0]&1)!=(spin2[pIn]&1)) width*=-1;
       double tf=  ndf[1]*m[pIn]*m[pIn]/m[0]/m[0];
       if(m[0]!=0)
       {
         if(m[pIn]/T > 20 ) tf*=exp(-m[pIn]/T+m[0]/T)*sqrt(m[0]/m[pIn]);   
         else 
         {  
           if(m[pIn]/T<0.01) tf/=m[pIn]/T; else tf*=bessK1(m[pIn]/T);
           if(m[0]/T<0.01)   tf*=m[0]/T  ; else tf/=bessK1(m[0]/T);
         }
         
//         printf("width(%s -> %s,%s)=%e  ft=%E tWidth=%E %E %E \n", p[pIn],p[0],p[pOut],width*ndf[0]/(double)ndf[pIn],tf,width*tf,m[pIn]/T,m[0]/T );

         width*=tf;
         if(!isfinite(width))
         { printf("tWidth(%s -> %s,%s T=%e)  is NaN\n", p[pIn],p[0],p[pOut],T);
          exit(0);
         } 
       }
       
       if(show)  printf(" %s %s -> %s %.2E\n",name, antiParticle(p[pIn]), p[pOut], width);
       tWidth+=width;
     }  
  }
  
  if(isFeeble(name)) return tWidth;
  char * ap=antiParticle(name);
  sprintf(proc,"%s,%s->1*x",name,ap);
  cc=newProcess(proc);
  if(cc)
  {
     procInfo1(cc,&ntot,&nin,&nout);  
     passParameters(cc);
     for(nsub=1;nsub<=ntot;nsub++)
     {
       REAL m[3];
       int ndf[3],spin2[3];
       char*p[3];

       for(int i=0;i<3;i++)
       {
         p[i]=cc->interface->pinf(nsub,i+1,m+i,NULL);
         cc->interface->pinfAux(nsub, i+1,spin2+i,NULL,NULL,ndf+i);
       }

       if(m[2]<=2*m[0]) continue;

       if( m[2]-m[0] > 20*T) continue;

       double md=0;
       double ms=2*m[1];
       double pRestOut=Sqrt(m[2]*m[2] - ms*ms)/2;
       double totcoef= pRestOut/(8. * M_PI * m[2]*m[2]);

       for(i=1;i<12;i++) pvect[i]=0;

       pvect[8]=m[2];

       pvect[0]=Sqrt(m[0]*m[0]+pRestOut*pRestOut);
       pvect[1]=pRestOut;
       
       pvect[4]=pvect[0];
       pvect[5]=-pRestOut;

       GG=sqrt(4*M_PI*alphaQCD(m[0]));

       width = totcoef * (cc->interface->sqme)(nsub,GG,pvect,NULL,&err_code);
       double tf=ndf[1]*m[2]*m[2]/m[0]/m[0]*bessK1(m[2]/T)/bessK1(m[0]/T);
       if(p[1]==p[2]) tf*=2;
       
       if(show) printf(" %s %s -> %s %.2E\n",p[0],p[1],p[2]);

       width*=tf;

       tWidth+=width;
     }       
  }
  if(show) printf("  sum=%.2E GeV\n",tWidth);
  return tWidth;
  
}
