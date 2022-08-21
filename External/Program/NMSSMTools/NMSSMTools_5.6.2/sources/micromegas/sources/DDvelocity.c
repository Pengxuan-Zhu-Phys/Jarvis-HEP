#include"micromegas_aux.h"
#include"micromegas_f.h"


/*===== Maxwell velosity distribution =====*/ 

double Maxwell(double v) 
{  double res,vsum,vdif,DV2;

   if(v<=0|| v>=vEsc+vEarth) return 0;

   double rv=vEsc/vRot;
   double NR= pow(M_PI,1.5)*vRot*vRot*vRot*(erf(rv)- 2*rv/sqrt(M_PI)*exp(-rv*rv));
             
   DV2=vRot*vRot;
   if(vEarth*v<0.001*DV2) res= 4*v*exp((-v*v-vEarth*vEarth)/(DV2))/(DV2);
   else 
   {   
     vsum=vEarth+v;
     if(vsum>vEsc) vsum=vEsc;
     vdif=vEarth-v;
     res=(exp(-vdif*vdif/DV2)-exp(-vsum*vsum/DV2))/vEarth;
   }
   return v*M_PI*vRot*vRot/NR*res;
}

static double erfi_int(double x){ return exp(x*x);} 
static double erfi( double x){ return 2/sqrt(M_PI)*simpson(erfi_int,0,x,1E-4,NULL);}

static double vR,vT,v_,sn; // integration parameters 
static double fiInt(double fi) { return exp( -pow(v_*sn*sin(fi)/vR,2) - pow(v_*sn*cos(fi)/vT,2));} // azimuth angle integrand 
                                          
static double cosInt(double cs) 
{  sn=sqrt(1-cs*cs);  return  simpson( fiInt,0,2*M_PI,1E-5,NULL)*exp(-pow((v_*cs-vEarth)/vT,2));} // polar anfle intergand
   

#define DIM 100

double SHMpp(double v) 
{ 
  static double vTab[DIM], SHMppTab[DIM] ;
  if(v<0 || v>vEarth+vEsc) return 0;
  double static beta=-1, eta=-1, vRo=-1,vEa=-1,vEs=-1;
  if( beta!=betaSHMpp || eta!=etaSHMpp ||   vEs!=vEsc || vRo!=vRot || vEa!=vEarth)
  {  beta=betaSHMpp;
     eta=etaSHMpp;
     vEs=vEsc;
     vRo=vRot;
     vEa=vEarth; 
     for(int i=0;i<DIM;i++) vTab[i]=i*(vEarth+vEsc)/(DIM-1);
     SHMppTab[0]=SHMppTab[DIM-1]=0;  
     double rv=vEsc/vRot;
     double NR= pow(M_PI,1.5)*vRot*vRot*vRot*(erf(rv)-     2*rv/sqrt(M_PI)*exp(-rv*rv));
     vR=vRot/sqrt(1-2./3.*beta), vT=vRot*sqrt((1-beta)/(1-2./3.*beta));
     double NS= pow(M_PI,1.5)*vR*vT*vT*(erf(vEsc/vR) - sqrt( (1/beta-1))*exp(-vEsc*vEsc/vT/vT)*erfi(vEsc/vR/sqrt(1/beta-1)));

     for(int i=1;i<DIM-1;i++)
     {  v_=vTab[i];  
       double vm= vEsc>v_+vEarth ? v_+vEarth : vEsc;
       double ediff; 
       if(vEarth*v_<0.001*vRot*vRot) ediff= 4*v*exp((-v_*v_-vEarth*vEarth)/vRot/vRot);
       else  ediff=exp(-pow((v_-vEarth)/vRot,2))-exp(-pow(vm/vRot,2));       
       SHMppTab[i]= (1-etaSHMpp)*M_PI*vRot*vRot/vEarth/NR*ediff;
              vm= vEsc<v_+vEarth ? v_+vEarth : vEsc;
       SHMppTab[i]+= etaSHMpp* v_*simpson(cosInt,-1+(vm*vm-vEsc*vEsc)/(2*vEarth*v_),1,1E-5,NULL)/NS;       
     }
  } 
  return v*polint3(v,DIM,vTab, SHMppTab);
}  
  