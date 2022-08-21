#define DUMMY
#ifdef DUMMY
   //  default dummy routine 
#include "micromegas.h"
#include "micromegas_aux.h"

void improveCrossSection(long n1,long n2,long n3,long n4 ,double PcmIn, 
double * res) { return; }

#else

//  Example for IDM. File should be disposed in IDM/lib
//  Loop improved cross sections depend on input parameters of IDM and can not be used in generic case. 
// It is only example presented here to help the user to construct his own version of loop improved calculations 


#include "../../include/micromegas.h"
#include "../../include/micromegas_aux.h"

 
void improveCrossSection(long n1,long n2,long n3,long n4 ,double PcmIn,  double * res) 
{

//  sorting PDG identifiers 
     if(n1>n2){ long m=n1;n1=n2;n2=m;}
     if(n3>n4){ long m=n3;n3=n4;n4=m;}
 
// For code flexibility it is better to organize reading of improved cross sections from file.
// But here for simplicity we use tabulated data. 
     
//   for ~X,~X->Z,Z       
    if(n1==35 &&  n2 == 35 &&   n3 == 23 &&  n4 == 23)
    { 
       double pcmTab[24]= {7.416E+00,2.097E+01,3.058E+01,4.064E+01,5.035E+01,6.033E+01,7.010E+01,8.008E+01,9.021E+01,1.004E+02,
                           1.100E+02,1.200E+02,1.301E+02,1.400E+02,1.500E+02,1.601E+02,1.700E+02,1.801E+02,1.900E+02,2.000E+02,
                           2.100E+02,2.200E+02,2.301E+02,2.397E+02};

// vcsTab contains cross section in GeV^{-2} units multiplied on v=2*Pcm/Mcdm. It allows to decrease number of points for interpolation.
// In case of Sommerfeld it is reasonable to tabulate cs*v^{1+eps} 
       double vcsTab[24]= {2.349e-09,2.343e-09,2.337e-09,2.327e-09,2.314e-09,2.299e-09,2.282e-09,2.263e-09,2.240e-09,2.216e-09,
                           2.191e-09,2.163e-09,2.133e-09,2.102e-09,2.070e-09,2.036e-09,2.002e-09,1.967e-09,1.932e-09,1.895e-09,
                           1.859e-09,1.822e-09,1.785e-09,1.749e-09};
       double v=2*PcmIn/Mcdm;    
       double cs;
       if(PcmIn>pcmTab[23])
       {  double PcmLast=pcmTab[23], vLast=2*PcmLast/Mcdm;   csLast=vcsTab[23]/vLast;
          cs=csLast*(Mcdm*Mcdm+PcmLast*PcmLast)/(Mcdm*Mcdm+PcmIn*PcmIn);
          printf("Warning from improveCrossSection: PcmIn=%.3E out of tabulated region. 1/s interpolation is applied \n",PcmIn);
       } 
       else  cs=polint3(PcmIn,24,pcmTab,vcsTab)/v;
       
// debugging print  
       printf(" improveCrossSection(~X,~X,Z,Z) PcmIn=%.3E   cs %.3E => %.3E\n", PcmIn, *res,cs);
       
       *res=cs;
       return;
    }     

// for ~X,~X->W-,W+

    if(n1==35 &&  n2 == 35 &&   n3 == -24 &&  n4 == 24)
    { 

       double pcmTab[24]={7.416E+00,2.097E+01,3.058E+01,4.064E+01,5.035E+01,6.033E+01,7.010E+01,8.008E+01,9.021E+01,1.001E+02,
                          1.110E+02,1.200E+02,1.301E+02,1.400E+02,1.500E+02,1.601E+02,1.700E+02,1.801E+02,1.900E+02,2.000E+02,
                          2.100E+02,2.200E+02,2.301E+02,2.397E+02};
                          
       double vcsTab[24]={3.062e-09,3.056e-09,3.047e-09,3.035e-09,3.020e-09,3.001e-09,2.980e-09,2.956e-09,2.928e-09,2.898e-09,
                          2.863e-09,2.832e-09,2.795e-09,2.757e-09,2.717e-09,2.674e-09,2.632e-09,2.588e-09,2.544e-09,2.499e-09,
                          2.453e-09,2.407e-09,2.361e-09,2.317e-09};

       double v=2*PcmIn/Mcdm;  
       double cs=polint3(PcmIn,24,pcmTab,vcsTab)/v;
       if(PcmIn>pcmTab[23]) printf("Warning from improveCrossSection: PcmIn=%.3E out of tabulated region\n",PcmIn); 
// debugging print  
       printf(" improveCrossSection(~X,~X,Z,Z) PcmIn=%.3E   cs %.3E => %.3E\n", PcmIn, *res,cs);
       
       *res=cs;
       return;
    }     
}
#endif 

