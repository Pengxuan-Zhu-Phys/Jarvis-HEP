#include "../include/micromegas_aux.h"


double uConversion(int u1,int u2)
{
   double inGeV[7]={ 1, 5.609588650E23,  8.6173303E-14, 1.973269788E-14,6.582119514E-25, 1.220910E19, 6.241509126E+02};
     
   if(u1<0||u2<0|| u1>6 || u2>6) return 0;
   
   if(u2==0) return inGeV[u1]; 
   if(u1==0) return 1/inGeV[u2];
   return inGeV[u1]/inGeV[u2];
}
 