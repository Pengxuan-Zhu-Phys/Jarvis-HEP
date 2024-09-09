
#include <stdio.h>
extern float upperlim_(float *CL,int *If,int *N,float*FC,float *muB,float*FB,int*Iflag);
float UpperLim(float CL,int If, int N, float* FC, float muB,float*FB,int *Iflag) { return upperlim_(&CL,&If,&N,FC,&muB,FB,Iflag); }  