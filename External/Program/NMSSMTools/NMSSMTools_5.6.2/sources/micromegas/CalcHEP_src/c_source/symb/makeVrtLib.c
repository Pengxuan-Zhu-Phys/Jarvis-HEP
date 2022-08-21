#include<stdio.h>
#include <string.h>
#include <math.h>
#include "model.h"
#include "read_mdl.h"
#include "procvar.h"
#include "reader_c.h"
#include "parser.h"
#include "symbolic.h"
#include "nType.h"

int menulevel; // for compilation 

int main(int argv, char**argc)
{
  int i,i10,nv,nLn,L,M,N;
  FILE*f,*fExt;  
  int nVar=0,nFunc=0,first;
  int mode,err; 
  char path[200];
  char * CalcHEP=NULL;
  
//  printf("argv=%d\n",  argv);
  if(argv<7 || argv>8) { printf("makeVrtLib: 6 or 7 arguments expected:  1) path To models 2) model number  3) label for library\n");
                         printf("4-7[,7] - particle names.\n"); return 1;}

  L=strlen(argc[0]);
  CalcHEP=malloc(L+10);
  strcpy(CalcHEP,argc[0]);
  CalcHEP[L-16]=0;  

  if(sscanf(argc[2],"%d",&M)!=1) { printf("Second argument should be a number\n"); return 1;}
  N=argv-4;
  char *field[4];
  for(i=0;i<N;i++) field[i]=argc[i+4];

//printf("Model=%d Nf=%d [",M,N); for(i=0;i<N;i++) printf(" %s ",field[i]); printf("]\n");

   
//   if(argv>=4)sscanf(argc[3],"%d",&mode); else mode=0;
   err= writeVertexCode(argc[1],M,N,field,argc[3]);  
   return err;
}

