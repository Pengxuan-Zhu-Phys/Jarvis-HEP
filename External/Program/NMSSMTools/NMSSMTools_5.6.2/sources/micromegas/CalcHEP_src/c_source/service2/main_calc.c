/*
Copyright (C) 2002, Alexander Pukhov
*/

#include <ctype.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "parser.h"
#include "calc.h"


int main(int argn, char** argc)
{ char s[500], format[20];
  double r;
  int i;

  if(argn<2) return 1;
  if(argc[1][0]=='%')
  { int d;
    sscanf(argc[1]+1,"%d",&d);
    sprintf(format,"%%.%dE\n",d);
    i=2;
  } else 
  { strcpy(format,"%.3E\n");
    i=1;
  }         
  for(s[0]=0; i < argn;i++) strcat(s,argc[i]);
  userFuncTxt=s;
  r = userFuncNumC();
  if(!rderrcode)  printf(format,r);
  return rderrcode;
}          

