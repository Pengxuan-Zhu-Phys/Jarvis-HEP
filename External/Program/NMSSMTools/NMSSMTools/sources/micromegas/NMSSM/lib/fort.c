#include"pmodel.h"
#include"../../include/micromegas_aux.h"
#include"pmodel_f.h"

int nmhwarn_(int *file)
{ 
  char fname[100];
  FILE*f;
  int r;

  if(*file==0) return slhaWarnings(NULL);
                                                                                
  sprintf(fname,"%d.tmptxt",getpid());
  f=fopen(fname,"w");
  r=slhaWarnings(f);
  fclose(f);
  fortreread_(file,fname,strlen(fname));
  unlink(fname);
  return r;
}


void o1contents_(int *Nch)
{
  char fname[20];
  FILE*f;

  sprintf(fname,"%d.tmptxt",getpid());
  f=fopen(fname,"w");
  o1Contents(f);
  fclose(f);
  fortreread_(Nch,fname,strlen(fname));
  unlink(fname);
}
int nmssmewsb_(void){ return nmssmEWSB();}

int  nmssmsugra_(double*m0, double*mhf, double*a0, double*tb,
double*sgn, double*Lambda, double*aLambda, double*aKappa, 
double*xif, double*xis, double*muP, double*MSPQ, double*M3HQ)
{
  return  nmssmSUGRA(*m0, *mhf, *a0, *tb, *sgn, *Lambda,*aLambda, *aKappa,
  *xif,*xis,*muP,*MSPQ,*M3HQ );
} 

int readslha_(char * fname, int *mode, int len)
{ char c_name[200];
  fName2c(fname, c_name,len);
  return readSLHA(c_name,*mode);
}

int  readvarnmssm_(char * f_name,int len)
{
  char c_name[100];
  fName2c(f_name,c_name,len);
  
  return readVarNMSSM(c_name);
}
