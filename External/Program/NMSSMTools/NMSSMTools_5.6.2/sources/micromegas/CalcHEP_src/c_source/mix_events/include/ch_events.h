#ifndef __CH_EVENT__
#define __CH_EVENT__

#include "model.h"

extern void ch2pythColors(int*CC,int nin, int nout, int * ch_color, int * pth_color);


typedef struct eventfile_info
{ 
  struct eventfile_info* next;
  FILE *F;  
  char * fileName;
  char pdf[2][100];
  long ltime;
  int Nin,Nout;
  double cs;
  double inMom[2];
  int inPID[2];
  int PIDs[10];
  char pName[10][P_NAME_SIZE];  
  double pmass[10];
  long nEvents;
  double wCoeff;
  long cEvent;
  long FirstEventPos;
  long CurrentEventPos;
  int  firstRd;  
} eventfile_info;


typedef struct decay_info
{ 
  struct decay_info* next;
  int ID;
  double totWidth;
  double  width;
  eventfile_info * List;
}  decay_info;

extern int nFiles;
extern eventfile_info * initEventFile(char* fname);
extern int readEvent(eventfile_info *Finfo, int *Nmom, double * mom, int * clr, double *Qf, double *alphaQCD, double * w);

extern eventfile_info * All;

extern decay_info * Decays;


#endif