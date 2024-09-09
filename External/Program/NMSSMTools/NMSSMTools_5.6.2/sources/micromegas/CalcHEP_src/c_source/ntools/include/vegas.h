#ifndef __VEGAS__
#define __VEGAS__

#include <pthread.h>
#include "n_proc.h"

#define MAX_DIM   15
#define MAX_NDMX  50


typedef  struct vegasGrid 
{
    int  dim;                           // dimension 
    double(*fxn)(double *,double);  // integrand
    int  ndmx;   
    double  x_grid[MAX_DIM][MAX_NDMX+1];
    int intCubes;
    int evnCubes;
    int NgI[MAX_DIM];
    int NgE[MAX_DIM];
    double  * fMax;
    pthread_mutex_t key;
} vegasGrid;

        
extern vegasGrid *  vegas_init
(  int dim,                        // dimension 
   double(*fxn)(double *,double),  // integrand
   int ndmx                        // size of grid 
);

extern void setEventCubes(vegasGrid*vegPtr, int nCubes);

extern void vegas_finish( vegasGrid * vegPtr);

extern int (*vegas_control)(double x); 

extern long vegas_int(vegasGrid * vegPtr, 
   long ncall0,                       /* number of integrand calls */
   double alph,                       /* rate of grid improvement  */
   int nCore,                         /* number of cores */ 
   double *ti,                        /* integral estimation */ 
   double *tsi                        /* standard deviation */
);

typedef struct 
{ double  eff;  /* efficiency */
  int     nexc;    /* number of points where max was improved */
  double  rmax;    /*  rate new/old max for such points */
  int     lmax;    /* number of subsequent  events from one cube */
  int     neg;     /* number of negative events */
  int     nan;      /* number of points with NaN */
  double  sumW;    /* sum of weights */
  double  sumW2;   /* sum of weights^2 */ 
} event_stat; 

extern long vegas_events(
vegasGrid * vegPtr, 
long  nEvents,
double gmax,  
void (*out)(double* ,double),
int recalc,    /* recalculate events in cube in case of new maximum */
int nCore,     /* number of cores */ 
event_stat * stat
);

#endif
