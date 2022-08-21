#ifndef __RANDOM__
#define __RANDOM__

#include <pthread.h>
  extern  double drandXX(void);
  extern  char * seedXX(char * init);
  extern  pthread_mutex_t drandXX_key;  
  
//  typedef struct {unsigned long Xlong,Xshort;} rndState;
  
//  extern double drandXX_(rndState * rnd );
  
#endif
