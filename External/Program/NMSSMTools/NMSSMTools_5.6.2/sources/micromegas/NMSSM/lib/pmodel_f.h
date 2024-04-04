#ifndef __NMSSN_F_ 
#define __NMSSN_F_
 
extern int nmhwarn_(int *file);                    

extern void o1contents_(int *file);
/*   
       subroutine o1Contents(file)
       integer file
*/ 

extern int nmssmewsb_(void);
extern int  nmssmsugra_(double *m0, double* mhf, double* a0, double* tb,double*sgn,
                        double*Lambda, double *aLambda, double*aKappa, 
                        double*xif,double *xis,double*muP,double*msPQ,double*M3HQ);
extern int readslha_(char * fname,int *mode, int len);
extern int readvarnmssm_(char *fname, int len);

#endif
