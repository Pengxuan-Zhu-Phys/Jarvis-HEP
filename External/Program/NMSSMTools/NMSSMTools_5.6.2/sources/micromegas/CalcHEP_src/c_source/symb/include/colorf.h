#ifndef __COLORS__
#define __COLORS__

#include"diagrams.h"

#define MAXE     3                 /*  Standard QCD                    */


typedef enum { zv=1,      // ZV     Colourless vertex        
               t88=4,     // 2*2    Transfer gluon vertex   { 8, 8, 0}
               vFabc=8,   // 2*2*2  Three gluon vertex      { 8, 8, 8}
               t33=15,    // 3*5    q-q_  tranfer vertex    { 0,-3, 3}
               v333=27,   // 3*3*3  q-q-q                   { 3, 3, 3}
               v833=30,   // 2*3*5  g-q-q_                  { 8,-3, 3}
               t66=77,    // 7*11   6-6_ transfer           { 0,-6, 6}
               V633=99,   // 11*3*3 6_-q-q                  {-6, 3, 3}
               V333=125,  // 5*5*5  q_-q_-q_                {-3,-3,-3}
               v633=175,  // 7*5*5  6-q_-q_                 { 6,-3,-3}
               v866=154   // 2*7*11 g-6-6_                  { 8,-6, 6}
              } vtype;

typedef struct cvertex
   {
      vtype        vt;           /* vertex type: zv, tv, g2, qg, vFabc */
      int        e[MAXE];        /* array of edges (linked vrtex numbers) */
   }  cvertex;


typedef struct factor
   { 
      long nc[4*MAXINOUT];
      int  len;
      int  pow2;
      long powN;
   }  factor;


extern factor * colorFactor (int nv, cvertex * vl);
extern void fct_num_calc(factor * fct,int Nc, long * n, long *d);
extern void fct_print(factor * fct, char *s);

extern int t2k2(vcsect* g, int * nv, cvertex * vl);
extern vtype typev(vert0 v,int valence);
#endif  
