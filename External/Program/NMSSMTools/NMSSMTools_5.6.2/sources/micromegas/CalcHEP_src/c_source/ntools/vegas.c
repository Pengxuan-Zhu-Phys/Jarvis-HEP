#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>

#include "vegas.h"
#include "drandXX.h"

const pthread_mutex_t keyInit=PTHREAD_MUTEX_INITIALIZER;

static void drand_arr(int dim, double * x)
{ int i;
  unsigned int umax= UINT_MAX;
  int rest=dim;
  unsigned int randpos=umax*drandXX();
  
  for(i=0;i<dim;i++) x[i]=-1; 
  
  for(i=0;i<dim;i++)
  {
     int pos=randpos%rest;
     int j;
     for(j=0;;j++) if(x[j]<0)
     {if(pos)pos--;else {x[j]=drandXX();break;}}
     randpos/=rest;
     umax/=rest;
     rest--;
     if(umax<rest){ umax= UINT_MAX;randpos=umax*drandXX();}
  }
}

static void  generateVegasCubes(int dim, int*nCubes, int*size) 
{  int i;
   long nC1=1,nCmem=*nCubes;
   double nCd=nCmem+0.001;
   for(i=0;i<dim;i++)
   {  int sz   =pow(nCd,1./(dim -i));
      nCd /= sz;
      nC1 *= sz;
      if(size) size[i]=sz;
   }
   if(nC1!=nCmem)  *nCubes=nC1;
}

void setEventCubes(vegasGrid*vegPtr, int nCubes)
{ 
  generateVegasCubes(vegPtr->dim,&nCubes, vegPtr->NgE);
  vegPtr->evnCubes=nCubes;
  free(vegPtr->fMax); 
  vegPtr->fMax=NULL;
}


static void Local2Global(vegasGrid *vegPtr, int*Kg,int*Ng, double*XLOC,double*XGLOB,double*JACOB,int*GRID_LOC)
{ int j,n;
  double xlj,xn,xn_;
  double dim=vegPtr->dim;
  int nd=vegPtr->ndmx, nd1=nd+1;
                                 
  *JACOB=1;     
  for (j = 0; j < dim; ++j)                       
  {  xlj = (Kg[j ] + XLOC[j])/Ng[j];
     n=(int)(xlj*nd);
     if (n) xn_=  vegPtr->x_grid[j][n]   ;else  xn_=0;               
            xn =  vegPtr->x_grid[j][n+1];
     XGLOB[j] = xn_ +(xn-xn_)*(xlj*nd-n);          
     *JACOB *= (xn-xn_)*nd;   
     if(GRID_LOC) GRID_LOC[j] = n;                 
  }
}


vegasGrid *  vegas_init(int dim, double(*fun)(double *,double), int nd)
{ 
   vegasGrid * vegPtr;
   
   if((dim>MAX_DIM)||(nd>MAX_NDMX)) return NULL;
   vegPtr=(vegasGrid * )malloc(sizeof(vegasGrid));
   if(vegPtr)
   {  int i,j;
      vegPtr->ndmx=nd;
      vegPtr->dim = dim;
      vegPtr->fxn=fun;
      for(j=0;j<dim;j++) for(i=0;i<=nd;i++)  vegPtr->x_grid[j][i]=i/(double)nd;         
      vegPtr->evnCubes=0;
      vegPtr->intCubes=0;
      for(i=0;i<dim;i++) { vegPtr->NgI[i]=vegPtr->NgE[i]=0;}
      vegPtr->fMax=NULL;
      vegPtr->key=keyInit;
   }
   return vegPtr;
}


void vegas_finish(vegasGrid * vegPtr) 
{   if(vegPtr)
    { 
       if(vegPtr->fMax) free(vegPtr->fMax);
       free(vegPtr); 
    }
}

/*     			*  VEGAS  *
      SUBROUTINE PERFORMS NDIM-DIMENSIONAL MONTE CARLO INTEG'N 
      - BY G.P. LEPAGE    SEPT 1976/(REV)AUG 1979 
      - ALGORITHM DESCRIBED IN J COMP PHYS 27,192(1978) 
*/

typedef struct { vegasGrid *vegPtr;
                 int npg; 
                 int * Kg;
                 int end;
                 int nCore;
                 pthread_mutex_t key;
               } cycle_str;   

typedef struct { double  ti;   
                 double  t2i;  
                 double dd[MAX_DIM][MAX_NDMX];
               } int_cycle_out; 

static void* vegas_cycle(void * par_)
{ 
  int i,j,l,k;
  double f,f2,fb,f2b;
  double x[MAX_DIM], xlocal[MAX_DIM];
  int    ia[MAX_DIM],kg[MAX_DIM],Kg_[MAX_DIM];
  int * ng,*Ng;

  cycle_str * par= (cycle_str *)par_; 
  
  int dim=par->vegPtr->dim;
  int Ndmx=par->vegPtr->ndmx;
  int npg=par->npg;
  int nCore=par->nCore;
  double*fMax=par->vegPtr->fMax; 
  int_cycle_out*result=malloc(sizeof(int_cycle_out));
  
  Ng=par->vegPtr->NgI;
  if(fMax) ng=par->vegPtr->NgE;

  for(i=0; i<dim; i++) for(j=0;j<Ndmx;j++)  result->dd[i][j]=0;  
  result->ti=0;
  result->t2i=0;
  
   for( ;  ;) 
   {
     if(nCore)pthread_mutex_lock(&(par->key)); 
       if(par->end) { if(nCore)pthread_mutex_unlock(&(par->key)); return result ;}
       for(i=0;i<dim;i++) Kg_[i]=par->Kg[i];
       for(i=dim-1; i>=0; i--){ if(++(par->Kg[i])< Ng[i]) break; else par->Kg[i]=0;}
       if(i<0) par->end=1;       
     if(nCore)pthread_mutex_unlock(&(par->key));                          
     
      fb  = 0;
      f2b = 0;
      for(i = 0; i<npg; i++) 
      {   double w;
          if(nCore)pthread_mutex_lock(&drandXX_key);  
            for(k=0;k<dim;k++)  xlocal[k]= drandXX();    
          if(nCore)pthread_mutex_unlock(&drandXX_key);            
          Local2Global(par->vegPtr ,Kg_, Ng, xlocal,x, &w, ia);
          f = par->vegPtr->fxn(x,w)*w;

          if(!isfinite(f)) 
          { if(nCore)pthread_mutex_lock(&(par->key));
            par->end=-1;
            if(nCore)pthread_mutex_unlock(&(par->key));
            return result;
          } 

          fb += f;
          f2= f*f;
          f2b += f2;
          for(j=0;j<dim;j++) result->dd[j][ia[j]] += f2;
          if(fMax)
          { 
            for(j=0;j<dim;j++) { double xj= (Kg_[j]+xlocal[j])/Ng[j];    kg[j]=xj*ng[j];}
            for(j=1, l=kg[0];j<dim;j++) l=l*ng[j]+kg[j];
            f=fabs(f);
            if(nCore)pthread_mutex_lock(&(par->vegPtr->key));
               if( fMax[l]<f ) fMax[l]=f;
            if(nCore)pthread_mutex_unlock(&(par->vegPtr->key)); 
          }  
      }
      f2b = sqrt(f2b/npg);
      fb /=npg; 
      f2b = (f2b - fb) * (f2b + fb)/(npg-1); 
      result->ti  += fb;
      result->t2i += f2b;   
    }
}

int (*vegas_control)(double x)=NULL;

static void* int_control(void * par_)
{  
  cycle_str * par=(cycle_str *)par_;
  int dim=par->vegPtr->dim;
  int nCore=par->nCore;
  if(!vegas_control) return NULL;
  vegas_control(0.);
  for(;;)
  {
   // usleep(2000);
    usleep(2000);
   if(nCore)pthread_mutex_lock(&(par->key));
    if(par->end) { if(nCore)pthread_mutex_unlock(&(par->key));return NULL;}
    {
      int i;
       double x=par->Kg[0];
        for(i=1;i<dim;i++) x=par->Kg[i]+par->vegPtr->NgI[i]*x;  
      x/=par->vegPtr->intCubes;
      if(vegas_control(x)) { par->end=2;if(nCore)pthread_mutex_unlock(&(par->key));return NULL;}
    }
   if(nCore)pthread_mutex_unlock(&(par->key));
  }   
}


long vegas_int(vegasGrid * vegPtr, long ncall0, double alph, int nCore,  double *ti, double *tsi)
{
   int dim= vegPtr->dim;
   int Ndmx=vegPtr->ndmx;
   double dd[MAX_DIM][MAX_NDMX];

   double x[MAX_DIM];
   int    Kg[MAX_DIM];
   int i,j,k;
   double  f2,fb,f2b;
   long l;
   cycle_str par;
   
   if(alph==0)
   {
     if(vegPtr->evnCubes && !vegPtr->fMax)
     { long cC;
       vegPtr->fMax=malloc(vegPtr->evnCubes*sizeof(double));
       for(cC=0;cC<vegPtr->evnCubes;cC++) vegPtr->fMax[cC]=0;
       
     }
   } else  if(vegPtr->fMax){free(vegPtr->fMax);vegPtr->fMax=NULL;}
   
   vegPtr->intCubes=ncall0/2;
   generateVegasCubes(vegPtr->dim,&(vegPtr->intCubes),vegPtr->NgI);
   par.nCore=nCore; 
   par.vegPtr=vegPtr;
   par.npg     = ncall0/vegPtr->intCubes;
   par.Kg      = Kg;
   par.end     = 0;
   par.key     =keyInit;
 

   for(i=0; i<dim;i++){ Kg[i]=0;}
/*    - MAIN INTEGRATION LOOP */
   {  
      int_cycle_out*result; 
      if(nCore)
      { pthread_t *threads = malloc(nCore*sizeof(pthread_t));
        pthread_t control;
        int_cycle_out*result; 
      
        for (k=0;k<nCore;k++) pthread_create(threads+k,NULL,vegas_cycle, &par);

        if(vegas_control) int_control(&par);

        for(i=0; i<dim;i++) for(j=0;j<Ndmx; j++) dd[i][j]=0;
       *ti=0; *tsi=0;
     
        for (k=0;k<nCore;k++) 
        {  pthread_join(threads[k],(void**) &result);
           for(i=0; i<dim;i++) for(j=0;j<Ndmx; j++) dd[i][j]+=result->dd[i][j];
           *ti+=result->ti;
           *tsi+=result->t2i;
           free(result);
        }          
        free(threads);
      } else
      {  
         result=vegas_cycle(&par);
         for(i=0; i<dim;i++) for(j=0;j<Ndmx; j++) dd[i][j]=result->dd[i][j];
         *ti=result->ti;
         *tsi=result->t2i;
         free(result);
      }        
   }
   if(par.end <0) return -1; else if(par.end >1) return 0; 
   { long nCubes=vegPtr->intCubes; 
     *ti/=nCubes;
     *tsi=sqrt( fabs(*tsi/(nCubes*nCubes)));
   }
    if(alph>0)  for(j=0; j<dim; j++)     /* REFINE GRID */
    {  double r[MAX_NDMX], rc = 0, dt=dd[j][0];
           
       for(i=1; i<Ndmx; i++) dt += dd[j][i];
       for(i=0; i<Ndmx; i++)
       {
           r[i] = 0;
           if( dd[j][i]  > 0)
           {  double xoln = log(dt/dd[j][i]);
              if(xoln<0.01)         r[i]=1;
              else if (xoln <= 70)  r[i] = pow( (1 - exp(-xoln))/xoln, alph);
              else                  r[i] = pow(  1/xoln,               alph);
           }
           rc += r[i];  
       }
       rc /= Ndmx;
       if(rc)
       {  double x[MAX_NDMX],dr=0,xo=0,xn=0;
          int k=0;
          for(i=1;i<Ndmx;) 
          {  for(;dr<rc;k++) dr+=r[k];
             xo=vegPtr->x_grid[j][k-1];
             xn=vegPtr->x_grid[j][k];
             for(;dr>=rc; i++) {dr-=rc; x[i] = xn-(xn-xo)*dr/r[k-1];}
          }
    	  for(i=1;i<Ndmx;i++) vegPtr->x_grid[j][i]=x[i];    
          vegPtr->x_grid[j][Ndmx] = 1;
       }
   }
   return  par.npg*vegPtr->intCubes;
} 


typedef struct { vegasGrid *vegPtr;
                 long  nEvents;
                 double gmax;
                 int recalc;
                 void(*out)( double *,double); 
                 long  Ntry;
                 long  cEvent;
                 int    end;
                 int   nCore;
                 pthread_mutex_t key;
               } event_str;   

static void* event_cycle(void * par_)
{
   event_str*par=(event_str*)par_;
   int nCore=par->nCore;
   int dim=par->vegPtr->dim;
   int Ndmx=par->vegPtr->ndmx;
   int nCubes=par->vegPtr->evnCubes;
   double (*fxn)(double*,double)=par->vegPtr->fxn;
   double gmax=par->gmax;
   double*smax=par->vegPtr->fMax;
      
   int  Ng[MAX_DIM];
   double xlocal[MAX_DIM],x[MAX_DIM];
   int i,k;
   long Ntry;
   struct  { double x[MAX_DIM]; double f; double rnd; int sgn; }  *events=NULL;
   int nRecT=0,nmax=0;
   event_stat * stat=malloc(sizeof(event_stat)); 
   
   stat->neg=0;
   stat->nexc=0;
   stat->lmax=1;
   stat->rmax=1;
   stat->nan=0;
   stat->sumW=0;
   stat->sumW2=0;
          
   generateVegasCubes(dim, &nCubes,Ng);
  
   for(;; ) 
   {  long L,L0,L1;
      double f,ds,w,sum,weight,rc1,rc2,oldMax,newMax;

      int n,sgn;   
      int Kg[MAX_DIM];
      
      if(nCore)pthread_mutex_lock(&(par->key));
        if(par->end) { if(nCore)pthread_mutex_unlock(&(par->key)); free(events); return stat;}  
        Ntry= ++(par->Ntry); 
      if(nCore)pthread_mutex_unlock(&(par->key));
      if(nCore)pthread_mutex_lock(&drandXX_key);     
        rc1=drandXX();
        rc2=drandXX();
        drand_arr(dim,xlocal);           
      if(nCore)pthread_mutex_unlock(&drandXX_key);
        
      L0=0;
      L1=nCubes-1;
      
      if(nCore)pthread_mutex_lock(&(par->vegPtr->key));
        sum=smax[nCubes-1];
        rc1*=sum;
        if(smax[0]>rc1) L1=0; else
        {
           while(L0+1 < L1)     
           {  L=(L0+L1)/2;
              if(smax[L]>rc1) L1=L;  else L0=L;  
           }
        } 
// L0<=L1          
//        if(smax[L0]>rc2) L=L0; else L=L1;
        oldMax= L1? smax[L1]-smax[L1-1] : smax[0];
      if(nCore)pthread_mutex_unlock(&par->vegPtr->key);
 
//      if(informline(cEvent,nEvents))  break;

      L=L1; for(i=dim-1;i>=0; i--) {Kg[i]=L%Ng[i]; L=L/Ng[i];}
      Local2Global(par->vegPtr,Kg,Ng, xlocal,x,&w,NULL);
      f=(par->vegPtr->fxn)(x,w)*w; if(f<0){f=-f;sgn=-1;} else sgn=1;      
      if(!isfinite(f)) {stat->nan++; continue;}
      if(f<=oldMax*gmax*rc2) continue;
      if(nCore)pthread_mutex_lock(&par->key);
       if(par->end){ if(nCore)pthread_mutex_unlock(&par->key);  continue;}
       par->cEvent++; 
       if(sgn<0) stat->neg++;
       weight=f/(oldMax*gmax);
       if(weight<1 || par->recalc) weight=1;
       weight*=sgn;
       stat->sumW+= weight;
       stat->sumW2+= weight*weight;
       (par->out)(x,sgn*weight);      
       if(par->cEvent>=par->nEvents) {par->end=1; if(nCore)pthread_mutex_unlock(&par->key);  continue;}
      if(nCore)pthread_mutex_unlock(&par->key);
      
      newMax=f;
      if(newMax>oldMax*gmax && par->recalc) 
      {
         long dN=0;
         int i,k,nRec=0;   
         for(dN=1;dN < Ntry*newMax/sum/gmax;dN++)
         { double rc;
             if(nCore)pthread_mutex_lock(&drandXX_key);                                            
              for(k=0;k<dim;k++) xlocal[k]=drandXX();
              rc=drandXX();
             if(nCore)pthread_mutex_unlock(&drandXX_key); 
              Local2Global(par->vegPtr,Kg,Ng, xlocal,x,&w,NULL);
              f=(*fxn)(x,w)*w; if(f<0){f=-f;sgn=-1;} else sgn=1;
              
              if(f>oldMax*gmax)  
              {  if(nRec >= nRecT) events=realloc(events, sizeof(*events)*(++nRecT)); 
                 events[nRec].f=f;
                 events[nRec].sgn=sgn;
                 events[nRec].rnd=rc;
                 for(i=0;i<dim;i++) events[nRec].x[i]=x[i];
                 nRec++;
                 if(f>newMax) newMax=f;
              }                 
         }
         if(nCore)pthread_mutex_lock(&par->vegPtr->key);
           oldMax= L1? smax[L1]-smax[L1-1] : smax[0];
         if(nCore)pthread_mutex_unlock(&par->vegPtr->key); 

         if(nCore)pthread_mutex_lock(&par->key); 
           for(n=0,k=0;k<nRec ;k++) if(events[k].f>oldMax &&
              events[k].f>newMax*events[k].rnd && par->cEvent<par->nEvents)
           {  par->cEvent++;
                 (par->out)(events[k].x,events[k].sgn);  
                 stat->sumW+=sgn;
                 stat->sumW2+=1;
                 n++;
           }                
           if(par->cEvent>=par->nEvents) par->end=1;            
            (par->Ntry)+=dN;
         if(nCore)pthread_mutex_unlock(&par->key);
         if(stat->lmax<n+1) stat->lmax=n+1; 
      }
      
      if(newMax>oldMax*gmax)
      {  ds=newMax-oldMax;
         if(nCore)pthread_mutex_lock(&(par->vegPtr->key));
           for(L=L1;L<nCubes;L++) smax[L]+=ds;
         if(nCore)pthread_mutex_unlock(&(par->vegPtr->key));
         stat->nexc++;
         if(stat->rmax < newMax/oldMax) stat->rmax=newMax/oldMax;
      }                                     
   }    
}

static void* event_control(void * par_)
{  
  event_str * par=(event_str *)par_;
  int nCore=par->nCore;
  if(!vegas_control) return NULL;
  vegas_control(0.);
  for(;;)
  {
   sleep(1);
   if(nCore)pthread_mutex_lock(&(par->key));
    if(par->end) { if(nCore)pthread_mutex_unlock(&(par->key));return NULL;}
    { double x=par->cEvent; x/=par->nEvents;
      if(vegas_control(x)) { par->end=2;if(nCore)pthread_mutex_unlock(&(par->key));return NULL;}
    }
   if(nCore)pthread_mutex_unlock(&(par->key));
  }   
}


long vegas_events(vegasGrid * vegPtr,  long  nEvents, double gmax, 
   void (*out)(double*,double),int recalc,  int nCore, event_stat * stat)
{
   int i, dim= vegPtr->dim;
   long cCube, nCubes=vegPtr->evnCubes;
   event_str par;
   event_stat*stat1;

   double * smax=vegPtr->fMax;
   
   if(!smax) return -1;
   for(cCube=1; cCube<nCubes; cCube++) smax[cCube]+=smax[cCube-1];
   par.nCore=nCore;   
   par.vegPtr=vegPtr;
   par.nEvents=nEvents;
   par.recalc=recalc;
   par.gmax=gmax;
   par.out=out;
   par.Ntry=0;    
   par.cEvent=0;
   par.end=0;
   par.key=keyInit;
  
   if(stat){ stat->neg=0; stat->lmax=0; stat->nexc=0; stat->rmax=1; stat->nan=0; stat->sumW=0; stat->sumW2=0;} 
   
/*    - MAIN INTEGRATION LOOP */
   if(nCore)
   { pthread_t *threads = malloc(nCore*sizeof(pthread_t));
     pthread_t control; 
     for (i=0;i<nCore;i++) pthread_create(threads+i,NULL,event_cycle, &par);
     event_control(&par);
     for (i=0;i<nCore;i++) 
     {  pthread_join(threads[i],(void**)&stat1);
        if(stat)
        { stat->neg+=stat1->neg;
          stat->nan+=stat1->nan;
          stat->sumW+=stat1->sumW;
          stat->sumW2+=stat1->sumW2;
          if(stat->lmax<stat1->lmax) stat->lmax=stat1->lmax;
          if(stat->nexc<stat1->nexc) stat->nexc=stat1->nexc;
          if(stat->rmax<stat1->rmax) stat->rmax=stat1->rmax;
        }  
        free(stat1); 
     }         
     free(threads);  
   } else 
   {  stat1=event_cycle(&par);
      if(stat)
      { stat->neg=stat1->neg;
        stat->nan=stat1->nan;
        stat->lmax=stat1->lmax;
        stat->nexc=stat1->nexc;
        stat->rmax=stat1->rmax;
      }
      free(stat1);  
   }                                                   
   
   { double mem=smax[0],mem_;
     for(cCube=1;cCube<nCubes;cCube++) { mem_=smax[cCube]; smax[cCube]-=mem; mem=mem_;} 
   } 
   if(stat) stat->eff=par.cEvent/(double)(par.Ntry);
   return par.cEvent;
}
