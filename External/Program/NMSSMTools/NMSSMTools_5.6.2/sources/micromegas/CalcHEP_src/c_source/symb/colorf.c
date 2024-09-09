/*
 Written by  Alexander Kryukov 1990-1999 
 Rewritten by Alexander Pukhov 2014 
    a) color 6 and  vertexes 333, 633, and 866  are added    
    b) general conception was saved, but implementation 
        was changed completely.  
*/


#ifdef CALCHEP
#include "syst2.h"
#define malloc m_alloc
#else
#include<stdio.h>
#include<stdlib.h>
#include <memory.h>
#include <string.h>

#endif
#include "syst.h"
#include "colorf.h"


int MAXGLEN, MAX_POW;
#define CERRLEV1 0    /*  Run time error level            */
#define CERRLEV2 0    /*  Halt level                      */ 

static int cerror(int n,char* s)
{ /*  - generate error message occur in color package - 08/01/90  */

   fprintf(stdout,"***** %s\n",s);
        if (n > CERRLEV1) fprintf(stderr,"Error  %u in  colorf.c \n",n), sortie(55);
   else if (n > CERRLEV2) sortie(56); 
   return 0;
}  

/*==============     FACTOR SESSION ==================*/

static void fct_init(factor * fct) 
{ int i;
  fct->len=0; fct->powN=0; fct->pow2=0; 
  for(i=0;i<4*MAXINOUT;i++) fct->nc[i]=0;
}

static void add_fct(factor *f, int sgn, int pow2, int powNm1,int powN,int powNp1)
{ int dd,i,k,Len,s[4*MAXINOUT];

  if(!sgn) return;

  if(pow2 > f->pow2)  sgn*=1<<(pow2 - f->pow2) ;   // pow2 matching
  else if(f->pow2 > pow2) 
  {  long q=1<<(f->pow2 - pow2);
     for(i=0;i<f->len;i++) f->nc[i]*=q;
     f->pow2=pow2;
  }       

  if(powN > f->powN ) powN-= f->powN;              // powN matching 
  else  
  {
    dd=f->powN-powN;
    f->powN=powN;
    powN=0;

    for(i=f->len-1;i>=0;i--) f->nc[i+dd]=f->nc[i];
    for(i=0;i<dd;i++) f->nc[i]=0;
    f->len+=dd;
  } 

  Len=1+powN+powNm1+powNp1;
  if(Len>f->len)
  { for(i=f->len;i<Len;i++) f->nc[i]=0;
    f->len=Len;
  } else Len=f->len;

  for(i=0;i<Len;i++) s[i]=0;
      
  s[0]=sgn; if(powNm1&1) s[0]*=-1; 


  for(k=1;k<=powNm1;k++)for(i=k;i>=1;i--) s[i]-=s[i-1];
  for(k=1;k<=powNp1;k++)for(i=k+powNm1;i>=1; i--) s[i]+=s[i-1];

  for(i=0;i<=powNm1+powNp1;i++) f->nc[i+powN]+=s[i];

  for(i=0; i<f->len && f->nc[i]==0 ;i++);
  if(i)
  {  dd=i;
     for( ;i<f->len; i++) f->nc[i-dd]= f->nc[i];
     f->len-=dd;
     f->powN+=dd;
  }
  for(;f->len && f->nc[f->len-1]==0;) f->len--; 
  
  if(!f->len) {f->pow2=0; f->powN=0; return;}
  for(;;)
  { 
    for(i=0;i<f->len;i++) if(f->nc[i]%2) break;
    if(i<f->len) break;
    f->pow2++; for(i=0;i<f->len;i++) f->nc[i]/=2;
  }  
}
/* ================= END OF FACTOR ==================== */

#define CDEBLEV  0    /*  Total debug level               */ 

/*################################################################*/
       
typedef struct cgraph 
   {  struct cgraph *    next;

      int          sgn, pow2, powN, powNm1,powNp1; /* factor */
      int          en;             /* Name of next edge */
      int          gl;             /* Number of vertecies (graph length) */
      int         mgl;             /* max Number of vertecies */
      cvertex  *   vl;             /* array of verticies */
   }  cgraph; 

/* ************************** Cross reference ************************* */ 
/* *                                                                  * */ 
/* *  GevV                                                            * */ 
/* *    +-----> CError                                                * */ 
/* *                                                                  * */ 
/* *  GetEN                                                           * */ 
/* *                                                                  * */ 
/* *  CError                                                          * */ 
/* *                                                                  * */ 
/* *  WrCG                                                            * */ 
/* *                                                                  * */ 
/* ******************************************************************** */ 


//#  if (CDEBLEV > DEBLEV) 

static int ltype(int n, int l, cgraph* cg)
{ 
/*
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
*/
  switch(cg->vl[n].vt) 
  {
    case  t88:  switch(l)
                { case 0: return 8;
                  case 1: return 8;
                  case 2: return 10;
                }   
    case  vFabc:switch(l)
                { case 0: return 8;
                  case 1: return 8;
                  case 2: return 8;
                }   
  
    case  t33: switch(l)
                { case 0: return 11;
                  case 1: return -3;
                  case 2: return  3;
                }   
    
    case  v333: switch(l)
                { case 0: return 3;
                  case 1: return 3;
                  case 2: return 3;
                }   
   
    case  v833:switch(l)
                { case 0: return 8;
                  case 1: return -3;
                  case 2: return 3;
                }   
   
    case  t66:switch(l)
                { case 0: return 12;
                  case 1: return -6;
                  case 2: return  6;
                }   
    
    case  V633:switch(l)
                { case 0: return -6;
                  case 1: return  3;
                  case 2: return  3;
                }   
   
    case  V333:switch(l)
                { case 0: return -3;
                  case 1: return -3;
                  case 2: return -3;
                }   
   
    case  v633:switch(l)
                { case 0: return 6;
                  case 1: return -3;
                  case 2: return -3;
                }   
   
    case  v866:switch(l)
                { case 0: return 8;
                  case 1: return -6;
                  case 2: return 6;
                }   
   default : return 13;    
  }	

}

static int findv2(vtype vt,cgraph* cg,int * n)
{ /* return True and number first vertex with type VT in C-graph - 06/01/90 */ 
  int  i=0; 
  for(i=0;i<cg->mgl;i++) if(cg->vl[i].vt == vt) {*n=i; return 1;} 
  return 0;   
}    

static int findl2(int e1,int n,cgraph* cg,int*l)
{ /* - find number of vertex contane edge e1 - 08/01/90  */ 
                                             /* 19/03/99 */
  int  i = 1,j, e = cg->vl[n].e[e1]; 

  for(i=0;i< cg->mgl;i++) if(cg->vl[i].vt!=zv)
  for(j=0;j<3;j++) if(cg->vl[i].e[j]==e && (i!=n || j!=e1) )
  {  if(l)*l=j; return i;}

printf(" no couple for line %d  \n", e);         
   cerror(253,"FindL(2): nonconnected edge"); 
} 


static int testGraph(cgraph* cg)
{ int n,l,n_,l_;
  for(n=0;n<cg->mgl;n++) if(cg->vl[n].vt != zv) for(l=0;l<3;l++)
  { int t=ltype(n,l,cg),t_;
    if(t>8) continue;
//printf("n=%d v=%d\n", n,cg->vl[n].vt);    
    n_=findl2(l,n,cg,&l_);
    t_=ltype(n_,l_,cg);
    if( !((t==8 && t_==8) || t+t_==0)){ printf(" type(%d %d)= %d    type(%d %d)=%d\n",n,l,t,n_,l_,t_);      return n;}   
  }
  return 0;
}  

static void wrcg(cgraph* cg)
{/*  Write C-graph on standard device - 04/01/90 */ 
   int      i; 
   printf("Kr: {(%d*2^%d) N^%d (N-1)^%d  (N+1)^%d\n",
	   cg->sgn,cg->pow2,cg->powN,cg->powNm1,cg->powNp1);
   for (i =0; i< cg->mgl; i++)
      if (cg->vl[i].vt != zv)
      {
	printf("[%d/",i);
        switch(cg->vl[i].vt) 
        {
          case  t88:    printf("t88");   break;
          case  vFabc:  printf("vFabc"); break;
          case  t33:    printf("t33");   break;
          case  v333:   printf("v333");  break; 
          case  v833:   printf("v833");  break;
          case  t66:    printf("t66");   break;
          case  V633:   printf("V633");  break; // 11*3*3 6_-q-q   
          case  V333:   printf("V333");  break; 
          case  v633:   printf("v633");  break;
          case  v866:   printf("v866");  break;  
          default  :    printf("%d",cg->vl[i].vt); break;
        }	
	printf(" (%d,%d,%d)]",
	 cg->vl[i].e[0],cg->vl[i].e[1],cg->vl[i].e[2]);
      }   /* if */
   printf("}");
   
   i=testGraph(cg);
   if(i) printf(" error in vertex %d\n",i); else printf("ok\n");   
} 
//#endif


static int getv2(cgraph* cg,int lbl)
{ /* return number of first free vertex in C-graph - 04/01/90 */ 
 
   int      n,i,j; 
   for(n=0;n<cg->mgl;n++) if(cg->vl[n].vt == zv) break;
   if(n==cg->mgl)
   { cg->vl=realloc(cg->vl, sizeof(cvertex)*(cg->mgl+4));
     for(i=cg->mgl; i<cg->mgl+4;i++)
     { cg->vl[i].vt=zv;
       for(j=0;j<3;j++) cg->vl[i].e[j]=0;
     }
     cg->mgl+=4;
   } 
   cg->vl[n].vt = lbl; 
   (cg->gl)++;  
   return n;
}     

#define SINGL   1    /* Colour singlet */ 
#define TRIPL  -3    /* Colour triplet */ 
#define ATRIPL  3    /* Colour antitriplet */ 
#define OCTET   8    /* Colour octet   */ 
#define DEBLEV 10    /*  Debug level   */ 

static cgraph * initcg(int nv, cvertex*vl)
{ 
  int i,j;  
  cgraph *cg=malloc(sizeof(cgraph));
  cg->sgn=1;  cg->pow2=cg->powN=cg->powNm1=cg->powNp1=0;
  cg->en = 0; 
  cg->gl = nv;
  cg->mgl= nv+2;
  cg->next=NULL;
  cg->vl= malloc( sizeof(cvertex)*cg->mgl); 
   for (i = 0; i < nv; i++) 
   { cg->vl[i]=vl[i];
     for(j=0;j<3;j++) if(vl[i].e[j]>cg->en) cg->en=vl[i].e[j];  
   }
   for (i = nv; i < cg->mgl; i++)
   {
      cg->vl[i].vt = zv; 
      cg->vl[i].e[0]=cg->vl[i].e[1] = cg->vl[i].e[2] = 0; 
   } 
   return cg;
}





static cgraph * addcg(cgraph ** cg) 
{ /*- add C-graph CG to weight structure - 16/10/99 -*/   
   cgraph  *pgl = (cgraph *) malloc(sizeof(cgraph));
   memcpy(pgl, *cg ,sizeof(cgraph)); 
   pgl->next = *cg; 
   *cg = pgl; 
   pgl->vl=malloc(sizeof(cvertex)*pgl->mgl);
   memcpy(pgl->vl,pgl->next->vl,sizeof(cvertex)*pgl->mgl);    
   return pgl;
}

static void rem833_833(int n0,int n1, cgraph ** c)
{ /* - remove gluon connected vertex n0 and n1 (see fugure) from first C-graph - 08/01/90  */ 
  int  l01,l02,l11,l12, n0_, n1_, i0_, i1_;    
  cgraph * cg=*c;   
#  if (CDEBLEV > DEBLEV) 
     printf("rem833_833(%d,%d)",n0,n1);  wrcg(cg);
#  endif 
   l01=cg->vl[n0].e[1]; l02=cg->vl[n0].e[2];
   l11=cg->vl[n1].e[1]; l12=cg->vl[n1].e[2]; 

   cg->pow2--;
   cg->vl[n0].vt = cg->vl[n1].vt = zv; cg->gl -= 2;
         
   if(l01==l12 && l02==l11) { cg->powNm1++; cg->powNp1++;}  
   else  if(l01 == l12)
   { 
     cg->powNm1++; cg->powNp1++; cg->powN--;                 
     n1_ = findl2(1,n1,cg,&i1_); 
     cg->vl[n1_].e[i1_] = l02;   
   } 
   else  if(l02 == l11)
   {  
      cg->powNm1++;  cg->powNp1++; cg->powN--; 
      n1_ = findl2(2,n1,cg,&i1_); 
      cg->vl[n1_].e[i1_] = l01;
   
   } else 
   {  
      n0_ = findl2(1,n0,cg,&i0_);
      cg->vl[n0_].e[i0_] = l12;
      n1_ = findl2(1,n1,cg,&i1_);        
      cg->vl[n1_].e[i1_] = l02;
                                             
      cg=addcg(c); 
      cg->sgn*=-1;
      cg->powN--; 
      cg->vl[n0_].e[i0_] =  l02; 
      cg->vl[n1_].e[i1_] =  l12;  
   }     
      
#  if (CDEBLEV > DEBLEV) 
      wrcg(cg); 
      if (cg->next) wrcg(cg->next); 
#  endif 
} 




static void rem833_Fabc(int nq,int ng, cgraph  ** c)
/* - remove gluon connected vertex n0 and n1 (see figure)
     from first C-graph - 08/01/90  */ 
{ int  i0,i1,i2,        n1, n2,lg1,lg2,lq1,lq2,lg;
   cgraph * cg = *c;
                        /*         v1        */ 
                        /*  v2.....*.....v3  */ 
                        /*         :         */ 
                        /*         :         */ 
                        /*  v3-->--*-->--v4  */ 
                        /*          v0       */ 

#  if (CDEBLEV > DEBLEV) 
   printf("rem833_Fabc(%d,%d)\n",nq,ng);
#  endif 

   findl2(0,nq,cg,&i0);
   i1=(i0+1)%3;
   i2=(i1+1)%3;

   lq1=cg->vl[nq].e[1];  lq2=cg->vl[nq].e[2];
   lg1=cg->vl[ng].e[i1]; lg2=cg->vl[ng].e[i2];
   lg=cg->vl[nq].e[0];

   cg->vl[ng].vt=v833; n1=ng;
   n2=nq;

   cg->vl[n1].e[1]=lq1;  cg->vl[n1].e[2]=cg->vl[n2].e[1]=lg; 
   cg->vl[n1].e[0]=lg1;  cg->vl[n2].e[0]=lg2;
                                          
   cg=addcg(c);
   cg->sgn*=-1;
   cg->vl[n1].e[0]= lg2; cg->vl[n2].e[0]=lg1;

#  if (CDEBLEV > DEBLEV) 
       wrcg(cg); 
       wrcg(cg->next); 
#  endif 
}  /* RemQG_3G */ 


static void rem866_Fabc(int n866,int n3g, cgraph  ** c)
{  int   i0,i1,i2,l0,l1,l2,N1,N2,n1, n2; 
   cgraph * cg = *c;    
#  if (CDEBLEV > DEBLEV) 
   printf("rem866_Fabc(%d,%d)\n", n866,n3g);
#  endif 

   cg->pow2++;
   findl2(0,n866,cg,&i0);
   i1=(i0+1)%3;
   i2=(i1+1)%3;
   l0=cg->vl[n3g].e[i0]; l1=cg->vl[n3g].e[i1]; l2=cg->vl[n3g].e[i2];
              
   N1=getv2(cg,V633); cg->vl[N1].e[0]=cg->vl[n866].e[1];
   N2=getv2(cg,v633); cg->vl[N2].e[0]=cg->vl[n866].e[2];

   cg->vl[N1].e[2]=cg->vl[N2].e[2]= ++(cg->en);
   cg->vl[N1].e[1]= ++(cg->en); cg->vl[N2].e[1]=++(cg->en);
   
   n1=n866; cg->vl[n1].vt=v833; 
   n2=n3g;  cg->vl[n2].vt=v833;
   
   cg->vl[n1].e[1]=cg->vl[N1].e[1];  cg->vl[n2].e[2]=cg->vl[N2].e[1];
   cg->vl[n1].e[2]=cg->vl[n2].e[1]=l0;   
   cg->vl[n1].e[0]= l1; cg->vl[n2].e[0]=l2;

   cg=addcg(c);
   cg->sgn*=-1;
   cg->vl[n1].e[0]= l2; cg->vl[n2].e[0]=l1;

#  if (CDEBLEV > DEBLEV) 
       wrcg(cg);  wrcg(cg->next); 
#  endif 
}  

static void exp866(int n866, cgraph  ** c)
{ int   N1,N2,n833;

  cgraph * cg = *c;

#  if (CDEBLEV > DEBLEV) 
   printf("exp866(%d)\n",n866);
#  endif 

   cg->pow2++;

   N1=getv2(cg,V633);  cg->vl[N1].e[0]=cg->vl[n866].e[1];
   N2=getv2(cg,v633);  cg->vl[N2].e[0]=cg->vl[n866].e[2];
   
   cg->vl[N1].e[2]=cg->vl[N2].e[2]= ++(cg->en);
   
   n833=n866; cg->vl[n833].vt=v833;
   
   cg->vl[n833].e[1]=cg->vl[N1].e[1]= ++(cg->en); 
   cg->vl[n833].e[2]=cg->vl[N2].e[1]= ++(cg->en);
#  if (CDEBLEV > DEBLEV) 
       wrcg(cg);
#  endif 
}  

static int istadpole2(int n,cgraph* cg)
{ /* return True if vertex n is teadpole - 08/01/90  */ 
   if(cg->vl[n].e[0] == cg->vl[n].e[1] || 
      cg->vl[n].e[1] == cg->vl[n].e[2] || 
      cg->vl[n].e[0] == cg->vl[n].e[2]    )
   { cg->sgn=0;return 1;} else return 0;   
  
}

static void remtv(cgraph * pgl )
{ /* Remove transfered vertex from all C-graphs - 06/01/90     */ 

 int         n, n1; 
 int         vt0;   /* Original type */   
                                           /*  -->--*-->--v1  */ 
#if (CDEBLEV > DEBLEV)                     /*                 */ 
  printf(".......RemTV........\n");/*                 */
#endif                                     /*       v0        */ 
                                           /*  .....*.....v1  */ 
                                           
   while (findv2(t33,pgl,&n) || findv2(t88,pgl,&n)||findv2(t66,pgl,&n) )
   { int l;
#  if (CDEBLEV > DEBLEV)                                   
         if (pgl != NULL) wrcg(pgl);                      
#  endif                                                   
      vt0 = pgl->vl[n].vt;                                     
      pgl->gl--; 
      n1=findl2(1,n,pgl,&l); 
      if (n1==n)  switch(vt0)
      { case t33 : pgl->powN++; break;
        case t88 : pgl->powNm1++; pgl->powNp1++; break;
        case t66: pgl->powN++; pgl->powNp1++; pgl->pow2--; break;
      }                                                            
      else  pgl->vl[n1].e[l]= (vt0==t88)?pgl->vl[n].e[0]:pgl->vl[n].e[2];
      pgl->vl[n].vt = zv; 
   } 
   #  if (CDEBLEV > DEBLEV) 
         wrcg(pgl);
   #  endif                                                       
}   


static void expFabc(int n0, cgraph ** c)
{ 
 int  n1,n2,l0,l1,l2,m0,m1,m2; 
 cgraph * cg=*c;
  
#  if (CDEBLEV > DEBLEV) 
      printf("expFabc(%d)\n",n0);
#  endif 

   l0=cg->vl[n0].e[0];
   l1=cg->vl[n0].e[1];
   l2=cg->vl[n0].e[2];
  
   n1 = getv2(cg,v833);
   n2 = getv2(cg,v833);
   cg->vl[n0].vt=v833;
    
   cg->vl[n0].e[0]=l0;  cg->vl[n1].e[0]=l1; cg->vl[n2].e[0]=l2;
   
   cg->vl[n0].e[2]=cg->vl[n1].e[1]=m0=++(cg->en);
   cg->vl[n1].e[2]=cg->vl[n2].e[1]=m1=++(cg->en);
   cg->vl[n2].e[2]=cg->vl[n0].e[1]=m2=++(cg->en);

   cg->pow2++;
   cg=addcg(c); 
   cg->vl[n0].e[1]=cg->vl[n1].e[2]=m0;
   cg->vl[n1].e[1]=cg->vl[n2].e[2]=m1;
   cg->vl[n2].e[1]=cg->vl[n0].e[2]=m2;
   cg->sgn*=-1;      
   
#  if (CDEBLEV > DEBLEV) 
      wrcg(cg); wrcg(cg->next); 
#  endif 
}  /* Exp3G */ 



static int exp6(int n6,  cgraph ** c)
{  int N6,n6_1,n6_2, l31,l32,L31,L32;
   cgraph * cg=*c;
   

   N6=findl2(0,n6,cg,NULL);
   if(cg->vl[N6].vt !=V633) return 0;

#  if (CDEBLEV > DEBLEV) 
      printf("exp6(%d)",n6);  wrcg(cg);
#  endif 

   
   n6_1=findl2(1,N6,cg,NULL); 
   n6_2=findl2(2,N6,cg,NULL); 
      
   if(n6_1==n6 && n6_2==n6) 
   { cg->powNp1++; cg->powN++; cg->pow2--;
     cg->vl[N6].vt=zv; cg->vl[n6].vt=zv; cg->gl-=2;
#  if (CDEBLEV > DEBLEV) 
      wrcg(cg); 
#  endif 
     return 1;
   }
   
   l31=  cg->vl[N6].e[1];
   l32=  cg->vl[N6].e[2];
   L31= cg->vl[n6].e[1];
   L32= cg->vl[n6].e[2];
   cg->pow2--;
   
   cg->vl[n6].vt=t33;     cg->vl[N6].vt=t33;
   cg->vl[n6].e[0]=0;     cg->vl[N6].e[0]=0;
   cg->vl[n6].e[1]=L31;   cg->vl[N6].e[1]=L32;       
   cg->vl[n6].e[2]=l31;   cg->vl[N6].e[2]=l32;
   
   cg=addcg(c); 
   cg->vl[n6].e[1]=L32;   cg->vl[N6].e[1]=L31; 
   cg->vl[n6].e[2]=l31;   cg->vl[N6].e[2]=l32;
   
#  if (CDEBLEV > DEBLEV) 
      wrcg(cg); 
      wrcg(cg->next); 
#  endif 
   return 1;
}  


static int rem333(int n0, cgraph ** c)
{ 
 int  n1,l01,l11,l02,l12,i; 
 cgraph * cg=*c;


#  if (CDEBLEV > DEBLEV) 
      printf("rem333(%d)",n0); wrcg(cg);
#  endif 

   n1 = findl2(0,n0,cg,&i);
   
   if(cg->vl[n1].vt!=V333) return 0;
   
   if(i)
   {   cg->vl[n1].e[i]=cg->vl[n1].e[0];
       cg->vl[n1].e[0]=cg->vl[n0].e[0];
       cg->sgn*=-1;
   }
   
   l01=cg->vl[n0].e[1];  l02=cg->vl[n0].e[2];
   l11=cg->vl[n1].e[1];  l12=cg->vl[n1].e[2];
   
   
   cg->vl[n0].vt = t33;  cg->vl[n0].e[0]=0; 
   cg->vl[n1].vt = t33;  cg->vl[n1].e[0]=0; 
   
   cg->vl[n0].e[1]=l11; cg->vl[n0].e[2]=l01;
   cg->vl[n1].e[1]=l12; cg->vl[n1].e[2]=l02; 
       
   cg=addcg(c); 
   cg->sgn*=-1;
   
   cg->vl[n0].e[1]=l12; 
   cg->vl[n1].e[1]=l11;

#  if (CDEBLEV > DEBLEV) 
      wrcg(cg); 
      wrcg(cg->next); 
#  endif 
   return 1;
}  /* Exp3Q */ 


factor *colorFactor (int nv, cvertex * vl)
{  /* - calculate color weight (two int n,d) - 08/01/90  */
   cgraph    * cg;
   factor *f=(factor *)malloc(sizeof(factor));
   fct_init(f);


   cg =initcg(nv,vl);

#  if (CDEBLEV > DEBLEV) 
    printf("INPUT:");   wrcg(cg); 
#  endif 
   while(cg)
   { 
      cgraph  * pgl;
      
      while (cg->gl && cg->sgn)
      {  int n0,ng, l;    
         remtv(cg); if(!cg->gl) break; 
if(findv2(vFabc,cg,&n0)) 
if(istadpole2(n0,cg)) break;  else {expFabc(n0,&cg); continue; }
         
         
                         
         if(findv2(v833,cg,&n0)) 
         {          
         if(istadpole2(n0,cg)) break;  
            ng=findl2(0,n0,cg,&l);
            if(istadpole2(ng,cg)) break;
            switch(cg->vl[ng].vt)
            { case v833:  rem833_833(n0,ng,&cg); continue;
              case vFabc: rem833_Fabc(n0,ng,&cg);continue;
            }  
         }
         
         if(findv2(v866,cg,&n0))
         { 
           if(istadpole2(n0,cg))  break;
              
           ng=findl2(0,n0,cg,&l);           
           if(istadpole2(ng,cg))  break;
            
           switch(cg->vl[ng].vt)
           {
             case vFabc:  rem866_Fabc(n0,ng,&cg); continue;  
//             case v866:   rem866_866(n0,ng,&cg); continue; 

             default: exp866(n0,&cg);  continue; 
           }  
         }
          
         if(findv2(v633,cg,&n0)) { if(exp6(n0,&cg)) continue;}
         
         if(findv2(vFabc,cg,&n0)) 
          if(istadpole2(n0,cg)) break;  else {expFabc(n0,&cg); continue; }
                   
         if(findv2(v333,cg,&n0)) { if(rem333(n0,&(cg))) continue; }

         wrcg(cg); cerror(251,"No rules to expend diagram");
      } 
      pgl = cg;
      if(pgl->sgn)  add_fct(f,pgl->sgn,pgl->pow2,pgl->powNm1,pgl->powN,pgl->powNp1);
      cg = cg->next; 
      free(pgl->vl);  free(pgl); 
   }  
   return f;
} 


static void rednd(long * n,long * d,int b)
{ /* - reduce N and D with respect to B - 08/01/90  */
   if (b != 1)
      while (*d != 1 && *n % b == 0 && *d % b == 0)
      {
         *n /= b;
         *d /= b;
      }  /* while */
} /* RedND */


void fct_num_calc(factor * fct,int Nc, long * n, long *d)
{ 
	int i;
	long p;
	for(i=0,*n=0,p=1;i<fct->len;i++) { *n+=p*fct->nc[i];p*=Nc;}
  
	*d=1;
	if(fct->pow2>0) *n*=1<<fct->pow2; else *d*=1<<(-fct->pow2);
	
	if(fct->powN>0) for(i=0;i<fct->powN;i++) *n*=Nc;
	else for(i=0;i<-fct->powN;i++) *d*=Nc;		
	rednd(n,d,Nc); 
	rednd(n,d,2);
}

void fct_print(factor *fct, char *s)
{ int i;

  if(!fct->len)
  {
	  strcpy(s,"(0)"); 
	  return;
  }

  sprintf(s,"(");
  for(i=0;i<fct->len;i++)  
	  if(fct->nc[i]) 
	  {  
		  sprintf(s+strlen(s),"%+ld",fct->nc[i]);
		  if(i) sprintf(s+strlen(s),"N^%d",i);
	  }
  sprintf(s+strlen(s),")");

  if(fct->powN)  sprintf(s+strlen(s),"*(N^%ld)",fct->powN);
  if(fct->pow2<0) sprintf(s+strlen(s),"*2^%d",fct->pow2);
}


vtype typev(vert0 v,int valence)
{ int n, color=1;
 
/* Return color type of vertex - 06/01/90  */ 
      
   for (n = 0; n < valence; n++) switch(prtclbase[v[n].partcl-1].cdim)
   {
      case  8: color*=2; break;
      case  3: color*=3; break;
      case -3: color*=5; break;
      case  6: color*=7; break;
      case -6: color*=11;break;  
   }
   
{ int i, permit[]={1,4,8,15,27,30,77,99,125,175,154};
  for(i=0;i<11;i++) if(color==permit[i]) break;
  if(i==11) { printf("unknown colour vertex\n");} 
}   
   return color;
   
}  /* TypeV */ 


int t2k2(vcsect* g, int * nv, cvertex * vl)
{
  int   i,k, ne, en=0, maxnv=*nv; 
  int   maptar[2*maxvert][MAXVALENCE];
    
/* Transfer Taranov's representation of graph to Kryukov's representation  */ 
/* - 07/01/90  */ 
    
   for(i=0;i<2*maxvert;i++) for(k=0;k<MAXVALENCE;k++) maptar[i][k]=0;
   *nv=0;
   for(i=0;i<g->sizet; i++)  if (typev(g->vertlist[i],g->valence[i]) != zv) 
   { 
      if(*nv>= maxnv) cerror(251,"To many vertices");
      vl[*nv].e[0]=0; vl[*nv].e[1]=0;  vl[*nv].e[2]=0;
      vl[*nv].vt = typev(g->vertlist[i],g->valence[i]); 

      for (k=0, ne=0; ne<g->valence[i] ;ne++)
      {  int dim= prtclbase[g->vertlist[i][ne].partcl-1].cdim;
         if(dim!=1)
         {  int l=maptar[i][ne];
            if(!l)
            { 
               l = ++en; 
               maptar[i][ne] = l;
               maptar[g->vertlist[i][ne].link.vno]
                     [g->vertlist[i][ne].link.edno]=l; 
            } 
            
            switch(vl[*nv].vt)
            { case  v833: case  t33: case v866: case t66: 
                switch (dim) 
                {
                  case  8:  vl[*nv].e[0]=l;  break; 
                  case -3:
                  case -6:  vl[*nv].e[1]=l;  break; 
                  case  3: 
                  case  6:  vl[*nv].e[2]=l;  break;
                 default: cerror(252,"t2k - invalid particle color");
                }
                break;
              case  V633: case v633:
                switch (dim) 
                {
                  case  6: case -6:  vl[*nv].e[0]=l;     break;
                  case  3: case -3: if(!k)k=1;  vl[*nv].e[k++] = l; break;
                }
                break;              
              default:  vl[*nv].e[k++] = l; break; // vFabc, t88 v333,V333
            } 
         }
      }  
      (*nv)++;
   } 
   int j,n, sgn=1;
   for(n=0;n<g->sizel; n++) if(vl[n].vt==vFabc || vl[n].vt==v333 || vl[n].vt==V333)
   for(i=0;i<g->valence[n]-1;i++) if(prtclbase1[g->vertlist[n][i].partcl].cdim!=1)
   for(j=i+1;j<g->valence[n];j++)  if(prtclbase1[g->vertlist[n][j].partcl].cdim!=1)
   {  int ni=g->vertlist[n][i].link.vno, nj=g->vertlist[n][j].link.vno;
      if(ni==nj && g->vertlist[n][i].link.edno > g->vertlist[n][j].link.edno) sgn*=-1;
   }
   return sgn; 
}
