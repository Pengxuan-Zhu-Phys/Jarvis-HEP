/*
 Copyright (C) 1997, Alexander Pukhov
*/

#include "physics.h"
#include "process.h"
#include "syst2.h"
#include "colorf.h"
#include "cweight.h"

int NcInfLimit=0;
int NcInfCC=0;

static int a2k(vampl * g, int nc, int * chains, int * nv, cvertex * vl)
{
  int   i,k, ne, en=0, maxnv=*nv; 
  int   maptar[maxvert][MAXVALENCE];
  int   endl[MAXINOUT];
  int   endL[MAXINOUT];
  int   endc[MAXINOUT];  

  for(i=0;i<MAXINOUT;i++){endl[i]=endL[i]=endc[i]=0;}
    
  for(i=0;i<maxvert;i++) for(k=0;k<MAXVALENCE;k++) maptar[i][k]=0;
  *nv=0;
  for(i=0;i<g->size; i++)  if (typev(g->vertlist[i],g->valence[i]) != zv) 
  { 
     if(*nv>= maxnv) return 1;
     vl[*nv].e[0]=0; vl[*nv].e[1]=0;  vl[*nv].e[2]=0;
     vl[*nv].vt = typev(g->vertlist[i],g->valence[i]); 
     for (k=0, ne=0; ne<g->valence[i] ;ne++)
     {  int dim= prtclbase[g->vertlist[i][ne].partcl-1].cdim;
        if(dim!=1)
        {  int l=maptar[i][ne];
           if(!l)
           { 
               l = ++en; 
               if(g->vertlist[i][ne].link.vno!= nullvert)
               {  maptar[i][ne] = l;
                  maptar[g->vertlist[i][ne].link.vno]
                     [g->vertlist[i][ne].link.edno]=l;
               }
               else 
               {  int  n=g->vertlist[i][ne].link.edno;
                  
                  if(dim==-3) endL[n]=l; else endl[n]=l;
                  endc[n]=dim; 
               }
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
                default: return 2;
              }
              break;
              case  V633: case v633:
                switch (dim) 
                {
                  case  6: case -6:  vl[*nv].e[0]=l;     break;
                  case  3: case -3: if(!k)k=1;  vl[*nv].e[k++] = l; break;
                }
                break;              
              default:  vl[*nv].e[k++] = l;  break; // vFabc, t88 v333,V333
           }  
         } 
      }  
      (*nv)++;   
  }

  for(i=0; i<MAXINOUT;i++) if(endc[i]==8 || abs(endc[i])==6)
  {  if(*nv>= maxnv) return 1; 
     vl[*nv].e[0]=endl[i];
     vl[*nv].e[1]=++en;
     vl[*nv].e[2]=++en;         

     switch(endc[i])
     { case 8:      
         vl[*nv].vt=v833;
         endl[i]=vl[*nv].e[2];
         endL[i]=vl[*nv].e[1];
       break;
    
       case 6:  
         vl[*nv].vt=V633;
         endl[i]=vl[*nv].e[1];
       break;

       case -6:
         vl[*nv].vt=v633;
         endL[i]=vl[*nv].e[1];
       break;
     }
     (*nv)++; 
  }


  for(i=0;i<nc;i++) if(abs(chains[4*i]==2))
  {  int from=chains[4*i+2]-1;
     int to  =chains[4*i+1]-1;  
     if(*nv>= maxnv) return 1;

     vl[*nv].vt=t33;
     vl[*nv].e[0]=0; 

     vl[*nv].e[1]=endl[from]++; 
     vl[*nv].e[2]=endL[to]++; 
     (*nv)++;
  } else 
  { int j;
     if(*nv>= maxnv) return 1;

     if(chains[4*i]==3)
     {  vl[*nv].vt=v333;
        for(j=0;j<3;j++) 
        { int m=chains[4*i+1+j]-1;
          vl[*nv].e[j]=endL[m]++; 
        }  
     } else 
     {
        vl[*nv].vt=V333;
        for(j=0;j<3;j++)
        { int m=chains[4*i+1+j]-1;
           vl[*nv].e[j]=endl[m]++; 
        }  
     }   
     (*nv)++; 
  }  
  
  return 0;
}  


static int maxNcPower(vcsect* g)
{
  int i,j;
  int n3=0,n_3=0;

  for(i=0;i<g->sizel; i++)
  for(j=0;j<g->valence[i];j++) if (g->vertlist[i][j].link.vno >=g->sizel)
  { int np=g->vertlist[i][j].partcl;
    switch(prtclbase[np-1].cdim)
    { case  3: n3++;       break; 
      case -3: n_3++;      break;
      case  6: n3+=2;      break;
      case -6: n_3+=2;     break;   
      case  8: n3++;n_3++; break;
    }
  }
   
  return (n3 + n_3 - abs(n3-n_3))/2 + 2*abs(n3-n_3)/3;
}

static void getLeadingTerm(factor * f, int  maxP, long *num, long *den)
{

//   if(maxP<f->len - f->dpow -1) {printf("BUG in my BRAIN \n"); sortie(58);}

   if( f->len==0 || maxP > f->len + f->powN -1) {*num=0; *den=1;} else
   {  *num=f->nc[f->len-1];
      if(f->pow2>0) {*num*=1<<f->pow2; *den=1;}  else  *den=1<<(-f->pow2);
      while (maxP--) (*num) *=3;
   }
}  


void c_basis_coef(vampl * g,int pow,int nc,int * chains,long * num,long * den)
{ int i,nv;
  cvertex   vl[3*MAXINOUT];
  factor * f;

  if(!pow) return;

  for(i=0; i<pow; i++)
  {  int maxP,j;
     nv=3*MAXINOUT;
     a2k(g, nc, chains +4*nc*i, &nv, vl);
     for(maxP=0,j=0;j<nc;j++) if(chains[4*(nc*i+j)]==2) maxP++; else maxP+=2;

     f=colorFactor(nv,vl);

     if(NcInfLimit)  getLeadingTerm(f,maxP,num+i,den+i);
     else fct_num_calc(f,3,num+i, den+i);
 
     free(f);
  }
}


void cwtarg(vcsect* g)
{  factor * f;
   cvertex   vl[2*MAXINOUT];
   int nv=2*MAXINOUT;
   int sgn=t2k2(g,&nv,vl);
   f=colorFactor(nv,vl);

   if(NcInfLimit) getLeadingTerm(f,maxNcPower(g),&(g->clrnum),&(g->clrdenum));
   else fct_num_calc(f,3,&(g->clrnum), &(g->clrdenum));
   
   g->clrnum*=sgn;
   free(f);
}  /* CWTarG */


static void  lreduce(long * l1, long  * l2)
{  long    c, i1, i2;

   i1 = *l1; i2 = *l2;
   if(i1<0) i1=-i1; if(i2<0) i2=-i2;

   if (i2 > i1) { c = i1; i1 = i2; i2 = c;}
   while (i2 != 0) { c = i2; i2 =i1%i2; i1 = c; }
   (*l1) /= i1;
   (*l2) /= i1;
}


int generateColorWeights(csdiagram*csdiagr,int cBasisPower,int nC,int*cChains,
     long * cCoefN,long * cCoefD)
{
   vcsect vcs;
   int NcInfLimit_tmp=NcInfLimit;
   int i;
   int inf;

   transfdiagr(csdiagr,&vcs);
   cwtarg(&vcs);

   if(NcInfLimit) inf=1; else inf=NcInfCC;
      
   NcInfLimit=inf;
      
   if(vcs.clrnum) 
   {  vampl left, right;
      long n=1,d=1;
      long * cCoefNr=malloc(cBasisPower*sizeof(long));
      long * cCoefDr=malloc(cBasisPower*sizeof(long));
 
      decompose(vcs,&left,&right);
      c_basis_coef(&left,cBasisPower,nC,cChains,cCoefN,cCoefD);
      c_basis_coef(&right,cBasisPower,nC,cChains,cCoefNr,cCoefDr);

      for(i=0;i<nC;i++) if(abs(cChains[4*i])==2)  d*=3; else {if(inf) d*=9; else d*=6;}
      for(i=0;i<nin+nout;i++) 
      { vertlink Q=left.outer[i];
        switch(prtclbase[left.vertlist[Q.vno][Q.edno].partcl-1].cdim)
        { case 8: n*=2; break;
          case 6: case -6: if(!inf){n*=3; d*=2;} break;
        }
      }
      
      for(i=0;i<right.size;i++)
      { int n8=0;
        int j;

        for(j=0;j<right.valence[i];j++)
        if(8==prtclbase[right.vertlist[i][j].partcl-1].cdim) n8++;
        if(n8==3) n*=-1;
      }
      
      for(i=0;i<cBasisPower;i++)
      { 
         cCoefN[i] *= cCoefNr[i]*n *vcs.clrdenum;
         cCoefD[i] *= cCoefDr[i]*d *vcs.clrnum;
         lreduce(cCoefN+i, cCoefD+i);
         if(cCoefN[i]<0){ cCoefN[i]*=-1;cCoefD[i]*=-1;}
      }
      free(cCoefNr); free(cCoefDr);       
   } else for(i=0;i<cBasisPower;i++){cCoefN[i]=0; cCoefD[i]=0;}

   NcInfLimit=NcInfLimit_tmp;
   return vcs.clrnum;
}


static int recurGen(int k,int rev,int np,int nap,int*pos3,int*pos_3,int*perm,int*used,int**w)
{
  int nc=np+(nap-np)/3; 
  int pow=0;
  int i,j0,j1,j2;
  
  if(k==nc) // writing 
  {    
     if(rev) for(i=0;i<np;i++) 
     { *((*w)++)=2; 
       *((*w)++)=pos3[i]; *((*w)++)=pos_3[perm[i]];
       *((*w)++)=0;
     } else 
     {int  f[MAXINOUT],s[MAXINOUT];
       for(i=0;i<np;i++){ f[i]=perm[i];s[i]=i;}
       i=0;
       while(i<np-1)
       { if(f[i]>f[i+1] || (f[i]==f[i+1] && s[i]>s[i+1]))
         { int b=f[i]; f[i]=f[i+1];f[i+1]=b;
               b=s[i]; s[i]=s[i+1];s[i+1]=b; 
           if(i) i--;else i++;
         }else i++;   
       }
       for(i=0;i<np;i++)
       { *((*w)++)=2;
         *((*w)++)=pos_3[f[i]]; *((*w)++)=pos3[s[i]];
         *((*w)++)=0;
       }  
     }    
     for(i=0;i<(nap-np)/3;i++) 
     { int j;
       if(rev) *((*w)++)=-3; else *((*w)++)=3;
       for(j=0;j<nap;j++) if(used[j]==-(i+1)) *((*w)++)=pos_3[j];
     }
     return 1; 
  }
  else if(k>=np)  // construction of eps chains
  { 
    int mark=np-k-1;

    for(j0=0;j0<nap-2;j0++) if(!used[j0])
    {  int prev[3]={0,0,0};
       if(pos_3[j0]<0 && mark<-1 &&  used[j0-1]==mark+1)
       { int k=0,l;
         for(l=0;l<nap;l++) if(used[l]==mark+1)prev[k++]=abs(pos_3[l]); 
       }
       if(abs(pos_3[j0])>prev[0]) {prev[1]=prev[2]=0;}        
       used[j0]= mark;
       for(j1=j0+1; j1<nap-1;j1++)if(!used[j1] &&(pos_3[j1]>0 || used[j1-1])&& abs(pos_3[j1])>=prev[1] )  
       { used[j1]= mark;
         if(abs(pos_3[j1])>prev[1])prev[2]=0;          
         for(j2=j1+1; j2<nap;j2++) if(!used[j2] &&(pos_3[j2]>0 || used[j2-1])&& abs(pos_3[j2])>=prev[2])
         if( abs(pos_3[j0])!=abs(pos_3[j1]) && abs(pos_3[j1])!=abs(pos_3[j2]))
         { used[j2]=mark;   
           pow+=recurGen(k+1,rev,np,nap,pos3,pos_3,perm,used,w);
           used[j2]=0;
         }   
         used[j1]=0;
       }
       used[j0]=0;
       break;    
    }
    return pow;
  }
  else  // construction of -3 3 chains
  { 
    if( pos3[k]<0) i=perm[k-1]+1; else i=0; 
    for(; i< nap; i++) if(!used[i] && pos_3[i]!=pos3[k]
        && (pos_3[i]>0 || used[i-1]!=0))
    {    
       used[i]=1;
       perm[k]=i;
       pow+=recurGen(k+1,rev,np,nap,pos3,pos_3,perm,used,w);
       used[i]=0;
    }
    return pow;
  }
}


static int Facto(int N)
{  int res=1,i;
   for(i=2;i<=N;i++) res*=i;
   return res;
} 

int infCbases(int np, int * cweight, int *nc, int ** chains)
{
   int n3=0, n_3=0, nc2, nc3;
   int i,nc_,rev,pow;
   int *used, *perm,*pos3,*pos_3,*wrt;
   pos3 =(int *)malloc(2*np*sizeof(int));
   pos_3=(int *)malloc(2*np*sizeof(int));

   for(i=1;i<=np;i++) switch (cweight[i-1])
   { 
      case  1:  break;
      case  3:  pos3[n3++]=i;   break; 
      case -3:  pos_3[n_3++]=i; break;
      case  6:  pos3[n3++]=i;   pos3[n3++]=-i;   break;
      case -6:  pos_3[n_3++]=i; pos_3[n_3++]=-i; break; 
      case  8:  pos3[n3++]=i;   pos_3[n_3++]=i;  break;
      default: {free(pos_3),    free(pos3);      return -1;} 
   }

   if(n3<=n_3) { nc2=n3;  nc3=(n_3-n3)/3; rev=0; }
   else        { nc2=n_3; nc3=(n3-n_3)/3; rev=1; }

   if(abs(n3-n_3) !=3*nc3) { free(pos_3), free(pos3); return -1;} 
      
   nc_=nc2+nc3;
   if(nc_)
   {  int pow_=Facto(nc2+3*nc3)/Facto(nc3);
      for(i=0;i<nc3;i++) pow_/=6;   
      wrt=(int *)malloc(4*(nc_)*pow_*sizeof(int));
      *chains=wrt;

      used= (int *) malloc((nc2+3*nc3)*sizeof(int)); for(i=0;i<nc2+3*nc3;i++) used[i]=0;
      perm= (int *) malloc(nc2*sizeof(int));   
      if(rev)  pow=recurGen(0,rev,n_3,n3,pos_3,pos3,perm,used,&wrt);
      else     pow=recurGen(0,rev,n3,n_3,pos3,pos_3,perm,used,&wrt);
if(pow_<pow) printf("ERROR: pow_=%d pow=%d  - wrong memmory allocation\n",pow_,pow);
      *nc=nc_;
      free(perm);  free(used); 
   } else { *nc=0; pow=0; *chains=NULL; pow=0;}

   
   for(i=0;i<nc_*pow;i++)
   {
     if((*chains)[4*i+1]<0) (*chains)[4*i+1]*=-1;
     if((*chains)[4*i+2]<0) (*chains)[4*i+2]*=-1; 
   }
   free(pos3), free(pos_3);   
   return pow;
}

#ifdef MAIN

int main(int argv, char**argc) 
{
  int np=argv-1;
  int cweight[100];
  int nc,pow,*chains;
  int i,j;

  if(!np) {printf("input expected\n"); exit(1);}
  for(i=0;i<np;i++)  if( 1!= sscanf(argc[i+1],"%d",cweight+i)) { printf("input mistyping\n"); exit(2);} 
    
  pow=infCbases(np, cweight, &nc, &chains);

  if(pow<0) { printf(" Wrong set of colors\n"); return 1;}
  printf("Basis power=%d.  Number of color chains %d. \n",pow,nc);
  
  for(i=0;i<pow;i++)
  { 
    for(j=0;j<nc;j++)
    { int *c=chains+4*(i*nc+j);
      if(c[0]==2)  printf("( %2d %2d)",c[1],c[2]);
      if(abs(c[0])==3)  printf("( %2d %2d %2d)",c[1],c[2],c[3]);
    }
    printf("\n");  
  }
  return 0;    
}

#endif