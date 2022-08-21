#include <stdio.h>
#include <stdlib.h>

static int recurGen(int k,int rev,int np,int nap,int*pos3,int*pos_3,int*perm,int*used,int**w)
{
  int nc=np+(nap-np)/3; 
  int pow=0;
  int i,j,j1,j2;
  
  if(k==nc) // writing 
  {    
     for(i=0;i<np;i++) 
     { *((*w)++)=2; 
       if(rev){ *((*w)++)=pos3[i]; *((*w)++)=pos_3[perm[i]];} 
       else {*((*w)++)=pos_3[perm[i]]; *((*w)++)=pos3[i];}
       *((*w)++)=0;
     }     
     for(i=0;i<(nap-np)/3;i++) 
     { if(rev) *((*w)++)=3; else *((*w)++)=-3; 
       for(j=0;j<nap;j++) if(used[j]==-(i+1)) *((*w)++)=pos_3[j];
     }
     return 1; 
  }
  else if(k>=np)  // construction of eps chains
  { 
    int mark=np-k-1;

    for(j=0;j<nap-2;j++) if(used[j]==0 &&(pos_3[j]>0 || used[j-1] ))
    {  
      used[j]= mark;
      for(j1=j+1; j1<nap-1;j1++)if(used[j1]==0 &&(pos_3[j1]>0 || used[j1-1]))  
      { used[j1]= mark;
        for(j2=j1+1; j2<nap;j2++)if(used[j2]==0 &&(pos_3[j2]>0 || used[j2-1])
        
        && abs(pos_3[j])!=abs(pos_3[j1]) 
        && abs(pos_3[j])!=abs(pos_3[j2]) && abs(pos_3[j1])!=abs(pos_3[j2]) )
        
        { used[j2]=mark;  
          pow+=recurGen(k+1,rev,np,nap,pos3,pos_3,perm,used,w);
          used[j2]=0;
        }  
        used[j1]=0;
      }
      used[j]=0;
      break;     
    }
    return pow;
  }
  else  // construction of -3 3 chains
  { 
    if( pos3[k]<0) i=perm[k-1]+1; else i=0;
    for(; i< nap; i++) if(!used[i] && pos_3[i]!=pos3[k]
        && (pos_3[perm[i]]>0 || used[perm[i]-1]!=0))
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

   
   for(i=0;i<4*nc_*pow;i++) if((*chains)[i]<0) (*chains)[i]*=-1; 
   free(pos3), free(pos_3);   
   return pow;
}

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
