/*
 Copyright (C) 1997, Slava Ilyin
*/
#include"interface.h"
#include"sf_epa.h"
#include "crt_util.h"
#include "vp.h"
#include "syst.h"
#include<math.h>

    static char   q_max_txt[2][20] = { "100.", "100." };
    static REAL q_max_val[2]={100,100};     
    static REAL*q_max_ptr[2]={q_max_val,q_max_val+1};
    static int pc[2]={1,-1};

/* ********************************************************* */
/*  Equivalent photon approximation structure function.   * */
/*     Improved Weizsaecker-Williams formula              * */
/*      C.F.Weizsaecker, Z.Phys. 88 (1934) 612            * */
/*      E.J.Williams,    Phys.Rev. 45 (1934) 729          * */
/*                                                        * */
/*   V.M.Budnev et al., Phys.Rep. 15C (1975) 181          * */
/* ********************************************************* */

int p_epa__(int *pNum){ if(pNum[0]==22 && pNum[1]==0 ) return 1; else return 0; } 
 

void n_epa__(int i,char * name)
{ 
  char * pname;
  i--;
  switch (pc[i])
  { case -2: pname="mu^-"; break;
    case -1: pname="e^-";  break;
    case  1: pname="e^+";  break;
    case  2: pname="mu^+"; break;
  }
  sprintf(name,"Equiv.Photon(particle= %s |Q|max= %s)", pname, q_max_txt[i]);
}  

int r_epa__(int i, char *name)
{
  char pname[10];
  double v;
  i--;
  if(2!= sscanf(name,
    "Equiv.Photon(particle= %s |Q|max= %[^)]",pname,q_max_txt[i])) return 0;
  trim(q_max_txt[i]);  
  if(1==sscanf(q_max_txt[i],"%lf",&v))
  {  
     if (v <= 0.) return 0;
     q_max_val[i]=v;
     q_max_ptr[i]=q_max_val+i;  
  } else 
  { q_max_ptr[i]=varAddress(q_max_txt[i]);
    if(!q_max_ptr[i]) return 0;
  }
       if(strcmp(pname,"mu^-")==0) pc[i]=-2;
  else if(strcmp(pname,"e^-") ==0) pc[i]=-1;
  else if(strcmp(pname,"e^+") ==0) pc[i]= 1;
  else if(strcmp(pname,"mu^+")==0) pc[i]= 2;
  else return 0;   
  return 1;
}

int i_epa__(int i, double *be, double * mass) 
{ i--; 
  *be=1.; 
  if(abs(pc[i])==1)*mass=5.11E-4; else *mass=0.10566;  
  return 1;
}

int mc_epa__(int i)
{ i--;
  switch(pc[i])
  { case -2: return  13;
    case -1: return  11;
    case  1: return -11;
    case  2: return -13; 
  }
}

int m_epa__(int i, int * pString)
{   
    void * pscr = NULL;
    i--;
    for(;;)
    { char mstr[]="\007"
                  "mu^-   "
                  "e^-    "
                  "e^+    "
                  "mu^+   ";


      char strmen[]="\050"
      " Incoming particle XXX                  "
      " |Q|max = ZZZ                           ";
      int mode=1;
      int n;
      if(pc[i]<0) n=3+pc[i]; else n=2+pc[i];
   
      improveStr(strmen,"XXX","%4.4s",mstr+1+(n-1)*mstr[0]);
      improveStr(strmen,"ZZZ","%s GeV",q_max_txt[i]);
      menu1(38,10,"",strmen, "n_sf_epa",&pscr,&mode);
    
      switch(mode)
      {
      case 0: { double v;
                if(1==sscanf(q_max_txt[i],"%lf",&v))
                {  
                   q_max_val[i]=v;
                   q_max_ptr[i]=q_max_val+i;
                   return 1;
                } else 
                { trim(q_max_txt[i]);
                  q_max_ptr[i]=varAddress(q_max_txt[i]);
                  if(!q_max_ptr[i]) 
                  { char buff[50];
                    sprintf(buff,"Can not find parameter '%s'",q_max_txt[i]);
                    continue;   
                  }
                }
                return 1;
              }      
      
      case 1: menu1(38,13,"",mstr,"",NULL,&n);
              if(n<3) pc[i]=n-3; else pc[i]=n-2; 
              break;
      case 2: correctStr (40,16,"Enter new value ",q_max_txt[i],19,1);
 break;
      }
    }
    return 1;     
}


double c_epa__(int i, double x, double q)
{
#define alpha  .0072992701
  i--;
  {  double delt = (abs(pc[i])>1)? 0.10566:5.11E-4;
     double f;   
     delt/=fabs(*(q_max_ptr[i]));
     delt*=delt;
     f=  alpha/(M_PI * 2) 
     * (log((1 - x) / (x * x * delt)) * ((1-x) * (1-x) + 1) / x 
      - (1 - x - delt * (x * x)) * 2/x);
     if(f<0)  f=0;
     return f; 
  }
}
