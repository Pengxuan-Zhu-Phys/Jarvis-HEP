#include"../../include/micromegas.h"
#include"../../include/micromegas_aux.h"
#include"pmodel.h"
#include"lpath.h"

#include<sys/wait.h>
#include<unistd.h>

#define FIN  "nngg.in"
#define FOUT "nngg.out"

int loopGamma(double * csAA, double *csAZ)
{
  double sigmav;
  char buff[2000];
  int err;
  FILE*f;
  int GI=0;
   
  *csAA=0,*csAZ=0; 

  if(!access(FOUT,R_OK)) unlink(FOUT);
  
  sprintf(buff, LPATH "/nngg-nm4b/lGamma.exe");
  if(access( buff,X_OK))
  { char buf[2000]; 
    sprintf(buf, "dir=%s/nngg-nm4b;  which  gmake; if(test $? -eq 0) then  gmake -C $dir; else make -C $dir; fi",LPATH);
    system(buf);
  } 
  if(access( buff,X_OK)) 
  {  
    printf("Cannot find/compile executable %s\n",buff);
    return 10;
  }  

  sprintf(buff+strlen(buff)," %E", vRot/299792.*1.5957691);

  err=System(buff);
  double ee=findValW("EE");
  double aFac=4*M_PI/ee/ee/137.036;
     

  if(err>=0) 
  {  err=slhaRead(FOUT,1);
     if(err) return err;
     *csAZ=slhaVal("Lgamma",0.,1,1)*2.9979E-26*aFac;
     *csAA=slhaVal("Lgamma",0.,1,2)*2.9979E-26*aFac*aFac;
  }  

//  if(!access(FOUT,R_OK)) unlink(FOUT);
//  if(!access(FIN,R_OK)) unlink(FIN);
  return err;
}  

extern int  loopgamma_( double * cs1, double *cs2); /* fortran */
int  loopgamma_(double * cs1, double *cs2){ return loopGamma(cs1,cs2); }
 
