#include "../include/micromegas.h"
#include "../include/micromegas_aux.h"

static char *VERSION="4.0";

int callSuperIsoSLHA(void)
{ 
   int err;
   char * command=malloc(strlen(micrO)+100);
   sprintf(command,"%s/Packages/superiso_v%s/slha.x",micrO,VERSION);

    
   if(access(command,X_OK))
   { 
     int sysTimeLimTmp=sysTimeLim;
     sysTimeLim=0;
     System("make -C %s/Packages/ -f SUPERISO.makef VERSION=%s",micrO,VERSION);
     sysTimeLim=sysTimeLimTmp;
   }

   if(access(command,X_OK)) 
   { printf("Can not compile superIso\n");
     free(command);
     return 1;
   }
   slhaWrite("slhaForSuperIso");
   System("%s  slhaForSuperIso >/dev/null  ",command);
//System("%s  slhaForSuperIso",command);
   err=slhaRead("output.flha",1);
//   printf("command-%s\n", command);
   free(command);
   return err;
}

int callsuperisoslha_(void){ return callSuperIsoSLHA();}
