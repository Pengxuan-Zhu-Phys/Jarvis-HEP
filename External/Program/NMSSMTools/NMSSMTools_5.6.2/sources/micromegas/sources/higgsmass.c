#include"micromegas_aux.h"
#include"micromegas.h"

void printHiggs(FILE* f)
{
   int i;
   fprintf(f,"Higgs masses and widths\n");
   for(i=0;i<nModelParticles;i++)
   { if (ModelPrtcls[i].name[0]!='~' && ModelPrtcls[i].spin2==0 && ModelPrtcls[i].cdim==1)
     { int j;
       double width=pWidth(ModelPrtcls[i].name, NULL);
       fprintf(f,"%7.7s  %7.2f %.2E",ModelPrtcls[i].name,
                    findValW(ModelPrtcls[i].mass), width);
       for(j=0;j<nModelVars;j++)if(strcmp(ModelPrtcls[i].width,varNames[j])==0)
       { printf("(%.2E)",width=varValues[j]); break;} 
        fprintf(f,"\n");
     }
   }  
}


void printPartialWidth(double width,txtList l,FILE*f)
{ 
   double br;
   char txt[50];  
   fprintf(f," Branching  Partial width   Channel\n");
   for(;l;l=l->next)
   {
       sscanf(l->txt, "%lf %[^\n]", &br,txt); 
       fprintf(f," %.3E  %.3E GeV   %s\n", br, br*width,txt); 
   }   
}   

