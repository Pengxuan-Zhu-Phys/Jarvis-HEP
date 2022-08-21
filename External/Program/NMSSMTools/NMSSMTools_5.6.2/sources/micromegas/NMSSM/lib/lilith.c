#include"../../include/micromegas.h"
#include"../../include/micromegas_aux.h"
#include "pmodel.h"

int LilithMDL(char*fname)
{
  unsigned int i, npart=0;

  double CU, Cb, Ctau, CV, Cgamma, Cg;

  char *parts[5]={"h1","h2","h3","ha","hb"};
  FILE*f; 
  f=fopen(fname,"w");
  
  fprintf(f, "<?xml version=\"1.0\"?>\n");
  fprintf(f, "<lilithinput>\n");

  for(i=0; i<5; i++) 
  {
    double mass = pMass(parts[i]);
    if(mass < 123 || mass > 128) continue;
    
    
    ++npart;

    // compute invisible and undetected branching ratios
    double invBR = 0., undBR = 0.;
    double w;
    txtList L;
    w=pWidth((char*)parts[i], &L);

    if(Mcdm1 < 0.5*mass) {
      char invdecay[50];
      char cdmName[50];
//      sortOddParticles(cdmName);
      strcpy(invdecay, CDM1);
      strcat(invdecay, ",");
      strcat(invdecay, CDM1);
      invBR = findBr(L, invdecay);
    }
    undBR = 1 - invBR - findBr(L, "b B") - findBr(L, "c C") - findBr(L, "l L") -
            findBr(L, "W+ W-") - findBr(L, "A A") - findBr(L, "Z Z") -
            findBr(L, "G G") - findBr(L, "m M") - findBr(L, "A Z") -
            findBr(L, "u U") - findBr(L, "d D") - findBr(L, "s S") - findBr(L, "t T");

    CU =   slhaVal("REDCOUP",0.,2,1+i,1);
    Cb =   slhaVal("REDCOUP",0.,2,1+i,3);
    Ctau = Cb;
    CV =   slhaVal("REDCOUP",0.,2,1+i,4);
    Cg=    slhaVal("REDCOUP",0.,2,1+i,5);
    Cgamma=slhaVal("REDCOUP",0.,2,1+i,6);
    
    
    fprintf(f, "  <reducedcouplings part=\"%s\">\n", parts[i]);
    fprintf(f, "    <mass>%f</mass>\n", mass);
    fprintf(f, "    <C to=\"uu\">%f</C>\n", CU);
    fprintf(f, "    <C to=\"bb\">%f</C>\n", Cb);
    fprintf(f, "    <C to=\"mumu\">%f</C>\n", Ctau);
    fprintf(f, "    <C to=\"tautau\">%f</C>\n", Ctau);
    fprintf(f, "    <C to=\"VV\">%f</C>\n", CV);
    fprintf(f, "    <C to=\"gammagamma\">%f</C>\n", Cgamma);
    fprintf(f, "    <C to=\"gg\">%f</C>\n", Cg);
//    fprintf(f, "    <C to=\"Zgamma\">%f</C>\n", 1.);
    fprintf(f, "    <precision>%s</precision>\n", "BEST-QCD");
    fprintf(f, "    <extraBR>\n");
    fprintf(f, "      <BR to=\"invisible\">%f</BR>\n", invBR);
    fprintf(f, "      <BR to=\"undetected\">%f</BR>\n", undBR);
    fprintf(f, "    </extraBR>\n");
    fprintf(f, "  </reducedcouplings>\n");
  }

  fprintf(f, "</lilithinput>\n");
  fclose(f);
  return npart;
}

int lilithmdl_(char*fname,int len)
{
   char * cname=malloc(len+2);
   int err;
   fName2c(fname,cname,len);
   err= LilithMDL(cname);
   free(cname);
   return err; 
}
