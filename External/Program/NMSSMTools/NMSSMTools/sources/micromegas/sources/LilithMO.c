#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"

int isSMP(int pdg)
{  
  int SMP[16]={1,2,3,4,5,6, 11,12,13,14,15,16, 21,22,23,24};
  int apdg=abs(pdg);
  for(int i=0;i<16;i++) if(apdg==SMP[i]) return 1;
  return 0;
}   
   
static int txt2plist(char*txt,double *br, int * ind)
{  
   if(0==sscanf(txt,"%lf",br)) return 0;
   char *chE, *chB=strstr(txt,"->");    
   if(!chB) return 0;
   chB+=2;
   int n;
   
   for(n=0, chE=chB; chE; n++, chE=strchr(chE+1,','))
   {  char buff[20];
         sscanf(chE+1,"%[^,]",buff);
         trim(buff); 
         ind[n]=pNum(buff);
   }
   return n;
}


int LilithMO(char * fname)  
{ int i,j,k;

//printf("LILITH_MO\n");
  
  double Mcp=1.478662E+00, Mbp=bPoleMass(), Mtp=tPoleMass();
  int VZdecay_save=VZdecay, VWdecay_save=VWdecay, changeVirtual=0;
      
  double coeff[10];

  char  *G =pdg2name(21);  if(!G) { printf("Gluon  is absent in the model");    return -1;} 
  char  *A =pdg2name(22);  if(!A) { printf("Photon is absent in the model");    return -1;}
  char  *Z = pdg2name(23); if(!Z) { printf("Z is absent in the model");         return -1;}
  char  *Wp= pdg2name(24); if(!Wp){ printf("W+ is absent in the model");        return -1;}
  char  *Wm=  antiParticle(Wp);
  double MW=pMass(Wp);
   
  lVert* WWA=getVertex(Wp,Wm,A,NULL); if(!WWA) { printf(" %s %s %s is absent\n", Wp,Wm,A); return -2;}
  getNumCoeff(WWA,coeff);
  for(k=0;k<WWA->nTerms;k++) if( strcmp(WWA->SymbVert[k],"p1.m1*m3.m2")==0) break;
  double EE=coeff[k];
//for(k=0;k<WWA->nTerms;k++) printf( "%e ", coeff[k]);
//printf("\n"); 
  
  lVert* WWZ=getVertex(Wp,Wm,Z,NULL); if(!WWZ) { printf(" %s %s %s is absent\n", Wp,Wm,Z); return -2;}
  getNumCoeff(WWZ,coeff);
  for(k=0;k<WWZ->nTerms;k++) if( strcmp(WWZ->SymbVert[k],"p1.m1*m3.m2")==0) break;
  double SW=sin(atan(EE/coeff[k]));
  double vev= 2*MW*SW/EE;

//  fprintf(f,"# SM basic parameters: ee=%E  MW=%.3E alphaEM^{-1}=%.3E SW^2=%.3E  VEV=%.3E \n",EE, MW,4*M_PI/EE/EE,SW*SW,vev);

// Find Higgses 

   char ** Higgs=NULL;
   int nHiggs=0;
   char ** chHiggs=NULL;
   int nHch=0;
   
   for(i=0;i<nModelParticles;i++) if(ModelPrtcls[i].name[0]!='~' && ModelPrtcls[i].spin2==0 
                                   && ModelPrtcls[i].cdim==1)
   { 
       double m=findValW(ModelPrtcls[i].mass);
       if(strcmp(ModelPrtcls[i].name,ModelPrtcls[i].aname)==0 && m>123 && m<128)
//       if(strcmp(ModelPrtcls[i].name,ModelPrtcls[i].aname)==0)
       {  Higgs=(char**)realloc(Higgs, (sizeof(char*))*(nHiggs+1));
          Higgs[nHiggs++]=ModelPrtcls[i].name;   
       }
   }
   if(nHiggs==0) return 0;

   FILE*f=fopen(fname,"w");
   if(VZdecay==0 || VWdecay==0) 
   { VZdecay=1; VWdecay=1; cleanDecayTable(); changeVirtual=1;}
      
   fprintf(f,"<?xml version=\"1.0\"?>\n<lilithinput>\n");
   for(i=0;i<nHiggs;i++)
   {  fprintf(f,"  <reducedcouplings part=\"%s\">\n",Higgs[i]);
      fprintf(f,"  <mass> %E </mass>\n",pMass(Higgs[i]));
      int ferm[4]={5,6,13,15};
      char* fn[4]={"bb","uu","mumu","tautau"};
      for(j=0;j<4;j++) 
      {  char *fe=pdg2name(ferm[j]); if(!fe) continue;
         char *Fe=antiParticle(fe);  if(!Fe) continue;
         double m=pMass(fe);
         if(m==0) continue; 
         lVert *hff=getVertex(Fe,fe,Higgs[i],NULL); 
         double c[2]={0,0};
         if(hff)
         {  getNumCoeff(hff,coeff); 
            for(k=0;k<hff->nTerms;k++)
            { if(strcmp(hff->SymbVert[k],"1")==0)    c[0]=-coeff[k]*vev/m; else 
              if(strcmp(hff->SymbVert[k],"G5*i")==0) c[1]=coeff[k]*vev/m;
            }
         }
//printf("Fe,fe,Higgs = %s %s %s\n", Fe,fe,Higgs[i]);
//printf("nTerms=%d\n", hff->nTerms);
//         for(int i=0;i<hff->nTerms;i++)   printf(" %s %e\n", hff->SymbVert[i], coeff[i]);  

         if(c[1]==0 ) fprintf(f,"       <C to=\"%s\" >  %E </C>\n", fn[j],c[0]);
         else 
         {  if(c[0]) fprintf(f,"       <C to=\"%s\" part=\"re\"> %E </C>\n",fn[j], c[0]); 
            if(c[1]) fprintf(f,"       <C to=\"%s\" part=\"im\"> %E </C>\n",fn[j], c[1]); 
         }       
      }   
 
      char *ve=pdg2name(24); if(!ve) continue;
      char *Ve=antiParticle(ve); if(!Ve) continue; 
      lVert *hvv=getVertex(Ve,ve,Higgs[i],NULL);          
      double c=0;
      if(hvv)
      {  getNumCoeff(hvv,coeff); 
         for(k=0;k<hvv->nTerms;k++) if(strcmp(hvv->SymbVert[k],"m2.m1")==0) c=0.5*coeff[k]*vev/pMass(ve)/pMass(ve);  
      } 
      fprintf(f, "    <C to=\"VV\">%f</C>\n", c);
         
      double complex ffE,faE,ffC,faC; 
      double mH=pMass(Higgs[i]); 
      HiggsLambdas(mH, Higgs[i], &ffE,&ffC, &faE, &faC);      
      double  LGGSM=lGGhSM(mH,alphaQCD(mH)/M_PI, Mcp,Mbp,Mtp,vev);
      double  LAASM=lAAhSM(mH,alphaQCD(mH)/M_PI, Mcp,Mbp,Mtp,vev);
      double  lamQGG=(ffC*conj(ffC) +4*faC*conj(faC));
      double  lamQAA=(ffE*conj(ffE) +4*faE*conj(faE));
      fprintf(f, "    <C to=\"gammagamma\">%E</C>\n", -sqrt(lamQAA)/LAASM);
      fprintf(f, "    <C to=\"gg\">%E</C>\n",  -sqrt(lamQGG)/LGGSM);

      double  BrAA=0, BrGG=0,BrDM=0,BrNoSM=0;
      double width,widthP;
      int idDM1=0,idDM2=0;
      if(CDM1) idDM1=abs(pNum(CDM1));
      if(CDM2) idDM2=abs(pNum(CDM2));
      txtList L,l;
      width=pWidth(Higgs[i],&L);
      for(l=L;l;l=l->next)
      {  double br;
         int id[10],nd;
         nd=txt2plist(l->txt,&br,id); 
         if(nd==2)
         { if(id[0]==21 && id[1]==21) BrGG=br;
           if(id[0]==22 && id[1]==22) BrAA=br;
           if((abs(id[0])==idDM1 &&  abs(id[1])==idDM1) || (abs(id[0])==idDM2 && abs(id[1])==idDM2)) BrDM+=br;
         }
         for(k=0;k<nd;k++) if(!isSMP(id[k])) { BrNoSM+=br; break;}   
      }
      
      double K=1;
      if(BrGG==0)
      { double a=alphaQCD(mH)/M_PI;
        double Rqcd=(1+a*(149./12.+a*(68.6482-a*212.447)));
        double wGG=2*Rqcd*pow(mH,3)*lamQGG/M_PI;
        K+=wGG/width;
      }
      if(BrAA==0) 
      {
        double  wAA=0.25*pow(mH,3)*lamQAA/M_PI;
        K+=wAA/width;
      }
      fprintf(f, "    <precision>%s</precision>\n", "BEST-QCD");
      fprintf(f, "    <extraBR>\n");
      fprintf(f, "      <BR to=\"invisible\">%f</BR>\n", BrDM/K);
      fprintf(f, "      <BR to=\"undetected\">%f</BR>\n",BrNoSM/K);
      fprintf(f, "    </extraBR>\n");
      fprintf(f, "  </reducedcouplings>\n");
    }                    
    fprintf(f, "</lilithinput>\n");
    fclose(f);
    if(changeVirtual)
    { VZdecay=VZdecay_save; VWdecay=VWdecay_save; cleanDecayTable();}
          
    return nHiggs;
}  

int lilithmo_(char*fname,int len)
{
   char * cname=malloc(len+2);
   fName2c(fname,cname,len);
   int nHiggs= LilithMO(cname);   
   free(cname);
   return nHiggs; 
}


