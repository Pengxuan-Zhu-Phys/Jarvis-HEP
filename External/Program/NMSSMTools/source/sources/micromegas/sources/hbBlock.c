#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"

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



int hbBlocksMO(char *fname, int *nHiggsCh)  
{ int i,j,k;

  FILE*f=fopen(fname,"w");
  
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
//for(k=0;k<WWA->nTerms;k++) printf(" %s  ", WWA->SymbVert[k]); printf("\n"); 
  for(k=0;k<WWA->nTerms;k++) if( strcmp(WWA->SymbVert[k],"m2.m1*p2.m3")==0) break;
  double EE=coeff[k];

//printf("EE=%E\n",EE);
  
//for(k=0;k<WWA->nTerms;k++) printf( "%e ", coeff[k]);
//printf("\n"); 
  
  lVert* WWZ=getVertex(Wp,Wm,Z,NULL); if(!WWZ) { printf(" %s %s %s is absent\n", Wp,Wm,Z); return -2;}
  getNumCoeff(WWZ,coeff);
  for(k=0;k<WWZ->nTerms;k++) if( strcmp(WWZ->SymbVert[k],"m2.m1*p2.m3")==0) break;
  double SW=sin(atan(EE/coeff[k]));
  double vev= fabs(2*MW*SW/EE);
//printf("SW=%E vev=%E\n",SW,vev);
  fprintf(f,"# SM basic parameters: MW=%.3E alphaEM^{-1}=%.3E SW^2=%.3E  VEV=%.3E \n", MW,4*M_PI/EE/EE,SW*SW,vev);

// Find Higgses 

   char ** Higgs=NULL;
   int nHiggs=0;
   char ** chHiggs=NULL;
   int nHch=0;

   int modsel=writeBlock("MODSEL",f);
   int nmssm=0;   
   for(i=0;i<nModelParticles;i++) if(/*ModelPrtcls[i].name[0]!='~' &&*/ ModelPrtcls[i].spin2==0 
                                   && ModelPrtcls[i].cdim==1)
   { 
       int hList[6]={25,35,45,36,46,37};
       for(j=0;j<6;j++) if( abs(ModelPrtcls[i].NPDG) == hList[j]) break;
       
       if(j==6) 
       { if(ModelPrtcls[i].name[0]!='~')  printf(" HB-Block Warning:  Colorless Scalar even particle  %s (pdg=%d) is out of HiggsBounds list\n",ModelPrtcls[i].name, ModelPrtcls[i].NPDG);  
         continue;
       }
       if(j==2 || j==4) nmssm=1; 
       if(strcmp(ModelPrtcls[i].name,ModelPrtcls[i].aname)==0)
       {  Higgs=(char**)realloc(Higgs, (sizeof(char*))*(nHiggs+1));
          Higgs[nHiggs++]=ModelPrtcls[i].name;
       } else
       {  chHiggs=(char**)realloc(chHiggs, (sizeof(char*))*(nHch+1));
          chHiggs[nHch++]=ModelPrtcls[i].name;
       }   
   }

   int DMList[7]= {1000022, 1000012, 1000014, 1000016, 1000017, 1000018, 1000019};
   if(CDM1)
   { int ndm=abs(pNum(CDM1));
     for(int i=0;i<7;i++) if(ndm==DMList[i]) break;
     if(i==7) printf(" HB-Block Warning: DM %s is out of HiggsBounds list\n",CDM1);
   }      
   if(CDM2)
   { int ndm=abs(pNum(CDM2));
     for(int i=0;i<7;i++) if(ndm==DMList[i]) break;
     if(i==7) printf(" HB-Block Warning: DM %s is out of HiggsBounds list\n",CDM2);
   }      
   
   if(modsel==0)
   { fprintf(f,"BLOCK MODSEL  # written by micrOMEGAs\n");
     fprintf(f," 3    %d\n\n",nmssm); 
   }
   
   fprintf(f,"BLOCK MASS\n");
   for(i=0;i<nHiggs;i++) fprintf(f," %d  %E # %s\n", pNum(Higgs[i]), pMass(Higgs[i]),Higgs[i]);    
   for(i=0;i<nHch;i++)   fprintf(f," %d  %E # %s\n", pNum(chHiggs[i]), pMass(chHiggs[i]),chHiggs[i]);
 
   fprintf(f,"Block HiggsCouplingsFermions\n");
   int ferm[3]={5,6,15};
   for(i=0;i<nHiggs;i++) for(j=0;j<3;j++)
   {  char *fe=pdg2name(ferm[j]); if(!fe) continue;
      char *Fe=antiParticle(fe);  if(!Fe) continue; 
      double mf=pMass(fe); if(mf==0) continue;
      lVert *hff=getVertex(Fe,fe,Higgs[i],NULL); 
      double c[2]={0,0};
      if(hff)
      {  getNumCoeff(hff,coeff); 
         for(k=0;k<hff->nTerms;k++)
         { if(strcmp(hff->SymbVert[k],"1")==0)    c[0]=-coeff[k]*vev/mf; else 
           if(strcmp(hff->SymbVert[k],"G5*i")==0) c[1]=-coeff[k]*vev/mf;
         }
      }       
      fprintf(f," %12.4E  %12.4E  3  %3d  %3d  %3d  # %s %s %s\n", c[0], c[1], pNum(Higgs[i]), ferm[j], ferm[j], Higgs[i],Fe,fe); 
   }   
   
   fprintf(f,"Block HiggsCouplingsBosons\n");
   int vbm[2]={23,24};
   for(i=0;i<nHiggs;i++) for(j=0;j<2;j++)
   {  char *ve=pdg2name(vbm[j]); if(!ve) continue;
      char *Ve=antiParticle(ve); if(!Ve) continue; 
      lVert *hvv=getVertex(Ve,ve,Higgs[i],NULL);
//printf("bosons: %s %s\n", ve,Ve);             
      double c=0;
      if(hvv)
      {  getNumCoeff(hvv,coeff); 
         for(k=0;k<hvv->nTerms;k++) if(strcmp(hvv->SymbVert[k],"m2.m1")==0) c=0.5*coeff[k]*vev/pMass(ve)/pMass(ve);  
      }    
      fprintf(f," %12.4E  3  %3d  %3d  %3d  # %s %s %s\n", c, pNum(Higgs[i]), vbm[j], vbm[j], Higgs[i],Ve,ve); 
   }

   for(i=0;i<nHiggs;i++) for(j=i+1;j<nHiggs;j++)
   {  char *ve=pdg2name(23); if(!ve) continue;
      lVert *hhv=getVertex(Higgs[i],Higgs[j],ve,NULL);
      double c=0;
      if(hhv)
      {  getNumCoeff(hhv,coeff); 
         for(k=0;k<hhv->nTerms;k++) if(strcmp(hhv->SymbVert[k],"p1.m3*i")==0) c=coeff[k]*vev/pMass(ve); 
      }    
      fprintf(f," %12.4E  3  %3d  %3d  %3d  # %s %s %s\n", c, pNum(Higgs[i]),pNum(Higgs[j]) , 23, Higgs[i],Higgs[j] ,ve); 
   }

   for(i=0;i<nHiggs;i++) for(j=0;j<nHch;j++)
   {  int charge3;
      qNumbers(chHiggs[j],NULL, &charge3, NULL); 
      char *ve=NULL;
      if(charge3==3) ve=pdg2name(-24);
      else if(charge3==-3) ve=pdg2name(24);
      if(!ve) continue;
      lVert *hhv=getVertex(Higgs[i],chHiggs[j],ve,NULL);
      double c[2]={0,0};
      if(hhv)
      {  getNumCoeff(hhv,coeff); 
         for(k=0;k<hhv->nTerms;k++) if(strcmp(hhv->SymbVert[k],"p1.m3*i")==0) c[0]=coeff[k]*vev/pMass(ve);
                               else if(strcmp(hhv->SymbVert[k],"p1.m3")==0) c[1]=coeff[k]*vev/pMass(ve); 
      }    
      fprintf(f," %12.4E  3  %3d  %3d  %3d  # %s %s %s\n", (double)sqrt(c[0]*c[0]+c[1]*c[1]), pNum(Higgs[i]),pNum(chHiggs[j]) ,pNum(ve), Higgs[i],chHiggs[j] ,ve); 
   }
   
   
// Gluon  and photon decays  
   double *lamQGG      =(double *)malloc(nHiggs*sizeof(double));
   double *lamQAA      =(double *)malloc(nHiggs*sizeof(double));          
   for(i=0;i< nHiggs;i++)
   {  double complex ffE,faE,ffC,faC; 
      double mH=pMass(Higgs[i]);
   
      HiggsLambdas(mH,Higgs[i], &ffE,&ffC, &faE, &faC);
           
      double  LGGSM=lGGhSM(mH,alphaQCD(mH)/M_PI, Mcp,Mbp,Mtp,vev);
      double  LAASM=lAAhSM(mH,alphaQCD(mH)/M_PI, Mcp,Mbp,Mtp,vev);
              lamQGG[i]=(ffC*conj(ffC) +4*faC*conj(faC));
              lamQAA[i]=(ffE*conj(ffE) +4*faE*conj(faE));      
      fprintf(f," %12.4E  3  %3d  %3d  %3d  # %s %s %s\n", sqrt(lamQGG[i])/fabs(LGGSM), pNum(Higgs[i]) ,pNum(G), pNum(G), Higgs[i] ,G ,G);
      fprintf(f," %12.4E  3  %3d  %3d  %3d  # %s %s %s\n", sqrt(lamQAA[i])/fabs(LAASM), pNum(Higgs[i]) ,pNum(A), pNum(A), Higgs[i] ,A ,A);
   }
  
   if(VZdecay==0 || VWdecay==0) 
   { VZdecay=1; VWdecay=1; cleanDecayTable(); changeVirtual=1;}

   char  *t =pdg2name(6);
   if(t) slhaDecayPrint(t, 0, f);

   for(i=0;i< nHiggs;i++)
   {  txtList L,l;
      double  BrAA=0, BrGG=0;
      double width,widthP;
      width=pWidth(Higgs[i],&L);
      for(l=L;l;l=l->next)
      {  double br;
         int id[10],nd;
         nd=txt2plist(l->txt,&br,id); 
         if(nd==2)
         { if(id[0]==21 && id[1]==21) BrGG=br;
           if(id[0]==22 && id[1]==22) BrAA=br;
         }   
      }
//BrGG=0;  BrAA=0;
      if(BrGG>0 && BrAA>0) slhaDecayPrint(Higgs[i], 0, f); else
      { double wGG,wAA,brGG_=0,brAA_=0;
        double mH=pMass(Higgs[i]);
        double a=alphaQCD(mH)/M_PI;
        double Rqcd=(1+a*(149./12.+a*(68.6482-a*212.447))); 
        wGG=2*Rqcd*pow(mH,3)*lamQGG[i]/M_PI;
        wAA=0.25*pow(mH,3)*lamQAA[i]/M_PI;        
        double K=1;
        if(BrGG<=0 && wGG>0) K+=wGG/width;
        if(BrAA<=0 && wAA>0) K+=wAA/width; 
        fprintf(f,"DECAY %d %E # %s\n", pNum(Higgs[i]), K*width, Higgs[i]);
        if(BrGG<=0 && wGG>0) fprintf(f," %12.4E 2 %4d %4d # %s %s\n", wGG/(K*width),pNum(G),pNum(G),G,G);
        if(BrAA<=0 && wAA>0) fprintf(f," %12.4E 2 %4d %4d # %s %s\n", wAA/(K*width),pNum(A),pNum(A),A,A); 
         
        for(l=L;l;l=l->next)
        { double br;
          int id[10],nd;
          nd=txt2plist(l->txt,&br,id);
          fprintf(f," %e %d ", br/K, nd);
          for(int j=0;j<nd;j++) fprintf(f," %4d", id[j]);
          fprintf(f," # %s\n", strstr(l->txt,"->")+2);
        }                                                                      
      }
   }
   
   for(i=0;i<nHch;i++) slhaDecayPrint(chHiggs[i], 0, f);
      
   free(Higgs); free(chHiggs); free(lamQGG); free(lamQAA);
   fclose(f);   
   if(nHiggsCh) *nHiggsCh=nHch;
   if(changeVirtual)
   { VZdecay=VZdecay_save; VWdecay=VWdecay_save; cleanDecayTable();}

   return nHiggs;
}  

int  hbblocksmo_(char *fname, int * nHch,int len)
{ 
  char * cname=malloc(len+2);
  int nHiggs;  
  fName2c(fname,cname,len);    
  nHiggs= hbBlocksMO(cname,nHch);      
  free(cname);        
  return nHiggs;          
}
  