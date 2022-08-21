//  SECTION  MODEL GENERATION 
//  SECTION  SCALAR COEFFICIENTS
//  SECTION  LOOP CORRECTIONS 
//  SECTION  nucleonAmplitudes
//  SECTION  NUCLEI  
//  SECTIONS Nuclear Form Factors for vector current

#include"micromegas_aux.h"
#include"micromegas_f.h"
#define  NEWMODEL
#define  GG_NLO


int DDLflag=1;
int QCDcorrections=1;
int Twist2On=1;

extern double DDmomSQ;
 
#define INMAX 100

//  SECTION  MODEL GENERATION 

static int  create2_th_model(char * disp,char * pname,char *ProcNameSI, char* ProcNameSD)
{  
   char * command;
   char buff[100];
   int line,i,k,err;
   int neutr;
   char aP_P[20];
   char * auxVert[10];
   int nAuxVerts;  
   int K;
   
   struct pProp
   { int num;
     char name[10];
     char aname[10];
     char mass[10];
     int spin2;
     int color;
     int line;
   } pInvolved[INMAX];
       
   for(i=0;i<INMAX;i++) pInvolved[i].num=0;   
   K=1;
   for(line=0;line<nModelParticles;line++)
   { 
     char * name1, *name2; 
     int num; 
     int color;
     name1=ModelPrtcls[line].name;
     name2=ModelPrtcls[line].aname;
     num=ModelPrtcls[line].NPDG;
     color=ModelPrtcls[line].cdim;
   
      if((strcmp(pname,name1)==0||strcmp(pname,name2)==0||
      (abs(color)!=1 && abs(num)!=21)))
      { 
        if(color!=1) k=K++; else k=0;
         pInvolved[k].num=num;
         strcpy(pInvolved[k].name,name1);
         strcpy(pInvolved[k].aname,name2);
         strcpy(pInvolved[k].mass,ModelPrtcls[line].mass);
         pInvolved[k].spin2=ModelPrtcls[line].spin2;
         pInvolved[k].color=ModelPrtcls[line].cdim;
         pInvolved[k].line=line+1;
      }
   }

   if(pInvolved[0].num==0)
   { printf("Direct detection module can not work because the model contains a particle\n");
     printf("with zero PDG code %s\n",  pInvolved[i].name);
     return 10;     
   }
      
   for(i=1;i<K;i++)
   if(pInvolved[i].num>=-6 && pInvolved[i].num<=6 && strcmp(pInvolved[i].mass,"0")==0)
   {  printf("\n Error! Direct detection module can not work because the model contains\n");
      printf("zero mass quark  '%s'\n",  pInvolved[i].name); 
      return 1;
   }
#ifdef NEWMODEL  
   command=malloc(strlen(WORK)+strlen(calchepDir)+7000);
   sprintf(command,"cd %s;  %s/bin/s_calchep -blind "
   "\"{[[{[{\\19\\24{{MDD__\\091\\09"    //   new  parameter MDD__ for mass of auxilarry particles  

   "{sqrt6D\\092.449489742783178\\09"   //     sqrt(6) constant
   
   "}}0\"",disp,calchepDir);
   err=system( command);
   if(err) {free(command);return 3;}

   sprintf(command,"cd %s;  %s/bin/s_calchep -blind \"[[{[[{"
   "\\19%d{\\04\\19\\${{}}0\";",disp,calchepDir,pInvolved[0].line);
   err=system( command);
   if(err) {free(command);return 3;}
      
 
    sprintf(command,"cd %s;  %s/bin/s_calchep -blind \"[[{[[{"   
   "{_S0_\\09_S0_\\09_S0_\\090\\090\\09MDD__\\090\\091\\09*\\09_S0_\\09_S0_\\09"   // new point-like scalar _S0_
   "{_V5_\\09_V5_\\09_V5_\\090\\092\\09MDD__\\090\\091\\09*\\09_V5_\\09_V5_\\09"   // new point-like vector _V5_
      
   "}"
   
   "}0\";",disp,calchepDir);
   err=system( command);
   if(err) {free(command);return 3;}

   sprintf(command,"cd %s; %s/bin/s_calchep -blind \"[[{[[[{",disp,calchepDir);
   for(i=1;i<K;i++) 
   { sprintf(aP_P,"%s\\09%s",pInvolved[i].aname,pInvolved[i].name); 
 
     if(pInvolved[i].spin2==1)
     {      sprintf(command+strlen(command),     
       "{%s\\09_S0_\\09\\09MDD__\\091\\09"                  //   (F f _S0_)   MDD__ interaction for fermions     
       "{%s\\09_V5_\\09\\09-MDD__\\09G5*G(m3)\\09"          //   (F f _V5_)  -MDD__*G5*G(m3) interaction
       ,aP_P,aP_P);        
     }
     else  if(pInvolved[i].spin2==0)
     {
       sprintf(command+strlen(command),
       "{%s\\09_S0_\\09\\09MDD__\\091\\09",aP_P);          // (s S _S0_) MDD__   interaction for scalars  
     }
   }
   sprintf(command+strlen(command),"}}0\""); 
   err=system( command);  
   
   if(err) {free(command);return 4;}

   sprintf(command,"cd %s; %s/bin/s_calchep -blind \"[[{[[[{",disp,calchepDir);

   neutr=strcmp(pInvolved[0].name,pInvolved[0].aname)==0; 
   sprintf(aP_P,"%s\\09%s", pInvolved[0].aname, pInvolved[0].name); 
   if (pInvolved[0].spin2==0)                                               // for scalar DM:
   {  sprintf(buff,"_S0_\\09\\092*MDD__*2*%s\\091\\09",pInvolved[0].mass);  //  (DM  dm _S0_)*2*MDD__*Mdm  interaction
      nAuxVerts=1;
      auxVert[0]=buff;
   }else 
   if(pInvolved[0].spin2==1)                             // for spinor DM:
   {  
      nAuxVerts=2;                                              
      auxVert[0]="_S0_\\09\\092*MDD__\\091\\09";         // (DM  dm _S0_)*2*MDD__ 
      auxVert[1]="_V5_\\09\\092*MDD__\\09G5*G(m3)\\09";  // (DM  dm _V5_)*2*MDD__*G5*G(m3)
   }else 
   if (pInvolved[0].spin2==2)                            // for vector DM:    
   { 
      nAuxVerts=2;
      sprintf(buff,"_S0_\\09\\092*MDD__*2*%s\\09m1.m2\\09",pInvolved[0].mass); //(DM dm _S0_)*2*MDD__*Mdm*m1.m2
      auxVert[0]=buff;
      auxVert[1]="_V5_\\09\\092*i*MDD__*sqrt6D\\09eps(m1,m2,(p1-p2),m3)\\09";  //(DM dm _V5_)*2*i*MDD__*sqrt(6)*eps(m1,m2,(p1-p2),m3)
   }   
   for(i=0;i<nAuxVerts;i++) 
   sprintf(command+strlen(command),"{%s\\09%s",aP_P,auxVert[i]);
    
   strcat(command, "}}0\"");

   err=system(command);
   free(command); 
    
   if(err) return 5; 
#endif                 
   sprintf(ProcNameSI,"QUARKS,%s->QUARKS,%s{",pname, pname);
   strcpy(ProcNameSD,ProcNameSI);
   {  char *lquarks=malloc(maxPlistLen), *allcol=malloc(maxPlistLen);
      for(i=1,lquarks[0]=0, allcol[0]=0;i<K;i++)
      { int forlQ=(abs(pInvolved[i].num)<=3);
        sprintf(allcol+strlen(allcol),"%s,",pInvolved[i].name);
        if(forlQ)sprintf(lquarks+strlen(lquarks),"%s,",pInvolved[i].name);
        if(strcmp(pInvolved[i].name,pInvolved[i].aname))
        {
          sprintf(allcol+strlen(allcol),"%s,",pInvolved[i].aname);
          if(forlQ)sprintf(lquarks+strlen(lquarks),"%s,",pInvolved[i].aname);
        }        
      }
      
      if(strlen(lquarks))lquarks[strlen(lquarks)-1]=0;
      if(strlen(allcol))  allcol[strlen( allcol)-1]=0; 
      if( allcol[0]) strcat(ProcNameSI, allcol); else ProcNameSI[0]=0; 
      if(lquarks[0]) strcat(ProcNameSD,lquarks); else ProcNameSD[0]=0;
      if(pInvolved[0].spin2==0) ProcNameSD[0]=0;
      free(lquarks), free(allcol);
   }    
   return 0;
}     

static int getAuxCodesForDD(char*pname,numout**ccSI,numout**ccSD)
{ 
  char libnameSD[40]="dirSD";
  char libnameSI[40]="dirSI";
  int newDir=0;
  char ProcNameSI[4000], ProcNameSD[4000];
  int err=0;
  char exclude[30];
  pname2lib(pname,libnameSD+5);
  pname2lib(pname,libnameSI+5);

  if(ccSI)*ccSI=getMEcode(0,1,"",NULL,NULL,libnameSI);
  if(ccSD)*ccSD=getMEcode(0,1,"",NULL,NULL,libnameSD);
  if( ( ccSI&&(!*ccSI))  || ( ccSD&&(!*ccSD) ) ) 
  { int line;
    newDir=prepareWorkPlace();
    err=create2_th_model(compDir,pname,ProcNameSI,ProcNameSD);
    if(err) return err;
  }  
  
   
  if(ccSI && !*ccSI && ProcNameSI[0])
  {  sprintf(exclude,"_S0_!=1,_V5_"); 
     *ccSI=getMEcode(0,1,ProcNameSI,exclude,NULL,libnameSI);
  }
  if(ccSD && !*ccSD && ProcNameSD[0])
  {  sprintf(exclude,"_V5_!=1,_S0_"); 
     *ccSD=getMEcode(0,1,ProcNameSD,exclude,NULL,libnameSD);
  }      
  
  if(ccSI && *ccSI) *((*ccSI)->interface->BWrange)=0;
  if(ccSD && *ccSD) *((*ccSD)->interface->BWrange)=0;
      
  if(newDir) cleanWorkPlace();
  return 0;
}

// SECTION SCALAR COEFFICIENTS

void calcScalarFF(double muDmd, double msDmd, double sigmaPiN, double sigma0)
{ 
  printf("This function is obsolete. Use calcScalarQuarkFF instead where\n"
         "the last argument is the sigma_s parameter\n");  

  double Mp=0.9383,Mn=0.9396;
  double y=1-sigma0/sigmaPiN, z=1.49,alpha=(2*z-(z-1)*y)/(2+(z-1)*y);
  ScalarFFPd=2*sigmaPiN/((1+ muDmd)*Mp*(1+alpha))/1000.;
  ScalarFFPu=muDmd*alpha*ScalarFFPd;
  ScalarFFPs=sigmaPiN*y*msDmd/(1+ muDmd)/Mp/1000.;

  ScalarFFNd=2*sigmaPiN*alpha/((1+ muDmd)*Mn*(1+alpha))/1000.;
  ScalarFFNu=muDmd*ScalarFFNd/alpha;
  ScalarFFNs=sigmaPiN*y*msDmd/(1+ muDmd)/Mn/1000.;
}


void calcScalarQuarkFF(double muDmd, double msDmd, double sigmaPiN, double sigmaS)
{ double Mp=0.9383,Mn=0.9396;
  double sigma0=sigmaPiN*(1-sigmaS/sigmaPiN*(1+muDmd)/msDmd);
  double y=1-sigma0/sigmaPiN, z=1.49,alpha=(2*z-(z-1)*y)/(2+(z-1)*y);
  ScalarFFPd=2*sigmaPiN/((1+ muDmd)*Mp*(1+alpha))/1000.;
  ScalarFFPu=muDmd*alpha*ScalarFFPd;
  ScalarFFPs=sigmaPiN*y*msDmd/(1+ muDmd)/Mp/1000.;

  ScalarFFNd=2*sigmaPiN*alpha/((1+ muDmd)*Mn*(1+alpha))/1000.;
  ScalarFFNu=muDmd*ScalarFFNd/alpha;
  ScalarFFNs=sigmaPiN*y*msDmd/(1+ muDmd)/Mn/1000.;   
}

//  SECTION  LOOP CORRECTIONS 

extern double (*loopFF__)(double,double,double,double);

static double zeroloopFactor(double sgn, double mq,double msq,double mne)
{ 
   if(msq>0.5)return 0;
  return   1/(msq*msq-mne*mne-mq*mq -2*sgn*mne*mq);
}

static double Delta(double M, double m1, double m2)
{ 
  double ms=m1*m1+m2*m2, md=m1*m1-m2*m2;
  return M*M*(M*M - 2*ms)+md*md;
}

static double Lambda(double M, double m1, double m2)
{
  double D=Delta(M,m1,m2);
  double Mm=m1*m1+m2*m2-M*M;
  double sD=sqrt(fabs(D));
//  if(sD<1E-5*Mm) return  2/Mm;
//  else 
  if(D>=0) return  log((Mm+sD)/(Mm-sD))/sD;
  else          return 2*atan( sD/Mm)/sD;  
}

// Fermionic Dark Matter  1502.02244   DreesNojiri
/*  Loop integrals I1 ... I5  */

#define SQ(x)  ((x)*(x))

double   LintIk(int II,double MSQ,double MQ,double MNE)
{  
  double LAM,SPPM,SPMM,R1,R2,R3,del,CMD;
  double msq2=MSQ*MSQ, mq2=MQ*MQ, mne2=MNE*MNE;

  SPPM=  msq2+mq2-mne2, SPMM= msq2-mq2-mne2;
  R1  =(msq2-mq2)/mne2;
  R2  =(mq2-mne2)/msq2;
  R3  =(msq2-mne2)/mq2;

  del =2.*mne2*(mq2+msq2)-mne2*mne2-SQ(msq2-mq2);
  
  if(del>0) LAM=2.*atan(sqrt(del)/SPPM)/sqrt(del);
  if(del<0) LAM=log((SPPM+sqrt(-del))/(SPPM-sqrt(-del)))/sqrt(-del);

  switch(II)
  { 
    case 1:  
      CMD=1./del*(R2/3-2/3.*R3-5/3.+ (2*msq2-2/3.*mne2)*LAM);	
    break; 
    case 2: 
      CMD=(log(msq2/mq2)-SPMM*LAM)/2./mne2/mne2+
         ( ((mq2*mq2-mq2*msq2)/mne2-7/3.*mq2+2/3.*(mne2-msq2))*LAM+R2/3+R1+2/3.
         )/del;
    break;
    case 3:
      CMD=-3/SQ(del)*SPPM+LAM/del*(-1+6*mq2*msq2/del);
    break;	
    case 4:
      CMD=((log(msq2/mq2) - SPMM*LAM)/2/mne2-1/msq2 -mq2*SPMM/del*LAM)/mne2/mne2
      
         +( mq2/mne2/mne2-SQ(1-mq2/mne2)/msq2+0.5/mne2
            +3*mq2/del*(1 +  R1 + (-R1*mq2-2*mq2-msq2+mne2)*LAM)
          )/del;
    break;
    case 5:
     CMD=(log(msq2/mq2)-SPMM*LAM)/(2*mne2*mne2)-(LAM*(2*(msq2-mne2)+3*mq2+R1*mq2)-3-R1)/del;
     break;
    default: CMD=0.; 
  }
  return CMD;
}

static double FeScLoop(double sgn, double mq,double msq,double mne)
{ 
/*  return   1/(msq*msq-mne*mne); */
return  1.5*mq*( 
       mq*(LintIk(1,msq,mq,mne)
  -2./3.*mne*mne*LintIk(3,msq,mq,mne))

 -sgn*mne*(LintIk(2,msq,mq,mne)-1./3.*LintIk(5,msq,mq,mne)-2./3.*mne*mne*LintIk(4,msq,mq,mne))
                ); 
}

//  Vectors  Dark Matter  1502.02244

static double FpVectA(double M, double m1, double m2)
{
  double D= Delta(M,m1,m2),L=Lambda(M,m1,m2), Mq=M*M, m1q=m1*m1, m2q=m2*m2;

  return  (D*(Mq*(m2q-m1q) +m1q*(m1q+5*m2q)) -6*m1q*m2q*( (m2q-m1q)*(m2q-m1q) -Mq*(m1q+3*m2q)))/(6*D*D*Mq) 
         -log(m1q/m2q)*m1q/(12*Mq*Mq)
         +(D*D*(m2q+m1q-Mq)+2*m2q*D*(5*m2q*m2q+20*m1q*m2q-m1q*m1q+Mq*(9*m2q+m1q))
         +12*m2q*m2q*(Mq*(m2q*m2q+10*m1q*m2q +5*m1q*m1q)-(m2q-m1q)*(m2q-m1q)*(m2q+3*m1q)))*m1q*L/(12*D*D*Mq*Mq);
}

static double FmVectA(double M, double m1, double m2)
{
   double D= Delta(M,m1,m2), L=Lambda(M,m1,m2), Mq=M*M, m1q=m1*m1, m2q=m2*m2;
   return -m2*(D*(2*m2q+m1q-2*Mq) +6*m1q*m2q*(m2q-m1q-Mq))/(6*m1*D*D) + m1*m2q*m2*(D+m1q*(m2q-m1q+Mq))*L/(D*D);
}
         
static double FpVectC(double M, double m1, double m2)
{
   double D= Delta(M,m1,m2), L=Lambda(M,m1,m2), Mq=M*M, m1q=m1*m1, m2q=m2*m2;
   return    -(2*Mq*Mq - 3*(m1q+m2q)*Mq+(m1q-m2q)*(m1q-m2q))/(6*D*Mq) +log(m1q/m2q)*(m1q-m2q)/(12*Mq*Mq)
             +L*(D*(Mq-m1q-m2q)*(m1q+m2q)+4*m1q*m2q*((m1q-m2q)*(m1q-m2q)-2*Mq*(m1q+m2q)))/(12*D*Mq*Mq); 
}  

static double FpVect(double M, double m1, double m2){ return  FpVectA(M, m1, m2)+ FpVectA(M, m2, m1)+FpVectC(M,m1,m2);}

static double FmVect(double M, double m1, double m2){ return  FmVectA(M, m1, m2)+ FmVectA(M, m2, m1);}

static double VectFermLoop(double sign, double m1,double m2, double M)
{ 
   return -3*m1/m2*(FmVect(M,m1,m2)-sign*(FpVect(m1,m2,M)*m2+FmVect(M,m1,m2)*m1)/M);
} 

// Scalar  Dark Matter
static double FpScalA(double M, double m1, double m2)
{
  double D= Delta(M,m1,m2),L=Lambda(M,m1,m2), Mq=M*M, m1q=m1*m1, m2q=m2*m2;

  return -L*m1q*m2q*m2q*(Mq+m1q-m2q)/(D*D) -((-Mq+m1q+2*m2q)*D +6*m1q*m2q*(Mq-m1q+m2q) )/(6*D*D);
}  
  
static double FmScalA(double M, double m1, double m2)
{
  double D= Delta(M,m1,m2),L=Lambda(M,m1,m2), Mq=M*M, m1q=m1*m1, m2q=m2*m2;
 
    return L*m1*m2q*m2*(D +m1q*(Mq-m1q+m2q))/(D*D) - m2*((-2*Mq+m1q+2*m2q)*D-6*m1q*m2q*(Mq+m1q-m2q))/(6*m1*D*D);  
}

static double FpScalC(double M, double m1, double m2)
{
  double D=Delta(M,m1,m2),L=Lambda(M,m1,m2), Mq=M*M, m1q=m1*m1, m2q=m2*m2; 
  return  (-Mq+m1q+m2q)/(2*D)-m1q*m2q*L/D;
}

static double FmScalC(double M, double m1, double m2)
{
  double D=Delta(M,m1,m2),L=Lambda(M,m1,m2), Mq=M*M, m1q=m1*m1, m2q=m2*m2; 
  return  2*m1*m2/D -m1*m2*(-Mq+m1q+m2q)*L/D;
}

static double FpScal(double M, double m1, double m2) { return FpScalA(M,m1,m2) + FpScalA(M,m2,m1) + FpScalC(M,m1,m2); }        
static double FmScal(double M, double m1, double m2) { return FmScalA(M,m1,m2) + FmScalA(M,m2,m1) + FmScalC(M,m1,m2); } 

static double ScalFermLoop(double sign, double m1,double m2,double M)
{ 
   return -3*m1/m2*(FmScal(M,m1,m2)+sign*(FpScal(M,m1,m2)*m2-FmScal(M,m1,m2)*m1)/M);
} 

static double pdfQnum,p1_p2;

static double twist2FF__(double sgn, double mq,double msq,double mne)
{ 
  double D=(msq*msq-mne*mne -mq*mq);
  double D2=D*D-4*mq*mq*mne*mne;
  
//return  parton_x(pdfQnum,1.2*(msq-mne*mne/msq))/D2*(D +sgn*2*p1_p2);
  return  parton_x(pdfQnum,msq-mne              )/D2*(D +sgn*2*p1_p2);
}

// Twist-2 subtraction

static double twist2_subtractionFF__(double sgn, double mq,double msq,double mne)
{ double D=(msq*msq-mne*mne -mq*mq);
  double D2=D*D-4*mq*mq*mne*mne;
  return   1/D2*(D +sgn*2*p1_p2);
}

// SECTION   nucleonAmplitudes

int nucleonAmplitudes(char * WIMP, double*pA0,double*pA5,double*nA0,double*nA5) 
{
  double wimpMass; 
  int wimpN,Qnum,aQnum,i,II,sgn,ntot,n;
  double  s,MN=0.939; 
  numout *ccSI,*ccSD,*cc; 
  double wGluP, wGluN,GG;
  REAL pvect[16]; 
  
  double wS0P__[6],wS0N__[6]; /*scalar */
  double wV5P__[3],wV5N__[3]; /*pseudo-vector*/
  double wSM0P[3],wSM0N[3]; /* sigma */
  double  wV0P[3]={1,2,0}, wV0N[3]={2,1,0}; /* vector current */

   for(i=0;i<3;i++) 
   { wS0P__[i]= *(&(ScalarFFPd)+i);
     wS0N__[i]= *(&(ScalarFFNd)+i);
     wV5P__[i]= *(&(pVectorFFPd)+i);    
     wV5N__[i]= *(&(pVectorFFNd)+i); 
     wSM0P[i]=  *(&(SigmaFFPd)+i);   
     wSM0N[i]=  *(&(SigmaFFNd)+i);
  }

  for(s=0,i=0;i<3;i++) s+= wS0P__[i];
  for(s=2./27.*(1-s),i=3;i<6;i++)wS0P__[i]=s;

  for(s=0,i=0;i<3;i++) s+= wS0N__[i];
  for(s=2./27.*(1-s),i=3;i<6;i++)wS0N__[i]=s;
                                   
  wGluP=wS0P__[4]*27./2.; wGluN=wS0N__[4]*27./2.;

  for(i=0;i<2;i++) {pA0[i]=0; pA5[i]=0; nA0[i]=0; nA5[i]=0;}
  
  wimpMass= pMass(WIMP);
  GG=sqrt(4*M_PI*alphaQCD(wimpMass));
  
  for(wimpN=0; wimpN<Nodd;wimpN++) if(strcmp(OddPrtcls[wimpN].name,WIMP)==0 ||
                           strcmp(OddPrtcls[wimpN].aname,WIMP)==0) break;

  if(OddPrtcls[wimpN].spin2) getAuxCodesForDD(WIMP,&ccSI,&ccSD);
  else { getAuxCodesForDD(WIMP,&ccSI,NULL); ccSD=NULL;}

  if(!ccSI && !ccSD)  return -1;
  
  for(II=0;II<2;II++)
  { if(II) cc=ccSD; else cc=ccSI;
    if(!cc) continue;
    procInfo1(cc,&ntot,NULL,NULL);
    for(i=1;i<=cc->interface->nvar;i++) 
    if(cc->link[i]) cc->interface->va[i]=*(cc->link[i]);
    else { printf("absence of link:%s\n", cc->interface->varName[i]);}
    
    if(cc->Q) *(cc->Q)=2*wimpMass;

    cc->interface->calcFunc();
/*    
    *(cc->interface->gtwidth)=0;
    *(cc->interface->twidth)=1;
    *(cc->interface->gswidth)=0;
*/
    for(n=1;n<=ntot;n++)
    { double cs0, ampl, qMass;
      char * names[4];
      REAL masses[4];
      int pdg[4];
      int l,err_code=0;
      int spin2, spin2Wimp, color,neutral,neutralWIMP;
      double alphaMq, qcdNLO;
      int lf;
      
      for(l=0;l<4;l++) names[l]= cc->interface->pinf(n,l+1,masses+l,pdg+l); 
      if(strcmp(names[0],names[2]) ) continue; 

      cc->interface->pinfAux(n,2,&spin2Wimp,NULL,&neutralWIMP,NULL);       
      cc->interface->pinfAux(n,1,&spin2,&color,&neutral,NULL);    
      Qnum=pdg[0];
      qMass=masses[0];
      if(!qMass) continue;
      sgn=Qnum>0?1:-1;
      aQnum=abs(Qnum);

      if(QCDcorrections && aQnum>3 )
      {  switch(aQnum)
         { case 4: alphaMq=0.39; break; 
           case 5: alphaMq=0.22; break; 
           default:alphaMq=alphaQCD(qMass);
         }  
      } else alphaMq=0;

      loopFF__=NULL;
      
      if(II==0 && DDLflag && aQnum>3)
      {  
         if(names[0][0]!='~')
         { 
           if(spin2==1) // even  color  fermion  particle
           switch(spin2Wimp)
           {  case 0:    loopFF__= ScalFermLoop; break;
              case 1:    loopFF__= FeScLoop;     break;
              case 2:    loopFF__= VectFermLoop; break;
              default:   loopFF__= NULL;         break;
           } else loopFF__=NULL;
         }  
         else
         { if((spin2==0 && spin2Wimp==1 ) || (spin2==1 && (spin2Wimp==0||spin2Wimp==2 )) ) loopFF__= zeroloopFactor;
           else loopFF__=NULL; 
         }             
      } else  loopFF__=NULL;
      
      for(i=0;i<16;i++) pvect[i]=0;     
      for(i=0;i<4;i++) pvect[4*i]=masses[i];
      err_code=0;
      cs0= (*cc->interface->sqme)(n,GG,pvect,NULL,&err_code);
      if(II==0 && loopFF__==NULL && spin2==1) //twist-2 trace  subtraction from SD amplitude
      {
          int k;
          double diff;
          double h=masses[0]/100;
          double d4[3]={0.5,-1,0.5};
          h=0.001; 
          loopFF__=twist2_subtractionFF__; 
          for(diff=0,k=0;k<3;k++)
          {   double E=masses[0]+h*k;
              double P= (k==0)?0: sqrt(E*E-masses[0]*masses[0]);
              pvect[0]=pvect[8]=E;
              pvect[3]=pvect[11]=P;
              p1_p2=E*masses[1];
              diff+=(*cc->interface->sqme)(n,GG,pvect,NULL,&err_code)*d4[k];                
          }
          diff/=pow(h*masses[1],2);            
          cs0-=3./4.*diff*masses[1]*masses[1]*masses[0]*masses[0];
      }
      loopFF__=NULL;
       
      if(spin2==1)      ampl= cs0/(128*masses[0]*masses[0]*masses[1]*masses[1]);
      else if(spin2==0) ampl= cs0/(32*masses[1]*masses[1]);
      else ampl=0;
/* Comment: 128=2*64.  2-because of we have calculated interference term. */       
/* if(ampl) printf("II=%d %E for %s %s -> %s %s\n",II,ampl, names[0],names[1],names[2],names[3]); */     
      if(II)  /* Spin dependent case */
      { if(aQnum<=3 ) 
        {  /* here (ampl/3) because of normalization of SD amplitude */ 
         pA5[0]+=(ampl/3)*wV5P__[aQnum-1]/2;
         nA5[0]+=(ampl/3)*wV5N__[aQnum-1]/2; 
         pA5[1]+=sgn*(ampl/3)*wSM0P[aQnum-1]/2;
         nA5[1]+=sgn*(ampl/3)*wSM0N[aQnum-1]/2;                                     
        }
      }
      else /* Spin independent case */
      { 
          /* scalar SI amplitude */
        if(aQnum>3)
        { int color,spin2,neutral;
          cc->interface->pinfAux(n,1,&spin2,&color,&neutral,NULL);
        
          switch(color)
          { case  3:
            case -3:
              switch(spin2)
              { case 0:
                  qcdNLO=1+(25./6.-16./9.) *alphaMq/M_PI;            
                  pA0[0]+=ampl*qcdNLO*wGluP*(1./108.)*MN/qMass/qMass/2;
                  nA0[0]+=ampl*qcdNLO*wGluN*(1./108.)*MN/qMass/qMass/2;
                break;
                case 1:
                  qcdNLO=1+(11./4. -16./9.)*alphaMq/M_PI;
                  pA0[0]+=ampl*qcdNLO*wGluP*(2./27.)*MN/qMass/2;
                  nA0[0]+=ampl*qcdNLO*wGluN*(2./27.)*MN/qMass/2;
                break;           
                default:
                printf("2*spin=%d, color=3 - not implemented \n",spin2);
              }
              break;
            case  8:
              switch(spin2)
              { case 0:
                  pA0[0]+=ampl*wGluP*(6./108.)*MN/qMass/qMass/2;
                  nA0[0]+=ampl*wGluN*(6./108.)*MN/qMass/qMass/2;
                break;
                case 1:

                  pA0[0]+=ampl*wGluP*(2.*6./27.)*MN/qMass/2 ;
                  nA0[0]+=ampl*wGluN*(2.*6./27.)*MN/qMass/2;
                break;           
                default:
                printf("2*spin=%d, color=8 - not implemented \n",spin2);
              }
              break;  
           }    
        }else
        {
            pA0[0]+=ampl*wS0P__[aQnum-1]*MN/qMass/2 ; 
            nA0[0]+=ampl*wS0N__[aQnum-1]*MN/qMass/2;
            pA0[1]+=sgn*ampl*wV0P[aQnum-1]/2;
            nA0[1]+=sgn*ampl*wV0N[aQnum-1]/2;
        }
                        
        if(aQnum<6 && Twist2On) /* twist-2  SI amplitude  */   
        { int k;
          double g=0,gp,gn;
          double h=masses[0]/100;
          double d4[3]={0.5,-1,0.5};
          h=0.001;
                     
          gp=0;
          gn=0;          
          loopFF__=twist2FF__;  /* twist-2 contribution  */
          for(k=0;k<3;k++)
          { double E=masses[0]+h*k;
            double P= (k==0)?0: sqrt(E*E-masses[0]*masses[0]);
            pvect[0]=pvect[8]=E;
            pvect[3]=pvect[11]=P;
            p1_p2=E*masses[1];
            pdfQnum=aQnum;
            gp+= (*cc->interface->sqme)(n,GG,pvect,NULL,&err_code)*d4[k];
            if(aQnum==1) 
            { pdfQnum=2;
              gn+= (*cc->interface->sqme)(n,GG,pvect,NULL,&err_code)*d4[k];
            }else if(pdfQnum==2)
            {
              pdfQnum=1;
              gn+= (*cc->interface->sqme)(n,GG,pvect,NULL,&err_code)*d4[k];
            }
          } 
          if(aQnum>2) gn=gp;         
          gp/=pow(h*masses[1],2)*128*masses[0]*masses[1]  * 2 ;
          gn/=pow(h*masses[1],2)*128*masses[0]*masses[1]  * 2 ;
 
          pA0[0]+=1.5*gp*masses[1]*MN/2;
          nA0[0]+=1.5*gn*masses[1]*MN/2;   
          loopFF__=NULL;
        }
      }              
    }    
  }
/*  OnlyTEQ0=0; */
  { double mem;  
    mem=pA0[1]; pA0[1]=pA0[0]-mem; pA0[0]+=mem;
    mem=pA5[1]; pA5[1]=pA5[0]-mem; pA5[0]+=mem;

    mem=nA0[1]; nA0[1]=nA0[0]-mem; nA0[0]+=mem;
    mem=nA5[1]; nA5[1]=nA5[0]-mem; nA5[0]+=mem;     
  }   
  return 0; 
}

// SECTION  NUCLEI  

double maxRecoil(double A) 
{ 
  double Mdm=0;
  if(CDM1!=NULL) Mdm=Mcdm1;
  if(CDM2!=NULL && Mcdm2>Mcdm1) Mdm=Mcdm2;
  if(Mdm==0) Mdm=Mcdm;        
  return 1E6*2*A*0.94* pow( Mdm*(vEsc+vEarth)/299792./(A*0.94+Mdm),2);
}


/*===== Intergrands for  Fermi nucleus density =====*/

static double R_, p_, a_, C_=1.23, B_=-0.6, A_=0.52;

static double FermiND0(double r){ return r*r/(1+exp((r-R_)/a_));}
static double FermiNDP(double r){ return sin(p_*r)*r/(1+exp((r-R_)/a_));}

double FermiFF(int A, double Qfermi)
{
  A_=Fermi_a;
  B_=Fermi_b;
  C_=Fermi_c;
  R_=C_*pow(A,1./3.)+B_;
  a_=A_;
  p_=Qfermi;
  return simpson(FermiNDP,0., R_+5, 1.E-4,NULL)/simpson(FermiND0,0., R_+5, 1.E-4,NULL)/p_;
}


/*==== nucleusRecoil: main functions ================*/

static double (*_fv_)(double v);
static double vfv(double v){ return _fv_(v)/v;}

static double nucleusRecoil_stat(double M_cdm, double(*fv)(double),
      int A, int Z, double J,
      
      void (*Sxx)(double, double *,double *, double *),
//      double(*S00)(double),double(*S01)(double),double(*S11)(double),
      
      double css,double csv00, double csv01,double csv11,
      double * dNdE)
{
  const double  Kg=1./1.782662E-27;          /* GeV  */
  const double  vC=299792;                   /* km/s */
  const double lDay=60.*60.*24.;             /* sec  */  

  const double step=1.E-6*RE_STEP;              /* 1KeV */
  _fv_=fv;
  int i;
  double MA,ffs,Rs,sum; 
  double FFs=1, s00,s01,s11;
  double E0,vmin,vmax;

  if(fv==Maxwell) vmax=vEsc+vEarth; 
  else vmax=1200.;
  
  MA=0.94*A;

  { 
    double Mr,SCcoeff;

    E0=M_cdm/(M_cdm+MA),E0=2*MA*E0*E0;
    Mr=MA*M_cdm/(MA+M_cdm);

    SCcoeff=4/M_PI*3.8937966E-28*Mr*Mr;

    css*=SCcoeff;
    csv00*=4*M_PI*SCcoeff/(2*J+1);    
    csv11*=4*M_PI*SCcoeff/(2*J+1);
    csv01*=4*M_PI*SCcoeff/(2*J+1);
  }
  
/*printf("css=%E, csv00=%E, csv11=%E, csv01=%E\n",css,csv00,csv11,csv01);*/

  if(A>1)
  {
    Rs=C_*pow(A,1./3.)+B_;
    R_=Rs; a_=A_; 
    ffs=simpson(FermiND0,0., R_+10*0.5, 1.E-4,NULL);
  }
  
  for(i=0,sum=0;i<RE_DIM;i++)
  { double E=RE_START*1E-6*pow(RE_STEP,i);

    vmin=sqrt(E/E0)*vC;
    double p=sqrt(E*MA*2)/0.197327;
    
    if(vmin>=vmax)  { dNdE[i]=0;continue;}
    if(A>1)
    { 
      p_=p;    
      R_=Rs; a_=A_; 
      FFs=simpson(FermiNDP,0., R_+10*0.5, 1.E-4,NULL)/p_/ffs;      
    } else FFs=1;

    dNdE[i]=FFs*FFs*css;

    if(J && Sxx)
    {  Sxx(p,&s00,&s01,&s11);
       dNdE[i]+=(csv00*s00+csv01*s01+csv11*s11);
    }   
    
    
    dNdE[i]*=simpson(vfv,vmin,vmax,1.E-4,NULL);
    dNdE[i]*=1/E0*1.E5*vC*vC*(rhoDM/M_cdm)*lDay*(Kg/MA)*1.E-6;
    if(i==0)            sum+=1E6*E*sqrt(RE_STEP)*dNdE[i];
    else if(i==RE_DIM-1) sum+=1E6*E/sqrt(RE_STEP)*dNdE[i];
    else                sum+=1E6*E*(sqrt(RE_STEP)+1/sqrt(RE_STEP))*dNdE[i];
  }
  return sum;
}

static double S00_0,S01_0,S11_0, RS_;

static void Sxx_(double p, double*s00,double*s01,double*s11)
{
  double  r=exp(-p*p*RS_*RS_/4);
  *s00=S00_0*r; *s01=S01_0*r; *s11=S11_0*r;
}

static double nucleusRecoil0_stat( double(*fv)(double),
int A, int Z, double J,double Sp,double Sn,
double css, double csv00, double csv01, double csv11,
double * dNdE)
{
  if(J)
  { double A3=pow(A,1./3.);
    RS_= 1.7*A3 -0.28 - 0.78*(A3-3.8 + sqrt((A3-3.8)*(A3-3.8) +0.2) );
      
    S00_0=  (Sp+Sn)*(Sp+Sn)*(2*J+1)*(J+1)/(4*M_PI*J);
    S11_0=  (Sp-Sn)*(Sp-Sn)*(2*J+1)*(J+1)/(4*M_PI*J);
    S01_0=2*(Sp+Sn)*(Sp-Sn)*(2*J+1)*(J+1)/(4*M_PI*J);
  }  
  return nucleusRecoil_stat(Mcdm,fv,A,Z,J,Sxx_, css,csv00,csv01,csv11,dNdE);
}

 
static double nucleusRecoil1(char *WINP, double(*fv)(double),
      int A, int Z, double J,
      void (*Sxx)(double,double*,double*,double*),
      double * dNdE)
{
  double css,csv00,csv01,csv11,M_cdm,MA,E; 
  double pA0[2],pA5[2],nA0[2],nA5[2],pA0_[2],pA5_[2],nA0_[2],nA5_[2];
  int i;
  _fv_=fv;

  M_cdm=pMass(WINP);
  MA=A*0.94;
  DDmomSQ=0; 
  nucleonAmplitudes(WINP,pA0,pA5,nA0,nA5);

  int smallMass=0;
  for(int i=0;i<2;i++) if(!isfinite(pA0[i]) || !isfinite(pA5[i]) || !isfinite(nA0[i]) || !isfinite(nA5[i])) smallMass=1;
  
  if(! smallMass)
  {
     E=RE_START*1E-6*pow(RE_STEP,RE_DIM-1);
     DDmomSQ=2*MA*E;
     nucleonAmplitudes(WINP,pA0_,pA5_,nA0_,nA5_);
     DDmomSQ=0;
     double eps=0.01;
     for(int i=0;i<2;i++) if( 
        fabs(pA0_[i]-pA0[i])>eps*(fabs(pA0_[i])+fabs(pA0[i]))||
        fabs(pA5_[i]-pA5[i])>eps*(fabs(pA5_[i])+fabs(pA5[i]))||
        fabs(nA0_[i]-nA0[i])>eps*(fabs(nA0_[i])+fabs(nA0[i]))||
        fabs(nA5_[i]-nA5[i])>eps*(fabs(nA5_[i])+fabs(nA5[i])) 
                            ) smallMass=1;
  }
  
  if(!smallMass)    
  {
     for(i=0,css=0,csv00=0,csv01=0,csv11=0;i<2;i++)
     { double AS=Z*pA0[i]+(A-Z)*nA0[i];
       double AVplus =pA5[i]+nA5[i];
       double AVminus=pA5[i]-nA5[i];
       double C=(1+dmAsymm*(1-2*i))/2;
       css+=AS*AS*C;
       csv00+=AVplus*AVplus*C;
       csv11+=AVminus*AVminus*C;
       csv01+=AVplus*AVminus*C;
     }  
     return nucleusRecoil_stat(M_cdm,fv,A,Z,J,Sxx,
                                css,csv00,csv01,csv11,dNdE);
  }
//  case of small Mass in t-channel


  const double  Kg=1./1.782662E-27;          /* GeV  */
  const double  vC=299792;                   /* km/s */
  const double lDay=60.*60.*24.;             /* sec  */  

  const double step=1.E-6*RE_STEP;              /* 1KeV */

  double ffs,Rs,sum; 
  double FFs=1, s00,s01,s11;
  double E0,vmin,vmax;
  double Mr, SCcoeff;
  
  if(fv==Maxwell) vmax=vEsc+vEarth; 
  else vmax=1200.;
  
  MA=0.94*A;

  
/*printf("css=%E, csv00=%E, csv11=%E, csv01=%E\n",css,csv00,csv11,csv01);*/

  if(A>1)
  {
    Rs=C_*pow(A,1./3.)+B_;
    R_=Rs; a_=A_; 
    ffs=simpson(FermiND0,0., R_+10*0.5, 1.E-4,NULL);
  }

  E0=M_cdm/(M_cdm+MA),E0=2*MA*E0*E0;
  Mr=MA*M_cdm/(MA+M_cdm);

  SCcoeff=4/M_PI*3.8937966E-28*Mr*Mr;
  
  for(i=0,sum=0;i<RE_DIM;i++)
  { double E=RE_START*1E-6*pow(RE_STEP,i);

    vmin=sqrt(E/E0)*vC;
    double p=sqrt(E*MA*2)/0.197327;
    
    if(vmin>=vmax)  { dNdE[i]=0;continue;}
    if(i && A>1)
    { 
      p_=p;    
      R_=Rs; a_=A_; 
      FFs=simpson(FermiNDP,0., R_+10*0.5, 1.E-4,NULL)/p_/ffs;      
    } else FFs=1;
    
    DDmomSQ=2*MA*E;
    nucleonAmplitudes(WINP,pA0,pA5,nA0,nA5);
    css=0,csv00=0,csv01=0,csv11=0;
    for(int k=0;k<2;k++)
    {  double AS=Z*pA0[k]+(A-Z)*nA0[k];
       double AVplus =pA5[k]+nA5[k];
       double AVminus=pA5[k]-nA5[k];
       double C=(1+dmAsymm*(1-2*k))/2;
       css+=AS*AS*C;
       csv00+=AVplus*AVplus*C;
       csv11+=AVminus*AVminus*C;
       csv01+=AVplus*AVminus*C;
    }

    css*=SCcoeff;
    csv00*=4*M_PI*SCcoeff/(2*J+1);    
    csv11*=4*M_PI*SCcoeff/(2*J+1);
    csv01*=4*M_PI*SCcoeff/(2*J+1);

    dNdE[i]=FFs*FFs*css;

    if(J && Sxx)
    {  Sxx(p,&s00,&s01,&s11);
       dNdE[i]+=(csv00*s00+csv01*s01+csv11*s11);
    }   
        
    dNdE[i]*=simpson(vfv,vmin,vmax,1.E-4,NULL);
    dNdE[i]*=1/E0*1.E5*vC*vC*(rhoDM/M_cdm)*lDay*(Kg/MA)*1.E-6;
    if(i==0)            sum+=1E6*E*sqrt(RE_STEP)*dNdE[i];
    else if(i==RE_DIM-1) sum+=1E6*E/sqrt(RE_STEP)*dNdE[i];
    else                sum+=1E6*E*(sqrt(RE_STEP)+1/sqrt(RE_STEP))*dNdE[i];
  }
   DDmomSQ=0;
  return sum;
}
 

double nucleusRecoil(double(*fv)(double),
int A, int Z, double J,
void (*Sxx)(double,double*,double*,double*),
double * dNdE)  
{  int i; 
   if(CDM1  && !CDM2)   return nucleusRecoil1(CDM1,fv,A,Z,J,Sxx,dNdE); 
   if(!CDM1 &&  CDM2)   return nucleusRecoil1(CDM2,fv,A,Z,J,Sxx,dNdE); 
   if(CDM1  &&  CDM2) 
   { double dNdE1[RE_DIM];
     double r1=0,r2=0;
     if(fracCDM2!=1) r1= nucleusRecoil1(CDM1,fv,A,Z,J,Sxx,dNdE1); 
     if(fracCDM2!=0) r2= nucleusRecoil1(CDM2,fv,A,Z,J,Sxx,dNdE);
     if(fracCDM2==1) {   return r2;}
     if(fracCDM2==0) { for(i=0;i<RE_DIM;i++) dNdE[i]=dNdE1[i];    return r1;}
      
     for(i=0;i<RE_DIM;i++) dNdE[i]=(1-fracCDM2)*dNdE1[i]+ fracCDM2*dNdE[i];
     return r1*(fracCDM2-1)+r2*fracCDM2;
   }
}

double nucleusRecoilCS( 
      double(*fv)(double),
      int A, int Z, double J,
      void(*Sxx)(double,double*,double*,double*),
      double  siP, double siN, double sdP,  double sdN,
      double * dNdE)
{

  double AS,AVplus, AVminus, css,csv00,csv01,csv11; 
  
  double LmbdP, XiP,LmbdN, XiN;
  double MN=0.939;
  double Mr=MN*Mcdm/(MN+Mcdm);
  double sCoeff= 2*sqrt(3.8937966E8/M_PI)*Mr;
   
  LmbdP=sqrt(fabs(siP))/sCoeff; if(siP<0) LmbdP*=-1;
  XiP  =sqrt(fabs(sdP)/3)/sCoeff; if(sdP<0) XiP*=-1;
  LmbdN=sqrt(fabs(siN))/sCoeff; if(siN<0) LmbdN*=-1;
  XiN  =sqrt(fabs(sdN)/3)/sCoeff; if(sdN<0) XiN*=-1;
   
  AS=Z*LmbdP+(A-Z)*LmbdN;
  AVplus =XiP+XiN;
  AVminus=XiP-XiN;
  css=AS*AS;
  csv00=AVplus*AVplus;
  csv11=AVminus*AVminus;
  csv01=AVplus*AVminus;

  return nucleusRecoil_stat(Mcdm,fv,A,Z,J,Sxx,
                                css,csv00,csv01,csv11,dNdE);
} 

double nucleusRecoil0CS( 
      double(*fv)(double),
      int A, int Z, double J,
      double Sp,double Sn, 
      double siP, double siN, double sdP,   double sdN,
      double * dNdE)
{
  double AS,AVplus, AVminus, css,csv00,csv01,csv11; 
  double LmbdP,XiP,LmbdN,XiN;
  double MN=0.939;
  double Mr=MN*Mcdm/(MN+Mcdm);
  double sCoeff= 2*sqrt(3.8937966E8/M_PI)*Mr;
 
  LmbdP=sqrt(fabs(siP))/sCoeff; if(siP<0) LmbdP*=-1;
  XiP  =sqrt(fabs(sdP)/3)/sCoeff; if(sdP<0) XiP*=-1;
  LmbdN=sqrt(fabs(siN))/sCoeff; if(siN<0) LmbdN*=-1;
  XiN  =sqrt(fabs(sdN)/3)/sCoeff; if(sdN<0) XiN*=-1;

  AS=Z*LmbdP+(A-Z)*LmbdN;
  AVplus =XiP+XiN;
  AVminus=XiP-XiN;
  css=AS*AS;
  csv00=AVplus*AVplus;
  csv11=AVminus*AVminus;
  csv01=AVplus*AVminus;

  return nucleusRecoil0_stat(fv,A,Z,J,Sp,Sn,
                                css,csv00,csv01,csv11,dNdE);
} 


double nucleusRecoil0(double(*fv)(double),
int A, int Z, double J,double Sp,double Sn,
double * dNdE)
{
  if(J)
  { double A3=pow(A,1./3.);
    RS_= 1.7*A3 -0.28 - 0.78*(A3-3.8 + sqrt((A3-3.8)*(A3-3.8) +0.2) );
      
    S00_0=  (Sp+Sn)*(Sp+Sn)*(2*J+1)*(J+1)/(4*M_PI*J);
    S11_0=  (Sp-Sn)*(Sp-Sn)*(2*J+1)*(J+1)/(4*M_PI*J);
    S01_0=2*(Sp+Sn)*(Sp-Sn)*(2*J+1)*(J+1)/(4*M_PI*J);
  }  
  return nucleusRecoil(fv,A,Z,J,Sxx_,dNdE);
}
 



double dNdERecoil(double E, double *tab)
{  double kE,alpha;
   int k;
   if(E<RE_START) return 0;
   kE= log(E/RE_START)/log(RE_STEP);
   if(kE<0 || kE>RE_DIM-1) return 0;
   k=kE;
   if(k>RE_DIM-2) k=RE_DIM-2;
   alpha= kE-k;
   if(tab[k]>0 && tab[k+1]>0) return pow(tab[k],1-alpha)*pow(tab[k+1],alpha);
   return (1-alpha)*tab[k]+alpha*tab[k+1];
}



double cutRecoilResult(double *tab, double E1, double E2)  // ????
{ int i,i1,i2;
  double sum,dx;
  if(E2>RE_STEP*(RE_DIM-1)) 
  { 
    printf("cutRecoilResult:  Maximal bound %.3E is larger than data limit %.3E\n",
       E2,RE_STEP*(RE_DIM-1));
    E2=RE_STEP*(RE_DIM-1); 
  } 
  if(E1<0)   E1=0;
  if(E1>=E2) return 0;
  
  i1=E1/RE_STEP;  i2=E2/RE_STEP;
  if(i1<E1/RE_STEP) i1++;
  
  for(i=i1,sum=0;i<=i2;i++) sum+=tab[i];
  sum-=(tab[i1]+tab[i2])/2;
  dx=i1 -E1/RE_STEP;
  if(dx>0) sum+=tab[i1]*dx + (tab[i1-1]-tab[i1])*dx*dx/2  ;
  dx=E2/RE_STEP-i2;
  if(dx) sum+=tab[i2]*dx  +(tab[i2+1]-tab[i2])*dx*dx/2;
  
  return sum*RE_STEP;
}









#ifdef INPROGRESS

static double nuSpect(double Enu/*MeV*/)  { return 1; /*1/MeV/cm^2/s*/}
static double nuSpectF(double Enu) { return 1000*nuSpect(Enu*1000)*(1-MAEr/(2*Enu*Emu));}

static double MAEr;


void neutrinoNucleusRecoil(
      int A, int Z, double J, void (*Sxx)(double, double *,double *, double *),
 
      double * dNdE)
{
  const double  Kg=1./1.782662E-27;          /* GeV  */
  const double  vC=299792;                   /* km/s */
  const double lDay=60.*60.*24.;             /* sec  */  

  const double step=1.E-6*RE_STEP;              /* 1KeV */

  int i;
  double MA,ffs,Rs,sum; 
  double FFs=1, s00,s01,s11;
  double E0,vmin,vmax;

  else vmax=1200.;
  
  MA=0.94*A;

  

  if(A>1)
  {
    Rs=C_*pow(A,1./3.)+B_;
    R_=Rs; a_=A_; 
    ffs=simpson(FermiND0,0., R_+10*0.5, 1.E-4,NULL);
  }
  
  for(i=0,sum=0;i<RE_DIM;i++)
  { double E=RE_START*1E-6*pow(RE_STEP,i);

    MAEr=MA*E;
    
    double r=simpson(nuSpectF, sqrt(MA*E/2),1/*GeV*/, 1E-3,NULL)*lDay; 
    double p=sqrt(E*MA*2)/0.197327;

    if(A>1)
    { 
      p_=p;    
      R_=Rs; a_=A_; 
      FFs=simpson(FermiNDP,0., R_+10*0.5, 1.E-4,NULL)/p_/ffs;      
    } else FFs=1;

    dNdE[i]=FFs*FFs*r*MA*Kg;

  }
}

#endif 