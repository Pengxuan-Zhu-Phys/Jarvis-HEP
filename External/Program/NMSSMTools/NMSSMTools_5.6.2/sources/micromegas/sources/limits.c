
#include"micromegas.h"
#include"micromegas_aux.h"

int Zinvisible(void)
{
  REAL mass[2]={Mcdm1,Mcdm2};
  char*  name[2]={CDM1,CDM2};
  char*  aname[2]={aCDM1,aCDM2};
  char* zname;
  txtList L;
  int i;
  double width,br;
  
  zname=pdg2name(23);
  if(!zname) { printf("Z boson  is absent\n"); return 0;}   
  width=pWidth(zname,&L);
  for(i=0,br=0;i<2;i++) if(name[i] && 2*mass[i]< pMass(zname))
  { char channel[40];
    sprintf(channel,"%s %s",name[i],aname[i]);
    br+=findBr(L,channel);   
  }
  
//  printf(" partial width Z->DM=%.1EMeV\n", br*1000*width);  
  if(br*1000*width>0.5)
  { printf(" partial width Z->DM,DM is %.2EMeV,  more than 0.5 MeV. See 1401.2447\n",br*1000*width);
    return 1;
  } else return 0;
}

extern int zinvisible_(void); //Fortran 
int zinvisible_(void) { return Zinvisible();}


int LspNlsp_LEP(double *cs_out)
{  
  char *e=pdg2name(11),*e_=pdg2name(-11),*z_=pdg2name(23);
  double res=0;
  double P=104;
  int VZ_tmp=VZdecay;
  int m;
  double resMax=0,csLim=0.5;

  if(!e || !e_ || !z_) return 0;
  if(!VZ_tmp) { VZdecay=1; cleanDecayTable();}

  for(m=0;m<2;m++) // cycle over two DM sectors
  {   
    char *CDM,*aCDM;
    double M,M_;    
    int n;
    
    if(m){ CDM=CDM2; M=Mcdm2;} else {CDM=CDM1;M=Mcdm1;}
    if(!CDM) continue;

    aCDM=pdg2name(-pNum(CDM));
    if(aCDM && strcmp(aCDM,CDM)==0) aCDM=NULL;
    
    for(n=1;;n++) // looking for next DM
    { 
      double width,brZ,brZu,brZd;
      char chan[50];
      txtList L;
      char *CDM_,*aCDM_;
      
      CDM_= nextOdd(n,&M_);
      
//printf("CDM_=%s\n", CDM_);     

      if(!CDM_ || M+M_>2*P)break; 
      if( (CDM[1]=='~' && m==0 ) || (CDM[1]!='~' && m!=0  ) ) continue;
      width=pWidth(CDM_,&L);
      sprintf(chan,"%s,%s", CDM,z_);
      brZ=findBr(L,chan);
      if(brZ==0 && aCDM)
      { sprintf(chan,"%s,%s", aCDM,z_);
        brZ=findBr(L,chan);
      } 
      if(brZ >0)
      {  double brZu=brZ*0.119, brZd=brZ*0.152;
         char process[100];
         char lib[40];
         char buff[20];
         int k,err;
         numout*cc;        

         brZ=2*brZd+brZu;
         if(M_-M > 2*1.5) brZ+=brZu; // c-quark
         if(M_-M > 2*5)   brZ+=brZd; // b-quark   
         sprintf(process,"%s,%s->_CDM_,_NCDM_{%s",e,e_,CDM);
         if(aCDM) sprintf(process+strlen(process),",%s",aCDM);
         sprintf(process+strlen(process),"{%s",CDM_);
         aCDM_=pdg2name(-pNum(CDM_));
         if(aCDM_ && strcmp(aCDM_,CDM_)==0) aCDM_=NULL;
         if(aCDM_) sprintf(process+strlen(process),",%s",aCDM_);
//printf("process=%s\n",process);
         sprintf(lib,"ee_");
         pname2lib(CDM,buff);
         strcat(lib,buff); 
         strcat(lib,"_");
         pname2lib(CDM_,buff);
         strcat(lib,buff); 
         cc= getMEcode(0,ForceUG, process, NULL,NULL,lib);
//printf(" cc->interface->nprc=%d\n",cc->interface->nprc);
         double dres=0;
         for(k=1;k<=cc->interface->nprc;k++) 
         { dres+=brZ*cs22(cc,k,P,-1,1 ,&err);} 
         res+=dres;
         if(dres>resMax) { resMax=dres; if(M_>100) csLim=0.1; else csLim=0.5; }
      }   
    }                 
  }
  if(VZ_tmp!=VZdecay) { VZdecay=VZ_tmp; cleanDecayTable();}
  if(cs_out) *cs_out=res;
  if( res>csLim) return 1; else return 0;    
}

int  lspnlsp_lep_(double *cs_out) { return  LspNlsp_LEP(cs_out);}