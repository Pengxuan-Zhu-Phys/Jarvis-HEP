#include "pmodel.h"
#include "pmodel_aux.h"
#include "../../include/micromegas.h"
#include "../../include/micromegas_aux.h"
/*#include"../../CalcHEP_src/c_source/model_aux/include/SLHAreader.h"*/

#include"lpath.h"


static int runTools(char * cmd, char * fout,int mode)
{  char * command;
   int err;
   double maxl;
   
   if(access(fout,F_OK)==0) unlink(fout);
   command=malloc(100+strlen(NMSSMTOOLS ));

   sprintf(command,
    "lPath=%s;EXPCON_PATH=$lPath/EXPCON;export EXPCON_PATH;$lPath/main/%s"
    ,NMSSMTOOLS ,cmd); 
//printf("command:%s\n", command);    
   err=System(command);
   free(command);
   if(err>=0) err=slhaRead(fout,0);

   return err; 
}
     
/* =====  end of header part ========= */

  double ewsbNMSSM(double tb,    double MG1, double MG2, double MG3,   
  double Ml2, double Ml3, double Mr2, double Mr3, double Mq2, double Mq3,   
  double Mu2, double Mu3, double Md2, double Md3, 
  double At, double Ab, double Al, double mu,     
  double LambdQ,double KappaQ,double aLmbdQ, double aKappQ,
  double mXiF,double mXiS, double muP, double msP,double m3h)
{ 

  int prec=2;
  if(mXiF || mXiS || muP|| msP||m3h) prec=0;
  
   int err,nw;
   FILE*  f=fopen("inp","w");
   if(f==NULL) return -1;
                                       
   fprintf(f,"Block MODSEL    # Select model\n"   
             "  1    0           # EWSB input\n"
             "  3    1           # NMSSM PARTICLE CONTENT\n"
             "  8    %d           # Precision for Higgs masses\n"
             "  9    0           # FLAG FOR MICROMEGAS (0=NO, 1=YES\n"
             " 10    0           # No scan, no ...\n"
             " 13    1           # Sparticles decays\n"
           "Block SMINPUTS    # Standard Model inputs\n",prec);
//   fprintf(f," 1   %.8E       # alpha_em^(-1)(MZ) SM MSbar\n",1/findValW("alfEMZ"));
//   fprintf(f," 2   %.8E       # G_Fermi \n",1.16637E-5); 
   fprintf(f," 3   %.8E       # alpha_s(MZ) SM MSbar\n",findValW("alfSMZ"));
   fprintf(f," 5   %.8E       # mb(mb) SM MSbar\n", findValW("MbMb"));
   fprintf(f," 6   %.8E       # mtop(pole)\n",      findValW("Mtp"));
   fprintf(f," 7   %.8E       #  Mtau     \n",      findValW("Ml"));

   fprintf(f,"Block MINPAR                 # Input parameters\n");
   fprintf(f," 3   %.8E       # tanb\n",     tb);

   fprintf(f,"Block EXTPAR\n");
   fprintf(f," 0   -1        # EWSB\n"); 
   fprintf(f," 1   %.8E      # MG1\n",MG1);
   fprintf(f," 2   %.8E      # MG2\n",MG2);
   fprintf(f," 3   %.8E      # MG3\n",MG3);

   fprintf(f," 11  %.8E      # At \n",At);
   fprintf(f," 12  %.8E      # Ab \n",Ab);
   fprintf(f," 13  %.8E      # Atau\n",Al);
/*   fprintf(f," 16  %.8E      # Am \n",Am);*/
/*   fprintf(f," 23  %.8E      # mu\n",mu); */
/*   fprintf(f," 24  %.8E      # MA\n",MH3);*/
/*   fprintf(f," 26  %.8E      # MA\n",MH3);*/
   fprintf(f," 31  %.8E      # Ml1\n",Ml2);
   fprintf(f," 32  %.8E      # Ml2\n",Ml2); 
   fprintf(f," 33  %.8E      # Ml3\n",Ml3);
   fprintf(f," 34  %.8E      # MR2\n",Mr2);
   fprintf(f," 35  %.8E      # MR2\n",Mr2); 
   fprintf(f," 36  %.8E      # MR3\n",Mr3);

   fprintf(f," 41  %.8E      # Mq1\n",Mq2);
   fprintf(f," 42  %.8E      # Mq2\n",Mq2); 
   fprintf(f," 43  %.8E      # Mq3\n",Mq3);
   fprintf(f," 44  %.8E      # Mu1\n",Mu2);
   fprintf(f," 45  %.8E      # Mu2\n",Mu2); 
   fprintf(f," 46  %.8E      # Mu3\n",Mu3);
   fprintf(f," 47  %.8E      # Md1\n",Md2);
   fprintf(f," 48  %.8E      # Md2\n",Md2); 
   fprintf(f," 49  %.8E      # Md3\n",Md3);

   fprintf(f," 61  %.8E      # L\n",  LambdQ);
   fprintf(f," 62  %.8E      # K\n",  KappaQ);
   fprintf(f," 63  %.8E      # AL\n", aLmbdQ);
   fprintf(f," 64  %.8E      # AK\n", aKappQ);
   fprintf(f," 65  %.8E      # MU\n", mu);

   fprintf(f," 66  %.8E        # XiF \n", (double)(mXiF*fabs(mXiF)));
   fprintf(f," 67  %.8E        # XiS \n", (double)(pow(mXiS,3)));
   fprintf(f," 68  %.8E        # muP \n",  muP);
   fprintf(f," 69  %.8E        # MS'^2 \n", (double)(msP*fabs(msP)));
   fprintf(f," 72  %.8E        # M3H^2\n", (double)(m3h*fabs(m3h)));
     
   fclose(f);

   err= runTools("nmhdecay","spectr",0);

   if(err) {FError=1;}
//   nw= slhaWarnings(NULL);
   return err;
}

double sugraNMSSM( double m0, double mhf, double a0, double tb, double sgn,  
     double Lambda,double aLambda, double aKappa, double mXiF, double mXiS,
     double muP, double msP,double m3h)
{
  int nw=0,err,prec=2;

  if(mXiF || mXiS || muP|| msP||m3h) prec=0;
  FILE*  f=fopen("inp","w");  
  if(f==NULL) return -1;
   
   fprintf(f,
           "Block MODSEL       # Select model\n"
           "  1    1            # SUGRA\n"   
           "  3    1            # NMSSM PARTICLE CONTENT\n"
           "  8    %d           # Precision for Higgs masses\n"
           "  9    0            # FLAG FOR MICROMEGAS (0=NO)\n" 
           " 10    0            # No scan, no ...\n"
           " 13    1           # Sparticles decays\n"
           "Block SMINPUTS               # Standard Model inputs\n",prec);
//   fprintf(f," 1   %.8E       # alpha_em^(-1)(MZ) SM MSbar\n",1/findValW("alfEMZ"));
//   fprintf(f," 2   %.8E       # G_Fermi \n",1.16637E-5); 
   fprintf(f," 3   %.8E       # alpha_s(MZ) SM MSbar\n",findValW("alfSMZ"));
   fprintf(f," 5   %.8E       # mb(mb) SM MSbar\n", findValW("MbMb"));
   fprintf(f," 6   %.8E       # mtop(pole)\n",      findValW("Mtp"));
   fprintf(f," 7   %.8E       #  Mtau     \n",      findValW("Ml"));
   fprintf(f,"Block MINPAR                 # Input parameters\n");
   
   fprintf(f," 1   %.8E       # m0\n",       m0); 
   fprintf(f," 2   %.8E       # m1/2\n",     mhf);   
   fprintf(f," 3   %.8E       # tanb\n",     tb);
   fprintf(f," 4   %.0f       # sign(mu)\n", sgn);
   fprintf(f," 5   %.8E       # A0\n",       a0);

   fprintf(f,"Block EXTPAR\n");
   fprintf(f," 61  %.8E        # L \n",  Lambda);
   fprintf(f," 63  %.8E        # A_LAMBDA\n", aLambda);
   fprintf(f," 64  %.8E        # A_K\n", aKappa);

   fprintf(f," 66  %.8E        # XiF \n",  mXiF*((double)fabs(mXiF)));
   fprintf(f," 67  %.8E        # XiS \n",  mXiS*mXiS*mXiS);
   fprintf(f," 68  %.8E        # muP \n", muP);
   fprintf(f," 69  %.8E        # MS'2 \n", msP*((double)fabs(msP)));
   fprintf(f," 72  %.8E        # M3H^2\n", m3h*((double)fabs(m3h)));
                           
   fclose(f);
     
  err= runTools("nmspec","spectr",0);
  if(err) {FError=1;} 
//  nw= slhaWarnings(NULL);
//  if(nw==-1) FError=1;
  return err;
}



