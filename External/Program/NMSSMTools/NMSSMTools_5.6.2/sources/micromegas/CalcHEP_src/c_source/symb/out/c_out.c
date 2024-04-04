/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include <unistd.h>
#include "chep_crt.h"
#include "getmem.h"
#include "syst2.h" 
#include "physics.h"
#include "s_files.h"
#include "procvar.h"
#include "pvars.h"
#include "diaprins.h"
#include "optimise.h"
#include "l_string.h"
#include "parser.h"
#include "reader_c.h"
#include "out_serv.h"
#include "saveres.h"
#include "denominators.h"
#include "process.h"
#include "cweight.h"
#include "c_out.h"
#include "writeF.h"
#include "sos.h"
#include "../../../include/version.h"
#include "nType.h"

/* ======================================================== */
typedef char prtclsarray[MAXINOUT][P_NAME_SIZE];

int noCChain=0;
int tWidths=0;
static int sumDiag=0;

static void  getprtcls(char* txt1,prtclsarray pnames)
{  int   j; 
   char   txt[STRSIZ]; 
   char*c0,*c1;
   
   strcpy(txt,txt1);
   c0=strstr(txt,"->"); 
   c0[0]=','; 
   c0[1]=' '; 

   for(j=nin+nout; j<MAXINOUT; j++) strcpy(pnames[j],"***"); 

   for(c0=txt,j=0; ;j++)
   { 
      c1=strchr(c0,',');      
      if(c1) c1[0]=0;
      trim(c0);
      strcpy(pnames[j],c0);
      if(!c1) return; else c0=c1+1; 
   }
} 


/*===========================================================*/


#define procinfptr struct procinfrec *
typedef struct procinfrec
   {
      procinfptr     next;
      int            tot;
      long           firstdiagpos;
      prtclsarray    p_name;
      int            p_masspos[MAXINOUT];
      int            p_code[MAXINOUT];
   }  procinfrec;
#undef procinfptr

typedef struct procinfrec *procinfptr;

static procinfptr   inf, inftmp;  /*  information about subProcess  */

       /*  statictics  */
static unsigned  ndiagrtot,  diagrcount;
static int  nvars,  nfunc;

static marktp   heapbeg;/*  for RELEASE  */

static int   nden_w, nden_s, nden_t, nden_0,  nsub1; /* from writesubprocess */


static int cBasisPower;
static int nC, *cChains=NULL;
static long  *cCoefN, *cCoefD;


static void clearstatistic(void)
{int  i; for (i = 17; i < 24; i++) { goto_xy(1,i); clr_eol();} }

static void init_stat(void)
{
   goto_xy(1,17);
   scrcolor(Yellow,Blue);
   print(" C Source Codes \n");
   scrcolor(Red,BGmain);
   print(" Process..........\n");
   print(" Total diagrams...\n");
   print(" Processed........\n");
   print(" Current..........\n");
   scrcolor(Yellow,Blue);
   print(" Press Esc to stop    ");
   scrcolor(Black,BGmain);
   goto_xy(20,18); print("%s",processch);
   goto_xy(20,19); print("%4u",ndiagrtot);
   goto_xy(20,20); print("   0");
   scrcolor(Yellow ,BGmain);
   goto_xy(20,21); print("   1");
   scrcolor(Yellow,BGmain);
}


static void writestatistic(void)
{
   scrcolor(Black ,BGmain);
   goto_xy(20,19); print("%4u",ndiagrtot);
   goto_xy(20,20);
   print("%2u (%%)",(((diagrcount - 1) * 100) / ndiagrtot));
   goto_xy(20,21); print("%4u",diagrcount);
}


static void writpict(unsigned ndiagr)
{  vcsect vcs;
   csdiagram  csdiagr;
   fseek(diagrq,ndiagr*sizeof(csdiagr),SEEK_SET);
   FREAD1(csdiagr,diagrq);
   transfdiagr(&csdiagr,&vcs);
   writeF("/*\n");
   DiagramToOutFile(&vcs,1,' ');
   writeF("*/\n");
}


static void labl(void)
{
  writeF("/*******************************\n");
  writeF("*    %s*\n",VERSION_);
  writeF("*******************************/\n");
}

  /* =========== Preliminary  calculations ================ */

static void calc_nvars_nfunc(void)
{ int   k;
  for(nvars=0, nfunc=0, k=1; k<=nmodelvar; k++) if(vararr[k].used)
  { if(k>nCommonVars && modelvars[k].pub==0) nfunc++; else  nvars++; }
} 


static void prepareprocinform(void)
{int ndel, ncalc, nrest;
 long recpos;
 char        txt[STRSIZ];
 int     i, k;
 csdiagram   csd;
 char * mass;
 int        nn;
 int nsubs;

   inf = NULL;
   fseek(menuq,0,SEEK_SET);
   for (nsubs=1;nsubs<=subproc_sq;nsubs++)
   {
      inftmp = inf;
      inf = (procinfptr)getmem_((unsigned)sizeof(procinfrec));
      inf->next = inftmp;
      rd_menu(2,nsubs,txt,&ndel,&ncalc,&nrest,&recpos);
      inf->firstdiagpos = recpos;
      getprtcls(txt,inf->p_name);
      for (i = 0; i < nin + nout; i++)
      {
         locateinbase(inf->p_name[i],&nn);
         mass=prtclbase[nn-1].massidnt; 
         if (strcmp(mass,"0")) {for(k=1; k<=nmodelvar; k++) if(strcmp(modelvars[k].varname,mass)==0) break;} 
         else k=0;
         if(k>nmodelvar){ printf("mass %s is not defined\n",mass); exit(125);} 
         inf->p_masspos[i] = k;
         inf->p_code[i]=prtclbase1[nn].N; 
      }
      for (i = nin + nout; i < MAXINOUT; i++)
      {
         strcpy(inf->p_name[i],"***");
         inf->p_masspos[i] = 0;
      }
      fseek(diagrq,recpos*sizeof(csdiagram),SEEK_SET);
      inf->tot = 0;
      for (i = 1; i <= ndel + ncalc + nrest; i++)
      {
         FREAD1(csd,diagrq);
         if (csd.status == 1) ++(inf->tot);
      }
      if(inf->tot==0) for (i = 0; i < nin + nout; i++) inf->p_masspos[i]=0;
   }
   nsubs--;
   revers((void **)&inf);
}


static void sortvars(void)
{
  int k;

  for(k=4;k<=nmodelvar;k++) if(vararr[k].used && (k<=nCommonVars||modelvars[k].pub)) 
            writeF("\n,\"%s\"",modelvars[k].varname);
  for(k=1+nCommonVars;k<=nmodelvar;k++)if (vararr[k].used && !modelvars[k].pub)
            writeF("\n,\"%s\"",modelvars[k].varname);     
}

/* ======= Information functions =========== */

static void geninf(char* name,int value)
{  writeF("const int %s = %d;\n\n",name, value); }


static void  writesubroutineinit(void)
{
   int         l;
   char        *ss;
   
   ext_h=NULL;   
   writeF("double (*aWidth_ext)(char*)=NULL;\n"); 
   writeF("static int calcFunc_stat(void)\n{\n");
   writeF(" REAL * V=va_ext;\n");
   writeF(" FError=0;\n");
   for(l=nCommonVars+1;l<=nmodelvar;l++)
   {
      if(vararr[l].used && ((modelvars[l].func && modelvars[l].pub==0) || modelvars[l].pwidth) )
      {  int num;
         checkNaN=0;
         if(modelvars[l].pwidth)
         {  writeF("   %s=aWidth_ext(\"%s\");\n",vararr[l].alias,
             prtclbase1[modelvars[l].pwidth].name); 
             checkNaN=1;
         } else
         { ss=(char *)readExpression(modelvars[l].func,rd_c,act_c,free);
/*	   writeF("   %s=%s;\n",vararr[l].alias,ss+3);*/
	   fprintf(outFile,"   %s=%s;\n",vararr[l].alias,ss+3);
	   free(ss);
	 }
	 if(checkNaN)
	 {
	 sscanf(vararr[l].alias,"V[%d]",&num); 
/*
if(!modelvars[l].pwidth)	 
           fprintf(outFile,"printf(\"%s= %%50.50s\\n\",\"%s\");\n",modelvars[l].varname,  modelvars[l].func);
writeF("printf(\"%s=%%E\\n\",%s);\n",vararr[l].alias,vararr[l].alias);
*/	
         writeF("   if(!isfinite(%s) || FError){ return %d;}\n",vararr[l].alias,num);
         }
      }
   }

   writeF("return 0;\n}\n"); 
}

static void  writeDenominators(deninforec* dendescript, int forAll)
{   
   int i,k;
   int nden_t2=nden_t, nden_s2=nden_s;
   char *gsw_txt,*gtw_txt;
   if(forAll) {gsw_txt="gswidth_ext"; gtw_txt="gtwidth_ext";} 
       else   {gsw_txt="gsw";         gtw_txt="gtw";} 

   for (i = 0; i < dendescript->tot_den; i++)
   {  int numm = dendescript->denarr[i].order_num; 
      if(dendescript->denarr[i].width==0)  numm += nden_w; else  
      { if(dendescript->denarr[i].stype)  nden_s2--; 
        else { numm += nden_s; nden_t2--;}
      }
      if(dendescript->denarr[i].power == 2) writeF("*Q2[%d]",numm); else 
      { 
        if(dendescript->denarr[i].stype) 
        {
          writeF("*(%s ? creal(Q1[%d]):",gsw_txt,numm);
          if (dendescript->denarr[i].power==1) writeF("Q1[%d])",numm);
          else  writeF("conj(Q1[%d]))",numm);
        } else
        { 
          writeF("*(%s ? creal(Q1[%d]):",gtw_txt,numm);
          if (dendescript->denarr[i].power==1) writeF("Q1[%d])",numm);
          else  writeF("conj(Q1[%d]))",numm);
        }        
      }  
   }

   writeF(";\n");

   if(nden_s2+nden_t2)
   {   
      if(nden_s2) 
      { 
        writeF("if(%s) Prop=Prop",gsw_txt);  
        for (k = 1; k <= nden_s; k++)					  	
        { int addpr = 1;							       
          for (i = 1; i <= dendescript->tot_den; i++)			       
          if (dendescript->denarr[i-1].width &&	dendescript->denarr[i-1].stype&&			       
   	     k == dendescript->denarr[i-1].order_num)  addpr = 0;	  	
          if (addpr)  writeF("*Q0[%d]",k);
        }
      writeF(";\n");
      }  
      if(nden_t2) 
      {  writeF("if(%s) Prop=Prop",gtw_txt); 
        for (k = nden_s+1; k <= nden_w; k++)					  	
        { int addpr = 1;							       
          for (i = 1; i <= dendescript->tot_den; i++)			       
          if (dendescript->denarr[i-1].width &&	(!dendescript->denarr[i-1].stype) &&			       
   	     k == dendescript->denarr[i-1].order_num+nden_s)  addpr = 0;	  	
          if (addpr)  writeF("*Q0[%d]",k);
        }
       writeF(";\n");
      }      					       
   }									       
}

static void calcColor(long diag)
{
   csdiagram  csdiagr;

   fseek(diagrq, (diag-1)*sizeof(csdiagr),SEEK_SET);
   FREAD1(csdiagr,diagrq);

   if(cBasisPower&&generateColorWeights(&csdiagr,cBasisPower,nC,cChains,cCoefN,cCoefD))
   {  
      int k;
      writeF(" if(cb_coeff)\n {\n");
           
      for(k=0; k<cBasisPower; k++)  if(cCoefN[k]) 
      { int i;
         writeF("  cb_coeff[%d] += (R*%d)/(%d); %s ",k, cCoefN[k],cCoefD[k],"/*");
         for(i=0;i<nC;i++) if(cChains[4*(nC*k+i)]==2)  writeF("(%d %d) ", cChains[4*(nC*k+i)+1], cChains[4*(nC*k+i)+2] );
         else writeF("(%d %d %d) ", cChains[4*(nC*k+i)+1], cChains[4*(nC*k+i)+2],cChains[4*(nC*k+i)+3] );
         writeF("*/\n");
      } 
      writeF(" }\n");
   }
}

static void  onediagram(deninforec* dendescript)
{  catrec      cr;
   marktp      bh;
   varptr      totnum, totdenum, rnum;
   long pos_c;
   int deg1,nConst;

   mark_(&bh);
   tmpNameMax=0;
   initinfo();
   initdegnames();
   
   fseek(catalog,dendescript->cr_pos,SEEK_SET);
   FREAD1(cr,catalog);
   ++(diagrcount);
   whichArchive(cr.nFile,'r');

   fseek(archiv,cr.factpos,SEEK_SET);
    
   readvardef(archiv);
   readpolynom(&totnum);
   readpolynom(&totdenum);
   clearvardef();
   
   fseek(archiv,cr.rnumpos,SEEK_SET);

   readvardef(archiv);
   readpolynom(&rnum);
   clearvardef();

   {  outFileOpen("%sresults%cf%d.c",pathtouser,f_slash,diagrcount);
      labl();
      writeF("#include\"num_out.h\"\n");
      writeF("#include\"num_in.h\"\n");   
   }
   writeF("extern FNN F%d_ext;\n",diagrcount);
   writeF("static void C%d(REAL*V,REAL * C)\n{\n",diagrcount);
   pos_c= ftell(outFile); writeF("%80s\n",""); 
   nConst=write_const();
   deg1=cleardegnames();       
   writeF("}\n"); 

   fseek(outFile,pos_c,SEEK_SET);
   if(deg1) writeF("REAL S[%d];",deg1);
   if(tmpNameMax) writeF("REAL tmp[%d];",tmpNameMax );                                

   fseek(outFile,0,SEEK_END);
   tmpNameMax=0;
   initdegnames();

   writeF("REAL F%d_ext(double GG, REAL*V, REAL*DP,REAL*Q0,COMPLEX*Q1,REAL*Q2,REAL*cb_coeff,int gsw,int gtw)\n{\n",diagrcount);

   if(!noPict) writpict(cr.ndiagr_ + inftmp->firstdiagpos - 1);

   writeF("REAL N,D,R; COMPLEX Prop;\n");
   pos_c= ftell(outFile); writeF("%80s\n","");
  
   writeF("if(!DP){C%d(V,C); return 0;} \n",diagrcount);
   if(nin==2) writeF("  REAL N_p1p2_=1/DP[0];\n");
   fortwriter("N",totnum);
   fortwriter("D",totdenum);
   fortwriter("R",rnum);
   
   writeF("R*=(N/D);\n");
   writeF("Prop=1");
   writeDenominators(dendescript,0);
   writeF("R*=creal(Prop);\n");
   if(!noCChain)calcColor(cr.ndiagr_+inftmp->firstdiagpos);

   writeF(" return R;\n");  
   writeF("}\n");

   deg1=cleardegnames();
   if(nConst==0) nConst=1;
   fseek(outFile,pos_c,SEEK_SET);
   writeF("static REAL C[%d];",nConst);
   if(deg1) writeF("REAL S[%d];",deg1);
   if(tmpNameMax) writeF("REAL tmp[%d];",tmpNameMax );
   fseek(outFile,0,SEEK_END);
   outFileClose();
   release_(&bh);
}


static int  alldiagrams(FILE * fd,  int nsub)
{  
   marktp     bh;
   varptr     totnum, totdenum, rnum;
   long       pos_c1,pos_c2; int deg1,deg2,tmpn1,tmpn2, nC;
   catrec     cr;
   deninforec dendescript;

   mark_(&bh); tmpNameMax=0; initinfo(); initdegnames();

   writeF("{\n");
   writeF("REAL N,D,R; COMPLEX Prop;\n");
   pos_c1= ftell(outFile); writeF("%70s\n","");   
   writeF("if(!momenta){ C%d(V,C); return 0;}\n",nsub);
   if(nin==2) writeF("  REAL N_p1p2_=1/DP[0];\n");
   while(FREAD1(dendescript,fd) == 1)
   {
      fseek(catalog,dendescript.cr_pos,SEEK_SET);
      FREAD1(cr,catalog); ++(diagrcount);
      if(!noPict)writpict(cr.ndiagr_ + inftmp->firstdiagpos - 1);
      whichArchive(cr.nFile,'r');
      fseek(archiv,cr.factpos,SEEK_SET);

      readvardef(archiv);
      readpolynom(&totnum);
      readpolynom(&totdenum);
      clearvardef();
   
      fseek(archiv,cr.rnumpos,SEEK_SET);

      readvardef(archiv);
      readpolynom(&rnum);
      clearvardef();
   
      fortwriter("N",totnum);
      fortwriter("D",totdenum);

      fortwriter("R",rnum);

      writeF("R*=(N/D);\n");
      if(nin+nout>3)
      {  writeF("Prop=1");
         writeDenominators(&dendescript,1);
         writeF("R*=creal(Prop);\n");
         writeF(" if(R>(*Fmax)) (*Fmax)=R; else if(R<-(*Fmax)) (*Fmax)=-R;\n");
      } else  writeF(";\n");
      if(!noCChain)calcColor(cr.ndiagr_+inftmp->firstdiagpos);
      writeF("ans+=R;\n");
   }   
   whichArchive(0,0);
   writeF("\n}\nreturn ans;\n}\n");

   deg1=cleardegnames();
   tmpn1=tmpNameMax;
   tmpNameMax=0;
   initdegnames();

   writeF("\nstatic void C%d(REAL*V,REAL*C)\n{\n",nsub); 
   pos_c2= ftell(outFile); writeF("%70s\n","");   

   nC=write_const(); 
   if(nC==0) nC=1; 
   writeF("}\n");

   fseek(outFile,pos_c1,SEEK_SET); 
   writeF("static REAL C[%d];",nC); 
   if(deg1) writeF("REAL S[%d];",deg1);
   if(tmpn1) writeF("REAL tmp[%d];",tmpn1); 

 
   fseek(outFile,pos_c2,SEEK_SET);
   deg2=cleardegnames();
   tmpn2=tmpNameMax;
   if(deg2) writeF("REAL S[%d];",deg2) ;
   if(tmpn2) writeF("REAL tmp[%d];",tmpn2 );
   fseek(outFile,0,SEEK_END);

   release_(&bh);
   if( escpressed()) return 1; else return 0;
}



static void  writesubprocess(int nsub,long firstDiag,long totDiag,int* breaker)
{  denlist    den_;
   int      i;
    
   deninforec   dendescript;
   FILE * fd;                /* file of (deninforec)  */
   char fd_name[STRSIZ];
   marktp mem_start;

   nsub1 = nsub;

   { outFileOpen("%sresults%cd%d.c",pathtouser,f_slash,nsub);
     labl();
     writeF("#include\"num_in.h\"\n");
     writeF("#include\"num_out.h\"\n");
   }

   if(totDiag==0) 
   { writeF("extern DNN S%d_ext;\n",nsub); 
     writeF("REAL S%d_ext(double GG,  REAL * momenta, REAL*cb_coeff, double*Fmax, int * err)\n{",nsub); 
     writeF("  return 0;\n}\n");
     outFileClose(); 
     return;
   }

   if(sumDiag) writeF("static void C%d(REAL*,REAL *);\n",nsub); else 
   {  writeF("extern FNN F%d_ext",firstDiag);
      for(i=1;i<totDiag;i++) writeF(",F%d_ext",i+firstDiag);
      writeF(";\n");
      writeF("static FNN *Farr[%d]={&F%d_ext",totDiag,firstDiag);
      for(i=1;i<totDiag;i++) writeF(",&F%d_ext",i+firstDiag);
      writeF("};\n");
   } 
   writeF("extern DNN S%d_ext;\n",nsub);
   writeF("REAL S%d_ext(double GG, REAL * momenta, REAL*cb_coeff, double*Fmax, int * err)\n{",nsub);
   writeF("REAL  ans=0;\n");

   sprintf(fd_name,"%stmp%cden.inf",pathtouser,f_slash);
   fd=fopen(fd_name,"wb"); 

   mark_(&mem_start);
   denominatorStatistic(nsub, &nden_s, &nden_t, &nden_0, &den_, fd); 
   fclose(fd);
   nden_w=nden_s+nden_t;
   writeF("REAL DP[%d];\n",((nin+nout)*(nin+nout-1))/2);
   writeF("REAL* V=va_ext;\n");
   if(nin+nout>3)
   {  int nden= nden_w+nden_0+1; 
      writeF("REAL mass[%d],width[%d];\n",nden,nden);
      writeF("char * Qtxt[%d];\n",nden);
      writeF("REAL Q0[%d]; COMPLEX Q1[%d]; REAL Q2[%d];\n",nden_w+nden_0+1, nden_w+nden_0+1,nden_w+nden_0+1);
//      writeF("sprod_(%d, momenta, DP);\n",nin+nout);
/*      writeF(" for(i=0;i<nin_ext;i++) s0max+=momenta[4*i];\n"); */

      if(sumDiag) writeF(" if(momenta)\n {"); else
      writeF("  if(!momenta) {int i; for(i=0;i<%d;i++) Farr[i](GG,va_ext,NULL,NULL,NULL,NULL,NULL,0,0); return 0;}\n",totDiag); 
                          
      for(;den_;den_ = den_->next)
      {  int m=0;
         i=den_->order_num;
         if(den_->width) 
         {
           if(den_->stype) fprintf(outFile,"width[%d]=%s; ",i,vararr[den_->width].alias);
           else 
           {  i+=nden_s;
             fprintf(outFile,"width[%d]=(twidth_ext)? %s : 0.; ",i,vararr[den_->width].alias); 
           }
         }else 
         { i+=nden_w;
           fprintf(outFile,"width[%d]=0.; ",i);
         }
         fprintf(outFile,"mass[%d]=%s; ",i,vararr[den_->mass].alias);
         fprintf(outFile," Qtxt[%d]=\"",i);       
/*         fprintf(outFile," Q[%d]=mass[%d]*mass[%d]-sqrMom(nin_ext,\"",i,i,i);*/

         while(den_->momStr[m]) fprintf(outFile,"\\%o",den_->momStr[m++]);
         fprintf(outFile,"\";\n");
/*         fprintf(outFile,"\",momenta);\n");   */ 
      }  
      writeF("*err=*err|prepDen(%d,nin_ext,BWrange_ext*BWrange_ext,mass,width,Qtxt,momenta,Q0,Q1,Q2);\n",
      nden_w+nden_0);
   } else
   {
      if(sumDiag) writeF(" if(momenta)\n {"); else
      writeF("  if(!momenta) {int i; for(i=0;i<%d;i++) Farr[i](GG,va_ext,NULL,NULL,NULL,NULL,NULL,0,0); return 0;}\n",totDiag); 
   }   
   writeF("sprod_(%d, momenta, DP);\n",nin+nout);
   if(sumDiag) writeF("}\n");
   release_(&mem_start);
   fd=fopen(fd_name,"rb"); 
   if(sumDiag)       
   {   
     *breaker = alldiagrams(fd,nsub); 
     writestatistic();   
     outFileClose();    
   } else 
   { 
      writeF("{int i; for(i=0;i<%d;i++) \n",totDiag);
      writeF(
      "{ REAL r=Farr[i](GG,va_ext,DP,Q0,Q1,Q2,cb_coeff,gswidth_ext,gtwidth_ext);\n"
      "  if(r>(*Fmax)) *Fmax=r;\n"
      "  ans+=r;\n"
      "}}\n"
      "return ans;\n}\n"
            );

      outFileClose();

      *breaker = 0;
      while(FREAD1(dendescript,fd) == 1)
      {
         if (escpressed())
         {  *breaker = 1;
            break;
         }
         onediagram(&dendescript);
         writestatistic();
      } 
   }
   fclose(fd);
   unlink(fd_name);
}  /*  WriteSubprocess  */



static void  make_pinf(void)
{
   int    i;

   writeF("char * pinf_ext(int nsub,int nprtcl,REAL* pmass,int * num)\n{\n");
   writeF("int n;\n");

   writeF(" static char *names[%d][%d] ={\n",subproc_sq,nin + nout);
   inftmp = inf;
   for (nsub = 1; nsub <= subproc_sq; nsub++)
   {  writeF("{");
      for (i = 1; i <= nin + nout; i++)
      {  if(i!=1) writeF(",");
         writeF("\"%s\"",inftmp->p_name[i-1]);
      }
      if (nsub== subproc_sq) writeF("}};\n"); else  writeF("},\n");
      inftmp = inftmp->next;
   }
   
   writeF("int const nvalue[%d][%d]={\n",subproc_sq,nin + nout);
   inftmp = inf;
   for (nsub = 1; nsub <= subproc_sq; nsub++)
   {  writeF("{");
      for (i = 1; i <= nin + nout; i++)
      {  int k=inftmp->p_masspos[i-1];      
         if(k) 
         {
            if(vararr[k].used) sscanf(vararr[k].alias,"V[%d]",&k); else k=-1;
         }
         if(i!=1) writeF(","); 
	 writeF("%d",k);
      }
      if (nsub== subproc_sq) writeF("}};\n"); else  writeF("},\n"); 
      inftmp = inftmp->next;
   }

   writeF("int const pcode[%d][%d]={\n",subproc_sq,nin + nout);
   inftmp = inf;
   for (nsub = 1; nsub <= subproc_sq; nsub++)
   {  writeF("{");
      for (i = 0; i < nin + nout; i++)
      { 
         if(i) writeF(",");  
	 writeF("%d",inftmp->p_code[i]);
      }
      if (nsub== subproc_sq) writeF("}};\n"); else  writeF("},\n"); 
      inftmp = inftmp->next;
   }

   writeF("if  (nsub<0 ||nsub>%d||nprtcl<0||nprtcl>%d) return NULL;\n",
   subproc_sq,nin + nout);
   writeF("if(pmass)\n{\n");
   writeF("  n=nvalue[nsub-1][nprtcl-1];\n");
   writeF("  if (n==0) *pmass=0; else *pmass=va_ext[n];\n"); 
   writeF("  if (*pmass<0) (*pmass)=-(*pmass);\n");  
   writeF("}\n");  
   writeF("if(num)*num=pcode[nsub-1][nprtcl-1];\n");
   writeF("return names[nsub-1][nprtcl-1];\n}\n");

   if(nin==1) writeF("char * polarized_ext[3]={\"\",\"\",\"\"};\n");
   else 
   {  writeF("char * polarized_ext[3]={\"\",\",");       
      for (i = 1; i <= nparticles; i++)
      if(polarized(1,i)) writeF("%s,",prtclbase1[i].name);
      writeF("\",\",");
      for (i = 1; i <= nparticles; i++)
      if(polarized(2,i)) writeF("%s,",prtclbase1[i].name);
      writeF("\"};\n");
   }

   writeF("int pinfAux_ext(int nsub,int nprtcl,int*spin2,int*color,int*neutral,int*ndf)\n{\n");
/*   writeF("int n;\n"); */

   writeF("int const pcode[%d][%d]={\n",subproc_sq,nin + nout);

   inftmp = inf;
   for (nsub = 1; nsub <= subproc_sq; nsub++)
   {  writeF("{");
      for (i = 0; i < nin + nout; i++)
      { 
         if(i) writeF(",");  
	 writeF("%d",inftmp->p_code[i]);
      }
      if (nsub== subproc_sq) writeF("}};\n"); else  writeF("},\n"); 
      inftmp = inftmp->next;
   }

   writeF("int const Spin2[%d][%d]={\n",subproc_sq,nin + nout);

   inftmp = inf;
   for (nsub = 1; nsub <= subproc_sq; nsub++)
   {  writeF("{");
      for (i = 0; i < nin + nout; i++)
      { int pos;
        if(i) writeF(",");  
	locateinbase(inftmp->p_name[i], &pos);
         writeF("%d",prtclbase1[pos].spin);	
      }
      if (nsub== subproc_sq) writeF("}};\n"); else  writeF("},\n"); 
      inftmp = inftmp->next;
   }

   writeF("int const Color[%d][%d]={\n",subproc_sq,nin + nout);

   inftmp = inf;
   for (nsub = 1; nsub <= subproc_sq; nsub++)
   {  writeF("{");
      for (i = 0; i < nin + nout; i++)
      { int pos;
        if(i) writeF(",");  
	locateinbase(inftmp->p_name[i], &pos);
         writeF("%d",prtclbase1[pos].cdim);	
      }
      if (nsub== subproc_sq) writeF("}};\n"); else  writeF("},\n"); 
      inftmp = inftmp->next;
   }

   writeF("int const Neutral[%d][%d]={\n",subproc_sq,nin + nout);

   inftmp = inf;
   for (nsub = 1; nsub <= subproc_sq; nsub++)
   {  writeF("{");
      for (i = 0; i < nin + nout; i++)
      { int pos;
        if(i) writeF(",");  
	locateinbase(inftmp->p_name[i], &pos);
	if(pos==prtclbase1[pos].anti)
         writeF("1");	else  writeF("0"); 
      }
      if (nsub== subproc_sq) writeF("}};\n"); else  writeF("},\n"); 
      inftmp = inftmp->next;
   }
   
   writeF("int const NDF[%d][%d]={\n",subproc_sq,nin + nout);

   inftmp = inf;
   for (nsub = 1; nsub <= subproc_sq; nsub++)
   {  writeF("{");
      for (i = 0; i < nin + nout; i++)
      { int pos;
        if(i) writeF(",");  
	locateinbase(inftmp->p_name[i], &pos);
	int ndf=fabs(prtclbase1[pos].cdim);
	int spin2=prtclbase1[pos].spin;
	if(strcmp(prtclbase1[pos].massidnt,"0")==0 && spin2==2    ) ndf*=2;
	else if(strchr("LR",prtclbase1[pos].hlp)==0)  ndf*=(spin2+1);
         writeF("%d",ndf);
      }
      if (nsub== subproc_sq) writeF("}};\n"); else  writeF("},\n"); 
      inftmp = inftmp->next;
   }
    
   writeF("if(nsub<0 ||nsub>%d||nprtcl<0||nprtcl>%d) return 0;\n",
   subproc_sq,nin + nout);
   writeF("if(spin2) *spin2=Spin2[nsub-1][nprtcl-1];\n");
   writeF("if(color) *color=Color[nsub-1][nprtcl-1];\n");
   writeF("if(neutral) *neutral=Neutral[nsub-1][nprtcl-1];\n");
   writeF("if(ndf) *ndf=NDF[nsub-1][nprtcl-1];\n");
   writeF("return pcode[nsub-1][nprtcl-1];\n}\n");   
}

static void  make_den_info(void)
{  int nden_s,nden_t,nden_0;

   writeF("\n char * den_info_ext(int nsub,int n, int * mass, int * width, int*pnum)\n{\n");
   writeF(" switch(nsub){\n");

   for (nsub = 1; nsub <= subproc_sq; nsub++)
   {  int n=1;
      marktp mem_start; 
      denlist    den_;
 
      writeF(" case %d: switch(n){",nsub);
      mark_(&mem_start);
      denominatorStatistic(nsub, &nden_s, &nden_t, &nden_0, &den_, NULL); 
      for(n=1 ;den_;den_ = den_->next,n++)
      { int m=0;
          writeF("\n    case %d: *mass=%d; *width=%d;  if(pnum) *pnum=%d; return \"",
          n, vararr[den_->mass].num, vararr[den_->width].num,den_->pnum);
         while(den_->momStr[m]) fprintf(outFile,"\\%o",den_->momStr[m++]);
         fprintf(outFile,"\";");    
      }  
      writeF("\n    default:*mass=0; *width=0; if(pnum) *pnum=0; return NULL;\n                  }\n");
 
      release_(&mem_start); 
   }
   writeF("   default: *mass=0; *width=0; return NULL;\n            }\n}\n");
}


static void  make_infbasis(void)
{
   int    i,j;
   int pcolor[MAXINOUT];
   int *cb_pow=malloc(subproc_sq*sizeof(int));
   int *cb_nc=malloc(subproc_sq*sizeof(int)); 
   
   for (nsub = 1, inftmp = inf; nsub <= subproc_sq; nsub++)
   {  

      for (i = 0; i < nin + nout; i++)
      {  int l;
         locateinbase(inftmp->p_name[i], &l);
         pcolor[i]=prtclbase[l-1].cdim;
         if(i<nin && pcolor[i]!=1 && pcolor[i]!=8) pcolor[i]*=-1; 
      }

      cBasisPower=infCbases(nin+nout,pcolor,&nC,&cChains);    
      cb_pow[nsub-1]=cBasisPower;
      cb_nc[nsub-1]=nC;
      if(nC*cBasisPower)
      {  writeF("\n static int cwb_%d[%d]=\n       {\n",nsub, 4*nC*cBasisPower);
         for(i=0;i<cBasisPower;i++)
         {  writeF("       ");
            for(j=0; j<nC; j++)   
            {
              writeF(" %d,%d,%d,%d ", cChains[4*(i*nC+j)],
              cChains[4*(i*nC+j)+1],cChains[4*(i*nC+j)+2],cChains[4*(i*nC+j)+3]);
              if(i==cBasisPower-1 && j==nC-1) writeF("\n       };"); else  writeF(",");
            }
            writeF("\n");
         }      
      }else  writeF("\n#define  cwb_%d NULL\n",nsub); 
          
      inftmp = inftmp->next;
      
   }
   
   writeF("colorBasis  cb_ext[%d]={\n",subproc_sq);
   for(i=0;i<subproc_sq;i++)
   { writeF(" { %d, %d,  cwb_%d}",cb_pow[i],cb_nc[i],i+1);
     if(i!=subproc_sq-1) writeF(",");
     writeF("\n");
   } 
   writeF("};\n");
   free(cb_pow), free(cb_nc);   
}

static int gt4(int*c1,int*c2)
{ int i;
  for(i=0;i<4;i++) if( abs(c1[i])> abs(c2[i])) return 1; else if( c1[i] < c2[i]) return 0;
  return 0;
}
    

static int* makeCperm(int ntot,int*np,int*perm)
{ int i,bPow,nc; 
  int pc[MAXINOUT];
  int *chains; 
  int * res;
  
  for(i=0;i<ntot;i++) 
  { pc[i]=prtclbase[np[i]-1].cdim;
    if(i<nin && pc[i]!=1 && pc[i]!=8) pc[i]*=-1;
  }     
  
  bPow=infCbases(ntot,pc,&nC,&chains);
//printf(" Ibasis: "); for(i=0;i<bPow*nC*4;i++) printf("%d ",chains[i]);  printf("\n");
  res=malloc(sizeof(int)*(bPow+1));
  res[0]=bPow;
  if(bPow>0)
  {  int*chains2=malloc(sizeof(int)*nC*bPow*4);
     memcpy(chains2,chains,nC*bPow*4*sizeof(int));
  {  int j,k,l; 

     for(k=0;k<4*nC*bPow;k++) if(k%4 && chains2[k]) chains2[k]=perm[chains2[k]-1]+1;
//printf(" Pbasis: "); for(i=0;i<bPow*nC*4;i++) printf("%d ",chains2[i]); printf("\n");

                        
     for(l=0;l<bPow*nC;l++)
     {                      
        int * c=chains2+4*l;
        if(abs(c[0])==3)
        { int b;
             
          if(c[1]>c[2]) { b=c[1];c[1]=c[2];c[2]=b;}
          if(c[2]>c[3]) { b=c[2];c[2]=c[3];c[3]=b;}
          if(c[1]>c[2]) { b=c[1];c[1]=c[2];c[2]=b;}
          if(c[2]>c[3]) { b=c[2];c[2]=c[3];c[3]=b;}
        }
     }


     for(l=0;l<bPow;l++)
     {  int l2;
        int * c=chains2+4*l*nC;                     
        k=0; 
        while(k<nC-1) if(gt4(c+4*k,c+4*k+4) )
           {  int buff[4];
              memcpy(buff,   c+4*k,  4*sizeof(int));                           
              memcpy(c+4*k,  c+4*k+4,4*sizeof(int));                           
              memcpy(c+4*k+4,buff,   4*sizeof(int));                           
              if(k>0) k--;else k++;                                            
           } else k++;                                                         

        for(l2=0;l2<bPow;l2++) if(memcmp(chains+4*nC*l2,c,4*nC*sizeof(int))==0)
        { 
           res[l+1]=l2+1;
           break;          
        }                  
                           
        if(l2==bPow) fprintf(stderr,"Can not construct permutation\n");
     }
//printf(" Sbasis: "); for(i=0;i<bPow*nC*4;i++) printf("%d ",chains2[i]); printf("\n");     
                                                                            
  }         
     free(chains2);
  }
 return res;
} 

static int** cPerm=NULL;
static int  ncPerm=0; 

static void permRec(int ntot,int* np, int i0, int *pused, int *perm)
{  
   int i,j;

   if(i0==ntot)
   { 
     for(j=0;j<ntot;j++) if(perm[j]!=j)
     { int i; 
       writeF(",\n{");
       for(i=0;i<ntot-1;i++)writeF("%d,",perm[i]+1);
       writeF("%d}",perm[i]+1);
       
       cPerm[ncPerm++]=makeCperm(ntot,np,perm);
       
       break; 
     } 
     return;
   }

   for(j=0;j<ntot;j++) if(np[j]==np[i0] && !pused[j] )
   { perm[i0]=j;
     pused[j]=1;
     permRec(ntot,np,i0+1,pused,perm);
     pused[j]=0;
   }
}


static void  make_perm(void)
{
   int  i,j,pos=0,cbPowMax;
   int * simMap=malloc(sizeof(int)*subproc_sq); 
   writeF("static  int permMap[%d][2]={\n",subproc_sq);
   

   for (nsub = 1, inftmp = inf; nsub <= subproc_sq; nsub++,inftmp = inftmp->next)
   {  
      int nP=1;
      int k=1;
      for (i = nin+1; i < nin + nout; i++)
      {   
         if(strcmp(inftmp->p_name[i],inftmp->p_name[i-1])==0){k++;nP*=k;} else k=1;  
      }
      simMap[nsub-1]=nP-1;
      writeF(" {%d,%d}", pos,nP-1);
      if(nP>1) pos+=nP-1; 
      if(nsub == subproc_sq) writeF("\n};\n"); else writeF(",");      
   }
   
   writeF("static int permP[%d][%d]={\n",pos,nin+nout);
//   writeF("{ "); for(i=0;i<nout-1;i++) writeF("0,"); writeF("0}");   

   cPerm=malloc(sizeof(int*)*pos);
   ncPerm=0;
   int first=1;   
   for (nsub = 1, inftmp = inf; nsub <= subproc_sq; nsub++,inftmp = inftmp->next)
   {  
      if(simMap[nsub-1])
      {  int used[MAXINOUT],buff[MAXINOUT],np[MAXINOUT];
         if(first)  first=0; else  writeF(",");
         writeF("\n // "); for(i=0;i<nin;i++)        writeF("%s ", inftmp->p_name[i]);
         writeF("-> "); for(i=nin;i<nin+nout;i++) writeF("%s ", inftmp->p_name[i]);
              
         for(i=0;i<nin+nout;i++) locateinbase(inftmp->p_name[i],np+i);
         for(i=0;i<nin;i++) { used[i]=1;buff[i]=i;}  for(i=nin;i<nin+nout;i++)used[i]=0;
         permRec(nin+nout,np,nin,used,buff);  
      }
   }
   writeF("\n};\n");

   for(i=0,cbPowMax=0 ;i<ncPerm;i++) if(cbPowMax<cPerm[i][0])cbPowMax=cPerm[i][0];
   writeF("static int permC[%d][%d]={\n",pos,cbPowMax);

   for(i=0;i<ncPerm;i++)
   {
      writeF("{"); 
        for(j=1;j<cPerm[i][0];j++) writeF("%d,",cPerm[i][j]);
        if(j==cPerm[i][0])  writeF("%d",cPerm[i][j]);
        writeF("}");
      if(i<ncPerm-1) writeF(",");
      writeF("\n");
   }
   writeF("};\n");
   
//      printf("bp=%d :",cPerm[i][0]);
//    for(j=1;j<=cPerm[i][0];j++) printf(" %d",cPerm[i][j]);  
//      printf("\n");
   
   
   for(i=0;i<ncPerm;i++) free(cPerm[i]);
   free(cPerm);
    
   free(simMap);
}




static void  make_vinf(void)
{  
  writeF("char * varName_ext[%d]={\"P(cms)\"",nvars+nfunc+1);
  sortvars();
  writeF("};\n");
}

static void zeroHeep(void)
{ goto_xy(1,1);print("Heep is empty!!!");inkey();
  sortie(70);
}


static int c_prog_int(void)
{
   int breaker;
   int i;
   long dfirst;
      
   if(nin+nout<=4) sumDiag=1; else sumDiag=0;   

   memerror=zeroHeep;
   mark_(&heapbeg);

   initvararray(0,'c',NULL);
  /* ======= Initialisation parth ======= */

   firstVar=nmodelvar;
   if(!strcmp( modelvars[firstVar].varname,strongconst))  firstVar--;
   prepareprocinform();
   calc_nvars_nfunc();
  /* ======= End of Initialisation ====== */

   {  outFileOpen("%sresults%cservice.c",pathtouser,f_slash); 
      labl();
      writeF("#include<math.h>\n");
      writeF("#include<complex.h>\n");                 
      writeF("#include\"num_out.h\"\n");
      writeF("#include\"num_in.h\"\n");

      writeF("double BWrange_ext=2.7;\n");
      writeF("int twidth_ext=0;\n");
      writeF("int gtwidth_ext=0;\n");
      writeF("int gswidth_ext=0;\n");
      writeF(" REAL va_ext[%d]={0};\n",nvars+nfunc+1); 
   }
   geninf("nin_ext",nin);
   geninf("nout_ext",nout);
   geninf("nprc_ext",subproc_sq);
   make_pinf();
   geninf("nvar_ext",nvars);
   geninf("nfunc_ext",nfunc);
   make_vinf();
   
   { 
      make_den_info();
      fprintf(outFile,"\nCalcHEP_interface interface_ext={ %d,\n\"%s\"\n,%d, %d, varName_ext,va_ext,"
          "%d, %d, %d, &pinf_ext, &pinfAux_ext, polarized_ext, &calcFunc_ext, &BWrange_ext,&twidth_ext,"
          "&gtwidth_ext,&gswidth_ext, &aWidth_ext, &sqme_ext,&den_info_ext,cb_ext};\n", 
      forceUG, pathtocalchep,nvars, nfunc, nin,nout,subproc_sq);

      writeF("\nCalcHEP_interface * PtrInterface_ext=&interface_ext;\n");

      outFileClose();
      outFileOpen("%sresults%csqme.c",pathtouser,f_slash); 
      labl();
      writeF("#include<stdio.h>\n");
      writeF("#include<math.h>\n");
      writeF("#include<complex.h>\n");
      writeF("#include\"num_out.h\"\n");
      writeF("#include\"num_in.h\"\n");
/*      
      writeF("extern double complex lAAhiggs(double  Q, char* hName);\n");
      writeF("extern double complex lGGhiggs(double  Q, char* hName);\n");
      writeF("extern double complex lAA5higgs(double Q, char* hName);\n");
      writeF("extern double complex lGG5higgs(double Q, char* hName);\n");
*/      
   }
   writeF("static int calcall[%d];\n",subproc_sq+1);
   {
   writeF("static int particles[%d]={0",1+nin+nout); 
   for(i=0;i<nin+nout;i++) writeF(",0");
   writeF("};\n");
   }
   writeF("extern DNN ");
   for(i=1;i<subproc_sq;i++)  writeF("S%d_ext,",i); 
   writeF("S%d_ext;\n",subproc_sq); 
   
   writeF("static  DNN * darr[%d]={",subproc_sq);
   for(i=1;i<subproc_sq;i++)  writeF("&S%d_ext,",i);
   writeF("&S%d_ext};\n",subproc_sq);
   
   
   fseek(catalog,0,SEEK_SET);
   { catrec  cr;
     ndiagrtot =0;
     while(FREAD1(cr,catalog)) if(cr.status==1) ndiagrtot++; 
   }
   writesubroutineinit();
   
   {  make_infbasis();
      make_perm();
      writeF("#include\"sqme.inc\"\n");
      outFileClose();
   }
   diagrcount = 0;
   inftmp = inf;
   init_stat();
   for (nsub = 1,dfirst=1; nsub <= subproc_sq; nsub++)
   {  int colors[MAXINOUT];

      if (inftmp->tot != 0)   /*  this subprocess IN archive  */
      {

         for(i=0;i<nin+nout;i++) 
         {  int l;
            locateinbase(inftmp->p_name[i], &l);
            colors[i]=prtclbase[l-1].cdim;
         }
         for(i=0;i<nin; i++) if(abs(colors[i])==3 || abs(colors[i])==6) colors[i]*=-1; 
         
         if(noCChain) for(i=0;i<nin+nout; i++) colors[i]=1; 
         
         cBasisPower=infCbases(nin+nout,colors,&nC,&cChains);
         if(cBasisPower>0)
         { 
            cCoefN=malloc(cBasisPower*sizeof(long));
            cCoefD=malloc(cBasisPower*sizeof(long));
         }
         writesubprocess(nsub,dfirst,inftmp->tot, &breaker);
         dfirst+=inftmp->tot;
         if (breaker) goto exi;
      
         if(cBasisPower)
         {
            if(cChains){free(cChains); cChains=NULL;} 
            free(cCoefN); free(cCoefD);
         }

      } else writesubprocess(nsub,dfirst,0, NULL);
      inftmp = inftmp->next;
   }
   
exi:
   clearstatistic();
   release_(&heapbeg);
   return !breaker;
}


int  c_prog(void)
{  
 int result;

 catalog=fopen(CATALOG_NAME,"rb"); 
 diagrq=fopen(DIAGRQ_NAME,"rb");
 menuq=fopen(MENUQ_NAME,"rb");

 result=c_prog_int();
 fclose(catalog);
 whichArchive(0,0);
 
 fclose(diagrq);
 fclose(menuq);

 return result; 
}

int vert_code(polyvars * vardef_ext)
{
   
//   outFile=stdout;
   labl();   
   writeF("#include<stdlib.h>\n");
   writeF("#include<math.h>\n");
   writeF("#include<complex.h>\n");                 
   writeF("#include\"nType.h\"\n");    
   initvararray(0,'c',vardef_ext);

   firstVar=nmodelvar;
   if(!strcmp( modelvars[firstVar].varname,strongconst))  firstVar--;

   calc_nvars_nfunc();

//   geninf(" nvar_ext",nvars);
//   geninf(" nfunc_ext",nfunc);
//   writeF(" static int nvars=%d;\n",nvars);
   
   writeF("static char*varName[%d]={\"zero\"",nvars+nfunc+1);
   sortvars();
   writeF("};\n");   

   writeF(" static double V[%d];\n",nvars+nfunc+1);
   

   writeF("static int vertexCoeff(double * coeff_out)\n{\n");

   int    l;  

        
   for(l=nCommonVars+1;l<=nmodelvar;l++)
   { 
      char *ss;
      if(vararr[l].used && ((modelvars[l].func && modelvars[l].pub==0) || modelvars[l].pwidth) )
      {  int num;
         checkNaN=0;
         if(modelvars[l].pwidth)
         {  writeF("   %s=aWidth_ext(\"%s\");\n",vararr[l].alias,
             prtclbase1[modelvars[l].pwidth].name); 
             checkNaN=1;
         } else
         {
           ss=(char *)readExpression(modelvars[l].func,rd_c,act_c,free);
/*	   writeF("   %s=%s;\n",vararr[l].alias,ss+3);*/
	   fprintf(outFile,"   %s=%s;\n",vararr[l].alias,ss+3);
	   free(ss);
	 }
	 if(checkNaN)
	 {
    	    sscanf(vararr[l].alias,"V[%d]",&num); 
            writeF("   if(!isfinite(%s)){ return %d;}\n",vararr[l].alias,num);
         }
      }
   }
   
   for(l=0;l<vardef_ext->nvar;l++) if(!strchr( vardef_ext->vars[l].name,'.'))  
         strcpy(vardef_ext->vars[l].name,vararr[vardef_ext->vars[l].num].alias);

   return nvars;
}

