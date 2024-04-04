
#include<math.h>
#include"micromegas.h"
#include"micromegas_aux.h"

#define SQR(x) (x)*(x)

// Routines to derive dilepton and dijet constraints on a new Z boson from LHC searches

static int Zprime_dilepton(double MZp, char * Zpname, txtList L)
{
// ATLAS Z' limits with 3.2 fb-1 at sqrt(s) = 13 TeV

// Two parts are purely UMSSM-dependent, the rest of the code can be used to apply in any model with a Z'
// Nevertheless caution between the difference in the computation of cross section with hCollider+specific PDF set and the theoretical curve derived by the experimental collaboration, a modified version of the second UMSSM-dependent part might be needed    
    
printf("\n-I- Using ATLAS Z' limits from dilepton final state with 3.2 fb-1 at sqrt(s) = 13 TeV\n");
printf("From Phys. Lett. B761 (2016) 372-392 (http://arxiv.org/abs/1607.03669)\n");
printf("----------------------------------\n");

// First UMSSM-dependent step when MZ2 is large and tE6 not in [-0.6864;0.6595] to avoid useless tests of very heavy Zprime :
// Derived lower bound on the Z2 mass and the associated tE6 angle for which the ATLAS limits exclude a UMSSM point in the case Br(Z2 -> SM)=1.
double AngleE6[44]={+1.5708,+1.5259,+1.4810,+1.4362,+1.3913,+1.3464,+1.3015,+1.2566,+1.2118,+1.1669,+1.1220,+1.0771,+1.0322,+0.9874,+0.9425,+0.8976,+0.8527,+0.8078,+0.7630,+0.7181,+0.6595,0.,-0.6864,-0.7181,-0.7630,-0.8078,-0.8527,-0.8976,-0.9425,-0.9874,-1.0322,-1.0771,-1.1220,-1.1669,-1.2118,-1.2566,-1.2925,-1.3015,-1.3464,-1.3913,-1.4362,-1.4810,-1.5259,-1.5708};
double ZpE6[44]=
{2750.4,2759.9,2769.1,2777.6,2785.7,2793.3,2800.5,2807.3,2813.9,2820.3,2826.6,2833.1,2839.7,2846.6,2853.8,2861.3,2869.2,2877.4,2885.9,2894.8,2906.9,3050.,2903.5,2891.4,2873.9,2856.0,2837.7,2819.0,2800.1,2781.4,2763.9,2747.7,2734.0,2723.1,2715.5,2711.2,2710.1,2710.2,2712.2,2716.9,2723.6,2731.8,2740.9,2750.4};
double tE6,Mlim=0;
int E6model=findVal("tE6",&tE6);
if(E6model!=0)
{
if(MZp>4000.){printf("POINT NOT TESTED - this subroutine cannot test %s mass above 4000 GeV.\n",Zpname);return 2;}
else printf("Pure UMSSM-dependent test skipped\n");
}
else
{
if(tE6>=0.6595 || tE6<=-0.6864) {Mlim=polint3(tE6,44,AngleE6,ZpE6);printf("Interpolated MZ2 limit in the case Br(Z2 -> SM)=1, for tE6=%+5.4f rad : %5.3f GeV\n",tE6,Mlim);}

if(MZp > Mlim)
{
printf("-> below the %s mass chosen\n",Zpname);
printf("==================================\n");
return 0;
}
}

if(Mlim==0.) printf("-> First derived lower bound not applicable here\n");
else printf("-> above the %s mass chosen\n",Zpname);
printf("Check whether this point is safe after including BSM branching fractions\n");
double BrZpll=(findBr(L,"e,E")+findBr(L,"m,M"))/2.;
printf("leptonic branching of %s : %4.2f %%\n",Zpname,200*BrZpll);


double cs,Pcm=6500,pTmin=0,Qren=pTmin,Qfact=MZp;
int nf=3;
printf("Computation of pp -> %s +jet(pt>%.2E GeV)  at %.2E GeV :\n",Zpname,pTmin,Pcm); 

  char oldPDF[50];
  strcpy(oldPDF,pdfName);
  setPDT("cteq6l1");  
cs=hCollider(Pcm,1,nf,Qren, Qfact,Zpname,NULL,pTmin,0);
printf("cs(pp->%s)=%.2E[fb]\n",Zpname,cs*pow(10,3));
printf("cs(pp->%s)*Br(%s->ll)=%.3E[fb]\n",Zpname,Zpname,cs*pow(10,3)*BrZpll);
restorePDF(oldPDF);

// Second UMSSM-dependent step to be compatible with the sigma*Br(Z' -> l+l-)_th shown in http://arxiv.org/abs/1607.03669 :
// Data from the theoretical curve sigma*Br(Z' -> l+l-) shown in http://arxiv.org/abs/1607.03669 (see Fig 2a) for the model psi :
double ZpPsi[43]={.5,.55,.6,.65,.7,.75,.8,.85,.9,.95,1.,1.05,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2.,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.,3.3,3.5,3.75};//Z_psi mass in TeV
double sigBPsiATLAS[43]={2333.,1645.,1172.,875.7,648.2,494.,391.4,301.2,243.3,194.7,155.7,128.3,86.18,72.37,60.19,51.04,42.87,36.7,30.53,26.14,22.38,18.98,16.57,14.19,12.39,10.71,9.351,8.085,7.128,6.163,4.743,3.651,2.837,2.205,1.730,1.358,1.076,.853,.669,.541,.272,.177,.104};//sigma*Br(Z' -> l+l-)_th
// Array of ratios between sigma*Br(Z' -> l+l-)_th and the sigma*Br(Z2 -> l+l-) obtained with hCollider using cteq6l1 for each Z' mass given in the array ZpPsi, all this in the UMSSM case Br(Z2 -> SM)=1 :
double RatiosigBPsi[43]={1.436,1.442,1.421,1.439,1.414,1.404,1.433,1.398,1.421,1.412,1.394,1.405,1.385,1.394,1.387,1.398,1.387,1.402,1.370,1.370,1.369,1.352,1.368,1.354,1.367,1.361,1.364,1.352,1.367,1.351,1.349,1.345,1.342,1.338,1.336,1.335,1.338,1.340,1.324,1.343,1.326,1.341,1.360};
// Get the rescaled sigma*Br(Z' -> l+l-) using a cubic interpolation of the array of ratios above : 
double sBth=cs*pow(10,3)*BrZpll*polint3(MZp*pow(10,-3),43,ZpPsi,RatiosigBPsi);
printf("Rescaled cs(pp->%s)*Br(%s->ll) with ATLAS data for (tE6,MZ2)=(%+5.4f rad, %5.3f GeV) = %.3E[fb]\n",Zpname,Zpname,tE6,MZp,sBth);

// Data from the observed limit curve sigma*Br(Z' -> l+l-) shown in http://arxiv.org/abs/1607.03669 (see Fig 2a) :
double Zpmass[36]={.5,.55,.6,.7,.75,.8,.85,.9,.95,1.,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2.,2.2,2.25,2.8,3.,3.3,4.};//Z' mass in TeV
double sigBexp[36]={13.27,8.437,12.65,4.635,4.788,4.275,2.466,3.695,4.416,2.852,2.311,2.63,2.761,2.806,3.143,2.899,3.636,3.353,2.899,2.03,1.645,1.422,1.333,1.29,1.493,1.542,1.469,1.29,1.098,1.029,0.98,0.949,0.964,0.98,1.012,1.171};//sigma*Br(Z' -> l+l-)_exp
double sBexp=polint3(MZp*pow(10,-3),36,Zpmass,sigBexp);
printf("Observed limit cs(pp->Z')*Br(Z'->ll) for MZprime=%5.3f GeV : %.3E[fb]\n",MZp,sBexp);


if(sBexp >= sBth)
{
printf("==================================\n");
}
else
{
printf("%s too light, parameter point is excluded by ATLAS dilepton searches\n",Zpname);
printf("==================================\n");
return 1;
}

return 0;
}




static int Zprime_dijet(double MZp, char * Zpname, txtList L)
{

// ATLAS and CMS dijet searches at sqrt(s) = 8 and 13 TeV
// recasted for a generic Zprime model by M. Fairbairn, J. Heal, F. Kahlhoefer and P. Tunney, "Constraints on Z' models from LHC dijet searches", http://arxiv.org/abs/1605.07940

printf("\n-II- Using recasted LHC dijet searches from M. Fairbairn et al., 'Constraints on Z' models from LHC dijet searches', http://arxiv.org/abs/1605.07940.\n");
printf("Done through a combination of ATLAS dijet searches at\n");
printf("sqrt(s) = 8 TeV : Phys. Rev. D91 (2015) 052007 (http://arxiv.org/abs/1407.1376)\n");
printf("sqrt(s) = 13 TeV : Phys. Lett. B754 (2016) 302-322 (http://arxiv.org/abs/1512.01530)\n");
printf("and CMS dijet searches at\n");
printf("sqrt(s) = 8 TeV : Phys. Rev. D91 (2015) 052009 (http://arxiv.org/abs/1501.04198) and Phys. Rev. Lett. 117 (2016) 031802 (http://arxiv.org/abs/1604.08907)\n");
printf("sqrt(s) = 13 TeV : Phys. Rev. Lett. 116 (2016) 071801 (http://arxiv.org/abs/1512.01224)\n");
printf("----------------------------------\n");

double gVZp_u=0,gAZp_u=0,gVZp_d=0,gAZp_d=0,gVZp_c=0,gAZp_c=0,gVZp_s=0,gAZp_s=0,gVZp_b=0,gAZp_b=0;
int Zpqqbar=findVal("gVZp_u",&gVZp_u)+findVal("gAZp_u",&gAZp_u)+findVal("gVZp_d",&gVZp_d)+findVal("gAZp_d",&gAZp_d)+findVal("gVZp_c",&gVZp_c)+findVal("gAZp_c",&gAZp_c)+findVal("gVZp_s",&gVZp_s)+findVal("gAZp_s",&gAZp_s)+findVal("gVZp_b",&gVZp_b)+findVal("gAZp_b",&gAZp_b);
if(Zpqqbar==20){printf("POINT NOT TESTED - there is no definition of any function corresponding to vectorial or axial Zprime-q-qbar coupling in this model in the form 'gVZp_q' or 'gAZp_q', with q=u,d,c,s or b.\n");return 2;}

double Gamma_q_inv=findBr(L,"u,U")+findBr(L,"d,D")+findBr(L,"c,C")+findBr(L,"s,S")+findBr(L,"b,B")+findBr(L,"b,B");
char*  name[8]={CDM1,CDM2,pdg2name(12),pdg2name(14),pdg2name(16),pdg2name(9912),pdg2name(9914),pdg2name(9916)};
char*  aname[8]={aCDM1,aCDM2,pdg2name(-12),pdg2name(-14),pdg2name(-16),pdg2name(-9912),pdg2name(-9914),pdg2name(-9916)};

for(int i=0;i<8;i++) if(name[i])
  { char channel[40];
    sprintf(channel,"%s %s",name[i],aname[i]);
    Gamma_q_inv+=findBr(L,channel);
  }

double Ratio=Gamma_q_inv/MZp;
printf("Ratio width/mass[%s] = %.2e, where only quarks (except top) and invisible particles are included in the width calculation\n",Zpname,Ratio);

printf("Compute jcalc=sum_q (gq^2_V + gq^2_A) * Br(Zprime->qqbar) for this point and compare the result to the interpolated upper limit at 95%% CL :\n");

double jcalc=(SQR(gVZp_u)+SQR(gAZp_u))*findBr(L,"u,U")+(SQR(gVZp_d)+SQR(gAZp_d))*findBr(L,"d,D")+(SQR(gVZp_c)+SQR(gAZp_c))*findBr(L,"c,C")+(SQR(gVZp_s)+SQR(gAZp_s))*findBr(L,"s,S")+(SQR(gVZp_b)+SQR(gAZp_b))*findBr(L,"b,B");

printf("sum_q (gq^2_V + gq^2_A)*Br(%s->qqbar) = %.3e\n",Zpname,jcalc);

int row, i, j, l, m;
double dijetgrid[4260][3],Zprimearray[71],Ratioarray[60],j95grid[71][60];
double ULj95=0;

char fname[300];
sprintf(fname, "%s/sources/data/dijetgrid.dat",micrO);
FILE *f=fopen(fname, "r");
if(!f) return 1;
for(i=0;i<5;i++) fscanf(f," %*[^\n]");
for(row=0; row<4260; row++)fscanf(f, "%lf%lf%lf", &dijetgrid[row][0], &dijetgrid[row][1], &dijetgrid[row][2]);

fclose(f);

for(i=0;i<71;i++) Zprimearray[i]=dijetgrid[60*i+1][0]; 
for(i=0;i<60;i++) Ratioarray[i]=dijetgrid[i][1];  

for(i=0;i<71;i++) { for(j=0;j<60;j++) j95grid[i][j]=dijetgrid[60*i+j][2]; }

if(MZp>Zprimearray[70]){printf("POINT NOT TESTED - this subroutine cannot test %s mass above %4.0f GeV.\n",Zpname,Zprimearray[70]);return 2;}
if(Ratio<0){printf("ERROR - this subroutine cannot test a ratio width/mass[%s] below 0.\n",Zpname);return 1;}
  
for(i=0;i<71;i++){ if(Zprimearray[i+1]>=MZp) break; }

for(j=0;j<60;j++){ if(Ratioarray[j+1]>=Ratio) break; }

// Computation of the Inverse-Distance-Weighted interpolation function

// Calculating the normalisation
double norm=0;
for(l=i;l<=i+1;l++){ for(m=j;m<=j+1;m++) norm=norm+1/sqrt(pow(Zprimearray[l]-MZp,2)+pow(Ratioarray[m]-Ratio,2)); }

for(l=i;l<=i+1;l++){ for(m=j;m<=j+1;m++) ULj95=ULj95+j95grid[l][m]/sqrt(pow(Zprimearray[l]-MZp,2)+pow(Ratioarray[m]-Ratio,2)); }
ULj95=ULj95/norm;

printf("The interpolated upper limit at 95%% CL on j=gq^2*Br(Z'->jj) for MZp=%5.3f GeV and ratio width/mass[Z'] = %.2e is : %.3e\n",MZp,Ratio,ULj95);

if(ULj95 >= jcalc)
{
printf("==================================\n");
}
else
{
printf("Parameter point is excluded at 95%% confidence level from LHC dijet searches\n");
printf("==================================\n");
return 1;
}
return 0;
  
}







int Zprimelimits(void)
{

if(!pdg2name(32)){printf("ERROR - there is no Zprime with PDG code = 32 in this model.\n");return 1;}

else
{
double MZp=pMass(pdg2name(32));
char * Zpname=pdg2name(32);
printf("%s mass : %5.3f GeV\n",Zpname,MZp);
printf("==================================\n");
printf("==== Limits on the %s boson : ====\n",Zpname);
printf("==================================\n");
if(MZp<500.) {printf("ERROR - this routine cannot test %s mass below 500 GeV.\n",Zpname);return 1;}

else
{
txtList L;
double width=pWidth(Zpname,&L);
printf("\n%s :   total width=%.2E[GeV]\n",Zpname,width);
if(width/MZp>.3) {printf("POINT NOT TESTED - this routine is only valid in the narrow-width approximation.\n");return 2;}

int flag_dilepton=Zprime_dilepton(MZp,Zpname,L);
if(flag_dilepton!=1)
{
int flag_dijet=Zprime_dijet(MZp,Zpname,L);
if(flag_dilepton+flag_dijet==0 || flag_dilepton+flag_dijet==2) return 0;
else if(flag_dilepton+flag_dijet==1 || flag_dilepton+flag_dijet==3) return 1;
//if(flag_dilepton+flag_dijet==4)
else {printf("POINT NOT TESTED - limits on the Zprime cannot apply to this scenario.\n");return 2;}
}
else return 1;
}

}
}

extern int zprimelimits_(void); //Fortran 
int  zprimelimits_(void) { return  Zprimelimits();}
