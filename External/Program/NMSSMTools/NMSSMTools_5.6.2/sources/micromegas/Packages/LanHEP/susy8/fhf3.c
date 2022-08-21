#include <math.h>
#include <stdio.h>


double fhf2f_ (double *, double *, double *, double *, double *,
	       double *, double *, double *, double *, double *,
	       double *, double *, int *);


static double tbsv,musv,Mq3sv,Mu3sv,Md3sv,Ml3sv,Mr3sv,Mq2sv,Mu2sv,Md2sv,Ml2sv,Mr2sv;
static double Atsv,Absv,Alsv,MH3sv,MG1sv,MG2sv,MG3sv,MWsv,MZsv,Mtsv,Mbsv;

static double res[6];

double feynhiggs1(double tb,double mu,double Mq3,double Mu3,double Md3,double Ml3,double Mr3,
			double Mq2,double Mu2,double Md2,double Ml2,double Mr2)
	{
	tbsv=tb;
	musv=mu;
	Mq3sv=Mq3;
	Mu3sv=Mu3;
	Md3sv=Md3;
	Ml3sv=Ml3;
	Mr3sv=Mr3;
	Mq2sv=Mq2;
	Mu2sv=Mu2;
	Md2sv=Md2;
	Ml2sv=Ml2;
	Mr2sv=Mr2;
	return 0.0;
	}

double feynhiggs2(double fh1,double At,double Ab,double Al,double MH3,double MG1,double MG2,
		double MG3,double MW,double MZ,double Mt,double Mb)
	{
	FILE *f;
	char cbuf[1024];
	int rh=0,rhh=0,rc=0,ra=0;
	Atsv=At;
	Absv=Ab;
	Alsv=Al;
	MH3sv=MH3;
	MG1sv=MG1;
	MG2sv=MG2;
	MG3sv=MG3;
	MWsv=MW;
	MZsv=MZ;
	Mtsv=Mt;
	Mbsv=Mb;
	res[1]=0.0;
	unlink("feynhiggs.out");
	f=fopen("feynhiggs.in","w");
	if(f==NULL)
		{
		puts("Error: can not open 'feynhiggs.in' for writing; unable to use FeynHiggs.");
		return 0.0;
		}
	fprintf(f,"MT           %f\n",Mtsv);
	fprintf(f,"MB           %f\n",Mbsv);
	fprintf(f,"MW           %f\n",MWsv);
	fprintf(f,"MZ           %f\n",MZsv);
	fprintf(f,"M3SQ         %f\n",Mq3sv);
	fprintf(f,"M3SU         %f\n",Mu3sv);
	fprintf(f,"M3SD         %f\n",Md3sv);
	fprintf(f,"M3SL         %f\n",Ml3sv);
	fprintf(f,"M3SE         %f\n",Mr3sv);
	fprintf(f,"M2SQ         %f\n",Mq2sv);
	fprintf(f,"M2SU         %f\n",Mu2sv);
	fprintf(f,"M2SD         %f\n",Md2sv);
	fprintf(f,"M2SL         %f\n",Ml2sv);
	fprintf(f,"M2SE         %f\n",Mr2sv);
	fprintf(f,"Abs(At)      %f\n",Atsv);
	fprintf(f,"Abs(Ab)      %f\n",Absv);
	fprintf(f,"Abs(Atau)    %f\n",Alsv);
	fprintf(f,"MA0          %f\n",MH3sv);
	fprintf(f,"TB           %f\n",tbsv);
	fprintf(f,"Abs(MUE)     %f\n",musv);
	fprintf(f,"Abs(M_2)     %f\n",MG2sv);
	fprintf(f,"Abs(M_1)     %f\n",MG1sv);
	fprintf(f,"Abs(M_3)     %f\n",MG3sv);
	fprintf(f,"Qt           -1\n");
	fprintf(f,"Qb           -1\n");
	fclose(f);
	if(system("FeynHiggs feynhiggs.in > feynhiggs.out 2>/dev/null")!=0 &&
	   system("./FeynHiggs feynhiggs.in > feynhiggs.out")!=0 &&
	   system("../FeynHiggs feynhiggs.in > feynhiggs.out")!=0 &&
	   system("../../FeynHiggs feynhiggs.in > feynhiggs.out")!=0)
	   {
	   	puts("Can not find FeynHiggs ./FeynHiggs ../FeynHiggs ../../FeynHiggs");
		return 0.0;
	   }
	f=fopen("feynhiggs.out","r");
	if(f==NULL)
		{
		puts("Error: can not open 'feynhiggs.out'; something went wrong with FeynHiggs.");
		return 0.0;
		}
	while(fgets(cbuf,1020,f))
	{
	rh+=sscanf(cbuf,"| Mh0 = %lf",&(res[1]));
	rhh+=sscanf(cbuf,"| MHH = %lf",res+2);
	rc+=sscanf(cbuf,"| MHp = %lf",res+5);
	ra+=sscanf(cbuf,"| SAeff = %lf",res+3);
	}
	fclose(f);
	
	if(rh==0 || rhh==0 || rc==0 || ra==0)
		{
		puts("Error: higgs masses are not calcuated by FeynHiggs!");
		system("tail feynhiggs.out");
		res[1]=0.0;
		return 0.0;
		}
	res[4]=sqrt(1.0-res[3]*res[3]);
	return 0.0;
	}


double 
fhf2 (double tb, double mu, double mq3, double mu3, double md3,
      double at, double ab, double mg2, double mg3,
      double ma, double mt, double mb, double tp)
{
  int sw;

  if (mq3 == 0.0 && mu3 == 0.0)
    return 0.0;

  sw = (int) floor (tp + 0.5);
  if (sw < 1 || sw > 5)
    {
      puts ("fhf1: wrong selector");
      return 0;
    }
	
  at=at-mu/tb;
  ab=ab-mu*tb;
  return fhf2f_ (&tb, &mu, &mq3, &mu3, &md3, &at, &ab, &mg2, &mg3, &ma, &mt, &mb, &sw);

}

double feynhiggs(double t, double tp)
	{
	if(res[1]==0.0)
		return fhf2(tbsv,musv,Mq3sv,Mu3sv,Md3sv,Atsv,Absv,MG2sv,MG3sv,MH3sv,Mtsv,Mbsv,tp);
	return res[(int)floor(tp+0.1)];
	}


/*
   void main()
   {
   printf("%f %f %f %f %f\n",
   fhf1(30.0,1000.0,1000.0,900.0,1500.0,500.0,175.0,1.0),
   fhf1(30.0,1000.0,1000.0,900.0,1500.0,500.0,175.0,2.0),
   fhf1(30.0,1000.0,1000.0,900.0,1500.0,500.0,175.0,3.0),
   fhf1(30.0,1000.0,1000.0,900.0,1500.0,500.0,175.0,4.0),
   fhf1(30.0,1000.0,1000.0,900.0,1500.0,500.0,175.0,5.0));
   }
 */
