#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>

float isasugf_ (float *, float *, float *, float *, float *,
		float *, int *, int *);

static double m0sv=0.0, m1sv=0.0, a0sv=0.0, tbsv=0.0, mtsv=0.0, smsv=0.0, MDsv=0.0; 
static double mssmpar[25];

double slha2fhf[6]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

static void readprm(void)
{
	int blt=0,i;
	char cbuf[1000];
	FILE *f=fopen("slhaOutput","r");
	
	slha2fhf[1]=0.0;
	
	if(f==0)
	{
		puts("Can not open slhaOutput!");
		return ;
	}
	
	while(fgets(cbuf,1000,f))
	{
		int i1,i2;
		double v;
		for(i=0;cbuf[i];i++)
			if(isalpha(cbuf[i]) && isupper(cbuf[i])) cbuf[i]=tolower(cbuf[i]);
		if(strncmp(cbuf,"block",5)==0)
		{
			char blnm[20];
			sscanf(cbuf+5,"%s",blnm);
			if(strcmp(blnm,"hmix")==0) blt=1;
			else if(strcmp(blnm,"msoft")==0) blt=2;
			else if(strcmp(blnm,"mass")==0) blt=3;
			else if(strcmp(blnm,"chepxisain")==0) blt=4;
			else if(strcmp(blnm,"au")==0) blt=11;
			else if(strcmp(blnm,"ad")==0) blt=12;
			else if(strcmp(blnm,"ae")==0) blt=13;
			else if(strcmp(blnm,"alpha")==0) blt=14;
			else blt=0;
			continue;
		}
		if(cbuf[0]!=' ')
			continue;
		if(blt==0)
			continue;
		
		if(blt==14)
		{
			if(sscanf(cbuf,"%lf",&v)!=1)
			continue;
		}
		else
		{
		if(blt<10 && sscanf(cbuf,"%d%lf",&i1,&v)!=2)
			continue;
		if(blt>10 && sscanf(cbuf,"%d%d%lf",&i1,&i2,&v)!=3)
			continue;
		}
		
		switch(blt)
		{
			case 1:
				if(i1==1) mssmpar[2]=v;
				/*if(i1==4) mssmpar[3]=sqrt(v);*/
				break;
			case 3:
				if(i1==36) mssmpar[3]=v;
				if(i1==1000021) mssmpar[1]=v;
				if(i1==25) slha2fhf[1]=v;
				if(i1==35) slha2fhf[2]=v;
				if(i1==37) slha2fhf[5]=v;
				break;
			case 4:
				if(i1>0 && i1<25) mssmpar[i1]=v;
				break;
			case 11:
				if(i1==3 && i2==3) mssmpar[20]=v;
				break;
			case 12:
				if(i1==3 && i2==3) mssmpar[21]=v;
				break;
			case 13:
				if(i1==3 && i2==3) mssmpar[22]=v;
				break;
			case 14:
					slha2fhf[3]=sin(v), slha2fhf[4]=cos(v);
				break;
			case 2:
				if(i1==1) mssmpar[23]=v;
				if(i1==2) mssmpar[24]=v;
				/*if(i1==3) mssmpar[1]=v;*/
				if(i1==32) mssmpar[13]=v;
				if(i1==33) mssmpar[18]=v;
				if(i1==35) mssmpar[14]=v;
				if(i1==36) mssmpar[19]=v;
				if(i1==42) mssmpar[10]=v;
				if(i1==43) mssmpar[15]=v;
				if(i1==45) mssmpar[12]=v;
				if(i1==46) mssmpar[17]=v;
				if(i1==48) mssmpar[11]=v;
				if(i1==49) mssmpar[16]=v;
			default:
				break;
		}
	}
	
	fclose(f);
	
}



double 
sugra (double am0, double am1, double aa0, double atb,
       double aasm, double MT, double MODE, double ch)
{
/*  float m0, m1, a0, tb, sm, mt, r;*/
  int p, q, i;
  FILE *f;

  if(am0!=m0sv || am1!=m1sv || aa0!=a0sv || atb!=tbsv || 
		  aasm!=smsv || MT!=mtsv || MODE!=MDsv)
  {
	  m0sv=am0;
	  m1sv=am1;
	  a0sv=aa0;
	  tbsv=atb;
	  smsv=aasm;
	  mtsv=MT;
	  MDsv=MODE;
	  q = (int) floor (MODE + 0.5);
	  
	  for(i=0;i<25;i++) mssmpar[i]=0.0;
	  
	  f=fopen("slhaInput","w");
	  /*fprintf(f,"Block SMINPUTS\n    6   %.15E   # Mtop\n",MT);*/
	  fprintf(f,"Block MODSEL\n    1   1  # sugra\nBlock MINPAR\n");
	  fprintf(f,"    1   %.15E  #   m0\n",am0);
	  fprintf(f,"    2   %.15E  #   m12\n",am1);
	  fprintf(f,"    3   %.15E  #   tb\n",atb);
	  fprintf(f,"    4   %.15E  #   sgn(mu)\n",aasm);
	  fprintf(f,"    5   %.15E  #   A0\n",aa0);
	  fclose(f);
	  
	  if(access("slhaScript",X_OK)==0)
		  system("./slhaScript");
	  else if(access("../slhaScript",X_OK)==0)
		  system("../slhaScript");
	  else
		  puts("slhaScript is not found!");
	  if(access("slhaOutput",R_OK)==0)
		  readprm();
	  else
		  puts("No slhaOutput file is found!");
  }

  p = (int) floor (ch + 0.5);
  
  if(p<30)  
  	return mssmpar[p];
  else
  	return slha2fhf[p-30];
  
/*
  
  m0 = am0;
  m1 = am1;
  a0 = aa0;
  tb = atb;
  sm = aasm;
  mt = MT;
  p = (int) floor (ch + 0.5);
  q = (int) floor (MODE + 0.5);

  printf("Hello world! parameter is %s\n",getcwd(0,0));
  
  r = isasugf_ (&m0, &m1, &a0, &tb, &sm, &mt, &q, &p);
  return (double) r;
*/
}

float isagmbf_ (float *, float *, float *, float *, float *,
		float *, float *, int *, int *);

static double gvsv=0.0;

double 
gmsb (double am0, double am1, double aa0, double atb,
      double aasm, double xcgv, double MT, double MODE, double ch)
{
/*  float m0, m1, a0, tb, sm, mt, r, cgv;
  int p, q;*/
		  
  int p, q, i;
  FILE *f;

  if(am0!=m0sv || am1!=m1sv || aa0!=a0sv || atb!=tbsv || 
		  aasm!=smsv || MT!=mtsv || xcgv!=gvsv || MODE!=MDsv)
  {
	  m0sv=am0;
	  m1sv=am1;
	  a0sv=aa0;
	  tbsv=atb;
	  smsv=aasm;
	  gvsv=xcgv;
	  MDsv=MODE;
	  mtsv=MT;
	  q = (int) floor (MODE + 0.5);
	  
	  for(i=0;i<25;i++) mssmpar[i]=0.0;
	  
	  f=fopen("slhaInput","w");
	  /*fprintf(f,"Block SMINPUTS\n    6   %.15E   # Mtop\n",MT);*/
	  fprintf(f,"Block MODSEL\n    1   2  # gmsb\nBlock MINPAR\n");
	  fprintf(f,"    1   %.15E  #   Lambda\n",am0);
	  fprintf(f,"    2   %.15E  #   Messenger scale\n",am1);
	  fprintf(f,"    3   %.15E  #   tb\n",atb);
	  fprintf(f,"    4   %.15E  #   sgn(mu)\n",aasm);
	  fprintf(f,"    5   %.15E  #   N5 messenger index\n",aa0);
	  fprintf(f,"    6   %.15E  #   Gravitino mass factor\n",xcgv);
	  fclose(f);
	  
	  if(access("slhaScript",X_OK)==0)
		  system("./slhaScript");
	  else if(access("../slhaScript",X_OK)==0)
		  system("../slhaScript");
	  else
		  puts("slhaScript is not found!");
	  if(access("slhaOutput",R_OK)==0)
		  readprm();
	  else
		  puts("No slhaOutput file is found!");
  }

  p = (int) floor (ch + 0.5);
  
  if(p<30)
	  return mssmpar[p];
  else
  	  return slha2fhf[p-30];
	  
/*
  m0 = am0;
  m1 = am1;
  a0 = aa0;
  tb = atb;
  sm = aasm;
  mt = MT;
  cgv = xcgv;
  p = (int) floor (ch + 0.5);
  q = (int) floor (MODE + 0.5);

  r = isagmbf_ (&m0, &m1, &a0, &tb, &sm, &cgv, &mt, &q, &p);
  return (double) r;
*/
}

double 
amsb (double am0, double am1,  double atb,
      double aasm,  double MT, double MODE, double ch)
{
/*  float m0, m1, a0, tb, sm, mt, r, cgv;
  int p, q;*/
		  
  int p, q, i;
  FILE *f;

  if(am0!=m0sv || am1!=m1sv || atb!=tbsv || 
		  aasm!=smsv || MT!=mtsv || MODE!=MDsv)
  {
	  m0sv=am0;
	  m1sv=am1;
	  tbsv=atb;
	  smsv=aasm;
	  MDsv=MODE;
	  mtsv=MT;
	  q = (int) floor (MODE + 0.5);
	  
	  for(i=0;i<25;i++) mssmpar[i]=0.0;
	  
	  f=fopen("slhaInput","w");
	  /*fprintf(f,"Block SMINPUTS\n    6   %.15E   # Mtop\n",MT);*/
	  fprintf(f,"Block MODSEL\n    1   3  # gmsb\nBlock MINPAR\n");
	  fprintf(f,"    1   %.15E  #   Scalar mass\n",am0);
	  fprintf(f,"    2   %.15E  #   Gravitino mass\n",am1);
	  fprintf(f,"    3   %.15E  #   tb\n",atb);
	  fprintf(f,"    4   %.15E  #   sgn(mu)\n",aasm);
	  fclose(f);
	  
	  if(access("slhaScript",X_OK)==0)
		  system("./slhaScript");
	  else if(access("../slhaScript",X_OK)==0)
		  system("../slhaScript");
	  else
		  puts("slhaScript is not found!");
	  if(access("slhaOutput",R_OK)==0)
		  readprm();
	  else
		  puts("No slhaOutput file is found!");
  }

  p = (int) floor (ch + 0.5);
  
  
  return mssmpar[p];
}

/*
   void main(void)
   {

   gmsb(1.0e4, 1.0e5, 10.0, 10.0, 1.0, -1.0, 175.0, 1.0, 1.0);
   }

 */
