/*
 Copyright (C) 1997, 1999  Alexander Pukhov 
*/

#include"../../include/nType.h"

#include"q_kin.h"
#include"decay.h"
#include"4_vector.h"

#include<math.h>
#include "interface.h"
#include"cut.h"
#include"kinaux.h"
#include"strfun.h"
#include"regfunal.h"
#include"subproc.h"
#include"q_kin.h"
#include"tools.h"
#include"drandXX.h"

//#define DEBUG

static void LorRot(REAL rapidity, int ntot,REAL*pvect)
{
  REAL sh,ch,ee, pp;
  int nv;
  
  if(rapidity)
  {   sh = Sinh(rapidity);
      ch = Sqrt(sh * sh + 1);
      for(nv=4*ntot-4; nv>=0; nv -=4)
      {  ee = pvect[nv];
         pp = pvect[nv+3];
         pvect[nv] = ee * ch + pp * sh;
         pvect[nv+3] = ee * sh + pp * ch;
      }
  }
} 

static void AzimuthRot(REAL fi, int ntot,REAL*pvect)
{
  REAL sn = Sin(fi);
  REAL cn = Cos(fi);
  int i;
    
  for(i=0;i<ntot;i++)
  {  int nx=4*i+1,ny=4*i+2;
     REAL X,Y;
     X=pvect[nx];
     Y=pvect[ny];
     pvect[nx]= X*cn+Y*sn;
     pvect[ny]=-X*sn+Y*cn;
  }
} 


static void Rot3D(double  cn0 , double fi, double psi, int ntot, REAL *pvect)
{
   REAL cn=cn0, sn=Sqrt(fabs(1-cn*cn));
   REAL fi_1=fi;
   REAL E[3]={ cn, sn*Sin(fi_1), sn*Cos(fi_1)};
   REAL EPS[3][3]= { {  0,  -E[2], E[1]},
                     { E[2],   0 ,-E[0]},
                     {-E[1], E[0],  0  }  }; 

   REAL P[3][3]=  {  {1-E[0]*E[0],  -E[0]*E[1], -E[0]*E[2] },
                     { -E[1]*E[0], 1-E[1]*E[1], -E[1]*E[2] },
                     { -E[2]*E[0],  -E[2]*E[1],1-E[2]*E[2] } };
    
  REAL ROT[3][3],cn2,sn2;
  double x1,x2,x3,f3; 
  int i,j,k;


  x1=psi;
  x2=M_PI;
  for(;;)
  { x3=0.5*(x1+x2);
    f3=x3-psi-Sin(x3);
    if(f3<0) x1=x3;else x2=x3;
    if(fabs(f3)<1.E-4) break;
  }
                             
  cn2=Cos(x3),sn2=Sin(x3);
  
  for(i=0;i<3;i++) for(j=0;j<3;j++) ROT[i][j]= (cn2-1)*P[i][j] + sn2*EPS[i][j]; 
      
  for(k=0;k<ntot;k++) 
  {  REAL * X= pvect+4*k+1;
     REAL Y[3]={X[0],X[1],X[2]}; 
     for(i=0;i<3;i++)  for(j=0;j<3;j++) Y[i]+=ROT[i][j]*X[j];  
     for(i=0;i<3;i++) X[i]=Y[i];
  }                       

}


#define DEPTH 10

typedef  struct
{  int ncsreg[2];
   int ncscut[2];
   int tcscut[2];
   char  lvpole[PLISTLEN];
   int  itypep;
   REAL sph_we;
}  sphereStr;




static REAL  cmfun( REAL E,REAL PM1,REAL PM2)
{
  REAL e2,s,d;
  e2=E*E;
  s=PM1+PM2;
  d=PM1-PM2;  
  return Sqrt((e2-s*s)*(e2-d*d))/(2*E);     
}


static REAL sqrt_S, rapidity, Pcm,E1,E2;
static double ssmin, ssmax;

static REAL tfact0;
static int nout1;

static REAL pm[PLISTLEN];

static int nvpos0, nvposx, nvposy;

static int nvout[DEPTH][2], nvin[DEPTH];
static int  lnkbab[DEPTH], lnkleg[DEPTH];
static REAL  summas[DEPTH][2];
static int nmsreg[DEPTH][2], nmscut[DEPTH][2], nss;	
static int nsph[DEPTH];

static REAL beta[2];
static  sphereStr  sph_inf[DEPTH][10]; 


int mkmom(double*x,double*tfact,double*xp1,double*xp2,REAL*pvect)
{
  int i,k,l;
  int nx=0;
  double  fct0;

  REAL  pIn[2][4]={{0,0,0,0},{0,0,0,0}};
  REAL  pXY[2][4]={{0,1,0,0},{0,0,1,0}};
    
  REAL  xcos, xfi, parfi;
  double cosmin, cosmax, parcos;
  REAL ytilda=0; 
  
  int  i__2;
  REAL d__1;

  REAL ff, al;
  REAL  xx, bes; 
  double fct;
 
  double stilda;
  REAL rstilda,  pcm, e1,e2; 
  int  ns;

  REAL psy1, psy2, x1,x2;

  REAL  amass[DEPTH][2];

  int nvpole;
  int nvpos;

  REAL hsum[2], hdif;

  REAL fct_1__;
  int nsing;
  
  iDecay memDecay;
    
  sing_struct singar[200];

  *tfact = tfact0;
/* **   MOMENTS */
    if (nin_int == 2) 
    {   REAL y1,y2;
        double be;
        double fct1, fct2;
	if (sf_num[0] || sf_num[1])   // integration over dstilda  
	{
           if(nout_int==1){ stilda=pm[2]*pm[2]; } else 
           {
	      nsing = 0;
	      getreg_(&nsing, singar, 0., 1., nss);
	      bes = 0; if(sf_num[0]) bes+=beta[0];  if(sf_num[1]) bes+=beta[1];
	      xx = x[nx++];
              if(bes< 1 && nsing)
              {
                  if(xx>0.5)
                  { xx=2*xx-1;
                    stilda=ssmax-pow(xx,1/bes)*(ssmax-ssmin);
                    fct2=ssmax/bes*pow(1-stilda/ssmax,1-bes)*pow(1-ssmin/ssmax,bes);
                    regfct_(2,nsing,singar,ssmin,ssmax,stilda,&fct1);
                  }else
                  { xx*=2; 
                    regfun_(2,nsing,singar,ssmin,ssmax,xx,&stilda,&fct1);
                    fct2=ssmax/bes*pow(1-stilda/ssmax,1-bes); 
                  }  
                  *tfact *= 2*fct1*fct2/(fct1+fct2);
              }
              else if(nsing)
              {
                regfun_(2,nsing,singar,ssmin,ssmax,xx,&stilda,&fct1);
                *tfact *=fct1;             
              }   
              else if(bes< 1)   // assuming  that integrand has singularity (1-stilda/ssmax)^(bes-1)
              {  
                 stilda=ssmax -pow(xx,1/bes)*(ssmax-ssmin);
                 *tfact *=ssmax/bes*pow(1-stilda/ssmax,1-bes) * pow(1-ssmin/ssmax,bes);
              }               
              else 
              {
                  stilda=xx*ssmax +(1-xx)*ssmin;
                  *tfact *= ssmax-ssmin;
              }              
           }
           
           rstilda=sqrt(stilda);
           
           pcm=cmfun(rstilda,pm[0],pm[1]);
           x1=(pcm+Sqrt(pcm*pcm+pm[0]*pm[0]))/(Pcm+E1);       // cm for partons == cm for hadrons
           x2=(pcm+Sqrt(pcm*pcm+pm[1]*pm[1]))/(Pcm+E2);
           *tfact *= x1*x2/(stilda);   // DxDs  indeed is x1*x2/(2*pcm*rstilda) 
           if(sf_num[0] && sf_num[1]) 
           {
                double x1max,x2max,x1min,x2min,Y;
                x1max=x1*(E1+E2)/rstilda;  //  s(boosted system + jet) < (E1+E2)^2  
                if(x1max>1) x1max=1;
                x2max=x2*(E1+E2)/rstilda;
                if(x2max>1) x2max=1;
                x1min= x1*x2/x2max;
                x2min= x1*x2/x1max;
                Y=-log(x1*x2); 
//if(pcm==0) { printf("pcm=%E   rstilda= %E x1*x2=%E x1min*x2min=%E  Y=%E \n",(double)pcm, (double) (x1*x2), (double)x1min,(double) x2min, (double)Y );}                
		xx = x[nx++];

		if (beta[0] < 1 && beta[1] < 1) 
		{
		
		    al = beta[1]/(beta[0]+beta[1]);
		    if (xx < al) 
		    {   xx /= al;
			be=beta[0];
			y1= pow((1-xx)*pow(-log(x1max),be)+xx*pow(-log(x1min),be),1/be);
			*tfact *= pow(-log(x1min),be) - pow(-log(x1max),be);
			ytilda=y1+log(x1);
			y2=Y-y1;
			if(y2<0) y2=0;	
		    } else 
		    {   xx = (1 - xx) / (1 - al);
                        be=beta[1];
                        y2= pow((1-xx)*pow(-log(x2max),be)+xx*pow(-log(x2min),be),1/be);
                        
                        *tfact *= pow(-log(x2min),be) - pow(-log(x2max),be);
                        ytilda=-y2-log(x2);
                        y1=Y-y2;
                        if(y1<0) y1=0;
		    }
		    *tfact *=(beta[0]+beta[1])*pow(divy_(y1),beta[0]-1)*pow(divy_(y2),beta[1]-1)/
		    	     (pow(y1,1-beta[0])+pow(y2,1-beta[1]));
		} 		
		else if (beta[0] < 1)   // addition factor not included in structure function   be*(1-x)^(be-1) =be*(1-exp(-y))^(be-1)=div_(y)*be*y^(be-1)
		{
		    be=beta[0];
                    y1= pow((1-xx)*pow(-log(x1max),be)+xx*pow(-log(x1min),be),1/be); 
                    ytilda=y1+log(x1);   
		    *tfact *= pow(divy_(y1),be-1)* (pow(-log(x1min),be) - pow(-log(x1max),be));                    
		} else if (beta[1] < 1) 
		{
		    be=beta[1];
                    y2= pow((1-xx)*pow(-log(x2max),be)+xx*pow(-log(x2min),be),1/be); 
                    ytilda=-y2-log(x2);   
		    *tfact *= pow(divy_(y2),be-1)* (pow(-log(x2min),be) - pow(-log(x2max),be));                    
		} else 
		{
                   ytilda= (1-xx)*log(x1max) + xx*log(x1min)-log(x1);    
                   *tfact *= log(x1max/x1min) ; 
		}

                x1*=Exp( ytilda);
                x2*=Exp(-ytilda);
                		
	    } else if (sf_num[0]) 
	    {   
	        x1*=x2;                
                if(beta[0]<1) *tfact *=beta[0]*pow(1-x1,beta[0]-1); 
	        ytilda=-log(x2);    
	    } else 
	    { 
	       x2*=x1;;
	       if(beta[1]<1) *tfact *=beta[1]*pow(1-x2,beta[1]-1);
	       ytilda=log(x1);
	    }
	} else 
	{
	    stilda = sqrt_S*sqrt_S;
	    rstilda=sqrt_S;
            pcm=Pcm;
            ytilda=0;
	}

        if(sf_num[0]==3 && x1>1-1E-9) { stilda*=(1-1E-9)/x1;  x1=1-1E-9;}
        if(sf_num[1]==3 && x2>1-1E-9) { stilda*=(1-1E-9)/x2;  x2=1-1E-9;}

        
	pIn[0][3]=pcm;
	pIn[1][3]=-pcm;
	     
    }
    
    
/* *  FILLING ZERO COMPONENTS FOR in-PARTICLES */

    for (k = 0; k < nin_int; ++k) pvFill(pm[k],pIn[k],k+1,pvect);    
	 
    if (nin_int == 2)   *tfact /= 4*pcm *rstilda; else 
    { rstilda = pm[0]; *tfact /= rstilda * 2; ytilda=0;} 
/* *    X & Y AXISES */

    pvFill(0,pXY[0],nvposx,pvect);
    pvFill(0,pXY[1],nvposy,pvect);  

    nvpos = nvpos0;
    nvpole =nvpos++;
    
    if(nout_int==1) { pvect[8]=pm[2],pvect[9]=0,pvect[10]=0,pvect[11]=0;}

/* *    MASS INTEGRATION */

    for (i = 0; i < nout1; ++i) 
    {   REAL sval= i? amass[lnkbab[i]][lnkleg[i]]: rstilda;   

       	for (k = 0; k < 2; ++k) 
	{   double  smin, smax;

	    if (kinmtc_1[i].lvout[k][1]) 
	    { double sqmass;
	       REAL  xx=x[nx++];
	        d__1=  k?  sval - amass[i][0] : sval - summas[i][1];
		         
	        smax = d__1 * d__1;
		
		d__1 = summas[i][k];
		smin = d__1 * d__1;

		if (nmscut[i][k]) rancor_(&smin, &smax,0., 1., nmscut[i][k]);
		
		if (smin >= smax)  {*tfact = 0; return 0;}
		
		if (nmsreg[i][k] )
		{
		    nsing = 0;
		    getreg_(&nsing, singar, 0., 1.,nmsreg[i][k]);
		    regfun_(2,nsing,singar,smin,smax,xx,&sqmass,&fct);
		} else
		{
		    sqmass = xx * smax + (1 - xx) * smin;
		    fct = smax - smin;
		} 
		amass[i][k] = Sqrt(sqmass);
		*tfact *= fct;
	    }
	    else amass[i][k]= summas[i][k];
	}
    }

    lvtonv(kinmtc_1[0].lvin, 0 , nvin[0],pvect); /*very stupid*/ 

        
    for (i = 0; i < nout1; ++i)  /*  MAIN CYCLE */
    {   int ns___=nsph[i]-1;
        double Emax[2];
	if (i == 0 && nin_int == 1)  xcos = 0.1 /* was fixed  0.1 */; 
                               else  xcos = x[nx++];
	al = 0;
	l = 0;
	if (i == 0 || (i == 1 && nin_int == 1)) 
	{
	    xfi = 0.1;  /* was fixed  0.1;    */
            for(;l<=ns___;l++)
	    {  al +=  sph_inf[i][l].sph_we;
	       if (xcos <= al) ns___ = l;
	    }
	    xcos = (al - xcos) / sph_inf[i][ns___].sph_we;
	} else 
	{
	    xfi = x[nx++];
            for(;l<=ns___;l++)
            {
	       al += sph_inf[i][l].sph_we;
	       if (xfi <= al) ns___ = l;
	    }
	    xfi = (al - xfi) /sph_inf[i][ns___].sph_we;
	}
	lvtonv( sph_inf[i][ns___].lvpole,nin_int, nvpole,pvect);

	decay_0(nvin[i], amass[i][0], amass[i][1], &fct0, Emax,&memDecay,pvect);
	if(fct0==0) {*tfact=0; return 0;}
	if(!isfinite(fct0))
	{ printf("decay_0: i=%d nvin[i]=%d   amass[i][0]=%E , amass[i][1]=%E\n", i, nvin[i],   (double)amass[i][0] ,(double) amass[i][1]); 
	  printf("mkmom+ infinite factor\n");
	  *tfact=0;
	  return 0;
	}  
        decay_1(nvpole, hsum, &hdif,&memDecay,pvect);
        
	cosmin = -1;
	cosmax = 1;
	nsing = 0;
	for (k = 0; k < 2; ++k) 
	{   int ncM=sph_inf[i][ns___].ncscut[k];
	    int ncT=sph_inf[i][ns___].tcscut[k];
	    d__1 = ((k << 1) - 1) / hdif;
            getreg_(&nsing,singar,hsum[k],d__1,sph_inf[i][ns___].ncsreg[k]);
            if(ncM) rancor_(&cosmin,&cosmax,hsum[k],d__1,ncM);
            if(ncT) rancor_t(&cosmax,hsum[k],d__1,Emax[k], pm[sph_inf[i][ns___].lvpole[0]-1],pcm, 
                            amass[i][k], invcut_1[ncT-1].cvmin );
	    if (cosmin >= cosmax) {*tfact = 0; return 0;}
	}
	regfun_(sph_inf[i][ns___].itypep,nsing,singar,cosmin,cosmax,xcos,&parcos,&fct);
	fct_1__ = sph_inf[i][ns___].sph_we / fct;
	parfi = (xfi * 2 - 1) * M_PI;
	decay_3(nvposy, parcos, parfi, nvout[i][0], nvout[i][1],&memDecay,pvect);
		
	i__2 = nsph[i];
	for (ns = 0; ns < i__2; ++ns)  if (ns != ns___)
	{
	    lvtonv(sph_inf[i][ns].lvpole, nin_int, nvpole,pvect);
	    decay_1(nvpole, hsum, &hdif,&memDecay,pvect);
	    decay_2(nvout[i][1], &parcos,&memDecay,pvect);
	    cosmin = -1;
	    cosmax = 1;
	    nsing = 0;
	    for (k = 0; k < 2; ++k) 
	    {   int ncM=sph_inf[i][ns].ncscut[k];
                int ncT=sph_inf[i][ns].tcscut[k];
		d__1 = ((k << 1) - 1) / hdif;
		getreg_(&nsing,singar,hsum[k],d__1,sph_inf[i][ns].ncsreg[k]);
		if(ncM) rancor_(&cosmin,&cosmax,hsum[k],d__1,ncM);
                if(ncT) rancor_t(&cosmax,hsum[k],d__1,Emax[k], pm[sph_inf[i][ns___].lvpole[0]-1],pcm, 
                            amass[i][k], invcut_1[ncT-1].cvmin );

		if (cosmin>=parcos || parcos>=cosmax){*tfact=0; return 0;}
	    }
	    regfct_(sph_inf[i][ns].itypep,nsing,singar,cosmin,cosmax,parcos, &fct);
	    fct_1__ += sph_inf[i][ns].sph_we/ fct;
	}
	*tfact = *tfact * fct0 / fct_1__;
    }
    
    
    if(nin_int==2)
    { *xp1=sf_num[0]? x1:1;
      *xp2=sf_num[1]? x2:1;
      LorRot(rapidity+ytilda,nin_int+nout_int,pvect);
      AzimuthRot(drandXX()*2*M_PI,nout_int, pvect+8);        
    } else 
    { 
        Rot3D(2*(drandXX()-0.5),drandXX()*2*M_PI,drandXX()*M_PI,nout_int,pvect+4); 
        LorRot(rapidity,nin_int+nout_int,pvect);
    }
    
    if(!isfinite(*tfact))
    { 
// fprintf(stderr,"mkmom: infinite factor\n");
       printf("mkmom!: infinite factor\n");
       *tfact=0;
      return 0;
    }

    for(i=0;i<(nin_int+nout_int)*4;i++) if(!isfinite(pvect[i])) {*tfact=0; return 0;}

    return 0;
}    


int imkmom(double P1, double P2)
{
    int i, j, k, l,ns;
    char lvbuf[PLISTLEN];
    int ndim;
    physValRec * pList;
    if(sf_num[0])  beta[0]=sf_be[0]; else beta[0]=1;
    if(sf_num[1])  beta[1]=sf_be[1]; else beta[1]=1;
    if(nin_int==2)  
    {  
       tfact0 = 2*M_PI*389379660.0;
       ndim = nout_int * 3 - 5 + ((nout_int==1)? 1:0);
       
       if (sf_num[0]) ndim++;
       if (sf_num[1]) ndim++;        
    }else { tfact0 = 2*M_PI; if(nout_int==2) ndim=1; else ndim = nout_int * 3 - 7;}

    for (i=0; i <  nin_int + nout_int; i++) pinf_int(Nsub,i+1,pm+i,NULL);

    nout1 = nout_int - 1; 
    if(nout1>DEPTH) return 0;
    nvposx = nin_int + nout_int + 1;
    nvposy = nvposx + 1;
    nvpos0 =  nvposy + 1;

/* *  NVOUT( , ) FILLING */
    for (i = 0; i < nout1; ++i)  for (k = 0; k < 2; ++k) 
    {
       if (kinmtc_1[i].lvout[k][1]) nvout[i][k] = nvpos0++;
       else                         nvout[i][k] = kinmtc_1[i].lvout[k][0];   
    }
    
    nvin[0] = nvpos0++;
    for (i = 1; i < nout1; ++i) 
    {   nvin[i]=0;
	for (j = 0; j < i; ++j) for (k = 0; k < 2; ++k) 
	{
	   if (eqvect_(kinmtc_1[i].lvin, kinmtc_1[j].lvout[k])) 
	   {
               nvin[i] = nvout[j][k];
	       lnkbab[i] = j;
	       lnkleg[i] = k;
	   }
	}
	if(!nvin[i]) { fprintf(stderr,"Error in kinematics \n"); sortie(52); }
    }


    for (i = 0; i < nout1; ++i) for (k = 0; k < 2; ++k) 
    {   
      REAL ss = 0; 
      int pn;
      
      for(j=0; (pn=kinmtc_1[i].lvout[k][j]);j++) ss += pm[pn - 1];
      summas[i][k] = ss;
    }

    ssmin=0; for(int i=nin_int;i<nin_int+nout_int;i++) ssmin+=pm[i];
    if(nin_int==1) { if(ssmin<pm[0]) ssmin=pm[0];} else { if(ssmin<pm[0]+pm[1]) ssmin=pm[0]+pm[1];}        
    ssmin*=ssmin;

    if (nin_int == 2) 
    {  REAL m1=sf_num[0]?sf_mass[0]:pm[0];
       REAL m2=sf_num[1]?sf_mass[1]:pm[1];

       incomkin(m1, m2, P1, P2,  &sqrt_S, &Pcm, &rapidity);
       E1=sqrt(Pcm*Pcm+m1*m1);
       E2=sqrt(Pcm*Pcm+m2*m2);

       double p1=0.5*( Pcm+E1 - pm[0]*pm[0]/(Pcm+E1));
       double e1=sqrt(p1*p1+pm[0]*pm[0]);
       double p2=0.5*( Pcm+E2 - pm[1]*pm[1]/(Pcm+E2));
       double e2=sqrt(p2*p2+pm[1]*pm[1]);
       double ss=pm[0]*pm[0]+pm[1]*pm[1]+2*(e1*e2+p1*p2);
       if(ss< sqrt_S*sqrt_S) ssmax=ss; else ssmax=sqrt_S*sqrt_S;
       
    } else 
    {   REAL m1=pm[0];

//        rapidity=Log((P1+Sqrt(P1*P1+m1*m1))/m1 ) ;
        rapidity=0;  
    }
    for(i = 0; i < nout1; ++i) 
    {
       nsph[i] = 0;
       for(k=0;k<2;k++)
       {
          for(ns=0; ns<10;ns++)
          {
             sph_inf[i][ns].ncsreg[k] = 0;
             sph_inf[i][ns].ncscut[k] = 0;
             sph_inf[i][ns].tcscut[k]=0;
          }
	  nmsreg[i][k] = 0;
	  nmscut[i][k] = 0;
       }
    }
    
    nss=0;
    for(l=0; invreg_1[l].lvinvr[0]; l++) 
    {   int orig=1, ll=0;
        for(;ll<l;ll++) {if( invreg_1[ll].nextrg == l+1) { orig=0; break;}}	
	if(orig) 
	{
	   sngpos_(invreg_1[l].lvinvr, &i, &k, lvbuf);
	   if (i==0)  nss = l+1; else
	   {  i--; k--;
	      if (lvbuf[0] == 0)  nmsreg[i][k] = l+1;
	      else 
	      {
		 for (ns = 0; ns <nsph[i]; ++ns) 
                 if (eqvect_(lvbuf,sph_inf[i][ns].lvpole)) 
		 {
		    sph_inf[i][ns].ncsreg[k] = l+1;
		    break;
		 }
		 if(ns==nsph[i] && ns<10)
		 {
		    nsph[i]++;
		    strcpy(sph_inf[i][ns].lvpole,lvbuf);
		    if (spole_(invreg_1[l].lvinvr))sph_inf[i][ns].itypep = -2; 
		    else sph_inf[i][ns].itypep = -1;
		         
	            sph_inf[i][ns].ncsreg[k] = l+1;
	         }
	      }
	   }
	}
    }

    for(l=0;l<nCuts;l++) if( invcut_1[l].key[0] == 'M')
    for(pList=invcut_1[l].pLists;pList;pList=pList->next)
    {   char buff[20];
        strcpy(buff,pList->pstr);
        coninv_(buff);
	sngpos_(buff, &i, &k, lvbuf);
	if (i == 0) rancor_(&ssmin, &ssmax, 0., 1., l+1);
	
	else if (lvbuf[0] == 0)  nmscut[i-1][k-1] = l+1;
	else 
	{   i--; k--;
	    for (ns = 0; ns <  nsph[i ]; ++ns) 
	    {
		if (eqvect_(lvbuf, sph_inf[i][ns].lvpole)) 
		{
		    sph_inf[i][ns].ncscut[k]  = l+1;
		    break;
		}
	    }
	    if(ns==nsph[i] && ns <10 )
	    {
	       nsph[i]++;
	       strcpy(sph_inf[i][ns].lvpole,lvbuf);
	       sph_inf[i][ns].itypep = 2;
               sph_inf[i][ns].ncscut[k] = l+1;
	    }
	}
    }
    if(nin_int==2 && ssmin>=ssmax) return 0;

    if(nin_int==2) for(l=0;l<nCuts;l++)
    {  int m;
       invcut_ tc=invcut_1[l];
/* 
       if( tc.key[0]=='T' && tc.key[1]!='^' && tc.minon 
                          && tc.pLists      && tc.pLists->pstr[1]==0 )  
       for(pList=invcut_1[l].pLists;pList;pList=pList->next)for(m=1;m<=2;m++)
       {  char str[4];
          strcpy(str+1,pList->pstr);
          str[0]=m;
          sngpos_(str, &i, &k, lvbuf);
	  {  i--; k--;
	     for (ns = 0; ns <  nsph[i ]; ++ns) 
	     {
		if(eqvect_(lvbuf, sph_inf[i][ns].lvpole)
                   &&strcmp(kinmtc_1[i].lvout[k],str+1)==0)
		{
		    sph_inf[i][ns].tcscut[k]  = l+1;
		    break;
		}
	     }
	  }
       }
*/ 
    }
        
    for(i = 0; i < nout1; ++i) 
    {
	if (nsph[i] == 0) 
	{
	   nsph[i] = 1;
           sph_inf[i][0].lvpole[0] = (i == 0 && nin_int == 1 )?  nvposx:1;
	   sph_inf[i][0].lvpole[1] = 0;
	   sph_inf[i][0].itypep = 1;
	} else 
	{
	   ns = nsph[i]-1;
	   if(sph_inf[i][ns].ncsreg[0] || sph_inf[i][ns].ncsreg[1])  
	   sph_inf[i][0].itypep *=-1;
	}
    }

    for(i = 0; i<nout1; ++i) 
    {   REAL wesum = 0;
	for (ns = 0; ns<nsph[i]; ++ns) 
	{
	    int   nwe = 0;
	    for (k = 0; k < 2; ++k) 
	    for (l = sph_inf[i][ns].ncsreg[k];l;l=invreg_1[l-1].nextrg) ++nwe;
		    
            if (sph_inf[i][ns].itypep >= 0) ++nwe;
	    sph_inf[i][ns].sph_we = nwe;
	    wesum += nwe;
	}
	for (ns = 0; ns <nsph[i]; ++ns)   sph_inf[i][ns].sph_we /= wesum;
    }

#ifdef DEBUG   
    for (i = 0; i < nout1; ++i) 
    {   
        printf("Decay number %d     nmscut= (%d,%d) nmsreg = (%d,%d)\n",
        i, nmscut[i][0], nmscut[i][1],  nmsreg[i][0], nmsreg[i][1]);
        
        {int  l,c;
          printf("kinematics= (");
          for (l=0;c= kinmtc_1[i].lvin[l];l++) printf("%d",c);
          printf(")->(");
          for (l=0;c=kinmtc_1[i].lvout[0][l];l++)  printf("%d",c);
          printf(")+(");
          for (l=0;c=kinmtc_1[i].lvout[1][l];l++)  printf("%d",c);
          printf(")\n");
        }        
        printf(" summas=(%f,%f)\n",summas[i][0],summas[i][1]);

        for (ns = 0; ns < nsph[i]; ++ns) 
	{   int c;
	    printf("   Sphere number = %d  weight=%f type=%d \n", 
	     ns, sph_inf[i][ns].sph_we, sph_inf[i][ns].itypep);
            printf("     pole vector(");
	    for(k=0; c=sph_inf[i][ns].lvpole[k]; k++) printf("%d",c);
            printf(")\n");
		
	    printf("    ncsreg=(%d,%d) ncscut=(%d,%d) tcscut=(%d,%d) \n",
	    sph_inf[i][ns].ncsreg[0], sph_inf[i][ns].ncsreg[1], 
            sph_inf[i][ns].ncscut[0], sph_inf[i][ns].ncscut[1],
            sph_inf[i][ns].tcscut[0], sph_inf[i][ns].tcscut[1]); 
	}
    }    
#endif 
    return ndim;
} /* mkmom_ */
