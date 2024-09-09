#include"micromegas.h"
#include"micromegas_aux.h"
#include"../CalcHEP_src/c_source/ntools/include/vegas.h"

/* Numerical recipes codes */

static double*dym,*dyt,*yt,*dysav,*ysav,*ytemp;
static int RKQCprnFlag=1;

static void rk4(double*y, double*dydx, int n, double x,double h,double * yout,
    void (*derivs)(double,double*,double*))
{
   int i;
   double hh=h/2, h6=h/6, xh=x+hh;

   for (i=0;i<n;i++) yt[i]=y[i]+hh*dydx[i];
   (*derivs)(xh,yt,dyt);
   for (i=0;i<n;i++) yt[i]=y[i]+hh*dyt[i];
   (*derivs)(xh,yt,dym);
   for (i=0;i<n;i++) { yt[i]=y[i]+h*dym[i]; dym[i] += dyt[i]; }
   (*derivs)(x+h,yt,dyt);
   for (i=0;i<n;i++) yout[i]=(y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]));
}


#define PGROW -0.20
#define PSHRNK -0.25
#define FCOR 0.06666666		/* 1.0/15.0 */
#define SAFETY 0.9
#define ERRCON 6.0e-4


static int rkqc(double * y, double * dydx, int n, double * x, double htry, 
     double eps,  double * yscal,  double* hdid, double* hnext,
     void (*derivs)(double,double *,double *))
{
   int i;
   double xsav=(*x),h=htry;

   for(i=0;i<n;i++) {ysav[i]=y[i]; dysav[i]=dydx[i];}

   for (;;) 
   {  double hh= 0.5*h, errmax=0;
      rk4(ysav,dysav,n,xsav,hh,ytemp,derivs);
      *x=xsav+hh; 
      (*derivs)(*x,ytemp,dydx);
      rk4(ytemp,dydx,n,*x,hh,y,derivs);
      *x=xsav+h;
      if (*x == xsav && RKQCprnFlag) 
      { printf("Step size too small in routine RKQC\n"); RKQCprnFlag=0;
        return 1;
      }
      rk4(ysav,dysav,n,xsav,h,ytemp,derivs);
      for (i=0;i<n;i++) 
      {  double temp;
         ytemp[i]=y[i]-ytemp[i];
         if(!isfinite( ytemp[i])) { errmax=ytemp[i]; break;} 
         temp= fabs(ytemp[i]/yscal[i]);
	 if (errmax < temp) errmax=temp;
      }
      if(!isfinite(errmax)){ h=h/10; continue;}
      errmax /= eps;

      if (errmax <= 1.0) 
      {
	 *hdid=h;
	 *hnext=((errmax > ERRCON ? SAFETY*h*exp(PGROW*log(errmax)) : 4*h));
	  break;
      }
      {  double h_=(SAFETY*h*exp(PSHRNK*log(errmax)));
//              if(h_/h>10 ) h=10*h;
//         else 
         if(h_/h<0.1) h=0.1*h; else  h=h_;
      }  
   }
   for (i=0;i<n;i++) y[i] += (double) (ytemp[i]*FCOR);
   return 0;
}

#define MAXSTP 10000000
#define TINY 1.0e-30


int  odeint(double * ystart, int nvar, double x1, double x2, double eps, 
         double h1, void (*derivs)(double,double *,double *))
{
   int nstp,i;
   double x,hnext,hdid,h;

   double *yscal,*y,*dydx;


   double ** allAlloc[9]={NULL,NULL,NULL,&dym,&dyt,&yt,&dysav,&ysav,&ytemp};
   allAlloc[0]=&yscal;allAlloc[1]=&y; allAlloc[2]=&dydx;
   for(i=0;i<9;i++) *allAlloc[i]=(double*)malloc(nvar*sizeof(double)); 

   RKQCprnFlag=1; 
   x=x1;
   h=((x2 > x1) ? fabs(h1) : -fabs(h1));
   for (i=0;i<nvar;i++) y[i]=ystart[i];
   for (nstp=1;nstp<=MAXSTP;nstp++) 
   {
      (*derivs)(x,y,dydx);
      for (i=0;i<nvar;i++) yscal[i]=(fabs(y[i])+fabs(dydx[i]*h)+TINY);

      if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
      if(rkqc(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs))break;
      
      if ((x-x2)*(x2-x1) >= 0.) 
      {
         for (i=0;i<nvar;i++) ystart[i]=y[i];
         for(i=0;i<9;i++) free(*allAlloc[i]);
         return 0;
      }
      h=hnext;
   }
   for(i=0;i<9;i++) free(*allAlloc[i]);
   return 1;
}


double bessI0(double  x)
{
	double ax,ans;
	double y;

	if ((ax=fabs(x)) < 3.75) {
		y=x/3.75;
		y*=y;
		ans= (1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
			+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2))))));
	} else {
		y=3.75/ax;
		ans= ((exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
			+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
			+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
			+y*0.392377e-2)))))))));
	}
	return ans;
}

static double bessI1(double x)
{
	double ax,ans;
	double y;

	if ((ax=fabs(x)) < 3.75) {
		y=x/3.75;
		y*=y;
		ans=(ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
			+y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3)))))));
	} else {
		y=3.75/ax;
		ans=(0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
			-y*0.420059e-2)));
		ans= (0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
			+y*(0.163801e-2+y*(-0.1031555e-1+y*ans)))));
		ans *= ((exp(ax)/sqrt(ax)));
	}
	return x < 0.0 ? -ans : ans;
}

double bessK0(double x)
{
/*
   M.Abramowitz and I.A.Stegun, Handbook of Mathematical Functions,
   Applied Mathematics Series vol. 55 (1964), Washington.
*/
         
	double y,ans;

	if (x <= 2.0) {
		y=x*x/4.0;
		ans=(-log(x/2.0)*bessI0(x))+(-0.57721566+y*(0.42278420
			+y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
			+y*(0.10750e-3+y*0.74e-5))))));
	} else {
		y=2.0/x;
		ans=(exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1
			+y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
			+y*(-0.251540e-2+y*0.53208e-3))))));
	}
	return  ans;
}


double bessK1(double x)
{
	double y,ans;

	if (x <= 2.0) {
		y=x*x/4.0;
		ans=(log(x/2.0)*bessI1(x))+(1.0/x)*(1.0+y*(0.15443144
			+y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
			+y*(-0.110404e-2+y*(-0.4686e-4)))))));
	} else {
		y=2.0/x;
		ans=(exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619
			+y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
			+y*(0.325614e-2+y*(-0.68245e-3)))))));
	}
	return  ans;
}

double bessK2(double x)
{
	double bk,bkm,bkp,tox;

	tox= 2.0/x;
	bkm=bessK0(x);
	bk=bessK1(x);
	bkp=bkm+tox*bk;
	bkm=bk;
	bk=bkp;
	return bk;
}

double K2pol(double x)
{
   if(x<0.1) return 1+ 1.875*x*(1+0.4375*x*(1-0.375*x));
   else      return bessK2(1/x)*exp(1/x)*sqrt(2/M_PI/x);
}

double K1pol(double x)
{
  if(x<0.1) return 1+ 0.375*x*(1-0.3125*x*(1+0.875*x));
  else      return bessK1(1/x)*exp(1/x)*sqrt(2/M_PI/x); 
}


static void del(int k, int * N, double *xa, double *ya)
{  int i;
    (*N)--;
    for(i=k;i<*N;i++) { xa[i]=xa[i+1]; ya[i]=ya[i+1];}    
} 

static void ins(int k,  double x, double y, int*N,double *xa,double *ya)
{ 
  int i;
  for(i=*N;i>k;i--) { xa[i]=xa[i-1]; ya[i]=ya[i-1];}
   xa[k]=x;ya[k]=y; (*N)++;        
}

int buildInterpolation(double (*Fun)(double), double x1,double x2, double eps,double delt, int*N_, double**xa_, double**ya_)
{  int i,cnt,N,k;
   double *xa,*ya,dx0;
   dx0=fabs(x2-x1)*delt;   
   N=5;
   xa=malloc(N*sizeof(double));
   ya=malloc(N*sizeof(double));
//printf("==================================\n");
   
   for(i=0;i<5;i++) {xa[i]=x1+ (x2-x1)/4*i; ya[i]=Fun(xa[i]); /*printf("x=%E y=%E\n",xa[i],ya[i]);*/ }  


   for(cnt=1;cnt;)
   { cnt=0; 
     for(i=0; i<N; i++)
     {  double x=xa[i], y=ya[i], yy;  
        if(i<N-1 && fabs(xa[i+1]-xa[i]) < dx0) continue; else
        if(i>0   && fabs(xa[i]-xa[i-1]) < dx0) continue;
         
        del(i,&N,xa,ya);
        yy=polint3(x, N, xa, ya);
        ins(i, x, y, &N,xa, ya);
        if( (eps>0 && fabs(yy-y) > eps) || (eps<0 && fabs(yy-y)> -eps*fabs(y)))  
        {
           cnt=1;
           xa=realloc(xa,sizeof(double)*(N+1));
           ya=realloc(ya,sizeof(double)*(N+1));
                if(i==0)   k=1;  
           else if(i==N-1) k=N-1;  
           else if(fabs(xa[i-1]-xa[i])< fabs(xa[i]-xa[i+1])) k=i+1;
           else  k=i;
           
           x=(xa[k-1]+xa[k])/2;
           y=Fun(x); 
//           printf("x=%E y=%E\n",x,y);
           
           ins(k,x,y,&N,xa,ya);
           i++;     
        }  
     }
   }   
   *N_=N;
   *xa_=xa;
   *ya_=ya;
   return 0;  
}



double MaxGapLim(double x, double mu) 
/* S.Yellin, Phys.Rev. D66,032005(2002)

   There is a theoretical model which predicts homogenious event distribution 
   with everage number of events mu. Let experiment gets a gap bitween points 
   where according to theory x point are expected. Then the theoretical model 
   is non-confirmed with probability MaxGap   
*/  
{
  int k;
  long double C0,lnkf;
  if(x>=mu) return 1;
  double lnMu[10]={0,           6.931472E-01,1.609438E+00,2.302585E+00,2.995732E+00,3.912023E+00,4.605170E+00,5.298317E+00,5.991465E+00,6.907755E+00};
  double xMin[10]={1.864675e-01,3.049546e-01,5.523649e-01,8.225020e-01,1.174237e+00,1.759748e+00,2.286058e+00,2.863048e+00,3.480203e+00,4.342722e+00};
  
  
  double xx=polint3(log(mu),10,lnMu,xMin);
  if(x<xx) return 1E-6;
  
  for(k=0,C0=0,lnkf=0;k<mu/x; k++,lnkf+=logl(k)) {long double dC=expl(-k*x +k*logl(mu-k*x) -lnkf)*(1+k/(mu-k*(long double)x));
                                                  if(k&1) C0-=dC; else C0+=dC;
//                                                  printf("k=%d dC=%E lnkf=%e  C0=%E \n", k, (double)dC,(double)lnkf,(double)C0);   
  }
  return C0;   
}


#define BUFFSIZE 500

int readTable(char * fileName, int *Ncolumn, double **tab)
{  FILE *f;
   char buff[BUFFSIZE];   
   int nRec=0,nCol=0,nCom=0;

   f=fopen(fileName,"r");
   
   
   if(!f) return 0;

   while(fgets(buff,BUFFSIZE,f))
   { int i;
     char*ch;
     for(i=0; buff[i] && buff[i]==' ';i++);
     if(buff[i]==0 || buff[i]=='#') {nCom++; continue;}
     ch=strtok(buff," \n");
     if(ch[0]=='#' || ch[0]==0) continue;
     for(i=0;ch;i++,ch=strtok(NULL," \n"))
     { 
       if(nRec==0) {tab[i]=malloc(sizeof(double)); nCol++;} 
       else
       { if(i==nCol){fclose(f); for(i=0;i<nCol;i++) free(tab[i]);  return -(nRec+1+nCom);}
         tab[i]=realloc(tab[i],(nRec+1)*sizeof(double));
       }
       if(1!=sscanf(ch,"%lf",tab[i]+nRec)) break;
     } 
     nRec++;
   }
   fclose(f);
   if(Ncolumn) *Ncolumn=nCol;
   return nRec;
}

static double amotry(double *p, double *y, int ndim,
	     double (*f)(double *), int ilo, double fac)
{
   int i,j;
   double  ytry, fac1=(1.0-fac)/ndim;
   double * p_buff=p+(ndim+1)*ndim;
   double * p_ilo =p+ilo*ndim;
   
   for(j=0;j<ndim;j++) p_buff[j]=p_ilo[j]*fac;
   for(i=0;i<=ndim;i++)  if(i!=ilo) 
     {double *p_i=p+i*ndim;  for(j=0;j<ndim;j++) p_buff[j] +=p_i[j]*fac1;} 
   ytry=-(*f)(p_buff);
   
   if (ytry > y[ilo]) 
   {  for(j=0;j<ndim;j++) p_ilo[j]=p_buff[j];
      y[ilo]=ytry;
   }
//printf("amotry returns fac=%f %E\n",fac,  ytry);   
   return ytry;
}

double amoeba(double *pp, double * dp, int ndim, double (*f)(double *), 
                                                    double eps, int *nCalls)
{
   int i,ilo,ihi,inlo,j;
   double ysave,ytry;
   double *p=(double*)malloc(ndim*(ndim+2)*sizeof(double));
   double *y=(double*)malloc((ndim+1)*sizeof(double));
   

   for(j=0;j<=ndim;j++) for(i=0;i<ndim;i++)  p[i+j*ndim]=pp[i];
   for(j=1;j<=ndim;j++) p[j-1+j*ndim]+=dp[j-1];
 
   for(j=0;j<=ndim;j++) y[j]=-f(p+j*ndim);    

//   for(j=0;j<=ndim;j++) printf(" %e %e %e   %e\n",p[0+j*ndim],p[1+j*ndim],p[2+j*ndim], -y[j]);
   
   
   for (;;) 
   {
      ihi=0;									     
      ilo = y[0]<y[1] ? (inlo=1,0) : (inlo=0,1);				     
      for (i=0;i<=ndim;i++)							     
      {										     
     	 if (y[i] >= y[ihi]) ihi=i;
     	 if (y[i] < y[ilo]) { inlo=ilo; ilo=i; } 
     	 else if (y[i] < y[inlo] && i != ilo) inlo=i;
      }										     

//printf("nCall=%d  ndim=%E\n",*nCalls,y[ilo]);   
     										     
      if( (nCalls &&  (*nCalls)<=0)||2*(y[ihi]-y[ilo])/(fabs(y[ilo])+fabs(y[ihi]))<eps)break;
     										     
      ytry=amotry(p,y,ndim,f,ilo,-1.0);  if(nCalls) (*nCalls)--;				     
      if (ytry >= y[ihi]) {ytry=amotry(p,y,ndim,f,ilo,2.);  if(nCalls) (*nCalls)--;}	     
      else if (ytry <= y[inlo])							     
      {										     
         ysave=y[ilo];								     
     	 ytry=amotry(p,y,ndim,f,ilo,0.5);  if(nCalls)(*nCalls)--;
     	 if (ytry <= ysave)
     	 {  
     	    for (i=0;i<=ndim;i++)
     	    {  double * p_ihi=p+ihi*ndim;
               if (i != ihi)
     	       {  double * p_i=p+i*ndim;
     		  for(j=0;j<ndim;j++) p_i[j]=0.5*(p_i[j]+p_ihi[j]);
     		  y[i]=-(*f)(p_i);
     	       }
            }
/*printf("srink\n");            */
     	    if(nCalls)(*nCalls) -= ndim;
         }									     
      }										     
   }
   for(i=0;i<ndim;i++){ pp[i]=p[i+ihi*ndim]; dp[i]=p[i+ilo*ndim]-p[i+ihi*ndim];}
   ysave=y[ihi];
   free(y); free(p);
   return -ysave;
}
/*========================== end of amoeba ================*/

#define MAXSTEP 15
double vegas_chain(int ndim, double (*Integrand)(double*, double),
int N0, double Nfact, double eps,double * dI,void (*clean)(void))   
{ vegasGrid *vegPtr=NULL;
  int k,l;
  double ti[MAXSTEP],dti[MAXSTEP];
  double ii,dii,chi2;

  vegPtr=vegas_init(ndim,Integrand,50);
  
  for(k=0;k<MAXSTEP;k++)
  { double s0=0,s1=0,s2=0; 
    vegas_int(vegPtr, N0 , 1.5, nPROCSS, ti+k, dti+k);
    if(clean) clean();
    printf("ti=%E dti=%E\n",ti[k], dti[k]);
    if(dti[k]==0) { ii=ti[k];  break;}
    for(l=k;l>=k/2;l--)
    { s0+=1/(dti[l]*dti[l]);
      s1+=ti[l]/(dti[l]*dti[l]);
      s2+=ti[l]*ti[l]/(dti[l]*dti[l]);
      if(l!=k)
      { 
        ii=s1/s0;
        dii=1/sqrt(s0);
        chi2=(s2-s1*s1/s0)/(k-l+1);
        if(chi2> 1 )dii*=sqrt(chi2);
        if(dii<eps*fabs(ii)) break;
      }  
    }
    if(k && dii<eps*fabs(ii)) break;
    N0*=Nfact;    
  }
  vegas_finish(vegPtr);
  if(dI) *dI=dii;
  return ii;
}  



void spline(double x[], double y[], int n, double y2[])
{
	int i,k;
	double p,qn,sig,un,*u;

	u=malloc(n*sizeof(double));
	
	y2[0]=u[0]=0.0;
	for (i=1;i<n-1;i++) 
	{
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2;
		y2[i]=(sig-1)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	
	qn=un=0;
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1);
	for (k=n-2;k>=0;k--) y2[k]=y2[k]*y2[k+1]+u[k];
	free(u);
}

void splint(double xa[], double ya[], double y2a[], int n, double x, double *y)
{
	int klo,khi,k;
	double h,b,a;

        if(xa[0]<xa[n-1]) {klo=0;  khi=n-1;} 
        else              {klo=n-1;khi=0;  }
	while (abs(khi-klo) > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}   
	h=xa[khi]-xa[klo];
	if (h == 0.0) printf("Bad xa input to routine splint");
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

#define FC_DIM 15000

static int findMinLimit(double mu,double b,double cl)
{
  double r[FC_DIM],p[FC_DIM],s,rMax;
  int N,n,n0,nMin,nMax;
  
  r[0]=exp(-mu);
  p[0]=exp(-mu-b);
  s=p[0];
  nMax=0;
  rMax=r[0];
   for(n=1;n<FC_DIM ;n++) p[n]=1;
  for(n=1; 1-s> cl/100 ;n++)
  {  double mu_best;
     if(n>b) mu_best=n-b; else mu_best=0;
     r[n]=pow((mu+b)/(mu_best+b),n)*exp(mu_best-mu);
     if(r[n]>rMax) {rMax=r[n];nMax=n;}
     p[n]=p[n-1]*(mu+b)/n;
     s+=p[n];
  }
  N=n;
  cl-=p[nMax]; nMin=nMax-1; nMax++; 
  for( ;cl>0;)
  {
    if(nMin>=0) { if(nMax<N && r[nMax]>r[nMin]) cl-=p[nMax++];  else  cl-=p[nMin--];}
    else  cl-=p[nMax++];
  }
//  if(nMin>0) return nMin; else return 0; 
   return nMin; 
} 


double FeldmanCousins(int n0, double b, double cl)
{
   double delta,mu;
   
   for(mu=0,delta=1; delta>0.001;  mu-=delta,delta/=2 )
   {  
     for(;findMinLimit(mu, b,cl)<n0;mu+=delta);
   }    
   return mu;
}

static double alpha=0.25,cc, p_exp;
static double ch2pva_int(double f)
{ 
  if(f<=0) return 0;
  return pow(f,(1-alpha)/alpha)* pow(-log(f)/alpha,p_exp)/cc;

}

double ch2pval(int nexp, double ch2obs)
{
   p_exp=0.5*nexp-1;  
   alpha=1/(1+p_exp/log(2));
  cc=exp(lgamma(nexp/2.))*alpha;
//   displayFunc(ch2pva_int,0,exp(-alpha*ch2obs/2),"ch2pva_int");   
   return simpson(ch2pva_int, 0, exp(-alpha*ch2obs/2), 1.E-3,NULL);///exp(lgamma(nexp/2.))/alpha;
}


double plr2pval(double l) { return  0.5*(1- erfl(sqrtl(-logl(l))));}


double pval2plr(double p) 
{  if(p>0.5) return 1;

   double l1=0,l2=1,p1=0,p2=0.5;
   
   for(;;)
   {  
      double l_=0.5*(l1+l2);   
      double p_= plr2pval(l_);
      if(fabs(p-p_) < 1E-3*p) return l_;
      
      if(p_<p) { l1=l_; p1=p_;} else {  l2=l_; p2=p_;}
   }
   
}

void  addErrorMess( char** All, char * one)
{
   if(*All==NULL || strstr(*All,one)==NULL)
   { int len; 
     if(*All) len=strlen(*All); else len=0;
     *All=realloc(*All, len+2+strlen(one)); 
     sprintf(*All+len,one);
   }
}

void delInterval(double x1,double x2,double **intervals, int*n)
{ 

//printf("delInterval [%e,%e]\n",x1,x2);
 
  int i1,i2;
  for(i1=0;i1<2*(*n);i1++) if(x1<=(*intervals)[i1]) break; i1--; 
  for(i2=0;i2<2*(*n);i2++) if(x2<(*intervals)[i2]) break; i2--;
  if(i1==i2)
  { 
//printf("i1=i2=%d  x1=%e x2=%e\n",i1,x1,x2);  
      if(i1&1==1) return; // this interval is already excluded
     *intervals=realloc(*intervals,2*(*n+1)*sizeof(double));
     for(int k=2*(*n)+1;k>=i1+3;k--) (*intervals)[k]=(*intervals)[k-2]; 
//printf(" %e %e %e %e\n", (*intervals)[0], (*intervals)[1], (*intervals)[2], (*intervals)[3]);      
     (*intervals)[i1+1]=x1; (*intervals)[i1+2]=x2;
//printf(" %e %e %e %e\n", (*intervals)[0], (*intervals)[1], (*intervals)[2], (*intervals)[3]);     
     (*n)++;  
     return;
  }
  if((i1&1)==0) (*intervals)[++i1]=x1;  
  if(!(i2&1)) (*intervals)[i2]=x2; else i2++;
  if(i2-i1>1)
  {
     for(int k=i1+1;k<2*(*n)-2;k++) (*intervals)[k]=(*intervals)[k-i1+i2-1];
     (*n)-=(i2-i1)/2;
  }  
}
