
/* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */

#include <math.h>
#include <stdlib.h>
#include<stdio.h>
#include"micromegas_aux.h"

//=====  ludcmp & lubksb auxiliary code
#define TINY 1.0e-30

static double*vv=NULL;

static int ludcmp(double*a, int n, int* indx, double* d)
{
   int i,imax,j,k;
   double big,dum,sum,temp;

   *d=1.0;
   for(i=0;i<n;i++) 
   {
      big=0.0;
      for(j=0;j<n;j++)if ((temp=fabs(a[i*n+j])) > big) big=temp;
      if (big == 0.0) return 1;
      vv[i]=1.0/big;
   }
   for(j=0;j<n;j++) 
   {
     for (i=0;i<j;i++) 
     { 
        sum=a[i*n+j];
        for (k=0;k<i;k++) sum -= a[i*n+k]*a[k*n+j];
        a[i*n+j]=sum;
     }
     big=0;
     for(i=j;i<n;i++) 
     {
        sum=a[i*n+j];
        for(k=0;k<j;k++)sum -= a[i*n+k]*a[k*n+j];
	a[i*n+j]=sum;
        if( (dum=vv[i]*fabs(sum)) >= big) {big=dum; imax=i;}
     }
     if(j != imax) 
     {   for (k=0;k<n;k++) 
         {   
	     dum=a[imax*n+k];
	     a[imax*n+k]=a[j*n+k];
	     a[j*n+k]=dum;
	}
	*d = -(*d);
	vv[imax]=vv[j];
     }
     indx[j]=imax;
     if (a[j*n+j] == 0.0) a[j*n+j]=(double) TINY;
     if (j != n-1) 
     {
        dum= (1.0/(a[j*n+j]));
        for (i=j+1;i<n;i++) a[i*n+j] *= dum;
     }
   }
}


static void lubksb(double* a, int n, int* indx, double b[])
{
   int i,ii=-1,ip,j;
   double sum;

   for (i=0;i<n;i++) 
   {
      ip=indx[i];     
      sum=b[ip];
      b[ip]=b[i];
      if(ii>=0) for (j=ii;j<=i-1;j++) sum -= a[i*n+j]*b[j]; else if (sum) ii=i;
      b[i]=sum;
   }
   for (i=n-1;i>=0;i--) 
   {
      sum=b[i];
      for (j=i+1;j<n;j++) sum -= a[i*n+j]*b[j];
      b[i]=sum/a[i*n+i];
   }
}

//======  common variables for stiff & stifbs

static double *a=NULL, *ysav=NULL,*yerr=NULL,*dydx=NULL,*dfdx=NULL,*dfdy=NULL;
static int *indx;

// ====== stiff() 

#define NRANSI
#define SAFETY 0.9
#define GROW 1.5
#define PGROW -0.25
#define SHRNK 0.5
#define PSHRNK (-1.0/3.0)
#define ERRCON 0.1296
#define MAXTRY 40
#define GAM (1.0/2.0)
#define A21 2.0
#define A31 (48.0/25.0)
#define A32 (6.0/25.0)
#define C21 -8.0
#define C31 (372.0/25.0)
#define C32 (12.0/5.0)
#define C41 (-112.0/125.0)
#define C42 (-54.0/125.0)
#define C43 (-2.0/5.0)
#define B1 (19.0/9.0)
#define B2 (1.0/2.0)
#define B3 (25.0/108.0)
#define B4 (125.0/108.0)
#define E1 (17.0/54.0)
#define E2 (7.0/36.0)
#define E3 0.0
#define E4 (125.0/108.0)
#define C1X (1.0/2.0)
#define C2X (-3.0/2.0)
#define C3X (121.0/50.0)
#define C4X (29.0/250.0)
#define A2X 1.0
#define A3X (3.0/5.0)


static double *dysav=NULL,*g1=NULL,*g2=NULL,*g3=NULL,*g4=NULL,*del=NULL,*ytemp=NULL;


static int  stiff_1(double y[], int n, double *x, double htry, double eps,
	double*yscal, double *hdid, double *hnext,
	void (*derivs)(double,double*,double*,double,double*,double*))
{
    int i,j,jtry;
    double d,errmax,h,xsav;

    xsav=(*x);    
    derivs(xsav,y,dydx,htry,dfdx,dfdy);

    for (i=0;i<n;i++) 
    {
       ysav[i]=y[i];
       dysav[i]=dydx[i];
    }
    h=htry;
    for (jtry=1;jtry<=MAXTRY;jtry++) 
    {    
    
      for(i=0;i<n;i++) 
      {
         for (j=0;j<n;j++) a[i*n+j] = -dfdy[i*n+j];
         a[i*n+i] += 1.0/(GAM*h);
      }
      ludcmp(a,n,indx,&d);
      for (i=0;i<n;i++) g1[i]=dysav[i]+h*C1X*dfdx[i];
      lubksb(a,n,indx,g1);
      for (i=0;i<n;i++) y[i]=ysav[i]+A21*g1[i];
      *x=xsav+A2X*h;
      derivs(*x,y,dydx,h,NULL,NULL);
      for(i=0;i<n;i++) g2[i]=dydx[i]+h*C2X*dfdx[i]+C21*g1[i]/h;
      lubksb(a,n,indx,g2);
      for (i=0;i<n;i++) y[i]=ysav[i]+A31*g1[i]+A32*g2[i];
      *x=xsav+A3X*h;
      derivs(*x,y,dydx,h,NULL,NULL);
      for (i=0;i<n;i++) g3[i]=dydx[i]+h*C3X*dfdx[i]+(C31*g1[i]+C32*g2[i])/h;
      lubksb(a,n,indx,g3);
      for(i=0;i<n;i++) g4[i]=dydx[i]+h*C4X*dfdx[i]+(C41*g1[i]+C42*g2[i]+C43*g3[i])/h;
      lubksb(a,n,indx,g4);
      for(i=0;i<n;i++) 
      {
        y[i]=ysav[i]+B1*g1[i]+B2*g2[i]+B3*g3[i]+B4*g4[i];
        yerr[i]=E1*g1[i]+E2*g2[i]+E3*g3[i]+E4*g4[i];
      }
      *x=xsav+h;
      if (*x == xsav) { printf("stepsize not significant in stiff\n"); return 1;}
      errmax=0.0;
      if(yscal) for(i=0;i<n;i++) if(errmax<fabs(yerr[i]/yscal[i])) errmax=fabs(yerr[i]/yscal[i]);
      else      for(i=0;i<n;i++) if(errmax<fabs(yerr[i])) errmax=fabs(yerr[i]);
      errmax /= eps;
       
      if (errmax <= 1) 
      {
        *hdid=h;
        *hnext=(errmax > ERRCON ? SAFETY*h*pow(errmax,PGROW) : GROW*h);
        return 0;
      } else 
      {
         *hnext=SAFETY*h*pow(errmax,PSHRNK);
         if (*hnext < SHRNK*h) *hnext=SHRNK*h;
         h=(*hnext); 
      }
    }
    printf("exceeded MAXTRY in stiff\n");
    return 2;
}

#undef SAFETY
#undef GROW
#undef PGROW
#undef SHRNK
#undef PSHRNK
#undef ERRCON
#undef MAXTRY
#undef GAM
#undef A21
#undef A31
#undef A32
#undef C21
#undef C31
#undef C32
#undef C41
#undef C42
#undef C43
#undef B1
#undef B2
#undef B3
#undef B4
#undef E1
#undef E2
#undef E3
#undef E4
#undef C1X
#undef C2X
#undef C3X
#undef C4X
#undef A2X
#undef A3X
#undef NRANSI

/* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */


int stiff(int first,double xstart, double xend, int n, double*y, double *yscal, double eps, double*htry, 
    void (*derivs)(double, double*, double*, double,double*,double*))
{ 
  double x=xstart;  
  double _h,h,h_;
  int cerr;
  
  
  h=*htry;
  if(first)
  { int mem=n*sizeof(double);

    vv   =realloc(vv   ,mem);
    indx =realloc(indx ,n*sizeof(int));

    dfdx =realloc(dfdx ,mem);
    dfdy =realloc(dfdy ,mem*n);
    yerr =realloc(yerr ,mem);
    dydx =realloc(dydx ,mem); 
    ysav =realloc(ysav ,mem);

    dysav=realloc(dysav,mem);
    a    =realloc(a    ,mem*n);
    g1   =realloc(g1   ,mem);
    g2   =realloc(g2   ,mem);
    g3   =realloc(g3   ,mem);
    g4   =realloc(g4   ,mem);
  } 
       
  for(; fabs(x-xend)>1.E-5*(fabs(xstart)+fabs(xend)); )
  { 
    h=fabs(h);
    if(h>fabs(x-xend)) h=fabs(x-xend);
    if(x>xend) h*=-1;
    cerr=stiff_1(y,n,&x,h,eps,yscal,&_h,&h_,derivs);
    if(cerr) return cerr;
    h=h_;
//     printf("T=%e Y= %e %e hdid=%e hnext %e\n",x, y[0],y[1],_h,h_);
  }
  *htry=h;
  return 0; 
}

//===== stifbs()

#define NRANSI
#define KMAXX 7
#define IMAXX (KMAXX+1)
#define SAFE1 0.25
#define SAFE2 0.7
#define REDMAX 1.0e-5
#define REDMIN 0.7
#define SCALMX 0.1

#define FMAX(x,y)  (x>y? x:y)
#define FMIN(x,y)  (x<y? x:y)


static double *yseq=NULL,*d=NULL,*cc=NULL;


static void simpr(double y[], double dydx[], double dfdx[], double *dfdy, int n,
	double xs, double htot, int nstep, double yout[],
	void (*derivs)(double, double*, double*, double,double*, double*))
{
	int i,j,nn;
	double d,h,x;

	h=htot/nstep;
	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) a[i*n+j] = -h*dfdy[i*n+j];
		++a[i*n+i];
	}
	ludcmp(a,n,indx,&d);
	for (i=0;i<n;i++) yout[i]=h*(dydx[i]+h*dfdx[i]);
	lubksb(a,n,indx,yout);
	for (i=0;i<n;i++) ytemp[i]=y[i]+(del[i]=yout[i]);
	x=xs+h;
	(*derivs)(x,ytemp,yout,h,NULL,NULL);
	for (nn=2;nn<=nstep;nn++) {
		for (i=0;i<n;i++) yout[i]=h*yout[i]-del[i];
		lubksb(a,n,indx,yout);
		for (i=0;i<n;i++) ytemp[i] += (del[i] += 2.0*yout[i]);
		x += h;
		(*derivs)(x,ytemp,yout,h,NULL,NULL);
	}
	for (i=0;i<n;i++) yout[i]=h*yout[i]-del[i];
	lubksb(a,n,indx,yout);
	for (i=0;i<n;i++) yout[i] += ytemp[i];
}

static double x[KMAXX];



static void pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv)
{
	int k1,j;
	double q,f2,f1,delta;

	x[iest]=xest;
	for (j=0;j<nv;j++) dy[j]=yz[j]=yest[j];
	if (iest == 0) 
	{
		for (j=0;j<nv;j++) d[j*KMAXX]=yest[j];
	} else {
		for (j=0;j<nv;j++) cc[j]=yest[j];
		for (k1=1;k1<=iest;k1++) {
			delta=1.0/(x[iest-k1]-xest);
			f1=xest*delta;
			f2=x[iest-k1]*delta;
			for (j=0;j<nv;j++) {
				q=d[j*KMAXX+k1-1];
				d[j*KMAXX+k1-1]=dy[j];
				delta=cc[j]-q;
				dy[j]=f1*delta;
				cc[j]=f2*delta;
				yz[j] += dy[j];
			}
		}
		for (j=0;j<nv;j++) d[j*KMAXX+iest]=dy[j];
	}
}



static int stifbs_1(int first, double y[], int nv, double *xx, double htry, double eps,
	double*yscal, double *hdid, double *hnext,
	void (*derivs)(double, double*,double*, double,double*, double*))
{
	int i,iq,k,kk,km;
//     static int first=1;	
        static int kmax,kopt,nvold = -1;
        static double epsold = -1.0,xnew;
        static double a[IMAXX+1],alf[KMAXX+1][KMAXX+1];
	double eps1,errmax,fact,h,red,scale,work,wrkmin,xest;
	double err[KMAXX];
	static int nseq[IMAXX]={2,6,10,14,22,34,50,70};
	int reduct,exitflag=0;
	

//	if(eps != epsold || nv != nvold) 
        if(first)
	{
		*hnext = xnew = -1.0e29;
		eps1=SAFE1*eps;
		a[0]=nseq[0]+1;
		for (k=0;k<KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
		for (iq=1;iq<KMAXX;iq++) 
		{
			for (k=0;k<iq;k++)
				alf[k][iq]=pow(eps1,((a[k+1]-a[iq+1])/
					((a[iq+1]-a[0]+1.0)*(2*k+3))));
		}
		epsold=eps;
		nvold=nv;
		a[0] += nv;
		for (k=0;k<KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
		for (kopt=1;kopt<KMAXX-1;kopt++)
			if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;
		kmax=kopt;
	}
	h=htry;
	for (i=0;i<nv;i++) ysav[i]=y[i];
	derivs(*xx,y,dydx, h,dfdx,dfdy);
	if (*xx != xnew || h != (*hnext)) 
	{
           first=1;
           kopt=kmax;
	}
	reduct=0;
	for (;;) 
	{
		for (k=0;k<kmax;k++) 
		{
			xnew=(*xx)+h;
			if (xnew == (*xx))
			{  printf("step size underflow in stifbs_1\n");
			   return 1;
                        }
			simpr(ysav,dydx,dfdx,dfdy,nv,*xx,h,nseq[k],yseq,derivs);
			xest=(h/nseq[k])*(h/nseq[k]);
			
			pzextr(k,xest,yseq,y,yerr,nv);

			
			if (k != 0) 
			{
				errmax=TINY;
				if(yscal) for(i=0;i<nv;i++) errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));
				else      for(i=0;i<nv;i++) errmax=FMAX(errmax,fabs(yerr[i]));
				errmax /= eps;
				km=k-1;
				err[km]=pow(errmax/SAFE1,1.0/(2*km+3));
			}
			if (k != 0 && (k >= kopt-1 || first)) 
			{
				if (errmax < 1.0) 
				{
					exitflag=1;
					break;
				}
				if (k == kmax || k == kopt+1) 
				{
					red=SAFE2/err[km];
					break;
				}
				else if (k == kopt && alf[kopt-1][kopt] < err[km]) 
				{
						red=1.0/err[km];
						break;
                           }
				else if (kopt == kmax && alf[km][kmax-1] < err[km]) 
				{
						red=alf[km][kmax-1]*SAFE2/err[km];
						break;
                           }
				else if (alf[km][kopt] < err[km]) 
				{
					red=alf[km][kopt-1]/err[km];
					break;
				}
			}
		}
				
		if (exitflag) break;
		red=FMIN(red,REDMIN);
		red=FMAX(red,REDMAX);
		h *= red;
		reduct=1;
	}	
	*xx=xnew;
	*hdid=h;
	first=0;
	wrkmin=1.0e35;
	for (kk=0;kk<=km;kk++) {
		fact=FMAX(err[kk],SCALMX);
		work=fact*a[kk+1];
		if (work < wrkmin) {
			scale=fact;
			wrkmin=work;
			kopt=kk+1;
		}
	}
	*hnext=h/scale;
	if (kopt >= k && kopt != kmax && !reduct) {
		fact=FMAX(scale/alf[kopt-1][kopt],SCALMX);
		if (a[kopt+1]*fact <= wrkmin) {
			*hnext=h/fact;
			kopt++;
		}
	}
    return 0;	
}

int stifbs(int first,double xstart, double xend, int nv, double*y, double *yscal, double eps, double*htry,
   void (*derivs)(double, double*,double*, double,double*,double*))
{   
   double x=xstart;
   double _h,h,h_;
   int cerr;
   
   if(first)
   { int mem=nv*sizeof(double);
	vv   =realloc(vv  ,mem);
        indx =realloc(indx,nv*sizeof(int));
	
	dfdx =realloc(dfdx,mem);
	dfdy =realloc(dfdy,mem*nv);
	ysav =realloc(ysav,mem);
        dydx =realloc(dydx,mem);
	yerr =realloc(yerr,mem);
   
	a    =realloc(a,   mem*nv);
	d    =realloc(d   ,KMAXX*mem);
	yseq =realloc(yseq,mem);
        cc   =realloc(cc  ,mem);
        ytemp=realloc(ytemp,mem);
        del  =realloc(del  ,mem);
   }   
   h=*htry;
   
   for(; fabs(x-xend)>1.E-5*(fabs(xstart)+fabs(xend));first=0 )
   { 
      h=fabs(h);
      if(h>fabs(x-xend)) h=fabs(x-xend);
      if(x>xend) h*=-1;
      cerr=stifbs_1(first,y,nv,&x,h,eps,yscal,&_h,&h_,derivs);
      if(cerr) return cerr;
      h=h_;
//      printf("x=%e Y= %e %e %e hdid=%e hnext %e\n",x, y[0],y[1],y[2],_h,h_);
   }
   *htry=h;
   return 0;
}


#undef KMAXX
#undef IMAXX
#undef SAFE1
#undef SAFE2
#undef REDMAX
#undef REDMIN
#undef SCALMX
#undef NRANSI

//#define MAIN
#ifdef MAIN

void derivs(double x,  double*y, double*f,double h,double*dfdx,double*dfdy)
{ 
  int i,n=3;
  f[0]= -0.013*y[0]-1000*y[0]*y[2];
  f[1]= -2500*y[1]*y[2];
  f[2]= -0.013*y[0]-1000*y[0]*y[2]-2500*y[1]*y[2];
  
  if(dfdx) for(i=0;i<3;i++) dfdx[i]=0;
  if(dfdy)
  {
    dfdy[0*n+0]= -0.013-1000*y[2];
    dfdy[0*n+1]=  0;
    dfdy[0*n+2]= -1000*y[0];
    dfdy[1*n+0]=  0;
    dfdy[1*n+1]= -2500*y[2];
    dfdy[1*n+2]= -2500*y[1];
    dfdy[2*n+0]= -0.013-1000*y[2];
    dfdy[2*n+1]= -2500*y[2];
    dfdy[2*n+2]= -1000*y[0]-2500*y[1];  
  }  
}


int main(void)
{
  int n=3;
  double x=0;
  double eps=1E-4;
  double h=2.9E-4;
  double y[3]={1,1,0};
  double yscal[3]={1,1,1};
  double dy[3];

 
  stifbs(1,0, 50, n, y,yscal,eps,&h,derivs);  
  
printf("y= %E %E %E\n", y[0],y[1],y[2]);
//  testing point stifbs: Y= 5.976547e-01 1.402343e+00 -1.893377e-06 
//                stiff : Y= 5.976613e-01 1.402337e+00 -1.890773e-06
                 
}


#endif
