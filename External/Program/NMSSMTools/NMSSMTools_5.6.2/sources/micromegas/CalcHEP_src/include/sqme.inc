
static REAL computer_eps=0.;

static void get_eps(void)
{  
   REAL z;
   for(computer_eps=1.E-8;;computer_eps/=1.5){z=1+computer_eps; if(z==1) return;}
}

int calcFunc_ext(void) 
{ int nsub,err,res;
  if(computer_eps==0.) get_eps();  
  res=calcFunc_stat();
  for(nsub=1;nsub<=nprc_ext;nsub++) sqme_ext(nsub,0,NULL,NULL,&err);
  return res;  
}

#define buffSize 100

double sqme_ext(int nsub, double GG,REAL*momenta,REAL*cb_coeff, int * err)
{
  int i,nperm,shift, ntot=nin_ext+nout_ext,cb_pow=cb_ext[nsub-1].pow;
  REAL result, momenta_[4*maxNp], cb_buf_[buffSize],*cb_buff=NULL;
  double Fmax=0; 
  if(!momenta){ darr[nsub-1](GG,NULL,NULL,NULL,err); return 0;}


  if(cb_coeff) for(i=0;i<cb_pow;i++) cb_coeff[i]=0;

  if(nin_ext==2)
  {  REAL p1p2=momenta[0]*momenta[4]-momenta[1]*momenta[5]-momenta[2]*momenta[6]-momenta[3]*momenta[7];
     REAL m1,m2; 
     pinf_ext(nsub,1,&m1,NULL); 
     pinf_ext(nsub,2,&m2,NULL);
     N_pol_11_= p1p2/sqrt(p1p2*p1p2-m1*m1*m2*m2)/m1;
//     printf(" p1p2=%E m1*m2=%E m1=%E \n", (double)p1p2, (double)(m1*m2), (double)(m1));
     N_pol_12_=-m1/sqrt(p1p2*p1p2-m1*m1*m2*m2);
     N_pol_22_= p1p2/sqrt(p1p2*p1p2-m1*m1*m2*m2)/m2;
     N_pol_21_=-m2/sqrt(p1p2*p1p2-m1*m1*m2*m2);
  }
   
  result=darr[nsub-1](GG,momenta,cb_coeff,&Fmax, err);
  nperm=permMap[nsub-1][1];
  if(nperm)
  { int n; 
     shift=permMap[nsub-1][0];
     if(cb_coeff && cb_pow)
     { if(cb_pow<=buffSize) cb_buff=cb_buf_; else cb_buff=malloc(sizeof(REAL)*cb_pow); }

     for(n=0;n<nperm;n++)
     { 
       if(cb_buff) for(i=0;i<cb_pow;i++) cb_buff[i]=0;
       for(i=0;i<ntot;i++) memcpy(momenta_+i*4,momenta+ 4*(permP[shift+n][i]-1),sizeof(REAL)*4);
       result+=darr[nsub-1](GG,momenta_,cb_coeff,&Fmax,err);
       if(cb_buff) for(i=0;i<cb_pow;i++) cb_coeff[permC[shift+n][i]-1]+=cb_buff[i];
     }
     result/=(1+nperm);
     if(cb_coeff) for(i=0;i<cb_pow;i++) cb_coeff[i]/=1+nperm;
     if(cb_buff && cb_buff !=cb_buf_) free(cb_buff);
  }

  if(!(*err) && Fmax*computer_eps*1000 > (result>0 ? result : -result) )
  { *err=1;
     if(result>0) { result-=(1000*Fmax*computer_eps -result)*0.1  ; if(result<0) result=0;}
     else { result+=(1000*Fmax*computer_eps+result)*0.05  ; if(result>0) result=0;}
  }


   return (double)result;
}
