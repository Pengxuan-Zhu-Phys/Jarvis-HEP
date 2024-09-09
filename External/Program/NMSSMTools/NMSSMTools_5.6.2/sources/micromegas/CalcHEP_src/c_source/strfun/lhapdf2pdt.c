#include<stdlib.h>
#include<string.h>
#include<stdio.h> 

#include "polint.h"

int main(int argc, char** argv)
{ 
   if(argc<2) { printf("One argument is expected,  name  of directory with  LHAPDF distribution.\n");
                printf("Second  argument if presented, corresponds to set number.\n");
                exit(1);
              }
              
   char *dirname=malloc(strlen(argv[1])+1);
   strcpy(dirname,argv[1]);
   char *ch=dirname+strlen(dirname)-1;
   if(ch[0]=='/') { ch[0]=0; ch--;}
   while(ch>dirname && ch[0]!='/') ch--;
   if(ch[0]=='/') ch++;
   char* finfo=malloc(strlen(dirname)+1+strlen(ch)+6);
   sprintf(finfo,"%s/%s.info",dirname,ch);
   char* fdat=malloc(strlen(dirname)+1+strlen(ch)+10);
   sprintf(fdat,"%s/%s_",dirname,ch);
   int memb=0;
   
   if(argc==2) strcat(fdat,"0000.dat");
   else  
   {  int i;
      for(i=0;i< 4-strlen(argv[2]);i++) strcat(fdat,"0");
       strcat(fdat,argv[2]);
       strcat(fdat,".dat");
       sscanf(argv[2],"%d", &memb); 
   }                
   FILE*Finf=fopen(finfo,"r");
   FILE*Fdat=fopen(fdat,"r");
   if(!Finf || !Fdat){ printf("Can not open files\n"); exit(2);}

   char key[50],val[20];
   

// Parameters which we read form .info file 
   int Index=0;
   char Ref[100]={""};
   char Format[100]={""};
   int  Particle=0;
   int  Adim=0;
   int  AdimV=0;
   double *Al=NULL;
   double *Qs=NULL;
//==============
   int i;

// Read info file    
   while(1==fscanf(Finf,"%[^:]:",key))
   {
//     printf("key='%s'\n", key);
     if(0==strcmp(key,"SetIndex"))        fscanf(Finf,"%d",&Index); 
     else if(0==strcmp(key,"Reference"))  fscanf(Finf,"%[^\n]",Ref);        
     else if(0==strcmp(key,"Format"))     
     {
       fscanf(Finf,"%s",Format); 
       if(0!=strcmp(Format,"lhagrid1")) { printf("Unknown format '%s'\n",Format); exit(3);}
     }
     else if(0==strcmp(key,"AlphaS_Qs")) 
     {  fscanf(Finf,"%*[^[][");
         Adim=0;
         double al;  
         while(1==fscanf(Finf,"%lf",&al))
         {  Qs=(double*) realloc(Qs, (Adim+1)*sizeof(double));
            Qs[Adim++]=al;
            fscanf(Finf," ,");
         } 
     }
     else if(0==strcmp(key,"AlphaS_Vals")) 
     {  fscanf(Finf,"%*[^[][");
        AdimV=0;
        double al;       
         while(1==fscanf(Finf,"%lf",&al))
         {  Al=(double*) realloc(Al, (AdimV+1)*sizeof(double));
            Al[AdimV++]=al;
            fscanf(Finf," ,");
         }
     }
     else if(0==strcmp(key,"Particle")) fscanf(Finf," %d",&Particle);
     char c;
     do  fscanf(Finf,"%c",&c); while(c!='\n');
   }
   fclose(Finf);
   
   printf(" Index=%d\n", Index);
   printf(" Set=%d\n",memb);
   printf(" Particle=%d\n", Particle);
   printf(" Format=%s\n", Format);
   printf(" Reference=%s\n", Ref);
   if(Particle!=2212){ printf("Proton is expected\n"); exit(4);} 
   if(Adim!=AdimV)   { printf("Error: different dimensions for 'AlphaS_Qs' and 'AlphaS_Vals'\n"); exit(5); }
   
   for(int i=1;i<AdimV;i++)
   {  if(Qs[i]==Qs[i-1]) { Adim--; continue;}
      if(Adim!=AdimV) { Qs[i-AdimV+Adim]=Qs[i]; Al[i-AdimV+Adim]=Al[i];} 
   }

//  Read data file        
   for(;;)
   { 
     fscanf(Fdat,"%s ",key);
     if(key[0]=='-') break;
     if(strcmp(key,"Format")==0)
     {  fscanf(Fdat,"%s ",val);
        if(strcmp(val,"lhagrid1")!=0)
        { printf("Unknown LHAPDF fprmat\n");
          exit(3);
        }
     }     
   }


   int Xdim=0,Qdim=0,Pdim=0;
   double *X=NULL,*Q=NULL,*Buff=NULL;
   int    *P=NULL;
   int ip,iq,ix;
   long posB,posE;
   
   
   for(;;)
   { 
     double *dX=NULL,*dQ=NULL;
     int *dP=NULL;
     int dXdim=0,dQdim=0,dPdim=0;
     
     fscanf(Fdat,"%*[ \t\n]");  posB=ftell(Fdat); fscanf(Fdat,"%*[^\n]"); posE=ftell(Fdat); fseek(Fdat,posB,SEEK_SET);
     if(posB==posE) break;  
     for(;ftell(Fdat)<posE;)
     {
       dX=(double*)realloc(dX,(dXdim+1)*sizeof(double));
       fscanf(Fdat,"%lf ", dX+dXdim);
       dXdim++;
     }  

     fscanf(Fdat,"%*[ \t\n]"); posB=ftell(Fdat); fscanf(Fdat,"%*[^\n]"); posE=ftell(Fdat);  fseek(Fdat,posB,SEEK_SET);  
     for(;ftell(Fdat)<posE;)
     { 
       dQ=(double*)realloc(dQ,(dQdim+1)*sizeof(double));
       fscanf(Fdat,"%lf ", dQ+dQdim);
       dQdim++;
     }
     
     fscanf(Fdat,"%*[ \t\n]"); posB=ftell(Fdat);  fscanf(Fdat,"%*[^\n]"); posE=ftell(Fdat); fseek(Fdat,posB,SEEK_SET);
     for(;ftell(Fdat)<posE;)
     { 
         dP=(int*)realloc(dP,(dPdim+1)*sizeof(int));
         fscanf(Fdat,"%d ", dP+dPdim);   
         dPdim++;
     }
     
     double *dBuff=(double*)malloc(dXdim*dQdim*dPdim*sizeof(double));
     double r;

     for(int i=0;i<dXdim*dQdim*dPdim;i++) 
     if(1!=fscanf(Fdat," %lf ",&r)) { printf(" Unexpected end of file\n"); exit(5); } else 
     {  ip=i%dPdim; 
        ix= i/(dPdim*dQdim); 
        iq=(i%(dPdim*dQdim))/dPdim;  
        dBuff[ ip+ ix*dPdim +iq*dPdim*dXdim] =r;  
//               ip+ iq*Pdim +ix*Pdim*Qdim 
     }                                                                                      

     if(!Buff)
     {  Buff=dBuff; X=dX; P=dP; Q=dQ;  
        Xdim=dXdim; Pdim=dPdim; Qdim=dQdim;        
     } else 
     {  if(dXdim!=Xdim) { printf("Unsupported option:  X-grid dimension  is changed\n"); return 3;}
        for(int i=0;i<Xdim;i++) if(X[i]!=dX[i]) { printf("Unsupported option:  X-grid is changed. position %d %E -> %E\n",i+1, X[i],dX[i]); return 3;}
        if(dPdim!=Pdim) { printf("Unsupported option:  parton list  is changed\n"); return 3;}
        for(int i=0;i<Pdim;i++) if(P[i]!=dP[i]) { printf("Unsupported option: parton list is changed\n"); return 3;}
        if(Q[Qdim-1]!=dQ[0]) { printf("Unsupported option: no matching for Q sets\n"); return 3;}
        Q=realloc(Q,(Qdim+dQdim-1)*sizeof(double)); 
        for(int i=0; i<dQdim-1;i++) Q[i+Qdim]=dQ[i+1];
        
        Buff=realloc(Buff, Pdim*Xdim*(Qdim+dQdim-1)*sizeof(double)); 
        for(int i=0;i<Pdim*Xdim*(dQdim-1);i++) Buff[Pdim*Xdim*Qdim+i]=dBuff[Pdim*Xdim+i];
        Qdim+=dQdim-1;
        free(dP),free(dQ); free(dX);free(dBuff);       
     }    
     fscanf(Fdat," %s ",key);
//     printf("last key=%s ftell=%d\n",key,ftell(Fdat));
   }
   fclose(Fdat);
   

//  Output 

//#ifdef QQQ
   if(argc>2)
   { char*ch_=malloc(strlen(ch)+strlen(argv[2])+2);
     sprintf(ch_,"%s:%s",ch,argv[2]);
     ch=ch_; 
   }
   
   char*fout=malloc(strlen(ch)+5);
   sprintf(fout,"%s.pdt",ch);   
   FILE*Fout=fopen(fout,"w");

   fprintf(Fout,"\n#distribution \"%s(proton)\"       2212 => ", ch);
   posB=ftell(Fout); for(int i=0;i<Pdim;i++)  fprintf(Fout,"    ");    
   fprintf(Fout,"\n#distribution \"%s(anti-proton)\" -2212 => ", ch);
   posE=ftell(Fout); for(int i=0;i<Pdim;i++)  fprintf(Fout,"    ");

   fprintf(Fout,"\n#Index %d",Index);
   fprintf(Fout,"\n#Memb %d",memb);
   fprintf(Fout,"\n#Source  LHAPDF6");
   fprintf(Fout,"\n#Reference %s",Ref); 
   fprintf(Fout,"\n#Interpolation biCubicLogXQ");
   
   fprintf(Fout,"\n#Q_grid\n");
   for(iq=0;iq<Qdim;iq++)
   { fprintf(Fout," %.4E",Q[iq]);
     if( (iq+1)%10==0) fprintf(Fout,"\n");
   }

   fprintf(Fout,"\n#Alpha\n");
   for(iq=0;iq<Qdim;iq++)
   { fprintf(Fout," %.4E",  polint3(Q[iq] , Adim, Qs, Al));
     if( (iq+1)%10==0) fprintf(Fout,"\n");
   }

   fprintf(Fout,"\n#X_grid\n");
   for(ix=0;ix<Xdim;ix++)
   { fprintf(Fout," %.4E",X[ix]);
     if( (ix+1)%10==0) fprintf(Fout,"\n");
   }
   
   char * Pstr=(char*)malloc(Pdim*4); Pstr[0]=0;
   char *aPstr=(char*)malloc(Pdim*4);aPstr[0]=0;
   int pk;
   for(ip=0,pk=0;ip<Pdim;ip++)
   {  int wrt=1;
      switch(P[ip])
      { case 21: case 22: case 23: 
            sprintf(Pstr+strlen(Pstr)," %d",P[ip]);
            sprintf(aPstr+strlen(aPstr)," %d",P[ip]);
            break;
        case 3: case 4: case 5: case 6: case 24:
            sprintf(Pstr+strlen(Pstr)," (%d %d)",P[ip],-P[ip]);
            sprintf(aPstr+strlen(aPstr)," (%d %d)",P[ip],-P[ip]);
            break;
        case -2: case -1: case 1: case 2: 
            sprintf(Pstr+strlen(Pstr)," %d",P[ip]);
            sprintf(aPstr+strlen(aPstr)," %d",-P[ip]);
            break;
        default : wrt=0;    
      }
      if(!wrt) continue;
      pk++; 
      fprintf(Fout,"\n#%d-parton\n",pk);
      for(iq=0;iq<Qdim;iq++)
      {  double d;
//         for(ix=0;ix<Xdim;ix++) { fprintf(Fout," %.4E",buff[ip+ iq*Pdim +ix*Pdim*Qdim] );} 
            for(ix=0;ix<Xdim;ix++) { fprintf(Fout," %.4E",Buff[ip+ ix*Pdim +iq*Pdim*Xdim] );} 
         fprintf(Fout,"\n");
      }    
   } 
   fseek(Fout,posB,SEEK_SET); fprintf(Fout,"%s",Pstr); 
   fseek(Fout,posE,SEEK_SET); fprintf(Fout,"%s",aPstr); 
        
   fclose(Fout);
   return 0;    
}
