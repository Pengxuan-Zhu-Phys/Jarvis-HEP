#define HBversion  "5.5.0" 
#define HSversion  "2.3.0"  
{
   system("rm -f HB.out HS.out");
   system("cp HB.in HS.in");
   system("cp HB.in hb.in");
   
   int len=strlen(micrO)+150;
   char*HB=(char*)malloc(len);
   char*HS=(char*)malloc(len);
   sprintf(HB,"%s/Packages/higgsbounds-%s/build/HiggsBounds", micrO,HBversion);
   sprintf(HS,"%s/Packages/higgssignals-%s/build/HiggsSignals",micrO,HSversion);    

#ifdef HIGGSBOUNDS

   if(access(HB ,X_OK ))
   {  char*command=(char*)malloc(len);
      sprintf(command,"cd %s/Packages; make -f HBOUNDS.makef VERSION=%s", micrO,HBversion);
      system(command);
      free(command);
   } 
   if(!access(HB ,X_OK ))
   {
     sprintf(HB+strlen(HB)," LandH SLHA %d  %d HB.in HB.out > hb.stdout",NH0,NHch);
     system(HB);

     slhaRead("HB.out",1+4);
     slhaSTRFormat("HiggsBoundsResults","0 %s %*s",HB_version);
     for(int i=0;i<3;i++) 
     {  char fmt[20];
        sprintf(fmt,"%d 1 %%lf",i+1);
        HB_id[i]=slhaValFormat("HiggsBoundsResults",0., fmt)+0.001; 
        sprintf(fmt,"%d 2 %%lf",i+1); 
        HB_result[i]=slhaValFormat("HiggsBoundsResults",0., fmt)+0.001;
        sprintf(fmt,"%d 3 %%lf",i+1);
        HB_obsratio[i]=slhaValFormat("HiggsBoundsResults",0.,fmt );
        sprintf(fmt,"%d 5 ||%%[^|]||",i+1);
        slhaSTRFormat("HiggsBoundsResults",fmt,HB_chan[i]);
     }
   }  
#endif 

#ifdef HIGGSSIGNALS
   if(access(HS ,X_OK ))
   {  char*command=(char*)malloc(len);
      sprintf(command,"cd %s/Packages; make -f HSIGNALS.makef VERSION=%s", micrO,HSversion);
      system(command);
      free(command);
   }   
   
   if(!access(HS ,X_OK ))
   {                           
#define DataSet " latestresults "
#define Method  " peak " 
//#define  Method " mass "
//#define  Method " both "
#define PDF  " 2 "  // Gaussian
//#define PDF " 1 "  // box 
//#define PDF " 3 "  // box+Gaussia

     sprintf(HS+strlen(HS)," " DataSet Method  PDF  " SLHA %d %d HS.in > hs.stdout",NH0,NHch); 
     system(HS);

     system("grep -A 10000  HiggsSignalsResults HS.in > HS.out");
     slhaRead("HS.out",1+4);
     slhaSTRFormat("HiggsSignalsResults"," 0 ||%[^|]||",HS_version);
     HS_observ=slhaVal("HiggsSignalsResults",0.,1,8);
     HS_chi2=slhaVal("HiggsSignalsResults",0.,1,17);
     HS_pval=slhaVal("HiggsSignalsResults",0.,1,21);
#undef PDF
#undef Method
#undef DataSet
   }
#endif
free(HB); free(HS);

}