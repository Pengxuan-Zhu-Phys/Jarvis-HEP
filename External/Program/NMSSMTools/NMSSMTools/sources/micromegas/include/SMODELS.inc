{
   char *VERSION="2.0.0";
   char smodels_database[100];
   char smodels_version[100];
   system("rm -f smodels.*");

   strcpy(smodelsInfo,VERSION);
   topology[0]=0;
   analysis[0]=0;
   status=0;
   Rvalue=0;
   
   if(VERSION[0]<'2' || VERSION[0]>'9')
   {
      char command[100];
      sprintf(command, "echo  \" SModelS '%s' is not supported\n\" > smodels.err",VERSION);
      system(command);
      status=-5;
   }else
   {        
      char*command=(char*)malloc(strlen(micrO)+300 );
      sprintf(command,"%s/Packages/smodels-%s/runSModelS.py",micrO,VERSION);
      if(access(command,X_OK ))
      {  printf("Downloading the program SModels\n");
         sprintf(command,"pwd=$PWD; cd  %s/Packages;  make -f  SMODELS.makef VERSION=%s 1>$pwd/smodels.err  2>$pwd/smodels.err", micrO,VERSION); 
         system(command);
         sprintf(command,"%s/Packages/smodels-%s/runSModelS.py",micrO,VERSION);
         if(access(command,X_OK )) status=-4; else system("rm -f smodels.err");
      } 

      if(status==0) if(smodels(LHCrun,5/*nf*/, 0.005/*fb*/, "smodels.slha",VERSION, 0)) status=-2;
      if(status==0)
      {    
         sprintf(command,"micrO=%s; python3  $micrO/Packages/smodels-%s/runSModelS.py -f smodels.slha -p $micrO/include/smodels_parameters.ini -o ./ -v error 2>smodels.err\n"
                             "if(test $? -eq 0) then rm -f smodels.err; fi \n",micrO,VERSION);
         system(command);
         if(!access("smodels.err",R_OK )) status=-3;
         else 
         {        
            err=slhaRead("smodels.slha.smodelsslha", 1);
            slhaSTRFormat("SModelS_Settings"," 0 %s",smodels_version);
            slhaSTRFormat("SModelS_Settings"," 1 %s",smodels_database);
            sprintf(smodelsInfo,"version %s with database %s",smodels_version, smodels_database);
            status=slhaVal("SModelS_Exclusion",0.,2,0,0); 
            if(status>=0)
            { 
               slhaSTRFormat("SModelS_Exclusion"," 1 0  %[^#]",topology);
               Rvalue=slhaVal("SModelS_Exclusion",0.,2,1,1);
               slhaSTRFormat("SModelS_Exclusion"," 1 4 %s",analysis);
            }
         }
      }
      free(command);
   }
   if(status>=-2) smodelsOK=1; else smodelsOK=0; 
}
