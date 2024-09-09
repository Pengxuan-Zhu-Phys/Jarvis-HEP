{ 
//   int n1,n2;
//   pythonversion_(&n1,&n2); 
//   if(n1==2)
   {
   
      char*command=(char*)malloc(strlen(micrO)+200 );
      sprintf(command, "python %s/Packages/Lilith-2-2.1db19.09patch1/run_lilith.py Lilith_in.xml  -s -r  Lilith_out.slha %d  %.2E ",
        micrO,  n_par, m2logL_reference);

      system(command);
      free(command);

      slhaRead("Lilith_out.slha", 1);
      m2logL = slhaVal("LilithResults",0.,1,0);
      exp_ndf = slhaVal("LilithResults",0.,1,1);
      m2logL_reference = slhaVal("LilithResults",0.,1,2);
      ndf = slhaVal("LilithResults",0.,1,3);
      pvalue = slhaVal("LilithResults",0.,1,4);
      slhaSTRFormat("LilithResults","5 %s ",Lilith_version);
    }
//     else 
//    {  printf("LILITH can not work with Python version > 2\n");
//       ndf=0;
//    }   
}
