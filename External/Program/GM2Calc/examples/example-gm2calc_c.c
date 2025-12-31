#include "gm2calc/gm2_1loop.h"
#include "gm2calc/gm2_2loop.h"
#include "gm2calc/gm2_uncertainty.h"
#include "gm2calc/MSSMNoFV_onshell.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void setup(MSSMNoFV_onshell* model) {
   const double pi = 3.14159265358979323846;

   /* fill SM parameters */
   gm2calc_mssmnofv_set_alpha_MZ(model,0.0077552);         /* 1L */
   gm2calc_mssmnofv_set_alpha_thompson(model,0.00729735);  /* 2L */
   gm2calc_mssmnofv_set_g3(model,sqrt(4 * pi * 0.1184));   /* 2L */
   gm2calc_mssmnofv_set_MT_pole(model,173.34);             /* 2L */
   gm2calc_mssmnofv_set_MB_running(model,4.18);            /* 2L, mb(mb) MS-bar */
   gm2calc_mssmnofv_set_MM_pole(model,0.1056583715);       /* 1L */
   gm2calc_mssmnofv_set_ML_pole(model,1.777);              /* 2L */
   gm2calc_mssmnofv_set_MW_pole(model,80.385);             /* 1L */
   gm2calc_mssmnofv_set_MZ_pole(model,91.1876);            /* 1L */

   /* fill DR-bar parameters */
   gm2calc_mssmnofv_set_TB(model, 10);      /* 1L */
   gm2calc_mssmnofv_set_Ae(model,1,1,0);    /* 1L */

   /* fill on-shell parameters */
   gm2calc_mssmnofv_set_Mu(model,350);      /* 1L */
   gm2calc_mssmnofv_set_MassB(model,150);   /* 1L */
   gm2calc_mssmnofv_set_MassWB(model,300);  /* 1L */
   gm2calc_mssmnofv_set_MassG(model,1000);  /* 2L */
   gm2calc_mssmnofv_set_Au(model,2,2,0);    /* 2L */
   gm2calc_mssmnofv_set_Ad(model,2,2,0);    /* 2L */
   gm2calc_mssmnofv_set_Ae(model,2,2,0);    /* 2L */
   gm2calc_mssmnofv_set_MAh_pole(model,1500);    /* 2L */
   gm2calc_mssmnofv_set_scale(model,454.7); /* 2L */

   for (unsigned i = 0; i < 3; i++) {
      gm2calc_mssmnofv_set_mq2(model,i,i,500*500); /* 2L */
      gm2calc_mssmnofv_set_ml2(model,i,i,500*500); /* 1L(smuon)/2L */
      gm2calc_mssmnofv_set_md2(model,i,i,500*500); /* 2L */
      gm2calc_mssmnofv_set_mu2(model,i,i,500*500); /* 2L */
      gm2calc_mssmnofv_set_me2(model,i,i,500*500); /* 1L(smuon)/2L */
   }

   /* calculate mass spectrum */
   const gm2calc_error error = gm2calc_mssmnofv_calculate_masses(model);

   if (gm2calc_mssmnofv_have_warning(model)) {
      char warning[400];
      gm2calc_mssmnofv_get_warnings(model, warning, sizeof(warning));
      printf("Warning: %s\n", warning);
   }

   if (error != gm2calc_NoError) {
      printf("Error: %s\n", gm2calc_error_str(error));
      abort();
   }
}

int main() {
   MSSMNoFV_onshell* model = gm2calc_mssmnofv_new();

   setup(model);

   const double amu =
      + gm2calc_mssmnofv_calculate_amu_1loop(model)
      + gm2calc_mssmnofv_calculate_amu_2loop(model);

   const double delta_amu =
      gm2calc_mssmnofv_calculate_uncertainty_amu_2loop(model);

   printf("amu = %e +- %e\n", amu, delta_amu);

   /* destroy model to prevent resource leak */
   gm2calc_mssmnofv_free(model);

   return 0;
}
