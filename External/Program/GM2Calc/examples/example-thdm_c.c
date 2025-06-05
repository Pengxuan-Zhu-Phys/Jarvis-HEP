#include "gm2calc/gm2_1loop.h"
#include "gm2calc/gm2_2loop.h"
#include "gm2calc/gm2_uncertainty.h"
#include "gm2calc/THDM.h"
#include "gm2calc/SM.h"

#include <stdio.h>

int main()
{
   const gm2calc_THDM_mass_basis basis = {
      .yukawa_type = gm2calc_THDM_type_2,
      .mh = 125,
      .mH = 400,
      .mA = 420,
      .mHp = 440,
      .sin_beta_minus_alpha = 0.999,
      .lambda_6 = 0,
      .lambda_7 = 0,
      .tan_beta = 3,
      .m122 = 40000,
      .zeta_u = 0,
      .zeta_d = 0,
      .zeta_l = 0,
      .Delta_u = { {0,0,0}, {0,0,0}, {0,0,0} },
      .Delta_d = { {0,0,0}, {0,0,0}, {0,0,0} },
      .Delta_l = { {0,0,0}, {0,0,0}, {0,0,0} },
      .Pi_u = { {0,0,0}, {0,0,0}, {0,0,0} },
      .Pi_d = { {0,0,0}, {0,0,0}, {0,0,0} },
      .Pi_l = { {0,0,0}, {0,0,0}, {0,0,0} }
   };

   gm2calc_SM sm;
   gm2calc_sm_set_to_default(&sm);
   sm.alpha_em_mz = 1.0/128.94579;
   sm.mu[2] = 173.34;
   sm.mu[1] = 1.28;
   sm.md[2] = 4.18;
   sm.ml[2] = 1.77684;

   gm2calc_THDM_config config;
   gm2calc_thdm_config_set_to_default(&config);

   gm2calc_THDM* model = 0;
   gm2calc_error error = gm2calc_thdm_new_with_mass_basis(&model, &basis, &sm, &config);

   if (error == gm2calc_NoError) {
      const double amu = gm2calc_thdm_calculate_amu_1loop(model)
                       + gm2calc_thdm_calculate_amu_2loop(model);

      const double delta_amu =
         gm2calc_thdm_calculate_uncertainty_amu_2loop(model);

      printf("amu = %g +- %g\n", amu, delta_amu);
   } else {
      printf("Error: %s\n", gm2calc_error_str(error));
   }

   gm2calc_thdm_free(model);

   return 0;
}
