#include "gm2calc/gm2_1loop.hpp"
#include "gm2calc/gm2_2loop.hpp"
#include "gm2calc/gm2_uncertainty.hpp"
#include "gm2calc/gm2_error.hpp"
#include "gm2calc/THDM.hpp"

#include <cstdio>

int main()
{
   gm2calc::thdm::Mass_basis basis;
   basis.yukawa_type = gm2calc::thdm::Yukawa_type::type_2;
   basis.mh = 125;
   basis.mH = 400;
   basis.mA = 420;
   basis.mHp = 440;
   basis.sin_beta_minus_alpha = 0.999;
   basis.lambda_6 = 0;
   basis.lambda_7 = 0;
   basis.tan_beta = 3;
   basis.m122 = 40000;
   basis.zeta_u = 0;
   basis.zeta_d = 0;
   basis.zeta_l = 0;
   basis.Delta_u << 0, 0, 0, 0, 0, 0, 0, 0, 0;
   basis.Delta_d << 0, 0, 0, 0, 0, 0, 0, 0, 0;
   basis.Delta_l << 0, 0, 0, 0, 0, 0, 0, 0, 0;
   basis.Pi_u << 0, 0, 0, 0, 0, 0, 0, 0, 0;
   basis.Pi_d << 0, 0, 0, 0, 0, 0, 0, 0, 0;
   basis.Pi_l << 0, 0, 0, 0, 0, 0, 0, 0, 0;

   gm2calc::SM sm;
   sm.set_alpha_em_mz(1.0/128.94579);
   sm.set_mu(2, 173.34);
   sm.set_mu(1, 1.28);
   sm.set_md(2, 4.18);
   sm.set_ml(2, 1.77684);

   gm2calc::thdm::Config config;
   config.force_output = false;
   config.running_couplings = true;

   try {
      const gm2calc::THDM model(basis, sm, config);

      const double amu = gm2calc::calculate_amu_1loop(model)
                       + gm2calc::calculate_amu_2loop(model);

      const double delta_amu =
         gm2calc::calculate_uncertainty_amu_2loop(model);

      std::printf("amu = %.5e +- %.5e\n", amu, delta_amu);
   } catch (const gm2calc::Error& e) {
      std::printf("%s\n", e.what());
   }

   return 0;
}
