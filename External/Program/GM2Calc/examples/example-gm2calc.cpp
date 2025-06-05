#include "gm2calc/gm2_1loop.hpp"
#include "gm2calc/gm2_2loop.hpp"
#include "gm2calc/gm2_uncertainty.hpp"
#include "gm2calc/gm2_error.hpp"
#include "gm2calc/MSSMNoFV_onshell.hpp"

#include <cstdio>
#include <cmath>

gm2calc::MSSMNoFV_onshell setup() {
   gm2calc::MSSMNoFV_onshell model;

   const double Pi = 3.141592653589793;
   const Eigen::Matrix<double,3,3> UnitMatrix
      = Eigen::Matrix<double,3,3>::Identity();

   // fill SM parameters
   model.set_alpha_MZ(0.0077552);               // 1L
   model.set_alpha_thompson(0.00729735);        // 2L
   model.set_g3(std::sqrt(4 * Pi * 0.1184));    // 2L
   model.get_physical().MFt   = 173.34;         // 2L
   model.get_physical().MFb   = 4.18;           // 2L, mb(mb) MS-bar
   model.get_physical().MFm   = 0.1056583715;   // 1L
   model.get_physical().MFtau = 1.777;          // 2L
   model.get_physical().MVWm  = 80.385;         // 1L
   model.get_physical().MVZ   = 91.1876;        // 1L

   // fill DR-bar parameters
   model.set_TB(10);                      // 1L
   model.set_Ae(1,1,0);                   // 1L

   // fill on-shell parameters
   model.set_Mu(350);                     // 1L
   model.set_MassB(150);                  // 1L
   model.set_MassWB(300);                 // 1L
   model.set_MassG(1000);                 // 2L
   model.set_mq2(500 * 500 * UnitMatrix); // 2L
   model.set_ml2(500 * 500 * UnitMatrix); // 1L(smuon)/2L
   model.set_md2(500 * 500 * UnitMatrix); // 2L
   model.set_mu2(500 * 500 * UnitMatrix); // 2L
   model.set_me2(500 * 500 * UnitMatrix); // 1L(smuon)/2L
   model.set_Au(2,2,0);                   // 2L
   model.set_Ad(2,2,0);                   // 2L
   model.set_Ae(2,2,0);                   // 2L
   model.set_MA0(1500);                   // 2L
   model.set_scale(454.7);                // 2L

   // calculate mass spectrum
   model.calculate_masses();

   return model;
}

int main() {
   try {
      gm2calc::MSSMNoFV_onshell model(setup());

      const double amu =
         + gm2calc::calculate_amu_1loop(model)
         + gm2calc::calculate_amu_2loop(model);

      const double delta_amu =
         gm2calc::calculate_uncertainty_amu_2loop(model);

      std::printf("amu = %.5e +- %.5e\n", amu, delta_amu);
   } catch (const gm2calc::Error& e) {
      std::printf("%s\n", e.what());
   }

   return 0;
}
