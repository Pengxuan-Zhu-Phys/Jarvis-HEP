// ====================================================================
// This file is part of GM2Calc.
//
// GM2Calc is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// GM2Calc is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GM2Calc.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#include "gm2calc/gm2_1loop.hpp"
#include "gm2calc/gm2_2loop.hpp"
#include "gm2calc/gm2_uncertainty.hpp"
#include "gm2calc/gm2_error.hpp"
#include "gm2calc/MSSMNoFV_onshell.hpp"

#include <cstdio>
#include <iostream>
#include <string>

gm2calc::MSSMNoFV_onshell setup()
{
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

   return model;
}

int main()
{
   const double tanb_start = 2.;
   const double tanb_stop = 100.;
   const unsigned nsteps = 100;

   printf("# %14s %16s %16s %16s\n",
          "tan(beta)", "amu", "uncertainty", "error");

   for (unsigned n = 0; n < nsteps; n++) {
      double amu{0.0}, delta_amu{0.0};
      const double tanb = tanb_start + (tanb_stop - tanb_start) * n / nsteps;
      std::string error;

      gm2calc::MSSMNoFV_onshell model(setup());
      model.set_TB(tanb);

      try {
         model.calculate_masses();
         amu = gm2calc::calculate_amu_1loop(model)
             + gm2calc::calculate_amu_2loop(model);
         delta_amu = gm2calc::calculate_uncertainty_amu_2loop(model);
      } catch (const gm2calc::Error& e) {
         error = "# " + std::string(e.what());
         amu = delta_amu = std::numeric_limits<double>::signaling_NaN();
      }

      printf("%16.8e %16.8e %16.8e %s\n",
             tanb, amu, delta_amu, error.c_str());
   }

   return 0;
}
