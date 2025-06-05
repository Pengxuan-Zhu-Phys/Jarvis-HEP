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

#include "gm2calc/SM.h"
#include "gm2calc/SM.hpp"
#include <complex>

extern "C"
{

/**
 * @brief Set SM prameters to default values
 */
void gm2calc_sm_set_to_default(gm2calc_SM* sm)
{
   if (sm == nullptr) {
      return;
   }

   gm2calc::SM def;
   sm->alpha_em_0 = def.get_alpha_em_0();
   sm->alpha_em_mz = def.get_alpha_em_mz();
   sm->alpha_s_mz = def.get_alpha_s_mz();
   sm->mh = def.get_mh();
   sm->mw = def.get_mw();
   sm->mz = def.get_mz();
   for (int i = 0; i < 3; i++) {
      sm->mu[i] = def.get_mu(i);
   }
   for (int i = 0; i < 3; i++) {
      sm->md[i] = def.get_md(i);
   }
   for (int i = 0; i < 3; i++) {
      sm->mv[i] = def.get_mv(i);
   }
   for (int i = 0; i < 3; i++) {
      sm->ml[i] = def.get_ml(i);
   }
   for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
         sm->ckm_real[i][k] = std::real(def.get_ckm(i, k));
         sm->ckm_imag[i][k] = std::imag(def.get_ckm(i, k));
      }
   }
}

} // extern "C"
