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

#include "gm2calc/gm2_2loop.hpp"
#include "gm2calc/THDM.hpp"
#include "THDM/gm2_2loop_helpers.hpp"

namespace gm2calc {

/**
 * Calculates 2-loop bosonic contribution to a_mu in the THDM.
 *
 * @param model THDM model parameters, masses and mixings
 * @return 2-loop contribution to a_mu
 */
double calculate_amu_2loop_bosonic(const THDM& model)
{
   thdm::THDM_B_parameters pars_b;
   pars_b.alpha_em = model.get_alpha_em();
   pars_b.mm = model.get_MFe(1);
   pars_b.mw = model.get_MVWm();
   pars_b.mz = model.get_MVZ();
   pars_b.mhSM = model.get_sm().get_mh();
   pars_b.mA = model.get_MAh(1);
   pars_b.mHp = model.get_MHm(1);
   pars_b.mh = model.get_Mhh();
   pars_b.tb = model.get_tan_beta();
   pars_b.zetal = model.get_zeta_l();
   pars_b.cos_beta_minus_alpha = model.get_cos_beta_minus_alpha();
   pars_b.lambda5 = model.get_LambdaFive();
   pars_b.lambda67 = model.get_LambdaSixSeven();

   return amu2L_B(pars_b);
}

/**
 * Calculates fermionic 2-loop contribution to a_mu in the THDM.
 *
 * @param model THDM model parameters, masses and mixings
 * @return 2-loop contribution to a_mu
 */
double calculate_amu_2loop_fermionic(const THDM& model)
{
   thdm::THDM_F_parameters pars_f;
   pars_f.alpha_em = model.get_alpha_em();
   pars_f.mm = model.get_MFe(1);
   pars_f.mw = model.get_MVWm();
   pars_f.mz = model.get_MVZ();
   pars_f.mhSM = model.get_sm().get_mh();
   pars_f.mA = model.get_MAh(1);
   pars_f.mHp = model.get_MHm(1);
   pars_f.mh = model.get_Mhh();
   pars_f.ml = model.get_MFe();
   pars_f.mu = model.get_MFu();
   pars_f.md = model.get_MFd();
   pars_f.yuh = model.get_yuh();
   pars_f.yuH = model.get_yuH();
   pars_f.yuA = model.get_yuA();
   pars_f.yuHp = model.get_yuHp();
   pars_f.ydh = model.get_ydh();
   pars_f.ydH = model.get_ydH();
   pars_f.ydA = model.get_ydA();
   pars_f.ydHp = model.get_ydHp();
   pars_f.ylh = model.get_ylh();
   pars_f.ylH = model.get_ylH();
   pars_f.ylA = model.get_ylA();
   pars_f.ylHp = model.get_ylHp();
   pars_f.vckm = model.get_sm().get_ckm();

   return amu2L_F(pars_f);
}

/**
 * Calculates full 2-loop contribution to a_mu in the general THDM.
 *
 * @param model THDM model parameters, masses and mixings
 * @return 2-loop contribution to a_mu
 */
double calculate_amu_2loop(const THDM& model)
{
   return calculate_amu_2loop_bosonic(model)
      + calculate_amu_2loop_fermionic(model);
}

} // namespace gm2calc
