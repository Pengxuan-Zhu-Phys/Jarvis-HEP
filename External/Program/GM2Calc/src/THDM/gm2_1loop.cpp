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
#include "gm2calc/THDM.hpp"
#include "THDM/gm2_1loop_helpers.hpp"
#include <cmath>

/**
 * \file gm2_1loop.cpp
 *
 * Contains functions necessary to calculate the THDM
 * contributions for g-2 at the 1-loop level.
 */

namespace gm2calc {

/**
 * Calculates full 1-loop contribution to a_mu in the general THDM.
 *
 * @param model THDM model parameters, masses and mixings
 * @return 1-loop contribution to a_mu
 */
double calculate_amu_1loop(const THDM& model)
{
   thdm::THDM_1L_parameters pars;
   pars.alpha_em = model.get_alpha_em();
   pars.mm = model.get_MFe(1);
   pars.mw = model.get_MVWm();
   pars.mz = model.get_MVZ();
   pars.mhSM = model.get_sm().get_mh();
   pars.mA = model.get_MAh(1);
   pars.mHp = model.get_MHm(1);
   pars.ml = model.get_MFe();
   pars.mv = model.get_MFv();
   pars.mh = model.get_Mhh();
   pars.ylh = model.get_ylh();
   pars.ylH = model.get_ylH();
   pars.ylA = model.get_ylA();
   pars.ylHp = model.get_ylHp();

   return thdm::amu1L(pars);
}

} // namespace gm2calc
