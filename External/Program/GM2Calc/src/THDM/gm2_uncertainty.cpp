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

#include "gm2calc/gm2_uncertainty.hpp"
#include "gm2calc/gm2_1loop.hpp"
#include "gm2calc/gm2_2loop.hpp"
#include "gm2calc/THDM.hpp"
#include "gm2_1loop_helpers.hpp"
#include "gm2_uncertainty_helpers.hpp"

#include <cmath>

/**
 * \file gm2_uncertainty.cpp
 *
 * Contains functions necessary to calculate theory uncertainty
 * associated with amu.
 */

namespace gm2calc {

/**
 * Calculates uncertainty associated with amu(0-loop)
 *
 * The estimated uncertainty is the magnitude amu(1-loop) plus
 * amu(2-loop).
 *
 * @param model model parameters (unused tag type)
 * @param amu_1L 1-loop contribution to amu
 * @param amu_2L 2-loop contribution to amu
 *
 * @return uncertainty for amu(0-loop)
 */
double calculate_uncertainty_amu_0loop(const THDM& /* model */, double amu_1L, double amu_2L)
{
   return std::abs(amu_1L) + std::abs(amu_2L);
}

/**
 * Calculates uncertainty associated with amu(0-loop)
 *
 * @param model model parameters
 *
 * @return uncertainty for amu(0-loop)
 */
double calculate_uncertainty_amu_0loop(const THDM& model)
{
   const double amu_1L = calculate_amu_1loop(model);
   const double amu_2L = calculate_amu_2loop(model);

   return calculate_uncertainty_amu_0loop(model, amu_1L, amu_2L);
}

/**
 * Calculates uncertainty associated with amu(1-loop)
 *
 * The estimated uncertainty is the sum of magnitude amu(2-loop) and
 * the 2-loop uncertainty, calculated by
 * calculate_uncertainty_amu_2loop().
 *
 * @param model model parameters
 * @param amu_1L 1-loop contribution to amu
 * @param amu_2L 2-loop contribution to amu
 *
 * @return uncertainty for amu(1-loop)
 */
double calculate_uncertainty_amu_1loop(const THDM& model, double amu_1L, double amu_2L)
{
   const double delta_amu_2L = calculate_uncertainty_amu_2loop(model, amu_1L, amu_2L);

   return std::abs(amu_2L) + std::abs(delta_amu_2L);
}

/**
 * Calculates uncertainty associated with amu(1-loop)
 *
 * @param model model parameters
 *
 * @return uncertainty for amu(1-loop)
 */
double calculate_uncertainty_amu_1loop(const THDM& model)
{
   const double amu_1L = calculate_amu_1loop(model);
   const double amu_2L = calculate_amu_2loop(model);

   return calculate_uncertainty_amu_1loop(model, amu_1L, amu_2L);
}

/**
 * Calculates uncertainty associated with amu(2-loop)
 *
 * Takes into account the neglected two-loop contribution
 * \f$a_\mu^{\Delta r\text{-shift}}\f$, Eq.(34) arxiv:1607.06292,
 * two-loop contributions of \f$O(m_\mu^4)\f$ and three-loop
 * contributions of \f$O(m_\mu^2)\f$
 *
 * @param model model parameters
 * @param amu_1L 1-loop contribution to amu
 * @param amu_2L 2-loop contribution to amu
 *
 * @return uncertainty for amu(2-loop)
 */
double calculate_uncertainty_amu_2loop(const THDM& model, double amu_1L, double amu_2L)
{
   const double pi = 3.1415926535897932;
   const double alpha_em = model.get_alpha_em();
   const double mm = model.get_MFe(1);
   const double mH = std::abs(model.get_Mhh(1));
   const double mA = std::abs(model.get_MAh(1));
   const double mHp = std::abs(model.get_MHm(1));
   const double mNP = std::fmin(mH,std::fmin(mA,mHp)); // new physics scale

   // universal 2-loop QED logarithmic correction from Eq.(51) hep-ph/9803384
   const double delta_alpha_em = -4*alpha_em/pi*std::log(std::abs(mNP/mm));

   // upper bound on neglected 2-loop contribution from \Delta r
   const double delta_amu_2L_delta_r = 2e-12;

   // estimate of 2-loop corrections O(mm^4)
   const double delta_amu_2L_mm4 = std::abs(amu_1L*delta_alpha_em);

   // estimate of 3-loop corrections O(mm^2)
   const double delta_amu_3L = std::abs(amu_2L*delta_alpha_em);

   return delta_amu_2L_delta_r + delta_amu_2L_mm4 + delta_amu_3L;
}

/**
 * Calculates uncertainty associated with amu(2-loop)
 *
 * @param model model parameters
 *
 * @return uncertainty for amu(2-loop)
 */
double calculate_uncertainty_amu_2loop(const THDM& model)
{
   const double amu_1L = calculate_amu_1loop(model);
   const double amu_2L = calculate_amu_2loop(model);

   return calculate_uncertainty_amu_2loop(model, amu_1L, amu_2L);
}

} // namespace gm2calc
