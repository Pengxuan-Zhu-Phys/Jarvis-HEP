/* ====================================================================
 * This file is part of GM2Calc.
 *
 * GM2Calc is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * GM2Calc is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GM2Calc.  If not, see
 * <http://www.gnu.org/licenses/>.
 * ==================================================================== */

#include "gm2calc/gm2_1loop.h"
#include "gm2calc/gm2_1loop.hpp"

#include <limits>

/**
 * @file gm2_1loop_c.cpp
 * @brief contains definitions of C interface functions for 1-loop calculation
 *
 * This file contains the definitions for the C interface functions
 * used to calculate \f$a_\mu\f$ at the 1-loop level.
 */

extern "C" {

/** calculates full 1-loop SUSY contributions to (g-2) in the MSSM (w/ tan(beta) resummation) */
double gm2calc_mssmnofv_calculate_amu_1loop(const MSSMNoFV_onshell* model)
{
   double amu = std::numeric_limits<double>::quiet_NaN();

   try {
      amu = gm2calc::calculate_amu_1loop(
         *reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model));
   } catch (...) {}

   return amu;
}

/** calculates full 1-loop SUSY contributions to (g-2) in the MSSM (no tan(beta) resummation) */
double gm2calc_mssmnofv_calculate_amu_1loop_non_tan_beta_resummed(const MSSMNoFV_onshell* model)
{
   double amu = std::numeric_limits<double>::quiet_NaN();

   try {
      amu = gm2calc::calculate_amu_1loop_non_tan_beta_resummed(
         *reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model));
   } catch (...) {}

   return amu;
}

/* === routines for individual 1-loop contributions === */

/** 1-loop neutralino contribution */
double gm2calc_mssmnofv_amu1LChi0(const MSSMNoFV_onshell* model)
{
   return gm2calc::amu1LChi0(
      *reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model));
}

/** 1-loop chargino contribution */
double gm2calc_mssmnofv_amu1LChipm(const MSSMNoFV_onshell* model)
{
   return gm2calc::amu1LChipm(
      *reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model));
}

} /* extern "C" */
