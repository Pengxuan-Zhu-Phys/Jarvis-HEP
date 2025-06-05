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

/** calculates full 1-loop BSM contributions to (g-2) in the THDM */
double gm2calc_thdm_calculate_amu_1loop(const gm2calc_THDM* model)
{
   double amu = std::numeric_limits<double>::quiet_NaN();

   try {
      amu = gm2calc::calculate_amu_1loop(
         *reinterpret_cast<const gm2calc::THDM*>(model));
   } catch (...) {}

   return amu;
}

} /* extern "C" */
