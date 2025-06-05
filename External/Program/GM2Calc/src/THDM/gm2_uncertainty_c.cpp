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

#include "gm2calc/gm2_uncertainty.h"
#include "gm2calc/gm2_uncertainty.hpp"
#include "gm2_uncertainty_helpers.h"
#include "gm2_uncertainty_helpers.hpp"

extern "C" {

double gm2calc_thdm_calculate_uncertainty_amu_0loop_amu1L_amu2L(const gm2calc_THDM* model, double amu_1L, double amu_2L)
{
   return gm2calc::calculate_uncertainty_amu_0loop(
      *reinterpret_cast<const gm2calc::THDM*>(model), amu_1L, amu_2L);
}

double gm2calc_thdm_calculate_uncertainty_amu_1loop_amu1L_amu2L(const gm2calc_THDM* model, double amu_1L, double amu_2L)
{
   return gm2calc::calculate_uncertainty_amu_1loop(
      *reinterpret_cast<const gm2calc::THDM*>(model), amu_1L, amu_2L);
}

double gm2calc_thdm_calculate_uncertainty_amu_2loop_amu1L_amu2L(const gm2calc_THDM* model, double amu_1L, double amu_2L)
{
   return gm2calc::calculate_uncertainty_amu_2loop(
      *reinterpret_cast<const gm2calc::THDM*>(model), amu_1L, amu_2L);
}

double gm2calc_thdm_calculate_uncertainty_amu_0loop(const gm2calc_THDM* model)
{
   return gm2calc::calculate_uncertainty_amu_0loop(
      *reinterpret_cast<const gm2calc::THDM*>(model));
}

double gm2calc_thdm_calculate_uncertainty_amu_1loop(const gm2calc_THDM* model)
{
   return gm2calc::calculate_uncertainty_amu_1loop(
      *reinterpret_cast<const gm2calc::THDM*>(model));
}

double gm2calc_thdm_calculate_uncertainty_amu_2loop(const gm2calc_THDM* model)
{
   return gm2calc::calculate_uncertainty_amu_2loop(
      *reinterpret_cast<const gm2calc::THDM*>(model));
}

} /* extern "C" */
