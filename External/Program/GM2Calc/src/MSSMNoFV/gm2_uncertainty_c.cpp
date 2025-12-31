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

double gm2calc_mssmnofv_calculate_uncertainty_amu_0loop_amu1L(const MSSMNoFV_onshell* model, double amu_1L)
{
   return gm2calc::calculate_uncertainty_amu_0loop(
      *reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model), amu_1L);
}

double gm2calc_mssmnofv_calculate_uncertainty_amu_0loop(const MSSMNoFV_onshell* model)
{
   return gm2calc::calculate_uncertainty_amu_0loop(
      *reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model));
}

double gm2calc_mssmnofv_calculate_uncertainty_amu_1loop_amu2L(const MSSMNoFV_onshell* model, double amu_2L)
{
   return gm2calc::calculate_uncertainty_amu_1loop(
      *reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model), amu_2L);
}

double gm2calc_mssmnofv_calculate_uncertainty_amu_1loop(const MSSMNoFV_onshell* model)
{
   return gm2calc::calculate_uncertainty_amu_1loop(
      *reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model));
}

double gm2calc_mssmnofv_calculate_uncertainty_amu_2loop(const MSSMNoFV_onshell* model)
{
   return gm2calc::calculate_uncertainty_amu_2loop(
      *reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model));
}

} /* extern "C" */
