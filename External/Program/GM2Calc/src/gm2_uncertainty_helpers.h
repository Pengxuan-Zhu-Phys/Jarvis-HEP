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

#ifndef GM2_UNCERTAINTY_HELPERS_H
#define GM2_UNCERTAINTY_HELPERS_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct gm2calc_THDM gm2calc_THDM;
typedef struct MSSMNoFV_onshell MSSMNoFV_onshell;

/** calculates uncertainty for amu(0-loop) */
double gm2calc_thdm_calculate_uncertainty_amu_0loop_amu1L_amu2L(const gm2calc_THDM* /* model */, double /* amu_1L */, double /* amu_2L */);

/** calculates uncertainty for amu(1-loop) */
double gm2calc_thdm_calculate_uncertainty_amu_1loop_amu1L_amu2L(const gm2calc_THDM* /* model */, double /* amu_1L */, double /* amu_2L */);

/** calculates uncertainty for amu(2-loop) */
double gm2calc_thdm_calculate_uncertainty_amu_2loop_amu1L_amu2L(const gm2calc_THDM* /* model */, double /* amu_1L */, double /* amu_2L */);

/** calculates uncertainty for amu(0-loop) w/ tan(beta) resummation */
double gm2calc_mssmnofv_calculate_uncertainty_amu_0loop_amu1L(const MSSMNoFV_onshell* /* model */, double /* amu_1L*/);

/** calculates uncertainty for amu(1-loop) w/ tan(beta) resummation */
double gm2calc_mssmnofv_calculate_uncertainty_amu_1loop_amu2L(const MSSMNoFV_onshell* /* model */, double /* amu_2L */);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
