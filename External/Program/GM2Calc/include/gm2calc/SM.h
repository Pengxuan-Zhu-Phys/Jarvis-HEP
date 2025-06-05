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

#ifndef GM2_SM_H
#define GM2_SM_H

/**
 * @file SM.h
 * @brief contains declarations of C interface functions for the SM
 *
 * This file contains the declarations for the C interface functions
 * used to set and retrieve the SM parameters and masses.
 */

#ifdef __cplusplus
extern "C" {
#endif

struct gm2calc_SM {
   double alpha_em_0;
   double alpha_em_mz;
   double alpha_s_mz;
   double mh;
   double mw;
   double mz;
   double mu[3];
   double md[3];
   double mv[3];
   double ml[3];
   double ckm_real[3][3];
   double ckm_imag[3][3];
};

typedef struct gm2calc_SM gm2calc_SM;

void gm2calc_sm_set_to_default(gm2calc_SM*);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
