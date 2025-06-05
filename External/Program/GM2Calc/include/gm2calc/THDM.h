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

#ifndef GM2_THDM_H
#define GM2_THDM_H

/**
 * @file THDM.h
 * @brief contains declarations of C interface functions for the model
 *
 * This file contains the declarations for the C interface functions
 * used to set and retrieve the model parameters and masses.
 */

#include "gm2calc/gm2_error.h"

#ifdef __cplusplus
extern "C" {
#endif

/** THDM configuration options */
struct gm2calc_THDM_config {
   int force_output;
   int running_couplings;
};
typedef struct gm2calc_THDM_config gm2calc_THDM_config;

/** create a new config */
void gm2calc_thdm_config_set_to_default(gm2calc_THDM_config*);

/** general THDM handle */
struct gm2calc_THDM;
typedef struct gm2calc_THDM gm2calc_THDM;

/** Yukawa schemes */
typedef enum {
   gm2calc_THDM_type_1 = 1,
   gm2calc_THDM_type_2,
   gm2calc_THDM_type_X,
   gm2calc_THDM_type_Y,
   gm2calc_THDM_aligned,
   gm2calc_THDM_general,
} gm2calc_THDM_yukawa_type;

/** convert int to gm2calc_THDM_yukawa_type */
gm2calc_THDM_yukawa_type int_to_c_yukawa_type(int);

/** general THDM general basis input */
struct gm2calc_THDM_gauge_basis {
   gm2calc_THDM_yukawa_type yukawa_type;
   double lambda[7];
   double tan_beta;
   double m122;
   double zeta_u;
   double zeta_d;
   double zeta_l;
   double Delta_u[3][3];
   double Delta_d[3][3];
   double Delta_l[3][3];
   double Pi_u[3][3];
   double Pi_d[3][3];
   double Pi_l[3][3];
};
typedef struct gm2calc_THDM_gauge_basis gm2calc_THDM_gauge_basis;

/** general THDM physical basis input */
struct gm2calc_THDM_mass_basis {
   gm2calc_THDM_yukawa_type yukawa_type;
   double mh;
   double mH;
   double mA;
   double mHp;
   double sin_beta_minus_alpha;
   double lambda_6;
   double lambda_7;
   double tan_beta;
   double m122;
   double zeta_u;
   double zeta_d;
   double zeta_l;
   double Delta_u[3][3];
   double Delta_d[3][3];
   double Delta_l[3][3];
   double Pi_u[3][3];
   double Pi_d[3][3];
   double Pi_l[3][3];
};
typedef struct gm2calc_THDM_mass_basis gm2calc_THDM_mass_basis;

struct gm2calc_SM;
typedef struct gm2calc_SM gm2calc_SM;

/** allocate new general THDM model with general basis input */
gm2calc_error gm2calc_thdm_new_with_gauge_basis(gm2calc_THDM**, const gm2calc_THDM_gauge_basis*, const gm2calc_SM*, const gm2calc_THDM_config*);

/** allocate new general THDM model with physical basis input */
gm2calc_error gm2calc_thdm_new_with_mass_basis(gm2calc_THDM**, const gm2calc_THDM_mass_basis*, const gm2calc_SM*, const gm2calc_THDM_config*);

/** delete general THDM model */
void gm2calc_thdm_free(gm2calc_THDM*);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
