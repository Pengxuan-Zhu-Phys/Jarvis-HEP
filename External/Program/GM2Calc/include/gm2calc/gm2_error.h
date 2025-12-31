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

#ifndef GM2_ERROR_H
#define GM2_ERROR_H

/**
 * @file gm2_error.h
 * @brief contains declarations of GM2Calc error codes for the C interface
 */

#ifdef __cplusplus
extern "C" {
#endif

/** error codes */
typedef enum {
   gm2calc_NoError = 0,
   gm2calc_InvalidInput,
   gm2calc_PhysicalProblem,
   gm2calc_UnknownError
} gm2calc_error;

/** translate error codes into a string */
const char* gm2calc_error_str(gm2calc_error);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
