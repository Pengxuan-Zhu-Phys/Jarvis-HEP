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

#include "gm2calc/gm2_error.h"

extern "C"
{

const char* gm2calc_error_str(gm2calc_error error)
{
   const char* error_str = "Unknown error";

   switch (error) {
   case gm2calc_NoError:
      error_str = "no error";
      break;
   case gm2calc_InvalidInput:
      error_str = "Input parameter set to invalid value";
      break;
   case gm2calc_PhysicalProblem:
      error_str = "Physical problem has occurred during calculation";
      break;
   default:
      break;
   }

   return error_str;
}

} /* extern "C" */
