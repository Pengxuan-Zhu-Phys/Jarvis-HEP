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

#include "gm2calc/THDM_parameters.hpp"

#include <iostream>

namespace gm2calc {

void THDM_parameters::print(std::ostream& ostr) const
{
   ostr << "----------------------------------------\n"
           "Parameters:\n"
           "----------------------------------------\n";
   ostr << "g1 = " << g1 << '\n';
   ostr << "g2 = " << g2 << '\n';
   ostr << "g3 = " << g3 << '\n';
   ostr << "lambda1 = " << lambda1 << '\n';
   ostr << "lambda2 = " << lambda2 << '\n';
   ostr << "lambda3 = " << lambda3 << '\n';
   ostr << "lambda4 = " << lambda4 << '\n';
   ostr << "lambda5 = " << lambda5 << '\n';
   ostr << "lambda6 = " << lambda6 << '\n';
   ostr << "lambda7 = " << lambda7 << '\n';
   ostr << "Gamma_u =\n" << Gamma_u << '\n';
   ostr << "Gamma_d =\n" << Gamma_d << '\n';
   ostr << "Gamma_l =\n" << Gamma_l << '\n';
   ostr << "Pi_u =\n" << Pi_u << '\n';
   ostr << "Pi_d =\n" << Pi_d << '\n';
   ostr << "Pi_l =\n" << Pi_l << '\n';
   ostr << "m122 = " << m122 << " GeV^2\n";
   ostr << "m112 = " << m112 << " GeV^2\n";
   ostr << "m222 = " << m222 << " GeV^2\n";
   ostr << "v1 = " << v1 << " GeV\n";
   ostr << "v2 = " << v2 << " GeV\n";
}

std::ostream& operator<<(std::ostream& ostr, const THDM_parameters& pars)
{
   pars.print(ostr);
   return ostr;
}

} // namespace gm2calc
