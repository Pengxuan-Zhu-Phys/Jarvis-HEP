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

#ifndef GM2_DILOG_HPP
#define GM2_DILOG_HPP

#include <complex>

namespace gm2calc {

/// real dilogarithm
double dilog(double) noexcept;

/// complex dilogarithm
std::complex<double> dilog(const std::complex<double>&) noexcept;

/// Clausen function Cl_2(x)
double clausen_2(double x) noexcept;

} // namespace gm2calc

#endif
