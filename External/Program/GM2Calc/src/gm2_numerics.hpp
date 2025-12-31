// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#ifndef GM2_NUMERICS_HPP
#define GM2_NUMERICS_HPP

#include <cmath>
#include <complex>
#include <limits>

namespace gm2calc {

/// returns number squared
template <typename T> T sqr(T x) noexcept { return x*x; }

/// returns number to the third power
template <typename T> T cube(T x) noexcept { return x*x*x; }

/// returns number to the third power
template <typename T> T pow3(T x) noexcept { return cube(x); }

/// returns number to the 4th power
template <typename T> T pow4(T x) noexcept { return sqr(sqr(x)); }

/// returns square root of absolute of number
double abs_sqrt(double) noexcept;

/// returns sign of real number
int sign(double) noexcept;

/// returns square root of absolute of number, times sign
double signed_abs_sqrt(double) noexcept;

/// returns square of number, times sign
double signed_sqr(double) noexcept;

template <typename T>
bool is_zero(T a, T eps) noexcept
{
   return std::fabs(a) < eps;
}

template <typename T>
bool is_equal(T a, T b, T eps) noexcept
{
   return is_zero(a - b, eps);
}

template <typename T>
bool is_equal_rel(T a, T b, T eps) noexcept
{
   const T max = std::max(std::abs(a), std::abs(b));
   return is_zero(a - b, eps*(1.0 + max));
}

} // namespace gm2calc

#endif
