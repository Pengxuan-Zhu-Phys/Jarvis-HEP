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

#ifndef GM2_MF_HPP
#define GM2_MF_HPP

namespace gm2calc {

/// calculates mb(Q) DR-bar
double calculate_mb_SM5_DRbar(double mb_mb, double alpha_s, double scale);

/// calculates mt(Q) MS-bar in the SM(6)
double calculate_mt_SM6_MSbar(double mt_pole, double alpha_s_mz, double mz, double scale) noexcept;

/// calculates mb(Q) MS-bar in the SM(6)
double calculate_mb_SM6_MSbar(double mb_mb, double mt_pole, double alpha_s_mz, double mz, double scale) noexcept;

/// calculates mtau(Q) MS-bar in the SM(6)
double calculate_mtau_SM6_MSbar(double mtau_pole, double alpha_em_mz, double scale) noexcept;

} // namespace gm2calc

#endif
