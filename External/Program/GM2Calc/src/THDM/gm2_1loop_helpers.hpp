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

#ifndef GM2_THDM_1LOOP_HELPERS_HPP
#define GM2_THDM_1LOOP_HELPERS_HPP

#include <Eigen/Core>

namespace gm2calc {

namespace thdm {

/// parameters to be passed to the 1-loop contribution functions
struct THDM_1L_parameters {
   double alpha_em{}; ///< electromagnetic coupling
   double mm{};       ///< muon mass for prefactor
   double mw{};       ///< W boson mass
   double mz{};       ///< Z boson mass
   double mhSM{};     ///< SM Higgs boson mass
   double mA{};       ///< CP-odd Higgs boson mass
   double mHp{};      ///< charged Higgs boson mass
   Eigen::Matrix<double,3,1> ml{Eigen::Matrix<double,3,1>::Zero()};  ///< down-type lepton masses
   Eigen::Matrix<double,3,1> mv{Eigen::Matrix<double,3,1>::Zero()};  ///< neutrino masses
   Eigen::Matrix<double,2,1> mh{Eigen::Matrix<double,2,1>::Zero()};  ///< CP-even Higgs bosons mass
   Eigen::Matrix<std::complex<double>,3,3> ylh{Eigen::Matrix<std::complex<double>,3,3>::Zero()}; ///< Y_l^h coefficients with l={e,m,τ}
   Eigen::Matrix<std::complex<double>,3,3> ylH{Eigen::Matrix<std::complex<double>,3,3>::Zero()}; ///< Y_l^H coefficients with l={e,m,τ}
   Eigen::Matrix<std::complex<double>,3,3> ylA{Eigen::Matrix<std::complex<double>,3,3>::Zero()}; ///< Y_l^A coefficients with l={e,m,τ}
   Eigen::Matrix<std::complex<double>,3,3> ylHp{Eigen::Matrix<std::complex<double>,3,3>::Zero()};///< Y_l^{H^\pm} coefficients with l={e,m,τ}
};

// === 1-loop contributions ===

double amu1L(const THDM_1L_parameters&) noexcept;

// === approximations ===

double amu1L_approx(const THDM_1L_parameters&) noexcept;

// === auxiliary functions ===

/// 1-loop THDM contribution to \f$\Delta\alpha\f$
double delta_alpha(double alpha, double mHp, double q) noexcept;

} // namespace thdm

} // namespace gm2calc

#endif
