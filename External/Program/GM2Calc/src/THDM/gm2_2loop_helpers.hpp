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

#ifndef GM2_THDM_2LOOP_HELPERS_HPP
#define GM2_THDM_2LOOP_HELPERS_HPP

#include <Eigen/Core>

namespace gm2calc {

namespace thdm {

/// parameters to be passed to the bosonic contribution functions
struct THDM_B_parameters {
   double alpha_em{};///< electromagnetic coupling
   double mm{};      ///< muon mass for prefactor
   double mw{};      ///< W boson mass
   double mz{};      ///< Z boson mass
   double mhSM{};    ///< SM Higgs boson mass
   double mA{};      ///< CP-odd Higgs boson mass
   double mHp{};     ///< charged Higgs boson mass
   Eigen::Matrix<double,2,1> mh{Eigen::Matrix<double,2,1>::Zero()}; ///< CP-even Higgs bosons mass
   double tb{};      ///< tan(beta)
   double zetal{};   ///< zeta_l
   double cos_beta_minus_alpha{}; ///< cos(beta - alpha_h)
   double lambda5{}; ///< Lambda_5
   double lambda67{}; ///< difference (Lambda_567 - Lambda_5)
};

/// parameters to be passed to the fermionic contribution functions
struct THDM_F_parameters {
   double alpha_em{}; ///< electromagnetic coupling
   double mm{};       ///< muon mass for prefactor
   double mw{};       ///< W boson mass
   double mz{};       ///< Z boson mass
   double mhSM{};     ///< SM Higgs boson mass
   double mA{};       ///< CP-odd Higgs boson mass
   double mHp{};      ///< charged Higgs boson mass
   Eigen::Matrix<double,2,1> mh{Eigen::Matrix<double,2,1>::Zero()};  ///< CP-even Higgs bosons mass
   Eigen::Matrix<double,3,1> ml{Eigen::Matrix<double,3,1>::Zero()};  ///< down-type lepton masses
   Eigen::Matrix<double,3,1> mu{Eigen::Matrix<double,3,1>::Zero()};  ///< up-type quark masses
   Eigen::Matrix<double,3,1> md{Eigen::Matrix<double,3,1>::Zero()};  ///< down-type quark masses
   Eigen::Matrix<std::complex<double>,3,3> yuh{Eigen::Matrix<std::complex<double>,3,3>::Zero()}; ///< y_f^S coefficients with f={u,c,t} and S=h
   Eigen::Matrix<std::complex<double>,3,3> yuH{Eigen::Matrix<std::complex<double>,3,3>::Zero()}; ///< y_f^S coefficients with f={u,c,t} and S=H
   Eigen::Matrix<std::complex<double>,3,3> yuA{Eigen::Matrix<std::complex<double>,3,3>::Zero()}; ///< y_f^S coefficients with f={u,c,t} and S=A
   Eigen::Matrix<std::complex<double>,3,3> yuHp{Eigen::Matrix<std::complex<double>,3,3>::Zero()};///< y_f^S coefficients with f={u,c,t} and S=H^+
   Eigen::Matrix<std::complex<double>,3,3> ydh{Eigen::Matrix<std::complex<double>,3,3>::Zero()}; ///< y_f^S coefficients with f={d,s,b} and S=h
   Eigen::Matrix<std::complex<double>,3,3> ydH{Eigen::Matrix<std::complex<double>,3,3>::Zero()}; ///< y_f^S coefficients with f={d,s,b} and S=H
   Eigen::Matrix<std::complex<double>,3,3> ydA{Eigen::Matrix<std::complex<double>,3,3>::Zero()}; ///< y_f^S coefficients with f={d,s,b} and S=A
   Eigen::Matrix<std::complex<double>,3,3> ydHp{Eigen::Matrix<std::complex<double>,3,3>::Zero()};///< y_f^S coefficients with f={d,s,b} and S=H^+
   Eigen::Matrix<std::complex<double>,3,3> ylh{Eigen::Matrix<std::complex<double>,3,3>::Zero()}; ///< y_f^S coefficients with f={e,m,τ} and S=h
   Eigen::Matrix<std::complex<double>,3,3> ylH{Eigen::Matrix<std::complex<double>,3,3>::Zero()}; ///< y_f^S coefficients with f={e,m,τ} and S=H
   Eigen::Matrix<std::complex<double>,3,3> ylA{Eigen::Matrix<std::complex<double>,3,3>::Zero()}; ///< y_f^S coefficients with f={e,m,τ} and S=A
   Eigen::Matrix<std::complex<double>,3,3> ylHp{Eigen::Matrix<std::complex<double>,3,3>::Zero()};///< y_f^S coefficients with f={e,m,τ} and S=H^+
   Eigen::Matrix<std::complex<double>,3,3> vckm{Eigen::Matrix<std::complex<double>,3,3>::Identity()};///< CKM matrix
};

// === 2-loop bosonic contributions ===

double amu2L_B(const THDM_B_parameters&) noexcept;

// routines for sub-expressions

double amu2L_B_EWadd(const THDM_B_parameters&) noexcept;
double amu2L_B_nonYuk(const THDM_B_parameters&) noexcept;
double amu2L_B_Yuk(const THDM_B_parameters&) noexcept;

// === 2-loop fermionic contributions ===

double amu2L_F(const THDM_F_parameters&) noexcept;

// routines for sub-expressions

double amu2L_F_charged(const THDM_F_parameters&) noexcept;
double amu2L_F_neutral(const THDM_F_parameters&) noexcept;

/// Eq (53), arxiv:1607.06292, f = u, S = h or H
double fuS(double ms2, double mu2, double mw2, double mz2) noexcept;
/// Eq (53), arxiv:1607.06292, f = d, S = h or H
double fdS(double ms2, double md2, double mw2, double mz2) noexcept;
/// Eq (53), arxiv:1607.06292, f = l, S = h or H
double flS(double ms2, double ml2, double mw2, double mz2) noexcept;
/// Eq (53), arxiv:1607.06292, f = u, S = A
double fuA(double ms2, double mu2, double mw2, double mz2) noexcept;
/// Eq (53), arxiv:1607.06292, f = d, S = A
double fdA(double ms2, double md2, double mw2, double mz2) noexcept;
/// Eq (53), arxiv:1607.06292, f = l, S = A
double flA(double ms2, double ml2, double mw2, double mz2) noexcept;
/// Eq (59), arxiv:1607.06292, S = H^\pm, f = u
double fuHp(double ms2, double md2, double mu2, double mw2, double mz2) noexcept;
/// Eq (59), arxiv:1607.06292, S = H^\pm, f = d
double fdHp(double ms2, double md2, double mu2, double mw2, double mz2) noexcept;
/// Eq (59), arxiv:1607.06292, S = H^\pm, f = l
double flHp(double ms2, double ml2, double mw2, double mz2) noexcept;

} // namespace thdm

} // namespace gm2calc

#endif
