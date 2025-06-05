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

#ifndef GM2_SM_HPP
#define GM2_SM_HPP

#include <Eigen/Core>
#include <complex>
#include <iosfwd>

namespace gm2calc {

class SM {
public:
   SM();

   void set_alpha_em_0(double a) noexcept { alpha_em_0 = a; }
   void set_alpha_em_mz(double a) noexcept { alpha_em_mz = a; }
   void set_alpha_s_mz(double a) noexcept { alpha_s_mz = a; }
   void set_mh(double m) noexcept { mh = m; }
   void set_mw(double m) noexcept { mw = m; }
   void set_mz(double m) noexcept { mz = m; }
   void set_mu(const Eigen::Matrix<double,3,1>& m) noexcept { mu = m; }
   void set_md(const Eigen::Matrix<double,3,1>& m) noexcept { md = m; }
   void set_mv(const Eigen::Matrix<double,3,1>& m) noexcept { mv = m; }
   void set_ml(const Eigen::Matrix<double,3,1>& m) noexcept { ml = m; }
   void set_mu(int i, double m) noexcept { mu(i) = m; }
   void set_md(int i, double m) noexcept { md(i) = m; }
   void set_mv(int i, double m) noexcept { mv(i) = m; }
   void set_ml(int i, double m) noexcept { ml(i) = m; }
   void set_ckm(const Eigen::Matrix<std::complex<double>,3,3>& m) { ckm = m; }
   void set_ckm(int i, int j, const std::complex<double>& m) { ckm(i, j) = m; }
   void set_ckm_from_wolfenstein(double lambdaW, double aCkm, double rhobar, double etabar);
   void set_ckm_from_angles(double theta_12, double theta_13, double theta_23, double delta);

   double get_alpha_em_0() const { return alpha_em_0; }
   double get_alpha_em_mz() const { return alpha_em_mz; }
   double get_alpha_s_mz() const { return alpha_s_mz; }
   double get_mh() const { return mh; }
   double get_mw() const { return mw; }
   double get_mz() const { return mz; }
   const Eigen::Matrix<double,3,1>& get_mu() const { return mu; }
   const Eigen::Matrix<double,3,1>& get_md() const { return md; }
   const Eigen::Matrix<double,3,1>& get_mv() const { return mv; }
   const Eigen::Matrix<double,3,1>& get_ml() const { return ml; }
   double get_mu(int i) const { return mu(i); }
   double get_md(int i) const { return md(i); }
   double get_mv(int i) const { return mv(i); }
   double get_ml(int i) const { return ml(i); }
   double get_e_0() const;
   double get_e_mz() const;
   double get_gY() const;
   double get_g2() const;
   double get_g3() const;
   double get_cw() const;
   double get_sw() const;
   double get_v() const;
   const Eigen::Matrix<std::complex<double>,3,3>& get_ckm() const { return ckm; }
   std::complex<double> get_ckm(int i, int j) const { return ckm(i, j); }

private:
   double alpha_em_0{0.0};  ///< electromagnetic coupling in Thompson limit
   double alpha_em_mz{0.0}; ///< electromagnetic coupling at Q = MZ
   double alpha_s_mz{0.0};  ///< strong coupling at Q = MZ
   double mh{0.0};          ///< Higgs boson pole mass
   double mw{0.0};          ///< W boson pole mass
   double mz{0.0};          ///< Z boson pole mass
   Eigen::Matrix<double,3,1> mu{Eigen::Matrix<double,3,1>::Zero()}; ///< up-type quark masses
   Eigen::Matrix<double,3,1> md{Eigen::Matrix<double,3,1>::Zero()}; ///< down-type quark masses
   Eigen::Matrix<double,3,1> mv{Eigen::Matrix<double,3,1>::Zero()}; ///< neutrino masses
   Eigen::Matrix<double,3,1> ml{Eigen::Matrix<double,3,1>::Zero()}; ///< down-type lepton pole masses
   Eigen::Matrix<std::complex<double>,3,3> ckm{Eigen::Matrix<std::complex<double>,3,3>::Identity()}; ///< CKM matrix
};

/// streaming operator
std::ostream& operator<<(std::ostream&, const SM&);

} // namespace gm2calc

#endif
