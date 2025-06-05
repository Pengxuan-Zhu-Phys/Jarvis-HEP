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

#include "gm2calc/SM.hpp"
#include "gm2calc/gm2_error.hpp"
#include "gm2_numerics.hpp"
#include "gm2_constants.hpp"
#include <cmath>
#include <complex>
#include <iostream>

namespace gm2calc {

namespace {

const double pi = 3.1415926535897932;

Eigen::Matrix<std::complex<double>,3,3> get_ckm_from_angles(
   double theta_12, double theta_13, double theta_23, double delta) noexcept
{
   const std::complex<double> eID(std::polar(1.0, delta));
   const double s12 = std::sin(theta_12);
   const double s13 = std::sin(theta_13);
   const double s23 = std::sin(theta_23);
   const double c12 = std::cos(theta_12);
   const double c13 = std::cos(theta_13);
   const double c23 = std::cos(theta_23);

   Eigen::Matrix<std::complex<double>,3,3> ckm;
   ckm(0, 0) = c12 * c13;
   ckm(0, 1) = s12 * c13;
   ckm(0, 2) = s13 / eID;
   ckm(1, 0) = -s12 * c23 - c12 * s23 * s13 * eID;
   ckm(1, 1) = c12 * c23 - s12 * s23 * s13 * eID;
   ckm(1, 2) = s23 * c13;
   ckm(2, 0) = s12 * s23 - c12 * c23 * s13 * eID;
   ckm(2, 1) = -c12 * s23 - s12 * c23 * s13 * eID;
   ckm(2, 2) = c23 * c13;

   return ckm;
}

Eigen::Matrix<std::complex<double>,3,3> get_ckm_from_wolfenstein(
   double lambdaW, double aCkm, double rhobar, double etabar)
{
   if (std::abs(lambdaW) > 1.0) { throw EInvalidInput("Error: Wolfenstein lambda out of range!"); }
   if (std::abs(aCkm)    > 1.0) { throw EInvalidInput("Error: Wolfenstein A parameter out of range!"); }
   if (std::abs(rhobar)  > 1.0) { throw EInvalidInput("Error: Wolfenstein rho-bar parameter out of range!"); }
   if (std::abs(etabar)  > 1.0) { throw EInvalidInput("Error: Wolfenstein eta-bar parameter out of range!"); }

   const double theta_12 = std::asin(lambdaW);
   const double theta_23 = std::asin(aCkm * sqr(lambdaW));
   const double lambdaW3 = cube(lambdaW);
   const double lambdaW4 = pow4(lambdaW);

   const std::complex<double> rpe(rhobar, etabar);
   const std::complex<double> V13conj = aCkm * lambdaW3 * rpe
      * std::sqrt(1.0 - sqr(aCkm) * lambdaW4) /
      std::sqrt(1.0 - sqr(lambdaW)) / (1.0 - sqr(aCkm) * lambdaW4 * rpe);

   double theta_13 = 0.0;
   double delta = 0.0;

   if (std::isfinite(std::real(V13conj)) && std::isfinite(std::imag(V13conj))) {
      theta_13 = std::asin(std::abs(V13conj));
      delta = std::arg(V13conj);
   }

   return get_ckm_from_angles(theta_12, theta_13, theta_23, delta);
}

} // anonymous namespace

SM::SM()
   : alpha_em_0(ALPHA_EM_THOMPSON)
   , alpha_em_mz(ALPHA_EM_MZ)
   , alpha_s_mz(ALPHA_S_MZ)
   , mh(MH)
   , mw(MW)
   , mz(MZ)
   , mu((Eigen::Matrix<double,3,1>() << MU, MC, MT).finished())
   , md((Eigen::Matrix<double,3,1>() << MD, MS, MBMB).finished())
   , ml((Eigen::Matrix<double,3,1>() << ME, MM, ML).finished())
   , ckm(get_ckm_from_angles(CKM_THETA12, CKM_THETA13, CKM_THETA23, CKM_DELTA))
{
}

double SM::get_e_0() const
{
   return std::sqrt(alpha_em_0*4*pi);
}

double SM::get_e_mz() const
{
   return std::sqrt(alpha_em_mz*4*pi);
}

double SM::get_gY() const
{
   return get_e_mz()/get_cw();
}

double SM::get_g2() const
{
   return get_e_mz()/get_sw();
}

double SM::get_g3() const
{
   return std::sqrt(alpha_s_mz*4*pi);
}

double SM::get_cw() const
{
   return std::abs(mw/mz);
}

double SM::get_sw() const
{
   return std::sqrt(1 - sqr(get_cw()));
}

double SM::get_v() const
{
   return 2*mw/get_g2();
}

void SM::set_ckm_from_wolfenstein(
   double lambdaW, double aCkm, double rhobar, double etabar)
{
   ckm = get_ckm_from_wolfenstein(lambdaW, aCkm, rhobar, etabar);
}

void SM::set_ckm_from_angles(
   double theta_12, double theta_13, double theta_23, double delta)
{
   ckm = get_ckm_from_angles(theta_12, theta_13, theta_23, delta);
}

std::ostream& operator<<(std::ostream& ostr, const SM& sm)
{
   ostr << "========================================\n"
           " SM\n"
           "========================================\n"
        << "alpha_em(MZ) = " << sm.get_alpha_em_mz() << '\n'
        << "alpha_em(0) = " << sm.get_alpha_em_0() << '\n'
        << "alpha_s(MZ) = " << sm.get_alpha_s_mz() << '\n'
        << "mw = " << sm.get_mw() << " GeV\n"
        << "mz = " << sm.get_mz() << " GeV\n"
        << "mh = " << sm.get_mh() << " GeV\n"
        << "mu = {" << sm.get_mu(0) << ", " << sm.get_mu(1) << ", " << sm.get_mu(2) << "} GeV\n"
        << "md = {" << sm.get_md(0) << ", " << sm.get_md(1) << ", " << sm.get_md(2) << "} GeV\n"
        << "mv = {" << sm.get_mv(0) << ", " << sm.get_mv(1) << ", " << sm.get_mv(2) << "} GeV\n"
        << "ml = {" << sm.get_ml(0) << ", " << sm.get_ml(1) << ", " << sm.get_ml(2) << "} GeV\n"
        << "CKM = {{" << sm.get_ckm(0,0) << ", " << sm.get_ckm(0,1) << ", " << sm.get_ckm(0,2) << "},\n"
           "       {" << sm.get_ckm(1,0) << ", " << sm.get_ckm(1,1) << ", " << sm.get_ckm(1,2) << "},\n"
           "       {" << sm.get_ckm(2,0) << ", " << sm.get_ckm(2,1) << ", " << sm.get_ckm(2,2) << "}}\n"
        ;
   return ostr;
}

} // namespace gm2calc
