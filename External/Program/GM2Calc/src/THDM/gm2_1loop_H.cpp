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

#include "THDM/gm2_1loop_helpers.hpp"
#include "gm2_ffunctions.hpp"
#include <cmath>
#include <complex>

/**
 * \file gm2_1loop_H.cpp
 *
 * Contains functions necessary to calculate the THDM
 * contributions for g-2 at the 1-loop level.
 */

namespace gm2calc {

namespace thdm {

namespace {

const double pi = 3.1415926535897932;
const double pi2 = 9.8696044010893586; // Pi^2

double sqr(double x) noexcept { return x*x; }

/// Eq.(28), arxiv:1607.06292
double Fh(double x) noexcept
{
   return F1C(x)/12 + F2C(x)/3;
}

/// Eq.(29), arxiv:1607.06292
double FA(double x) noexcept
{
   return F1C(x)/12 - F2C(x)/3;
}

/// Eq.(30), arxiv:1607.06292
double FHp(double x) noexcept
{
   return -F1N(x)/12;
}

double AS(int gen, const Eigen::Matrix<double,3,1>& ml, double mS2, const Eigen::Matrix<std::complex<double>,3,3>& y) noexcept
{
   const auto x = sqr(ml(gen))/mS2;
   const auto y2 = std::conj(y(gen, 1))*std::conj(y(1, gen));

   return
      + (std::norm(y(gen, 1)) + std::norm(y(1, gen)))*F1C(x)/24
      + std::real(y2)*ml(gen)/ml(1)*F2C(x)/3;
}

double AA(int gen, const Eigen::Matrix<double,3,1>& ml, double mS2, const Eigen::Matrix<std::complex<double>,3,3>& y) noexcept
{
   const auto x = sqr(ml(gen))/mS2;
   const auto y2 = std::conj(y(gen, 1))*std::conj(y(1, gen));

   return
      + (std::norm(y(gen, 1)) + std::norm(y(1, gen)))*F1C(x)/24
      - std::real(y2)*ml(gen)/ml(1)*F2C(x)/3;
}

double AHp(int gen, const Eigen::Matrix<double,3,1>& mv, double mS2, const Eigen::Matrix<std::complex<double>,3,3>& y) noexcept
{
   return -std::norm(y(gen, 1))/48*(
      F1N(sqr(mv(1))/mS2) + F1N(sqr(mv(gen))/mS2));
}


} // anonymous namespace

/**
 * Approximation for 1-loop contribution
 * Eq (27) from arxiv:1607.06292
 *
 * @note The factor 1/2 in front of the charged Higgs contribution
 * stems from the fact that \f$Y_l^{H^\pm} = \sqrt{2}\; Y_l^A\f$.  In
 * Eq.(27) the prefactor of the charged Higgs contribution is
 * \f$(Y_l^A)^2\f$, not \f$(Y_l^{H^\pm})^2\f$.
 */
double amu1L_approx(const THDM_1L_parameters& pars) noexcept
{
   const auto mm2 = sqr(pars.mm);
   const auto mw2 = sqr(pars.mw);
   const auto mz2 = sqr(pars.mz);
   const auto mhSM2 = sqr(pars.mhSM);
   const auto mh2 = sqr(pars.mh(0));
   const auto mH2 = sqr(pars.mh(1));
   const auto mA2 = sqr(pars.mA);
   const auto mHp2 = sqr(pars.mHp);
   const auto sw2 = 1 - mw2/mz2;
   const auto e2 = 4*pi*pars.alpha_em;
   const auto g22 = e2/sw2;
   const auto v2 = 4*mw2/g22;

   // Eq.(27), arxiv:1607.06292
   const auto res =
      + std::norm(pars.ylh(1,1))/mh2*Fh(mm2/mh2)
      + std::norm(pars.ylH(1,1))/mH2*Fh(mm2/mH2)
      + std::norm(pars.ylA(1,1))/mA2*FA(mm2/mA2)
      + 0.5 * std::norm(pars.ylHp(1,1))/mHp2*FHp(mm2/mHp2)
      // subtract SM contribution
      - mm2/(v2*mhSM2)*Fh(mm2/mhSM2);

   return mm2*res/(8*pi2);
}

/**
 * Full (CP-conserving) 1-loop contribution
 */
double amu1L(const THDM_1L_parameters& pars) noexcept
{
   const auto mm2 = sqr(pars.mm);
   const auto mw2 = sqr(pars.mw);
   const auto mz2 = sqr(pars.mz);
   const auto mhSM2 = sqr(pars.mhSM);
   const auto mh2 = sqr(pars.mh(0));
   const auto mH2 = sqr(pars.mh(1));
   const auto mA2 = sqr(pars.mA);
   const auto mHp2 = sqr(pars.mHp);
   const auto sw2 = 1 - mw2/mz2;
   const auto e2 = 4*pi*pars.alpha_em;
   const auto g2 = std::sqrt(e2/sw2);
   const auto v = 2*pars.mw/g2;

   const Eigen::Matrix<double, 3, 3> ylhSM{
      (Eigen::Matrix<double, 3, 3>()
       << 0.0, 0.0, 0.0,
          0.0, pars.mm/v, 0.0,
          0.0, 0.0, 0.0).finished()};

   double res = 0.0;

   for (int g = 0; g < 3; ++g) {
      res += AS(g, pars.ml, mh2, pars.ylh)/mh2;
      res += AS(g, pars.ml, mH2, pars.ylH)/mH2;
      res += AA(g, pars.ml, mA2, pars.ylA)/mA2;
      res += AHp(g, pars.mv, mHp2, pars.ylHp)/mHp2;
   }

   // subtract SM contribution
   res -= AS(1, pars.ml, mhSM2, ylhSM)/mhSM2;

   return mm2*res/(8*pi2);
}

/**
 * Calculates the 1-loop THDM contribution to \f$\Delta\alpha\f$.
 * \f$\alpha^{\text{THDM}} = \alpha^{\text{SM}}/(1 - \Delta\alpha)\f$
 *
 * @param alpha electromagnetic coupling
 * @param mHp charged Higgs mass
 * @param q renormalization scale
 *
 * @return 1-loop THDM contribution to \f$\Delta\alpha\f$
 */
double delta_alpha(double alpha, double mHp, double q) noexcept
{
   return -alpha/(6*pi)*std::log(std::abs(mHp/q));
}

} // namespace thdm

} // namespace gm2calc
