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

#include "gm2_mf.hpp"
#include "gm2_log.hpp"
#include "gm2_numerics.hpp"

#include <cmath>
#include <boost/math/tools/roots.hpp>

/**
 * \file gm2_mf.cpp
 *
 * Contains functions necessary to calculate the running fermion masse
 * in the MS-bar and DR-bar scheme.
 */

namespace gm2calc {

namespace {

const double pi = 3.14159265358979323846;

/**
 * Calculates the strong coupling constant \f$\alpha_s(Q)\f$ in the
 * Standard Model with 5 active quark flavours in the MS-bar scheme at
 * the scale \f$Q\f$ using Eq (9) from arxiv:hep-ph/0207126 .
 *
 * @param scale renormalization scale \f$Q\f$
 * @param lambda_qcd \f$\Lambda_{\text{QCD}}\f$
 *
 * @return \f$\alpha_s(Q)\f$ MS-bar in the SM w/ 5 active flavours
 */
double calculate_alpha_s_SM5_at(double scale, double lambda_qcd) noexcept
{
   const double t = std::log(sqr(scale/lambda_qcd));
   const double it = 1/t;
   const double logt = std::log(t);

   return 12.*pi/23*it * (
      1 + it*(- 348./529*logt
              + sqr(348./529)*it*(sqr(logt - 0.5) - 78073./242208))
      );
}

/**
 * Calculates \f$\Lambda_{\text{QCD}}\f$ from a given
 * \f$\alpha_s(Q)\f$, where \f$\alpha_s(Q)\f$ is defined in the MS-bar
 * scheme in the SM with 5 active quark flavours.  A root finding
 * algorithm is used to solve Eq. (9) from arxiv:hep-ph/0207126
 * numerically for \f$\Lambda_{\text{QCD}}\f$.
 *
 * @param alpha \f$\alpha_s(Q)\f$ MS-bar, SM w/ 5 active quark flavours
 * @param scale renormalization scale Q
 * @param lambda_qcd_min minimum \f$\Lambda_{\text{QCD}}\f$
 * @param lambda_qcd_max maximum \f$\Lambda_{\text{QCD}}\f$
 * @param precision_goal accuracy goal
 * @param max_iterations maximum number of iterations
 *
 * @return \f$\Lambda_{\text{QCD}}\f$
 */
double calculate_lambda_qcd(double alpha, double scale,
                            double lambda_qcd_min = 0.001,
                            double lambda_qcd_max = 10.,
                            double precision_goal = 1e-10,
                            unsigned max_iterations = 1000) noexcept
{
   boost::uintmax_t it = max_iterations;

   // difference between goal alpha and alpha calculated using the
   // current lambda_qcd
   auto Difference_alpha = [alpha, scale](double lambda) -> double {
      const double alpha_current = calculate_alpha_s_SM5_at(scale, lambda);
      return alpha - alpha_current;
   };

   // stopping criterion, given two brackets a, b
   auto Stop_crit = [precision_goal](double a, double b) -> bool {
      return std::abs(a - b) < precision_goal;
   };

   double lambda_qcd = 0.217; // Nf = 5, PDG

   // find the root
   try {
      const std::pair<double,double> root =
         boost::math::tools::toms748_solve(Difference_alpha, lambda_qcd_min,
                                           lambda_qcd_max, Stop_crit, it);

      lambda_qcd = 0.5 * (root.first + root.second);
   } catch (const std::exception& e) {
      WARNING("Could not determine lambda_QCD: " << e.what()
              << ".  Using lambda_QCD = " << lambda_qcd);
   }

   if (it >= max_iterations) {
      const double precision = std::abs(Difference_alpha(lambda_qcd));
      WARNING("Calculation of Lambda_QCD did not converge"
              " (reached accuracy: " << precision <<
              ", accuracy goal: " << precision_goal <<
              ", max. iterations: " << max_iterations << ")");
   }

   return lambda_qcd;
}

/**
 * Calculates \f$F_b(\mu)\f$, Eq (5) of arxiv:hep-ph/0207126 .
 *
 * @param alpha strong coupling constant at the scale \f$\mu\f$
 *
 * @note Eq (5) of arxiv:hep-ph/0207126 assumes 5 active quark
 * flavours (nf = 5).  For this reason this function can only be used
 * for renormalization group running between Q = mb and Q = mt.
 *
 * @return \f$F_b(\mu)\f$, Eq (5) of arxiv:hep-ph/0207126
 */
double Fb(double alpha) noexcept
{
   const double as = alpha / pi;

   return std::pow(23./6*as, 12./23) * (
      1 + as*(3731./3174 + 1.500706*as)
      );
}

/**
 * MS-bar to DR-bar conversion for mb, Eq (11) of arxiv:hep-ph/0207126 .
 *
 * @param alpha strong coupling constant
 *
 * @return MS-bar to DR-bar conversion factor Eq (11)
 */
double conversion_mb_MSbar_to_DRbar(double alpha) noexcept
{
   const double as = alpha / pi;

   return 1 + as*(-1./3 - 29./72*as);
}

/**
 * Calculates the MS-bar strong coupling alpha_s(Q=mt_pole) from
 * alpha_s(Q=mz) using the 1-loop QCD beta function (nf = 5).
 *
 * @note Taken from SOFTSUSY
 *
 * @param mt_pole top quark pole mass
 * @param alpha_s_mz strong coupling at Q = mz
 * @param mz Z boson pole mass
 *
 * @return alpha_s(mt_pole)
 */
double calculate_alpha_s_SM6_MSbar_at_mt(
   double mt_pole, double alpha_s_mz, double mz) noexcept
{
   return alpha_s_mz /(1 - 23/(6*pi)*alpha_s_mz*std::log(mz/mt_pole));
}

/**
 * Calculates top quark MS-bar mass in the SM mt(MS-bar,Q=mt_pole)
 * from the top quark pole mass, using the 1-loop QCD contribution.
 *
 * @note Taken from SOFTSUSY
 *
 * @param mt_pole top quark pole mass
 * @param as strong coupling at Q = mt_pole
 *
 * @return mt(MS-bar,Q=mt_pole)
 */
double calculate_mt_SM6_MSbar_at(double mt_pole, double as) noexcept
{
   return mt_pole/(1 + 4/(3*pi)*as);
}

} // anonymous namespace

/**
 * Calculates mb(Q) in the DR-bar scheme in the SM w/ 5 active quark
 * flavours using the approach described in arxiv:hep-ph/0207126 .
 *
 * @param mb_mb mb(mb) MS-bar in SM w/ 5 active quark flavours
 * @param alpha_s alpha_s MS-bar in SM w/ 5 quark flavours at scale Q
 * @param scale renormalization scale Q
 *
 * @return mb(Q) DR-bar in the SM w/ 5 quarks
 */
double calculate_mb_SM5_DRbar(
   double mb_mb, double alpha_s, double scale)
{
   // determine Lambda_QCD
   const double lambda_qcd = calculate_lambda_qcd(alpha_s, scale);

   // calculate alpha_s(mb)
   const double alpha_s_mb = calculate_alpha_s_SM5_at(mb_mb, lambda_qcd);

   // run mb to destination scale
   // Here alpha_s must be given at the destination scale `scale'.
   const double mb = mb_mb * Fb(alpha_s) / Fb(alpha_s_mb);

   // DR-bar conversion
   const double mb_DRbar = mb * conversion_mb_MSbar_to_DRbar(alpha_s);

   return mb_DRbar;
}

/**
 * Calculates the running top quark MS-bar mass mt(SM(6),Q) at the
 * scale Q.
 *
 * @param mt_pole top quark pole mass
 * @param alpha_s_at_mz strong coupling at the scale Q = mz
 * @param mz Z boson pole mass
 * @param scale renormalization scale
 *
 * @return mt(SM(6),MS-bar,Q)
 */
double calculate_mt_SM6_MSbar(
   double mt_pole, double alpha_s_mz, double mz, double scale) noexcept
{
   // alpha_s(SM(6), MS-bar, Q = mt_pole)
   const double alpha_s_mt = calculate_alpha_s_SM6_MSbar_at_mt(mt_pole, alpha_s_mz, mz);

   // mt(SM(6), MS-bar, Q = mt_pole)
   const double mt_mt = calculate_mt_SM6_MSbar_at(mt_pole, alpha_s_mt);

   // run from Q = mt_pole to Q = scale
   const double mt_scale = mt_mt * std::pow(scale/mt_pole, -2/pi*alpha_s_mt);

   // Note: QCD beta function of the quark mass parameter:
   // dm/d(log(Q)) = -2/pi*alpha_s*m

   return mt_scale;
}

/**
 * Calculates the running bottom quark MS-bar mass mb(SM(6),Q) in the
 * SM(6) at the scale Q.
 *
 * @param mb_mb bottom quark MS-bar mass mb(mb) in the SM(5)
 * @param mt_pole top quark pole mass
 * @param alpha_s_mz strong coupling at the scale mz
 * @param mz Z boson pole mass
 * @param scale renormalization scale
 *
 * @return mb(MS-bar,SM(6),Q)
 */
double calculate_mb_SM6_MSbar(
   double mb_mb, double mt_pole, double alpha_s_mz, double mz, double scale) noexcept
{
   // determine Lambda_QCD
   const double lambda_qcd = calculate_lambda_qcd(alpha_s_mz, mz);

   // calculate alpha_s(mb)
   const double alpha_s_mb = calculate_alpha_s_SM5_at(mb_mb, lambda_qcd);

   // calculate alpha_s(mt)
   const double alpha_s_mt = calculate_alpha_s_SM5_at(mt_pole, lambda_qcd);

   // run mb(mb) to Q = mt_pole
   const double mb_mt = mb_mb * Fb(alpha_s_mt) / Fb(alpha_s_mb);

   // run mb(mt) to Q = scale
   const double mb_scale = mb_mt * std::pow(scale/mt_pole, -2/pi*alpha_s_mt);

   return mb_scale;
}

/**
 * Calculates the running tau lepton MS-bar mass mtau(SM(6),Q) in the
 * SM(6) at the scale Q.
 *
 * @param mtau_pole tau lepton pole mass
 * @param alpha_em_mz electromagnetic coupling at the scale Q = MZ
 * @param scale renormalization scale
 *
 * @return mtau(MS-bar,SM(6),Q)
 */
double calculate_mtau_SM6_MSbar(
   double mtau_pole, double alpha_em_mz, double scale) noexcept
{
   // calculate mtau(mtau)
   const double mtau_mtau = mtau_pole; // neglecting loop corrections

   // Note: QED beta function of the lepton mass parameter:
   // dm/d(log(Q)) = -3/(2*pi)*alpha_em*m

   // run mtau(mtau) to Q = scale
   const double mtau_scale = mtau_mtau * std::pow(scale/mtau_mtau, -3/(2*pi)*alpha_em_mz);

   return mtau_scale;
}

} // namespace gm2calc
