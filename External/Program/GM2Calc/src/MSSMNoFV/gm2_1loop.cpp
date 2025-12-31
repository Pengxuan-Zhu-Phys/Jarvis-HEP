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

#include "gm2calc/gm2_1loop.hpp"
#include "gm2calc/MSSMNoFV_onshell.hpp"

#include "MSSMNoFV/gm2_1loop_helpers.hpp"
#include "gm2_ffunctions.hpp"
#include "gm2_numerics.hpp"

#include <cmath>
#include <complex>
#include <limits>

/**
 * \file gm2_1loop.cpp
 *
 * Contains functions necessary to calculate the SUSY contributions
 * for g-2 at the 1-loop level.
 */

namespace gm2calc {

namespace {

constexpr double Pi = 3.141592653589793;
constexpr double oneOver16PiSqr = 6.332573977646111e-3; // 1/(4 Pi)^2
constexpr double root2 = 1.414213562373095; // Sqrt[2]
constexpr double eps = std::numeric_limits<double>::epsilon();

} // anonymous namespace

/**
 * Calculates full 1-loop SUSY contribution to (g-2), Eq (45) of
 * arXiv:hep-ph/0609168.
 *
 * This function re-defines the muon Yukawa coupling in terms of the
 * tree-level relation with the muon pole mass, i.e. \f$y_\mu =
 * \frac{\sqrt{2} m_\mu^\text{pole}}{v_d}\f$.  Therefore, this
 * function does not use tan(beta) resummation.
 */
double calculate_amu_1loop_non_tan_beta_resummed(const MSSMNoFV_onshell& model)
{
   // create new model and reset its muon Yukawa coupling to tree-level
   MSSMNoFV_onshell model_ytree(model);
   model_ytree.convert_to_non_tan_beta_resummed();

   return amu1LChi0(model_ytree) + amu1LChipm(model_ytree);
}

/**
 * Calculates full 1-loop SUSY contribution to (g-2), Eq (45) of
 * arXiv:hep-ph/0609168.
 *
 * This function assumes that the Yukawa coupling is defined according
 * to arXiv:1504.05500, Eq. (13) and footnote 2.  Therefore, this
 * function uses tan(beta) resummation in the Yukawa coupling.
 */
double calculate_amu_1loop(const MSSMNoFV_onshell& model)
{
   return amu1LChi0(model) + amu1LChipm(model);
}

/**
 * Calculates 1-loop neutralino contribution to (g-2), Eq (2.11a) of
 * arXiv:1311.1775.
 */
double amu1LChi0(const MSSMNoFV_onshell& model)
{
   const Eigen::Array<double,2,1>& m_smu(model.get_MSm());
   const Eigen::Array<double,4,2> AAN_(AAN(model));
   const Eigen::Array<double,4,2> BBN_(BBN(model));
   const Eigen::Array<double,4,2> x(x_im(model));
   const Eigen::Array<double,4,1>& m_chi(model.get_MChi());

   double result = 0;

   for (int i = 0; i < 4; ++i) {
      for (int m = 0; m < 2; ++m) {
         result += - AAN_(i, m) * F1N(x(i, m)) / (12 * sqr(m_smu(m)))
                   - m_chi(i) * BBN_(i, m) * F2N(x(i, m))
                      / (6 * model.get_MM() * sqr(m_smu(m)));
      }
   }

   return result * sqr(model.get_MM()) * oneOver16PiSqr;
}

/**
 * Calculates 1-loop chargino contribution to (g-2), Eq (2.11b) of
 * arXiv:1311.1775.
 */
double amu1LChipm(const MSSMNoFV_onshell& model)
{
   const Eigen::Array<double,2,1> x(x_k(model));
   const double MSvm(model.get_MSvmL());
   const Eigen::Array<double,2,1> AAC_(AAC(model));
   const Eigen::Array<double,2,1> BBC_(BBC(model));
   const Eigen::Array<double,2,1>& m_cha(model.get_MCha());

   double result = 0;

   for (int k = 0; k < 2; ++k) {
      result +=
         AAC_(k) * F1C(x(k)) / 12 +
         m_cha(k) * BBC_(k) * F2C(x(k)) / (3 * model.get_MM());
   }

   return result * sqr(model.get_MM()/MSvm) * oneOver16PiSqr;
}

/**
 * Calculates \f$n^L_{i\tilde{l}_k}\f$, Eq (2.5o) of arXiv:1311.1775,
 * for \f$l=\mu\f$
 */
Eigen::Array<std::complex<double>,4,2> n_L(const MSSMNoFV_onshell& model)
{
   const double gY(model.get_gY());
   const double g2(model.get_g2());
   const double ymu(model.get_Ye(1, 1));
   const Eigen::Matrix<std::complex<double>,4,4>& ZN(model.get_ZN());
   const Eigen::Matrix<double,2,2>& u_smu(model.get_USm());

   Eigen::Array<std::complex<double>,4,2> result;

   for (int i = 0; i < 4; ++i) {
      for (int m = 0; m < 2; ++m) {
         result(i, m) =
            1 / root2 * (+ gY * std::conj(ZN(i,0))
                         + g2 * std::conj(ZN(i,1))) * u_smu(m,0)
            - ymu * std::conj(ZN(i,2)) * u_smu(m,1);
      }
   }

   return result;
}

/**
 * Calculates \f$n^R_{i\tilde{l}_k}\f$, Eq (2.5p) of arXiv:1311.1775,
 * with \f$l=\mu\f$.
 */
Eigen::Array<std::complex<double>,4,2> n_R(const MSSMNoFV_onshell& model)
{
   const double gY(model.get_gY());
   const double ymu(model.get_Ye(1, 1));
   const Eigen::Matrix<std::complex<double>,4,4>& ZN(model.get_ZN());
   const Eigen::Matrix<double,2,2>& u_smu(model.get_USm());

   Eigen::Array<std::complex<double>,4,2> result;

   for (int i = 0; i < 4; ++i) {
      for (int m = 0; m < 2; ++m) {
         result(i, m) = -(root2 * gY * ZN(i, 0) * u_smu(m, 1) +
                          ymu * ZN(i, 2) * u_smu(m, 0));
      }
   }

   return result;
}

/**
 * Calculates \f$c^L_{i\tilde{\nu}_\mu}\f$, Eq (2.5a) of
 * arXiv:1311.1775 for \f$l=\mu\f$.
 *
 * This expression is the complex conjugate of Eq. (50) of
 * arXiv:hep-ph/0609168.
 */
Eigen::Array<std::complex<double>,2,1> c_L(const MSSMNoFV_onshell& model)
{
   return -model.get_g2() * model.get_UP().col(0).conjugate();
}

/**
 * Calculates \f$c^R_{i\tilde{\nu}_\mu}\f$, Eq (2.5b) of
 * arXiv:1311.1775 for \f$l=\mu\f$.
 *
 * This expression is equal to Eq. (51) of arXiv:hep-ph/0609168.
 */
Eigen::Array<std::complex<double>,2,1> c_R(const MSSMNoFV_onshell& model)
{
   return model.get_Ye(1, 1) * model.get_UM().col(1);
}

/**
 * Calculates \f$\mathcal{A}^{c+}_{ii\tilde{\nu}_\mu} =
 * |c^L_{i\tilde{\nu}_\mu}|^2 + |c^R_{i\tilde{\nu}_\mu}|^2\f$,
 * Eqs (2.11b), (2.7a) in arXiv:1311.1775.
 *
 * This expression is identical to \f$|c^L_{im}|^2 + |c^R_{im}|^2\f$,
 * which appears in Eq. (47) of arXiv:hep-ph/0609168.
 */
Eigen::Array<double,2,1> AAC(const MSSMNoFV_onshell& model)
{
   return c_L(model).cwiseAbs2() + c_R(model).cwiseAbs2();
}

/**
 * Calculates \f$\mathcal{A}^{n+}_{ii\tilde{\mu}_m} =
 * |n^L_{i\tilde{\mu}_m}|^2 + |n^R_{i\tilde{\mu}_m}|^2\f$,
 * Eqs (2.11a), (2.7a) in arXiv:1311.1775.
 *
 * This expression is identical to \f$|n^L_{im}|^2 + |n^R_{im}|^2\f$,
 * which appears in Eq. (46) of arXiv:hep-ph/0609168.
 */
Eigen::Array<double,4,2> AAN(const MSSMNoFV_onshell& model)
{
   return n_L(model).cwiseAbs2() + n_R(model).cwiseAbs2();
}

/**
 * Calculates \f$\mathcal{B}^{c+}_{ii\tilde{\nu}_\mu} = 2 \text{Re}
 * [(c^L_{i\tilde{\nu}_\mu})^* c^R_{i\tilde{\nu}_\mu}]\f$,
 * Eqs (2.11b), (2.7b) in arXiv:1311.1775.
 */
Eigen::Array<double,2,1> BBC(const MSSMNoFV_onshell& model)
{
   return 2 * (c_L(model).conjugate() * c_R(model)).real();
}

/**
 * Calculates \f$\mathcal{B}^{n+}_{ii\tilde{\mu}_m} = 2 \text{Re}
 * [(n^L_{i\tilde{\mu}_m})^* n^R_{i\tilde{\mu}_m}]\f$, Eqs (2.11a),
 * (2.7b) in arXiv:1311.1775.
 */
Eigen::Array<double,4,2> BBN(const MSSMNoFV_onshell& model)
{
   return 2 * (n_L(model).conjugate() * n_R(model)).real();
}

/**
 * Calculates \f$x_{im} = m^2_{\chi_i^0} / m^2_{\tilde{\mu}_m}\f$,
 * which appears in Eq (46) of arXiv:hep-ph/0609168.
 */
Eigen::Array<double,4,2> x_im(const MSSMNoFV_onshell& model)
{
   const Eigen::Matrix<double,4,1> m_chi(model.get_MChi());
   const Eigen::Matrix<double,2,1> m_sm(model.get_MSm());

   return (m_chi * m_sm.cwiseInverse().transpose()).array().square();
}

/**
 * Calculates \f$x_k = m^2_{\chi_k^\pm} / m^2_{\tilde{\nu}_\mu}\f$,
 * which appears in Eq (47) of arXiv:hep-ph/0609168.
 */
Eigen::Array<double,2,1> x_k(const MSSMNoFV_onshell& model)
{
   return (model.get_MCha() / model.get_MSvmL()).square();
}

// === approximations ===

/**
 * Calculates the 1-loop leading log approximation: Wino--Higgsino,
 * muon-sneutrino, Eq (6.2a) arXiv:1311.1775
 */
double amu1LWHnu(const MSSMNoFV_onshell& model)
{
   const double tan_beta = model.get_TB();
   const double M2 = model.get_MassWB();
   const double mu = model.get_Mu();
   const double msv_2 = model.get_MSvmL();

   return sqr(model.get_g2()) * 2 * oneOver16PiSqr
      * (sqr(model.get_MM()) * M2 * mu * tan_beta) / sqr(sqr(msv_2))
      * Fa(sqr(M2 / msv_2), sqr(mu / msv_2));
}

/**
 * Calculates the 1-loop leading log approximation: Wino--Higgsino,
 * left-handed smuon, Eq (6.2b) arXiv:1311.1775
 */
double amu1LWHmuL(const MSSMNoFV_onshell& model)
{
   const double tan_beta = model.get_TB();
   const double M2 = model.get_MassWB();
   const double mu = model.get_Mu();
   const double msl_2 = std::sqrt(model.get_ml2(1, 1));

   return - sqr(model.get_g2()) * oneOver16PiSqr
      * (sqr(model.get_MM()) * M2 * mu * tan_beta) / sqr(sqr(msl_2))
      * Fb(sqr(M2 / msl_2), sqr(mu / msl_2));
}

/**
 * Calculates the 1-loop leading log approximation: Bino--Higgsino,
 * left-handed smuon, Eq (6.2c) arXiv:1311.1775
 */
double amu1LBHmuL(const MSSMNoFV_onshell& model)
{
   const double tan_beta = model.get_TB();
   const double M1 = model.get_MassB();
   const double mu = model.get_Mu();
   const double msl_2 = std::sqrt(model.get_ml2(1, 1));
   const double gY = model.get_gY();

   return sqr(gY) * oneOver16PiSqr
      * (sqr(model.get_MM()) * M1 * mu * tan_beta) / sqr(sqr(msl_2))
      * Fb(sqr(M1 / msl_2), sqr(mu / msl_2));
}

/**
 * Calculates the 1-loop leading log approximation: Bino--Higgsino,
 * right-handed smuon, Eq (6.2d) arXiv:1311.1775
 */
double amu1LBHmuR(const MSSMNoFV_onshell& model)
{
   const double tan_beta = model.get_TB();
   const double M1 = model.get_MassB();
   const double mu = model.get_Mu();
   const double mse_2 = std::sqrt(model.get_me2(1, 1));
   const double gY = model.get_gY();

   return - sqr(gY) * 2 * oneOver16PiSqr
      * (sqr(model.get_MM()) * M1 * mu * tan_beta) / sqr(sqr(mse_2))
      * Fb(sqr(M1 / mse_2), sqr(mu / mse_2));
}

/**
 * Calculates the 1-loop leading log approximation: Bino, left-handed
 * smuon, right-handed smuon, Eq (6.2e) arXiv:1311.1775
 */
double amu1LBmuLmuR(const MSSMNoFV_onshell& model)
{
   const double tan_beta = model.get_TB();
   const double M1 = model.get_MassB();
   const double mu = model.get_Mu();
   const double msl_2 = std::sqrt(model.get_ml2(1, 1));
   const double mse_2 = std::sqrt(model.get_me2(1, 1));
   const double gY = model.get_gY();

   if (is_zero(M1, eps)) {
      return 0;
   }

   return sqr(gY) * 2 * oneOver16PiSqr
      * (sqr(model.get_MM()) * mu * tan_beta) / (M1 * sqr(M1))
      * Fb(sqr(msl_2 / M1), sqr(mse_2 / M1));
}

/**
 * Calculates the full 1-loop leading log approximation, Eq (6.1)
 * arXiv:1311.1775
 * as it stands, without tan(beta) resummation
 */
double amu1Lapprox_non_tan_beta_resummed(const MSSMNoFV_onshell& model)
{
   return amu1LWHnu(model) + amu1LWHmuL(model) + amu1LBHmuL(model) +
          amu1LBHmuR(model) + amu1LBmuLmuR(model);
}

/**
 * Calculates the full 1-loop leading log approximation, Eq (6.1)
 * arXiv:1311.1775
 * but include tan(beta) resummation
 */
double amu1Lapprox(const MSSMNoFV_onshell& model)
{
   return amu1Lapprox_non_tan_beta_resummed(model) * tan_beta_cor(model);
}

// tan(beta) corrections

/**
 * Calculates \f$\frac{1}{1 + \Delta_{\mu}}\f$
 *
 * @param model model parameters
 *
 * @return \f$\frac{1}{1 + \Delta_{\mu}}\f$
 */
double tan_beta_cor(const MSSMNoFV_onshell& model)
{
   const double delta_mu = delta_mu_correction(model);

   return 1.0 / (1.0 + delta_mu);
}

/**
 * Calculates \f$\Delta_{\ell}\f$ using a generalized form of Eq (8)
 * arxiv:0808.1530.
 *
 * @param model model parameters
 * @param gen lepton generation (0 = electron, 1 = muon, 2 = tau)
 *
 * @return \f$\Delta_{\ell}\f$
 */
double delta_down_lepton_correction(const MSSMNoFV_onshell& model, int gen)
{
   const double mu = model.get_Mu();
   const double tb = model.get_TB();
   const double g2 = model.get_g2();
   const double gY = model.get_gY();
   const double M1 = model.get_MassB();
   const double M2 = model.get_MassWB();
   const double mw2 = sqr(model.get_MW());
   const double mz2 = sqr(model.get_MZ());
   const double sw = std::sqrt(1 - mw2/mz2);
   const double sw2 = sw*sw;
   const double x1 = sqr(M2) + sqr(mu) + 2*mw2;
   const double x2 = abs_sqrt(sqr(x1) - sqr(2*M2*mu));
   const double m1 = abs_sqrt(0.5*(x1 - x2));
   const double m2 = abs_sqrt(0.5*(x1 + x2));
   const double m_sneu_lep = abs_sqrt(model.get_ml2(gen, gen) - 0.5*mz2);
   const double m_slep_L = abs_sqrt(model.get_ml2(gen, gen) - mz2*(sw2 - 0.5));
   const double m_slep_R = abs_sqrt(model.get_me2(gen, gen) + mz2*sw2);

   const double delta_lep =
      - mu*tb*oneOver16PiSqr
       *(+ sqr(g2)*M2*(Iabc(m1, m2, m_sneu_lep) + 0.5*Iabc(m1, m2, m_slep_L))
         + sqr(gY)*M1*(Iabc(mu, M1, m_slep_R) - 0.5*Iabc(mu, M1, m_slep_L)
                       - Iabc(M1, m_slep_L, m_slep_R)));

   return delta_lep;
}

/**
 * Calculates \f$\Delta_{\mu}\f$ using delta_down_lepton_correction()
 *
 * @param model model parameters
 *
 * @return \f$\Delta_{\mu}\f$
 */
double delta_mu_correction(const MSSMNoFV_onshell& model)
{
   return delta_down_lepton_correction(model, 1);
}

/**
 * Calculates \f$\Delta_{\tau}\f$ using delta_down_lepton_correction()
 *
 * @param model model parameters
 *
 * @return \f$\Delta_{\tau}\f$
 */
double delta_tau_correction(const MSSMNoFV_onshell& model)
{
   return delta_down_lepton_correction(model, 2);
}

/**
 * Returns the \f$\Delta_b\f$ corrections from arxiv:0901.2065,
 * Eq (103) and Eqs (31)-(35)
 */
double delta_bottom_correction(const MSSMNoFV_onshell& model)
{
   const double tb = model.get_TB();
   const double At = model.get_Au(2,2);
   const double yt = model.get_Yu(2,2);
   const double gY = model.get_gY();
   const double g2 = model.get_g2();
   const double alpha_S = sqr(model.get_g3()) / (4*Pi);
   const double mu = model.get_Mu();
   const double M1 = model.get_MassB();
   const double M2 = model.get_MassWB();
   const double M3 = model.get_MassG();
   const double mstL = abs_sqrt(model.get_mq2(2,2));
   const double mstR = abs_sqrt(model.get_mu2(2,2));
   const double msbL = abs_sqrt(model.get_mq2(2,2));
   const double msbR = abs_sqrt(model.get_md2(2,2));

   // Note:
   // 1/z^2 H_2(x^2/z^2, y^2/z^2) = - Iabc(x,y,z)

   // Eq.(31)
   const double eps_0 =
        2/(3*Pi) * alpha_S * mu * M3
        * Iabc(msbL, msbR, std::abs(M3))

      - 1/(96*sqr(Pi)) * sqr(gY) * mu * M1
        * (Iabc(msbL, std::abs(mu), std::abs(M1))
           + 2 * Iabc(msbR, std::abs(mu), std::abs(M1)))

      - 1/(144*sqr(Pi)) * sqr(gY) * mu * M1
        * Iabc(msbL, msbR, std::abs(M1))

      - 3/(32*sqr(Pi)) * sqr(g2) * mu * M2
        * Iabc(msbL, std::abs(mu), std::abs(M2))
      ;

   const double Ittm = Iabc(mstL, mstR, std::abs(mu));

   // Eq.(35), term ~ self-energy neglected
   const double eps_Y_vM = - oneOver16PiSqr * At * mu * Ittm;

   // Eq.(32)
   const double eps_Y = oneOver16PiSqr * At * mu * Ittm + eps_Y_vM;

   // Eq.(33)
   const double eps_3_tilde = eps_0 + sqr(yt) * eps_Y;

   // Eq.(103)
   const double delta_b = tb * eps_3_tilde;

   return delta_b;
}

} // namespace gm2calc
