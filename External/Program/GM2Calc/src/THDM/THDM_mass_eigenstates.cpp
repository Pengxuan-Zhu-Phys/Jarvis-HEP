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

/**
 * @file THDM_mass_eigenstates.cpp
 * @brief implementation of the THDM model class
 *
 * Contains the definition of the THDM model class methods
 * which solve EWSB and calculate masses and mixings from MSbar
 * parameters.
 *
 * This file was generated with FlexibleSUSY 2.6.0 and SARAH 4.14.3.
 */

#include "gm2calc/THDM_mass_eigenstates.hpp"

#include "gm2_eigen_utils.hpp"
#include "gm2_linalg.hpp"
#include "gm2_numerics.hpp"
#include "gm2_raii.hpp"

#include <cmath>
#include <iostream>

namespace gm2calc {

namespace {

const double pi = 3.1415926535897932;
const double sqrt2_inv = 0.70710678118654752; // 1/Sqrt[2]
const double gut_normalization = 0.77459666924148338; // Sqrt[3/5]

} // anonymous namespace

#define CLASSNAME THDM_mass_eigenstates

void CLASSNAME::do_force_output(bool flag)
{
   force_output = flag;
}

bool CLASSNAME::do_force_output() const
{
   return force_output;
}

const THDM_problems& CLASSNAME::get_problems() const
{
   return problems;
}

int CLASSNAME::solve_ewsb_tree_level()
{
   int error = 0;

   const double old_m112 = m112;
   const double old_m222 = m222;

   m112 = (0.25*(2*m122*v2 + 2*v2*m122 - 2*lambda1*cube(v1) - lambda7*
      cube(v2) - lambda7*cube(v2) - 3*lambda6*v2*sqr(v1) - 3*v2*lambda6*sqr(v1)
      - 2*lambda3*v1*sqr(v2) - 2*lambda4*v1*sqr(v2) - lambda5*v1*sqr(v2)
      - v1*lambda5*sqr(v2)))/v1;
   m222 = (0.25*(2*m122*v1 + 2*v1*m122 - lambda6*cube(v1) - lambda6
      *cube(v1) - 2*lambda2*cube(v2) - 2*lambda3*v2*sqr(v1) - 2*lambda4*v2*sqr(v1)
      - lambda5*v2*sqr(v1) - v2*lambda5*sqr(v1) - 3*lambda7*v1*sqr(v2) - 3*
      v1*lambda7*sqr(v2)))/v2;

   const bool is_finite = std::isfinite(m112) && std::isfinite(m222);

   if (!is_finite) {
      m112 = old_m112;
      m222 = old_m222;
      error = 1;
   }

   return error;
}

int CLASSNAME::solve_ewsb()
{
   return solve_ewsb_tree_level();
}

void CLASSNAME::print(std::ostream& ostr) const
{
   ostr << "========================================\n"
           "THDM\n"
           "========================================\n";
   THDM_parameters::print(ostr);
   ostr << "----------------------------------------\n"
           "Masses:\n"
           "----------------------------------------\n";
   ostr << "Mhh = {" << Mhh(0) << ", " << Mhh(1) << "} GeV\n";
   ostr << "MAh = {" << MAh(0) << ", " << MAh(1) << "} GeV\n";
   ostr << "MHm = {" << MHm(0) << ", " << MHm(1) << "} GeV\n";
   ostr << "MFu = {" << MFu(0) << ", " << MFu(1) << ", " << MFu(2) << "} GeV\n";
   ostr << "MFd = {" << MFd(0) << ", " << MFd(1) << ", " << MFd(2) << "} GeV\n";
   ostr << "MFv = {" << MFv(0) << ", " << MFv(1) << ", " << MFv(2) << "} GeV\n";
   ostr << "MFe = {" << MFe(0) << ", " << MFe(1) << ", " << MFe(2) << "} GeV\n";
   ostr << "MVWm = " << MVWm << " GeV\n";
   ostr << "MVZ = " << MVZ << " GeV\n";

   ostr << "----------------------------------------\n"
           "Mixing matrices:\n"
           "----------------------------------------\n";
   ostr << "ZH =\n" << ZH << '\n';
   ostr << "ZA =\n" << ZA << '\n';
   ostr << "ZP =\n" << ZP << '\n';
   ostr << "Vd =\n" << Vd << '\n';
   ostr << "Ud =\n" << Ud << '\n';
   ostr << "Vu =\n" << Vu << '\n';
   ostr << "Uu =\n" << Uu << '\n';
   ostr << "Ve =\n" << Ve << '\n';
   ostr << "Ue =\n" << Ue << '\n';

   ostr << "----------------------------------------\n"
           "Derived parameters:\n"
           "----------------------------------------\n";
   ostr << "v = " << get_v() << " GeV\n";
   ostr << "theta_w = " << ThetaW() << '\n';
   ostr << "alpha_h = " << get_alpha_h() << '\n';
   ostr << "beta = " << get_beta() << '\n';
   ostr << "sin(beta - alpha_h) = " << get_sin_beta_minus_alpha() << '\n';
   ostr << "cos(beta - alpha_h) = " << get_cos_beta_minus_alpha() << '\n';
   ostr << "eta = " << get_eta() << '\n';
   ostr << "tan(beta) = " << get_tan_beta() << '\n';
}

/**
 * routine which finds the MSbar mass eigenstates and mixings.
 */
void CLASSNAME::calculate_MSbar_masses()
{
   calculate_boson_masses();
   calculate_fermion_masses();
}

/**
 * routine which finds the boson mass eigenstates and mixings.
 */
void CLASSNAME::calculate_boson_masses()
{
   const auto save_m112_raii = make_raii_save(m112);
   const auto save_m222_raii = make_raii_save(m222);

   solve_ewsb_tree_level();

   calculate_MVZ();
   calculate_MVWm();
   calculate_MHm();
   calculate_MAh();
   calculate_Mhh();

   reorder_MSbar_masses();
}

/**
 * routine which finds the fermion mass eigenstates and mixings.
 */
void CLASSNAME::calculate_fermion_masses()
{
   calculate_MFu();
   calculate_MFd();
   calculate_MFv();
   calculate_MFe();
}

/**
 * reorders MSbar masses so that golstones are placed at the index
 * specified in the model files definition of the associated
 * gauge boson (see Z-boson definition in default particles.m file
 * in the Models directory of your SARAH distribution for example)
 */
void CLASSNAME::reorder_MSbar_masses()
{
   move_goldstone_to(0, MVZ, MAh, ZA);
   move_goldstone_to(0, MVWm, MHm, ZP);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fv() const
{
   return Eigen::Matrix<double,3,3>::Zero();
}

void CLASSNAME::calculate_MFv()
{
   MFv.setZero();
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_hh() const
{
   Eigen::Matrix<double,2,2> mass_matrix_hh;

   mass_matrix_hh(0,0) = m112 + 1.5*lambda6*v1*v2 + 1.5*v1*v2*lambda6 +
      1.5*lambda1*sqr(v1) + 0.5*lambda3*sqr(v2) + 0.5*lambda4*sqr(v2) + 0.25*
      lambda5*sqr(v2) + 0.25*lambda5*sqr(v2);
   mass_matrix_hh(0,1) = -0.5*m122 + lambda3*v1*v2 + lambda4*v1*v2 + 0.5*
      lambda5*v1*v2 + 0.5*v1*v2*lambda5 - 0.5*m122 + 0.75*lambda6*
      sqr(v1) + 0.75*lambda6*sqr(v1) + 0.75*lambda7*sqr(v2) + 0.75*
      lambda7*sqr(v2);
   mass_matrix_hh(1,1) = m222 + 1.5*lambda7*v1*v2 + 1.5*v1*v2*lambda7 +
      0.5*lambda3*sqr(v1) + 0.5*lambda4*sqr(v1) + 0.25*lambda5*sqr(v1) + 0.25*
      lambda5*sqr(v1) + 1.5*lambda2*sqr(v2);

   symmetrize(mass_matrix_hh);

   return mass_matrix_hh;
}

void CLASSNAME::calculate_Mhh()
{
   const auto mass_matrix_hh(get_mass_matrix_hh());
   fs_diagonalize_hermitian<double,double,2>(mass_matrix_hh, Mhh, ZH);
   normalize_to_interval<2,2>(ZH);

   if (Mhh.minCoeff() < 0.) {
      problems.flag_tachyon("hh");
   }

   Mhh = sqrt(Mhh.cwiseAbs());
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Ah() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Ah;

   const double tw = gut_normalization*g1/g2;
   const double rt = std::sqrt(1 + tw*tw);
   const double sw = tw/rt;
   const double cw = 1/rt;

   mass_matrix_Ah(0,0) = m112 + 0.5*lambda6*v1*v2 + 0.5*v1*v2*lambda6 +
      0.5*lambda1*sqr(v1) + 0.3872983346207417*g1*g2*cw*sw*sqr(v1) +
      0.5*lambda3*sqr(v2) + 0.5*lambda4*sqr(v2) - 0.25*lambda5*sqr(
      v2) - 0.25*lambda5*sqr(v2) + 0.25*sqr(g2)*sqr(v1)*sqr(cw
      ) + 0.15*sqr(g1)*sqr(v1)*sqr(sw);
   mass_matrix_Ah(0,1) = -0.5*m122 + 0.5*lambda5*v1*v2 + 0.5*v1*v2*lambda5
      - 0.5*m122 + 0.3872983346207417*g1*g2*v1*v2*cw*sw + 0.25*lambda6*sqr(v1)
      + 0.25*lambda6*sqr(v1) + 0.25*
      lambda7*sqr(v2) + 0.25*lambda7*sqr(v2) + 0.25*v1*v2*sqr(g2)*sqr(cw) +
      0.15*v1*v2*sqr(g1)*sqr(sw);
   mass_matrix_Ah(1,1) = m222 + 0.5*lambda7*v1*v2 + 0.5*v1*v2*lambda7 +
      0.5*lambda3*sqr(v1) + 0.5*lambda4*sqr(v1) - 0.25*lambda5*sqr(v1) - 0.25*
      lambda5*sqr(v1) + 0.5*lambda2*sqr(v2) + 0.3872983346207417*g1*g2*
      cw*sw*sqr(v2) + 0.25*sqr(g2)*sqr(v2)*sqr(cw) + 0.15*sqr(g1)*sqr(v2)*sqr(sw);

   symmetrize(mass_matrix_Ah);

   return mass_matrix_Ah;
}

void CLASSNAME::calculate_MAh()
{
   const auto mass_matrix_Ah(get_mass_matrix_Ah());
   fs_diagonalize_hermitian<double,double,2>(mass_matrix_Ah, MAh, ZA);
   normalize_to_interval<2,2>(ZA);

   if (MAh.minCoeff() < 0.) {
      problems.flag_tachyon("Ah");
   }

   MAh = sqrt(MAh.cwiseAbs());
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Hm() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Hm;

   mass_matrix_Hm(0,0) = m112 + 0.5*lambda6*v1*v2 + 0.5*v1*v2*lambda6 +
      0.5*lambda1*sqr(v1) + 0.25*sqr(g2)*sqr(v1) + 0.5*lambda3*sqr(v2);
   mass_matrix_Hm(0,1) = 0.5*lambda4*v1*v2 + 0.5*lambda5*v1*v2 - m122 +
      0.25*v1*v2*sqr(g2) + 0.5*lambda6*sqr(v1) + 0.5*lambda7*sqr(v2
      );
   mass_matrix_Hm(1,0) = -m122 + 0.5*lambda4*v1*v2 + 0.5*v1*v2*lambda5 +
      0.25*v1*v2*sqr(g2) + 0.5*lambda6*sqr(v1) + 0.5*lambda7*sqr(v2);
   mass_matrix_Hm(1,1) = m222 + 0.5*lambda7*v1*v2 + 0.5*v1*v2*lambda7 +
      0.5*lambda3*sqr(v1) + 0.5*lambda2*sqr(v2) + 0.25*sqr(g2)*sqr(v2);

   return mass_matrix_Hm;
}

void CLASSNAME::calculate_MHm()
{
   const auto mass_matrix_Hm(get_mass_matrix_Hm());
   fs_diagonalize_hermitian<double,double,2>(mass_matrix_Hm, MHm, ZP);
   normalize_to_interval<2,2>(ZP);

   if (MHm.minCoeff() < 0.) {
      problems.flag_tachyon("Hm");
   }

   MHm = sqrt(MHm.cwiseAbs());
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::get_mass_matrix_Fd() const
{
   return sqrt2_inv*(v1*Gamma_d + v2*Pi_d);
}

void CLASSNAME::calculate_MFd()
{
   const auto mass_matrix_Fd(get_mass_matrix_Fd());
   fs_svd<double,std::complex<double>,3,3>(mass_matrix_Fd, MFd, Vd, Ud);
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::get_mass_matrix_Fu() const
{
   return sqrt2_inv*(v1*Gamma_u + v2*Pi_u);
}

void CLASSNAME::calculate_MFu()
{
   const auto mass_matrix_Fu(get_mass_matrix_Fu());
   fs_svd<double,std::complex<double>,3,3>(mass_matrix_Fu, MFu, Vu, Uu);
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::get_mass_matrix_Fe() const
{
   return sqrt2_inv*(v1*Gamma_l + v2*Pi_l);
}

void CLASSNAME::calculate_MFe()
{
   const auto mass_matrix_Fe(get_mass_matrix_Fe());
   fs_svd<double,std::complex<double>,3,3>(mass_matrix_Fe, MFe, Ve, Ue);
}

double CLASSNAME::get_mass_matrix_VWm() const
{
   const double mass_matrix_VWm = 0.25*sqr(g2)*(sqr(v1) + sqr(v2));
   return mass_matrix_VWm;
}

void CLASSNAME::calculate_MVWm()
{
   const auto mass_matrix_VWm = get_mass_matrix_VWm();
   MVWm = abs_sqrt(mass_matrix_VWm);
}

double CLASSNAME::get_mass_matrix_VZ() const
{
   const double tw = gut_normalization*g1/g2;
   const double rt = std::sqrt(1 + tw*tw);
   const double sw = tw/rt;
   const double cw = 1/rt;

   const double mass_matrix_VZ = 0.25*(sqr(v1) + sqr(v2))
      *sqr(g2*cw + gut_normalization*g1*sw);

   return mass_matrix_VZ;
}

void CLASSNAME::calculate_MVZ()
{
   const auto mass_matrix_VZ = get_mass_matrix_VZ();
   MVZ = abs_sqrt(mass_matrix_VZ);
}

double CLASSNAME::get_ewsb_eq_hh_1() const
{
   double result = m112*v1 - 0.5*m122*v2 - 0.5*v2*m122 + 0.5*lambda1*cube
      (v1) + 0.25*lambda7*cube(v2) + 0.25*lambda7*cube(v2) + 0.75*lambda6*v2
      *sqr(v1) + 0.75*v2*lambda6*sqr(v1) + 0.5*lambda3*v1*sqr(v2) + 0.5*
      lambda4*v1*sqr(v2) + 0.25*lambda5*v1*sqr(v2) + 0.25*v1*lambda5*sqr(v2);

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_2() const
{
   double result = -0.5*m122*v1 + m222*v2 - 0.5*v1*m122 + 0.25*lambda6*
      cube(v1) + 0.25*lambda6*cube(v1) + 0.5*lambda2*cube(v2) + 0.5*lambda3*
      v2*sqr(v1) + 0.5*lambda4*v2*sqr(v1) + 0.25*lambda5*v2*sqr(v1) + 0.25*v2*
      lambda5*sqr(v1) + 0.75*lambda7*v1*sqr(v2) + 0.75*v1*lambda7*sqr(v2);

   return result;
}

/**
 * Returns the CP-odd Higgs mixing angle \f$\beta\f$.
 *
 * The mixing angle \f$\beta\f$ is chosen such that
 * \f$0 \leq \beta \leq \pi/2\f$.
 *
 * @return CP-odd Higgs mixing angle \f$\beta\f$
 */
double CLASSNAME::get_beta() const
{
   return std::atan(get_tan_beta());
}

double CLASSNAME::get_sin_beta() const
{
   const double tb = get_tan_beta();
   return tb/std::sqrt(1.0 + sqr(tb));
}

double CLASSNAME::get_cos_beta() const
{
   const double tb = get_tan_beta();
   return 1.0/std::sqrt(1.0 + sqr(tb));
}

double CLASSNAME::get_tan_beta() const
{
   return v2/v1;
}

/**
 * Returns the CP-even Higgs mixing angle \f$\alpha_h\f$.
 *
 * The mixing angle \f$\alpha_h\f$ is chosen such that
 * \f$-\pi/2 \leq \beta - \alpha_h \leq \pi/2\f$.
 *
 * @return CP-even Higgs mixing angle \f$\alpha_h\f$
 */
double CLASSNAME::get_alpha_h() const
{
   double alpha_h = std::asin(ZH(1,1));

   const double bma = get_beta() - alpha_h;
   const double eps = 10*std::numeric_limits<double>::epsilon();

   if (bma < -pi/2 - eps) { alpha_h -= pi; }
   if (bma >  pi/2 + eps) { alpha_h += pi; }

   return alpha_h;
}

double CLASSNAME::get_sin_beta_minus_alpha() const
{
   return std::sin(get_beta() - get_alpha_h());
}

double CLASSNAME::get_cos_beta_minus_alpha() const
{
   return std::cos(get_beta() - get_alpha_h());
}

double CLASSNAME::get_alpha_em() const
{
   const double gY = g1*gut_normalization;
   const double e2 = sqr(gY)*sqr(g2)/(sqr(gY) + sqr(g2));
   return e2/12.566370614359173; // e^2/(4 Pi)
}

/**
 * Returns the CP-even Higgs mixing angle contribution \f$\eta\f$ that
 * describes the deviation from the SM limit.
 *
 * The \f$\eta\f$ parameter is defined as
 * \f$\beta - \alpha = \pi/2 - \eta\f$.
 *
 * @return \f$\eta\f$
 */
double CLASSNAME::get_eta() const
{
   return pi/2 + get_alpha_h() - get_beta();
}

/**
 * Returns \f$\Lambda_5\f$, Eq (14) arxiv:1607.06292
 * @return \f$\Lambda_5\f$
 */
double CLASSNAME::get_LambdaFive() const
{
   const double tb = get_tan_beta();
   const double sbcb = tb/(1.0 + sqr(tb)); // Sin[beta] Cos[beta]
   return 2*get_m122()/(get_v_sqr() * sbcb);
}

/**
 * Returns \f$\left(\Lambda_{567} - \Lambda_{5}\right)\left(\tan\beta - \frac{1}{\tan\beta}\right) = \frac{\lambda_6}{\sin^2\beta} - \frac{\lambda_7}{\cos^2\beta}\f$
 */
double CLASSNAME::get_LambdaSixSeven() const
{
   const double sb = get_sin_beta();
   const double cb = get_cos_beta();
   return lambda6/sqr(sb) - lambda7/sqr(cb);
}

double CLASSNAME::ThetaW() const
{
   return std::atan(gut_normalization*g1/g2);
}

double CLASSNAME::get_v() const
{
   return std::sqrt(get_v_sqr());
}

double CLASSNAME::get_v_sqr() const
{
   return sqr(v1) + sqr(v2);
}

void CLASSNAME::set_alpha_em_and_cw(double alpha_em, double cw)
{
   const double e = std::sqrt(alpha_em*4*pi);
   const double sw = std::sqrt(1 - sqr(cw));
   g1 = e/(cw*gut_normalization);
   g2 = e/sw;
}

void CLASSNAME::set_tan_beta_and_v(double tan_beta, double v)
{
   const double sq = std::sqrt(1.0 + sqr(tan_beta));
   const double sb = tan_beta/sq;
   const double cb = 1/sq;
   v1 = v*cb;
   v2 = v*sb;
}

std::ostream& operator<<(std::ostream& ostr, const CLASSNAME& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace gm2calc
