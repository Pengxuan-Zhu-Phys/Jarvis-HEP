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
 * @file THDM_mass_eigenstates.hpp
 *
 * @brief contains class for general THDM model with routines needed
 *        to solve EWSB and determine the masses and mixings.
 *
 * This file was generated with FlexibleSUSY 2.6.0 and SARAH 4.14.3 .
 */

#ifndef GM2_THDM_MASS_EIGENSTATES_H
#define GM2_THDM_MASS_EIGENSTATES_H

#include "THDM_problems.hpp"
#include "THDM_parameters.hpp"

#include <iosfwd>

#include <Eigen/Core>

namespace gm2calc {

/**
 * @class THDM_mass_eigenstates
 * @brief model class with routines for determing masses and mixinga and EWSB
 */
class THDM_mass_eigenstates : public THDM_parameters
{
public:
   void print(std::ostream&) const;

   void calculate_MSbar_masses();

   void do_force_output(bool);
   bool do_force_output() const;
   void reorder_MSbar_masses();
   const THDM_problems& get_problems() const;
   int solve_ewsb_tree_level();
   int solve_ewsb();

   /// gluon mass
   double get_MVG() const { return MVG; }
   /// photon mass
   double get_MVP() const { return MVP; }
   /// W boson mass
   double get_MVWm() const { return MVWm; }
   /// Z boson mass
   double get_MVZ() const { return MVZ; }
   /// CP-even Higgs boson masses
   const Eigen::Array<double,2,1>& get_Mhh() const { return Mhh; }
   /// CP-even Higgs boson i mass
   double get_Mhh(int i) const { return Mhh(i); }
   /// Goldstone and CP-odd Higgs boson masses (in that order)
   const Eigen::Array<double,2,1>& get_MAh() const { return MAh; }
   /// Goldstone (i = 0) or CP-odd Higgs boson (i = 1) mass
   double get_MAh(int i) const { return MAh(i); }
   /// Goldstone and charged Higgs boson masses (in that order)
   const Eigen::Array<double,2,1>& get_MHm() const { return MHm; }
   /// Goldstone (i = 0) or charged Higgs boson (i = 1) mass
   double get_MHm(int i) const { return MHm(i); }
   /// up-type quark masses
   const Eigen::Array<double,3,1>& get_MFu() const { return MFu; }
   /// up-type quark i mass
   double get_MFu(int i) const { return MFu(i); }
   /// down-type quark masses
   const Eigen::Array<double,3,1>& get_MFd() const { return MFd; }
   /// down-type quark i mass
   double get_MFd(int i) const { return MFd(i); }
   /// neutrino masses
   const Eigen::Array<double,3,1>& get_MFv() const { return MFv; }
   /// neutrino i mass
   double get_MFv(int i) const { return MFv(i); }
   /// charged lepton masses
   const Eigen::Array<double,3,1>& get_MFe() const { return MFe; }
   /// charged lepton i mass
   double get_MFe(int i) const { return MFe(i); }

   /// CP-even Higgs boson mixing matrix
   const Eigen::Matrix<double,2,2>& get_ZH() const { return ZH; }
   /// CP-even Higgs boson mixing matrix element
   double get_ZH(int i, int k) const { return ZH(i,k); }
   /// CP-odd Higgs boson mixing matrix
   const Eigen::Matrix<double,2,2>& get_ZA() const { return ZA; }
   /// CP-odd Higgs boson mixing matrix element
   double get_ZA(int i, int k) const { return ZA(i,k); }
   /// charged Higgs boson mixing matrix
   const Eigen::Matrix<double,2,2>& get_ZP() const { return ZP; }
   /// charged Higgs boson mixing matrix element
   double get_ZP(int i, int k) const { return ZP(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Vd() const { return Vd; }
   std::complex<double> get_Vd(int i, int k) const { return Vd(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ud() const { return Ud; }
   std::complex<double> get_Ud(int i, int k) const { return Ud(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Vu() const { return Vu; }
   std::complex<double> get_Vu(int i, int k) const { return Vu(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Uu() const { return Uu; }
   std::complex<double> get_Uu(int i, int k) const { return Uu(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ve() const { return Ve; }
   std::complex<double> get_Ve(int i, int k) const { return Ve(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ue() const { return Ue; }
   std::complex<double> get_Ue(int i, int k) const { return Ue(i,k); }

   double get_ewsb_eq_hh_1() const;
   double get_ewsb_eq_hh_2() const;

   double get_sin_beta() const;   ///< sin(beta)
   double get_cos_beta() const;   ///< cos(beta)
   double get_tan_beta() const;   ///< tan(beta) = ratio of VEVs v2/v1
   double get_beta() const;       ///< CP-odd and charged Higgs mixing angle
   double get_alpha_h() const;    ///< CP-even Higgs mixing angle
   double get_sin_beta_minus_alpha() const; ///< sin(beta - alpha_h)
   double get_cos_beta_minus_alpha() const; ///< cos(beta - alpha_h)
   double get_alpha_em() const;   ///< electromagnetic coupling
   double get_eta() const;        ///< deviation of CP-even Higgs mixing angle from SM limit
   double get_LambdaFive() const; ///< capital Lambda5, Eq (14) arxiv:1607.06292
   double get_LambdaSixSeven() const; ///< (Lambda_{567} - Lambda_{5})(tan(b) - 1/tan(b))
   double ThetaW() const;         ///< weak mixing angle
   double get_v() const;          ///< SM-like VEV
   double get_v_sqr() const;      ///< squared SM-like VEV

   /// set tan(beta) and vacuum expectation value
   void set_tan_beta_and_v(double, double);

   /// set alpha_em and cos(theta_w)
   void set_alpha_em_and_cw(double, double);

protected:
   void calculate_boson_masses();
   void calculate_fermion_masses();

   double get_mass_matrix_VZ() const;
   void calculate_MVZ();
   double get_mass_matrix_VWm() const;
   void calculate_MVWm();
   Eigen::Matrix<double,3,3> get_mass_matrix_Fv() const;
   void calculate_MFv();
   Eigen::Matrix<double,2,2> get_mass_matrix_hh() const;
   void calculate_Mhh();
   Eigen::Matrix<double,2,2> get_mass_matrix_Ah() const;
   void calculate_MAh();
   Eigen::Matrix<double,2,2> get_mass_matrix_Hm() const;
   void calculate_MHm();
   Eigen::Matrix<std::complex<double>,3,3> get_mass_matrix_Fd() const;
   void calculate_MFd();
   Eigen::Matrix<std::complex<double>,3,3> get_mass_matrix_Fu() const;
   void calculate_MFu();
   Eigen::Matrix<std::complex<double>,3,3> get_mass_matrix_Fe() const;
   void calculate_MFe();

private:
   bool force_output{false};        ///< switch to force output of pole masses
   THDM_problems problems{}; ///< problems

   // masses
   double MVG{0.0};
   double MVWm{0.0};
   double MVP{0.0};
   double MVZ{0.0};
   Eigen::Array<double,3,1> MFv{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,2,1> Mhh{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MAh{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MHm{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,3,1> MFd{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> MFu{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> MFe{Eigen::Array<double,3,1>::Zero()};

   // mixing matrices
   Eigen::Matrix<double,2,2> ZH{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZA{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZP{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Vd{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ud{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Vu{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Uu{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ve{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ue{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
};

std::ostream& operator<<(std::ostream&, const THDM_mass_eigenstates&);

} // namespace gm2calc

#endif
