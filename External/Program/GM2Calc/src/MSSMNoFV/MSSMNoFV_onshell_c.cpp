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

#include "gm2calc/MSSMNoFV_onshell.h"
#include "gm2calc/MSSMNoFV_onshell.hpp"
#include "gm2calc/gm2_error.hpp"

#include "gm2_log.hpp"

#include <complex>
#include <string>
#include <iostream>

/**
 * @file MSSMNoFV_onshell_c.cpp
 * @brief contains definitions of C interface functions for the model
 *
 * This file contains the definitions for the C interface functions
 * used to set and retrieve the model parameters and masses.
 */

extern "C"
{

/**
 * @brief Allocate a new MSSMNoFV model.
 *
 * This function allocates a new MSSMNoFV model and returns a pointer
 * to the created object.  To prevent a resource leak, the model
 * should be destroyed using gm2calc_mssmnofv_free() .
 *
 * @return pointer to model object
 */
MSSMNoFV_onshell* gm2calc_mssmnofv_new()
{
   return reinterpret_cast<MSSMNoFV_onshell*>(new gm2calc::MSSMNoFV_onshell());
}

/**
 * @brief Deletes a MSSMNoFV model.
 *
 * This function deletes a MSSMNoFV model object, which has been
 * created using gm2calc_mssmnofv_new() .
 *
 * @param model pointer to model object
 */
void gm2calc_mssmnofv_free(MSSMNoFV_onshell* model)
{
   delete reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model);
}

void gm2calc_mssmnofv_set_alpha_MZ(MSSMNoFV_onshell* model, double alpha_MZ)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_alpha_MZ(alpha_MZ);
}

void gm2calc_mssmnofv_set_alpha_thompson(MSSMNoFV_onshell* model, double alpha_0)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_alpha_thompson(alpha_0);
}

void gm2calc_mssmnofv_set_Ae(MSSMNoFV_onshell* model, unsigned i, unsigned k, double a)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_Ae(i,k,a);
}

void gm2calc_mssmnofv_set_Au(MSSMNoFV_onshell* model, unsigned i, unsigned k, double a)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_Au(i,k,a);
}

void gm2calc_mssmnofv_set_Ad(MSSMNoFV_onshell* model, unsigned i, unsigned k, double a)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_Ad(i,k,a);
}

void gm2calc_mssmnofv_set_g3(MSSMNoFV_onshell* model, double g3)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_g3(g3);
}

void gm2calc_mssmnofv_set_MassB(MSSMNoFV_onshell* model, double mass_b)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_MassB(mass_b);
}

void gm2calc_mssmnofv_set_MassWB(MSSMNoFV_onshell* model, double mass_wb)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_MassWB(mass_wb);
}

void gm2calc_mssmnofv_set_MassG(MSSMNoFV_onshell* model, double mass_g)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_MassG(mass_g);
}

void gm2calc_mssmnofv_set_mq2(MSSMNoFV_onshell* model, unsigned i, unsigned k, double mq2)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_mq2(i,k,mq2);
}

void gm2calc_mssmnofv_set_mu2(MSSMNoFV_onshell* model, unsigned i, unsigned k, double mu2)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_mu2(i,k,mu2);
}

void gm2calc_mssmnofv_set_md2(MSSMNoFV_onshell* model, unsigned i, unsigned k, double md2)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_md2(i,k,md2);
}

void gm2calc_mssmnofv_set_ml2(MSSMNoFV_onshell* model, unsigned i, unsigned k, double ml2)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_ml2(i,k,ml2);
}

void gm2calc_mssmnofv_set_me2(MSSMNoFV_onshell* model, unsigned i, unsigned k, double me2)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_me2(i,k,me2);
}

void gm2calc_mssmnofv_set_Mu(MSSMNoFV_onshell* model, double mu)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_Mu(mu);
}

void gm2calc_mssmnofv_set_TB(MSSMNoFV_onshell* model, double tan_beta)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_TB(tan_beta);
}

void gm2calc_mssmnofv_set_scale(MSSMNoFV_onshell* model, double scale)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_scale(scale);
}

void gm2calc_mssmnofv_set_MAh_pole(MSSMNoFV_onshell* model, double MA0)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_MA0(MA0);
}

void gm2calc_mssmnofv_set_MZ_pole(MSSMNoFV_onshell* model, double MZ)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->get_physical().MVZ = MZ;
}

void gm2calc_mssmnofv_set_MW_pole(MSSMNoFV_onshell* model, double MW)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->get_physical().MVWm = MW;
}

void gm2calc_mssmnofv_set_MT_pole(MSSMNoFV_onshell* model, double MFt)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->get_physical().MFt = MFt;
}

void gm2calc_mssmnofv_set_MB_running(MSSMNoFV_onshell* model, double MFb)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->get_physical().MFb = MFb;
}

void gm2calc_mssmnofv_set_ML_pole(MSSMNoFV_onshell* model, double MFtau)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->get_physical().MFtau = MFtau;
}

void gm2calc_mssmnofv_set_MM_pole(MSSMNoFV_onshell* model, double MFm)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->get_physical().MFm = MFm;
}

void gm2calc_mssmnofv_set_MSm_pole(MSSMNoFV_onshell* model, unsigned i, double MSm)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->get_physical().MSm(i) = MSm;
}

void gm2calc_mssmnofv_set_MSvmL_pole(MSSMNoFV_onshell* model, double MSvmL)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->get_physical().MSvmL = MSvmL;
}

void gm2calc_mssmnofv_set_MCha_pole(MSSMNoFV_onshell* model, unsigned i, double MCha)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->get_physical().MCha(i) = MCha;
}

void gm2calc_mssmnofv_set_MChi_pole(MSSMNoFV_onshell* model, unsigned i, double MChi)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->get_physical().MChi(i) = MChi;
}

void gm2calc_mssmnofv_set_verbose_output(MSSMNoFV_onshell* model, int verbose_output)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_verbose_output(verbose_output != 0);
}

double gm2calc_mssmnofv_get_Ae(const MSSMNoFV_onshell* model, unsigned i, unsigned k)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_Ae(i,k);
}

double gm2calc_mssmnofv_get_Ad(const MSSMNoFV_onshell* model, unsigned i, unsigned k)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_Ad(i,k);
}

double gm2calc_mssmnofv_get_Au(const MSSMNoFV_onshell* model, unsigned i, unsigned k)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_Au(i,k);
}

double gm2calc_mssmnofv_get_EL(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_EL();
}

double gm2calc_mssmnofv_get_EL0(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_EL0();
}

double gm2calc_mssmnofv_get_gY(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_gY();
}

double gm2calc_mssmnofv_get_g1(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_g1();
}

double gm2calc_mssmnofv_get_g2(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_g2();
}

double gm2calc_mssmnofv_get_g3(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_g3();
}

double gm2calc_mssmnofv_get_TB(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_TB();
}

double gm2calc_mssmnofv_get_MassB(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_MassB();
}

double gm2calc_mssmnofv_get_MassWB(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_MassWB();
}

double gm2calc_mssmnofv_get_MassG(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_MassG();
}

double gm2calc_mssmnofv_get_Mu(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_Mu();
}

double gm2calc_mssmnofv_get_mq2(const MSSMNoFV_onshell* model, unsigned i, unsigned k)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_mq2(i,k);
}

double gm2calc_mssmnofv_get_md2(const MSSMNoFV_onshell* model, unsigned i, unsigned k)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_md2(i,k);
}

double gm2calc_mssmnofv_get_mu2(const MSSMNoFV_onshell* model, unsigned i, unsigned k)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_mu2(i,k);
}

double gm2calc_mssmnofv_get_ml2(const MSSMNoFV_onshell* model, unsigned i, unsigned k)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_ml2(i,k);
}

double gm2calc_mssmnofv_get_me2(const MSSMNoFV_onshell* model, unsigned i, unsigned k)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_me2(i,k);
}

double gm2calc_mssmnofv_get_scale(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_scale();
}

double gm2calc_mssmnofv_get_vev(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_vev();
}

double gm2calc_mssmnofv_get_MW(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_MW();
}

double gm2calc_mssmnofv_get_MZ(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_MZ();
}

double gm2calc_mssmnofv_get_ME(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_ME();
}

double gm2calc_mssmnofv_get_MM(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_MM();
}

double gm2calc_mssmnofv_get_ML(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_ML();
}

double gm2calc_mssmnofv_get_MU(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_MU();
}

double gm2calc_mssmnofv_get_MC(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_MC();
}

double gm2calc_mssmnofv_get_MT(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_MT();
}

double gm2calc_mssmnofv_get_MD(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_MD();
}

double gm2calc_mssmnofv_get_MS(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_MS();
}

double gm2calc_mssmnofv_get_MB(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_MB();
}

double gm2calc_mssmnofv_get_MBMB(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_MBMB();
}

double gm2calc_mssmnofv_get_MAh(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_MAh(1);
}

double gm2calc_mssmnofv_get_Mhh(const MSSMNoFV_onshell* model, unsigned i)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_Mhh(i);
}

double gm2calc_mssmnofv_get_MCha(const MSSMNoFV_onshell* model, unsigned i)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_MCha(i);
}

double gm2calc_mssmnofv_get_UM(const MSSMNoFV_onshell* model, unsigned i, unsigned j, double* u_imag)
{
   const std::complex<double> U_ij =
      reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_UM(i,j);

   if (u_imag != nullptr) {
      *u_imag = std::imag(U_ij);
   }

   return std::real(U_ij);
}

double gm2calc_mssmnofv_get_UP(const MSSMNoFV_onshell* model, unsigned i, unsigned j, double* u_imag)
{
   const std::complex<double> U_ij =
      reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_UP(i,j);

   if (u_imag != nullptr) {
      *u_imag = std::imag(U_ij);
   }

   return std::real(U_ij);
}

double gm2calc_mssmnofv_get_MChi(const MSSMNoFV_onshell* model, unsigned i)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_MChi(i);
}

double gm2calc_mssmnofv_get_ZN(const MSSMNoFV_onshell* model, unsigned i, unsigned j, double* u_imag)
{
   const std::complex<double> U_ij =
      reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_ZN(i,j);

   if (u_imag != nullptr) {
      *u_imag = std::imag(U_ij);
   }

   return std::real(U_ij);
}

double gm2calc_mssmnofv_get_MSe(const MSSMNoFV_onshell* model, unsigned i)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_MSe(i);
}

double gm2calc_mssmnofv_get_MSveL(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_MSveL();
}

double gm2calc_mssmnofv_get_MSm(const MSSMNoFV_onshell* model, unsigned i)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_MSm(i);
}

double gm2calc_mssmnofv_get_MSvmL(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_MSvmL();
}

double gm2calc_mssmnofv_get_MStau(const MSSMNoFV_onshell* model, unsigned i)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_MStau(i);
}

double gm2calc_mssmnofv_get_MSvtL(const MSSMNoFV_onshell* model)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_MSvtL();
}

double gm2calc_mssmnofv_get_MSu(const MSSMNoFV_onshell* model, unsigned i)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_MSu(i);
}

double gm2calc_mssmnofv_get_MSd(const MSSMNoFV_onshell* model, unsigned i)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_MSd(i);
}

double gm2calc_mssmnofv_get_MSc(const MSSMNoFV_onshell* model, unsigned i)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_MSc(i);
}

double gm2calc_mssmnofv_get_MSs(const MSSMNoFV_onshell* model, unsigned i)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_MSs(i);
}

double gm2calc_mssmnofv_get_MSt(const MSSMNoFV_onshell* model, unsigned i)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_MSt(i);
}

double gm2calc_mssmnofv_get_MSb(const MSSMNoFV_onshell* model, unsigned i)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_MSb(i);
}

double gm2calc_mssmnofv_get_USe(const MSSMNoFV_onshell* model, unsigned i, unsigned j)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_ZE(i,j);
}

double gm2calc_mssmnofv_get_USm(const MSSMNoFV_onshell* model, unsigned i, unsigned j)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_ZM(i,j);
}

double gm2calc_mssmnofv_get_UStau(const MSSMNoFV_onshell* model, unsigned i, unsigned j)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_ZTau(i,j);
}

double gm2calc_mssmnofv_get_USu(const MSSMNoFV_onshell* model, unsigned i, unsigned j)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_ZU(i,j);
}

double gm2calc_mssmnofv_get_USd(const MSSMNoFV_onshell* model, unsigned i, unsigned j)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_ZD(i,j);
}

double gm2calc_mssmnofv_get_USc(const MSSMNoFV_onshell* model, unsigned i, unsigned j)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_ZC(i,j);
}

double gm2calc_mssmnofv_get_USs(const MSSMNoFV_onshell* model, unsigned i, unsigned j)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_ZS(i,j);
}

double gm2calc_mssmnofv_get_USt(const MSSMNoFV_onshell* model, unsigned i, unsigned j)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_ZT(i,j);
}

double gm2calc_mssmnofv_get_USb(const MSSMNoFV_onshell* model, unsigned i, unsigned j)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_ZB(i,j);
}

double gm2calc_mssmnofv_get_Ye(const MSSMNoFV_onshell* model, unsigned i, unsigned k)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_Ye(i,k);
}

double gm2calc_mssmnofv_get_Yd(const MSSMNoFV_onshell* model, unsigned i, unsigned k)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_Yd(i,k);
}

double gm2calc_mssmnofv_get_Yu(const MSSMNoFV_onshell* model, unsigned i, unsigned k)
{
   return reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model)->get_Yu(i,k);
}

/**
 * This function converts the model parameters to a mixed
 * on-shell/DR-bar scheme, used to calculate \f$a_\mu\f$.  The
 * function uses default values for the conversion precision goal and
 * the maximum number of iterations.
 *
 * @param model pointer to model object
 *
 * @return error code gm2calc_error
 */
gm2calc_error gm2calc_mssmnofv_convert_to_onshell(MSSMNoFV_onshell* model)
{
   gm2calc_error error = gm2calc_NoError;

   try {
      reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->convert_to_onshell();
   } catch (const gm2calc::EInvalidInput&) {
      error = gm2calc_InvalidInput;
   } catch (const gm2calc::EPhysicalProblem&) {
      error = gm2calc_PhysicalProblem;
   } catch (...) {
      error = gm2calc_UnknownError;
   }

   return error;
}

/**
 * This function converts the model parameters to a mixed
 * on-shell/DR-bar scheme, used to calculate \f$a_\mu\f$.
 *
 * @param model pointer to model object
 * @param precision precision goal of the conversion
 * @param max_iterations maximum number of iterations
 *
 * @return error code gm2calc_error
 */
gm2calc_error gm2calc_mssmnofv_convert_to_onshell_params(
   MSSMNoFV_onshell* model, double precision, unsigned max_iterations)
{
   gm2calc_error error = gm2calc_NoError;

   try {
      reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->convert_to_onshell(precision, max_iterations);
   } catch (const gm2calc::EInvalidInput&) {
      error = gm2calc_InvalidInput;
   } catch (const gm2calc::EPhysicalProblem&) {
      error = gm2calc_PhysicalProblem;
   } catch (...) {
      error = gm2calc_UnknownError;
   }

   return error;
}

/**
 * This function calculates the masses of the particles in the model.
 *
 * @param model pointer to model object
 *
 * @return error code gm2calc_error
 */
gm2calc_error gm2calc_mssmnofv_calculate_masses(MSSMNoFV_onshell* model)
{
   gm2calc_error error = gm2calc_NoError;

   try {
      reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->calculate_masses();
   } catch (const gm2calc::EInvalidInput&) {
      error = gm2calc_InvalidInput;
   } catch (const gm2calc::EPhysicalProblem&) {
      error = gm2calc_PhysicalProblem;
   } catch (...) {
      error = gm2calc_UnknownError;
   }

   return error;
}

/**
 * Returns true if there are problems
 *
 * @param model pointer to model object
 *
 * @return true if there are problems, false otherwise
 */
int gm2calc_mssmnofv_have_problem(MSSMNoFV_onshell* model)
{
   return static_cast<int>(reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->get_problems().have_problem());
}

/**
 * Returns true if there are warnings
 *
 * @param model pointer to model object
 *
 * @return true if there are warnings, false otherwise
 */
int gm2calc_mssmnofv_have_warning(MSSMNoFV_onshell* model)
{
   return static_cast<int>(reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->get_problems().have_warning());
}

/**
 * Fills string with problem descriptions
 *
 * @param model pointer to model object
 * @param msg buffer for message string
 * @param len available length of message string
 *
 */
void gm2calc_mssmnofv_get_problems(MSSMNoFV_onshell* model, char* msg, unsigned len)
{
   if (msg == nullptr) {
      return;
   }

   const std::string str(reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->get_problems().get_problems());
   msg[str.copy(msg, len - 1)] = '\0';
}

/**
 * Fills string with warning descriptions
 *
 * @param model pointer to model object
 * @param msg buffer for message string
 * @param len available length of message string
 *
 */
void gm2calc_mssmnofv_get_warnings(MSSMNoFV_onshell* model, char* msg, unsigned len)
{
   if (msg == nullptr) {
      return;
   }

   const std::string str(reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->get_problems().get_warnings());
   msg[str.copy(msg, len - 1)] = '\0';
}

void print_mssmnofv(const MSSMNoFV_onshell* model)
{
   VERBOSE(*reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model));
}

} // extern "C"
