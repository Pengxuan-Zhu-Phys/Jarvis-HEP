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

#include "THDM/gm2_2loop_helpers.hpp"
#include "gm2_dilog.hpp"
#include "gm2_ffunctions.hpp"

#include <cmath>
#include <limits>

/**
 * \file gm2_2loop_F.cpp
 *
 * Contains functions necessary to calculate the fermionic THDM
 * contributions for g-2 at the 2-loop level.
 */

namespace gm2calc {

namespace thdm {

namespace {

struct F_sm_pars {
   double mw2{}; ///< squared W boson mass
   double mz2{}; ///< squared Z boson mass
};

struct F_neut_pars {
   double qf{};  ///< electromagnetic charge of fermion f
   double ql{};  ///< electromagnetic charge of fermion l
   double t3f{}; ///< SU(2)_L charge of ferimon f
   double t3l{}; ///< SU(2)_L charge of ferimon l
   double nc{};  ///< number of colors of fermion f
};

struct F_char_pars {
   double qd{};  ///< electromagnetic charge of up-type fermion
   double qu{};  ///< electromagnetic charge of down-type fermion
   double nc{};  ///< number of colors
};

const double pi = 3.1415926535897932;
const double pi2 = 9.8696044010893586; // Pi^2
const double q_u =  2.0/3.0;  ///< electric charge of up-type quarks
const double q_d = -1.0/3.0;  ///< electric charge of down-type quarks
const double q_v =  0.0;      ///< electric charge of up-type leptons
const double q_l = -1.0;      ///< electric charge of down-type leptons
const double t3_u = +0.5;     ///< SU(2)_L charge of up-type quark
const double t3_d = -0.5;     ///< SU(2)_L charge of down-type quark
const double t3_l = -0.5;     ///< SU(2)_L charge of charged lepton

double sqr(double x) noexcept { return x*x; }

double calc_v2(const THDM_F_parameters& thdm) noexcept
{
   const double mw2 = sqr(thdm.mw);
   const double mz2 = sqr(thdm.mz);
   const double sw2 = 1.0 - mw2/mz2;
   const double e2 = 4*pi*thdm.alpha_em;
   const double g22 = e2/sw2;
   return 4*mw2/g22;
}

/// Eq (56), arxiv:1607.06292, S = h or H
double FH(double ms2, double mf2) noexcept
{
   const double x = mf2/ms2;
   return 0.5/x*f_S(x);
}

/// Eq (57), arxiv:1607.06292, S = A
double FA(double ms2, double mf2) noexcept
{
   const double x = mf2/ms2;
   return 0.5/x*f_PS(x);
}

double FHZ(double mf2, double ms2, double mz2) noexcept
{
   const double x = mf2/ms2;
   const double y = mf2/mz2;
   return FSZ(x, y);
}

double FAZ(double mf2, double ms2, double mz2) noexcept
{
   const double x = mf2/ms2;
   const double y = mf2/mz2;
   return FPZ(x, y);
}

/// Eq (54), arxiv:1607.06292, S = h or H
template <typename F>
double fSgamma(double ms2, double mf2, const F_neut_pars& pars, F FH) noexcept
{
   const double qf2 = sqr(pars.qf);
   const double nc = pars.nc;

   return qf2*nc * mf2/ms2 * FH(ms2, mf2);
}

/// Eq (55), arxiv:1607.06292, S = h or H
template <typename F>
double fSZ(double ms2, double mf2, const F_neut_pars& pars, const F_sm_pars& sm, F FHZ) noexcept
{
   const double mw2 = sm.mw2;
   const double mz2 = sm.mz2;
   const double cw2 = mw2/mz2;
   const double sw2 = 1.0 - cw2;
   const double qf = pars.qf;
   const double ql = pars.ql;
   const double nc = pars.nc;
   const double gvf = 0.5*pars.t3f - qf*sw2;
   const double gvl = 0.5*pars.t3l - ql*sw2;

   // Note: 0.5*FHZ(mf2, ms2, mz2) = -mf2/(ms2 - mz2) * (FH(ms2, mf2) - FH(mz2, mf2))
   // with FH = {f_S, f_PS}

   return 0.5*nc*qf*gvl*gvf/(sw2*cw2)*FHZ(mf2, ms2, mz2);
}

/// Eq (53), arxiv:1607.06292, S = h or H
template <typename F, typename FZ>
double ffS(double ms2, double mf2, const F_neut_pars& pars, const F_sm_pars& sm, F FH, FZ FHZ) noexcept
{
   return fSgamma(ms2, mf2, pars, FH) + fSZ(ms2, mf2, pars, sm, FHZ);
}

/// Eq (59), arxiv:1607.06292, S = H^\pm, f = l
double flHp(double ms2, double ml2, const F_char_pars& pars, const F_sm_pars& sm) noexcept
{
   const double mw2 = sm.mw2;
   const double nc = pars.nc;
   const double x = ml2/ms2;
   const double y = ml2/mw2;

   return -nc*FCWl(x, y);
}

/// Eq (59), arxiv:1607.06292, S = H^\pm, f = u
double fuHp(double ms2, double md2, double mu2, const F_char_pars& pars, const F_sm_pars& sm) noexcept
{
   const double mw2 = sm.mw2;
   const double qd = pars.qd;
   const double qu = pars.qu;
   const double nc = pars.nc;
   const double xu = mu2/ms2;
   const double xd = md2/ms2;
   const double yu = mu2/mw2;
   const double yd = md2/mw2;

   return -nc*FCWu(xu, xd, yu, yd, qu, qd);
}

/// Eq (59), arxiv:1607.06292, S = H^\pm, f = d
double fdHp(double ms2, double md2, double mu2, const F_char_pars& pars, const F_sm_pars& sm) noexcept
{
   const double mw2 = sm.mw2;
   const double qd = pars.qd;
   const double qu = pars.qu;
   const double nc = pars.nc;
   const double xu = mu2/ms2;
   const double xd = md2/ms2;
   const double yu = mu2/mw2;
   const double yd = md2/mw2;

   return -nc*FCWd(xu, xd, yu, yd, qu, qd);
}

} // anonymous namespace

/// Eq (53), arxiv:1607.06292, f = u, S = h or H
double fuS(double ms2, double mu2, double mw2, double mz2) noexcept
{
   const F_neut_pars pars{q_u, q_l, t3_u, t3_l, 3.0};
   const F_sm_pars sm{ mw2, mz2 };
   const auto lFH = [] (double ms2, double mu2) { return FH(ms2, mu2); };
   const auto lFHZ = [] (double mf2, double ms2, double mz2) { return FHZ(mf2, ms2, mz2); };

   return ffS(ms2, mu2, pars, sm, lFH, lFHZ);
}

/// Eq (53), arxiv:1607.06292, f = d, S = h or H
double fdS(double ms2, double md2, double mw2, double mz2) noexcept
{
   const F_neut_pars pars{q_d, q_l, t3_d, t3_l, 3.0};
   const F_sm_pars sm{ mw2, mz2 };
   const auto lFH = [] (double ms2, double md2) { return FH(ms2, md2); };
   const auto lFHZ = [] (double mf2, double ms2, double mz2) { return FHZ(mf2, ms2, mz2); };

   return ffS(ms2, md2, pars, sm, lFH, lFHZ);
}

/// Eq (53), arxiv:1607.06292, f = l, S = h or H
double flS(double ms2, double ml2, double mw2, double mz2) noexcept
{
   const F_neut_pars pars{q_l, q_l, t3_l, t3_l, 1.0};
   const F_sm_pars sm{ mw2, mz2 };
   const auto lFH = [] (double ms2, double ml2) { return FH(ms2, ml2); };
   const auto lFHZ = [] (double mf2, double ms2, double mz2) { return FHZ(mf2, ms2, mz2); };

   return ffS(ms2, ml2, pars, sm, lFH, lFHZ);
}

/// Eq (53), arxiv:1607.06292, f = u, S = A
double fuA(double ms2, double mu2, double mw2, double mz2) noexcept
{
   const F_neut_pars pars{q_u, q_l, t3_u, t3_l, 3.0};
   const F_sm_pars sm{ mw2, mz2 };
   const auto lFA = [] (double ms2, double mu2) { return FA(ms2, mu2); };
   const auto lFAZ = [] (double mf2, double ms2, double mz2) { return FAZ(mf2, ms2, mz2); };

   return ffS(ms2, mu2, pars, sm, lFA, lFAZ);
}

/// Eq (53), arxiv:1607.06292, f = d, S = A
double fdA(double ms2, double md2, double mw2, double mz2) noexcept
{
   const F_neut_pars pars{q_d, q_l, t3_d, t3_l, 3.0};
   const F_sm_pars sm{ mw2, mz2 };
   const auto lFA = [] (double ms2, double md2) { return FA(ms2, md2); };
   const auto lFAZ = [] (double mf2, double ms2, double mz2) { return FAZ(mf2, ms2, mz2); };

   return ffS(ms2, md2, pars, sm, lFA, lFAZ);
}

/// Eq (53), arxiv:1607.06292, f = l, S = A
double flA(double ms2, double ml2, double mw2, double mz2) noexcept
{
   const F_neut_pars pars{q_l, q_l, t3_l, t3_l, 1.0};
   const F_sm_pars sm{mw2, mz2};
   const auto lFA = [] (double ms2, double ml2) { return FA(ms2, ml2); };
   const auto lFAZ = [] (double mf2, double ms2, double mz2) { return FAZ(mf2, ms2, mz2); };

   return ffS(ms2, ml2, pars, sm, lFA, lFAZ);
}

/// Eq (59), arxiv:1607.06292, S = H^\pm, f = u
double fuHp(double ms2, double md2, double mu2, double mw2, double mz2) noexcept
{
   const F_char_pars pars{q_d, q_u, 3.0};
   const F_sm_pars sm{mw2, mz2};
   return fuHp(ms2, md2, mu2, pars, sm);
}

/// Eq (59), arxiv:1607.06292, S = H^\pm, f = d
double fdHp(double ms2, double md2, double mu2, double mw2, double mz2) noexcept
{
   const F_char_pars pars{q_d, q_u, 3.0};
   const F_sm_pars sm{mw2, mz2};
   return fdHp(ms2, md2, mu2, pars, sm);
}

/// Eq (59), arxiv:1607.06292, S = H^\pm, f = l
double flHp(double ms2, double ml2, double mw2, double mz2) noexcept
{
   const F_char_pars pars{q_l, q_v, 1.0};
   const F_sm_pars sm{mw2, mz2};
   return flHp(ms2, ml2, pars, sm);
}

/**
 * Calculates 2-loop fermionic contributions with charged Higgs
 * boson.
 *
 * Eq (52), arXiv:2110.13238
 */
double amu2L_F_charged(const THDM_F_parameters& thdm) noexcept
{
   const double mHp2 = sqr(thdm.mHp);
   const double mw2 = sqr(thdm.mw);
   const double mz2 = sqr(thdm.mz);
   const double v2 = calc_v2(thdm);
   const double sw2 = 1 - mw2/mz2;
   const double pref = sqr(thdm.alpha_em*thdm.mm/(8*pi*thdm.mw*sw2))*v2/thdm.ml(1);

   double res = 0;

   // loop over generations
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         // u
         res += fuHp(mHp2, sqr(thdm.md(j)), sqr(thdm.mu(i)), mw2, mz2)
                *std::real(std::conj(thdm.yuHp(i,j))*thdm.vckm(i,j)*thdm.ylHp(1,1))/thdm.mu(i);
         // d
         res += fdHp(mHp2, sqr(thdm.md(j)), sqr(thdm.mu(i)), mw2, mz2)
                *std::real(std::conj(thdm.ydHp(i,j))*thdm.vckm(i,j)*thdm.ylHp(1,1))/thdm.md(j);
      }
      // l
      res += flHp(mHp2, sqr(thdm.ml(i)), mw2, mz2)
             *std::real(std::conj(thdm.ylHp(i,i))*thdm.ylHp(1,1))/thdm.ml(i);
   }

   return pref*res;
}

/**
 * Calculates 2-loop fermionic contributions with neutral Higgs
 * bosons.
 *
 * Eq (53), arxiv:1607:06292
 */
double amu2L_F_neutral(const THDM_F_parameters& thdm) noexcept
{
   const double mh2 = sqr(thdm.mh(0));
   const double mH2 = sqr(thdm.mh(1));
   const double mA2 = sqr(thdm.mA);
   const double mhSM2 = sqr(thdm.mhSM);
   const double mw2 = sqr(thdm.mw);
   const double mz2 = sqr(thdm.mz);
   const double v2 = calc_v2(thdm);
   const double sw2 = 1 - mw2/mz2;
   const auto pref = sqr(thdm.alpha_em*thdm.mm/(2*pi*thdm.mw))/sw2;

   double res = 0;

   // loop over generations
   for (int i = 0; i < 3; ++i) {
      // h
      res += fuS(mh2, sqr(thdm.mu(i)), mw2, mz2)*std::real(std::conj(thdm.yuh(i,i))*thdm.ylh(1,1))*v2/(thdm.mu(i)*thdm.ml(1));
      res += fdS(mh2, sqr(thdm.md(i)), mw2, mz2)*std::real(std::conj(thdm.ydh(i,i))*thdm.ylh(1,1))*v2/(thdm.md(i)*thdm.ml(1));
      res += flS(mh2, sqr(thdm.ml(i)), mw2, mz2)*std::real(std::conj(thdm.ylh(i,i))*thdm.ylh(1,1))*v2/(thdm.ml(i)*thdm.ml(1));

      // H
      res += fuS(mH2, sqr(thdm.mu(i)), mw2, mz2)*std::real(std::conj(thdm.yuH(i,i))*thdm.ylH(1,1))*v2/(thdm.mu(i)*thdm.ml(1));
      res += fdS(mH2, sqr(thdm.md(i)), mw2, mz2)*std::real(std::conj(thdm.ydH(i,i))*thdm.ylH(1,1))*v2/(thdm.md(i)*thdm.ml(1));
      res += flS(mH2, sqr(thdm.ml(i)), mw2, mz2)*std::real(std::conj(thdm.ylH(i,i))*thdm.ylH(1,1))*v2/(thdm.ml(i)*thdm.ml(1));

      // A
      res += fuA(mA2, sqr(thdm.mu(i)), mw2, mz2)*std::real(std::conj(thdm.yuA(i,i))*thdm.ylA(1,1))*v2/(thdm.mu(i)*thdm.ml(1));
      res += fdA(mA2, sqr(thdm.md(i)), mw2, mz2)*std::real(std::conj(thdm.ydA(i,i))*thdm.ylA(1,1))*v2/(thdm.md(i)*thdm.ml(1));
      res += flA(mA2, sqr(thdm.ml(i)), mw2, mz2)*std::real(std::conj(thdm.ylA(i,i))*thdm.ylA(1,1))*v2/(thdm.ml(i)*thdm.ml(1));

      // subtract hSM
      res -= fuS(mhSM2, sqr(thdm.mu(i)), mw2, mz2);
      res -= fdS(mhSM2, sqr(thdm.md(i)), mw2, mz2);
      res -= flS(mhSM2, sqr(thdm.ml(i)), mw2, mz2);
   }

   return pref*res;
}

/**
 * Calculates the sum of the 2-loop fermionic contributions with
 * neutral and charged Higgs bosons.
 *
 * Eq (63), arxiv:1607:06292
 */
double amu2L_F(const THDM_F_parameters& thdm) noexcept
{
   return amu2L_F_neutral(thdm) + amu2L_F_charged(thdm);
}

} // namespace thdm

} // namespace gm2calc
