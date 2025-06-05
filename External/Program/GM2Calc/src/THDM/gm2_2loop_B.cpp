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
#include "gm2_numerics.hpp"

/**
 * \file gm2_2loop_B.cpp
 *
 * Contains functions necessary to calculate the bosonic THDM
 * contributions for g-2 at the 2-loop level.
 */

namespace gm2calc {

namespace thdm {

namespace {

const double pi = 3.1415926535897932;
const double pi2 = 9.8696044010893586; // Pi^2

const double eps_shift = 1e-8; // parameter shift to avoid spurious divergences

/// shift value away from limit, if it is close to the limit
void shift(double& val, double limit, double eps) noexcept
{
   if (is_equal_rel(val, limit, eps)) {
      val = (1 + eps)*limit;
   }
}

/// Eq.(102), arxiv:1607.06292
double YF1(double u, double w, double cw2) noexcept
{
   shift(u, 1.0, eps_shift);

   const auto cw4 = cw2*cw2;
   const auto c0 = cw2*(-1 + cw2)*(u + 2*w)/u;

   // Note: Phi(w,w,1) == 0.5/w*f_PS(w)*(1 - 4*w)
   // Note: Phi(u,w,w) == 0.5*u/w*f_PS(w/u)*(u - 4*w)

   return
      - 72*c0 - 36*c0*std::log(w)
      + 9*(-8*cw4 - 3*u + 2*cw2*(4 + u))*(u + 2*w)/(2*(u - 1)*u)*std::log(u)
      + 9./2*(3 - 10*cw2 + 8*cw4)*(u + 2*w)/(u - 1)*f_PS(w)
      - 9./2*(8*cw4 + 3*u - 2*cw2*(4 + u))*(u + 2*w)/((u - 1)*u)*f_PS(w/u)
      ;
}

/// Eq.(121), arxiv:1607.06292
double YFZ(double u, double cw2) noexcept
{
   const auto cw4 = cw2*cw2;
   const auto u2 = u*u;
   const auto lu = std::log(u);
   const auto li = dilog(1.0 - u);
   const auto phi = 0.5*u*f_PS(1/u)*(u - 4); // Phi(u,1,1);

   const auto z1 = 3*(17 - 48*cw2 + 32*cw4); // Eq.(122)
   const auto z2 = 5 - 12*cw2 + 8*cw4;       // Eq.(123)
   const auto z3 = 3*(1 - 3*cw2 + 2*cw4);    // Eq.(124)

   const double res =
      + z1*u*li
      + z2/(2*u2)*(6*(-4 + u)*u + pi2*(4 + 3*u) + 6*u*(4 + u)*lu
                   - 6*(4 + 3*u)*li + 6*u*(2 + u)*phi)
      + z3*u*(6 + pi2*(-4 + u)*u + 3*lu*(4 + (-4 + u)*u*lu)
              + 12*(-4 + u)*u*li + 6*(-2 + u)*phi);

   return res;
}

/// Eq.(125), arxiv:1607.06292
double YFW(double u, double cw2) noexcept
{
   const auto cw4 = cw2*cw2;
   const auto cw6 = cw4*cw2;
   const auto u2 = u*u;
   const auto u3 = u2*u;

   // Note: Phi(u,cw2,cw2) == 0.5*u/cw2*f_PS(cw2/u)*(u - 4*cw2)

   const double res =
      - 57.0/2*cw2 - 4*cw6*pi2/u2 + 3./4*cw4*(32 - 3*pi2)/u
      + 3./2*(16*cw6 + 9*cw4*u + 12*cw2*u2 - 19*u3)/u2*dilog(1.0 - u/cw2)
      + 3./2*cw2*(16*cw2 + 19*u)/u*std::log(cw2/u)
      - 3./4*(4*cw4 - 50*cw2*u + 19*u2)/cw2*f_PS(cw2/u);

   return res;
}

/// Eq.(105), arxiv:1607.06292
double YF2(double u, double cw2) noexcept
{
   shift(u, 4*cw2, eps_shift);
   shift(u, 1.0, eps_shift);

   const auto cw4 = cw2*cw2;
   const auto cw6 = cw4*cw2;
   const auto cw8 = cw4*cw4;
   const auto u2 = u*u;
   const auto u3 = u2*u;

   const double f0 = 3.0/4*cw4*(-640 + 576*cw2 + 7*pi2);     // Eq.(106)
   const double f1 = 96*cw6*(11 - 53*cw2 + 36*cw4);          // Eq.(107)
   const double f2 = -3.0/4*cw2*(-66*cw2 - 48*cw4 + 672*cw6);// Eq.(108)
   const double f3 = -3.0/4*cw2*(109 - 430*cw2 + 120*cw4);   // Eq.(109)
   const double f4 = 96*cw6*(-11 + 9*cw2);                   // Eq.(110)
   const double f5 = 45.0/2*cw4 + 192*cw6;                   // Eq.(111)
   const double f6 = 3.0/4*cw2*(157 + 90*cw2);               // Eq.(112)
   const double f7 = -3.0/4*(18 + 61*cw2);                   // Eq.(113)
   const double f8 = -7 + 61*cw2 - 162*cw4 + 96*cw6;         // Eq.(114)
   const double f9 = 1 - 5*cw2 + 10*cw4;                     // Eq.(115)
   const double f10 = -1728*cw8*(-1 + cw2);                  // Eq.(116)
   const double f11 = 3*cw6*(-899 + 768*cw2);                // Eq.(117)
   const double f12 = 387*cw4 - 363*cw6;                     // Eq.(118)
   const double f13 = 9.0/2*cw2*(57 + 106*cw2);              // Eq.(119)
   const double f14 = -15.0/2*(7 + 45*cw2);                  // Eq.(120)

   // Note: Phi(cw2,cw2,1) == 0.5/cw2*f_PS(cw2)*(1 - 4*cw2)
   // Note: Phi(u,cw2,cw2) == 0.5*u/cw2*f_PS(cw2/u)*(u - 4*cw2)

   const double res =
      + 8*cw6*pi2/u2 + f0/u + 393.0/8*cw2
      + (f1/u + f2 + f3*u)*std::log(cw2)/((4*cw2-1)*(4*cw2-u))
      + (f4/u + f5 + f6*u + f7*u2)*std::log(u)/((u-1)*(4*cw2-u))
      - 3.0/2*(32*cw6/u2 + 21*cw4/u + 15*cw2 - 35*u)*dilog(1.0 - u/cw2)
      + (f8 + f9*u)*9./4*(-3 + 4*cw2)*f_PS(cw2)/((1-4*cw2)*(u-1))
      + 0.5*(f10/u + f11 + f12*u + f13*u2 + f14*u3 + 105.0/2*u2*u2)*f_PS(cw2/u)
        /(cw2*(u-4*cw2)*(u-1))
      ;

   return YFW(u, cw2) + YFZ(u, cw2) + res;
}

/// Eq.(126), arxiv:1607.06292
double YF3(double u, double w, double cw2) noexcept
{
   shift(u, 4*cw2, eps_shift);
   shift(w, cw2, eps_shift);

   const auto cw4 = cw2*cw2;
   const auto cw6 = cw4*cw2;
   const auto cw8 = cw4*cw4;
   const auto u2 = u*u;
   const auto u3 = u2*u;
   const auto u4 = u2*u2;
   const auto u5 = u4*u;
   const auto w2 = w*w;
   const auto lc = std::log(cw2);
   const auto lu = std::log(u);
   const auto lw = std::log(w);

   // Eq.(127)
   const auto a1 = -9*cw2*u3 + 9*cw2*u2*(3*cw2+w) + 27*cw4*u*(w-cw2)
      + 9*(cw8 - 4*cw6*w + 3*cw4*w2);
   // Eq.(128)
   const auto a2 = 9*cw4*w/2 - 9*u2*(5*cw2 + w) + u*(36*cw4 + 153*cw2*w/4)
      + 9*u3;
   // Eq.(129)
   const auto a3 = 9*cw2*u2 - 9.0/2*cw2*u*(4*cw2 + w);
   // Eq.(130)
   const auto a4 = -9.0/2*u2*w*(2*cw4 + 9*cw2*w + 2*w2)
      + 9.0/8*u*w*(32*cw6 + 13*cw4*w + 35*cw2*w2) + 9*u3*w2;
   // Eq.(131)
   const auto a5 = -9*u3*(cw2 + w) - 9*u*(3*cw6 + 2*cw2*w2)
      + 9*u2*(3*cw4 + 4*cw2*w + w2) + 9.0/2*cw4*(2*cw4 - 6*cw2*w + w2);
   // Eq.(132)
   const auto a6 = -9*u4*(9*cw2 + w) + u*(81*cw6*w - 225*cw8) + 9*cw8*(w - cw2)
      - 9.0/2*u2*(3*cw6 + 37*cw4*w) + u3*(198*cw4 + 72*cw2*w) + 9*u5;
   // Eq.(133)
   const auto a7 = -9*cw2*u4 + 18*cw2*u3*(2*cw2 + w) + 36*u*(cw8 - 2*cw6*w)
      - 9*cw2*u2*(6*cw4 - cw2*w + w2) - 9*cw2*(cw2 - 3*w)*(cw6 - 2*cw4*w + cw2*w2);

   // Note: Phi(u,cw2,cw2) == 0.5*u/cw2*f_PS(cw2/u)*(u - 4*cw2)

   const double res =
      + 9*u*(2*cw2 - u + w)/w
      + (a1*(lu - lc) + 9*cw4*(cw4 - 4*cw2*w + 3*w2)*lc)*(lw - lc)/(2*w2*(cw2-w))
      + a2*lu/(w*(4*cw2-u))
      + a3*lw/(w*(cw2-w))
      + a4*lc/(w2*(4*cw2-u)*(cw2-w))
      + a5/(cw2*w2)*dilog(1.0 - u/cw2)
      + a6/(sqr(cw2)*(u-4*cw2)*(cw2-w))*0.5*f_PS(cw2/u)
      + a7/(w2*(cw2-w)*(cw4-2*cw2*(u+w)+sqr(u-w)))*Phi(u,w,cw2)
      ;

   return res;
}

/// Eq.(72), arxiv:1607.06292
double T0(double u, double w, double cw2) noexcept
{
   const auto cw4 = cw2*cw2;

   return
      9/cw4*(u - w)*(cw2*(u - w)*(u + 2*w) - cube(u - w) + cw4*w)
      /(cw4 + sqr(u - w) - 2*cw2*(u + w))*Phi(u,w,cw2);
}

/// Eq.(73), arxiv:1607.06292
double T1(double u, double w, double cw2) noexcept
{
   const auto cw4 = cw2*cw2;

   return 9/cw4*(u - w)*(cw2*w - sqr(u - w))*dilog(1.0 - u/w);
}

/// calculates potentially divergent (a^2 Log[a] - b^2 Log[b])/(a - b)
double dxlog(double a, double b) noexcept
{
   const double lb = std::log(b);

   if (is_equal_rel(a, b, 1e-4)) {
      return b*(1 + 2*lb) + (a - b)*(1.5 + lb) + sqr(a - b)/(3*b);
   }

   return (sqr(a)*std::log(a) - sqr(b)*lb)/(a - b);
}

/**
 * Calculates the following combination of T2 loop functions, Eq.(74),
 * arxiv:1607.06292:
 *
 * (xA - xH)/(xA - xHp)*T2p[xA, xH] + T2m[xH, xHp] + T2p[xHp, xH] + T2p[xHp, xA]
 */
double TX(double xH, double xA, double xHp, double cw2) noexcept
{
   const auto cw4 = sqr(cw2);
   const auto sw2 = 1.0 - cw2;
   const auto f6 = (7 - 14*cw2 + 4*cw4)/(4*cw2*sw2);
   const auto f7 = 1 - 6*cw2 + 4*cw4;
   const auto f8 = (13 - 20*cw2 + 4*cw4)/(cw2*sw2);
   const auto f9 = 7 - 12*cw2 + 8*cw4;
   const auto lH = std::log(xH);
   const auto lA = std::log(xA);
   const auto xAH  = dxlog(xA, xH);
   const auto xAHp = dxlog(xA, xHp);
   const auto xHHp = dxlog(xH, xHp);

   return
      (cw2 + 2*cw4 - 3*f9*xA)*lA/2
      + (3*cw2*f6 - (3*f8*xA)/2)*(xA - xHp)*lA
      + 3*f6*sqr(xA - xHp)*lA
      + (f6*cube(xA - xHp)*lA)/cw2
      + (cw2 + 2*cw4 - 3*f9*xH)*lH/2
      + (3*cw2*f6 - (3*f8*xH)/2)*(xH - xHp)*lH
      + 3*f6*sqr(xH - xHp)*lH
      + (f6*cube(xH - xHp)*lH)/cw2
      + 3*(f7*xAH + xAHp + xHHp);
}

/// Eq.(75), arxiv:1607.06292, prefactor (u-v) has been pulled out
double T4(double u, double cw2, double xH, double xA) noexcept
{
   const auto cw4 = cw2*cw2;
   const auto sw2 = 1.0 - cw2;
   const auto f5 = cw2*(5 - 16*cw2 + 8*cw4)/sw2;

   return std::log(u)/4*f5*(xA*(3 + 2*xH) - sqr(xA) + 3*xH - sqr(xH) - 3);
}

/// Eq.(76), arxiv:1607.06292
double T5(double u, double w, double cw2) noexcept
{
   const auto cw4 = cw2*cw2;
   const auto sw2 = 1.0 - cw2;
   const auto f6 = (7 - 14*cw2 + 4*cw4)/(4*cw2*sw2);
   const auto f8 = (13 - 20*cw2 + 4*cw4)/(cw2*sw2);

   return std::log(u)*(
      3.0/2*u + f6/cw2*(cube(u - w) + 3*cw2*sqr(u - w) + 3*cw4*(u - w))
      - 3.0/2*f8*u*(u - w) - cw2/2 - cw4);
}

/// Eq.(77), arxiv:1607.06292
double T6(double u, double w, double cw2) noexcept
{
   const auto cw4 = cw2*cw2;
   const auto u2 = u*u;

   return 9.0/2*(
      (u - w)*(u2 - 2*u*w + w*(w - cw2))/cw4*std::log(u/w)*std::log(w/cw2)
      + std::log(cw2)/cw2*(2*u2 + u*(cw2 - 4*w) - w*(cw2 - 2*w)));
}

/**
 * Eq (78), arxiv:1607.06292
 *
 * @note overall factor -1/2 is missing in arxiv:1607.06292v2
 */
double T7(double u, double w, double cw2) noexcept
{
   const auto cw4 = cw2*cw2;
   const auto sw2 = 1.0 - cw2;
   const auto f5 = cw2*(5 - 16*cw2 + 8*cw4)/sw2;
   const auto ra = std::complex<double>(1 + sqr(u - w) - 2*(u + w), 0.0);
   const auto s1 = u + w - 1.0 + std::sqrt(ra); // Eq.(79)

   const auto res =
      -0.5*f5*(2*(u + w) - sqr(u - w) - 1)*std::log(s1/(2*std::sqrt(u*w)))
      *(u + w - 1 - 4*u*w/s1);

   return std::real(res);
}

/// Eq.(80), arxiv:1607.06292
double T8(double u, double w, double cw2) noexcept
{
   const auto cw4 = cw2*cw2;
   const auto sw2 = 1.0 - cw2;
   const auto f6 = (7 - 14*cw2 + 4*cw4)/(4*cw2*sw2);
   const auto ra = std::complex<double>(sqr(u + w - cw2) - 4*u*w, 0.0);
   const auto s2 = u + w - cw2 + std::sqrt(ra); // Eq.(81)

   const auto res =
      2.0*f6*(4*u*w - sqr(u + w - cw2))*std::log(s2/(2*std::sqrt(u*w)))
      *((u + w)/cw2 - 4.0*u*w/(cw2*s2) - 1.0);

   return std::real(res);
}

/// Eq.(103), arxiv:1607.06292
double T9(double u, double w, double cw2) noexcept
{
   shift(w, cw2, eps_shift);

   const auto cw4 = cw2*cw2;
   const auto u2 = u*u;
   const auto w2 = w*w;

   // Note: Phi(u,w,w) == 0.5*u/w*f_PS(w/u)*(u - 4*w)

   return
      - 2*(cw4*w + cw2*(u2 + u*w - 2*w2) - cube(u-w))*Phi(u,w,cw2)
        /((cw2 - w)*(cw4 - 2*cw2*(u+w) + sqr(u-w)))
      + cw4*(u2 - 4*u*w + 2*w2)*u*f_PS(w/u)/(w*w2*(w-cw2))
      - 2*(cw2*u*(u-2*w) + w*sqr(u-w))*dilog(1.0 - u/w)/w2;
}

/// Eq.(104), arxiv:1607.06292
double T10(double u, double w, double cw2) noexcept
{
   shift(w, cw2, eps_shift);

   const auto u2 = u*u;
   const auto w2 = w*w;
   const auto lwu = std::log(w/u);
   const auto lwc = std::log(w/cw2);

   return
      (u2 - cw2*w - 2*u*w + w2)/(2*(cw2-w))*lwu*lwc
      + cw2*(cw2 + 2*u - 2*w)/(2*(cw2-w))*lwc
      + cw2*(u*lwu + w - u)/w;
}

/// Eq.(99), arxiv:1607.06292
double fb(double u, double w, double al, double cw2) noexcept
{
   return al*pi/(cw2*(-1.0 + cw2))*(u + 2*w);
}

/// Eq.(100), arxiv:1607.06292
double Fm0(double u, double w, double al, double cw2) noexcept
{
   return 1.0/(al*pi) * cw2*(-1 + cw2)/(u + 2*w) * YF1(u,w,cw2);
}

/// Eq.(101), arxiv:1607.06292
double Fmp(double u, double w, double al, double cw2) noexcept
{
   return (-9*(-1 + cw2))/(al*pi) * (T9(u,w,cw2)/2 + T10(u,w,cw2));
}

} // anonymous namespace

/**
 * Calculates 2-loop bosonic pure electroweak contributions.
 *
 * Eq (49), arxiv:1607.06292
 */
double amu2L_B_EWadd(const THDM_B_parameters& thdm) noexcept
{
   const double zeta2 = 1.6449340668482264; // Zeta[2]
   const double mw2 = sqr(thdm.mw);
   const double mz2 = sqr(thdm.mz);
   const double cw2 = mw2/mz2;
   const double cw4  = cw2*cw2;
   const double cw6  = cw4*cw2;
   const double cw8  = cw4*cw4;
   const double cw10 = cw8*cw2;
   const double cw12 = cw8*cw4;
   const double cw14 = cw8*cw6;
   const double mh2 = sqr(thdm.mh(0));
   const double xh_ = mh2/mz2; // temporary value
   const double xh = is_equal_rel(1.0, xh_, eps_shift) ? xh_*(1 + eps_shift) : xh_;
   const double xw = xh/cw2;
   const double lh = std::log(xh);
   const double lw = std::log(cw2);
   const double liw = dilog(1 - xw);
   const double lih = dilog(1 - xh);
   const double lh2 = lh*lh;
   const double phi1 = 3*xh*f_PS(1/xh)*(xh - 4); // Phi(xh,1,1) == 0.5*xh*f_PS(1/xh)*(xh - 4)
   const double phi3 = 6*(-liw + zeta2);
   const double phi4 = 6*(lh2/2 + 2*lih + zeta2);
   const double phi5 = 3/cw2*f_PS(cw2); // Phi(cw2,cw2,1) == 0.5/cw2*f_PS(cw2)*(1 - 4*cw2)
   const double phi6 = 6*(-lih + zeta2);
   const double phi7 = 3*cw2*xw*f_PS(1/xw)*(xw - 4); // Phi(xw,1,1) == 0.5*xw*f_PS(1/xw)*(xw - 4)

   const double xm2 = 256*(32*cw10 - 5*cw4 + 32*cw6 - 56*cw8)*phi6
      - 2304*(5*cw10 - 4*cw12 - cw8)*phi7 - 512*(cw10 - 4*cw12)*phi3;

   const double xm1 = 64*(24*(144*cw12 + 5*cw4 - 32*cw6 + 94*cw8) -
      2*((5*cw4 - 32*cw6)*(12*lh + phi1) +
         8*cw10*(330 + 147*lw - 195*lh - 4*phi1 - phi3) +
         16*cw12*(-54*lw + 54*lh + phi3) +
         cw8*(-240*lw + 912*lh + 56*phi1 + phi3)) +
      (-32*cw10 + 10*cw2 - 59*cw4 + 80*cw6 - 8*cw8)*phi6) -
      4*(3072*cw10 + 907*cw6 - 4396*cw8)*phi7;

   const double x0 = 96*(730*cw10 - 936*cw12 + 384*cw14 + 21*cw6 - 211*cw8)*phi5 +
    4*(231*cw4 - 1037*cw6 + 452*cw8)*phi7 +
    16*(-837*cw6 + 852*cw8 + 6672*cw10*lw - 6912*cw12*(2 + lw) +
       20*cw2*(-12 + 12*lh + phi1) - 5*phi6 + 22*cw2*phi6 -
       32*cw10*(-468 + 42*lh + 4*phi1 + phi3 + 12*phi6) +
       2*cw4*(468 - 588*lh - 54*phi1 + 34*phi6) +
       4*cw8*(306*lw - 141*lh + 24*phi1 - 2*phi3 + 184*phi6 - 6*pi2) +
       cw6*(-561*lw + 1089*lh + 96*phi1 + 4*phi3 - 464*phi6 + 6*pi2));

   const double x1 = 32*cw2*(54 + 6*lh + 3*phi1 - 19*phi6) + 90*cw2*phi7 +
    64*cw10*(1824 + 684*lw + 384*lh - 128*phi1 + 627*phi5 - 128*phi6 +
       128*pi2) - 4*(120*lh + 10*phi1 + 3648*cw12*phi5 - 5*(24 + phi6) +
       16*cw8*(2853 + 768*lw + 291*lh - 240*phi1 - 4*phi3 + 519*phi5 -
          272*phi6 + 218*pi2) + cw4*
        (5190 + 345*lw - 3081*lh - 416*phi1 - 134*phi3 + 252*phi5 - 1096*phi6 +
          33*phi7 + 412*pi2) - 4*cw6*
        (3939 + 780*lw - 1158*lh - 560*phi1 - 138*phi3 + 615*phi5 - 808*phi6 -
          57*phi7 + 598*pi2));

   const double x2 = cw2*(3867 + 426*lw - 2106*lh - 1392*lh2 - 672*phi1 - 262*phi3 +
       464*phi4 + 126*phi5 - 70*phi7 - 586*pi2) -
    128*cw10*(144 + 288*lh - 96*lh2 - 72*phi1 + 128*phi4 - 3*phi5 - 32*pi2) +
    2*cw4*(-5322 - 144*lw + 1434*lh + 2472*lh2 + 1392*phi1 + 256*phi3 -
       56*phi4 - 561*phi5 + 396*phi7 + 916*pi2) -
    4*(-5*(-30 + 18*lh + phi1 + 3*phi6) + 8*phi7 +
       4*cw8*(72*lw - 3264*lh + 960*lh2 + 752*phi1 - 1664*phi4 + 201*phi5 +
          320*(-3 + pi2)) + cw6*(-4104 + 696*lw + 2892*lh + 48*lh2 - 192*phi1 -
          536*phi3 + 2672*phi4 - 867*phi5 + 672*pi2));

   const double x3 = -24 - 60*lh + 102*lh2 + 68*phi1 - 3072*cw10*phi1 + 32*phi3 -
    34*phi4 + 15360*cw10*phi4 - 6*(3*cw2 - 19*cw4 + 50*cw6 - 40*cw8)*phi5 +
    32*phi7 - 16*cw6*(45*lw + 8*(123 + 240*lh - 78*lh2 - 38*phi1 + 5*phi4 -
          26*pi2)) + 4*cw4*(1683 + 417*lw + 3222*lh - 1056*lh2 - 688*phi1 -
       262*phi3 + 1216*phi4 - 90*pi2) +
    256*cw8*(36 + 72*lh - 24*lh2 + 3*phi1 - 73*phi4 - 8*pi2) + 2*pi2 -
    cw2*(747 + 426*lw + 1350*lh - 120*lh2 - 112*phi1 - 134*phi3 + 808*phi4 +
       128*phi7 + 94*pi2);

   const double x4 = -2*(-72 - 144*lh + 51*lh2 + 36*phi1 + 16*phi3 - 65*phi4 +
      1536*cw10*phi4 + 32*(-24*cw8*phi1 + 36*cw8*phi4 +
         cw6*(18 + 36*lh - 12*lh2 + 33*phi1 - 152*phi4 - 4*pi2)) -
      8*cw4*(126 + 252*lh - 84*lh2 + 21*phi1 - 284*phi4 - 28*pi2) +
      4*cw2*(126 + 252*lh - 87*lh2 - 39*phi1 - 16*phi3 - 7*phi4 - 13*pi2) + pi2);

   const double x5 = 24*((-5 + 27*cw2 - 14*cw4 - 72*cw6 + 64*cw8)*phi4
      + (1 - 7*cw2 + 14*cw4 - 8*cw6)*phi1);

   const double x6 = -24*(-1 + 7*cw2 - 14*cw4 + 8*cw6)*phi4;

   const double res = (xm2/xh + xm1)/xh + x0 + xh*(x1 + xh*(x2 + xh*(x3 + xh*(x4 + xh*(x5 + xh*x6)))));

   const double pref = sqr(thdm.alpha_em*thdm.mm/(48*thdm.mz*pi*cw2*(4*cw2 - xh)*(1 - cw2)))
      /(2*(-1 + 4*cw2)*(-1 + xh));

   return pref * res * thdm.cos_beta_minus_alpha * thdm.zetal;
}

/**
 * Calculates 2-loop bosonic non-Yukawa contributions.
 *
 * Eq (71), arxiv:1607.06292
 */
double amu2L_B_nonYuk(const THDM_B_parameters& thdm) noexcept
{
   const auto mw2 = sqr(thdm.mw);
   const auto mz2 = sqr(thdm.mz);
   const auto cw2 = mw2/mz2;
   const auto cw4 = cw2*cw2;
   const auto cw6 = cw4*cw2;
   const auto cw8 = cw4*cw4;
   const auto sw2 = 1.0 - cw2;
   const auto sw4 = sw2*sw2;
   const auto mA2 = sqr(thdm.mA);
   const auto mH2 = sqr(thdm.mh(1));
   const auto mHp2 = sqr(thdm.mHp);
   const auto xH = mH2/mz2;
   const auto xA = mA2/mz2;
   const auto xHp = mHp2/mz2;

   const auto f1 = 7.0/2 - 25.0/(2*cw2) + 4*cw2 - 4*cw4;
   const auto f2 = 2*(17 - 24*cw2 + 56*cw4 - 128*cw6 + 64*cw8);
   const auto f3 = (25 - 32*cw2 + 4*cw4)/(cw2*sw2);
   const auto f4 = 13.0/2 - 15*cw2 + 10*cw4;
   const auto f5 = cw2*(5 - 16*cw2 + 8*cw4)/sw2;

   const double res =
      + TX(xH, xA, xHp, cw2)
      + (xA - xH)*(T4(xA, cw2, xH, xA) - T4(xH, cw2, xH, xA))
      + T5(xHp, xH, cw2)
      + T5(xHp, xA, cw2)
      + T6(xA, xHp, cw2)
      + T6(xH, xHp, cw2)
      + T7(xA, xH, cw2)
      + T7(xHp, xHp, cw2)*sqr(1 - 2*cw2)
      + T8(xA, xHp, cw2)
      + T8(xH, xHp, cw2)
      - 16.0/3*cw2*sw2*(1 + 8*cw2 - 8*cw4)
      + 8*cw4*sw4/(5*xHp)
      + f2*xHp
      - f3*sqr(xHp)
      + f1*(sqr(xA) + sqr(xH))
      + f3*xHp*(xA + xH)
      + f4*(xA + xH)
      - f5*xA*xH
      + T1(xA, xHp, cw2)
      + T1(xH, xHp, cw2)
      + T0(xA, xHp, cw2)
      + T0(xH, xHp, cw2)
      ;

   const auto pref = sqr(thdm.alpha_em/(24*pi*cw2*sw2) * thdm.mm/thdm.mz);

   return pref*res;
}

/**
 * Calculates 2-loop bosonic Yukawa contributions.
 *
 * Eq (52), arxiv:1607.06292
 */
double amu2L_B_Yuk(const THDM_B_parameters& thdm) noexcept
{
   const auto tb = thdm.tb;
   const auto sc = tb - 1.0/tb;
   const auto zetal = thdm.zetal;
   const auto lambda5 = thdm.lambda5;
   const auto lambda67 = thdm.lambda67;
   const auto cos_beta_minus_alpha = thdm.cos_beta_minus_alpha;
   const auto al = thdm.alpha_em;
   const auto mhSM2 = sqr(thdm.mhSM);
   const auto mH2 = sqr(thdm.mh(1));
   const auto mHp2 = sqr(thdm.mHp);
   const auto mw2 = sqr(thdm.mw);
   const auto mz2 = sqr(thdm.mz);

   const auto cw2 = mw2/mz2;
   const auto xhSM = mhSM2/mz2;
   const auto xHp = mHp2/mz2;
   const auto xH = mH2/mz2;

   // Eq.(91), arxiv:1607.06292
   const auto a000 = fb(xhSM, xHp, al, cw2)*Fm0(xhSM, xHp, al, cw2);
   // Eq.(92), arxiv:1607.06292
   const auto a0z0 = -fb(xH, 0.0, al, cw2)*(
      Fm0(xH, xHp, al, cw2) + Fmp(xH, xHp, al, cw2));
   // Eq.(93), arxiv:1607.06292
   const auto a500 = Fm0(xhSM, xHp, al, cw2);
   // Eq.(94), arxiv:1607.06292
   const auto a5z0 = -0.5*(
      Fm0(xH, xHp, al, cw2) + Fmp(xH, xHp, al, cw2));
   // Eq.(95), arxiv:1607.06292
   const auto a001 =
      + fb(xH, 0.0, al, cw2)*Fm0(xH, xHp, al, cw2)
      - fb(xhSM, 0.0, al, cw2)*Fm0(xhSM, xHp, al, cw2);
   // Eq.(96), arxiv:1607.06292
   const auto a0z1 =
      - (
         + fb(xH, xHp, al, cw2)*(
            + Fm0(xH, xHp, al, cw2)
            + Fmp(xH, xHp, al, cw2)
         )
         - YF3(xH, xHp, cw2)
         // SM contributions with opposite sign:
         - fb(xhSM, xHp, al, cw2)*(
            + Fm0(xhSM, xHp, al, cw2)
            + Fmp(xhSM, xHp, al, cw2)
         )
         + YF3(xhSM, xHp, cw2)
      )
      + YF2(xH, cw2);
   // Eq.(97), arxiv:1607.06292
   const auto a501 =
      + Fm0(xH, xHp, al, cw2)/2
      - Fm0(xhSM, xHp, al, cw2)/2;
   // Eq.(98), arxiv:1607.06292
   const auto a5z1 =
      - Fm0(xH, xHp, al, cw2)
      - Fmp(xH, xHp, al, cw2)
      + Fm0(xhSM, xHp, al, cw2)
      + Fmp(xhSM, xHp, al, cw2);

   const double res =
      + a000
      + a0z0*sc*zetal
      + a500*lambda5
      + a5z0*(sc*lambda5 + lambda67)*zetal
      + (
         + a001*sc
         + a0z1*zetal
         + a501*(sc*lambda5 + lambda67)
         + a5z1*lambda5*zetal
         )*cos_beta_minus_alpha;

   const auto pref = sqr(al/(24*pi*cw2*(1.0 - cw2)) * thdm.mm/thdm.mz);

   return pref*res;
}

/**
 * Calculates the sum of the 2-loop bosonic contributions.
 *
 * Eq (48), arxiv:1607.06292
 */
double amu2L_B(const THDM_B_parameters& thdm) noexcept
{
   return amu2L_B_EWadd(thdm) + amu2L_B_nonYuk(thdm) + amu2L_B_Yuk(thdm);
}

} // namespace thdm

} // namespace gm2calc
