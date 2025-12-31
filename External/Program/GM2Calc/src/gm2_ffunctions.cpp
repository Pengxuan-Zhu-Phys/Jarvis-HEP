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

#include "gm2_ffunctions.hpp"
#include "gm2_dilog.hpp"
#include "gm2_log.hpp"
#include "gm2_numerics.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>

namespace gm2calc {

namespace {
   constexpr double eps = 10.0*std::numeric_limits<double>::epsilon();
   const double qdrt_eps = std::pow(eps, 0.25);

   /// shift values symmetrically away from equality, if they are close
   void shift(double& x, double& y, double rel_diff) noexcept
   {
      if (is_equal_rel(x, y, rel_diff)) {
         const double mid = 0.5*std::abs(y + x);
         if (x < y) {
            x = (1 - rel_diff)*mid;
            y = (1 + rel_diff)*mid;
         } else {
            x = (1 + rel_diff)*mid;
            y = (1 - rel_diff)*mid;
         }
      }
   }

   void sort(double& x, double& y) noexcept
   {
      if (x > y) { std::swap(x, y); }
   }

   void sort(double& x, double& y, double& z) noexcept
   {
      if (x > y) { std::swap(x, y); }
      if (y > z) { std::swap(y, z); }
      if (x > y) { std::swap(x, y); }
   }

   /// calculates phi(xd, xu, 1)/y with y = (xu - xd)^2 - 2*(xu + xd) + 1,
   /// properly handle the case y = 0
   double phi_over_y(double xu, double xd) noexcept
   {
      const double sqrtxd = std::sqrt(xd);
      const double ixd = 1/xd;
      constexpr double eps = 1e-8;

      // test two cases where y == 0
      if (std::abs((xu - 1)*ixd + 2/sqrtxd - 1) < eps) {
         return -std::log(std::abs(-1 + sqrtxd))/sqrtxd + std::log(xd)/(2*(-1 + sqrtxd));
      } else if (std::abs((xu - 1)*ixd - 2/sqrtxd - 1) < eps) {
         return std::log(1 + sqrtxd)/sqrtxd - std::log(xd)/(2*(1 + sqrtxd));
      }

      const double y = sqr(xu - xd) - 2*(xu + xd) + 1;
      const double phi = Phi(xd, xu, 1);

      return phi/y;
   }

   /// lambda^2(u,v)
   double lambda_2(double u, double v) noexcept
   {
      return sqr(1 - u - v) - 4*u*v;
   }

   /// expansion of (1 - lambda + u - v)/2 for u ~ v ~ 0 up to including O(u^3 v^3)
   double l00(double u, double v) noexcept
   {
      return v*(1 + u*(1 + u*(1 + u)) + v*(u*(1 + u*(3 + 6*u)) + u*(1 + u*(6 + 20*u))*v));
   }

   /// expansion of (1 - lambda + u - v)/2 for u ~ 0 and v < 1 up to including O(u^3 v^3)
   double l0v(double u, double v) noexcept
   {
      const double a = 1 - v;
      const double a2 = a*a;
      const double a3 = a2*a;
      return u*(0.5*(1 + (1 + v)/a) + u*(v + u*v*(1 + v)/a2)/a3);
   }

   /// expansion of (1 - lambda - u + v)/2 for u ~ 0 and v < 1 up to including O(u^3 v^3)
   double lv0(double u, double v) noexcept
   {
      const double a = 1 - v;
      const double a2 = a*a;
      const double a3 = a2*a;
      return v + u*(0.5*(-1 + (1 + v)/a) + u*(v + u*v*(1 + v)/a2)/a3);
   }

   /// returns tuple (0.5*(1 - lambda + u - v), 0.5*(1 - lambda - u + v))
   std::tuple<double,double> luv(double lambda, double u, double v) noexcept
   {
      if (v < qdrt_eps) {
         return std::make_tuple(l00(u, v), l00(v, u));
      } else if (u < qdrt_eps) {
         return std::make_tuple(l0v(u, v), lv0(u, v));
      }
      return std::make_tuple(0.5*(1 - lambda + u - v),
                             0.5*(1 - lambda - u + v));
   }

   /// u < 1 && v < 1, lambda^2(u,v) > 0; note: phi_pos(u,v) = phi_pos(v,u)
   double phi_pos(double u, double v) noexcept
   {
      if (is_equal_rel(u, 1.0, eps) && is_equal_rel(v, 1.0, eps)) {
         return 2.343907238689459;
      }

      const double pi23 = 3.2898681336964529; // Pi^2/3
      const auto lambda = std::sqrt(lambda_2(u,v));

      if (is_equal_rel(u, v, eps)) {
         const double x = u < qdrt_eps ? u*(1 + u*(1 + u*(2 + 5*u))) : 0.5*(1 - lambda);

         return (- sqr(std::log(u)) + 2*sqr(std::log(x))
                 - 4*dilog(x) + pi23)/lambda;
      }

      double x = 0, y = 0;
      std::tie(x, y) = luv(lambda, u, v);

      return (- std::log(u)*std::log(v) + 2*std::log(x)*std::log(y)
              - 2*dilog(x) - 2*dilog(y) + pi23)/lambda;
   }

   /// clausen_2(2*acos(x))
   double cl2acos(double x) noexcept
   {
      return clausen_2(2*std::acos(x));
   }

   /// lambda^2(u,v) < 0, u = 1
   double phi_neg_1v(double v) noexcept
   {
      return 2*(cl2acos(1 - 0.5*v) + 2*cl2acos(0.5*std::sqrt(v)));
   }

   /// lambda^2(u,v) < 0; note: phi_neg(u,v) = phi_neg(v,u)
   double phi_neg(double u, double v) noexcept
   {
      if (is_equal_rel(u, 1.0, eps) && is_equal_rel(v, 1.0, eps)) {
         // -I/9 (Pi^2 - 36 PolyLog[2, (1 - I Sqrt[3])/2])/Sqrt[3]
         return 2.343907238689459;
      }

      const auto lambda = std::sqrt(-lambda_2(u,v));

      if (is_equal_rel(u, v, eps)) {
         return 4*clausen_2(2*std::asin(std::sqrt(0.25/u)))/lambda;
      }

      if (is_equal_rel(u, 1.0, eps)) {
         return phi_neg_1v(v)/lambda;
      }

      if (is_equal_rel(v, 1.0, eps)) {
         return phi_neg_1v(u)/lambda;
      }

      const auto sqrtu = std::sqrt(u);
      const auto sqrtv = std::sqrt(v);

      return 2*(+ cl2acos(0.5*(1 + u - v)/sqrtu)
                + cl2acos(0.5*(1 - u + v)/sqrtv)
                + cl2acos(0.5*(-1 + u + v)/(sqrtu*sqrtv)))/lambda;
   }

   /**
    * Phi(u,v) with u = x/z, v = y/z.
    *
    * The following identities hold:
    * Phi(u,v) = Phi(v,u) = Phi(1/u,v/u)/u = Phi(1/v,u/v)/v
    */
   double phi_uv(double u, double v) noexcept
   {
      const auto lambda = lambda_2(u,v);

      if (is_zero(lambda, eps)) {
         // phi_uv is always multiplied by lambda.  So, in order to
         // avoid nans if lambda == 0, we simply return 0
         return 0.0;
      }

      if (lambda > 0.) {
         if (u <= 1 && v <= 1) {
            return phi_pos(u,v);
         }
         const auto vou = v/u;
         if (u >= 1 && vou <= 1) {
            const auto oou = 1/u;
            return phi_pos(oou,vou)*oou;
         }
         // v >= 1 && u/v <= 1
         const auto oov = 1/v;
         return phi_pos(oov,1/vou)*oov;
      }

      return phi_neg(u,v);
   }

} // anonymous namespace

double F1C(double x) noexcept {
   if (is_zero(x, eps)) {
      return 4.0;
   }

   const double d = x - 1.0;

   if (is_equal_rel(x, 1.0, 0.03)) {
      return 1.0 + d*(-0.6 + d*(0.4 + d*(-2.0/7.0
         + d*(3.0/14.0 + d*(-1.0/6.0
         + 2.0/15.0*d)))));
   }

   return 2.0/pow4(d)*(2.0 + x*(3.0 + 6.0*std::log(x) + x*(-6.0 + x)));
}

double F2C(double x) noexcept {
   if (is_zero(x, eps)) {
      return 0.0;
   }

   if (is_equal_rel(x, 1.0, 0.03)) {
      const double d = x - 1.0;

      return 1.0 + d*(-0.75 + d*(0.6 + d*(-0.5 + d*(3.0/7.0
         + d*(-0.375 + 1.0/3.0*d)))));
   }

   return 3.0/(2.0*pow3(1.0 - x))*(-3.0 - 2.0*std::log(x) + x*(4.0 - x));
}

double F3C(double x) noexcept {
   const double d = x - 1.0;

   if (is_equal_rel(x, 1.0, 0.03)) {
      return 1.0
         + d*(1059.0/1175.0
         + d*(-4313.0/3525.0
         + d*(70701.0/57575.0
         + d*(-265541.0/230300.0
         + d*(+48919.0/46060.0
         - 80755.0/82908.0*d)))));
   }

   const double lx = std::log(x);
   const double x2 = sqr(x);

   return 4.0/(141.0*pow4(d)) * (
      + (1.0 - x) * (151.0 * x2 - 335.0 * x + 592.0)
      + 6.0 * (21.0 * pow3(x) - 108.0 * x2 - 93.0 * x + 50.0) * lx
      - 54.0 * x * (x2 - 2.0 * x - 2.0) * sqr(lx)
      - 108.0 * x * (x2 - 2.0 * x + 12.0) * dilog(1.0 - x)
      );
}

double F4C(double x) noexcept {
   if (is_zero(x, eps)) {
      return 0.0;
   }

   if (is_equal_rel(x, 1.0, 0.03)) {
      const double d = x - 1.0;

      return 1.0
         + d*(-45.0/122.0
         + d*(941.0/6100.0
         + d*(-17.0/305.0
         + d*(+282.0/74725.0
         + d*(+177.0/6832.0 - 47021.0/1076040.0*d)))));
   }

   const double lx = std::log(x);
   const double x2 = sqr(x);

   return -9.0/(122.0 * pow3(1.0 - x)) * (
      + 8.0 * (x2 - 3.0 * x + 2.0)
      + (11.0 * x2 - 40.0 * x + 5.0) * lx
      - 2.0 * (x2 - 2.0 * x - 2.0) * sqr(lx)
      - 4.0 * (x2 - 2.0 * x + 9.0) * dilog(1.0 - x)
      );
}

double F1N(double x) noexcept {
   if (is_zero(x, eps)) {
      return 2.0;
   }

   const double d = x - 1.0;

   if (is_equal_rel(x, 1.0, 0.03)) {
      return 1.0 + d*(-0.4 + d*(0.2 + d*(-4.0/35.0
         + d*(1.0/14.0 + d*(-1.0/21.0 + 1.0/30.0*d)))));
   }

   return 2.0/pow4(d)*(1.0 + x*(-6.0 + x*(+3.0 - 6.0 * std::log(x) + 2.0 * x)));
}

double F2N(double x) noexcept {
   if (is_zero(x, eps)) {
      return 3.0;
   }

   if (is_equal_rel(x, 1.0, 0.04)) {
      const double d = x - 1.0;

      return 1. + d*(-0.5 + d*(0.3 + d*(-0.2
         + d*(1.0/7.0 + d*(-3.0/28.0 + 1.0/12.0*d)))));
   }

   return 3.0/pow3(1.0 - x) * (1.0 + x*(2.0 * std::log(x) - x));
}

double F3N(double x) noexcept {
   if (is_zero(x, eps)) {
      return 8.0/105.0;
   }

   const double d = x - 1.0;

   if (is_equal_rel(x, 1.0, 0.03)) {
      return 1.0 + d*(76/875.0 + d*(-431/2625.0 + d*(5858/42875.0
         + d*(-3561/34300.0 + d*(23/294.0 - 4381/73500.0*d)))));
   }

   const double x2 = sqr(x);

   return 4.0/105.0/pow4(d) * (
      + (1.0 - x) * (-97.0 * x2 - 529.0 * x + 2.0)
      + 6.0 * x2 * (13.0 * x + 81.0) * std::log(x)
      + 108.0 * x * (7.0 * x + 4.0) * dilog(1.0 - x)
      );
}

double F4N(double x) noexcept {
   const double PI2 = 9.8696044010893586; // Pi^2

   if (is_zero(x, eps)) {
      return -3.0/4.0*(-9.0 + PI2);
   }

   if (is_equal_rel(x, 1.0, 0.03)) {
      const double d = x - 1.0;

      return 1.0 + sqr(d)*(-111.0/800.0 + d*(59.0/400.0 + d*(-129.0/980.0
         + d*(177.0/1568.0 - 775.0/8064.0*d))));
   }

   return -2.25/pow3(1.0 - x) * (
      + (x + 3.0) * (x * std::log(x) + x - 1.0)
      + (6.0 * x + 2.0) * dilog(1.0 - x)
      );
}

namespace {

/// expansion of Fb(x,y) around x ~ 1 and y ~ 1 up to including O((x-1)^2 (y-1)^2)
double Fb11(double x, double y) noexcept {
   const double x1 = x - 1;
   const double y1 = y - 1;

   return
      + 1.0/12 + (-0.05 + 1.0/30*y1)*y1
      + x1*(-0.05 + (1.0/30 - 1.0/42*y1)*y1
      + x1*(1.0/30 + (-1.0/42 + 1.0/56*y1)*y1));
}

/// expansion of Fb(x,y) around y ~ x, x != 0
double Fbx(double x, double y) noexcept {
   if (is_equal_rel(x, 1.0, 1e-2)) {
      const double d = x - 1;
      return 1.0/12 + d*(-0.1 + d*(0.1 + d*(-2.0/21 + d*(5.0/56 + d*(-1.0/12 + d*(7.0/90 - 4.0/55*d))))));
   }

   const double x1 = x - 1.0;
   const double d = y - x;
   const double lx = std::log(x);
   const double x14 = pow4(x1);
   const double x15 = x14*x1;
   const double x16 = x15*x1;

   return (-5 - 2*lx + x*(4 - 4*lx + x))/(2*x14)
      - d*(-1 + x*(-9 - 6*lx + x*(9 - 6*lx + x)))/(2*x15*x)
      - sqr(d)*(-1 + x*(12 + x*(36 + 36*lx + x*(-44 + 24*lx - 3*x))))/(6*x16*sqr(x));
}

} // anonymous namespace

double Fb(double x, double y) noexcept {
   if (x < 0 || y < 0) {
      ERROR("Fb: x and y must not be negative!");
      return std::numeric_limits<double>::quiet_NaN();
   }

   sort(x, y);

   if (is_zero(y, eps)) {
      return 0;
   }

   if (is_equal_rel(x, 1.0, 1e-4) && is_equal_rel(y, 1.0, 1e-4)) {
      return Fb11(x, y);
   }

   if (is_equal_rel(x, y, 1e-5)) {
      return Fbx(x, y);
   }

   return (G4(y) - G4(x))/(x - y);
}

namespace {

/// expansion of Fa(x,y) around x ~ 1 and y ~ 1 up to including O((x-1)^2 (y-1)^2)
double Fa11(double x, double y) noexcept {
   const double x1 = x - 1;
   const double y1 = y - 1;

   return
      0.25 + (-0.2 + 1.0/6*y1)*y1
      + x1*(-0.2 + (1.0/6 - 1.0/7*y1)*y1
      + x1*(1.0/6 + (-1.0/7 + 1.0/8*y1)*y1));
}

/// expansion of Fa(x,y) around y ~ x, x != 0
double Fax(double x, double y) noexcept {
   if (is_equal_rel(x, 1.0, 1e-2)) {
      const double d = x - 1;
      return 0.25 + d*(-0.4 + d*(0.5 + d*(-4.0/7 + d*(5.0/8 + d*(-2./3 + d*(0.7 - 8.0/11*d))))));
   }

   const double x1 = x - 1.0;
   const double d = y - x;
   const double lx = std::log(x);
   const double x14 = pow4(x1);
   const double x15 = x14*x1;
   const double x16 = x15*x1;
   const double x2 = sqr(x);
   const double x3 = x2*x;

   return (2 + x*(3 + 6*lx + x*(-6 + x)))/(2*x14*x)
      - d*(-1 + x*(8 + x*(12*lx + x*(-8 + x))))/(2*x15*x2)
      - sqr(d)*(-2 + x*(15 + x*(-60 + x*(20 - 60*lx + x*(30 - 3*x)))))/(6*x16*x3);
}

} // anonymous namespace

double Fa(double x, double y) noexcept {
   if (x < 0 || y < 0) {
      ERROR("Fa: x and y must not be negative!");
      return std::numeric_limits<double>::quiet_NaN();
   }

   sort(x, y);

   if (is_zero(y, eps)) {
      return 0;
   }

   if (is_equal_rel(x, 1.0, 1e-4) && is_equal_rel(y, 1.0, 1e-4)) {
      return Fa11(x,y);
   }

   if (is_equal_rel(x, y, 1e-5)) {
      return Fax(x,y);
   }

   return (G3(y) - G3(x))/(x - y);
}

double G3(double x) noexcept {
   if (is_equal_rel(x, 1.0, 1e-2)) {
      const double d = x - 1;
      return 1.0/3 + d*(-0.25 + d*(0.2 + d*(-1.0/6 + d*(1.0/7 + d*(-1.0/8
         + d*(1.0/9 + d*(-0.1 + 1.0/11*d)))))));
   }

   return ((x - 1)*(x - 3) + 2*std::log(x))/(2*pow3(x - 1));
}

double G4(double x) noexcept {
   if (is_equal_rel(x, 1.0, 1e-2)) {
      const double d = x - 1;
      return 1.0/6 + d*(-1.0/12 + d*(0.05 + d*(-1.0/30 + d*(1.0/42
         + d*(-1.0/56 + d*(1.0/72 + d*(-1.0/90 + 1.0/110*d)))))));
   }

   return ((x - 1)*(x + 1) - 2*x*std::log(x))/(2*pow3(x - 1));
}

namespace {

/// Ixy(0,y), squared arguments, y != 0
double I0y(double y) noexcept {
   if (is_equal_rel(y, 1.0, eps)) {
      const double d = y - 1;
      return 1 + d*(-0.5 + 1./3*d);
   }

   return std::log(y)/(y - 1);
}

/// I(x,y), squared arguments, x == 1, y != 0
double I1y(double x, double y) noexcept {
   const double dy = y - 1;
   const double dy2 = sqr(dy);
   const double dx = (x - 1)/dy2;
   const double y2 = sqr(y);
   const double yly = y*std::log(y);

   return (1 - y + yly)/dy2
      + dx*(0.5 - 0.5*y2 + yly)/dy
      + sqr(dx)*(1./3 + 0.5*y + yly + y2*(1./6*y - 1));
}

/// I(x,y), squared arguments, x == y, x != 0, y != 0
double Ixx(double x, double y) noexcept {
   const double eps_eq = 0.0001;

   if (is_equal_rel(y, 1.0, eps_eq)) {
      const double dx = x - 1;
      const double dy = y - 1;
      const double dy2 = sqr(dy);

      return 0.5 + dx*(-1./6 + 1./12*dy - 1./20*dy2)
         + sqr(dx)*(1./12 - 1./20*dy + 1./30*dy2)
         - 1./6*dy + 1./12*dy2;
   }

   const double y2 = sqr(y);
   const double dy = y - 1;
   const double dy2 = sqr(dy);
   const double dxy = (x - y)/dy2;
   const double ly = std::log(y);

   return (dy - ly)/dy2
      + dxy*(0.5 - 0.5*y2 + y*ly)/(dy*y)
      + sqr(dxy)*(1./6 - y + y2*(0.5 + 1./3*y - ly))/y2;
}

/// I(x,y), x < y, x and y are squared arguments
double Ixy(double x, double y) noexcept {
   const double eps_eq = 0.0001;

   if (is_zero(y, eps)) {
      return 0;
   }

   if (is_zero(x, eps)) {
      return I0y(y);
   }

   if (is_equal_rel(x/y, 1.0, eps_eq)) {
      return Ixx(x, y);
   }

   if (is_equal_rel(x, 1.0, eps_eq)) {
      return I1y(x, y);
   }

   if (is_equal_rel(y, 1.0, eps_eq)) {
      return I1y(y, x);
   }

   const double lx = std::log(x);
   const double ly = std::log(y);

   return (x*(y - 1)*lx - y*(x - 1)*ly)/((x - 1)*(x - y)*(y - 1));
}

/// I(x,y,z), x, y and z are squared arguments
double Ixyz(double x, double y, double z) noexcept {
   sort(x, y, z);

   if (is_zero(z, eps)) {
      return 0;
   }

   return Ixy(x/z, y/z)/z;
}

} // anonymous namespace

double Iabc(double a, double b, double c) noexcept {
   return Ixyz(sqr(a), sqr(b), sqr(c));
}

/**
 * Calculates \f$f_{PS}(z)\f$, Eq (70) arXiv:hep-ph/0609168
 * @author Alexander Voigt
 */
double f_PS(double z) noexcept {
   if (z < 0.0) {
      ERROR("f_PS: z must not be negative!");
      return std::numeric_limits<double>::quiet_NaN();
   } else if (z == 0.0) {
      return 0.0;
   } else if (z < std::numeric_limits<double>::epsilon()) {
      const double pi23 = 3.2898681336964529; // Pi^2/3
      const double lz = std::log(z);
      return z*(pi23 + lz*lz);
   } else if (z < 0.25) {
      const double y = std::sqrt(1 - 4*z); // 0 < y < 1
      const double c = -9.8696044010893586; // -Pi^2
      const double q = (1 + y)/(1 - y);
      const double lq = std::log(q);
      return z/y*(4*dilog(1 + q) - lq*(2*std::log(z) - lq) + c);
   } else if (z == 0.25) {
      return 1.3862943611198906; // Log[4]
   }

   // z > 0.25
   const double y = std::sqrt(-1 + 4*z);
   const double theta = std::atan2(y, 2*z - 1);
   return 4*z/y*clausen_2(theta);
}

/**
 * Calculates \f$f_S(z)\f$, Eq (71) arXiv:hep-ph/0609168
 */
double f_S(double z) noexcept {
   if (z < 0.0) {
      ERROR("f_S: z must not be negative!");
      return std::numeric_limits<double>::quiet_NaN();
   } else if (z == 0.0) {
      return 0.0;
   } if (z > 1e2) {
      const double lz = std::log(z);
      const double iz = 1/z;
      return (-13./9 - 2./3*lz) + iz*(-26./150 - 15./150*lz
         + iz*(-673./22050 - 420./22050*lz + iz*(-971./158760 - 630./158760*lz)));
   }

   return (2*z - 1)*f_PS(z) - 2*z*(2 + std::log(z));
}

/**
 * Calculates \f$f_{\tilde{f}}(z)\f$, Eq (72) arXiv:hep-ph/0609168
 */
double f_sferm(double z) noexcept {
   if (z < 0.0) {
      ERROR("f_sferm: z must not be negative!");
      return std::numeric_limits<double>::quiet_NaN();
   } else if (z == 0.0) {
      return 0.0;
   }

   return 0.5*z*(2 + std::log(z) - f_PS(z));
}

/**
 * Calculates Barr-Zee 2-loop function for diagram with lepton loop
 * and charged Higgs and W boson mediators, Eq (60), arxiv:1607.06292,
 * with extra global prefactor z.
 */
double f_CSl(double z) noexcept {
   if (z < 0.0) {
      ERROR("f_CSl: z must not be negative!");
      return std::numeric_limits<double>::quiet_NaN();
   } else if (z == 0) {
      return 0.0;
   }

   constexpr double pi26 = 1.6449340668482264;

   return z*(z + z*(z - 1)*(dilog(1 - 1/z) - pi26) + (z - 0.5)*std::log(z));
}

/**
 * Eq (61), arxiv:1607.06292, with extra global prefactor xd
 *
 * @note There is a misprint in Eq (61), arxiv:1607.06292v2: There
 * should be no Phi function in the 2nd line of (61).
 */
double f_CSd(double xu, double xd, double qu, double qd) noexcept
{
   if (xd < 0.0 || xu < 0.0) {
      ERROR("f_CSd: xu and xd must not be negative!");
      return std::numeric_limits<double>::quiet_NaN();
   } else if (xd == 0.0) {
      return 0.0;
   }

   const double s = 0.25*(qu + qd);
   const double c = sqr(xu - xd) - qu*xu + qd*xd;
   const double cbar = (xu - qu)*xu - (xd + qd)*xd;
   const double lxu = std::log(xu);
   const double lxd = std::log(xd);
   const double phiy = phi_over_y(xu, xd);

   return xd*(-(xu - xd) + (cbar - c*(xu - xd)) * phiy
              + c*(dilog(1.0 - xd/xu) - 0.5*lxu*(lxd - lxu))
              + (s + xd)*lxd + (s - xu)*lxu);
}

/// Eq (62), arxiv:1607.06292, with extra global prefactor xu
double f_CSu(double xu, double xd, double qu, double qd) noexcept
{
   if (xd < 0.0 || xu < 0.0) {
      ERROR("f_CSu: xu and xd must not be negative!");
      return std::numeric_limits<double>::quiet_NaN();
   }

   const double s = 1 + 0.25*(qu + qd);
   const double c = sqr(xu - xd) - (qu + 2)*xu + (qd + 2)*xd;
   const double cbar = (xu - qu - 2)*xu - (xd + qd + 2)*xd;
   const double lxu = std::log(xu);
   const double lxd = std::log(xd);
   const double phiy = phi_over_y(xu, xd);
   const double fCSd = -(xu - xd) + (cbar - c*(xu - xd)) * phiy
      + c*(dilog(1.0 - xd/xu) - 0.5*lxu*(lxd - lxu))
      + (s + xd)*lxd + (s - xu)*lxu;

   return xu*(fCSd - 4.0/3*(xu - xd - 1)*phiy
              - 1.0/3*(lxd + lxu)*(lxd - lxu));
}

/**
 * \f$\mathcal{F}_1(\omega)\f$, Eq (25) arxiv:1502.04199
 */
double F1(double w) noexcept {
   if (w < 0.0) {
      ERROR("F1: w must not be negative!");
      return std::numeric_limits<double>::quiet_NaN();
   } else if (w == 0.0) {
      return 0.0;
   } else if (w == 0.25) {
      return -0.5;
   }

   return (w - 0.5)*f_PS(w) - w*(2 + std::log(w));
}

/**
 * \f$\tilde{\mathcal{F}}_1(\omega)\f$, Eq (26) arxiv:1502.04199
 */
double F1t(double w) noexcept {
   return 0.5*f_PS(w);
}

/**
 * \f$\mathcal{F}_2(\omega)\f$, Eq (27) arxiv:1502.04199
 */
double F2(double w) noexcept {
   if (w < 0.0) {
      ERROR("F2: w must not be negative!");
      return std::numeric_limits<double>::quiet_NaN();
   } else if (w == 0.25) {
      return -0.38629436111989062; // 1 - Log[4]
   }

   return 1 + 0.5*(std::log(w) - f_PS(w));
}

/**
 * \f$\mathcal{F}_3(\omega)\f$, Eq (28) arxiv:1502.04199
 * @author Alexander Voigt
 */
double F3(double w) noexcept {
   if (w < 0.0) {
      ERROR("F3: w must not be negative!");
      return std::numeric_limits<double>::quiet_NaN();
   } else if (w == 0.25) {
      return 19.0/4.0;
   } if (w >= 1e2) {
      const double lw = std::log(w);
      const double iw = 1/w;
      return 89./12 + 42./12*lw + iw*(284./360 + 165./360*lw
         + iw*(6199./44100 + 3885./44100*lw
         + iw*(30017./1.0584e6 + 19530./1.0584e6*lw
         + iw*(83351./1.37214e7 + 55440./1.37214e7*lw
         + iw*(34978051./2.597186592e10 + 23603580./2.597186592e10*lw)))));
   }

   return (0.5 + 7.5*w)*(2 + std::log(w)) + (4.25 - 7.5*w)*f_PS(w);
}

/**
 * Barr-Zee 2-loop function with fermion loop and pseudoscalar and Z
 * boson mediators.
 *
 * @param x squared mass ratio (mf/ms)^2.
 * @param y squared mass ratio (mf/mz)^2.
 */
double FPZ(double x, double y) noexcept
{
   if (x < 0 || y < 0) {
      ERROR("FPZ: arguments must not be negative.");
      return std::numeric_limits<double>::quiet_NaN();
   }

   sort(x, y);

   constexpr double eps = 1e-8;

   if (x == 0 || y == 0) {
      return 0;
   } else if (std::abs(1 - x/y) < eps) {
      if (std::abs(x - 0.25) < eps) {
         // -(1 + 2*Log[2])/3 + O(x - 1/4)
         return -0.79543145370663021 - 1.7453806518612167*(x - 0.25);
      }
      return 2*x*(f_PS(x) + std::log(x))/(1 - 4*x);
   }

   return (y*f_PS(x) - x*f_PS(y))/(x - y);
}

/**
 * Barr-Zee 2-loop function with fermion loop and scalar and Z boson
 * mediators.
 *
 * @param x squared mass ratio (mf/ms)^2.
 * @param y squared mass ratio (mf/mz)^2.
 */
double FSZ(double x, double y) noexcept
{
   if (x < 0 || y < 0) {
      ERROR("FSZ: arguments must not be negative.");
      return std::numeric_limits<double>::quiet_NaN();
   }

   sort(x, y);

   constexpr double eps = 1e-8;

   if (x == 0 || y == 0) {
      return 0;
   } else if (std::abs(1 - x/y) < eps) {
      if (std::abs(x - 0.25) < eps) {
         // (-1 + 4*Log[2])/3 + O(x - 1/4)
         return 0.59086290741326041 + 1.2361419555836500*(x - 0.25);
      } else if (x >= 1e3) {
         const double ix = 1/x;
         const double lx = std::log(x);
         return 7./9 + 2./3*lx
            + ix*(37./150 + 1./5*lx
            + ix*(533./7350 + 2./35*lx
            + ix*(1627./79380 + 1./63*lx
            + ix*(18107./3201660 + 1./231*lx))));
      }
      return 2*x*(1 - 4*x + 2*x*f_PS(x) + std::log(x)*(1 - 2*x))/(4*x - 1);
   }

   return (y*f_S(x) - x*f_S(y))/(x - y);
}

/**
 * Barr-Zee 2-loop function with lepton loop and charge scalar and W
 * boson mediators.
 *
 * @param x squared mass ratio (mf/ms)^2.
 * @param y squared mass ratio (mf/mw)^2.
 */
double FCWl(double x, double y) noexcept
{
   if (x < 0 || y < 0) {
      ERROR("FCWl: arguments must not be negative.");
      return std::numeric_limits<double>::quiet_NaN();
   }

   sort(x, y);

   constexpr double eps = 1e-8;

   if (x == 0 || y == 0) {
      return 0;
   } else if (std::abs(1 - x/y) < eps) {
      const double pi26 = 1.6449340668482264;
      return -f_CSl(x) + x*(-0.5 + x*(3 + (3*x - 2)*(dilog(1 - 1/x) - pi26))
         + (3*x - 0.5)*std::log(x));
   }

   return (y*f_CSl(x) - x*f_CSl(y))/(x - y);
}

/**
 * Barr-Zee 2-loop function with up-type quark loop and charge scalar
 * and W boson mediators.
 *
 * @param xu squared mass ratio (mu/ms)^2.
 * @param xd squared mass ratio (md/ms)^2.
 * @param yu squared mass ratio (mu/mw)^2.
 * @param yd squared mass ratio (md/mw)^2.
 * @param qu electric charge count of up-type quark
 * @param qd electric charge count of down-type quark
 */
double FCWu(double xu, double xd, double yu, double yd, double qu, double qd) noexcept
{
   if (xu < 0 || xd < 0 || yu < 0 || yd < 0) {
      ERROR("FCWu: arguments must not be negative.");
      return std::numeric_limits<double>::quiet_NaN();
   }

   constexpr double eps = 1e-8;

   // Note: xd == yd  <=>  xu == yu, per definition
   if (std::abs(1 - xu/yu) < eps) {
      shift(xu, yu, eps);
      shift(xd, yd, eps);
   }

   return (yu*f_CSu(xu, xd, qu, qd) - xu*f_CSu(yu, yd, qu, qd))/(xu - yu);
}

/**
 * Barr-Zee 2-loop function with down-type quark loop and charge
 * scalar and W boson mediators.
 *
 * @param xu squared mass ratio (mu/ms)^2.
 * @param xd squared mass ratio (md/ms)^2.
 * @param yu squared mass ratio (mu/mw)^2.
 * @param yd squared mass ratio (md/mw)^2.
 * @param qu electric charge count of up-type quark
 * @param qd electric charge count of down-type quark
 */
double FCWd(double xu, double xd, double yu, double yd, double qu, double qd) noexcept
{
   if (xu < 0 || xd < 0 || yu < 0 || yd < 0) {
      ERROR("FCWd: arguments must not be negative.");
      return std::numeric_limits<double>::quiet_NaN();
   }

   constexpr double eps = 1e-8;

   // Note: xd == yd  <=>  xu == yu, per definition
   if (std::abs(1 - xu/yu) < eps) {
      shift(xu, yu, eps);
      shift(xd, yd, eps);
   }

   return (yd*f_CSd(xu, xd, qu, qd) - xd*f_CSd(yu, yd, qu, qd))/(xd - yd);
}

/**
 * Källén lambda function \f$\lambda^2(x,y,z) = x^2 + y^2 + z^2 - 2xy - 2yz - 2xz\f$.
 * The arguments u and v are interpreted as squared masses.
 *
 * @param x squared mass
 * @param y squared mass
 * @param z squared mass
 *
 * @return \f$\lambda^2(x,y,z)\f$
 */
double lambda_2(double x, double y, double z) noexcept
{
   return z*z*lambda_2(x/z, y/z);
}

/**
 * \f$\Phi(x,y,z)\f$ function from arxiv:1607.06292 Eq.(68).

 * @note The arguments x, y and z are interpreted as squared masses.
 *
 * @note Proportional to Phi from Davydychev and Tausk, Nucl. Phys. B397 (1993) 23
 *
 * @param x squared mass
 * @param y squared mass
 * @param z squared mass
 *
 * @return \f$\Phi(x,y,z)\f$
 */
double Phi(double x, double y, double z) noexcept
{
   sort(x, y, z);
   const auto u = x/z, v = y/z;
   return phi_uv(u,v)*z*lambda_2(u, v)/2;
}

} // namespace gm2calc
