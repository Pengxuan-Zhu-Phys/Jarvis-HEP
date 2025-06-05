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

#ifndef GM2_EIGEN_UTILS_HPP
#define GM2_EIGEN_UTILS_HPP

#include <algorithm>
#include <limits>
#include <Eigen/Core>

namespace gm2calc {

template <typename Derived>
unsigned closest_index(double mass, const Eigen::ArrayBase<Derived>& v)
{
   unsigned pos = 0;
   typename Derived::PlainObject tmp;
   tmp.setConstant(mass);

   (v - tmp).abs().minCoeff(&pos);

   return pos;
}

template <class Derived>
bool is_equal(const Eigen::ArrayBase<Derived>& a,
              const Eigen::ArrayBase<Derived>& b,
              double precision_goal)
{
   return (a - b).cwiseAbs().maxCoeff() < precision_goal;
}

template <class Derived>
bool is_zero(const Eigen::ArrayBase<Derived>& a, double eps)
{
   return a.cwiseAbs().maxCoeff() < eps;
}

/**
 * Normalize each element of the given real matrix to be within the
 * interval [min, max].  Values < min are set to min.  Values > max
 * are set to max.
 *
 * @param m matrix
 * @param min minimum
 * @param max maximum
 */
template <int M, int N>
void normalize_to_interval(Eigen::Matrix<double,M,N>& m, double min = -1., double max = 1.)
{
   auto data = m.data();
   const auto size = m.size();

   for (int i = 0; i < size; i++) {
      if (data[i] < min) {
         data[i] = min;
      } else if (data[i] > max) {
         data[i] = max;
      }
   }
}

/**
 * The element of v, which is closest to mass, is moved to the
 * position idx.
 *
 * @param idx new index of the mass eigenvalue
 * @param mass mass to compare against
 * @param v vector of masses
 * @param z corresponding mixing matrix
 */
template <typename DerivedArray, typename DerivedMatrix>
void move_goldstone_to(int idx, double mass, Eigen::ArrayBase<DerivedArray>& v,
                       Eigen::MatrixBase<DerivedMatrix>& z)
{
   int pos = closest_index(mass, v);
   if (pos == idx) {
      return;
   }

   const int sign = (idx - pos) < 0 ? -1 : 1;
   int steps = std::abs(idx - pos);

   // now we shuffle the states
   while (steps--) {
      const int new_pos = pos + sign;
      v.row(new_pos).swap(v.row(pos));
      z.row(new_pos).swap(z.row(pos));
      pos = new_pos;
   }
}

/**
 * Returns all elements from src, which are not close to the elements
 * in cmp.  The returned vector will have the length (src.size() -
 * cmp.size()).
 *
 * @param src source vector
 * @param cmp vector with elements to compare against
 * @return vector with elements of src not close to cmp
 */
template<class Real, int Nsrc, int Ncmp>
Eigen::Array<Real,Nsrc - Ncmp,1> remove_if_equal(
   const Eigen::Array<Real,Nsrc,1>& src,
   const Eigen::Array<Real,Ncmp,1>& cmp)
{
   static_assert(Nsrc > Ncmp, "Error: size of source vector is not greater "
                              "than size of comparison vector.");

   Eigen::Array<Real,Nsrc,1> non_equal(src);
   Eigen::Array<Real,Nsrc - Ncmp,1> dst;

   for (int i = 0; i < Ncmp; i++) {
      const int idx = closest_index(cmp(i), non_equal);
      non_equal(idx) = std::numeric_limits<double>::infinity();
   }

   std::remove_copy_if(non_equal.data(), non_equal.data() + Nsrc,
                       dst.data(), [] (auto x) { return !std::isfinite(x); });

   return dst;
}

/**
 * @brief reorders vector v according to ordering in vector v2
 * @param v vector with elementes to be reordered
 * @param v2 vector with reference ordering
 */
template<class Real, int N>
void reorder_vector(
   Eigen::Array<Real,N,1>& v,
   const Eigen::Array<Real,N,1>& v2)
{
   Eigen::PermutationMatrix<N> p;
   p.setIdentity();
   std::sort(p.indices().data(), p.indices().data() + p.indices().size(),
             [&v2] (int i, int j) { return std::abs(v2[i]) < std::abs(v2[j]); });

#if EIGEN_VERSION_AT_LEAST(3,1,4)
   v.matrix().transpose() *= p.inverse();
#else
   Eigen::Map<Eigen::Matrix<Real,N,1> >(v.data()).transpose() *= p.inverse();
#endif
}

/**
 * @brief reorders vector v according to ordering of diagonal elements in mass_matrix
 * @param v vector with elementes to be reordered
 * @param matrix matrix with diagonal elements with reference ordering
 */
template<class Derived>
void reorder_vector(
   Eigen::Array<double,Eigen::MatrixBase<Derived>::RowsAtCompileTime,1>& v,
   const Eigen::MatrixBase<Derived>& matrix)
{
   reorder_vector(v, matrix.diagonal().array().eval());
}

template <typename Derived>
void symmetrize(Eigen::MatrixBase<Derived>& m)
{
   static_assert(Eigen::MatrixBase<Derived>::RowsAtCompileTime ==
                 Eigen::MatrixBase<Derived>::ColsAtCompileTime,
                 "symmetrize is only defined for squared matrices");

   for (int i = 0; i < Eigen::MatrixBase<Derived>::RowsAtCompileTime; i++) {
      for (int k = 0; k < i; k++) {
         m(i,k) = m(k,i);
      }
   }
}

} // namespace gm2calc

#endif
