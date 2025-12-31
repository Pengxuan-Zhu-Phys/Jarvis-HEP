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

#ifndef GM2_SLHA_IO_HPP
#define GM2_SLHA_IO_HPP

#include "gm2calc/gm2_error.hpp"

#include "slhaea.h"

#include <cmath>
#include <iosfwd>
#include <string>
#include <type_traits>

#include <Eigen/Core>

namespace gm2calc {

struct Config_options;
class SM;
class MSSMNoFV_onshell;
struct MSSMNoFV_onshell_physical;

namespace thdm {
struct Config;
struct Gauge_basis;
struct Mass_basis;
}

/**
 * @class GM2_slha_io
 * @brief class for reading input files and writing SLHA output files
 */
class GM2_slha_io {
public:
   using Tuple_processor = std::function<void(int, double)>;

   void clear();

   // reading functions
   void read_from_file(const std::string&);
   void read_from_source(const std::string&);
   void read_from_stream(std::istream&);
   void read_block(const std::string&, const Tuple_processor&, double scale = 0) const;
   template <class Derived>
   void read_block(const std::string&, Eigen::MatrixBase<Derived>&, double scale = 0) const;
   double read_scale(const std::string&) const;

   // writing functions
   void write_to_file(const std::string&);
   void write_to_stream(std::ostream&);
   void fill_block_entry(const std::string&, unsigned, double, const std::string&);
   void fill_block_entry(const std::string&, unsigned, const std::string&);

   /// read model parameters (GM2Calc input format)
   void fill_gm2calc(MSSMNoFV_onshell&) const;

   /// read model parameters (SLHA input format)
   void fill_slha(MSSMNoFV_onshell&) const;

   /// read SM parameters
   void fill(SM&) const;

   /// read THDM gauge basis parameters
   void fill(thdm::Gauge_basis&) const;

   /// read THDM mass basis parameters
   void fill(thdm::Mass_basis&) const;

   /// read configuration
   void fill(Config_options&) const;

private:
   SLHAea::Coll data;          ///< SHLA data

   /// convert string to number
   template <class Scalar>
   static Scalar convert_to(const std::string&);
   /// compare block scale
   static bool is_at_scale(const SLHAea::Block&, double, double eps = 0.01);
   /// read scale from block
   static double read_scale(const SLHAea::Block&);
   /// read block with tuple processor
   static void read_block(const SLHAea::Block&, const Tuple_processor&);
   /// read block into Eigen::MatrixBase
   template <class Derived>
   static void read_block(const SLHAea::Block&, Eigen::MatrixBase<Derived>&);
   /// read block into matrix
   template <class Derived>
   static void read_matrix(const SLHAea::Block&, Eigen::MatrixBase<Derived>&);
   /// read block into vector
   template <class Derived>
   static void read_vector(const SLHAea::Block&, Eigen::MatrixBase<Derived>&);

   void fill_scale(MSSMNoFV_onshell&) const;
   void fill_alpha_from_gm2calcinput(MSSMNoFV_onshell&) const;
   void fill_from_A(MSSMNoFV_onshell&) const;
   void fill_from_gm2calcinput(MSSMNoFV_onshell&) const;
   void fill_from_hmix(MSSMNoFV_onshell&) const;
   void fill_from_mass(MSSMNoFV_onshell_physical&) const;
   void fill_from_msoft(MSSMNoFV_onshell&) const;
   void fill_from_sminputs(MSSMNoFV_onshell&) const;
};

template <class Scalar>
Scalar GM2_slha_io::convert_to(const std::string& str)
{
   Scalar value;
   try {
      if (std::is_same<Scalar, int>::value) {
         value = std::stoi(str);
      } else if (std::is_same<Scalar, long>::value) {
         value = std::stol(str);
      } else if (std::is_same<Scalar, long long>::value) {
         value = std::stoll(str);
      } else if (std::is_same<Scalar, unsigned long>::value) {
         value = std::stoul(str);
      } else if (std::is_same<Scalar, unsigned long long>::value) {
         value = std::stoull(str);
      } else if (std::is_same<Scalar, float>::value) {
         value = std::stof(str);
      } else if (std::is_same<Scalar, double>::value) {
         value = std::stod(str);
      } else if (std::is_same<Scalar, long double>::value) {
         value = std::stold(str);
      } else {
         value = SLHAea::to<Scalar>(str);
      }
      if (!std::isfinite(static_cast<double>(value))) {
         throw 1;
      }
   }  catch (...) {
      throw EReadError("non-numeric input");
   }
   return value;
}

/**
 * Fills a matrix from an SLHA block
 *
 * @param block the block
 * @param matrix matrix to be filled
 */
template <class Derived>
void GM2_slha_io::read_matrix(const SLHAea::Block& block, Eigen::MatrixBase<Derived>& matrix)
{
   using Index_t = Eigen::Index;
   static_assert(std::is_signed<Index_t>::value, "Eigen::Index must be a signed integer type.");

   const Index_t cols = matrix.cols(), rows = matrix.rows();

   for (const auto& line : block) {
      if (line.is_data_line() && line.size() >= 3) {
         const Index_t i = convert_to<Index_t>(line[0]) - 1;
         const Index_t k = convert_to<Index_t>(line[1]) - 1;
         if (0 <= i && i < rows && 0 <= k && k < cols) {
            matrix(i, k) = convert_to<double>(line[2]);
         }
      }
   }
}

/**
 * Fills a vector from an SLHA block
 *
 * @param block the block
 * @param vector vector to be filled
 */
template <class Derived>
void GM2_slha_io::read_vector(const SLHAea::Block& block, Eigen::MatrixBase<Derived>& vector)
{
   using Index_t = Eigen::Index;
   static_assert(std::is_signed<Index_t>::value, "Eigen::Index must be a signed integer type.");

   const Index_t rows = vector.rows();

   for (const auto& line : block) {
      if (line.is_data_line() && line.size() >= 2) {
         const Index_t i = convert_to<Index_t>(line[0]) - 1;
         if (0 <= i && i < rows) {
            vector(i) = convert_to<double>(line[1]);
         }
      }
   }
}

/**
 * Fills a matrix from an SLHA block
 *
 * @param block the block
 * @param matrix matrix to be filled
 */
template <class Derived>
void GM2_slha_io::read_block(const SLHAea::Block& block, Eigen::MatrixBase<Derived>& matrix)
{
   if (matrix.cols() == 1) {
      GM2_slha_io::read_vector(block, matrix);
   } else {
      GM2_slha_io::read_matrix(block, matrix);
   }
}

/**
 * Fills a matrix from a SLHA block
 *
 * @param block_name block name
 * @param matrix matrix to be filled
 * @param scale (or 0 if scale should be ignored)
 */
template <class Derived>
void GM2_slha_io::read_block(const std::string& block_name,
                             Eigen::MatrixBase<Derived>& matrix,
                             double scale) const
{
   auto block = data.find(block_name);

   while (block != data.cend()) {
      if (is_at_scale(*block, scale)) {
         GM2_slha_io::read_block(*block, matrix);
      }

      ++block;
      block = data.find(block, data.end(), block_name);
   }
}

} // namespace gm2calc

#endif
