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

#ifndef GM2_THDM_PROBLEMS_HPP
#define GM2_THDM_PROBLEMS_HPP

#include <iosfwd>
#include <string>
#include <vector>

namespace gm2calc {

/**
 * @class THDM_problems
 * @brief contains problem and warning flags
 *
 * If a problem has occurred, the physical particle spectrum cannot be
 * trusted (for example tachyons are present).  A warning means an
 * imprecision has occurred and care must be taken when interpreting
 * the particle spectrum.
 */
class THDM_problems {
public:
   void clear();              ///< delete all problems and warnings
   void clear_problems();     ///< delete all problems
   void clear_warnings();     ///< delete all warnings
   void flag_tachyon(const std::string&);
   bool have_tachyon() const; ///< returns true if tachyon exists
   bool have_problem() const; ///< returns true if problem has occurred
   bool have_warning() const; ///< returns true if there is a warning
   std::string get_warnings() const; ///< get warnings as string
   std::string get_problems() const; ///< get problems as string
   void print(std::ostream&) const; ///< print problems and warnings to stream
   void print_problems(std::ostream&) const; ///< print problems to stream
   void print_warnings(std::ostream&) const; ///< print warnings to stream

private:
   std::vector<std::string> tachyons;
};

std::ostream& operator<<(std::ostream&, const THDM_problems&);

} // namespace gm2calc

#endif
