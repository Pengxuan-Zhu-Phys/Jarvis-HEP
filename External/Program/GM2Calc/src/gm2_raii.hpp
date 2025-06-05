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

#ifndef GM2_RAII_HPP
#define GM2_RAII_HPP

namespace gm2calc {

/**
 * @class RAII_save
 * @brief Saves value of variable and restores it at destruction
 */
template <typename T>
class RAII_save {
public:
   explicit RAII_save(T& var_) noexcept : var(var_), value(var_) {}
   RAII_save(const RAII_save&) = delete;
   RAII_save(RAII_save&&) noexcept = default;
   ~RAII_save() { var = value; }
   RAII_save& operator=(const RAII_save&) = delete;
   RAII_save& operator=(RAII_save&& other) noexcept = default;

private:
   T& var;
   T value{};
};

template <typename T>
constexpr RAII_save<T> make_raii_save(T& var)
{
   return RAII_save<T>(var);
}

} // namespace gm2calc

#endif
