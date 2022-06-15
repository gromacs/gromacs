/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \inpublicapi \file
 * \brief
 * Implements nblib utilities
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */

#ifndef NBLIB_UTIL_UTIL_HPP
#define NBLIB_UTIL_UTIL_HPP

#include <cassert>

#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

namespace nblib
{

/*! \brief A template to create structs as a type-safe alternative to using declarations
 *
 * \inpublicapi
 * \ingroup nblib
 *
 * Used in public API functions where a distinction between different
 * arguments of the same underlying type is desired. This provides a type-safe
 * version to using declarations. Instead of naming a type alias, the name
 * is used to define a struct that inherits from StrongType<T>, where T is
 * the underlying type. For example:
 *
 * struct C6 : StrongType<real>
 * {
 *     using StrongType::StrongType;
 *     using StrongType::operator=;
 * };
 *
 * Due to the T() conversion and assignment from T,
 * an instance of the resulting C6 struct behaves essentially like a real, while construction
 * from real is disabled. This makes it impossible to pass a real as a function parameter
 * of type C6.
 */
template<class T, class Phantom>
struct StrongType
{
    //! default ctor
    StrongType() : value_{} {}
    //! construction from the underlying type T, implicit conversions disabled
    explicit StrongType(T v) : value_(std::move(v)) {}

    //! assignment from T
    StrongType& operator=(T v)
    {
        value_ = std::move(v);
        return *this;
    }

    //! conversion to T
    operator T() const { return value_; }

    //! access the underlying value
    T value() const { return value_; }

private:
    T value_;
};

/*! \brief StrongType equality comparison
 *
 * Requires that both T and Phantom template parameters match.
 * For the case where a comparison between StrongTypes with matching T, but differing Phantom
 * parameters is desired, the underlying value attribute should be compared instead
 */
template<class T, class Phantom>
[[maybe_unused]] inline bool operator==(const StrongType<T, Phantom>& lhs, const StrongType<T, Phantom>& rhs)
{
    return lhs.value() == rhs.value();
}

//! comparison function <
template<class T, class Phantom>
inline bool operator<(const StrongType<T, Phantom>& lhs, const StrongType<T, Phantom>& rhs)
{
    return lhs.value() < rhs.value();
}

//! comparison function >
template<class T, class Phantom>
inline bool operator>(const StrongType<T, Phantom>& lhs, const StrongType<T, Phantom>& rhs)
{
    return lhs.value() > rhs.value();
}

//! \brief Utility to call function with each element in tuple_
template<class F, class... Ts>
void for_each_tuple(F&& func, std::tuple<Ts...>& tuple_)
{
    std::apply(
            [f = func](auto&... args) {
                [[maybe_unused]] auto list = std::initializer_list<int>{ (f(args), 0)... };
            },
            tuple_);
}

//! \brief Utility to call function with each element in tuple_ with const guarantee
template<class F, class... Ts>
void for_each_tuple(F&& func, const std::tuple<Ts...>& tuple_)
{
    std::apply(
            [f = func](auto&... args) {
                [[maybe_unused]] auto list = std::initializer_list<int>{ (f(args), 0)... };
            },
            tuple_);
}

//! \brief Format strings for use in error messages
template<class... Args>
std::string formatString(std::string fmt, Args... args)
{
    std::ostringstream os;
    std::string        delimiter = "{}";

    auto next_token = [](std::string& s, const std::string& delimiter_)
    {
        std::string token = s.substr(0, s.find(delimiter_));

        std::size_t next = s.find(delimiter_);
        if (next == std::string::npos)
            s.clear();
        else
            s.erase(0, next + delimiter_.length());

        return token;
    };

    [[maybe_unused]]
    std::initializer_list<int> unused{ 0, (os << next_token(fmt, delimiter) << args, 0)... };

    os << next_token(fmt, delimiter);

    return os.str();
}

} // namespace nblib

#endif // NBLIB_UTIL_UTIL_HPP
