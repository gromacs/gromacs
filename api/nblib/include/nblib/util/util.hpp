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

#include <limits>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>

#include "nblib/util/annotation.hpp"
#include "nblib/util/tuple.hpp"

namespace util
{

//! \brief like util::integral_constant, but device-qualified
template<typename T, T v>
struct integral_constant
{
    static constexpr T              value = v;
    typedef T                       value_type;
    typedef integral_constant<T, v> type;

    HOST_DEVICE_FUN
    constexpr operator value_type() const noexcept { return value; } // NOLINT
};

template<class... Ts, size_t... Is>
HOST_DEVICE_FUN constexpr auto discardFirstImpl(const tuple<Ts...>& tuple, std::index_sequence<Is...>)
{
    return make_tuple(util::get<Is + 1>(tuple)...);
}

template<class... Ts>
HOST_DEVICE_FUN constexpr auto discardFirstElement(const tuple<Ts...>& tuple)
{
    constexpr int tupleSize = std::tuple_size<util::tuple<Ts...>>::value;
    static_assert(tupleSize > 1);

    using Seq = std::make_index_sequence<tupleSize - 1>;
    return discardFirstImpl(tuple, Seq{});
}

} // namespace util

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
    StrongType() = default;
    //! construction from the underlying type T, implicit conversions disabled
    HOST_DEVICE_FUN explicit StrongType(T v) : value_(std::move(v)) {}

    //! assignment from T
    HOST_DEVICE_FUN StrongType& operator=(T v)
    {
        value_ = std::move(v);
        return *this;
    }

    //! conversion to T
    HOST_DEVICE_FUN operator T() const { return value_; } // NOLINT

    //! access the underlying value
    HOST_DEVICE_FUN const T& value() const { return value_; }
    HOST_DEVICE_FUN T& value() { return value_; }

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
[[maybe_unused]] HOST_DEVICE_FUN inline bool operator==(const StrongType<T, Phantom>& lhs,
                                                        const StrongType<T, Phantom>& rhs)
{
    return lhs.value() == rhs.value();
}

//! comparison function <
template<class T, class Phantom>
HOST_DEVICE_FUN inline bool operator<(const StrongType<T, Phantom>& lhs, const StrongType<T, Phantom>& rhs)
{
    return lhs.value() < rhs.value();
}

//! comparison function >
template<class T, class Phantom>
HOST_DEVICE_FUN inline bool operator>(const StrongType<T, Phantom>& lhs, const StrongType<T, Phantom>& rhs)
{
    return lhs.value() > rhs.value();
}

/*! \brief create a switch on \p runtimeIndex with CompileTimeIndices as the possible cases
 *
 * \tparam CompileTimeIndices  integer sequence of possible switch case values
 * \tparam F                   callable object taking a compile time integer, e.g
 * integral_constant<int, N> \param runtimeIndex         runtime value to call f with \param f
 * callable of type F \return                     f(integral_constant<int, runtimeIndex>{})
 *
 * f is called with the compile time constant value of runtimeIndex.
 * This function (createSwitch) is equivalent to:
 *
 * switch(runtimeIndex)
 * {
 *   case 0: f(integral_constant<int, 0>{}); break;
 *   case 1: f(integral_constant<int, 1>{}); break;
 *   ...
 * }
 */
template<int... CompileTimeIndices, class F>
HOST_DEVICE_FUN auto createSwitch(int runtimeIndex, std::integer_sequence<int, CompileTimeIndices...>, F&& f)
{
    using ReturnType =
            std::common_type_t<decltype(f(util::integral_constant<int, CompileTimeIndices>{}))...>;

    ReturnType                                  ret;
    [[maybe_unused]] std::initializer_list<int> list{ (
            runtimeIndex == CompileTimeIndices
            ? (ret = f(util::integral_constant<int, CompileTimeIndices>{})),
            0
            : 0)... };

    return ret;
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

    auto next_token = [](std::string& s, const std::string& delimiter_) {
        std::string token = s.substr(0, s.find(delimiter_));

        std::size_t next = s.find(delimiter_);
        if (next == std::string::npos)
            s.clear();
        else
            s.erase(0, next + delimiter_.length());

        return token;
    };

    [[maybe_unused]] std::initializer_list<int> unused{ 0, (os << next_token(fmt, delimiter) << args, 0)... };

    os << next_token(fmt, delimiter);

    return os.str();
}

//! wrapper for std::numeric_limits<T> to avoid calling the host std::epsilon function in device code
template<class T>
struct MachineEpsilon
{
    static constexpr T value = std::numeric_limits<T>::epsilon();
};

//! \brief return ceil(dividend/divisor) as integer
template<class I, std::enable_if_t<std::is_integral_v<I>, int> = 0>
I iceil(I dividend, I divisor)
{
    return (dividend + divisor - 1) / divisor;
}

} // namespace nblib

#endif // NBLIB_UTIL_UTIL_HPP
