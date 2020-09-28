/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
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

#ifndef NBLIB_UTIL_INTERNAL_H
#define NBLIB_UTIL_INTERNAL_H

#include <cassert>

#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include "nblib/basicdefinitions.h"
#include "nblib/vector.h"

namespace nblib
{

namespace detail
{
//! Format strings for use in error messages
std::string next_token(std::string& s, const std::string& delimiter);

// Like std::void_t but for values
template<auto...>
using void_value_t = void;

template<class T, class = void>
struct HasValueMember : std::false_type
{
};

template<class T>
struct HasValueMember<T, void_value_t<T::value>> : std::true_type
{
};

template<class T, class = void>
struct AccessTypeMemberIfPresent
{
    typedef T type;
};

template<class T>
struct AccessTypeMemberIfPresent<T, typename std::void_t<typename T::type>>
{
    typedef typename T::type type;
};

template<class T>
using AccessTypeMemberIfPresent_t = typename AccessTypeMemberIfPresent<T>::type;

//! this trait evaluates to std::true_type if T is the same as Tuple[N]
//! OR if T is the same as the type member of Tuple[N]
template<int N, typename T, typename Tuple>
struct MatchTypeOrTypeMember :
    std::disjunction<std::is_same<T, std::tuple_element_t<N, Tuple>>,
                     std::is_same<T, AccessTypeMemberIfPresent_t<std::tuple_element_t<N, Tuple>>>>
{
};

//! recursion to check the next field N+1
template<int N, class T, class Tuple, template<int, class, class> class Comparison, bool Match = false>
struct MatchField_ :
    std::integral_constant<size_t, MatchField_<N + 1, T, Tuple, Comparison, Comparison<N + 1, T, Tuple>{}>{}>
{
};

//! recursion stop when Comparison<N, T, Tuple>::value is true
template<int N, class T, class Tuple, template<int, class, class> class Comparison>
struct MatchField_<N, T, Tuple, Comparison, true> : std::integral_constant<size_t, N>
{
};

} // namespace detail

/*! \brief The value member of this struct evaluates to the integral constant N for which
 *  the value member of Comparison<N, T, Tuple> is true
 *  and generates a compiler error if there is no such N
 */
template<class T, class Tuple, template<int, class, class> class Comparison>
struct MatchField : detail::MatchField_<0, T, Tuple, Comparison, Comparison<0, T, Tuple>{}>
{
};

/*! \brief Function to return the index in Tuple whose type matches T
 *  - If there are more than one, the first occurrence will be returned
 *  - If there is no such type, a compiler error from accessing a tuple out of range is generated
 *  Note that the default comparison operation supplied here also matches if the type member of Tuple[N] matches T
 */
template<typename T, typename Tuple, template<int, class, class> class Comparison = detail::MatchTypeOrTypeMember>
struct FindIndex : std::integral_constant<size_t, MatchField<T, Tuple, Comparison>{}>
{
};

//! Function to return the element in Tuple whose type matches T
//! Note: if there are more than one, the first occurrence will be returned
template<typename T, typename Tuple>
decltype(auto) pickType(Tuple& tup)
{
    return std::get<FindIndex<T, Tuple>{}>(tup);
}

//! Utility to call function with each element in tuple_
template<class F, class... Ts>
void for_each_tuple(F&& func, std::tuple<Ts...>& tuple_)
{
    std::apply(
            [f = func](auto&... args) {
                [[maybe_unused]] auto list = std::initializer_list<int>{ (f(args), 0)... };
            },
            tuple_);
}

//! Utility to call function with each element in tuple_ with const guarantee
template<class F, class... Ts>
void for_each_tuple(F&& func, const std::tuple<Ts...>& tuple_)
{
    std::apply(
            [f = func](auto&... args) {
                [[maybe_unused]] auto list = std::initializer_list<int>{ (f(args), 0)... };
            },
            tuple_);
}

//! Format strings for use in error messages
template<class... Args>
std::string formatString(std::string fmt, Args... args)
{
    std::ostringstream os;
    std::string        delimiter = "{}";

    std::initializer_list<int> unused{ 0, (os << detail::next_token(fmt, delimiter) << args, 0)... };
    static_cast<void>(unused); // unused is not actually used

    os << detail::next_token(fmt, delimiter);

    return os.str();
}

} // namespace nblib

#endif // NBLIB_UTIL_INTERNAL_H
