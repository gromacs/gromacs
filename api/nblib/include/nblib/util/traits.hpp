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
 * Implements general purpose STL-like type traits
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */

#ifndef NBLIB_UTIL_TRAITS_HPP
#define NBLIB_UTIL_TRAITS_HPP

#include <cassert>

#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>


namespace nblib
{

//! \brief Base template for a holder of entries of different data types
template<class... Ts>
struct TypeList
{
};

namespace detail
{
//! \brief unimplemented base template
template<template<class...> class P, class L>
struct [[maybe_unused]] Map_;

/*! \brief Implementation of Map_
 *
 * This is a specialization of the Map_ base template
 * for the case that the L template parameter itself has template parameters
 * in this case, the template parameters of L are caught in Ts...
 *
 */
template<template<class...> class P, template<class...> class L, class... Ts>
struct Map_<P, L<Ts...>>
{
    // resulting type is a TypeList of the P-template instantiated
    // with all template parameters of L
    typedef TypeList<P<Ts>...> type;
};

//! \brief unimplemented base template
template<template<class...> class P, class L>
struct [[maybe_unused]] Reduce_;

//! \brief Implementation of Reduce_
template<template<class...> class P, template<class...> class L, class... Ts>
struct Reduce_<P, L<Ts...>>
{
    typedef P<Ts...> type;
};

//! \brief unimplemented base template
template<class L1, class L2>
struct [[maybe_unused]] FuseTwo_;

//! \brief implementation of FuseTwo_
template<template<class...> class L1, template<class...> class L2, class... Ts1, class... Ts2>
struct FuseTwo_<L1<Ts1...>, L2<Ts2...>>
{
    typedef TypeList<Ts1..., Ts2...> type;
};

//! \brief unimplemented base template
template<class... Ls>
struct [[maybe_unused]] Fuse_;

//! \brief recursion endpoint
template<class L>
struct Fuse_<L>
{
    typedef L type;
};

//! \brief recurse until only one type is left
template<class L1, class L2, class... Ls>
struct Fuse_<L1, L2, Ls...>
{
    typedef typename Fuse_<typename FuseTwo_<L1, L2>::type, Ls...>::type type;
};


//! \brief keep adding the template parameter pack to the type list
template<class L, int N, class... Ts>
struct RepeatHelper_
{
    typedef typename RepeatHelper_<typename FuseTwo_<L, TypeList<Ts...>>::type, N - 1, Ts...>::type type;
};

//! \brief stop recurision
template<class L, class... Ts>
struct RepeatHelper_<L, 1, Ts...>
{
    typedef L type;
};

//! \brief base case
template<class L, int N, class = void>
struct Repeat_
{
};

//! \brief capture original template parameter pack, protect against N < 1
template<template<class...> class L, int N, class... Ts>
struct Repeat_<L<Ts...>, N, std::enable_if_t<N >= 1>>
{
    typedef typename RepeatHelper_<L<Ts...>, N, Ts...>::type type;
};


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

/*! \brief Comparison meta function that compares T to Tuple[N]
 *
 * This trait evaluates to std::true_type if T is the same as Tuple[N]
 * OR if T is the same as the type member of Tuple[N]
 */
template<int N, typename T, typename Tuple>
struct MatchTypeOrTypeMember :
    std::disjunction<std::is_same<T, std::tuple_element_t<N, Tuple>>,
                     std::is_same<T, AccessTypeMemberIfPresent_t<std::tuple_element_t<N, Tuple>>>>
{
};

//! \brief Recursion to check the next field N+1
template<int N, class T, class Tuple, template<int, class, class> class Comparison, class Match = void>
struct MatchField_ : std::integral_constant<size_t, MatchField_<N + 1, T, Tuple, Comparison>{}>
{
};

//! \brief recursion stop when Comparison<N, T, Tuple>::value is true
template<int N, class T, class Tuple, template<int, class, class> class Comparison>
struct MatchField_<N, T, Tuple, Comparison, std::enable_if_t<Comparison<N, T, Tuple>{}>> :
    std::integral_constant<size_t, N>
{
};

} // namespace detail

/*! \brief Create a TypeList of P instantiated with each template parameter of L
 *
 * returns TypeList<P<Ts>...>, with Ts... = template parameters of L
 * does not compile if L has no template parameters
 */
template<template<class...> class P, class L>
using Map = typename detail::Map_<P, L>::type;

/*! \brief Base template for expressing a datatype P templated with all the entries in type list L
 *
 * The result is P instantiated with all the template parameters of L
 */
template<template<class...> class P, class L>
using Reduce = typename detail::Reduce_<P, L>::type;

//! \brief Concatenates template parameters of two variadic templates into a TypeList
template<class... Ls>
using FuseTwo = typename detail::FuseTwo_<Ls...>::type;

/*! \brief This traits concatenates an arbitrary number of variadic templates into a single TypeList
 *
 * For clarity reasons, the fuse operation to fuse two lists into one has been decoupled
 * into a separate trait from the handling of the recursion over the variadic arguments.
 */
template<class... Ls>
using Fuse = typename detail::Fuse_<Ls...>::type;

/*! \brief Repeat the template parameters of L N times
 *
 * L must have template parameters
 * N must be bigger than 0
 * Repeated types are put in a TypeList
 */
template<class L, int N>
using Repeat = typename detail::Repeat_<L, N>::type;

/*! \brief Meta function to return the first index in Tuple whose type matches T
 *
 *  If there are more than one, the first occurrence will be returned.
 *  If there is no such type, the size of Tuple will be returned.
 *  Note that the default comparison operation supplied here also matches if the type member Tuple[N]::type matches T
 */
template<typename T, class TL, template<int, class, class> class Comparison = detail::MatchTypeOrTypeMember>
struct FindIndex
{
};

/*! \brief Specialization to only enable this trait if TL has template parameters
 *
 * \tparam T          a type to look for in the template parameters of TL
 * \tparam TL         a template template parameter, e.g. std::tuple or nblib::TypeList
 * \tparam Ts         template parameters of TL
 * \tparam Comparison comparison operation
 *
 *  Note that \a T is added to \a TL as a sentinel to terminate the recursion
 *  and prevent an out of bounds tuple access compiler error.
 */
template<typename T, template<class...> class TL, class... Ts, template<int, class, class> class Comparison>
struct FindIndex<T, TL<Ts...>, Comparison> : detail::MatchField_<0, T, std::tuple<Ts..., T>, Comparison>
{
};

/*! \brief Meta function to return the element in Tuple whose type matches T
 *
 * If there are more than one, the first occurrence will be returned
 * If there is no such that, a compiler error is generated due to accessing
 * the tuple out of bounds
 */
template<typename T, typename Tuple>
decltype(auto) pickType(Tuple& tup)
{
    return std::get<FindIndex<T, std::decay_t<Tuple>>{}>(tup);
}

//! \brief template meta function to determine whether T is contained in TL
template<class T, class TL>
struct Contains
{
};

/*! \brief implementation of the Contains trait to look for T in TL
 *
 * \tparam T   type to look for in TL
 * \tparam TL  a variadic type, such as std::tuple or TypeList
 * \tparam Ts  the template parameters of TL
 *
 * Note that this clang-format enforced formatting is unfortunate, it should be:
 * struct Contains<T, TL<Ts...>> : std::bool_constant<FindIndex<T, TL<Ts...>>{} < sizeof...(Ts)>
 */
template<class T, template<class...> class TL, class... Ts>
        struct Contains<T, TL<Ts...>> : std::bool_constant < FindIndex<T, TL<Ts...>>{}<sizeof...(Ts)>
{
};

} // namespace nblib

#endif // NBLIB_UTIL_TRAITS_HPP
