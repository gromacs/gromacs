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

/*! \file
 * \brief  Wrapper around different types of tuples compatible with device code
 *
 * \author Sebastian Keller <sebastian.f.keller@gmail.com>
 */

#ifndef NBLIB_LISTEDFORCES_TUPLE_HPP
#define NBLIB_LISTEDFORCES_TUPLE_HPP

#include <tuple>

#if defined(__CUDACC__)
#    include <thrust/tuple.h>
#endif


#if defined(__CUDACC__)

namespace util
{

template<class... Ts>
using tuple = thrust::tuple<Ts...>;

template<size_t N, class T>
constexpr __host__ __device__ decltype(auto) get(T&& tup) noexcept
{
    return thrust::get<N>(tup);
}

template<class... Ts>
constexpr __host__ __device__ tuple<Ts&...> tie(Ts&... args) noexcept
{
    return thrust::tuple<Ts&...>(args...);
}

template<typename... Ts>
constexpr __host__ __device__ tuple<typename std::decay_t<Ts>...> make_tuple(Ts&&... args)
{
    typedef tuple<typename std::decay_t<Ts>...> TupleType;
    return TupleType(std::forward<Ts>(args)...);
}

} // namespace util

//! \brief specializations of tuple traits in std:: namespace to make structured binding work with thrust tuples
namespace std
{

template<size_t N, class... Ts>
struct tuple_element<N, thrust::tuple<Ts...>>
{
    typedef typename thrust::tuple_element<N, thrust::tuple<Ts...>>::type type;
};

template<class... Ts>
struct tuple_size<thrust::tuple<Ts...>>
{
    static const int value = thrust::tuple_size<thrust::tuple<Ts...>>::value;
};

} // namespace std

#else

namespace util
{

template<class... Ts>
using tuple = std::tuple<Ts...>;


template<size_t N, class T>
constexpr decltype(auto) get(T&& tup) noexcept
{
    return std::get<N>(tup);
}

template<class... Ts>
constexpr tuple<Ts&...> tie(Ts&... args) noexcept
{
    return std::tuple<Ts&...>(args...);
}

template<typename... Ts>
constexpr tuple<typename std::decay_t<Ts>...> make_tuple(Ts&&... args)
{
    typedef tuple<typename std::decay_t<Ts>...> TupleType;
    return TupleType(std::forward<Ts>(args)...);
}

} // namespace util

#endif

#endif // NBLIB_LISTEDFORCES_TUPLE_HPP
