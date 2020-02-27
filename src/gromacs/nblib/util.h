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
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#ifndef GROMACS_UTIL_H
#define GROMACS_UTIL_H

#include <functional>
#include <tuple>
#include <type_traits>
#include <vector>

#include "gromacs/math/vectypes.h"

namespace nblib
{

//! generate Velocites from a Maxwell Boltzmann distro, masses should be the
//! same as the ones specified for the Topology object
std::vector<gmx::RVec> generateVelocity(real Temperature, unsigned int seed, std::vector<real> const& masses);

bool checkNumericValues(const std::vector<gmx::RVec>& values);

template<class T>
inline void ignore_unused(T& x)
{
    static_cast<void>(x);
}

namespace detail
{

template<class T, class U>
struct disjunction : std::integral_constant<bool, T::value || U::value>
{
};

template<class T>
struct void_t
{
    typedef void type;
};

template<class T, class Enable = void>
struct HasTypeMember
{
    typedef T type;
};

template<class T>
struct HasTypeMember<T, typename void_t<typename T::type>::type>
{
    typedef typename T::type type;
};

template<int N, typename T, typename Tuple>
struct CompareField :
    disjunction<std::is_same<T, typename std::tuple_element<N, Tuple>::type>,
                std::is_same<T, typename HasTypeMember<typename std::tuple_element<N, Tuple>::type>::type>>
{
};

template<int N, class T, class Tuple, bool Match = false>
struct MatchingField
{
    static decltype(auto) get(Tuple& tp)
    {
        // check next element
        return MatchingField<N + 1, T, Tuple, CompareField<N + 1, T, Tuple>::value>::get(tp);
    }
};

template<int N, class T, class Tuple>
struct MatchingField<N, T, Tuple, true>
{
    static decltype(auto) get(Tuple& tp) { return std::get<N>(tp); }
};

} // namespace detail

//! Function to return the element in Tuple whose type matches T
//! Note: if there are more than one, the first occurrence will be returned
template<typename T, typename Tuple>
decltype(auto) pickType(Tuple& tup)
{
    return detail::MatchingField<0, T, Tuple, detail::CompareField<0, T, Tuple>::value>::get(tup);
}


} // namespace nblib

#endif // GROMACS_UTIL_H
