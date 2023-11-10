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
 * The functions in this file are equivalent to the standard library implementations,
 * but qualified for usage in device code. If one would compile device code with
 * the --expt-relaxed-constexpr flag, functions declared as constexpr could be directly
 * used from the standard library.
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */

#ifndef NBLIB_UTIL_FUNCTIONS_HPP
#define NBLIB_UTIL_FUNCTIONS_HPP

#include <cassert>
#include <cmath>

#include "nblib/util/annotation.hpp"

namespace util
{

//! \brief This does what you think it does
template<class T>
HOST_DEVICE_FUN constexpr const T& min(const T& a, const T& b)
{
    if (b < a)
        return b;
    return a;
}

//! \brief This does what you think it does
template<class T>
HOST_DEVICE_FUN constexpr const T& max(const T& a, const T& b)
{
    if (a < b)
        return b;
    return a;
}

//! @brief a simplified version of std::upper_bound that can be compiled as device code
template<class ForwardIt, class T>
HOST_DEVICE_FUN ForwardIt upper_bound(ForwardIt first, ForwardIt last, const T& value)
{
    ForwardIt     it;
    long long int step;
    long long int count = last - first;

    while (count > 0)
    {
        it   = first;
        step = count / 2;
        it += step;
        if (!(value < *it)) // NOLINT
        {
            first = ++it;
            count -= step + 1;
        }
        else
            count = step;
    }
    return first;
}

template<class T>
HOST_DEVICE_FUN constexpr T invsqrt(T x)
{
#ifdef __CUDA_ARCH__
    return rsqrt(x);
#else
    return T(1.0) / std::sqrt(x);
#endif
}

} // namespace util

#endif // NBLIB_UTIL_FUNCTIONS_HPP
