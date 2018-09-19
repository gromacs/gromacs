/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
#ifndef GMX_MATH_DATAVECTORMEASURE_H
#define GMX_MATH_DATAVECTORMEASURE_H

#include <algorithm>

#include "gromacs/math/container/containermeasure.h"
#include "gromacs/utility/real.h"

namespace gmx
{

namespace containeroperation
{

template <typename T>
void multiply(T * vecData, typename T::value_type factor)
{
    std::for_each(std::begin(*vecData), std::end(*vecData), [factor](typename T::value_type &v) { v *= factor; });
}

template <typename T>
void divide(T * vecData, typename T::value_type divisor)
{
    std::for_each(std::begin(*vecData), std::end(*vecData), [divisor](typename T::value_type &v) { v /= divisor; });
}

template <typename T>
void addScalar(T * vecData, typename T::value_type summand)
{
    std::for_each(std::begin(*vecData), std::end(*vecData), [summand](typename T::value_type &datum) { datum += summand; });
}

template <typename T, typename U>
void add(T * vecData, const U &summand)
{
    std::transform(std::begin(*vecData), std::end(*vecData), std::begin(summand), std::begin(*vecData), [](typename T::value_type a, typename U::value_type b) { return a+b; });
}

template <typename T>
void normalize(T * vecData)
{
    divide(vecData, containermeasure::mean(*vecData));
}


template <typename T>
void setZero(T * vecData)
{
    std::fill(vecData->begin(), vecData->end(), 0);
}

} // namespace containeroperation
} // namespace gmx

#endif
