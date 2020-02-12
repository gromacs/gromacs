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

#include <tuple>
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

//! Allow creation of vector<tuple<T, T>> from two vector<T>
template<typename T>
std::vector<std::tuple<T, T>> operator+(std::vector<T>&& lhs, std::vector<T>&& rhs)
{

    std::vector<std::tuple<T, T>> ret(lhs.size());

    std::transform(std::make_move_iterator(lhs.cbegin()), std::make_move_iterator(lhs.cend()),
                   std::make_move_iterator(rhs.cbegin()), ret.begin(), [](auto&& lhs_val, auto&& rhs_val) {
                       return std::make_tuple(std::move(lhs_val), std::move(rhs_val));
                   });
    return ret;
}


} // namespace nblib

#endif // GROMACS_UTIL_H
