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
/*! \libinternal \file
 * \brief
 * Declares utilities for template metaprogramming
 *
 * \author Roland Schulz <roland.schulz@intel.com>
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_TEMPLATE_MP_H
#define GMX_UTILITY_TEMPLATE_MP_H

#include <cassert>
#include <cstddef>

#include <utility>

#include "gromacs/compat/mp11.h"

namespace gmx
{

template<class Function>
auto dispatchTemplatedFunction(Function&& f)
{
    return std::forward<Function>(f)();
}

/** \internal \brief Helper function to select appropriate template based on runtime values.
 *
 * Can only use enums for template parameters.
 * These enums must have a member \c Count indicating the total number of valid values.
 *
 * Example usage:
 * \code
    enum class Options {
        Op1 = 0,
        Op2 = 1,
        Count = 2
    };

    template<Options p1, Options p2>
    bool foo(int i);

    bool bar(Options p1, Options p2, int i) {
        return dispatchTemplatedFunction(
            [=](auto p1, auto p2) {
                return foo<p1, p2>(i);
            },
            p1, p2);
    }
 * \endcode
 */
template<class Function, class Enum, class... Enums>
auto dispatchTemplatedFunction(Function&& f, Enum e, Enums... es)
{
    return dispatchTemplatedFunction(
            [&](auto... es_) {
                return compat::mp_with_index<size_t(Enum::Count)>(size_t(e), [&](auto e_) {
                    return std::forward<Function>(f)(
                            std::integral_constant<Enum, static_cast<Enum>(size_t(e_))>(), es_...);
                });
            },
            es...);
}

} // namespace gmx

#endif // GMX_UTILITY_TEMPLATE_MP_H
