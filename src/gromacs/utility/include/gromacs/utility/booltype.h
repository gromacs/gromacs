/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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
/*! \internal \file
 *
 * \brief Define a boolean datatype that can be stored in a std::vector and
 *        have a view on it.
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_utility
 */

#ifndef GMX_BOOLTYPE_H
#define GMX_BOOLTYPE_H

#include <vector>

namespace gmx
{

template<typename>
class ArrayRef;

/*! \brief A clone of a bool as a workaround on the template specialization
 *         of std::vector<bool> that is incompatible with ArrayRef.
 *
 * Use when you need to create an ArrayRef on a vector of boolean values.
 *
 * \note In contrast to bool this type is always initialized to false.
 *
 */
struct BoolType
{
    BoolType() = default;

    /*! \brief Allow implicit construction from plain bool.*/
    BoolType(bool value);

    /*! \brief Conversion to bool. */
    constexpr operator bool() const { return value_; }

    bool value_ = false;
};

/*! \brief
 * Create ArrayRef to bool from reference to std::vector<BoolType>.
 *
 * Allow to easily make views of bool from vectors of BoolType.
 *
 * \see ArrayRef
 */
// NOLINTNEXTLINE(google-runtime-references)
ArrayRef<bool> makeArrayRef(std::vector<BoolType>& boolVector);

/*! \brief
 * Create ArrayRef to const bool from reference to std::vector<BoolType>.
 *
 * Allow to easily make views of const bool from vectors of BoolType.
 *
 * \see ArrayRef
 */
ArrayRef<const bool> makeConstArrayRef(const std::vector<BoolType>& boolVector);

} // namespace gmx
#endif
