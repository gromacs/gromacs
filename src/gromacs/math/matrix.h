/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
/*! \libinternal
 * \file
 * \brief Declares special case of 3x3 matrix frequently used.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_math
 */

#ifndef GMX_MATH_MATRIX_H_
#define GMX_MATH_MATRIX_H_

#include <array>

#include "gromacs/utility/real.h"

#include "multidimarray.h"

namespace gmx
{

/*! \brief Three-by-three matrix of ElementType.
 * \tparam ElementType type of element to be stored in matrix
 */
template <class ElementType>
using BasicMatrix3x3 = MultiDimArray<std::array<ElementType, 3*3>, extents <3, 3> >;

/*! \brief Three-by-three real number matrix.
 * \note will replace the C-style real[3][3] "matrix"
 */
using Matrix3x3 = BasicMatrix3x3<real>;

} // namespace gmx

#endif
