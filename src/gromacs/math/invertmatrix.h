/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * Declares routines to invert 3x3 matrices
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_math
 * \inlibraryapi
 */
#ifndef GMX_MATH_INVERTMATRIX_H
#define GMX_MATH_INVERTMATRIX_H

#include "gromacs/math/vec.h"
#include "gromacs/utility/basedefinitions.h"

namespace gmx
{

/*! \brief Invert a simulation-box matrix in \c src, return in \c dest
 *
 * This routine assumes that src is a simulation-box matrix, i.e. has
 * zeroes in the upper-right triangle. A fatal error occurs if the
 * product of the leading diagonal is too small. The inversion can be
 * done "in place", i.e \c src and \c dest can be the same matrix.
 */
void invertBoxMatrix(const matrix src, matrix dest);

/*! \brief Invert a general 3x3 matrix in \c src, return in \c dest
 *
 * A fatal error occurs if the determinant is too small. \c src and
 * \c dest cannot be the same matrix.
 */
void invertMatrix(const matrix src, matrix dest);

} // namespace gmx

#endif
