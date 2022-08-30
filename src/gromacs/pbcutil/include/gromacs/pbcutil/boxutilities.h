/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
#ifndef GMX_PBCUTIL_BOXUTILITIES_H
#define GMX_PBCUTIL_BOXUTILITIES_H

#include <stdio.h>

#include "gromacs/math/vectypes.h"

/*! \brief Change box components to preserve the relative box shape
 *
 * Typically used with bInit set to false, in which case it changes
 * box components to b[XX][XX]*box_rel to preserve the relative box
 * shape. If bInit is true, then the values in b are used to set
 * those in box_rel so that subsquent calls can use that box_rel to
 * adjust b to maintain a consistent box.
 */
void do_box_rel(int ndim, const matrix deform, matrix box_rel, matrix b, bool bInit);

namespace gmx
{

/*! \brief
 * Returns whether two boxes are of equal size and shape (within reasonable
 * tolerance).
 */
bool boxesAreEqual(const matrix box1, const matrix box2);

/*! \brief
 * Returns whether a box is only initialised to zero or not.
 */
bool boxIsZero(const matrix box);

} // namespace gmx

#endif
