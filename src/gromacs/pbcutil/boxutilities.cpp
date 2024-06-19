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
/*! \internal \file
 * \brief
 * Implements routines in boxutilities.h.
 *
 * Utility functions for handling boxes.
 */
#include "gmxpre.h"

#include "gromacs/pbcutil/boxutilities.h"

#include <cmath>

#include <algorithm>

#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

void do_box_rel(int ndim, const matrix deform, matrix box_rel, matrix b, bool bInit)
{
    for (int d = YY; d <= ZZ; ++d)
    {
        for (int d2 = XX; d2 < ndim; ++d2)
        {
            /* We need to check if this box component is deformed
             * or if deformation of another component might cause
             * changes in this component due to box corrections.
             */
            if (deform[d][d2] == 0
                && !(d == ZZ && d2 == XX && deform[d][YY] != 0 && (b[YY][d2] != 0 || deform[YY][d2] != 0)))
            {
                if (bInit)
                {
                    box_rel[d][d2] = b[d][d2] / b[XX][XX];
                }
                else
                {
                    b[d][d2] = b[XX][XX] * box_rel[d][d2];
                }
            }
        }
    }
}

namespace gmx
{

namespace
{

//! Whether two box elements are equal (with a tolerance).
bool boxElementEqual(real element1, real element2)
{
    // Compare with a relative tolerance (for big boxes) and with
    // an absolute tolerance (small boxes are generally not specified with very
    // high number of decimals).
    return gmx_within_tol(element1, element2, 10 * GMX_REAL_EPS) || std::fabs(element1 - element2) < 1e-3;
}

} // namespace

bool boxesAreEqual(const matrix box1, const matrix box2)
{
    return boxElementEqual(box1[XX][XX], box2[XX][XX]) && boxElementEqual(box1[YY][XX], box2[YY][XX])
           && boxElementEqual(box1[YY][YY], box2[YY][YY])
           && boxElementEqual(box1[ZZ][XX], box2[ZZ][XX]) && boxElementEqual(box1[ZZ][YY], box2[ZZ][YY])
           && boxElementEqual(box1[ZZ][ZZ], box2[ZZ][ZZ]);
}

bool boxIsZero(const matrix box)
{
    return boxElementEqual(box[XX][XX], 0.0) && boxElementEqual(box[YY][XX], 0.0)
           && boxElementEqual(box[YY][YY], 0.0) && boxElementEqual(box[ZZ][XX], 0.0)
           && boxElementEqual(box[ZZ][YY], 0.0) && boxElementEqual(box[ZZ][ZZ], 0.0);
}

} // namespace gmx
