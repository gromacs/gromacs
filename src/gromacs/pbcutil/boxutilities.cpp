/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Implements routines in boxutilities.h.
 *
 * Utility functions for handling boxes.
 */
#include "gmxpre.h"

#include "boxutilities.h"

#include "gromacs/math/vectypes.h"

void do_box_rel(int ndim, const matrix deform, matrix box_rel,
                matrix b, bool bInit)
{
    for (int d = YY; d <= ZZ; ++d)
    {
        for (int d2 = XX; d2 < ndim; ++d2)
        {
            /* We need to check if this box component is deformed
             * or if deformation of another component might cause
             * changes in this component due to box corrections.
             */
            if (deform[d][d2] == 0 &&
                !(d == ZZ && d2 == XX && deform[d][YY] != 0 &&
                  (b[YY][d2] != 0 || deform[YY][d2] != 0)))
            {
                if (bInit)
                {
                    box_rel[d][d2] = b[d][d2]/b[XX][XX];
                }
                else
                {
                    b[d][d2] = b[XX][XX]*box_rel[d][d2];
                }
            }
        }
    }
}
