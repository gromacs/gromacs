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
/*! \internal \file
 *
 * \brief Implements routines in matrix.h .
 *
 * \author Christian Blau <blau@kth.se>
 */

#include "gmxpre.h"

#include "gromacs/math/matrix.h"

#include "gromacs/math/vectypes.h"
#include "gromacs/mdspan/layouts.h"
#include "gromacs/mdspan/mdspan.h"
#include "gromacs/utility/real.h"

namespace gmx
{


Matrix3x3 transpose(Matrix3x3ConstSpan matrixView)
{

    return Matrix3x3({ matrixView(0, 0),
                       matrixView(1, 0),
                       matrixView(2, 0),
                       matrixView(0, 1),
                       matrixView(1, 1),
                       matrixView(2, 1),
                       matrixView(0, 2),
                       matrixView(1, 2),
                       matrixView(2, 2) });
}

void matrixVectorMultiply(Matrix3x3ConstSpan matrix, RVec* v)
{
    const real resultXX =
            matrix(XX, XX) * (*v)[XX] + matrix(XX, YY) * (*v)[YY] + matrix(XX, ZZ) * (*v)[ZZ];
    const real resultYY =
            matrix(YY, XX) * (*v)[XX] + matrix(YY, YY) * (*v)[YY] + matrix(YY, ZZ) * (*v)[ZZ];
    (*v)[ZZ] = matrix(ZZ, XX) * (*v)[XX] + matrix(ZZ, YY) * (*v)[YY] + matrix(ZZ, ZZ) * (*v)[ZZ];
    (*v)[XX] = resultXX;
    (*v)[YY] = resultYY;
}


} // namespace gmx
