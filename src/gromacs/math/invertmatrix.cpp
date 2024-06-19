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
 * Routines to invert 3x3 matrices
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/invertmatrix.h"

#include <cmath>

#include <filesystem>

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

namespace gmx
{

void invertMatrix(const matrix src, matrix dest)
{
    const real smallreal = 1.0e-24_real;
    const real largereal = 1.0e24_real;

    real determinant = det(src);
    real c           = 1.0_real / determinant;
    real fc          = std::fabs(c);

    if ((fc <= smallreal) || (fc >= largereal))
    {
        gmx_fatal(FARGS, "Can not invert matrix, determinant = %e", determinant);
    }
    GMX_ASSERT(dest != src, "Cannot do in-place inversion of matrix");

    dest[XX][XX] = c * (src[YY][YY] * src[ZZ][ZZ] - src[ZZ][YY] * src[YY][ZZ]);
    dest[XX][YY] = -c * (src[XX][YY] * src[ZZ][ZZ] - src[ZZ][YY] * src[XX][ZZ]);
    dest[XX][ZZ] = c * (src[XX][YY] * src[YY][ZZ] - src[YY][YY] * src[XX][ZZ]);
    dest[YY][XX] = -c * (src[YY][XX] * src[ZZ][ZZ] - src[ZZ][XX] * src[YY][ZZ]);
    dest[YY][YY] = c * (src[XX][XX] * src[ZZ][ZZ] - src[ZZ][XX] * src[XX][ZZ]);
    dest[YY][ZZ] = -c * (src[XX][XX] * src[YY][ZZ] - src[YY][XX] * src[XX][ZZ]);
    dest[ZZ][XX] = c * (src[YY][XX] * src[ZZ][YY] - src[ZZ][XX] * src[YY][YY]);
    dest[ZZ][YY] = -c * (src[XX][XX] * src[ZZ][YY] - src[ZZ][XX] * src[XX][YY]);
    dest[ZZ][ZZ] = c * (src[XX][XX] * src[YY][YY] - src[YY][XX] * src[XX][YY]);
}

} // namespace gmx
