/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2021, by the GROMACS development team, led by
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
 *
 * \brief Basic routines to handle periodic boundary conditions with SYCL.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \inlibraryapi
 * \ingroup module_pbcutil
 */
#ifndef GMX_PBCUTIL_PBC_AIUC_SYCL_H
#define GMX_PBCUTIL_PBC_AIUC_SYCL_H

#include "gromacs/pbcutil/pbc_aiuc.h"

/*! \brief Computes the vector between two points taking PBC into account.
 *
 * Computes the vector dr between points r2 and r1, taking into account the
 * periodic boundary conditions, described in pbcAiuc object. Note that this
 * routine always does the PBC arithmetic for all directions, multiplying the
 * displacements by zeroes if the corresponding direction is not periodic.
 * For triclinic boxes only distances up to half the smallest box diagonal
 * element are guaranteed to be the shortest. This means that distances from
 * 0.5/sqrt(2) times a box vector length (e.g. for a rhombic dodecahedron)
 * can use a more distant periodic image.
 *
 * \todo This routine operates on rvec types and uses PbcAiuc to define
 *       periodic box, but essentially does the same thing as SIMD and GPU
 *       version. These will have to be unified in future to avoid code
 *       duplication. See Issue #2863:
 *       https://gitlab.com/gromacs/gromacs/-/issues/2863
 *
 * \param[in]  pbcAiuc  PBC object.
 * \param[in]  r1       Coordinates of the first point.
 * \param[in]  r2       Coordinates of the second point.
 * \param[out]    dr       Resulting distance.
 */
static inline void pbcDxAiucSycl(const PbcAiuc& pbcAiuc, const rvec& r1, const rvec& r2, rvec dr)
{
    dr[XX] = r1[XX] - r2[XX];
    dr[YY] = r1[YY] - r2[YY];
    dr[ZZ] = r1[ZZ] - r2[ZZ];

    float shz = cl::sycl::rint(dr[ZZ] * pbcAiuc.invBoxDiagZ);
    dr[XX] -= shz * pbcAiuc.boxZX;
    dr[YY] -= shz * pbcAiuc.boxZY;
    dr[ZZ] -= shz * pbcAiuc.boxZZ;

    float shy = cl::sycl::rint(dr[YY] * pbcAiuc.invBoxDiagY);
    dr[XX] -= shy * pbcAiuc.boxYX;
    dr[YY] -= shy * pbcAiuc.boxYY;

    float shx = cl::sycl::rint(dr[XX] * pbcAiuc.invBoxDiagX);
    dr[XX] -= shx * pbcAiuc.boxXX;
}

#endif // GMX_PBCUTIL_PBC_AIUC_SYCL_H
