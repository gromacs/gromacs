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
#ifndef GMX_PBCUTIL_PBC_AIUC_H
#define GMX_PBCUTIL_PBC_AIUC_H

#include "gromacs/pbcutil/ishift.h"

/*! \brief Compact and ordered version of the PBC matrix.
 *
 * The structure contains all the dimensions of the periodic box,
 * arranged in a convenient order. This duplicates the information,
 * stored in PBC 'box' matrix object. The structure can be set by
 * setPbcAiuc( ... ) routine below.
 *
 */
struct PbcAiuc
{
    float invBoxDiagZ;
    float boxZX;
    float boxZY;
    float boxZZ;
    float invBoxDiagY;
    float boxYX;
    float boxYY;
    float invBoxDiagX;
    float boxXX;
};

/*! \brief Set the PBC data to use in GPU kernels.
 *
 * \param[in]   numPbcDim   Number of periodic dimensions:
 *                          0 - no periodicity.
 *                          1 - periodicity along X-axis.
 *                          2 - periodicity in XY plane.
 *                          3 - periodicity along all dimensions.
 * \param[in]   box         Matrix, describing the periodic cell.
 * \param[out]  pbcAiuc     Pointer to PbcAiuc structure to be initialized.
 *
 */
static void setPbcAiuc(int           numPbcDim,
                       const matrix  box,
                       PbcAiuc      *pbcAiuc)
{

    pbcAiuc->invBoxDiagZ = 0.0f;
    pbcAiuc->boxZX       = 0.0f;
    pbcAiuc->boxZY       = 0.0f;
    pbcAiuc->boxZZ       = 0.0f;
    pbcAiuc->invBoxDiagY = 0.0f;
    pbcAiuc->boxYX       = 0.0f;
    pbcAiuc->boxYY       = 0.0f;
    pbcAiuc->invBoxDiagX = 0.0f;
    pbcAiuc->boxXX       = 0.0f;

    if (numPbcDim > ZZ)
    {
        pbcAiuc->invBoxDiagZ = 1.0f/box[ZZ][ZZ];
        pbcAiuc->boxZX       = box[ZZ][XX];
        pbcAiuc->boxZY       = box[ZZ][YY];
        pbcAiuc->boxZZ       = box[ZZ][ZZ];
    }
    if (numPbcDim > YY)
    {
        pbcAiuc->invBoxDiagY = 1.0f/box[YY][YY];
        pbcAiuc->boxYX       = box[YY][XX];
        pbcAiuc->boxYY       = box[YY][YY];
    }
    if (numPbcDim > XX)
    {
        pbcAiuc->invBoxDiagX = 1.0f/box[XX][XX];
        pbcAiuc->boxXX       = box[XX][XX];
    }
}

/*! \brief Computes the vector between two points taking PBC into account.
 *
 * Computes the vector dr between points x2 and x1, taking into account the
 * periodic boundary conditions, described in pbcAiuc object. Note that this
 * routine always does the PBC arithmetic for all directions, multiplying the
 * displacements by zeroes if the corresponding direction is not periodic.
 * For triclinic boxes only distances up to half the smallest box diagonal
 * element are guaranteed to be the shortest. This means that distances from
 * 0.5/sqrt(2) times a box vector length (e.g. for a rhombic dodecahedron)
 * can use a more distant periodic image.
 *
 * \param[in]  pbcAiuc  PBC object.
 * \param[in]  r1       Coordinates of the first point.
 * \param[in]  r2       Coordinates of the second point.
 * \param[out] dr       Resulting distance.
 */
static inline
void pbcDxAiuc(const PbcAiuc  &pbcAiuc,
               const rvec     &r1,
               const rvec     &r2,
               rvec            dr)
{
    dr[XX] = r1[XX] - r2[XX];
    dr[YY] = r1[YY] - r2[YY];
    dr[ZZ] = r1[ZZ] - r2[ZZ];

    float shz  = rintf(dr[ZZ]*pbcAiuc.invBoxDiagZ);
    dr[XX]    -= shz*pbcAiuc.boxZX;
    dr[YY]    -= shz*pbcAiuc.boxZY;
    dr[ZZ]    -= shz*pbcAiuc.boxZZ;

    float shy  = rintf(dr[YY]*pbcAiuc.invBoxDiagY);
    dr[XX]    -= shy*pbcAiuc.boxYX;
    dr[YY]    -= shy*pbcAiuc.boxYY;

    float shx  = rintf(dr[XX]*pbcAiuc.invBoxDiagX);
    dr[XX]    -= shx*pbcAiuc.boxXX;
}


#endif //GMX_PBCUTIL_PBC_AIUC_H
