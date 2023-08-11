/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
/*! \file
 *
 * \brief Structure and basic routines to handle periodic boundary conditions.
 *
 * This file contains CPU
 *
 * \todo CPU, GPU and SIMD routines essentially do the same operations on
 *       different data-types. Currently this leads to code duplication,
 *       which has to be resolved. For details, see Issue #2863
 *       https://gitlab.com/gromacs/gromacs/-/issues/2863
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Berk Hess <hess@kth.se>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \inlibraryapi
 * \ingroup module_pbcutil
 */
#ifndef GMX_PBCUTIL_PBC_AIUC_H
#define GMX_PBCUTIL_PBC_AIUC_H

#include "gromacs/pbcutil/ishift.h"

/*! \brief Compact and ordered version of the PBC matrix.
 *
 * The structure contains all the dimensions of the periodic box,
 * arranged so that the memory access pattern is more efficient.
 * This duplicates the information, stored in PBC 'box' matrix object,
 * but without duplicating off-diagonal members of the matrix.
 * The structure can be set by setPbcAiuc( ... ) routine below.
 */
struct PbcAiuc
{
    //! 1/box[ZZ][ZZ]
    float invBoxDiagZ;
    //! box[ZZ][XX]
    float boxZX;
    //! box[ZZ][YY]
    float boxZY;
    //! box[ZZ][ZZ]
    float boxZZ;
    //! 1/box[YY][YY]
    float invBoxDiagY;
    //! box[YY][XX]
    float boxYX;
    //! box[YY][YY]
    float boxYY;
    //! 1/box[XX][XX]
    float invBoxDiagX;
    //! box[XX][XX]
    float boxXX;
};

/*! \brief Set the PBC data-structure.
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
static inline void setPbcAiuc(int numPbcDim, const matrix box, PbcAiuc* pbcAiuc)
{

    pbcAiuc->invBoxDiagZ = 0.0F;
    pbcAiuc->boxZX       = 0.0F;
    pbcAiuc->boxZY       = 0.0F;
    pbcAiuc->boxZZ       = 0.0F;
    pbcAiuc->invBoxDiagY = 0.0F;
    pbcAiuc->boxYX       = 0.0F;
    pbcAiuc->boxYY       = 0.0F;
    pbcAiuc->invBoxDiagX = 0.0F;
    pbcAiuc->boxXX       = 0.0F;

    if (numPbcDim > ZZ)
    {
        pbcAiuc->invBoxDiagZ = 1.0F / box[ZZ][ZZ];
        pbcAiuc->boxZX       = box[ZZ][XX];
        pbcAiuc->boxZY       = box[ZZ][YY];
        pbcAiuc->boxZZ       = box[ZZ][ZZ];
    }
    if (numPbcDim > YY)
    {
        pbcAiuc->invBoxDiagY = 1.0F / box[YY][YY];
        pbcAiuc->boxYX       = box[YY][XX];
        pbcAiuc->boxYY       = box[YY][YY];
    }
    if (numPbcDim > XX)
    {
        pbcAiuc->invBoxDiagX = 1.0F / box[XX][XX];
        pbcAiuc->boxXX       = box[XX][XX];
    }
}

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
 * \param[out] dr       Resulting distance.
 */
static inline void pbcDxAiuc(const PbcAiuc& pbcAiuc, const rvec& r1, const rvec& r2, rvec dr)
{
    dr[XX] = r1[XX] - r2[XX];
    dr[YY] = r1[YY] - r2[YY];
    dr[ZZ] = r1[ZZ] - r2[ZZ];

    float shz = std::rintf(dr[ZZ] * pbcAiuc.invBoxDiagZ);
    dr[XX] -= shz * pbcAiuc.boxZX;
    dr[YY] -= shz * pbcAiuc.boxZY;
    dr[ZZ] -= shz * pbcAiuc.boxZZ;

    float shy = std::rintf(dr[YY] * pbcAiuc.invBoxDiagY);
    dr[XX] -= shy * pbcAiuc.boxYX;
    dr[YY] -= shy * pbcAiuc.boxYY;

    float shx = std::rintf(dr[XX] * pbcAiuc.invBoxDiagX);
    dr[XX] -= shx * pbcAiuc.boxXX;
}


#endif // GMX_PBCUTIL_PBC_AIUC_H
