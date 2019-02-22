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
#ifndef GMX_PBCUTIL_PBC_AIUC_CUDA_CUH
#define GMX_PBCUTIL_PBC_AIUC_CUDA_CUH

#include "gromacs/pbcutil/pbc_aiuc.h"

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
template <bool returnShift>
static __forceinline__ __device__
int pbcDxAiuc(const PbcAiuc &pbcAiuc,
              const float4  &r1,
              const float4  &r2,
              fvec           dr)
{
    dr[XX] = r1.x - r2.x;
    dr[YY] = r1.y - r2.y;
    dr[ZZ] = r1.z - r2.z;

    float shz  = rintf(dr[ZZ]*pbcAiuc.invBoxDiagZ);
    dr[XX]    -= shz*pbcAiuc.boxZX;
    dr[YY]    -= shz*pbcAiuc.boxZY;
    dr[ZZ]    -= shz*pbcAiuc.boxZZ;

    float shy  = rintf(dr[YY]*pbcAiuc.invBoxDiagY);
    dr[XX]    -= shy*pbcAiuc.boxYX;
    dr[YY]    -= shy*pbcAiuc.boxYY;

    float shx  = rintf(dr[XX]*pbcAiuc.invBoxDiagX);
    dr[XX]    -= shx*pbcAiuc.boxXX;

    if (returnShift)
    {
        ivec ishift;

        ishift[XX] = -__float2int_rn(shx);
        ishift[YY] = -__float2int_rn(shy);
        ishift[ZZ] = -__float2int_rn(shz);

        return IVEC2IS(ishift);
    }
    else
    {
        return 0;
    }
}

/*! \brief Computes the vector between two points taking PBC into account.
 *
 * Computes the vector dr between points x2 and x1, taking into account the
 * periodic boundary conditions, described in pbcAiuc object. Same as above,
 * only takes and returns data in float3 format. Does not return shifts.
 *
 * \param[in]  pbcAiuc  PBC object.
 * \param[in]  r1       Coordinates of the first point.
 * \param[in]  r2       Coordinates of the second point.
 * \returns    dr       Resulting distance.
 */
static __forceinline__ __host__ __device__
float3 pbcDxAiuc(const PbcAiuc &pbcAiuc,
                 const float3  &r1,
                 const float3  &r2)
{
    float3 dr;
    dr.x = r1.x - r2.x;
    dr.y = r1.y - r2.y;
    dr.z = r1.z - r2.z;

    float shz  = rintf(dr.z*pbcAiuc.invBoxDiagZ);
    dr.x    -= shz*pbcAiuc.boxZX;
    dr.y    -= shz*pbcAiuc.boxZY;
    dr.z    -= shz*pbcAiuc.boxZZ;

    float shy  = rintf(dr.y*pbcAiuc.invBoxDiagY);
    dr.x    -= shy*pbcAiuc.boxYX;
    dr.y    -= shy*pbcAiuc.boxYY;

    float shx  = rintf(dr.x*pbcAiuc.invBoxDiagX);
    dr.x    -= shx*pbcAiuc.boxXX;

    return dr;
}

#endif //GMX_PBCUTIL_PBC_AIUC_CUDA_CUH
