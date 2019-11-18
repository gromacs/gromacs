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
/*! \libinternal \file
 *
 * \brief Basic routines to handle periodic boundary conditions with CUDA.
 *
 * This file contains GPU implementation of the PBC-aware vector evaluation.
 *
 * \todo CPU, GPU and SIMD routines essentially do the same operations on
 *       different data-types. Currently this leads to code duplication,
 *       which has to be resolved. For details, see Redmine task #2863
 *       https://redmine.gromacs.org/issues/2863
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Berk Hess <hess@kth.se>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \inlibraryapi
 * \ingroup module_pbcutil
 */
#ifndef GMX_PBCUTIL_PBC_AIUC_CUDA_CUH
#define GMX_PBCUTIL_PBC_AIUC_CUDA_CUH

#include "gromacs/gpu_utils/gpu_vec.cuh"
#include "gromacs/gpu_utils/vectype_ops.cuh"
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
 * \todo This routine uses CUDA float4 types for input coordinates and
 *       returns in rvec data-type. Other than that, it does essentially
 *       the same thing as the version below, as well as SIMD and CPU
 *       versions. This routine is used in gpubonded module.
 *       To avoid code duplication, these implementations should be
 *       unified. See Redmine task #2863:
 *       https://redmine.gromacs.org/issues/2863
 *
 * \param[in]  pbcAiuc  PBC object.
 * \param[in]  r1       Coordinates of the first point.
 * \param[in]  r2       Coordinates of the second point.
 * \param[out] dr       Resulting distance.
 */
template<bool returnShift>
static __forceinline__ __device__ int
                       pbcDxAiuc(const PbcAiuc& pbcAiuc, const float4& r1, const float4& r2, fvec dr)
{
    dr[XX] = r1.x - r2.x;
    dr[YY] = r1.y - r2.y;
    dr[ZZ] = r1.z - r2.z;

    float shz = rintf(dr[ZZ] * pbcAiuc.invBoxDiagZ);
    dr[XX] -= shz * pbcAiuc.boxZX;
    dr[YY] -= shz * pbcAiuc.boxZY;
    dr[ZZ] -= shz * pbcAiuc.boxZZ;

    float shy = rintf(dr[YY] * pbcAiuc.invBoxDiagY);
    dr[XX] -= shy * pbcAiuc.boxYX;
    dr[YY] -= shy * pbcAiuc.boxYY;

    float shx = rintf(dr[XX] * pbcAiuc.invBoxDiagX);
    dr[XX] -= shx * pbcAiuc.boxXX;

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
 * Computes the vector dr between points r2 and r1, taking into account the
 * periodic boundary conditions, described in pbcAiuc object. Same as above,
 * only takes and returns data in float3 format. Does not return shifts.
 *
 * \todo This routine uses CUDA float3 types for both input and returns
 *       values. Other than that, it does essentially the same thing as the
 *       version above, as well as SIMD and CPU versions. This routine is
 *       used in GPU-based constraints.
 *       To avoid code duplication, these implementations should be
 *       unified. See Redmine task #2863:
 *       https://redmine.gromacs.org/issues/2863
 *
 * \param[in]  pbcAiuc  PBC object.
 * \param[in]  r1       Coordinates of the first point.
 * \param[in]  r2       Coordinates of the second point.
 * \returns    dr       Resulting distance.
 */
static __forceinline__ __host__ __device__ float3 pbcDxAiuc(const PbcAiuc& pbcAiuc,
                                                            const float3&  r1,
                                                            const float3&  r2)
{
    float3 dr = r1 - r2;

    float shz = rintf(dr.z * pbcAiuc.invBoxDiagZ);
    dr.x -= shz * pbcAiuc.boxZX;
    dr.y -= shz * pbcAiuc.boxZY;
    dr.z -= shz * pbcAiuc.boxZZ;

    float shy = rintf(dr.y * pbcAiuc.invBoxDiagY);
    dr.x -= shy * pbcAiuc.boxYX;
    dr.y -= shy * pbcAiuc.boxYY;

    float shx = rintf(dr.x * pbcAiuc.invBoxDiagX);
    dr.x -= shx * pbcAiuc.boxXX;

    return dr;
}

#endif // GMX_PBCUTIL_PBC_AIUC_CUDA_CUH
