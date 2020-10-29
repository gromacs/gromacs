/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
 *
 * \brief
 * CUDA kernels for GPU versions of copy_rvec_to_nbat_real and add_nbat_f_to_f.
 *
 *  \author Alan Gray <alang@nvidia.com>
 *  \author Jon Vincent <jvincent@nvidia.com>
 */

#include "gromacs/gpu_utils/vectype_ops.cuh"
#include "gromacs/nbnxm/nbnxm.h"

/*! \brief CUDA kernel for transforming position coordinates from rvec to nbnxm layout.
 *
 * TODO:
 *  - rename kernel so naming matches with the other NBNXM kernels;
 *  - enable separate compilation unit

 * \param[in]     numColumns          Extent of cell-level parallelism.
 * \param[out]    gm_xq               Coordinates buffer in nbnxm layout.
 * \tparam        setFillerCoords     Whether to set the coordinates of the filler particles.
 * \param[in]     gm_x                Coordinates buffer.
 * \param[in]     gm_atomIndex        Atom index mapping.
 * \param[in]     gm_numAtoms         Array of number of atoms.
 * \param[in]     gm_cellIndex        Array of cell indices.
 * \param[in]     cellOffset          First cell.
 * \param[in]     numAtomsPerCell     Number of atoms per cell.
 */
template<bool setFillerCoords>
static __global__ void nbnxn_gpu_x_to_nbat_x_kernel(int numColumns,
                                                    float4* __restrict__ gm_xq,
                                                    const float3* __restrict__ gm_x,
                                                    const int* __restrict__ gm_atomIndex,
                                                    const int* __restrict__ gm_numAtoms,
                                                    const int* __restrict__ gm_cellIndex,
                                                    int cellOffset,
                                                    int numAtomsPerCell)
{


    const float farAway = -1000000.0f;

    // Map cell-level parallelism to y component of CUDA block index.
    int cxy = blockIdx.y;

    if (cxy < numColumns)
    {

        const int numAtoms = gm_numAtoms[cxy];
        const int offset   = (cellOffset + gm_cellIndex[cxy]) * numAtomsPerCell;
        int       numAtomsRounded;
        if (setFillerCoords)
        {
            // TODO: This can be done more efficiently
            numAtomsRounded = (gm_cellIndex[cxy + 1] - gm_cellIndex[cxy]) * numAtomsPerCell;
        }
        else
        {
            // We fill only the real particle locations.
            // We assume the filling entries at the end have been
            // properly set before during pair-list generation.
            numAtomsRounded = numAtoms;
        }

        const int threadIndex = blockIdx.x * blockDim.x + threadIdx.x;

        // Destination address where x should be stored in nbnxm layout. We use this cast here to
        // save only x, y and z components, not touching the w (q) component, which is pre-defined.
        float3* gm_xqDest = (float3*)&gm_xq[threadIndex + offset];

        // Perform layout conversion of each element.
        if (threadIndex < numAtomsRounded)
        {
            if (threadIndex < numAtoms)
            {
                *gm_xqDest = gm_x[gm_atomIndex[threadIndex + offset]];
            }
            else
            {
                *gm_xqDest = make_float3(farAway);
            }
        }
    }
}
