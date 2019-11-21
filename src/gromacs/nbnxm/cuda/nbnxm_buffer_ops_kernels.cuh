/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
 *  - improve/simplify/document use of cxy_na and na_round
 *  - rename kernel so naming matches with the other NBNXM kernels;
 *  - enable separate compilation unit

 * \param[in]     numColumns          extent of cell-level parallelism
 * \param[out]    gm_coordinatesNbnxm coordinates buffer in nbnxm layout
 * \param[in]     setFillerCoords     tells whether to set the coordinates of the filler particles
 * \param[in]     gm_coordinatesRvec  coordinates buffer in rvec format
 * \param[in]     gm_atomIndex        atom index mapping
 * \param[in]     gm_numAtoms         array of number of atoms
 * \param[in]     gm_cellIndex        array of cell indices
 * \param[in]     cellOffset          first cell
 * \param[in]     numAtomsPerCell     number of atoms per cell
 */
__global__ void nbnxn_gpu_x_to_nbat_x_kernel(int numColumns,
                                             float* __restrict__ gm_coordinatesNbnxm,
                                             bool setFillerCoords,
                                             const rvec* __restrict__ gm_coordinatesRvec,
                                             const int* __restrict__ gm_atomIndex,
                                             const int* __restrict__ gm_numAtoms,
                                             const int* __restrict__ gm_cellIndex,
                                             int cellOffset,
                                             int numAtomsPerCell);


__global__ void nbnxn_gpu_x_to_nbat_x_kernel(int numColumns,
                                             float* __restrict__ gm_coordinatesNbnxm,
                                             bool setFillerCoords,
                                             const rvec* __restrict__ gm_coordinatesRvec,
                                             const int* __restrict__ gm_atomIndex,
                                             const int* __restrict__ gm_numAtoms,
                                             const int* __restrict__ gm_cellIndex,
                                             int cellOffset,
                                             int numAtomsPerCell)
{


    const float farAway = -1000000.0f;

    /* map cell-level parallelism to y component of CUDA block index */
    int cxy = blockIdx.y;

    if (cxy < numColumns)
    {

        int na = gm_numAtoms[cxy];
        int a0 = (cellOffset + gm_cellIndex[cxy]) * numAtomsPerCell;
        int na_round;
        if (setFillerCoords)
        {
            // TODO: This can be done more efficiently
            na_round = (gm_cellIndex[cxy + 1] - gm_cellIndex[cxy]) * numAtomsPerCell;
        }
        else
        {
            /* We fill only the real particle locations.
             * We assume the filling entries at the end have been
             * properly set before during pair-list generation.
             */
            na_round = na;
        }

        /* map parallelism within a cell to x component of CUDA block index linearized
         * with threads within a block */
        int i, j0;
        i = blockIdx.x * blockDim.x + threadIdx.x;

        j0 = a0 * STRIDE_XYZQ;

        // destination address where x shoud be stored in nbnxm layout
        float3* gm_coordinatesDest = (float3*)&gm_coordinatesNbnxm[j0 + 4 * i];

        /* perform conversion of each element */
        if (i < na_round)
        {
            if (i < na)
            {
                *gm_coordinatesDest = *((float3*)gm_coordinatesRvec[gm_atomIndex[a0 + i]]);
            }
            else
            {
                *gm_coordinatesDest = make_float3(farAway);
            }
        }
    }
}

/*! \brief CUDA kernel to sum up the force components
 *
 * \tparam        accumulateForce  If the initial forces in \p gm_fTotal should be saved.
 * \tparam        addPmeForce      Whether the PME force should be added to the total.
 *
 * \param[in]     gm_forcesNbnxm   Non-bonded forces in nbnxm format.
 * \param[in]     gm_forcesPme     PME forces.
 * \param[in,out] gm_forcesTotal   Force buffer to be reduced into.
 * \param[in]     cell             Cell index mapping.
 * \param[in]     atomStart        Start atom index.
 * \param[in]     numAtoms         Number of atoms.
 */
template<bool accumulateForce, bool addPmeForce>
__global__ void nbnxn_gpu_add_nbat_f_to_f_kernel(const float3* __restrict__ gb_forcesNbnxm,
                                                 const float3* __restrict__ gm_forcesPme,
                                                 float3* gm_forcesTotal,
                                                 const int* __restrict__ gm_cell,
                                                 const int atomStart,
                                                 const int numAtoms);
template<bool accumulateForce, bool addPmeForce>
__global__ void nbnxn_gpu_add_nbat_f_to_f_kernel(const float3* __restrict__ gb_forcesNbnxm,
                                                 const float3* __restrict__ gm_forcesPme,
                                                 float3* gm_forcesTotal,
                                                 const int* __restrict__ gm_cell,
                                                 const int atomStart,
                                                 const int numAtoms)
{

    /* map particle-level parallelism to 1D CUDA thread and block index */
    int threadIndex = blockIdx.x * blockDim.x + threadIdx.x;

    /* perform addition for each particle*/
    if (threadIndex < numAtoms)
    {

        int     i             = gm_cell[atomStart + threadIndex];
        float3* gm_forcesDest = (float3*)&gm_forcesTotal[atomStart + threadIndex];
        float3  temp;

        if (accumulateForce)
        {
            temp = *gm_forcesDest;
            temp += gb_forcesNbnxm[i];
        }
        else
        {
            temp = gb_forcesNbnxm[i];
        }
        if (addPmeForce)
        {
            temp += gm_forcesPme[atomStart + threadIndex];
        }
        *gm_forcesDest = temp;
    }
    return;
}
