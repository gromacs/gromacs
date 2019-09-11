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

 * \param[in]     numColumns        extent of cell-level parallelism
 * \param[out]    xnb               position buffer in nbnxm layout
 * \param[in]     setFillerCoords   tells whether to set the coordinates of the filler particles
 * \param[in]     x                 position buffer
 * \param[in]     a                 atom index mapping stride between atoms in memory
 * \param[in]     cxy_na            array of extents
 * \param[in]     cxy_ind           array of cell indices
 * \param[in]     cellOffset        first cell
 * \param[in]     numAtomsPerCell   number of atoms per cell
 */
__global__ void nbnxn_gpu_x_to_nbat_x_kernel(int                         numColumns,
                                             float *  __restrict__       xnb,
                                             bool                        setFillerCoords,
                                             const rvec *  __restrict__  x,
                                             const int *  __restrict__   a,
                                             const int *  __restrict__   cxy_na,
                                             const int *  __restrict__   cxy_ind,
                                             int                         cellOffset,
                                             int                         numAtomsPerCell);


__global__ void nbnxn_gpu_x_to_nbat_x_kernel(int                         numColumns,
                                             float *  __restrict__       xnb,
                                             bool                        setFillerCoords,
                                             const rvec *  __restrict__  x,
                                             const int *  __restrict__   a,
                                             const int *  __restrict__   cxy_na,
                                             const int *  __restrict__   cxy_ind,
                                             int                         cellOffset,
                                             int                         numAtomsPerCell)
{


    const float farAway = -1000000.0f;

    /* map cell-level parallelism to y component of CUDA block index */
    int cxy = blockIdx.y;

    if (cxy < numColumns)
    {

        int na = cxy_na[cxy];
        int a0 = (cellOffset + cxy_ind[cxy])*numAtomsPerCell;
        int na_round;
        if (setFillerCoords)
        {
            // TODO: This can be done more efficiently
            na_round =
                (cxy_ind[cxy+1] - cxy_ind[cxy])*numAtomsPerCell;
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
        i = blockIdx.x*blockDim.x+threadIdx.x;

        j0 = a0*STRIDE_XYZQ;

        // destination address where x shoud be stored in nbnxm layout
        float3 *x_dest = (float3 *)&xnb[j0 + 4*i];

        /* perform conversion of each element */
        if (i < na_round)
        {
            if (i < na)
            {
                *x_dest = *((float3 *)x[a[a0 + i]]);
            }
            else
            {
                *x_dest = make_float3(farAway);
            }
        }
    }

}

/*! \brief CUDA kernel to sum up the force components
 *
 * \tparam        accumulateForce  If the initial forces in \p d_fTotal should be saved.
 * \tparam        addPmeForce      Whether the PME force should be added to the total.
 *
 * \param[in]     d_fNB            Non-bonded forces in nbat format.
 * \param[in]     d_fPme           PME forces.
 * \param[in,out] d_fTotal         Force buffer to be reduced into.
 * \param[in]     cell             Cell index mapping.
 * \param[in]     atomStart        Start atom index.
 * \param[in]     numAtoms         Number of atoms.
 */
template <bool accumulateForce, bool addPmeForce>
__global__ void
nbnxn_gpu_add_nbat_f_to_f_kernel(const float3 *__restrict__  d_fNB,
                                 const float3 *__restrict__  d_fPme,
                                 float3                     *d_fTotal,
                                 const int *__restrict__     d_cell,
                                 const int                   atomStart,
                                 const int                   numAtoms);
template <bool accumulateForce, bool addPmeForce>
__global__ void
nbnxn_gpu_add_nbat_f_to_f_kernel(const float3 *__restrict__  d_fNB,
                                 const float3 *__restrict__  d_fPme,
                                 float3                     *d_fTotal,
                                 const int *__restrict__     d_cell,
                                 const int                   atomStart,
                                 const int                   numAtoms)
{

    /* map particle-level parallelism to 1D CUDA thread and block index */
    int threadIndex = blockIdx.x*blockDim.x+threadIdx.x;

    /* perform addition for each particle*/
    if (threadIndex < numAtoms)
    {

        int     i        = d_cell[atomStart+threadIndex];
        float3 *fDest    = (float3 *)&d_fTotal[atomStart+threadIndex];
        float3  temp;

        if (accumulateForce)
        {
            temp  = *fDest;
            temp += d_fNB[i];
        }
        else
        {
            temp = d_fNB[i];
        }
        if (addPmeForce)
        {
            temp += d_fPme[atomStart+threadIndex];
        }
        *fDest = temp;

    }
    return;
}
