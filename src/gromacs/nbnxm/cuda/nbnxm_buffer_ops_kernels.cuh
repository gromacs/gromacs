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
 * CUDA kernel for GPU version of copy_rvec_to_nbat_real.
 * Converts coordinate data from rvec to nb format.
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
 * \param[in]     gridIndex         grid index
 * \param[in]     FillLocal         boolean to specify if Fill Local is true
 * \param[in]     x                 position buffer
 * \param[in]     a                 atom index mapping stride between atoms in memory
 * \param[in]     cxy_na            array of extents
 * \param[in]     cxy_ind           array of cell indices
 * \param[in]     cellOffset        first cell
 * \param[in]     numAtomsPerCell   number of atoms per cell
 */
__global__ void nbnxn_gpu_x_to_nbat_x_kernel(int                         numColumns,
                                             float *  __restrict__       xnb,
                                             int                         gridIndex,
                                             bool                        FillLocal,
                                             const rvec *  __restrict__  x,
                                             const int *  __restrict__   a,
                                             const int *  __restrict__   cxy_na,
                                             const int *  __restrict__   cxy_ind,
                                             int                         cellOffset,
                                             int                         numAtomsPerCell);


__global__ void nbnxn_gpu_x_to_nbat_x_kernel(int                         numColumns,
                                             float *  __restrict__       xnb,
                                             int                         gridIndex,
                                             bool                        FillLocal,
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
        if (gridIndex == 0 && FillLocal)
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
