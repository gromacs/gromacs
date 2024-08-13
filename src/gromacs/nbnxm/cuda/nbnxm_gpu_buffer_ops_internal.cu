/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
 *  \brief Define CUDA kernel (and its wrapper) for transforming position coordinates from rvec to nbnxm layout.
 *
 *  \author Alan Gray <alang@nvidia.com>
 *  \author Jon Vincent <jvincent@nvidia.com>
 *  \author Szilard Pall <pall.szilard@gmail.com>
 */

#include "gmxpre.h"

#include "gromacs/nbnxm/nbnxm_gpu_buffer_ops_internal.h"

#include "gromacs/gpu_utils/typecasts_cuda_hip.h"
#include "gromacs/gpu_utils/vectype_ops.cuh"
#include "gromacs/nbnxm/cuda/nbnxm_cuda_types.h"
#include "gromacs/nbnxm/grid.h"

namespace gmx
{

/*! \brief CUDA kernel for transforming position coordinates from rvec to nbnxm layout.
 *
 * TODO:
 *  - rename kernel so naming matches with the other NBNXM kernels;
 *  - enable separate compilation unit

 * \param[in]     numColumns          Extent of cell-level parallelism.
 * \param[out]    gm_xq               Coordinates buffer in nbnxm layout.
 * \param[in]     gm_x                Coordinates buffer.
 * \param[in]     gm_atomIndex        Atom index mapping.
 * \param[in]     gm_numAtoms         Array of number of atoms.
 * \param[in]     gm_cellIndex        Array of cell indices.
 * \param[in]     cellOffset          First cell.
 * \param[in]     numAtomsPerCell     Number of atoms per cell.
 */
static __global__ void nbnxn_gpu_x_to_nbat_x_kernel(int numColumns,
                                                    float4* __restrict__ gm_xq,
                                                    const float3* __restrict__ gm_x,
                                                    const int* __restrict__ gm_atomIndex,
                                                    const int* __restrict__ gm_numAtoms,
                                                    const int* __restrict__ gm_cellIndex,
                                                    int cellOffset,
                                                    int numAtomsPerCell)
{

    const float farAway = -1000000.0F;

    // Map cell-level parallelism to y component of CUDA block index.
    int cxy = blockIdx.y;

    if (cxy < numColumns)
    {

        const int numAtoms = gm_numAtoms[cxy];
        const int offset   = (cellOffset + gm_cellIndex[cxy]) * numAtomsPerCell;

        const int threadIndex = blockIdx.x * blockDim.x + threadIdx.x;

        // Destination address where x should be stored in nbnxm layout. We use this cast here to
        // save only x, y and z components, not touching the w (q) component, which is pre-defined.
        float3* gm_xqDest = reinterpret_cast<float3*>(&gm_xq[threadIndex + offset]);

        // Perform layout conversion of each element.
        if (threadIndex < numAtoms)
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


//! Number of CUDA threads in a block
// TODO Optimize this through experimentation
constexpr static int c_bufOpsThreadsPerBlock = 128;

void launchNbnxmKernelTransformXToXq(const Grid&          grid,
                                     NbnxmGpu*            nb,
                                     DeviceBuffer<Float3> d_x,
                                     const DeviceStream&  deviceStream,
                                     const unsigned int   numColumnsMax,
                                     const int            gridId)
{
    const int numColumns      = grid.numColumns();
    const int cellOffset      = grid.cellOffset();
    const int numAtomsPerCell = grid.numAtomsPerCell();

    KernelLaunchConfig config;
    config.blockSize[0] = c_bufOpsThreadsPerBlock;
    config.blockSize[1] = 1;
    config.blockSize[2] = 1;
    config.gridSize[0] =
            gmx::divideRoundUp(grid.numCellsColumnMax() * numAtomsPerCell, c_bufOpsThreadsPerBlock);
    config.gridSize[1] = numColumns;
    config.gridSize[2] = 1;
    GMX_ASSERT(config.gridSize[0] > 0, "Can not have empty grid, early return above avoids this");
    config.sharedMemorySize = 0;

    auto       kernelFn      = nbnxn_gpu_x_to_nbat_x_kernel;
    float3*    d_xFloat3     = asFloat3(d_x);
    float4*    d_xq          = nb->atdat->xq;
    const int* d_atomIndices = nb->atomIndices;
    const int* d_cxy_na      = &nb->cxy_na[numColumnsMax * gridId];
    const int* d_cxy_ind     = &nb->cxy_ind[numColumnsMax * gridId];
    const auto kernelArgs    = prepareGpuKernelArguments(
            kernelFn, config, &numColumns, &d_xq, &d_xFloat3, &d_atomIndices, &d_cxy_na, &d_cxy_ind, &cellOffset, &numAtomsPerCell);
    launchGpuKernel(kernelFn, config, deviceStream, nullptr, "XbufferOps", kernelArgs);
}

} // namespace gmx
