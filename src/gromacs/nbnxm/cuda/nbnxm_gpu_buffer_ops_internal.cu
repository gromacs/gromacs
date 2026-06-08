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
#include "gromacs/gpu_utils/vectype_ops_cuda.h"
#include "gromacs/nbnxm/cuda/nbnxm_cuda_types.h"
#include "gromacs/nbnxm/grid.h"

namespace gmx
{

/*! \brief CUDA kernel for transforming position coordinates from rvec to nbnxm layout.
 *
 * Processes columns from multiple grids in a single kernel launch.
 * blockIdx.x covers atoms within a column, blockIdx.y maps to a global
 * column index across all grids. The grid membership of each column is
 * determined using \p gridParams: columnsPrefix[] holds the prefix sum
 * of per-grid column counts and binOffset[] the starting bin for each
 * grid.
 *
 * \param[in]     gridBegin           Index of first grid in gridset.
 * \param[in]     numColumnsMax       Max columns per grid (stride for per-grid arrays).
 * \param[in]     numAtomsPerBin      Number of atoms per bin (same for all grids).
 * \param[in]     gridParams          Per-grid parameters (prefix sums and bin offsets).
 * \param[out]    gm_xq               Coordinates buffer in nbnxm layout.
 * \param[in]     gm_x                Coordinates buffer.
 * \param[in]     gm_atomIndex        Atom index mapping.
 * \param[in]     gm_numAtoms         Array of number of atoms per column (all grids, stride numColumnsMax).
 * \param[in]     gm_binIndex         Array of bin indices per column (all grids, stride numColumnsMax).
 */
static __global__ void nbnxn_gpu_x_to_nbat_x_kernel(int                  gridBegin,
                                                    int                  numColumnsMax,
                                                    int                  numAtomsPerBin,
                                                    FusedXToXqGridParams gridParams,
                                                    float4* __restrict__ gm_xq,
                                                    const float3* __restrict__ gm_x,
                                                    const int* __restrict__ gm_atomIndex,
                                                    const int* __restrict__ gm_numAtoms,
                                                    const int* __restrict__ gm_binIndex)
{
    const int globalCol = blockIdx.y;

    // Start assuming the column belongs to grid 0, using constant [0] indices.
    int localCol  = globalCol;
    int gridIdx   = gridBegin;
    int binOffset = gridParams.binOffset[0];

    // Last-write-wins scan over strictly increasing columnsPrefix[]: the highest matching g is
    // the correct grid. Full unrolling with constant indices is required to avoid local memory
    // spills; unused entries use INT_MAX sentinels so they never match.
#pragma unroll
    for (int g = 1; g < c_maxGridsPerKernelLaunch; g++)
    {
        if (globalCol >= gridParams.columnsPrefix[g])
        {
            localCol  = globalCol - gridParams.columnsPrefix[g];
            gridIdx   = g + gridBegin;
            binOffset = gridParams.binOffset[g];
        }
    }

    // Access per-grid device arrays using numColumnsMax stride
    const int colDataOffset = numColumnsMax * gridIdx + localCol;
    const int numAtoms      = gm_numAtoms[colDataOffset];
    const int offset        = (binOffset + gm_binIndex[colDataOffset]) * numAtomsPerBin;

    const int threadIndex = blockIdx.x * blockDim.x + threadIdx.x;

    // Destination address where x should be stored in nbnxm layout. We use this cast here to
    // save only x, y and z components, not touching the w (q) component, which is pre-defined.
    float3* gm_xqDest = reinterpret_cast<float3*>(&gm_xq[threadIndex + offset]);

    // Perform layout conversion of each element.
    if (threadIndex < numAtoms)
    {
        *gm_xqDest = gm_x[gm_atomIndex[threadIndex + offset]];
    }
}


//! Number of CUDA threads in a block
// TODO Optimize this through experimentation
constexpr static int c_bufOpsThreadsPerBlock = 128;

void launchNbnxmKernelTransformXToXq(const FusedXToXqLaunchParams& launchParams,
                                     NbnxmGpu*                     nb,
                                     DeviceBuffer<Float3>          d_x,
                                     const DeviceStream&           deviceStream)
{
    KernelLaunchConfig config;
    config.blockSize[0] = c_bufOpsThreadsPerBlock;
    config.blockSize[1] = 1;
    config.blockSize[2] = 1;
    config.gridSize[0] = gmx::divideRoundUp(launchParams.maxNumAtomsPerColumn, c_bufOpsThreadsPerBlock);
    config.gridSize[1] = launchParams.totalNumColumns;
    config.gridSize[2] = 1;
    GMX_ASSERT(config.gridSize[0] > 0, "Can not have empty grid, caller should skip empty grids");
    config.sharedMemorySize = 0;

    auto        kernelFn            = nbnxn_gpu_x_to_nbat_x_kernel;
    float3*     d_xFloat3           = asFloat3(d_x);
    float4*     d_xq                = nb->atdat->xq;
    const int*  d_atomIndices       = nb->atomIndices;
    const int*  d_numAtomsPerColumn = nb->numAtomsPerCell;
    const int*  d_columnToBin       = nb->cellToBin;
    const auto& gridParams          = launchParams.gridParams;
    const int   numAtomsPerBin      = launchParams.numAtomsPerBin;
    const int   gridBegin           = launchParams.gridBegin;
    const int   numColumnsMaxInt    = launchParams.numColumnsMax;
    const auto  kernelArgs          = prepareGpuKernelArguments(kernelFn,
                                                      config,
                                                      &gridBegin,
                                                      &numColumnsMaxInt,
                                                      &numAtomsPerBin,
                                                      &gridParams,
                                                      &d_xq,
                                                      &d_xFloat3,
                                                      &d_atomIndices,
                                                      &d_numAtomsPerColumn,
                                                      &d_columnToBin);
    launchGpuKernel(kernelFn, config, deviceStream, nullptr, "XbufferOps", kernelArgs);
}

} // namespace gmx
