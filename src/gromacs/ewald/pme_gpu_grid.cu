/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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

/*! \internal \file
 * \brief Implements PME GPU halo exchange and PME GPU - Host FFT grid conversion
 * functions. These functions are used for PME decomposition in mixed-mode
 *
 * \author Gaurav Garg <gaugarg@nvidia.com>
 *
 * \ingroup module_ewald
 */

#include "gmxpre.h"

#include "pme_gpu_grid.h"

#include "config.h"

#include <cstdlib>

#include "gromacs/fft/parallel_3dfft.h"
#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/devicebuffer.cuh"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/timing/wallcycle.h"

#include "pme.cuh"
#include "pme_gpu_types.h"
#include "pme_gpu_types_host.h"
#include "pme_gpu_types_host_impl.h"

/*! \brief
 * A CUDA kernel which packs non-contiguous overlap data in all 8 neighboring directions
 *
 * \param[in] gm_realGrid          PME device grid
 * \param[out] gm_transferGrid*    device arrays used to pack data in 8-neighboring directions
 * \param[in] overlapSize*         halo size in 4 directions
 * \param[in] myGrid*              local domain size in X and Y dimension
 * \param[in] pmeSize              Local PME grid size
 *
 */
static __global__ void pmeGpuPackHaloExternal(const float* __restrict__ gm_realGrid,
                                              float* __restrict__ gm_transferGridUp,
                                              float* __restrict__ gm_transferGridDown,
                                              float* __restrict__ gm_transferGridLeft,
                                              float* __restrict__ gm_transferGridRight,
                                              float* __restrict__ gm_transferGridUpLeft,
                                              float* __restrict__ gm_transferGridDownLeft,
                                              float* __restrict__ gm_transferGridUpRight,
                                              float* __restrict__ gm_transferGridDownRight,
                                              int  overlapSizeUp,
                                              int  overlapSizeDown,
                                              int  overlapSizeLeft,
                                              int  overlapSizeRight,
                                              int  myGridX,
                                              int  myGridY,
                                              int3 pmeSize)
{
    int iz = threadIdx.x + blockIdx.x * blockDim.x;
    int iy = threadIdx.y + blockIdx.y * blockDim.y;
    int ix = threadIdx.z + blockIdx.z * blockDim.z;

    // we might get iz greater than pmeSize.z when pmeSize.z is not multiple of
    // threadsAlongZDim(see below), same for iy when it's not multiple of threadsAlongYDim
    if (iz >= pmeSize.z || iy >= myGridY)
    {
        return;
    }

    // up
    if (ix < overlapSizeUp)
    {
        int pmeIndex = (ix + pmeSize.x - overlapSizeUp) * pmeSize.y * pmeSize.z + iy * pmeSize.z + iz;
        int packedIndex                = ix * myGridY * pmeSize.z + iy * pmeSize.z + iz;
        gm_transferGridUp[packedIndex] = gm_realGrid[pmeIndex];
    }

    // down
    if (ix >= myGridX - overlapSizeDown)
    {
        int pmeIndex = (ix + overlapSizeDown) * pmeSize.y * pmeSize.z + iy * pmeSize.z + iz;
        int packedIndex = (ix - (myGridX - overlapSizeDown)) * myGridY * pmeSize.z + iy * pmeSize.z + iz;
        gm_transferGridDown[packedIndex] = gm_realGrid[pmeIndex];
    }

    // left
    if (iy < overlapSizeLeft)
    {
        int pmeIndex = ix * pmeSize.y * pmeSize.z + (iy + pmeSize.y - overlapSizeLeft) * pmeSize.z + iz;
        int packedIndex                  = ix * overlapSizeLeft * pmeSize.z + iy * pmeSize.z + iz;
        gm_transferGridLeft[packedIndex] = gm_realGrid[pmeIndex];
    }

    // right
    if (iy >= myGridY - overlapSizeRight)
    {
        int pmeIndex    = ix * pmeSize.y * pmeSize.z + (iy + overlapSizeRight) * pmeSize.z + iz;
        int packedIndex = ix * overlapSizeRight * pmeSize.z
                          + (iy - (myGridY - overlapSizeRight)) * pmeSize.z + iz;
        gm_transferGridRight[packedIndex] = gm_realGrid[pmeIndex];
    }

    // up left
    if (ix < overlapSizeUp && iy < overlapSizeLeft)
    {
        int pmeIndex = (ix + pmeSize.x - overlapSizeUp) * pmeSize.y * pmeSize.z
                       + (iy + pmeSize.y - overlapSizeLeft) * pmeSize.z + iz;
        int packedIndex                    = ix * overlapSizeLeft * pmeSize.z + iy * pmeSize.z + iz;
        gm_transferGridUpLeft[packedIndex] = gm_realGrid[pmeIndex];
    }

    // down left
    if (ix >= myGridX - overlapSizeDown && iy < overlapSizeLeft)
    {
        int pmeIndex = (ix + overlapSizeDown) * pmeSize.y * pmeSize.z
                       + (iy + pmeSize.y - overlapSizeLeft) * pmeSize.z + iz;
        int packedIndex =
                (ix - (myGridX - overlapSizeDown)) * overlapSizeLeft * pmeSize.z + iy * pmeSize.z + iz;
        gm_transferGridDownLeft[packedIndex] = gm_realGrid[pmeIndex];
    }

    // up right
    if (ix < overlapSizeUp && iy >= myGridY - overlapSizeRight)
    {
        int pmeIndex = (ix + pmeSize.x - overlapSizeUp) * pmeSize.y * pmeSize.z
                       + (iy + overlapSizeRight) * pmeSize.z + iz;
        int packedIndex = ix * overlapSizeRight * pmeSize.z
                          + (iy - (myGridY - overlapSizeRight)) * pmeSize.z + iz;
        gm_transferGridUpRight[packedIndex] = gm_realGrid[pmeIndex];
    }

    // down right
    if (ix >= myGridX - overlapSizeDown && iy >= myGridY - overlapSizeRight)
    {
        int pmeIndex = (ix + overlapSizeDown) * pmeSize.y * pmeSize.z
                       + (iy + overlapSizeRight) * pmeSize.z + iz;
        int packedIndex = (ix - (myGridX - overlapSizeDown)) * overlapSizeRight * pmeSize.z
                          + (iy - (myGridY - overlapSizeRight)) * pmeSize.z + iz;
        gm_transferGridDownRight[packedIndex] = gm_realGrid[pmeIndex];
    }
}

/*! \brief
 * A CUDA kernel which assigns data in halo region in all 8 neighboring directions
 *
 * \param[in] gm_realGrid          PME device grid
 * \param[out] gm_transferGrid*    packed data in 8-neighboring directions
 * \param[in] overlapSize*         halo size in 4 directions
 * \param[in] myGrid*              local domain size in X and Y dimension
 * \param[in] pmeSize              Local PME grid size
 */
static __global__ void pmeGpuUnpackHaloExternal(float* __restrict__ gm_realGrid,
                                                const float* __restrict__ gm_transferGridUp,
                                                const float* __restrict__ gm_transferGridDown,
                                                const float* __restrict__ gm_transferGridLeft,
                                                const float* __restrict__ gm_transferGridRight,
                                                const float* __restrict__ gm_transferGridUpLeft,
                                                const float* __restrict__ gm_transferGridDownLeft,
                                                const float* __restrict__ gm_transferGridUpRight,
                                                const float* __restrict__ gm_transferGridDownRight,
                                                int  overlapSizeUp,
                                                int  overlapSizeDown,
                                                int  overlapSizeLeft,
                                                int  overlapSizeRight,
                                                int  myGridX,
                                                int  myGridY,
                                                int3 pmeSize)
{
    int iz = threadIdx.x + blockIdx.x * blockDim.x;
    int iy = threadIdx.y + blockIdx.y * blockDim.y;
    int ix = threadIdx.z + blockIdx.z * blockDim.z;

    // we might get iz greater than pmeSize.z when pmeSize.z is not multiple of
    // threadsAlongZDim(see below), same for iy when it's not multiple of threadsAlongYDim
    if (iz >= pmeSize.z || iy >= myGridY)
    {
        return;
    }

    // up
    if (ix < overlapSizeUp)
    {
        int pmeIndex = (ix + pmeSize.x - overlapSizeUp) * pmeSize.y * pmeSize.z + iy * pmeSize.z + iz;
        int packedIndex       = ix * myGridY * pmeSize.z + iy * pmeSize.z + iz;
        gm_realGrid[pmeIndex] = gm_transferGridUp[packedIndex];
    }

    // down
    if (ix >= myGridX - overlapSizeDown)
    {
        int pmeIndex = (ix + overlapSizeDown) * pmeSize.y * pmeSize.z + iy * pmeSize.z + iz;
        int packedIndex = (ix - (myGridX - overlapSizeDown)) * myGridY * pmeSize.z + iy * pmeSize.z + iz;
        gm_realGrid[pmeIndex] = gm_transferGridDown[packedIndex];
    }

    // left
    if (iy < overlapSizeLeft)
    {
        int pmeIndex = ix * pmeSize.y * pmeSize.z + (iy + pmeSize.y - overlapSizeLeft) * pmeSize.z + iz;
        int packedIndex       = ix * overlapSizeLeft * pmeSize.z + iy * pmeSize.z + iz;
        gm_realGrid[pmeIndex] = gm_transferGridLeft[packedIndex];
    }

    // right
    if (iy >= myGridY - overlapSizeRight)
    {
        int pmeIndex    = ix * pmeSize.y * pmeSize.z + (iy + overlapSizeRight) * pmeSize.z + iz;
        int packedIndex = ix * overlapSizeRight * pmeSize.z
                          + (iy - (myGridY - overlapSizeRight)) * pmeSize.z + iz;
        gm_realGrid[pmeIndex] = gm_transferGridRight[packedIndex];
    }

    // up left
    if (ix < overlapSizeUp && iy < overlapSizeLeft)
    {
        int pmeIndex = (ix + pmeSize.x - overlapSizeUp) * pmeSize.y * pmeSize.z
                       + (iy + pmeSize.y - overlapSizeLeft) * pmeSize.z + iz;
        int packedIndex       = ix * overlapSizeLeft * pmeSize.z + iy * pmeSize.z + iz;
        gm_realGrid[pmeIndex] = gm_transferGridUpLeft[packedIndex];
    }

    // down left
    if (ix >= myGridX - overlapSizeDown && iy < overlapSizeLeft)
    {
        int pmeIndex = (ix + overlapSizeDown) * pmeSize.y * pmeSize.z
                       + (iy + pmeSize.y - overlapSizeLeft) * pmeSize.z + iz;
        int packedIndex =
                (ix - (myGridX - overlapSizeDown)) * overlapSizeLeft * pmeSize.z + iy * pmeSize.z + iz;
        gm_realGrid[pmeIndex] = gm_transferGridDownLeft[packedIndex];
    }

    // up right
    if (ix < overlapSizeUp && iy >= myGridY - overlapSizeRight)
    {
        int pmeIndex = (ix + pmeSize.x - overlapSizeUp) * pmeSize.y * pmeSize.z
                       + (iy + overlapSizeRight) * pmeSize.z + iz;
        int packedIndex = ix * overlapSizeRight * pmeSize.z
                          + (iy - (myGridY - overlapSizeRight)) * pmeSize.z + iz;
        gm_realGrid[pmeIndex] = gm_transferGridUpRight[packedIndex];
    }

    // down right
    if (ix >= myGridX - overlapSizeDown && iy >= myGridY - overlapSizeRight)
    {
        int pmeIndex = (ix + overlapSizeDown) * pmeSize.y * pmeSize.z
                       + (iy + overlapSizeRight) * pmeSize.z + iz;
        int packedIndex = (ix - (myGridX - overlapSizeDown)) * overlapSizeRight * pmeSize.z
                          + (iy - (myGridY - overlapSizeRight)) * pmeSize.z + iz;
        gm_realGrid[pmeIndex] = gm_transferGridDownRight[packedIndex];
    }
}

/*! \brief
 * A CUDA kernel which adds grid overlap data received from neighboring ranks
 *
 * \param[in] gm_realGrid          PME device grid
 * \param[out] gm_transferGrid*    packed data in 8-neighboring directions
 * \param[in] overlapSize*         halo size in 4 directions
 * \param[in] myGrid*              local domain size in X and Y dimension
 * \param[in] pmeSize              Local PME grid size
 */

static __global__ void pmeGpuUnpackAndAddHaloInternal(float* __restrict__ gm_realGrid,
                                                      const float* __restrict__ gm_transferGridUp,
                                                      const float* __restrict__ gm_transferGridDown,
                                                      const float* __restrict__ gm_transferGridLeft,
                                                      const float* __restrict__ gm_transferGridRight,
                                                      const float* __restrict__ gm_transferGridUpLeft,
                                                      const float* __restrict__ gm_transferGridDownLeft,
                                                      const float* __restrict__ gm_transferGridUpRight,
                                                      const float* __restrict__ gm_transferGridDownRight,
                                                      int  overlapSizeX,
                                                      int  overlapSizeY,
                                                      int  overlapUp,
                                                      int  overlapLeft,
                                                      int  myGridX,
                                                      int  myGridY,
                                                      int3 pmeSize)
{
    int iz = threadIdx.x + blockIdx.x * blockDim.x;
    int iy = threadIdx.y + blockIdx.y * blockDim.y;
    int ix = threadIdx.z + blockIdx.z * blockDim.z;

    // we might get iz greater than pmeSize.z when pmeSize.z is not multiple of
    // threadsAlongZDim(see below), same for iy when it's not multiple of threadsAlongYDim
    if (iz >= pmeSize.z || iy >= myGridY)
    {
        return;
    }

    int pmeIndex = ix * pmeSize.y * pmeSize.z + iy * pmeSize.z + iz;

    float val = gm_realGrid[pmeIndex];

    // up rank
    if (ix < overlapSizeX)
    {
        int packedIndex = ix * myGridY * pmeSize.z + iy * pmeSize.z + iz;
        val += gm_transferGridUp[packedIndex];
    }

    // down rank
    if (ix >= myGridX - overlapSizeX && overlapUp > 0)
    {
        int packedIndex = (ix - (myGridX - overlapSizeX)) * myGridY * pmeSize.z + iy * pmeSize.z + iz;
        val += gm_transferGridDown[packedIndex];
    }

    // left rank
    if (iy < overlapSizeY)
    {
        int packedIndex = ix * overlapSizeY * pmeSize.z + iy * pmeSize.z + iz;
        val += gm_transferGridLeft[packedIndex];
    }

    // right rank
    if (iy >= myGridY - overlapSizeY && overlapLeft > 0)
    {
        int packedIndex = ix * overlapSizeY * pmeSize.z + (iy - (myGridY - overlapSizeY)) * pmeSize.z + iz;
        val += gm_transferGridRight[packedIndex];
    }

    // up left rank
    if (ix < overlapSizeX && iy < overlapSizeY)
    {
        int packedIndex = ix * overlapSizeY * pmeSize.z + iy * pmeSize.z + iz;
        val += gm_transferGridUpLeft[packedIndex];
    }

    // up right rank
    if (ix < overlapSizeX && iy >= myGridY - overlapSizeY && overlapLeft > 0)
    {
        int packedIndex = ix * overlapSizeY * pmeSize.z + (iy - (myGridY - overlapSizeY)) * pmeSize.z + iz;
        val += gm_transferGridUpRight[packedIndex];
    }

    // down left rank
    if (ix >= myGridX - overlapSizeX && overlapUp > 0 && iy < overlapSizeY)
    {
        int packedIndex = (ix - (myGridX - overlapSizeX)) * overlapSizeY * pmeSize.z + iy * pmeSize.z + iz;
        val += gm_transferGridDownLeft[packedIndex];
    }

    // down right rank
    if (ix >= myGridX - overlapSizeX && overlapUp > 0 && iy >= myGridY - overlapSizeY && overlapLeft > 0)
    {
        int packedIndex = (ix - (myGridX - overlapSizeX)) * overlapSizeY * pmeSize.z
                          + (iy - (myGridY - overlapSizeY)) * pmeSize.z + iz;
        val += gm_transferGridDownRight[packedIndex];
    }

    gm_realGrid[pmeIndex] = val;
}

/*! \brief
 * A CUDA kernel which packs non-contiguous overlap data in all 8 neighboring directions
 *
 * \param[in] gm_realGrid          PME device grid
 * \param[out] gm_transferGrid*    packed data in 8-neighboring directions
 * \param[in] overlapSize*         halo size in 4 directions
 * \param[in] myGrid*              local domain size in X and Y dimension
 * \param[in] pmeSize              Local PME grid size
 */
static __global__ void pmeGpuPackHaloInternal(const float* __restrict__ gm_realGrid,
                                              float* __restrict__ gm_transferGridUp,
                                              float* __restrict__ gm_transferGridDown,
                                              float* __restrict__ gm_transferGridLeft,
                                              float* __restrict__ gm_transferGridRight,
                                              float* __restrict__ gm_transferGridUpLeft,
                                              float* __restrict__ gm_transferGridDownLeft,
                                              float* __restrict__ gm_transferGridUpRight,
                                              float* __restrict__ gm_transferGridDownRight,
                                              int  overlapSizeX,
                                              int  overlapSizeY,
                                              int  overlapUp,
                                              int  overlapLeft,
                                              int  myGridX,
                                              int  myGridY,
                                              int3 pmeSize)
{
    int iz = threadIdx.x + blockIdx.x * blockDim.x;
    int iy = threadIdx.y + blockIdx.y * blockDim.y;
    int ix = threadIdx.z + blockIdx.z * blockDim.z;

    // we might get iz greater than pmeSize.z when pmeSize.z is not multiple of
    // threadsAlongZDim(see below), same for iy when it's not multiple of threadsAlongYDim
    if (iz >= pmeSize.z || iy >= myGridY)
    {
        return;
    }

    int pmeIndex = ix * pmeSize.y * pmeSize.z + iy * pmeSize.z + iz;

    float val = gm_realGrid[pmeIndex];

    // up rank
    if (ix < overlapSizeX)
    {
        int packedIndex                = ix * myGridY * pmeSize.z + iy * pmeSize.z + iz;
        gm_transferGridUp[packedIndex] = val;
    }

    // down rank
    if (ix >= myGridX - overlapSizeX && overlapUp > 0)
    {
        int packedIndex = (ix - (myGridX - overlapSizeX)) * myGridY * pmeSize.z + iy * pmeSize.z + iz;
        gm_transferGridDown[packedIndex] = val;
    }

    // left rank
    if (iy < overlapSizeY)
    {
        int packedIndex                  = ix * overlapSizeY * pmeSize.z + iy * pmeSize.z + iz;
        gm_transferGridLeft[packedIndex] = val;
    }

    // right rank
    if (iy >= myGridY - overlapSizeY && overlapLeft > 0)
    {
        int packedIndex = ix * overlapSizeY * pmeSize.z + (iy - (myGridY - overlapSizeY)) * pmeSize.z + iz;
        gm_transferGridRight[packedIndex] = val;
    }

    // up left rank
    if (ix < overlapSizeX && iy < overlapSizeY)
    {
        int packedIndex                    = ix * overlapSizeY * pmeSize.z + iy * pmeSize.z + iz;
        gm_transferGridUpLeft[packedIndex] = val;
    }

    // down left rank
    if (ix >= myGridX - overlapSizeX && overlapUp > 0 && iy < overlapSizeY)
    {
        int packedIndex = (ix - (myGridX - overlapSizeX)) * overlapSizeY * pmeSize.z + iy * pmeSize.z + iz;
        gm_transferGridDownLeft[packedIndex] = val;
    }

    // up right rank
    if (ix < overlapSizeX && iy >= myGridY - overlapSizeY && overlapLeft > 0)
    {
        int packedIndex = ix * overlapSizeY * pmeSize.z + (iy - (myGridY - overlapSizeY)) * pmeSize.z + iz;
        gm_transferGridUpRight[packedIndex] = val;
    }

    // down right rank
    if (ix >= myGridX - overlapSizeX && overlapUp > 0 && iy >= myGridY - overlapSizeY && overlapLeft > 0)
    {
        int packedIndex = (ix - (myGridX - overlapSizeX)) * overlapSizeY * pmeSize.z
                          + (iy - (myGridY - overlapSizeY)) * pmeSize.z + iz;
        gm_transferGridDownRight[packedIndex] = val;
    }
}

/*! \brief
 * A CUDA kernel which copies data from pme grid to FFT grid and back
 *
 * \param[in] gm_pmeGrid          local PME grid
 * \param[in] gm_fftGrid          local FFT grid
 * \param[in] fftNData           local FFT grid size without padding
 * \param[in] fftSize            local FFT grid padded size
 * \param[in] pmeSize            local PME grid padded size
 *
 * \tparam  pmeToFft               A boolean which tells if this is conversion from PME grid to FFT grid or reverse
 */
template<bool pmeToFft>
static __global__ void pmegrid_to_fftgrid(float* __restrict__ gm_realGrid,
                                          float* __restrict__ gm_fftGrid,
                                          int3 fftNData,
                                          int3 fftSize,
                                          int3 pmeSize)
{
    int iz = threadIdx.x + blockIdx.x * blockDim.x;
    int iy = threadIdx.y + blockIdx.y * blockDim.y;
    int ix = threadIdx.z + blockIdx.z * blockDim.z;

    if (ix >= fftNData.x || iy >= fftNData.y || iz >= fftNData.z)
    {
        return;
    }

    int fftidx   = ix * fftSize.y * fftSize.z + iy * fftSize.z + iz;
    int pmeIndex = ix * pmeSize.y * pmeSize.z + iy * pmeSize.z + iz;

    if (pmeToFft)
    {
        gm_fftGrid[fftidx] = gm_realGrid[pmeIndex];
    }
    else
    {
        gm_realGrid[pmeIndex] = gm_fftGrid[fftidx];
    }
}

/*! \brief
 * Launches CUDA kernel to pack non-contiguous external halo data
 */
static void packHaloDataExternal(const PmeGpu*       pmeGpu,
                                 int                 overlapUp,
                                 int                 overlapDown,
                                 int                 overlapLeft,
                                 int                 overlapRight,
                                 int                 myGridX,
                                 int                 myGridY,
                                 const ivec&         pmeSize,
                                 DeviceBuffer<float> realGrid,
                                 DeviceBuffer<float> packedGridUp,
                                 DeviceBuffer<float> packedGridDown,
                                 DeviceBuffer<float> packedGridLeft,
                                 DeviceBuffer<float> packedGridRight,
                                 DeviceBuffer<float> packedGridUpLeft,
                                 DeviceBuffer<float> packedGridDownLeft,
                                 DeviceBuffer<float> packedGridUpRight,
                                 DeviceBuffer<float> packedGridDownRight)
{
    // Keeping threadsAlongZDim same as warp size for better coalescing.
    // Not keeping to higher value such as 64 to avoid high masked out
    // inactive threads as FFT grid sizes tend to be quite small
    const int threadsAlongZDim = 32;
    const int threadsAlongYDim = 4;

    // right grid
    KernelLaunchConfig config;
    config.blockSize[0]     = threadsAlongZDim;
    config.blockSize[1]     = threadsAlongYDim;
    config.blockSize[2]     = 1;
    config.gridSize[0]      = gmx::divideRoundUp(pmeSize[ZZ], threadsAlongZDim);
    config.gridSize[1]      = gmx::divideRoundUp(myGridY, threadsAlongYDim);
    config.gridSize[2]      = myGridX;
    config.sharedMemorySize = 0;

    auto kernelFn   = pmeGpuPackHaloExternal;
    auto kernelArgs = prepareGpuKernelArguments(kernelFn,
                                                config,
                                                &realGrid,
                                                &packedGridUp,
                                                &packedGridDown,
                                                &packedGridLeft,
                                                &packedGridRight,
                                                &packedGridUpLeft,
                                                &packedGridDownLeft,
                                                &packedGridUpRight,
                                                &packedGridDownRight,
                                                &overlapUp,
                                                &overlapDown,
                                                &overlapLeft,
                                                &overlapRight,
                                                &myGridX,
                                                &myGridY,
                                                &pmeSize);

    launchGpuKernel(kernelFn,
                    config,
                    pmeGpu->archSpecific->pmeStream_,
                    nullptr,
                    "PME Domdec GPU Pack Grid Halo Exchange",
                    kernelArgs);
}

/*! \brief
 * Launches CUDA kernel to pack non-contiguous internal halo data
 */
static void packHaloDataInternal(const PmeGpu*       pmeGpu,
                                 int                 overlapSizeX,
                                 int                 overlapSizeY,
                                 int                 overlapUp,
                                 int                 overlapLeft,
                                 int                 myGridX,
                                 int                 myGridY,
                                 const ivec&         pmeSize,
                                 DeviceBuffer<float> realGrid,
                                 DeviceBuffer<float> packedGridUp,
                                 DeviceBuffer<float> packedGridDown,
                                 DeviceBuffer<float> packedGridLeft,
                                 DeviceBuffer<float> packedGridRight,
                                 DeviceBuffer<float> packedGridUpLeft,
                                 DeviceBuffer<float> packedGridDownLeft,
                                 DeviceBuffer<float> packedGridUpRight,
                                 DeviceBuffer<float> packedGridDownRight)
{
    // Keeping threadsAlongZDim same as warp size for better coalescing,
    // Not keeping to higher value such as 64 to avoid high masked out
    // inactive threads as FFT grid sizes tend to be quite small
    const int threadsAlongZDim = 32;
    const int threadsAlongYDim = 4;

    // right grid
    KernelLaunchConfig config;
    config.blockSize[0]     = threadsAlongZDim;
    config.blockSize[1]     = threadsAlongYDim;
    config.blockSize[2]     = 1;
    config.gridSize[0]      = gmx::divideRoundUp(pmeSize[ZZ], threadsAlongZDim);
    config.gridSize[1]      = gmx::divideRoundUp(myGridY, threadsAlongYDim);
    config.gridSize[2]      = myGridX;
    config.sharedMemorySize = 0;

    auto kernelFn   = pmeGpuPackHaloInternal;
    auto kernelArgs = prepareGpuKernelArguments(kernelFn,
                                                config,
                                                &realGrid,
                                                &packedGridUp,
                                                &packedGridDown,
                                                &packedGridLeft,
                                                &packedGridRight,
                                                &packedGridUpLeft,
                                                &packedGridDownLeft,
                                                &packedGridUpRight,
                                                &packedGridDownRight,
                                                &overlapSizeX,
                                                &overlapSizeY,
                                                &overlapUp,
                                                &overlapLeft,
                                                &myGridX,
                                                &myGridY,
                                                &pmeSize);

    launchGpuKernel(kernelFn,
                    config,
                    pmeGpu->archSpecific->pmeStream_,
                    nullptr,
                    "PME Domdec GPU Pack Grid Halo Exchange",
                    kernelArgs);
}


/*! \brief
 * Launches CUDA kernel to unpack and reduce overlap data
 */
static void unpackAndAddHaloDataInternal(const PmeGpu*       pmeGpu,
                                         int                 overlapSizeX,
                                         int                 overlapSizeY,
                                         int                 overlapUp,
                                         int                 overlapLeft,
                                         int                 myGridX,
                                         int                 myGridY,
                                         const ivec&         pmeSize,
                                         DeviceBuffer<float> realGrid,
                                         DeviceBuffer<float> packedGridUp,
                                         DeviceBuffer<float> packedGridDown,
                                         DeviceBuffer<float> packedGridLeft,
                                         DeviceBuffer<float> packedGridRight,
                                         DeviceBuffer<float> packedGridUpLeft,
                                         DeviceBuffer<float> packedGridDownLeft,
                                         DeviceBuffer<float> packedGridUpRight,
                                         DeviceBuffer<float> packedGridDownRight)
{
    // Keeping threadsAlongZDim same as warp size for better coalescing,
    // Not keeping to higher value such as 64 to avoid high masked out
    // inactive threads as FFT grid sizes tend to be quite small
    const int threadsAlongZDim = 32;
    const int threadsAlongYDim = 4;

    // right grid
    KernelLaunchConfig config;
    config.blockSize[0]     = threadsAlongZDim;
    config.blockSize[1]     = threadsAlongYDim;
    config.blockSize[2]     = 1;
    config.gridSize[0]      = gmx::divideRoundUp(pmeSize[ZZ], threadsAlongZDim);
    config.gridSize[1]      = gmx::divideRoundUp(myGridY, threadsAlongYDim);
    config.gridSize[2]      = myGridX;
    config.sharedMemorySize = 0;

    auto kernelFn = pmeGpuUnpackAndAddHaloInternal;

    auto kernelArgs = prepareGpuKernelArguments(kernelFn,
                                                config,
                                                &realGrid,
                                                &packedGridUp,
                                                &packedGridDown,
                                                &packedGridLeft,
                                                &packedGridRight,
                                                &packedGridUpLeft,
                                                &packedGridDownLeft,
                                                &packedGridUpRight,
                                                &packedGridDownRight,
                                                &overlapSizeX,
                                                &overlapSizeY,
                                                &overlapUp,
                                                &overlapLeft,
                                                &myGridX,
                                                &myGridY,
                                                &pmeSize);


    launchGpuKernel(kernelFn,
                    config,
                    pmeGpu->archSpecific->pmeStream_,
                    nullptr,
                    "PME Domdec GPU Pack Grid Halo Exchange",
                    kernelArgs);
}

/*! \brief
 * Launches CUDA kernel to initialize overlap data
 */
static void unpackHaloDataExternal(const PmeGpu*       pmeGpu,
                                   int                 overlapUp,
                                   int                 overlapDown,
                                   int                 overlapLeft,
                                   int                 overlapRight,
                                   int                 myGridX,
                                   int                 myGridY,
                                   const ivec&         pmeSize,
                                   DeviceBuffer<float> realGrid,
                                   DeviceBuffer<float> packedGridUp,
                                   DeviceBuffer<float> packedGridDown,
                                   DeviceBuffer<float> packedGridLeft,
                                   DeviceBuffer<float> packedGridRight,
                                   DeviceBuffer<float> packedGridUpLeft,
                                   DeviceBuffer<float> packedGridDownLeft,
                                   DeviceBuffer<float> packedGridUpRight,
                                   DeviceBuffer<float> packedGridDownRight)
{
    // Keeping threadsAlongZDim same as warp size for better coalescing,
    // Not keeping to higher value such as 64 to avoid high masked out
    // inactive threads as FFT grid sizes tend to be quite small
    const int threadsAlongZDim = 32;
    const int threadsAlongYDim = 4;

    // right grid
    KernelLaunchConfig config;
    config.blockSize[0]     = threadsAlongZDim;
    config.blockSize[1]     = threadsAlongYDim;
    config.blockSize[2]     = 1;
    config.gridSize[0]      = gmx::divideRoundUp(pmeSize[ZZ], threadsAlongZDim);
    config.gridSize[1]      = gmx::divideRoundUp(myGridY, threadsAlongYDim);
    config.gridSize[2]      = myGridX;
    config.sharedMemorySize = 0;

    auto kernelFn   = pmeGpuUnpackHaloExternal;
    auto kernelArgs = prepareGpuKernelArguments(kernelFn,
                                                config,
                                                &realGrid,
                                                &packedGridUp,
                                                &packedGridDown,
                                                &packedGridLeft,
                                                &packedGridRight,
                                                &packedGridUpLeft,
                                                &packedGridDownLeft,
                                                &packedGridUpRight,
                                                &packedGridDownRight,
                                                &overlapUp,
                                                &overlapDown,
                                                &overlapLeft,
                                                &overlapRight,
                                                &myGridX,
                                                &myGridY,
                                                &pmeSize);


    launchGpuKernel(kernelFn,
                    config,
                    pmeGpu->archSpecific->pmeStream_,
                    nullptr,
                    "PME Domdec GPU Pack Grid Halo Exchange",
                    kernelArgs);
}

/*! \brief
 * utility function to send and recv halo data from neighboring ranks
 */
static void receiveAndSend(DeviceBuffer<float> sendBuf,
                           int                 sendCount,
                           int                 dest,
                           MPI_Request*        sendRequest,
                           DeviceBuffer<float> recvBuf,
                           int                 recvCount,
                           int                 src,
                           MPI_Request*        recvRequest,
                           int                 tag,
                           MPI_Comm            comm)
{
    // send data to dest rank and recv from src rank
    MPI_Irecv(recvBuf, recvCount, MPI_FLOAT, src, tag, comm, recvRequest);

    MPI_Isend(sendBuf, sendCount, MPI_FLOAT, dest, tag, comm, sendRequest);
}

void pmeGpuGridHaloExchange(const PmeGpu* pmeGpu, gmx_wallcycle* wcycle)
{
#if GMX_MPI
    // Note here we are assuming that width of the chunks is not so small that we need to
    // transfer to/from multiple ranks i.e. that the distributed grid contains chunks at least order-1 points wide.

    auto* kernelParamsPtr = pmeGpu->kernelParams.get();
    ivec  localPmeSize;
    localPmeSize[XX] = kernelParamsPtr->grid.realGridSizePadded[XX];
    localPmeSize[YY] = kernelParamsPtr->grid.realGridSizePadded[YY];
    localPmeSize[ZZ] = kernelParamsPtr->grid.realGridSizePadded[ZZ];

    int overlapX = pmeGpu->haloExchange->haloSizeX[gmx::DirectionX::Center];
    int overlapY = pmeGpu->haloExchange->haloSizeY[gmx::DirectionY::Center];

    int overlapDown = pmeGpu->haloExchange->haloSizeX[gmx::DirectionX::Down];
    int overlapUp   = pmeGpu->haloExchange->haloSizeX[gmx::DirectionX::Up];

    int overlapRight = pmeGpu->haloExchange->haloSizeY[gmx::DirectionY::Right];
    int overlapLeft  = pmeGpu->haloExchange->haloSizeY[gmx::DirectionY::Left];

    int myGridX = pmeGpu->haloExchange->gridSizeX;
    int myGridY = pmeGpu->haloExchange->gridSizeY;

    int sizeX = pmeGpu->common->nnodesX;
    int down  = pmeGpu->haloExchange->ranksX[gmx::DirectionX::Down];
    int up    = pmeGpu->haloExchange->ranksX[gmx::DirectionX::Up];

    int sizeY = pmeGpu->common->nnodesY;
    int right = pmeGpu->haloExchange->ranksY[gmx::DirectionY::Right];
    int left  = pmeGpu->haloExchange->ranksY[gmx::DirectionY::Left];

    for (int gridIndex = 0; gridIndex < pmeGpu->common->ngrids; gridIndex++)
    {
        MPI_Request req[16];
        int         reqCount = 0;
        float*      realGrid = pmeGpu->kernelParams->grid.d_realGrid[gridIndex];

        float* sendGridUp =
                pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Up][gmx::DirectionY::Center];
        float* sendGridDown =
                pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Down][gmx::DirectionY::Center];

        // no need to pack if slab-decomposition in X-dimension as data is already contiguous
        if (pmeGpu->common->nnodesY == 1)
        {
            int sendOffsetDown = myGridX * localPmeSize[YY] * localPmeSize[ZZ];
            int sendOffsetUp = (localPmeSize[XX] - overlapUp) * localPmeSize[YY] * localPmeSize[ZZ];
            sendGridUp       = &realGrid[sendOffsetUp];
            sendGridDown     = &realGrid[sendOffsetDown];
        }
        else
        {
            wallcycle_start(wcycle, WallCycleCounter::LaunchGpuPme);

            // launch packing kernel
            packHaloDataExternal(
                    pmeGpu,
                    overlapUp,
                    overlapDown,
                    overlapLeft,
                    overlapRight,
                    myGridX,
                    myGridY,
                    localPmeSize,
                    realGrid,
                    sendGridUp,
                    sendGridDown,
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Center][gmx::DirectionY::Left],
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Center][gmx::DirectionY::Right],
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Up][gmx::DirectionY::Left],
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Down][gmx::DirectionY::Left],
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Up][gmx::DirectionY::Right],
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Down][gmx::DirectionY::Right]);

            wallcycle_stop(wcycle, WallCycleCounter::LaunchGpuPme);
        }

        wallcycle_start(wcycle, WallCycleCounter::WaitGpuPmeSpread);

        // Make sure data is ready on GPU before MPI communication.
        // Wait for spread to finish in case of slab decomposition along X-dimension and
        // wait for packing to finish otherwise.
        // Todo: Consider using events to create dependcy on spread
        pmeGpu->archSpecific->pmeStream_.synchronize();

        wallcycle_stop(wcycle, WallCycleCounter::WaitGpuPmeSpread);


        wallcycle_start(wcycle, WallCycleCounter::PmeHaloExchangeComm);

        // major dimension
        if (sizeX > 1)
        {
            constexpr int mpiTag = 403; // Arbitrarily chosen

            // send data to down rank and recv from up rank
            receiveAndSend(sendGridDown,
                           overlapDown * myGridY * localPmeSize[ZZ],
                           down,
                           &req[reqCount],
                           pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Up][gmx::DirectionY::Center],
                           overlapX * myGridY * localPmeSize[ZZ],
                           up,
                           &req[reqCount + 1],
                           mpiTag,
                           pmeGpu->common->mpiCommX);
            reqCount += 2;

            if (overlapUp > 0)
            {
                // send data to up rank and recv from down rank
                receiveAndSend(
                        sendGridUp,
                        overlapUp * myGridY * localPmeSize[ZZ],
                        up,
                        &req[reqCount],
                        pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Down][gmx::DirectionY::Center],
                        overlapX * myGridY * localPmeSize[ZZ],
                        down,
                        &req[reqCount + 1],
                        mpiTag,
                        pmeGpu->common->mpiCommX);
                reqCount += 2;
            }
        }

        // minor dimension
        if (sizeY > 1)
        {
            constexpr int mpiTag = 404; // Arbitrarily chosen

            // recv from left rank and send to right rank
            receiveAndSend(
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Center][gmx::DirectionY::Right],
                    overlapRight * myGridX * localPmeSize[ZZ],
                    right,
                    &req[reqCount],
                    pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Center][gmx::DirectionY::Left],
                    overlapY * myGridX * localPmeSize[ZZ],
                    left,
                    &req[reqCount + 1],
                    mpiTag,
                    pmeGpu->common->mpiCommY);
            reqCount += 2;

            if (overlapLeft > 0)
            {
                // recv from right rank and send data to left rank
                receiveAndSend(
                        pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Center][gmx::DirectionY::Left],
                        overlapLeft * myGridX * localPmeSize[ZZ],
                        left,
                        &req[reqCount],
                        pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Center][gmx::DirectionY::Right],
                        overlapY * myGridX * localPmeSize[ZZ],
                        right,
                        &req[reqCount + 1],
                        mpiTag,
                        pmeGpu->common->mpiCommY);
                reqCount += 2;
            }
        }

        if (sizeX > 1 && sizeY > 1)
        {
            int rankUpLeft   = up * sizeY + left;
            int rankDownLeft = down * sizeY + left;

            int rankUpRight   = up * sizeY + right;
            int rankDownRight = down * sizeY + right;

            constexpr int mpiTag = 405; // Arbitrarily chosen

            // send data to down rank and recv from up rank
            receiveAndSend(
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Down][gmx::DirectionY::Right],
                    overlapDown * overlapRight * localPmeSize[ZZ],
                    rankDownRight,
                    &req[reqCount],
                    pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Up][gmx::DirectionY::Left],
                    overlapX * overlapY * localPmeSize[ZZ],
                    rankUpLeft,
                    &req[reqCount + 1],
                    mpiTag,
                    pmeGpu->common->mpiComm);
            reqCount += 2;

            if (overlapLeft > 0)
            {
                // send data to down left rank and recv from up right rank
                receiveAndSend(
                        pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Down][gmx::DirectionY::Left],
                        overlapDown * overlapLeft * localPmeSize[ZZ],
                        rankDownLeft,
                        &req[reqCount],
                        pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Up][gmx::DirectionY::Right],
                        overlapX * overlapY * localPmeSize[ZZ],
                        rankUpRight,
                        &req[reqCount + 1],
                        mpiTag,
                        pmeGpu->common->mpiComm);
                reqCount += 2;
            }

            if (overlapUp > 0)
            {
                // send data to up right rank and recv from down left rank
                receiveAndSend(
                        pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Up][gmx::DirectionY::Right],
                        overlapUp * overlapRight * localPmeSize[ZZ],
                        rankUpRight,
                        &req[reqCount],
                        pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Down][gmx::DirectionY::Left],
                        overlapX * overlapY * localPmeSize[ZZ],
                        rankDownLeft,
                        &req[reqCount + 1],
                        mpiTag,
                        pmeGpu->common->mpiComm);
                reqCount += 2;
            }

            if (overlapUp > 0 && overlapLeft > 0)
            {
                // send data to up left rank and recv from down right rank
                receiveAndSend(
                        pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Up][gmx::DirectionY::Left],
                        overlapUp * overlapLeft * localPmeSize[ZZ],
                        rankUpLeft,
                        &req[reqCount],
                        pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Down][gmx::DirectionY::Right],
                        overlapX * overlapY * localPmeSize[ZZ],
                        rankDownRight,
                        &req[reqCount + 1],
                        mpiTag,
                        pmeGpu->common->mpiComm);
                reqCount += 2;
            }
        }

        MPI_Waitall(reqCount, req, MPI_STATUSES_IGNORE);

        wallcycle_stop(wcycle, WallCycleCounter::PmeHaloExchangeComm);

        wallcycle_start(wcycle, WallCycleCounter::LaunchGpuPme);

        // reduce halo data
        unpackAndAddHaloDataInternal(
                pmeGpu,
                overlapX,
                overlapY,
                overlapUp,
                overlapLeft,
                myGridX,
                myGridY,
                localPmeSize,
                realGrid,
                pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Up][gmx::DirectionY::Center],
                pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Down][gmx::DirectionY::Center],
                pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Center][gmx::DirectionY::Left],
                pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Center][gmx::DirectionY::Right],
                pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Up][gmx::DirectionY::Left],
                pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Down][gmx::DirectionY::Left],
                pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Up][gmx::DirectionY::Right],
                pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Down][gmx::DirectionY::Right]);

        wallcycle_stop(wcycle, WallCycleCounter::LaunchGpuPme);
    }
#else
    GMX_UNUSED_VALUE(pmeGpu);
#endif
}

void pmeGpuGridHaloExchangeReverse(const PmeGpu* pmeGpu, gmx_wallcycle* wcycle)
{
#if GMX_MPI
    // Note here we are assuming that width of the chunks is not so small that we need to
    // transfer to/from multiple ranks i.e. that the distributed grid contains chunks at least order-1 points wide.

    auto* kernelParamsPtr = pmeGpu->kernelParams.get();
    ivec  localPmeSize;
    localPmeSize[XX] = kernelParamsPtr->grid.realGridSizePadded[XX];
    localPmeSize[YY] = kernelParamsPtr->grid.realGridSizePadded[YY];
    localPmeSize[ZZ] = kernelParamsPtr->grid.realGridSizePadded[ZZ];

    int overlapX = pmeGpu->haloExchange->haloSizeX[gmx::DirectionX::Center];
    int overlapY = pmeGpu->haloExchange->haloSizeY[gmx::DirectionY::Center];

    int overlapDown = pmeGpu->haloExchange->haloSizeX[gmx::DirectionX::Down];
    int overlapUp   = pmeGpu->haloExchange->haloSizeX[gmx::DirectionX::Up];

    int overlapRight = pmeGpu->haloExchange->haloSizeY[gmx::DirectionY::Right];
    int overlapLeft  = pmeGpu->haloExchange->haloSizeY[gmx::DirectionY::Left];

    int myGridX = pmeGpu->haloExchange->gridSizeX;
    int myGridY = pmeGpu->haloExchange->gridSizeY;

    int sizeX = pmeGpu->common->nnodesX;
    int down  = pmeGpu->haloExchange->ranksX[gmx::DirectionX::Down];
    int up    = pmeGpu->haloExchange->ranksX[gmx::DirectionX::Up];

    int sizeY = pmeGpu->common->nnodesY;
    int right = pmeGpu->haloExchange->ranksY[gmx::DirectionY::Right];
    int left  = pmeGpu->haloExchange->ranksY[gmx::DirectionY::Left];

    for (int gridIndex = 0; gridIndex < pmeGpu->common->ngrids; gridIndex++)
    {
        MPI_Request req[16];
        int         reqCount = 0;

        float* realGrid = pmeGpu->kernelParams->grid.d_realGrid[gridIndex];

        float* sendGridUp =
                pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Up][gmx::DirectionY::Center];
        float* sendGridDown =
                pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Down][gmx::DirectionY::Center];
        float* recvGridUp =
                pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Up][gmx::DirectionY::Center];
        float* recvGridDown =
                pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Down][gmx::DirectionY::Center];

        // no need to pack if slab-decomposition in X-dimension as data is already contiguous
        if (sizeY == 1)
        {
            int sendOffsetUp   = 0;
            int sendOffsetDown = (myGridX - overlapX) * localPmeSize[YY] * localPmeSize[ZZ];
            int recvOffsetUp = (localPmeSize[XX] - overlapUp) * localPmeSize[YY] * localPmeSize[ZZ];
            int recvOffsetDown = myGridX * localPmeSize[YY] * localPmeSize[ZZ];
            sendGridUp         = &realGrid[sendOffsetUp];
            sendGridDown       = &realGrid[sendOffsetDown];
            recvGridUp         = &realGrid[recvOffsetUp];
            recvGridDown       = &realGrid[recvOffsetDown];
        }
        else
        {
            wallcycle_start(wcycle, WallCycleCounter::LaunchGpuPme);

            // launch packing kernel
            packHaloDataInternal(
                    pmeGpu,
                    overlapX,
                    overlapY,
                    overlapUp,
                    overlapLeft,
                    myGridX,
                    myGridY,
                    localPmeSize,
                    realGrid,
                    sendGridUp,
                    sendGridDown,
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Center][gmx::DirectionY::Left],
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Center][gmx::DirectionY::Right],
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Up][gmx::DirectionY::Left],
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Down][gmx::DirectionY::Left],
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Up][gmx::DirectionY::Right],
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Down][gmx::DirectionY::Right]);

            wallcycle_stop(wcycle, WallCycleCounter::LaunchGpuPme);
        }

        wallcycle_start(wcycle, WallCycleCounter::WaitGpuFftToPmeGrid);

        // Make sure data is ready on GPU before MPI communication.
        // Wait for FFT to PME grid conversion to finish in case of slab decomposition along X-dimension and
        // wait for packing to finish otherwise.
        // Todo: Consider using events to create dependcy on FFT->PME grid operation
        pmeGpu->archSpecific->pmeStream_.synchronize();

        wallcycle_stop(wcycle, WallCycleCounter::WaitGpuFftToPmeGrid);

        wallcycle_start(wcycle, WallCycleCounter::PmeHaloExchangeComm);

        // major dimension
        if (sizeX > 1)
        {
            constexpr int mpiTag = 406; // Arbitrarily chosen

            // send data to up rank and recv from down rank
            receiveAndSend(sendGridUp,
                           overlapX * myGridY * localPmeSize[ZZ],
                           up,
                           &req[reqCount],
                           recvGridDown,
                           overlapDown * myGridY * localPmeSize[ZZ],
                           down,
                           &req[reqCount + 1],
                           mpiTag,
                           pmeGpu->common->mpiCommX);
            reqCount += 2;

            if (overlapUp > 0)
            {
                // send data to down rank and recv from up rank
                receiveAndSend(sendGridDown,
                               overlapX * myGridY * localPmeSize[ZZ],
                               down,
                               &req[reqCount],
                               recvGridUp,
                               overlapUp * myGridY * localPmeSize[ZZ],
                               up,
                               &req[reqCount + 1],
                               mpiTag,
                               pmeGpu->common->mpiCommX);
                reqCount += 2;
            }
        }

        // minor dimension
        if (sizeY > 1)
        {
            constexpr int mpiTag = 407; // Arbitrarily chosen

            // recv from right rank and send data to left rank
            receiveAndSend(
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Center][gmx::DirectionY::Left],
                    overlapY * myGridX * localPmeSize[ZZ],
                    left,
                    &req[reqCount],
                    pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Center][gmx::DirectionY::Right],
                    overlapRight * myGridX * localPmeSize[ZZ],
                    right,
                    &req[reqCount + 1],
                    mpiTag,
                    pmeGpu->common->mpiCommY);
            reqCount += 2;

            if (overlapLeft > 0)
            {
                // recv from left rank and send data to right rank
                receiveAndSend(
                        pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Center][gmx::DirectionY::Right],
                        overlapY * myGridX * localPmeSize[ZZ],
                        right,
                        &req[reqCount],
                        pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Center][gmx::DirectionY::Left],
                        overlapLeft * myGridX * localPmeSize[ZZ],
                        left,
                        &req[reqCount + 1],
                        mpiTag,
                        pmeGpu->common->mpiCommY);
                reqCount += 2;
            }
        }

        if (sizeX > 1 && sizeY > 1)
        {
            int rankUpLeft   = up * sizeY + left;
            int rankDownLeft = down * sizeY + left;

            int rankUpRight   = up * sizeY + right;
            int rankDownRight = down * sizeY + right;

            constexpr int mpiTag = 408; // Arbitrarily chosen

            // send data to up left and recv from down right rank
            receiveAndSend(
                    pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Up][gmx::DirectionY::Left],
                    overlapX * overlapY * localPmeSize[ZZ],
                    rankUpLeft,
                    &req[reqCount],
                    pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Down][gmx::DirectionY::Right],
                    overlapDown * overlapRight * localPmeSize[ZZ],
                    rankDownRight,
                    &req[reqCount + 1],
                    mpiTag,
                    pmeGpu->common->mpiComm);
            reqCount += 2;

            if (overlapLeft > 0)
            {
                // send data to up right rank and recv from down left rank
                receiveAndSend(
                        pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Up][gmx::DirectionY::Right],
                        overlapX * overlapY * localPmeSize[ZZ],
                        rankUpRight,
                        &req[reqCount],
                        pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Down][gmx::DirectionY::Left],
                        overlapDown * overlapLeft * localPmeSize[ZZ],
                        rankDownLeft,
                        &req[reqCount + 1],
                        mpiTag,
                        pmeGpu->common->mpiComm);
                reqCount += 2;
            }

            if (overlapUp > 0)
            {
                // send data to down left rank and recv from up right rank
                receiveAndSend(
                        pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Down][gmx::DirectionY::Left],
                        overlapX * overlapY * localPmeSize[ZZ],
                        rankDownLeft,
                        &req[reqCount],
                        pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Up][gmx::DirectionY::Right],
                        overlapUp * overlapRight * localPmeSize[ZZ],
                        rankUpRight,
                        &req[reqCount + 1],
                        mpiTag,
                        pmeGpu->common->mpiComm);
                reqCount += 2;
            }

            if (overlapUp > 0 && overlapLeft > 0)
            {
                // send data to down right rank and recv from up left rank
                receiveAndSend(
                        pmeGpu->haloExchange->d_sendGrids[gmx::DirectionX::Down][gmx::DirectionY::Right],
                        overlapX * overlapY * localPmeSize[ZZ],
                        rankDownRight,
                        &req[reqCount],
                        pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Up][gmx::DirectionY::Left],
                        overlapUp * overlapLeft * localPmeSize[ZZ],
                        rankUpLeft,
                        &req[reqCount + 1],
                        mpiTag,
                        pmeGpu->common->mpiComm);
                reqCount += 2;
            }
        }

        MPI_Waitall(reqCount, req, MPI_STATUSES_IGNORE);

        wallcycle_stop(wcycle, WallCycleCounter::PmeHaloExchangeComm);

        // data is written at the right place as part of MPI communication if slab-decomposition is
        // used in X-dimension, but we need to unpack if decomposition happens (also) along Y
        if (sizeY > 1)
        {
            wallcycle_start(wcycle, WallCycleCounter::LaunchGpuPme);

            // assign halo data
            unpackHaloDataExternal(
                    pmeGpu,
                    overlapUp,
                    overlapDown,
                    overlapLeft,
                    overlapRight,
                    myGridX,
                    myGridY,
                    localPmeSize,
                    realGrid,
                    recvGridUp,
                    recvGridDown,
                    pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Center][gmx::DirectionY::Left],
                    pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Center][gmx::DirectionY::Right],
                    pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Up][gmx::DirectionY::Left],
                    pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Down][gmx::DirectionY::Left],
                    pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Up][gmx::DirectionY::Right],
                    pmeGpu->haloExchange->d_recvGrids[gmx::DirectionX::Down][gmx::DirectionY::Right]);

            wallcycle_stop(wcycle, WallCycleCounter::LaunchGpuPme);
        }
    }
#else
    GMX_UNUSED_VALUE(pmeGpu);
#endif
}

template<bool pmeToFft>
void convertPmeGridToFftGrid(const PmeGpu* pmeGpu, float* h_fftRealGrid, gmx_parallel_3dfft* fftSetup, const int gridIndex)
{
    ivec localFftNData, localFftOffset, localFftSize;
    ivec localPmeSize;

    gmx_parallel_3dfft_real_limits(fftSetup, localFftNData, localFftOffset, localFftSize);

    localPmeSize[XX] = pmeGpu->kernelParams->grid.realGridSizePadded[XX];
    localPmeSize[YY] = pmeGpu->kernelParams->grid.realGridSizePadded[YY];
    localPmeSize[ZZ] = pmeGpu->kernelParams->grid.realGridSizePadded[ZZ];

    // this is true in case of slab decomposition
    if (localPmeSize[ZZ] == localFftSize[ZZ] && localPmeSize[YY] == localFftSize[YY])
    {
        int fftSize = localFftSize[ZZ] * localFftSize[YY] * localFftNData[XX];
        if (pmeToFft)
        {
            copyFromDeviceBuffer(h_fftRealGrid,
                                 &pmeGpu->kernelParams->grid.d_realGrid[gridIndex],
                                 0,
                                 fftSize,
                                 pmeGpu->archSpecific->pmeStream_,
                                 pmeGpu->settings.transferKind,
                                 nullptr);
        }
        else
        {
            copyToDeviceBuffer(&pmeGpu->kernelParams->grid.d_realGrid[gridIndex],
                               h_fftRealGrid,
                               0,
                               fftSize,
                               pmeGpu->archSpecific->pmeStream_,
                               pmeGpu->settings.transferKind,
                               nullptr);
        }
    }
    else
    {
        // launch copy kernel
        // Keeping threadsAlongZDim same as warp size for better coalescing,
        // Not keeping to higher value such as 64 to avoid high masked out
        // inactive threads as FFT grid sizes tend to be quite small
        const int threadsAlongZDim = 32;

        KernelLaunchConfig config;
        config.blockSize[0] = threadsAlongZDim;
        config.blockSize[1] = 4;
        config.blockSize[2] = 1;
        config.gridSize[0]  = gmx::divideRoundUp<size_t>(localFftNData[ZZ], config.blockSize[0]);
        config.gridSize[1]  = gmx::divideRoundUp<size_t>(localFftNData[YY], config.blockSize[1]);
        config.gridSize[2]  = localFftNData[XX];
        config.sharedMemorySize = 0;

        auto kernelFn = pmegrid_to_fftgrid<pmeToFft>;

        const auto kernelArgs =
                prepareGpuKernelArguments(kernelFn,
                                          config,
                                          &pmeGpu->kernelParams->grid.d_realGrid[gridIndex],
                                          &h_fftRealGrid,
                                          &localFftNData,
                                          &localFftSize,
                                          &localPmeSize);

        launchGpuKernel(kernelFn,
                        config,
                        pmeGpu->archSpecific->pmeStream_,
                        nullptr,
                        "Convert PME grid to FFT grid",
                        kernelArgs);
    }

    if (pmeToFft)
    {
        pmeGpu->archSpecific->syncSpreadGridD2H.markEvent(pmeGpu->archSpecific->pmeStream_);
    }
}

template<bool pmeToFft>
void convertPmeGridToFftGrid(const PmeGpu* pmeGpu, DeviceBuffer<float>* d_fftRealGrid, const int gridIndex)
{
    ivec localPmeSize;

    ivec localFftNData, localFftSize;

    localPmeSize[XX] = pmeGpu->kernelParams->grid.realGridSizePadded[XX];
    localPmeSize[YY] = pmeGpu->kernelParams->grid.realGridSizePadded[YY];
    localPmeSize[ZZ] = pmeGpu->kernelParams->grid.realGridSizePadded[ZZ];

    localFftNData[XX] = pmeGpu->archSpecific->localRealGridSize[XX];
    localFftNData[YY] = pmeGpu->archSpecific->localRealGridSize[YY];
    localFftNData[ZZ] = pmeGpu->archSpecific->localRealGridSize[ZZ];

    localFftSize[XX] = pmeGpu->archSpecific->localRealGridSizePadded[XX];
    localFftSize[YY] = pmeGpu->archSpecific->localRealGridSizePadded[YY];
    localFftSize[ZZ] = pmeGpu->archSpecific->localRealGridSizePadded[ZZ];

    // this is true in case of slab decomposition
    if (localPmeSize[ZZ] == localFftSize[ZZ] && localPmeSize[YY] == localFftSize[YY])
    {
        int fftSize = localFftSize[ZZ] * localFftSize[YY] * localFftNData[XX];
        if (pmeToFft)
        {
            copyBetweenDeviceBuffers(d_fftRealGrid,
                                     &pmeGpu->kernelParams->grid.d_realGrid[gridIndex],
                                     fftSize,
                                     pmeGpu->archSpecific->pmeStream_,
                                     pmeGpu->settings.transferKind,
                                     nullptr);
        }
        else
        {
            copyBetweenDeviceBuffers(&pmeGpu->kernelParams->grid.d_realGrid[gridIndex],
                                     d_fftRealGrid,
                                     fftSize,
                                     pmeGpu->archSpecific->pmeStream_,
                                     pmeGpu->settings.transferKind,
                                     nullptr);
        }
    }
    else
    {
        // launch copy kernel
        // keeping same as warp size for better coalescing
        // Not keeping to higher value such as 64 to avoid high masked out
        // inactive threads as FFT grid sizes tend to be quite small
        const int threadsAlongZDim = 32;

        KernelLaunchConfig config;
        config.blockSize[0] = threadsAlongZDim;
        config.blockSize[1] = 4;
        config.blockSize[2] = 1;
        config.gridSize[0]  = gmx::divideRoundUp<size_t>(localFftNData[ZZ], config.blockSize[0]);
        config.gridSize[1]  = gmx::divideRoundUp<size_t>(localFftNData[YY], config.blockSize[1]);
        config.gridSize[2]  = localFftNData[XX];
        config.sharedMemorySize = 0;

        auto kernelFn = pmegrid_to_fftgrid<pmeToFft>;

        const auto kernelArgs =
                prepareGpuKernelArguments(kernelFn,
                                          config,
                                          &pmeGpu->kernelParams->grid.d_realGrid[gridIndex],
                                          d_fftRealGrid,
                                          &localFftNData,
                                          &localFftSize,
                                          &localPmeSize);

        launchGpuKernel(kernelFn,
                        config,
                        pmeGpu->archSpecific->pmeStream_,
                        nullptr,
                        "Convert PME grid to FFT grid",
                        kernelArgs);
    }
}

template void convertPmeGridToFftGrid<true>(const PmeGpu*       pmeGpu,
                                            float*              h_fftRealGrid,
                                            gmx_parallel_3dfft* fftSetup,
                                            const int           gridIndex);

template void convertPmeGridToFftGrid<false>(const PmeGpu*       pmeGpu,
                                             float*              h_fftRealGrid,
                                             gmx_parallel_3dfft* fftSetup,
                                             const int           gridIndex);

template void convertPmeGridToFftGrid<true>(const PmeGpu*        pmeGpu,
                                            DeviceBuffer<float>* d_fftRealGrid,
                                            const int            gridIndex);

template void convertPmeGridToFftGrid<false>(const PmeGpu*        pmeGpu,
                                             DeviceBuffer<float>* d_fftRealGrid,
                                             const int            gridIndex);
