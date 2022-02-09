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

#include "gromacs/math/vec.h"
#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/devicebuffer.cuh"

#include "pme.cuh"
#include "pme_gpu_types_host.h"
#include "pme_gpu_types.h"
#include "pme_gpu_types_host_impl.h"
#include "gromacs/fft/parallel_3dfft.h"

/*! \brief
 * A CUDA kernel which packs non-contiguous overlap data in Y-dimension
 *
 * \param[in] gm_realGrid          local grid
 * \param[in] gm_transferGrid      device array used to pack data
 * \param[in] offset               offset of y-overlap region
 * \param[in] overlapSize          overlap Size in y-overlap region
 * \param[in] pmeSize              Local PME grid size
 */
static __global__ void pmeGpuPackHaloY(const float* __restrict__ gm_realGrid,
                                       float* __restrict__ gm_transferGrid,
                                       int  offset,
                                       int  overlapSize,
                                       int3 pmeSize)
{
    int iz = threadIdx.x + blockIdx.x * blockDim.x;
    int iy = threadIdx.y + blockIdx.y * blockDim.y;
    int ix = threadIdx.z + blockIdx.z * blockDim.z;

    // we might get iz greather than pmeSize.z when pmeSize.z is not
    // multiple of threadsAlongZDim(see below)
    if (iz >= pmeSize.z)
    {
        return;
    }

    int pmeIndex    = ix * pmeSize.y * pmeSize.z + (iy + offset) * pmeSize.z + iz;
    int packedIndex = ix * overlapSize * pmeSize.z + iy * pmeSize.z + iz;

    gm_transferGrid[packedIndex] = gm_realGrid[pmeIndex];
}

/*! \brief
 * A CUDA kernel which adds/puts grid overlap data received from neighboring rank in Y-dim
 *
 * \param[in] gm_realGrid          local grid
 * \param[in] gm_transferGrid      overlapping region from neighboring rank
 * \param[in] starty               offset of y-overlap region
 * \param[in] overlapSize          overlap Size in y-overlap region
 * \param[in] pmeSize              Local PME grid size
 *
 * \tparam  reduce                A boolean which tells whether to reduce values or just assign
 */
template<bool reduce>
static __global__ void pmeGpuAddHaloY(float* __restrict__ gm_realGrid,
                                      const float* __restrict__ gm_transferGrid,
                                      int  offset,
                                      int  overlapSize,
                                      int3 pmeSize)
{
    int iz = threadIdx.x + blockIdx.x * blockDim.x;
    int iy = threadIdx.y + blockIdx.y * blockDim.y;
    int ix = threadIdx.z + blockIdx.z * blockDim.z;

    // we might get iz greather than pmeSize.z when pmeSize.z is not
    // multiple of threadsAlongZDim(see below)
    if (iz >= pmeSize.z)
    {
        return;
    }

    int pmeIndex    = ix * pmeSize.y * pmeSize.z + (iy + offset) * pmeSize.z + iz;
    int packedIndex = ix * overlapSize * pmeSize.z + iy * pmeSize.z + iz;

    if (reduce)
    {
        gm_realGrid[pmeIndex] += gm_transferGrid[packedIndex];
    }
    else
    {
        gm_realGrid[pmeIndex] = gm_transferGrid[packedIndex];
    }
}

/*! \brief
 * A CUDA kernel which adds grid overlap data received from neighboring rank
 *
 * \param[in] gm_realGrid          local grid
 * \param[in] gm_transferGrid      overlapping region from neighboring rank
 * \param[in] size                 Number of elements in overlap region
 */
static __global__ void pmeGpuAddHalo(float* __restrict__ gm_realGrid,
                                     const float* __restrict__ gm_transferGrid,
                                     int size)
{
    int val = threadIdx.x + blockIdx.x * blockDim.x;
    if (val < size)
    {
        gm_realGrid[val] += gm_transferGrid[val];
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
 * Launches CUDA kernel to pack non-contiguous overlap data in Y-dimension
 *
 * \param[in]  pmeGpu              The PME GPU structure.
 * \param[in] overlapSize          overlap Size in y-overlap region
 * \param[in] yOffset              offset of y-overlap region
 * \param[in] localXSize           Local x size
 * \param[in] pmeSize              PME grid size
 * \param[in] realGrid             local grid
 * \param[in] packrdGrid           device array used to pack data
 */
static void packYData(const PmeGpu* pmeGpu,
                      int           overlapSize,
                      int           yOffset,
                      int           localXSize,
                      const ivec&   pmeSize,
                      float*        realGrid,
                      float*        packrdGrid)
{
    // keeping same as warp size for better coalescing
    // Not keeping to higher value such as 64 to avoid high masked out
    // inactive threads as FFT grid sizes tend to be quite small
    const int threadsAlongZDim = 32;

    // right grid
    KernelLaunchConfig config;
    config.blockSize[0]     = threadsAlongZDim;
    config.blockSize[1]     = overlapSize;
    config.blockSize[2]     = 1;
    config.gridSize[0]      = (pmeSize[ZZ] + threadsAlongZDim - 1) / threadsAlongZDim;
    config.gridSize[1]      = 1;
    config.gridSize[2]      = localXSize;
    config.sharedMemorySize = 0;


    auto kernelFn = pmeGpuPackHaloY;

    auto kernelArgs = prepareGpuKernelArguments(
            kernelFn, config, &realGrid, &packrdGrid, &yOffset, &overlapSize, &pmeSize);

    launchGpuKernel(kernelFn,
                    config,
                    pmeGpu->archSpecific->pmeStream_,
                    nullptr,
                    "PME Domdec GPU Pack Grid Halo Exchange",
                    kernelArgs);
}

/*! \brief
 * Launches CUDA kernel to reduce/unpack overlap data in Y-dimension
 *
 * \param[in]  pmeGpu              The PME GPU structure.
 * \param[in] overlapSize          overlap Size in y-overlap region
 * \param[in] yOffset              offset of y-overlap region
 * \param[in] localXSize           Local x size
 * \param[in] pmeSize              PME grid size
 * \param[in] realGrid             local grid
 * \param[in] packrdGrid           device array used to pack data
 *
 * \tparam  reduce                A boolean which tells whether to reduce values or just assign
 */
template<bool reduce>
static void reduceYData(const PmeGpu* pmeGpu,
                        int           overlapSize,
                        int           yOffset,
                        int           localXSize,
                        const ivec&   pmeSize,
                        float*        realGrid,
                        float*        packrdGrid)
{
    // keeping same as warp size for better coalescing
    // Not keeping to higher value such as 64 to avoid high masked out
    // inactive threads as FFT grid sizes tend to be quite small
    const int threadsAlongZDim = 32;

    // right grid
    KernelLaunchConfig config;
    config.blockSize[0]     = threadsAlongZDim;
    config.blockSize[1]     = overlapSize;
    config.blockSize[2]     = 1;
    config.gridSize[0]      = (pmeSize[ZZ] + threadsAlongZDim - 1) / threadsAlongZDim;
    config.gridSize[1]      = 1;
    config.gridSize[2]      = localXSize;
    config.sharedMemorySize = 0;

    auto kernelFn = pmeGpuAddHaloY<reduce>;

    auto kernelArgs = prepareGpuKernelArguments(
            kernelFn, config, &realGrid, &packrdGrid, &yOffset, &overlapSize, &pmeSize);

    launchGpuKernel(kernelFn,
                    config,
                    pmeGpu->archSpecific->pmeStream_,
                    nullptr,
                    "PME Domdec GPU Pack Grid Halo Exchange",
                    kernelArgs);
}

/*! \brief
 * Launches CUDA kernel to reduce overlap data in X-dimension
 *
 * \param[in]  pmeGpu              The PME GPU structure.
 * \param[in] overlapSize          overlap Size in y-overlap region
 * \param[in] realGrid             local grid
 * \param[in] packrdGrid           device array used to pack data
 */
static void reduceXData(const PmeGpu* pmeGpu, int overlapSize, float* realGrid, float* packrdGrid)
{
    // launch reduction kernel
    // ToDo: Experiment with different block size and decide
    const int threadsPerBlock = 256;

    KernelLaunchConfig config;
    config.blockSize[0]     = threadsPerBlock;
    config.blockSize[1]     = 1;
    config.blockSize[2]     = 1;
    config.gridSize[0]      = (overlapSize + threadsPerBlock - 1) / threadsPerBlock;
    config.gridSize[1]      = 1;
    config.gridSize[2]      = 1;
    config.sharedMemorySize = 0;

    auto kernelFn = pmeGpuAddHalo;

    auto kernelArgs = prepareGpuKernelArguments(kernelFn, config, &realGrid, &packrdGrid, &overlapSize);

    launchGpuKernel(kernelFn,
                    config,
                    pmeGpu->archSpecific->pmeStream_,
                    nullptr,
                    "PME Domdec GPU Apply Grid Halo Exchange",
                    kernelArgs);
}

void pmeGpuGridHaloExchange(const PmeGpu* pmeGpu)
{
#if GMX_MPI
    // Note here we are assuming that width of the chunks is not so small that we need to
    // transfer to/from multiple ranks i.e. that the distributed grid contains chunks at least order-1 points wide.

    auto* kernelParamsPtr = pmeGpu->kernelParams.get();
    ivec  localPmeSize;
    localPmeSize[XX] = kernelParamsPtr->grid.realGridSizePadded[XX];
    localPmeSize[YY] = kernelParamsPtr->grid.realGridSizePadded[YY];
    localPmeSize[ZZ] = kernelParamsPtr->grid.realGridSizePadded[ZZ];

    int overlapSize = pmeGpu->common->gridHalo;

    // minor dimension
    if (pmeGpu->common->nnodesY > 1)
    {
        int rank  = pmeGpu->common->nodeidY;
        int size  = pmeGpu->common->nnodesY;
        int right = (rank + 1) % size;
        int left  = (rank + size - 1) % size;

        // Note that s2g0[size] is the grid size (array is allocated to size+1)
        int myGrid    = pmeGpu->common->s2g0Y[rank + 1] - pmeGpu->common->s2g0Y[rank];
        int rightGrid = pmeGpu->common->s2g0Y[right + 1] - pmeGpu->common->s2g0Y[right];
        int leftGrid  = pmeGpu->common->s2g0Y[left + 1] - pmeGpu->common->s2g0Y[left];

        // current implementation transfers from/to only immediate neighbours
        GMX_ASSERT(overlapSize <= myGrid && overlapSize <= rightGrid && overlapSize <= leftGrid,
                   "Exchange supported only with immediate neighbor");

        int overlapRecv  = std::min(overlapSize, myGrid);
        int overlapRight = std::min(overlapSize, rightGrid);
        int overlapLeft  = std::min(overlapSize, leftGrid);

        // if only 2 PME ranks in Y-domain and overlap width more than slab width
        // just transfer all grid points from neighbor
        if (right == left && overlapRight + overlapLeft >= rightGrid)
        {
            overlapRecv  = myGrid;
            overlapRight = rightGrid;
            overlapLeft  = 0;
        }

        int pmegridNx = pmeGpu->common->pmegridNk[XX];

        for (int gridIndex = 0; gridIndex < pmeGpu->common->ngrids; gridIndex++)
        {
            // launch packing kernel
            float* realGrid = pmeGpu->kernelParams->grid.d_realGrid[gridIndex];

            // Pack data that needs to be sent to right rank
            packYData(pmeGpu,
                      overlapRight,
                      myGrid,
                      pmegridNx,
                      localPmeSize,
                      realGrid,
                      pmeGpu->archSpecific->d_sendGridRightY);

            if (overlapLeft > 0)
            {
                // Pack data that needs to be sent to left rank
                packYData(pmeGpu,
                          overlapLeft,
                          localPmeSize[YY] - overlapLeft,
                          pmegridNx,
                          localPmeSize,
                          realGrid,
                          pmeGpu->archSpecific->d_sendGridLeftY);
            }

            // synchronize before starting halo exchange
            pme_gpu_synchronize(pmeGpu);

            constexpr int mpiTag = 403; // Arbitrarily chosen

            // send data to right rank and recv from left rank
            MPI_Sendrecv(pmeGpu->archSpecific->d_sendGridRightY,
                         overlapRight * pmegridNx * localPmeSize[ZZ],
                         MPI_FLOAT,
                         right,
                         mpiTag,
                         pmeGpu->archSpecific->d_recvGridLeftY,
                         overlapRecv * pmegridNx * localPmeSize[ZZ],
                         MPI_FLOAT,
                         left,
                         mpiTag,
                         pmeGpu->common->mpiCommY,
                         MPI_STATUS_IGNORE);

            if (overlapLeft > 0)
            {
                // send data to left rank and recv from right rank
                MPI_Sendrecv(pmeGpu->archSpecific->d_sendGridLeftY,
                             overlapLeft * pmegridNx * localPmeSize[ZZ],
                             MPI_FLOAT,
                             left,
                             mpiTag,
                             pmeGpu->archSpecific->d_recvGridRightY,
                             overlapRecv * pmegridNx * localPmeSize[ZZ],
                             MPI_FLOAT,
                             right,
                             mpiTag,
                             pmeGpu->common->mpiCommY,
                             MPI_STATUS_IGNORE);
            }

            // reduce data received from left rank
            reduceYData<true>(
                    pmeGpu, overlapRecv, 0, pmegridNx, localPmeSize, realGrid, pmeGpu->archSpecific->d_recvGridLeftY);

            if (overlapLeft > 0)
            {
                // reduce data received from right rank
                reduceYData<true>(pmeGpu,
                                  overlapRecv,
                                  myGrid - overlapRecv,
                                  pmegridNx,
                                  localPmeSize,
                                  realGrid,
                                  pmeGpu->archSpecific->d_recvGridRightY);
            }
        }
    }

    // wait for spread to finish before starting halo exchange
    pmeGpu->archSpecific->spreadCompleted.waitForEvent();

    // major dimension
    if (pmeGpu->common->nnodesX > 1)
    {
        int rank  = pmeGpu->common->nodeidX;
        int size  = pmeGpu->common->nnodesX;
        int right = (rank + 1) % size;
        int left  = (rank + size - 1) % size;

        // Note that s2g0[size] is the grid size (array is allocated to size+1)
        int myGrid    = pmeGpu->common->s2g0X[rank + 1] - pmeGpu->common->s2g0X[rank];
        int rightGrid = pmeGpu->common->s2g0X[right + 1] - pmeGpu->common->s2g0X[right];
        int leftGrid  = pmeGpu->common->s2g0X[left + 1] - pmeGpu->common->s2g0X[left];

        // current implementation transfers from/to only immediate neighbours
        GMX_ASSERT(overlapSize <= myGrid && overlapSize <= rightGrid && overlapSize <= leftGrid,
                   "Exchange supported only with immediate neighbor");

        int overlapRecv  = std::min(overlapSize, myGrid);
        int overlapRight = std::min(overlapSize, rightGrid);
        int overlapLeft  = std::min(overlapSize, leftGrid);

        // if only 2 PME ranks in X-domain and overlap width more than slab width
        // just transfer all grid points from neighbor
        if (right == left && overlapRight + overlapLeft >= rightGrid)
        {
            overlapRecv  = myGrid;
            overlapRight = rightGrid;
            overlapLeft  = 0;
        }

        int transferStartRight = myGrid * localPmeSize[YY] * localPmeSize[ZZ];
        int transferStartLeft = (localPmeSize[XX] - overlapLeft) * localPmeSize[YY] * localPmeSize[ZZ];

        // Current implementation transfers the whole grid along y, an optimization is
        // possible where only local y-length can be transferred
        // But, this will require executing packing kernel
        int transferSizeSendRight = overlapRight * localPmeSize[YY] * localPmeSize[ZZ];
        int transferSizeSendLeft  = overlapLeft * localPmeSize[YY] * localPmeSize[ZZ];
        int transferSizeRecv      = overlapRecv * localPmeSize[YY] * localPmeSize[ZZ];

        for (int gridIndex = 0; gridIndex < pmeGpu->common->ngrids; gridIndex++)
        {
            float* realGrid = pmeGpu->kernelParams->grid.d_realGrid[gridIndex];

            constexpr int mpiTag = 403; // Arbitrarily chosen

            // send data to right rank and recv from left rank
            MPI_Sendrecv(&realGrid[transferStartRight],
                         transferSizeSendRight,
                         MPI_FLOAT,
                         right,
                         mpiTag,
                         pmeGpu->archSpecific->d_recvGridLeftX,
                         transferSizeRecv,
                         MPI_FLOAT,
                         left,
                         mpiTag,
                         pmeGpu->common->mpiCommX,
                         MPI_STATUS_IGNORE);

            if (overlapLeft > 0)
            {
                // send data to left rank and recv from right rank
                MPI_Sendrecv(&realGrid[transferStartLeft],
                             transferSizeSendLeft,
                             MPI_FLOAT,
                             left,
                             mpiTag,
                             pmeGpu->archSpecific->d_recvGridRightX,
                             transferSizeRecv,
                             MPI_FLOAT,
                             right,
                             mpiTag,
                             pmeGpu->common->mpiCommX,
                             MPI_STATUS_IGNORE);
            }

            // reduce data received from left rank
            reduceXData(pmeGpu, transferSizeRecv, realGrid, pmeGpu->archSpecific->d_recvGridLeftX);

            if (overlapLeft > 0)
            {
                // reduce data received from right rank
                int    offset       = (myGrid - overlapRecv) * localPmeSize[YY] * localPmeSize[ZZ];
                float* offsetedGrid = realGrid + offset;
                reduceXData(pmeGpu, transferSizeRecv, offsetedGrid, pmeGpu->archSpecific->d_recvGridRightX);
            }
        }
    }
#else
    GMX_UNUSED_VALUE(packYData);
    GMX_UNUSED_VALUE(pmeGpu);
    GMX_UNUSED_VALUE(reduceXData);
#endif
}

void pmeGpuGridHaloExchangeReverse(const PmeGpu* pmeGpu)
{
#if GMX_MPI
    auto* kernelParamsPtr = pmeGpu->kernelParams.get();
    ivec  localPmeSize;
    localPmeSize[XX] = kernelParamsPtr->grid.realGridSizePadded[XX];
    localPmeSize[YY] = kernelParamsPtr->grid.realGridSizePadded[YY];
    localPmeSize[ZZ] = kernelParamsPtr->grid.realGridSizePadded[ZZ];

    int overlapSize = pmeGpu->common->gridHalo;

    // minor dimension
    if (pmeGpu->common->nnodesY > 1)
    {
        int rank  = pmeGpu->common->nodeidY;
        int size  = pmeGpu->common->nnodesY;
        int right = (rank + 1) % size;
        int left  = (rank + size - 1) % size;

        int myGrid    = pmeGpu->common->s2g0Y[rank + 1] - pmeGpu->common->s2g0Y[rank];
        int rightGrid = pmeGpu->common->s2g0Y[right + 1] - pmeGpu->common->s2g0Y[right];
        int leftGrid  = pmeGpu->common->s2g0Y[left + 1] - pmeGpu->common->s2g0Y[left];

        // current implementation transfers from/to only immediate neighbours
        GMX_ASSERT(overlapSize <= myGrid && overlapSize <= rightGrid && overlapSize <= leftGrid,
                   "Exchange supported only with immediate neighbor");
        int overlapSend  = std::min(overlapSize, myGrid);
        int overlapRight = std::min(overlapSize, rightGrid);
        int overlapLeft  = std::min(overlapSize, leftGrid);

        // if only 2 PME ranks in Y-domain and overlap width more than slab width
        // just transfer all grid points from neighbor
        if (right == left && overlapRight + overlapLeft >= rightGrid)
        {
            overlapSend  = myGrid;
            overlapRight = rightGrid;
            overlapLeft  = 0;
        }

        int pmegridNx = pmeGpu->common->pmegridNk[XX];

        for (int gridIndex = 0; gridIndex < pmeGpu->common->ngrids; gridIndex++)
        {
            // launch packing kernel
            float* realGrid = pmeGpu->kernelParams->grid.d_realGrid[gridIndex];

            // Pack data that needs to be sent to left rank
            packYData(
                    pmeGpu, overlapSend, 0, pmegridNx, localPmeSize, realGrid, pmeGpu->archSpecific->d_sendGridLeftY);

            if (overlapLeft > 0)
            {
                // Pack data that needs to be sent to right rank
                packYData(pmeGpu,
                          overlapSend,
                          (myGrid - overlapSend),
                          pmegridNx,
                          localPmeSize,
                          realGrid,
                          pmeGpu->archSpecific->d_sendGridRightY);
            }

            // synchronize before starting halo exchange
            pme_gpu_synchronize(pmeGpu);

            constexpr int mpiTag = 403; // Arbitrarily chosen

            // send data to left rank and recv from right rank
            MPI_Sendrecv(pmeGpu->archSpecific->d_sendGridLeftY,
                         overlapSend * pmegridNx * localPmeSize[ZZ],
                         MPI_FLOAT,
                         left,
                         mpiTag,
                         pmeGpu->archSpecific->d_recvGridRightY,
                         overlapRight * pmegridNx * localPmeSize[ZZ],
                         MPI_FLOAT,
                         right,
                         mpiTag,
                         pmeGpu->common->mpiCommY,
                         MPI_STATUS_IGNORE);

            if (overlapLeft > 0)
            {
                // send data to right rank and recv from left rank
                MPI_Sendrecv(pmeGpu->archSpecific->d_sendGridRightY,
                             overlapSend * pmegridNx * localPmeSize[ZZ],
                             MPI_FLOAT,
                             right,
                             mpiTag,
                             pmeGpu->archSpecific->d_recvGridLeftY,
                             overlapLeft * pmegridNx * localPmeSize[ZZ],
                             MPI_FLOAT,
                             left,
                             mpiTag,
                             pmeGpu->common->mpiCommY,
                             MPI_STATUS_IGNORE);
            }

            // unpack data received from right rank
            reduceYData<false>(pmeGpu,
                               overlapRight,
                               myGrid,
                               pmegridNx,
                               localPmeSize,
                               realGrid,
                               pmeGpu->archSpecific->d_recvGridRightY);

            if (overlapLeft > 0)
            {
                // unpack data received from left rank
                reduceYData<false>(pmeGpu,
                                   overlapLeft,
                                   localPmeSize[YY] - overlapLeft,
                                   pmegridNx,
                                   localPmeSize,
                                   realGrid,
                                   pmeGpu->archSpecific->d_recvGridLeftY);
            }
        }
    }

    // Wait for conversion from FFT to Pme grid to finish before reverse halo exchange
    pmeGpu->archSpecific->syncFftToPmeGrid.waitForEvent();

    // major dimension
    if (pmeGpu->common->nnodesX > 1)
    {
        int rank  = pmeGpu->common->nodeidX;
        int size  = pmeGpu->common->nnodesX;
        int right = (rank + 1) % size;
        int left  = (rank + size - 1) % size;

        int myGrid    = pmeGpu->common->s2g0X[rank + 1] - pmeGpu->common->s2g0X[rank];
        int rightGrid = pmeGpu->common->s2g0X[right + 1] - pmeGpu->common->s2g0X[right];
        int leftGrid  = pmeGpu->common->s2g0X[left + 1] - pmeGpu->common->s2g0X[left];

        // current implementation transfers from/to only immediate neighbours
        GMX_ASSERT(overlapSize <= myGrid && overlapSize <= rightGrid && overlapSize <= leftGrid,
                   "Exchange supported only with immediate neighbor");
        int overlapSend  = std::min(overlapSize, myGrid);
        int overlapRight = std::min(overlapSize, rightGrid);
        int overlapLeft  = std::min(overlapSize, leftGrid);

        // if only 2 PME ranks in X-domain and overlap width more than slab width
        // just transfer all grid points from neighbor
        if (right == left && overlapRight + overlapLeft >= rightGrid)
        {
            overlapSend  = myGrid;
            overlapRight = rightGrid;
            overlapLeft  = 0;
        }

        int transferstartRight = myGrid * localPmeSize[YY] * localPmeSize[ZZ];
        int transferstartLeft = (localPmeSize[XX] - overlapLeft) * localPmeSize[YY] * localPmeSize[ZZ];

        // Current implementation transfers the whole grid along y, an optimization is
        // possible where only local y-length can be trasnferred
        // But, this will require executing packing kernel
        int transferSizeSend      = overlapSend * localPmeSize[YY] * localPmeSize[ZZ];
        int transferSizeRecvRight = overlapRight * localPmeSize[YY] * localPmeSize[ZZ];
        int transferSizeRecvLeft  = overlapLeft * localPmeSize[YY] * localPmeSize[ZZ];

        for (int gridIndex = 0; gridIndex < pmeGpu->common->ngrids; gridIndex++)
        {
            float* realGrid = pmeGpu->kernelParams->grid.d_realGrid[gridIndex];

            constexpr int mpiTag = 403; // Arbitrarily chosen

            // send data to left rank and recv from right rank
            MPI_Sendrecv(&realGrid[0],
                         transferSizeSend,
                         MPI_FLOAT,
                         left,
                         mpiTag,
                         &realGrid[transferstartRight],
                         transferSizeRecvRight,
                         MPI_FLOAT,
                         right,
                         mpiTag,
                         pmeGpu->common->mpiCommX,
                         MPI_STATUS_IGNORE);

            if (overlapLeft > 0)
            {
                // send data to right rank and recv from left rank
                int offset = (myGrid - overlapSend) * localPmeSize[YY] * localPmeSize[ZZ];
                MPI_Sendrecv(&realGrid[offset],
                             transferSizeSend,
                             MPI_FLOAT,
                             right,
                             mpiTag,
                             &realGrid[transferstartLeft],
                             transferSizeRecvLeft,
                             MPI_FLOAT,
                             left,
                             mpiTag,
                             pmeGpu->common->mpiCommX,
                             MPI_STATUS_IGNORE);
            }
        }
    }
#else
    GMX_UNUSED_VALUE(pmeGpu);
#endif
}

template<bool pmeToFft>
void convertPmeGridToFftGrid(const PmeGpu* pmeGpu, float* h_grid, gmx_parallel_3dfft_t* fftSetup, const int gridIndex)
{
    ivec localFftNData, localFftOffset, localFftSize;
    ivec localPmeSize;

    gmx_parallel_3dfft_real_limits(fftSetup[gridIndex], localFftNData, localFftOffset, localFftSize);

    localPmeSize[XX] = pmeGpu->kernelParams->grid.realGridSizePadded[XX];
    localPmeSize[YY] = pmeGpu->kernelParams->grid.realGridSizePadded[YY];
    localPmeSize[ZZ] = pmeGpu->kernelParams->grid.realGridSizePadded[ZZ];

    // this is true in case of slab decomposition
    if (localPmeSize[ZZ] == localFftSize[ZZ] && localPmeSize[YY] == localFftSize[YY])
    {
        int fftSize = localFftSize[ZZ] * localFftSize[YY] * localFftNData[XX];
        if (pmeToFft)
        {
            copyFromDeviceBuffer(h_grid,
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
                               h_grid,
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
        // ToDo: Experiment with different block size and decide on optimal configuration

        // keeping same as warp size for better coalescing
        // Not keeping to higher value such as 64 to avoid high masked out
        // inactive threads as FFT grid sizes tend to be quite small
        const int threadsAlongZDim = 32;

        KernelLaunchConfig config;
        config.blockSize[0] = threadsAlongZDim;
        config.blockSize[1] = 4;
        config.blockSize[2] = 1;
        config.gridSize[0]  = (localFftNData[ZZ] + config.blockSize[0] - 1) / config.blockSize[0];
        config.gridSize[1]  = (localFftNData[YY] + config.blockSize[1] - 1) / config.blockSize[1];
        config.gridSize[2]  = localFftNData[XX];
        config.sharedMemorySize = 0;

        auto kernelFn = pmegrid_to_fftgrid<pmeToFft>;

        const auto kernelArgs =
                prepareGpuKernelArguments(kernelFn,
                                          config,
                                          &pmeGpu->kernelParams->grid.d_realGrid[gridIndex],
                                          &h_grid,
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
    else
    {
        pmeGpu->archSpecific->syncFftToPmeGrid.markEvent(pmeGpu->archSpecific->pmeStream_);
    }
}

template void convertPmeGridToFftGrid<true>(const PmeGpu*         pmeGpu,
                                            float*                h_grid,
                                            gmx_parallel_3dfft_t* fftSetup,
                                            const int             gridIndex);

template void convertPmeGridToFftGrid<false>(const PmeGpu*         pmeGpu,
                                             float*                h_grid,
                                             gmx_parallel_3dfft_t* fftSetup,
                                             const int             gridIndex);
