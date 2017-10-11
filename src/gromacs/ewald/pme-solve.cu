/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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
 *  \brief Implements PME GPU Fourier grid solving in CUDA.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 */

#include "gmxpre.h"

#include "config.h"

#include "gromacs/gpu_utils/cuda_arch_utils.cuh"
#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

#include "pme.cuh"
#include "pme-timings.cuh"

//! Solving kernel max block width in warps picked among powers of 2 (2, 4, 8, 16) for max. occupancy and min. runtime
//! (560Ti (CC2.1), 660Ti (CC3.0) and 750 (CC5.0)))
constexpr int c_solveMaxWarpsPerBlock = 8;
//! Solving kernel max block size in threads
constexpr int c_solveMaxThreadsPerBlock = (c_solveMaxWarpsPerBlock * warp_size);

// CUDA 6.5 can not compile enum class as a template kernel parameter,
// so we replace it with a duplicate simple enum
#if GMX_CUDA_VERSION >= 7000
using GridOrderingInternal = GridOrdering;
#else
enum GridOrderingInternal
{
    YZX,
    XYZ
};
#endif

/*! \brief
 * PME complex grid solver kernel function.
 *
 * \tparam[in] gridOrdering             Specifies the dimension ordering of the complex grid.
 * \tparam[in] computeEnergyAndVirial   Tells if the reciprocal energy and virial should be computed.
 * \param[in]  kernelParams             Input PME CUDA data in constant memory.
 */
template<
    GridOrderingInternal gridOrdering,
    bool computeEnergyAndVirial
    >
__launch_bounds__(c_solveMaxThreadsPerBlock)
__global__ void pme_solve_kernel(const struct pme_gpu_cuda_kernel_params_t kernelParams)
{
    /* This kernel supports 2 different grid dimension orderings: YZX and XYZ */
    int majorDim, middleDim, minorDim;
    switch (gridOrdering)
    {
        case GridOrderingInternal::YZX:
            majorDim  = YY;
            middleDim = ZZ;
            minorDim  = XX;
            break;

        case GridOrderingInternal::XYZ:
            majorDim  = XX;
            middleDim = YY;
            minorDim  = ZZ;
            break;

        default:
            assert(false);
    }

    /* Global memory pointers */
    const float * __restrict__ gm_splineValueMajor    = kernelParams.grid.d_splineModuli + kernelParams.grid.splineValuesOffset[majorDim];
    const float * __restrict__ gm_splineValueMiddle   = kernelParams.grid.d_splineModuli + kernelParams.grid.splineValuesOffset[middleDim];
    const float * __restrict__ gm_splineValueMinor    = kernelParams.grid.d_splineModuli + kernelParams.grid.splineValuesOffset[minorDim];
    float * __restrict__       gm_virialAndEnergy     = kernelParams.constants.d_virialAndEnergy;
    float2 * __restrict__      gm_grid                = (float2 *)kernelParams.grid.d_fourierGrid;

    /* Various grid sizes and indices */
    const int localOffsetMinor = 0, localOffsetMajor = 0, localOffsetMiddle = 0; //unused
    const int localSizeMinor   = kernelParams.grid.complexGridSizePadded[minorDim];
    const int localSizeMiddle  = kernelParams.grid.complexGridSizePadded[middleDim];
    const int localCountMiddle = kernelParams.grid.complexGridSize[middleDim];
    const int localCountMinor  = kernelParams.grid.complexGridSize[minorDim];
    const int nMajor           = kernelParams.grid.realGridSize[majorDim];
    const int nMiddle          = kernelParams.grid.realGridSize[middleDim];
    const int nMinor           = kernelParams.grid.realGridSize[minorDim];
    const int maxkMajor        = (nMajor + 1) / 2;  // X or Y
    const int maxkMiddle       = (nMiddle + 1) / 2; // Y OR Z => only check for !YZX
    const int maxkMinor        = (nMinor + 1) / 2;  // Z or X => only check for YZX

    /* Each thread works on one cell of the Fourier space complex 3D grid (gm_grid).
     * Each block handles up to c_solveMaxThreadsPerBlock cells -
     * depending on the grid contiguous dimension size,
     * that can range from a part of a single gridline to several complete gridlines.
     */
    const int threadLocalId     = threadIdx.x;
    const int gridLineSize      = localCountMinor;
    const int gridLineIndex     = threadLocalId / gridLineSize;
    const int gridLineCellIndex = threadLocalId - gridLineSize * gridLineIndex;
    const int gridLinesPerBlock = blockDim.x / gridLineSize;
    const int activeWarps       = (blockDim.x / warp_size);
    const int indexMinor        = blockIdx.x * blockDim.x + gridLineCellIndex;
    const int indexMiddle       = blockIdx.y * gridLinesPerBlock + gridLineIndex;
    const int indexMajor        = blockIdx.z;

    /* Optional outputs */
    float energy = 0.0f;
    float virxx  = 0.0f;
    float virxy  = 0.0f;
    float virxz  = 0.0f;
    float viryy  = 0.0f;
    float viryz  = 0.0f;
    float virzz  = 0.0f;

    assert(indexMajor < kernelParams.grid.complexGridSize[majorDim]);
    if ((indexMiddle < localCountMiddle) & (indexMinor < localCountMinor) & (gridLineIndex < gridLinesPerBlock))
    {
        /* The offset should be equal to the global thread index for coalesced access */
        const int             gridIndex     = (indexMajor * localSizeMiddle + indexMiddle) * localSizeMinor + indexMinor;
        float2 * __restrict__ gm_gridCell   = gm_grid + gridIndex;

        const int             kMajor  = indexMajor + localOffsetMajor;
        /* Checking either X in XYZ, or Y in YZX cases */
        const float           mMajor  = (kMajor < maxkMajor) ? kMajor : (kMajor - nMajor);

        const int             kMiddle = indexMiddle + localOffsetMiddle;
        float                 mMiddle = kMiddle;
        /* Checking Y in XYZ case */
        if (gridOrdering == GridOrderingInternal::XYZ)
        {
            mMiddle = (kMiddle < maxkMiddle) ? kMiddle : (kMiddle - nMiddle);
        }
        const int             kMinor  = localOffsetMinor + indexMinor;
        float                 mMinor  = kMinor;
        /* Checking X in YZX case */
        if (gridOrdering == GridOrderingInternal::YZX)
        {
            mMinor = (kMinor < maxkMinor) ? kMinor : (kMinor - nMinor);
        }
        /* We should skip the k-space point (0,0,0) */
        const bool notZeroPoint  = (kMinor > 0) | (kMajor > 0) | (kMiddle > 0);

        float      mX, mY, mZ;
        switch (gridOrdering)
        {
            case GridOrderingInternal::YZX:
                mX = mMinor;
                mY = mMajor;
                mZ = mMiddle;
                break;

            case GridOrderingInternal::XYZ:
                mX = mMajor;
                mY = mMiddle;
                mZ = mMinor;
                break;

            default:
                assert(false);
        }

        /* 0.5 correction factor for the first and last components of a Z dimension */
        float corner_fac = 1.0f;
        switch (gridOrdering)
        {
            case GridOrderingInternal::YZX:
                if ((kMiddle == 0) | (kMiddle == maxkMiddle))
                {
                    corner_fac = 0.5f;
                }
                break;

            case GridOrderingInternal::XYZ:
                if ((kMinor == 0) | (kMinor == maxkMinor))
                {
                    corner_fac = 0.5f;
                }
                break;

            default:
                assert(false);
        }

        if (notZeroPoint)
        {
            const float mhxk = mX * kernelParams.step.recipBox[XX][XX];
            const float mhyk = mX * kernelParams.step.recipBox[XX][YY] + mY * kernelParams.step.recipBox[YY][YY];
            const float mhzk = mX * kernelParams.step.recipBox[XX][ZZ] + mY * kernelParams.step.recipBox[YY][ZZ] + mZ * kernelParams.step.recipBox[ZZ][ZZ];

            const float m2k        = mhxk * mhxk + mhyk * mhyk + mhzk * mhzk;
            assert(m2k != 0.0f);
            //TODO: use LDG/textures for gm_splineValue
            float       denom = m2k * float(M_PI) * kernelParams.step.boxVolume * gm_splineValueMajor[kMajor] * gm_splineValueMiddle[kMiddle] * gm_splineValueMinor[kMinor];
            assert(isfinite(denom));
            assert(denom != 0.0f);
            const float   tmp1   = expf(-kernelParams.grid.ewaldFactor * m2k);
            const float   etermk = kernelParams.constants.elFactor * tmp1 / denom;

            float2        gridValue    = *gm_gridCell;
            const float2  oldGridValue = gridValue;
            gridValue.x   *= etermk;
            gridValue.y   *= etermk;
            *gm_gridCell   = gridValue;

            if (computeEnergyAndVirial)
            {
                const float tmp1k = 2.0f * (gridValue.x * oldGridValue.x + gridValue.y * oldGridValue.y);

                float       vfactor = (kernelParams.grid.ewaldFactor + 1.0f / m2k) * 2.0f;
                float       ets2    = corner_fac * tmp1k;
                energy = ets2;

                float ets2vf  = ets2 * vfactor;

                virxx   = ets2vf * mhxk * mhxk - ets2;
                virxy   = ets2vf * mhxk * mhyk;
                virxz   = ets2vf * mhxk * mhzk;
                viryy   = ets2vf * mhyk * mhyk - ets2;
                viryz   = ets2vf * mhyk * mhzk;
                virzz   = ets2vf * mhzk * mhzk - ets2;
            }
        }
    }

    /* Optional energy/virial reduction */
    if (computeEnergyAndVirial)
    {
#if (GMX_PTX_ARCH >= 300)
        /* A tricky shuffle reduction inspired by reduce_force_j_warp_shfl.
         * The idea is to reduce 7 energy/virial components into a single variable (aligned by 8).
         * We will reduce everything into virxx.
         */

        /* We can only reduce warp-wise */
        const int          width      = warp_size;
        const unsigned int activeMask = c_fullWarpMask;

        /* Making pair sums */
        virxx  += gmx_shfl_down_sync(activeMask, virxx, 1, width);
        viryy  += gmx_shfl_up_sync  (activeMask, viryy, 1, width);
        virzz  += gmx_shfl_down_sync(activeMask, virzz, 1, width);
        virxy  += gmx_shfl_up_sync  (activeMask, virxy, 1, width);
        virxz  += gmx_shfl_down_sync(activeMask, virxz, 1, width);
        viryz  += gmx_shfl_up_sync  (activeMask, viryz, 1, width);
        energy += gmx_shfl_down_sync(activeMask, energy, 1, width);
        if (threadLocalId & 1)
        {
            virxx = viryy; // virxx now holds virxx and viryy pair sums
            virzz = virxy; // virzz now holds virzz and virxy pair sums
            virxz = viryz; // virxz now holds virxz and viryz pair sums
        }

        /* Making quad sums */
        virxx  += gmx_shfl_down_sync(activeMask, virxx, 2, width);
        virzz  += gmx_shfl_up_sync  (activeMask, virzz, 2, width);
        virxz  += gmx_shfl_down_sync(activeMask, virxz, 2, width);
        energy += gmx_shfl_up_sync  (activeMask, energy, 2, width);
        if (threadLocalId & 2)
        {
            virxx = virzz;  // virxx now holds quad sums of virxx, virxy, virzz and virxy
            virxz = energy; // virxz now holds quad sums of virxz, viryz, energy and unused paddings
        }

        /* Making octet sums */
        virxx += gmx_shfl_down_sync(activeMask, virxx, 4, width);
        virxz += gmx_shfl_up_sync  (activeMask, virxz, 4, width);
        if (threadLocalId & 4)
        {
            virxx = virxz; // virxx now holds all 7 components' octet sums + unused paddings
        }

        /* We only need to reduce virxx now */
#pragma unroll
        for (int delta = 8; delta < width; delta <<= 1)
        {
            virxx += gmx_shfl_down_sync(activeMask, virxx, delta, width);
        }
        /* Now first 7 threads of each warp have the full output contributions in virxx */

        const int        componentIndex      = threadLocalId & (warp_size - 1);
        const bool       validComponentIndex = (componentIndex < c_virialAndEnergyCount);
        /* Reduce 7 outputs per warp in the shared memory */
        const int        stride              = 8; // this is c_virialAndEnergyCount==7 rounded up to power of 2 for convenience, hence the assert
        assert(c_virialAndEnergyCount == 7);
        const int        reductionBufferSize = (c_solveMaxThreadsPerBlock / warp_size) * stride;
        __shared__ float sm_virialAndEnergy[reductionBufferSize];

        if (validComponentIndex)
        {
            const int warpIndex = threadLocalId / warp_size;
            sm_virialAndEnergy[warpIndex * stride + componentIndex] = virxx;
        }
        __syncthreads();

        /* Reduce to the single warp size */
        const int targetIndex = threadLocalId;
#pragma unroll
        for (int reductionStride = reductionBufferSize >> 1; reductionStride >= warp_size; reductionStride >>= 1)
        {
            const int sourceIndex = targetIndex + reductionStride;
            if ((targetIndex < reductionStride) & (sourceIndex < activeWarps * stride))
            {
                // TODO: the second conditional is only needed on first iteration, actually - see if compiler eliminates it!
                sm_virialAndEnergy[targetIndex] += sm_virialAndEnergy[sourceIndex];
            }
            __syncthreads();
        }

        /* Now use shuffle again */
        if (threadLocalId < warp_size)
        {
            float output = sm_virialAndEnergy[threadLocalId];
#pragma unroll
            for (int delta = stride; delta < warp_size; delta <<= 1)
            {
                output += gmx_shfl_down_sync(activeMask, output, delta, warp_size);
            }
            /* Final output */
            if (validComponentIndex)
            {
                assert(isfinite(output));
                atomicAdd(gm_virialAndEnergy + componentIndex, output);
            }
        }
#else
        /* Shared memory reduction with atomics for compute capability < 3.0.
         * Each component is first reduced into warp_size positions in the shared memory;
         * Then first c_virialAndEnergyCount warps reduce everything further and add to the global memory.
         * This can likely be improved, but is anyway faster than the previous straightforward reduction,
         * which was using too much shared memory (for storing all 7 floats on each thread).
         * [48KB (shared mem limit per SM on CC2.x) / sizeof(float) (4) / c_solveMaxThreadsPerBlock (256) / c_virialAndEnergyCount (7) ==
         * 6 blocks per SM instead of 16 which is maximum on CC2.x].
         */

        const int        lane      = threadLocalId & (warp_size - 1);
        const int        warpIndex = threadLocalId / warp_size;
        const bool       firstWarp = (warpIndex == 0);
        __shared__ float sm_virialAndEnergy[c_virialAndEnergyCount * warp_size];
        if (firstWarp)
        {
            sm_virialAndEnergy[0 * warp_size + lane] = virxx;
            sm_virialAndEnergy[1 * warp_size + lane] = viryy;
            sm_virialAndEnergy[2 * warp_size + lane] = virzz;
            sm_virialAndEnergy[3 * warp_size + lane] = virxy;
            sm_virialAndEnergy[4 * warp_size + lane] = virxz;
            sm_virialAndEnergy[5 * warp_size + lane] = viryz;
            sm_virialAndEnergy[6 * warp_size + lane] = energy;
        }
        __syncthreads();
        if (!firstWarp)
        {
            atomicAdd(sm_virialAndEnergy + 0 * warp_size + lane, virxx);
            atomicAdd(sm_virialAndEnergy + 1 * warp_size + lane, viryy);
            atomicAdd(sm_virialAndEnergy + 2 * warp_size + lane, virzz);
            atomicAdd(sm_virialAndEnergy + 3 * warp_size + lane, virxy);
            atomicAdd(sm_virialAndEnergy + 4 * warp_size + lane, virxz);
            atomicAdd(sm_virialAndEnergy + 5 * warp_size + lane, viryz);
            atomicAdd(sm_virialAndEnergy + 6 * warp_size + lane, energy);
        }
        __syncthreads();

        GMX_UNUSED_VALUE(activeWarps);
        assert(activeWarps >= c_virialAndEnergyCount); // we need to cover all components, or have multiple iterations otherwise
        const int componentIndex = warpIndex;
        if (componentIndex < c_virialAndEnergyCount)
        {
            const int targetIndex = threadLocalId;
#pragma unroll
            for (int reductionStride = warp_size >> 1; reductionStride >= 1; reductionStride >>= 1)
            {
                if (lane < reductionStride)
                {
                    sm_virialAndEnergy[targetIndex] += sm_virialAndEnergy[targetIndex + reductionStride];
                }
            }
            if (lane == 0)
            {
                atomicAdd(gm_virialAndEnergy + componentIndex, sm_virialAndEnergy[targetIndex]);
            }
        }
#endif
    }
}

void pme_gpu_solve(const pme_gpu_t *pmeGpu, t_complex *h_grid,
                   GridOrdering gridOrdering, bool computeEnergyAndVirial)
{
    const bool   copyInputAndOutputGrid = pme_gpu_is_testing(pmeGpu) || !pme_gpu_performs_FFT(pmeGpu);

    cudaStream_t stream          = pmeGpu->archSpecific->pmeStream;
    const auto  *kernelParamsPtr = pmeGpu->kernelParams.get();

    if (copyInputAndOutputGrid)
    {
        cu_copy_H2D_async(kernelParamsPtr->grid.d_fourierGrid, h_grid, pmeGpu->archSpecific->complexGridSize * sizeof(float), stream);
    }

    int majorDim = -1, middleDim = -1, minorDim = -1;
    switch (gridOrdering)
    {
        case GridOrdering::YZX:
            majorDim  = YY;
            middleDim = ZZ;
            minorDim  = XX;
            break;

        case GridOrdering::XYZ:
            majorDim  = XX;
            middleDim = YY;
            minorDim  = ZZ;
            break;

        default:
            GMX_ASSERT(false, "Implement grid ordering here and below for the kernel launch");
    }

    const int maxBlockSize      = c_solveMaxThreadsPerBlock;
    const int gridLineSize      = pmeGpu->kernelParams->grid.complexGridSize[minorDim];
    const int gridLinesPerBlock = std::max(maxBlockSize / gridLineSize, 1);
    const int blocksPerGridLine = (gridLineSize + maxBlockSize - 1) / maxBlockSize;
    const int cellsPerBlock     = gridLineSize * gridLinesPerBlock;
    const int blockSize         = (cellsPerBlock + warp_size - 1) / warp_size * warp_size;
    // rounding up to full warps so that shuffle operations produce defined results
    dim3 threads(blockSize);
    dim3 blocks(blocksPerGridLine,
                (pmeGpu->kernelParams->grid.complexGridSize[middleDim] + gridLinesPerBlock - 1) / gridLinesPerBlock,
                pmeGpu->kernelParams->grid.complexGridSize[majorDim]);

    pme_gpu_start_timing(pmeGpu, gtPME_SOLVE);
    if (gridOrdering == GridOrdering::YZX)
    {
        if (computeEnergyAndVirial)
        {
            pme_solve_kernel<GridOrderingInternal::YZX, true> <<< blocks, threads, 0, stream>>> (*kernelParamsPtr);
        }
        else
        {
            pme_solve_kernel<GridOrderingInternal::YZX, false> <<< blocks, threads, 0, stream>>> (*kernelParamsPtr);
        }
    }
    else if (gridOrdering == GridOrdering::XYZ)
    {
        if (computeEnergyAndVirial)
        {
            pme_solve_kernel<GridOrderingInternal::XYZ, true> <<< blocks, threads, 0, stream>>> (*kernelParamsPtr);
        }
        else
        {
            pme_solve_kernel<GridOrderingInternal::XYZ, false> <<< blocks, threads, 0, stream>>> (*kernelParamsPtr);
        }
    }
    CU_LAUNCH_ERR("pme_solve_kernel");
    pme_gpu_stop_timing(pmeGpu, gtPME_SOLVE);

    if (computeEnergyAndVirial)
    {
        cu_copy_D2H_async(pmeGpu->staging.h_virialAndEnergy, kernelParamsPtr->constants.d_virialAndEnergy,
                          c_virialAndEnergyCount * sizeof(float), stream);
        cudaError_t stat = cudaEventRecord(pmeGpu->archSpecific->syncEnerVirD2H, stream);
        CU_RET_ERR(stat, "PME solve energy/virial event record failure");
    }

    if (copyInputAndOutputGrid)
    {
        cu_copy_D2H_async(h_grid, kernelParamsPtr->grid.d_fourierGrid, pmeGpu->archSpecific->complexGridSize * sizeof(float), stream);
        cudaError_t stat = cudaEventRecord(pmeGpu->archSpecific->syncSolveGridD2H, stream);
        CU_RET_ERR(stat, "PME solve grid sync event record failure");
    }
}
