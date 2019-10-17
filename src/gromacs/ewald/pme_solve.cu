/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018,2019, by the GROMACS development team, led by
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

#include <cassert>

#include <math_constants.h>

#include "gromacs/gpu_utils/cuda_arch_utils.cuh"

#include "pme.cuh"

/*! \brief
 * PME complex grid solver kernel function.
 *
 * \tparam[in] gridOrdering             Specifies the dimension ordering of the complex grid.
 * \tparam[in] computeEnergyAndVirial   Tells if the reciprocal energy and virial should be
 * computed. \param[in]  kernelParams             Input PME CUDA data in constant memory.
 */
template<GridOrdering gridOrdering, bool computeEnergyAndVirial>
__launch_bounds__(c_solveMaxThreadsPerBlock) CLANG_DISABLE_OPTIMIZATION_ATTRIBUTE __global__
        void pme_solve_kernel(const struct PmeGpuCudaKernelParams kernelParams)
{
    /* This kernel supports 2 different grid dimension orderings: YZX and XYZ */
    int majorDim, middleDim, minorDim;
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

        default: assert(false);
    }

    /* Global memory pointers */
    const float* __restrict__ gm_splineValueMajor =
            kernelParams.grid.d_splineModuli + kernelParams.grid.splineValuesOffset[majorDim];
    const float* __restrict__ gm_splineValueMiddle =
            kernelParams.grid.d_splineModuli + kernelParams.grid.splineValuesOffset[middleDim];
    const float* __restrict__ gm_splineValueMinor =
            kernelParams.grid.d_splineModuli + kernelParams.grid.splineValuesOffset[minorDim];
    float* __restrict__ gm_virialAndEnergy = kernelParams.constants.d_virialAndEnergy;
    float2* __restrict__ gm_grid           = (float2*)kernelParams.grid.d_fourierGrid;

    /* Various grid sizes and indices */
    const int localOffsetMinor = 0, localOffsetMajor = 0, localOffsetMiddle = 0; // unused
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
    const int gridLinesPerBlock = max(blockDim.x / gridLineSize, 1);
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
    if ((indexMiddle < localCountMiddle) & (indexMinor < localCountMinor)
        & (gridLineIndex < gridLinesPerBlock))
    {
        /* The offset should be equal to the global thread index for coalesced access */
        const int gridIndex = (indexMajor * localSizeMiddle + indexMiddle) * localSizeMinor + indexMinor;
        float2* __restrict__ gm_gridCell = gm_grid + gridIndex;

        const int kMajor = indexMajor + localOffsetMajor;
        /* Checking either X in XYZ, or Y in YZX cases */
        const float mMajor = (kMajor < maxkMajor) ? kMajor : (kMajor - nMajor);

        const int kMiddle = indexMiddle + localOffsetMiddle;
        float     mMiddle = kMiddle;
        /* Checking Y in XYZ case */
        if (gridOrdering == GridOrdering::XYZ)
        {
            mMiddle = (kMiddle < maxkMiddle) ? kMiddle : (kMiddle - nMiddle);
        }
        const int kMinor = localOffsetMinor + indexMinor;
        float     mMinor = kMinor;
        /* Checking X in YZX case */
        if (gridOrdering == GridOrdering::YZX)
        {
            mMinor = (kMinor < maxkMinor) ? kMinor : (kMinor - nMinor);
        }
        /* We should skip the k-space point (0,0,0) */
        const bool notZeroPoint = (kMinor > 0) | (kMajor > 0) | (kMiddle > 0);

        float mX, mY, mZ;
        switch (gridOrdering)
        {
            case GridOrdering::YZX:
                mX = mMinor;
                mY = mMajor;
                mZ = mMiddle;
                break;

            case GridOrdering::XYZ:
                mX = mMajor;
                mY = mMiddle;
                mZ = mMinor;
                break;

            default: assert(false);
        }

        /* 0.5 correction factor for the first and last components of a Z dimension */
        float corner_fac = 1.0f;
        switch (gridOrdering)
        {
            case GridOrdering::YZX:
                if ((kMiddle == 0) | (kMiddle == maxkMiddle))
                {
                    corner_fac = 0.5f;
                }
                break;

            case GridOrdering::XYZ:
                if ((kMinor == 0) | (kMinor == maxkMinor))
                {
                    corner_fac = 0.5f;
                }
                break;

            default: assert(false);
        }

        if (notZeroPoint)
        {
            const float mhxk = mX * kernelParams.current.recipBox[XX][XX];
            const float mhyk = mX * kernelParams.current.recipBox[XX][YY]
                               + mY * kernelParams.current.recipBox[YY][YY];
            const float mhzk = mX * kernelParams.current.recipBox[XX][ZZ]
                               + mY * kernelParams.current.recipBox[YY][ZZ]
                               + mZ * kernelParams.current.recipBox[ZZ][ZZ];

            const float m2k = mhxk * mhxk + mhyk * mhyk + mhzk * mhzk;
            assert(m2k != 0.0f);
            // TODO: use LDG/textures for gm_splineValue
            float denom = m2k * float(CUDART_PI_F) * kernelParams.current.boxVolume
                          * gm_splineValueMajor[kMajor] * gm_splineValueMiddle[kMiddle]
                          * gm_splineValueMinor[kMinor];
            assert(isfinite(denom));
            assert(denom != 0.0f);

            const float tmp1   = expf(-kernelParams.grid.ewaldFactor * m2k);
            const float etermk = kernelParams.constants.elFactor * tmp1 / denom;

            float2       gridValue    = *gm_gridCell;
            const float2 oldGridValue = gridValue;
            gridValue.x *= etermk;
            gridValue.y *= etermk;
            *gm_gridCell = gridValue;

            if (computeEnergyAndVirial)
            {
                const float tmp1k =
                        2.0f * (gridValue.x * oldGridValue.x + gridValue.y * oldGridValue.y);

                float vfactor = (kernelParams.grid.ewaldFactor + 1.0f / m2k) * 2.0f;
                float ets2    = corner_fac * tmp1k;
                energy        = ets2;

                float ets2vf = ets2 * vfactor;

                virxx = ets2vf * mhxk * mhxk - ets2;
                virxy = ets2vf * mhxk * mhyk;
                virxz = ets2vf * mhxk * mhzk;
                viryy = ets2vf * mhyk * mhyk - ets2;
                viryz = ets2vf * mhyk * mhzk;
                virzz = ets2vf * mhzk * mhzk - ets2;
            }
        }
    }

    /* Optional energy/virial reduction */
    if (computeEnergyAndVirial)
    {
        /* A tricky shuffle reduction inspired by reduce_force_j_warp_shfl.
         * The idea is to reduce 7 energy/virial components into a single variable (aligned by 8).
         * We will reduce everything into virxx.
         */

        /* We can only reduce warp-wise */
        const int          width      = warp_size;
        const unsigned int activeMask = c_fullWarpMask;

        /* Making pair sums */
        virxx += __shfl_down_sync(activeMask, virxx, 1, width);
        viryy += __shfl_up_sync(activeMask, viryy, 1, width);
        virzz += __shfl_down_sync(activeMask, virzz, 1, width);
        virxy += __shfl_up_sync(activeMask, virxy, 1, width);
        virxz += __shfl_down_sync(activeMask, virxz, 1, width);
        viryz += __shfl_up_sync(activeMask, viryz, 1, width);
        energy += __shfl_down_sync(activeMask, energy, 1, width);
        if (threadLocalId & 1)
        {
            virxx = viryy; // virxx now holds virxx and viryy pair sums
            virzz = virxy; // virzz now holds virzz and virxy pair sums
            virxz = viryz; // virxz now holds virxz and viryz pair sums
        }

        /* Making quad sums */
        virxx += __shfl_down_sync(activeMask, virxx, 2, width);
        virzz += __shfl_up_sync(activeMask, virzz, 2, width);
        virxz += __shfl_down_sync(activeMask, virxz, 2, width);
        energy += __shfl_up_sync(activeMask, energy, 2, width);
        if (threadLocalId & 2)
        {
            virxx = virzz;  // virxx now holds quad sums of virxx, virxy, virzz and virxy
            virxz = energy; // virxz now holds quad sums of virxz, viryz, energy and unused paddings
        }

        /* Making octet sums */
        virxx += __shfl_down_sync(activeMask, virxx, 4, width);
        virxz += __shfl_up_sync(activeMask, virxz, 4, width);
        if (threadLocalId & 4)
        {
            virxx = virxz; // virxx now holds all 7 components' octet sums + unused paddings
        }

        /* We only need to reduce virxx now */
#pragma unroll
        for (int delta = 8; delta < width; delta <<= 1)
        {
            virxx += __shfl_down_sync(activeMask, virxx, delta, width);
        }
        /* Now first 7 threads of each warp have the full output contributions in virxx */

        const int  componentIndex      = threadLocalId & (warp_size - 1);
        const bool validComponentIndex = (componentIndex < c_virialAndEnergyCount);
        /* Reduce 7 outputs per warp in the shared memory */
        const int stride =
                8; // this is c_virialAndEnergyCount==7 rounded up to power of 2 for convenience, hence the assert
        assert(c_virialAndEnergyCount == 7);
        const int        reductionBufferSize = (c_solveMaxThreadsPerBlock / warp_size) * stride;
        __shared__ float sm_virialAndEnergy[reductionBufferSize];

        if (validComponentIndex)
        {
            const int warpIndex                                     = threadLocalId / warp_size;
            sm_virialAndEnergy[warpIndex * stride + componentIndex] = virxx;
        }
        __syncthreads();

        /* Reduce to the single warp size */
        const int targetIndex = threadLocalId;
#pragma unroll
        for (int reductionStride = reductionBufferSize >> 1; reductionStride >= warp_size;
             reductionStride >>= 1)
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
        /* NOTE: This reduction assumes there are at least 4 warps (asserted).
         *       To use fewer warps, add to the conditional:
         *       && threadLocalId < activeWarps * stride
         */
        assert(activeWarps * stride >= warp_size);
        if (threadLocalId < warp_size)
        {
            float output = sm_virialAndEnergy[threadLocalId];
#pragma unroll
            for (int delta = stride; delta < warp_size; delta <<= 1)
            {
                output += __shfl_down_sync(activeMask, output, delta, warp_size);
            }
            /* Final output */
            if (validComponentIndex)
            {
                assert(isfinite(output));
                atomicAdd(gm_virialAndEnergy + componentIndex, output);
            }
        }
    }
}

//! Kernel instantiations
template __global__ void pme_solve_kernel<GridOrdering::YZX, true>(const PmeGpuCudaKernelParams);
template __global__ void pme_solve_kernel<GridOrdering::YZX, false>(const PmeGpuCudaKernelParams);
template __global__ void pme_solve_kernel<GridOrdering::XYZ, true>(const PmeGpuCudaKernelParams);
template __global__ void pme_solve_kernel<GridOrdering::XYZ, false>(const PmeGpuCudaKernelParams);
