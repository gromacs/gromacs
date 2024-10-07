/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
 *  \brief Implements PME GPU Fourier grid solving in HIP.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 */

#include "gmxpre.h"

#include <cassert>

#include <hip/hip_math_constants.h>

#include "gromacs/gpu_utils/device_utils_hip_sycl.h"
#include "gromacs/gpu_utils/gputraits.h"
#include "gromacs/gpu_utils/hip_kernel_utils.h"
#include "gromacs/gpu_utils/vectype_ops_hip.h"

#include "pme_gpu_constants.h"
#include "pme_gpu_internal.h"
#include "pme_gpu_types.h"

template<int parallelExecutionWidth>
static constexpr int sc_solveMaxThreadsPerBlock = c_solveMaxWarpsPerBlock* parallelExecutionWidth;

/*! \brief
 * PME complex grid solver kernel function.
 *
 * \tparam     gridOrdering             Specifies the dimension ordering of the complex grid.
 * \tparam     computeEnergyAndVirial   Tells if the reciprocal energy and virial should be computed.
 * \tparam     gridIndex                The index of the grid to use in the kernel.
 * \param[in]  kernelParams             Input PME HIP data in constant memory.
 */
template<GridOrdering gridOrdering, bool computeEnergyAndVirial, const int gridIndex, int parallelExecutionWidth>
LAUNCH_BOUNDS_EXACT_SINGLE(sc_solveMaxThreadsPerBlock<parallelExecutionWidth>)
__global__ void pmeSolveKernel(const struct PmeGpuKernelParamsBase kernelParams)
{
    /* This kernel supports 2 different grid dimension orderings: YZX and XYZ */
    constexpr int majorDim  = gridOrdering == GridOrdering::YZX ? YY : XX;
    constexpr int middleDim = gridOrdering == GridOrdering::YZX ? ZZ : YY;
    constexpr int minorDim  = gridOrdering == GridOrdering::YZX ? XX : ZZ;

    /* Global memory pointers */
    const float* __restrict__ gm_splineValueMajor = kernelParams.grid.d_splineModuli[gridIndex]
                                                    + kernelParams.grid.splineValuesOffset[majorDim];
    const float* __restrict__ gm_splineValueMiddle = kernelParams.grid.d_splineModuli[gridIndex]
                                                     + kernelParams.grid.splineValuesOffset[middleDim];
    const float* __restrict__ gm_splineValueMinor = kernelParams.grid.d_splineModuli[gridIndex]
                                                    + kernelParams.grid.splineValuesOffset[minorDim];
    float* __restrict__ gm_virialAndEnergy = kernelParams.constants.d_virialAndEnergy[gridIndex];
    float2* __restrict__ gm_grid =
            reinterpret_cast<float2*>(kernelParams.grid.d_fftComplexGrid[gridIndex]);

    /* Various grid sizes and indices */
    const int localOffsetMinor  = kernelParams.grid.kOffsets[minorDim];
    const int localOffsetMiddle = kernelParams.grid.kOffsets[middleDim];
    const int localOffsetMajor  = kernelParams.grid.kOffsets[majorDim];
    const int localSizeMinor    = kernelParams.grid.localComplexGridSizePadded[minorDim];
    const int localSizeMiddle   = kernelParams.grid.localComplexGridSizePadded[middleDim];
    const int localCountMiddle  = kernelParams.grid.localComplexGridSize[middleDim];
    const int localCountMinor   = kernelParams.grid.localComplexGridSize[minorDim];
    const int nMajor            = kernelParams.grid.realGridSize[majorDim];
    const int nMiddle           = kernelParams.grid.realGridSize[middleDim];
    const int nMinor            = kernelParams.grid.realGridSize[minorDim];
    const int maxkMajor         = (nMajor + 1) / 2;  // X or Y
    const int maxkMiddle        = (nMiddle + 1) / 2; // Y OR Z => only check for !YZX
    const int maxkMinor         = (nMinor + 1) / 2;  // Z or X => only check for YZX

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
    const int activeWarps       = blockDim.x / parallelExecutionWidth;
    const int indexMinor        = blockIdx.x * blockDim.x + gridLineCellIndex;
    const int indexMiddle       = blockIdx.y * gridLinesPerBlock + gridLineIndex;
    const int indexMajor        = blockIdx.z;

    /* Optional outputs */
    float energy = 0.0F;
    float virxx  = 0.0F;
    float virxy  = 0.0F;
    float virxz  = 0.0F;
    float viryy  = 0.0F;
    float viryz  = 0.0F;
    float virzz  = 0.0F;

    assert(indexMajor < kernelParams.grid.localComplexGridSize[majorDim]);
    if ((indexMiddle < localCountMiddle) & (indexMinor < localCountMinor)
        & (gridLineIndex < gridLinesPerBlock))
    {
        /* The offset should be equal to the global thread index for coalesced access */
        const int gridThreadIndex =
                (indexMajor * localSizeMiddle + indexMiddle) * localSizeMinor + indexMinor;
        float2* __restrict__ gm_gridCell = gm_grid + gridThreadIndex;

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

        const AmdPackedFloat3 mm = (gridOrdering == GridOrdering::YZX)
                                           ? AmdPackedFloat3(mMinor, mMajor, mMiddle)
                                           : AmdPackedFloat3(mMajor, mMiddle, mMinor);

        const bool isLastComponent = gridOrdering == GridOrdering::YZX
                                             ? ((kMiddle == 0) | (kMiddle == maxkMiddle))
                                             : ((kMinor == 0) | (kMinor == maxkMinor));
        /* 0.5 correction factor for the first and last components of a Z dimension */
        const float corner_fac = isLastComponent ? 0.5F : 1.0F;

        // TODO: use textures for gm_splineValue
        float vMajor  = LDG(gm_splineValueMajor + kMajor);
        float vMiddle = LDG(gm_splineValueMiddle + kMiddle);
        float vMinor  = LDG(gm_splineValueMinor + kMinor);

        if (notZeroPoint)
        {
            AmdPackedFloat3 mhk, recip1, recip2, recip3;
            recip1.x_         = kernelParams.current.recipBox[XX][XX];
            recip1.y_         = kernelParams.current.recipBox[XX][YY];
            recip1.z_         = kernelParams.current.recipBox[XX][ZZ];
            recip2.x_         = 0.0f;
            recip2.y_         = kernelParams.current.recipBox[YY][YY];
            recip2.z_         = kernelParams.current.recipBox[YY][ZZ];
            recip3.x_         = 0.0f;
            recip3.y_         = 0.0f;
            recip3.z_         = kernelParams.current.recipBox[ZZ][ZZ];
            AmdPackedFloat3 a = mm.x() * recip1;
            AmdPackedFloat3 b = mm.y() * recip2;
            AmdPackedFloat3 c = mm.z() * recip3;
            mhk               = a + b + c;

            const float m2k = mhk.norm2();
            assert(m2k != 0.0F);

            float denom = m2k * float(HIP_PI) * kernelParams.current.boxVolume * vMajor * vMiddle * vMinor;
            assert(isfinite(denom));
            assert(denom != 0.0F);

            const float tmp1   = __expf(-kernelParams.grid.ewaldFactor * m2k);
            const float etermk = kernelParams.constants.elFactor * tmp1 / denom;

            float2       gridValue    = *gm_gridCell;
            const float2 oldGridValue = gridValue;
            gridValue                 = gridValue * etermk;
            *gm_gridCell              = gridValue;

            if constexpr (computeEnergyAndVirial)
            {
                const float tmp1k =
                        2.0F * (gridValue.x * oldGridValue.x + gridValue.y * oldGridValue.y);

                float vfactor = (kernelParams.grid.ewaldFactor + 1.0F / m2k) * 2.0F;
                float ets2    = corner_fac * tmp1k;
                energy        = ets2;

                float ets2vf = ets2 * vfactor;

                virxx = ets2vf * mhk.x() * mhk.x() - ets2;
                virxy = ets2vf * mhk.x() * mhk.y();
                virxz = ets2vf * mhk.x() * mhk.z();
                viryy = ets2vf * mhk.y() * mhk.y() - ets2;
                viryz = ets2vf * mhk.y() * mhk.z();
                virzz = ets2vf * mhk.z() * mhk.z() - ets2;
            }
        }
    }

    /* Optional energy/virial reduction */
    if constexpr (computeEnergyAndVirial)
    {

        virxx += amdDppUpdateShfl<float, /* row_shl:1 */ 0x101>(virxx);
        viryy += amdDppUpdateShfl<float, /* row_shr:1 */ 0x111>(viryy);
        virzz += amdDppUpdateShfl<float, /* row_shl:1 */ 0x101>(virzz);
        virxy += amdDppUpdateShfl<float, /* row_shr:1 */ 0x111>(virxy);
        virxz += amdDppUpdateShfl<float, /* row_shl:1 */ 0x101>(virxz);
        viryz += amdDppUpdateShfl<float, /* row_shr:1 */ 0x111>(viryz);
        energy += amdDppUpdateShfl<float, /* row_shl:1 */ 0x101>(energy);
        if (threadLocalId & 1)
        {
            virxx = viryy; // virxx now holds virxx and viryy pair sums
            virzz = virxy; // virzz now holds virzz and virxy pair sums
            virxz = viryz; // virxz now holds virxz and viryz pair sums
        }

        virxx += amdDppUpdateShfl<float, /* row_shl:2 */ 0x102>(virxx);
        virzz += amdDppUpdateShfl<float, /* row_shr:2 */ 0x112>(virzz);
        virxz += amdDppUpdateShfl<float, /* row_shl:2 */ 0x102>(virxz);
        energy += amdDppUpdateShfl<float, /* row_shr:2 */ 0x112>(energy);
        if (threadLocalId & 2)
        {
            virxx = virzz;  // virxx now holds quad sums of virxx, virxy, virzz and virxy
            virxz = energy; // virxz now holds quad sums of virxz, viryz, energy and unused paddings
        }

        virxx += amdDppUpdateShfl<float, /* row_shl:4 */ 0x104>(virxx);
        virxz += amdDppUpdateShfl<float, /* row_shr:4 */ 0x114>(virxz);
        if (threadLocalId & 4)
        {
            virxx = virxz; // virxx now holds all 7 components' octet sums + unused paddings
        }

        /* We only need to reduce virxx now */
#pragma unroll
        for (int delta = 8; delta < parallelExecutionWidth; delta <<= 1)
        {
            virxx += __shfl_down(virxx, delta, parallelExecutionWidth);
        }

        /* Now first 7 threads of each warp have the full output contributions in virxx */

        const int  componentIndex      = threadLocalId & (parallelExecutionWidth - 1);
        const bool validComponentIndex = (componentIndex < c_virialAndEnergyCount);
        /* Reduce 7 outputs per warp in the shared memory */
        constexpr int stride =
                8; // this is c_virialAndEnergyCount==7 rounded up to power of 2 for convenience, hence the assert
        static_assert(c_virialAndEnergyCount == 7);
        const int reductionBufferSize = sc_solveMaxThreadsPerBlock<parallelExecutionWidth> * stride;
        __shared__ float sm_virialAndEnergy[reductionBufferSize];

        if (validComponentIndex)
        {
            const int warpIndex = threadLocalId / parallelExecutionWidth;
            sm_virialAndEnergy[warpIndex * stride + componentIndex] = virxx;
        }
        __syncthreads();

        /* Reduce to the single warp size */
        const int targetIndex = threadLocalId;
#pragma unroll
        for (int reductionStride = reductionBufferSize >> 1; reductionStride >= parallelExecutionWidth;
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
        assert(activeWarps * stride >= parallelExecutionWidth);
        if (threadLocalId < parallelExecutionWidth)
        {
            float output = sm_virialAndEnergy[threadLocalId];
#pragma unroll
            for (int delta = stride; delta < parallelExecutionWidth; delta <<= 1)
            {
                output += __shfl_down(output, delta, parallelExecutionWidth);
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

//! Kernel class instantiations
/* Disable the "explicit template instantiation 'PmeSplineAndSpreadKernel<...>' will emit a vtable in every
 * translation unit [-Wweak-template-vtables]" warning.
 * It is only explicitly instantiated in this translation unit, so we should be safe.
 */
CLANG_DIAGNOSTIC_IGNORE("-Wweak-template-vtables")

#define INSTANTIATE(parallelExecutionWidth)                                                       \
    template __global__ void pmeSolveKernel<GridOrdering::XYZ, false, 0, parallelExecutionWidth>( \
            PmeGpuKernelParamsBase kernelParams);                                                 \
    template __global__ void pmeSolveKernel<GridOrdering::XYZ, true, 0, parallelExecutionWidth>(  \
            PmeGpuKernelParamsBase kernelParams);                                                 \
    template __global__ void pmeSolveKernel<GridOrdering::YZX, false, 0, parallelExecutionWidth>( \
            PmeGpuKernelParamsBase kernelParams);                                                 \
    template __global__ void pmeSolveKernel<GridOrdering::YZX, true, 0, parallelExecutionWidth>(  \
            PmeGpuKernelParamsBase kernelParams);                                                 \
    template __global__ void pmeSolveKernel<GridOrdering::XYZ, false, 1, parallelExecutionWidth>( \
            PmeGpuKernelParamsBase kernelParams);                                                 \
    template __global__ void pmeSolveKernel<GridOrdering::XYZ, true, 1, parallelExecutionWidth>(  \
            PmeGpuKernelParamsBase kernelParams);                                                 \
    template __global__ void pmeSolveKernel<GridOrdering::YZX, false, 1, parallelExecutionWidth>( \
            PmeGpuKernelParamsBase kernelParams);                                                 \
    template __global__ void pmeSolveKernel<GridOrdering::YZX, true, 1, parallelExecutionWidth>(  \
            PmeGpuKernelParamsBase kernelParams);

INSTANTIATE(32);
INSTANTIATE(64);

CLANG_DIAGNOSTIC_RESET
