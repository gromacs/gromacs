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
 *  \brief Implements PME GPU Fourier grid solving in SYCL.
 *
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 */

#include "gmxpre.h"

#include "pme_solve_sycl.h"

#include <cassert>

#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/gpu_utils/sycl_kernel_utils.h"
#include "gromacs/math/units.h"
#include "gromacs/utility/basedefinitions.h"

#include "pme_gpu_constants.h"

using mode = sycl::access_mode;

/*! \brief
 * PME complex grid solver kernel function.
 *
 * \tparam     gridOrdering             Specifies the dimension ordering of the complex grid.
 * \tparam     computeEnergyAndVirial   Tells if the reciprocal energy and virial should be
 *                                        computed.
 * \tparam     subGroupSize             Describes the width of a SYCL subgroup
 */
template<GridOrdering gridOrdering, bool computeEnergyAndVirial, int subGroupSize>
auto makeSolveKernel(sycl::handler& cgh,
                     const float* __restrict__ gm_splineModuli,
                     SolveKernelParams solveKernelParams,
                     float* __restrict__ gm_virialAndEnergy,
                     float* __restrict__ gm_fourierGrid_)
{
    /* Reduce 7 outputs per warp in the shared memory */
    const int stride =
            8; // this is c_virialAndEnergyCount==7 rounded up to power of 2 for convenience, hence the assert
    static_assert(c_virialAndEnergyCount == 7);
    const int reductionBufferSize = c_solveMaxWarpsPerBlock * stride;

    // Help compiler eliminate local buffer when it is unused.
    auto sm_virialAndEnergy = [&]() {
        if constexpr (computeEnergyAndVirial)
        {
            return sycl::local_accessor<float, 1>(sycl::range<1>(reductionBufferSize), cgh);
        }
        else
        {
            return nullptr;
        }
    }();

    /* Each thread works on one cell of the Fourier space complex 3D grid (gm_grid).
     * Each block handles up to c_solveMaxWarpsPerBlock * subGroupSize cells -
     * depending on the grid contiguous dimension size,
     * that can range from a part of a single gridline to several complete gridlines.
     */
    return [=](sycl::nd_item<3> itemIdx) [[intel::reqd_sub_group_size(subGroupSize)]]
    {
        if constexpr (skipKernelCompilation<subGroupSize>())
        {
            return;
        }
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

            default: SYCL_ASSERT(false);
        }

        /* Global memory pointers */
        const float* __restrict__ gm_splineValueMajor =
                gm_splineModuli + solveKernelParams.splineValuesOffset[majorDim];
        const float* __restrict__ gm_splineValueMiddle =
                gm_splineModuli + solveKernelParams.splineValuesOffset[middleDim];
        const float* __restrict__ gm_splineValueMinor =
                gm_splineModuli + solveKernelParams.splineValuesOffset[minorDim];
        // The Fourier grid is allocated as float values, even though
        // it logically contains complex values. (It also can be
        // the same memory as the real grid for in-place transforms.)
        // The buffer underlying the accessor may have a size that is
        // larger than the active grid, because it is allocated with
        // reallocateDeviceBuffer. The size of that larger-than-needed
        // grid can be an odd number of floats, even though actual
        // grid code only accesses up to an even number of floats. If
        // we would use the reinterpet method of the accessor to
        // convert from float to float2, runtime boundary checks can
        // fail because of this mismatch. So, we extract the
        // underlying global_ptr and use that to construct
        // sycl::float2 values when needed.
        sycl::global_ptr<float> gm_fourierGrid = gm_fourierGrid_;

        /* Various grid sizes and indices */
        const int localOffsetMinor  = solveKernelParams.kOffsets[minorDim];
        const int localOffsetMiddle = solveKernelParams.kOffsets[middleDim];
        const int localOffsetMajor  = solveKernelParams.kOffsets[majorDim];
        const int localSizeMinor    = solveKernelParams.complexGridSizePadded[minorDim];
        const int localSizeMiddle   = solveKernelParams.complexGridSizePadded[middleDim];
        const int localCountMiddle  = solveKernelParams.complexGridSize[middleDim];
        const int localCountMinor   = solveKernelParams.complexGridSize[minorDim];
        const int nMajor            = solveKernelParams.realGridSize[majorDim];
        const int nMiddle           = solveKernelParams.realGridSize[middleDim];
        const int nMinor            = solveKernelParams.realGridSize[minorDim];
        const int maxkMajor         = (nMajor + 1) / 2;  // X or Y
        const int maxkMiddle        = (nMiddle + 1) / 2; // Y OR Z => only check for !YZX
        const int maxkMinor         = (nMinor + 1) / 2;  // Z or X => only check for YZX

        const int threadLocalId     = itemIdx.get_local_linear_id();
        const int gridLineSize      = localCountMinor;
        const int gridLineIndex     = threadLocalId / gridLineSize;
        const int gridLineCellIndex = threadLocalId - gridLineSize * gridLineIndex;
        const int gridLinesPerBlock =
                sycl::max(itemIdx.get_local_range(2) / size_t(gridLineSize), size_t(1));
        const int activeWarps = (itemIdx.get_local_range(2) / subGroupSize);
        const int indexMinor = itemIdx.get_group(2) * itemIdx.get_local_range(2) + gridLineCellIndex;
        const int indexMiddle = itemIdx.get_group(1) * gridLinesPerBlock + gridLineIndex;
        const int indexMajor  = itemIdx.get_group(0);

        /* Optional outputs */
        float energy = 0.0F;
        float virxx  = 0.0F;
        float virxy  = 0.0F;
        float virxz  = 0.0F;
        float viryy  = 0.0F;
        float viryz  = 0.0F;
        float virzz  = 0.0F;

        SYCL_ASSERT(indexMajor < solveKernelParams.complexGridSize[majorDim]);
        if ((indexMiddle < localCountMiddle) & (indexMinor < localCountMinor)
            & (gridLineIndex < gridLinesPerBlock))
        {
            /* The offset should be equal to the global thread index for coalesced access */
            const int gridThreadIndex =
                    (indexMajor * localSizeMiddle + indexMiddle) * localSizeMinor + indexMinor;

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

                default: SYCL_ASSERT(false);
            }

            /* 0.5 correction factor for the first and last components of a Z dimension */
            float corner_fac = 1.0F;
            switch (gridOrdering)
            {
                case GridOrdering::YZX:
                    if ((kMiddle == 0) | (kMiddle == maxkMiddle))
                    {
                        corner_fac = 0.5F;
                    }
                    break;

                case GridOrdering::XYZ:
                    if ((kMinor == 0) | (kMinor == maxkMinor))
                    {
                        corner_fac = 0.5F;
                    }
                    break;

                default: SYCL_ASSERT(false);
            }

            if (notZeroPoint)
            {
                const float mhxk = mX * solveKernelParams.recipBox[XX][XX];
                const float mhyk = mX * solveKernelParams.recipBox[XX][YY]
                                   + mY * solveKernelParams.recipBox[YY][YY];
                const float mhzk = mX * solveKernelParams.recipBox[XX][ZZ]
                                   + mY * solveKernelParams.recipBox[YY][ZZ]
                                   + mZ * solveKernelParams.recipBox[ZZ][ZZ];

                const float m2k = mhxk * mhxk + mhyk * mhyk + mhzk * mhzk;
                SYCL_ASSERT(m2k != 0.0F);
                float denom = m2k * float(M_PI) * solveKernelParams.boxVolume * gm_splineValueMajor[kMajor]
                              * gm_splineValueMiddle[kMiddle] * gm_splineValueMinor[kMinor];
                SYCL_ASSERT(sycl::isfinite(denom));
                SYCL_ASSERT(denom != 0.0F);

                const float tmp1   = sycl::exp(-solveKernelParams.ewaldFactor * m2k);
                const float etermk = solveKernelParams.elFactor * tmp1 / denom;

                sycl::float2 gridValue;
                gridValue.load(gridThreadIndex, sycl::global_ptr<const float>(gm_fourierGrid));
                const sycl::float2 oldGridValue = gridValue;
                gridValue *= etermk;
                gridValue.store(gridThreadIndex, gm_fourierGrid);

                if (computeEnergyAndVirial)
                {
                    const float tmp1k = 2.0F * sycl::dot(gridValue, oldGridValue);

                    float vfactor = (solveKernelParams.ewaldFactor + 1.0F / m2k) * 2.0F;
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
        if constexpr (computeEnergyAndVirial)
        {
            /* A tricky shuffle reduction inspired by reduce_force_j_warp_shfl.
             * The idea is to reduce 7 energy/virial components into a single variable (aligned by
             * 8). We will reduce everything into virxx.
             */

            /* We can only reduce warp-wise */
            const int width = subGroupSize;
            static_assert(subGroupSize >= 8);

            sycl::sub_group sg = itemIdx.get_sub_group();

            /* Making pair sums */
            virxx += sycl::shift_group_left(sg, virxx, 1);
            viryy += sycl::shift_group_right(sg, viryy, 1);
            virzz += sycl::shift_group_left(sg, virzz, 1);
            virxy += sycl::shift_group_right(sg, virxy, 1);
            virxz += sycl::shift_group_left(sg, virxz, 1);
            viryz += sycl::shift_group_right(sg, viryz, 1);
            energy += sycl::shift_group_left(sg, energy, 1);
            if (threadLocalId & 1)
            {
                virxx = viryy; // virxx now holds virxx and viryy pair sums
                virzz = virxy; // virzz now holds virzz and virxy pair sums
                virxz = viryz; // virxz now holds virxz and viryz pair sums
            }

            /* Making quad sums */
            virxx += sycl::shift_group_left(sg, virxx, 2);
            virzz += sycl::shift_group_right(sg, virzz, 2);
            virxz += sycl::shift_group_left(sg, virxz, 2);
            energy += sycl::shift_group_right(sg, energy, 2);
            if (threadLocalId & 2)
            {
                virxx = virzz; // virxx now holds quad sums of virxx, virxy, virzz and virxy
                virxz = energy; // virxz now holds quad sums of virxz, viryz, energy and unused paddings
            }

            /* Making octet sums */
            virxx += sycl::shift_group_left(sg, virxx, 4);
            virxz += sycl::shift_group_right(sg, virxz, 4);
            if (threadLocalId & 4)
            {
                virxx = virxz; // virxx now holds all 7 components' octet sums + unused paddings
            }

            /* We only need to reduce virxx now */
#pragma unroll
            for (int delta = 8; delta < width; delta <<= 1)
            {
                virxx += sycl::shift_group_left(sg, virxx, delta);
            }
            /* Now first 7 threads of each warp have the full output contributions in virxx */

            const int  componentIndex      = threadLocalId & (subGroupSize - 1);
            const bool validComponentIndex = (componentIndex < c_virialAndEnergyCount);

            if (validComponentIndex)
            {
                const int warpIndex = threadLocalId / subGroupSize;
                sm_virialAndEnergy[warpIndex * stride + componentIndex] = virxx;
            }
            itemIdx.barrier(sycl::access::fence_space::local_space);

            /* Reduce to the single warp size */
            const int targetIndex = threadLocalId;
#pragma unroll
            for (int reductionStride = reductionBufferSize >> 1; reductionStride >= subGroupSize;
                 reductionStride >>= 1)
            {
                const int sourceIndex = targetIndex + reductionStride;
                if ((targetIndex < reductionStride) & (sourceIndex < activeWarps * stride))
                {
                    sm_virialAndEnergy[targetIndex] += sm_virialAndEnergy[sourceIndex];
                }
                itemIdx.barrier(sycl::access::fence_space::local_space);
            }

            /* Now use shuffle again */
            /* NOTE: This reduction assumes there are at least 4 warps (asserted).
             *       To use fewer warps, add to the conditional:
             *       && threadLocalId < activeWarps * stride
             */
            SYCL_ASSERT(activeWarps * stride >= subGroupSize);
            if (threadLocalId < subGroupSize)
            {
                float output = sm_virialAndEnergy[threadLocalId];
#pragma unroll
                for (int delta = stride; delta < subGroupSize; delta <<= 1)
                {
                    output += sycl::shift_group_left(sg, output, delta);
                }
                /* Final output */
                if (validComponentIndex)
                {
                    SYCL_ASSERT(sycl::isfinite(output));
                    atomicFetchAdd(gm_virialAndEnergy[componentIndex], output);
                }
            }
        }
    };
}

template<GridOrdering gridOrdering, bool computeEnergyAndVirial, int gridIndex, int subGroupSize>
PmeSolveKernel<gridOrdering, computeEnergyAndVirial, gridIndex, subGroupSize>::PmeSolveKernel()
{
    reset();
}

template<GridOrdering gridOrdering, bool computeEnergyAndVirial, int gridIndex, int subGroupSize>
void PmeSolveKernel<gridOrdering, computeEnergyAndVirial, gridIndex, subGroupSize>::setArg(size_t argIndex,
                                                                                           void* arg)
{
    if (argIndex == 0)
    {
        auto* params = reinterpret_cast<PmeGpuKernelParams*>(arg);

        constParams_                             = &params->constants;
        gridParams_                              = &params->grid;
        solveKernelParams_.ewaldFactor           = params->grid.ewaldFactor;
        solveKernelParams_.realGridSize          = params->grid.realGridSize;
        solveKernelParams_.kOffsets              = params->grid.kOffsets;
        solveKernelParams_.complexGridSize       = params->grid.localComplexGridSize;
        solveKernelParams_.complexGridSizePadded = params->grid.localComplexGridSizePadded;
        solveKernelParams_.splineValuesOffset    = params->grid.splineValuesOffset;
        solveKernelParams_.recipBox[XX]          = params->current.recipBox[XX];
        solveKernelParams_.recipBox[YY]          = params->current.recipBox[YY];
        solveKernelParams_.recipBox[ZZ]          = params->current.recipBox[ZZ];
        solveKernelParams_.boxVolume             = params->current.boxVolume;
        solveKernelParams_.elFactor              = params->constants.elFactor;
    }
    else
    {
        GMX_RELEASE_ASSERT(argIndex == 0, "Trying to pass too many args to the solve kernel");
    }
}

template<GridOrdering gridOrdering, bool computeEnergyAndVirial, int gridIndex, int subGroupSize>
void PmeSolveKernel<gridOrdering, computeEnergyAndVirial, gridIndex, subGroupSize>::launch(
        const KernelLaunchConfig& config,
        const DeviceStream&       deviceStream)
{
    GMX_RELEASE_ASSERT(gridParams_, "Can not launch the kernel before setting its args");
    GMX_RELEASE_ASSERT(constParams_, "Can not launch the kernel before setting its args");

    using KernelNameType = PmeSolveKernel<gridOrdering, computeEnergyAndVirial, gridIndex, subGroupSize>;

    // SYCL has different multidimensional layout than OpenCL/CUDA.
    const sycl::range<3> localSize{ config.blockSize[2], config.blockSize[1], config.blockSize[0] };
    const sycl::range<3> groupRange{ config.gridSize[2], config.gridSize[1], config.gridSize[0] };
    const sycl::nd_range<3> range{ groupRange * localSize, localSize };

    sycl::queue q = deviceStream.stream();

    q.submit(GMX_SYCL_DISCARD_EVENT[&](sycl::handler & cgh) {
        auto kernel = makeSolveKernel<gridOrdering, computeEnergyAndVirial, subGroupSize>(
                cgh,
                gridParams_->d_splineModuli[gridIndex].get_pointer(),
                solveKernelParams_,
                constParams_->d_virialAndEnergy[gridIndex].get_pointer(),
                gridParams_->d_fftComplexGrid[gridIndex].get_pointer());
        cgh.parallel_for<KernelNameType>(range, kernel);
    });

    // Delete set args, so we don't forget to set them before the next launch.
    reset();
}

template<GridOrdering gridOrdering, bool computeEnergyAndVirial, int gridIndex, int subGroupSize>
void PmeSolveKernel<gridOrdering, computeEnergyAndVirial, gridIndex, subGroupSize>::reset()
{
    gridParams_  = nullptr;
    constParams_ = nullptr;
}

//! Kernel class instantiations
/* Disable the "explicit template instantiation 'PmeSplineAndSpreadKernel<...>' will emit a vtable in every
 * translation unit [-Wweak-template-vtables]" warning.
 * It is only explicitly instantiated in this translation unit, so we should be safe.
 */
CLANG_DIAGNOSTIC_IGNORE("-Wweak-template-vtables")

#define INSTANTIATE(subGroupSize)                                             \
    template class PmeSolveKernel<GridOrdering::XYZ, false, 0, subGroupSize>; \
    template class PmeSolveKernel<GridOrdering::XYZ, true, 0, subGroupSize>;  \
    template class PmeSolveKernel<GridOrdering::YZX, false, 0, subGroupSize>; \
    template class PmeSolveKernel<GridOrdering::YZX, true, 0, subGroupSize>;  \
    template class PmeSolveKernel<GridOrdering::XYZ, false, 1, subGroupSize>; \
    template class PmeSolveKernel<GridOrdering::XYZ, true, 1, subGroupSize>;  \
    template class PmeSolveKernel<GridOrdering::YZX, false, 1, subGroupSize>; \
    template class PmeSolveKernel<GridOrdering::YZX, true, 1, subGroupSize>;

#if GMX_SYCL_DPCPP
INSTANTIATE(16);
#endif
INSTANTIATE(32);
INSTANTIATE(64);

CLANG_DIAGNOSTIC_RESET
