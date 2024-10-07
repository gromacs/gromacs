/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
 *  \brief Implements PME GPU spline calculation and charge spreading in HIP.
 *  TODO: consider always pre-sorting particles (as in DD case).
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 */

#include "gmxpre.h"

#include <cassert>

#include "gromacs/gpu_utils/hip_kernel_utils.h"
#include "gromacs/gpu_utils/typecasts_cuda_hip.h"

#include "pme_gpu_calculate_splines_hip.h"
#include "pme_gpu_constants.h"
#include "pme_gpu_types.h"

template<int parallelExecutionWidth>
static constexpr int sc_spreadHipMaxWarpsPerBlock = (parallelExecutionWidth == 64) ? 8 : 4;

template<int parallelExecutionWidth>
static constexpr int                                          sc_spreadMaxThreadsPerBlock =
        sc_spreadHipMaxWarpsPerBlock<parallelExecutionWidth>* parallelExecutionWidth;

/*! \brief
 * Charge spreading onto the grid.
 * This corresponds to the CPU function spread_coefficients_bsplines_thread().
 * Optional second stage of the spline_and_spread_kernel.
 *
 * \tparam     order                PME interpolation order.
 * \tparam     wrapX                Whether the grid overlap in dimension X should be wrapped.
 * \tparam     wrapY                Whether the grid overlap in dimension Y should be wrapped.
 * \tparam     gridIndex            The index of the grid to use in the kernel.
 * \tparam     threadsPerAtom       How many threads work on each atom
 *
 * \param[in]  kernelParams         Input PME HIP data in constant memory.
 * \param[in]  atomCharge           Atom charge/coefficient of atom processed by thread.
 * \param[in]  sm_gridlineIndices   Atom gridline indices in the shared memory.
 * \param[in]  sm_theta             Atom spline values in the shared memory.
 */
template<int order, bool wrapX, bool wrapY, int gridIndex, ThreadsPerAtom threadsPerAtom, int parallelExecutionWidth>
__device__ __forceinline__ void spreadChargeKernel(const PmeGpuKernelParamsBase kernelParams,
                                                   const float*                 atomCharge,
                                                   const int* __restrict__ sm_gridlineIndices,
                                                   const float* __restrict__ sm_theta)
{
    /* Global memory pointer to the output grid */
    float* __restrict__ gm_grid = kernelParams.grid.d_realGrid[gridIndex];

    // Number of atoms processed by a single warp in spread and gather
    constexpr int threadsPerAtomValue = (threadsPerAtom == ThreadsPerAtom::Order) ? order : order * order;
    constexpr int atomsPerWarp        = parallelExecutionWidth / threadsPerAtomValue;

    const int nx  = kernelParams.grid.realGridSize[XX];
    const int ny  = kernelParams.grid.realGridSize[YY];
    const int nz  = kernelParams.grid.realGridSize[ZZ];
    const int pny = kernelParams.grid.realGridSizePadded[YY];
    const int pnz = kernelParams.grid.realGridSizePadded[ZZ];

    const int offx = 0, offy = 0, offz = 0; // unused for now

    const int atomIndexLocal = threadIdx.z;

    const bool chargeCheck = pme_gpu_check_atom_charge(*atomCharge);
    if (chargeCheck)
    {
        // Spline Z coordinates
        const int ithz = threadIdx.x;

        const int ixBase = sm_gridlineIndices[atomIndexLocal * DIM + XX] - offx;
        const int iyBase = sm_gridlineIndices[atomIndexLocal * DIM + YY] - offy;
        int       iz     = sm_gridlineIndices[atomIndexLocal * DIM + ZZ] - offz + ithz;
        if (iz >= nz)
        {
            iz -= nz;
        }
        /* Atom index w.r.t. warp - alternating 0 1 0 1 .. */
        const int atomWarpIndex = atomIndexLocal % atomsPerWarp;
        /* Warp index w.r.t. block - could probably be obtained easier? */
        const int warpIndex = atomIndexLocal / atomsPerWarp;

        const int splineIndexBase = getSplineParamIndexBase<order, atomsPerWarp>(warpIndex, atomWarpIndex);
        const int splineIndexZ = getSplineParamIndex<order, atomsPerWarp>(splineIndexBase, ZZ, ithz);
        const float thetaZ     = sm_theta[splineIndexZ];

        /* loop not used if order*order threads per atom */
        const int ithyMin = (threadsPerAtom == ThreadsPerAtom::Order) ? 0 : threadIdx.y;
        const int ithyMax = (threadsPerAtom == ThreadsPerAtom::Order) ? order : threadIdx.y + 1;
        for (int ithy = ithyMin; ithy < ithyMax; ithy++)
        {
            int iy = iyBase + ithy;
            if (wrapY & (iy >= ny))
            {
                iy -= ny;
            }

            const int splineIndexY = getSplineParamIndex<order, atomsPerWarp>(splineIndexBase, YY, ithy);
            float       thetaY = sm_theta[splineIndexY];
            const float Val    = thetaZ * thetaY * (*atomCharge);
            assert(isfinite(Val));
            const int offset = iy * pnz + iz;

#pragma unroll
            for (int ithx = 0; (ithx < order); ithx++)
            {
                int ix = ixBase + ithx;
                if (wrapX & (ix >= nx))
                {
                    ix -= nx;
                }
                const int gridIndexGlobal = ix * pny * pnz + offset;
                const int splineIndexX =
                        getSplineParamIndex<order, atomsPerWarp>(splineIndexBase, XX, ithx);
                const float thetaX = sm_theta[splineIndexX];
                assert(isfinite(thetaX));
                assert(isfinite(gm_grid[gridIndexGlobal]));
                atomicAdd(gm_grid + gridIndexGlobal, thetaX * Val);
            }
        }
    }
}

/*! \brief
 * A spline computation and charge spreading kernel function.
 *
 * Two tuning parameters can be used for additional performance. For small systems and for debugging
 * writeGlobal should be used removing the need to recalculate the theta values in the gather kernel.
 * Similarly for useOrderThreads large systems order threads per atom gives higher performance than order*order threads
 *
 * \tparam     order                PME interpolation order.
 * \tparam     computeSplines       A boolean which tells if the spline parameter and
 *                                  gridline indices' computation should be performed.
 * \tparam     spreadCharges        A boolean which tells if the charge spreading should be performed.
 * \tparam     wrapX                A boolean which tells if the grid overlap in dimension X should be wrapped.
 * \tparam     wrapY                A boolean which tells if the grid overlap in dimension Y should be wrapped.
 * \tparam     numGrids             The number of grids to use in the kernel. Can be 1 or 2.
 * \tparam     writeGlobal          A boolean which tells if the theta values and gridlines should be written to global memory.
 * \tparam     threadsPerAtom       How many threads work on each atom
 * param[in]  kernelParams         Input PME HIP data in constant memory.
 */
template<int order, bool computeSplines, bool spreadCharges, bool wrapX, bool wrapY, int numGrids, bool writeGlobal, ThreadsPerAtom threadsPerAtom, int parallelExecutionWidth>
LAUNCH_BOUNDS_EXACT_SINGLE(sc_spreadMaxThreadsPerBlock<parallelExecutionWidth>)
__global__ void pmeSplineAndSpreadKernel(const PmeGpuKernelParamsBase kernelParams)
{
    constexpr int threadsPerAtomValue = (threadsPerAtom == ThreadsPerAtom::Order) ? order : order * order;
    constexpr int atomsPerBlock = sc_spreadMaxThreadsPerBlock<parallelExecutionWidth> / threadsPerAtomValue;
    // Number of atoms processed by a single warp in spread and gather
    constexpr int atomsPerWarp = parallelExecutionWidth / threadsPerAtomValue;
    // Gridline indices, ivec
    __shared__ int sm_gridlineIndices[atomsPerBlock * DIM];
    // Charges
    __shared__ float sm_coefficients[atomsPerBlock];
    // Spline values
    __shared__ float sm_theta[atomsPerBlock * DIM * order];
    float            dtheta;

    float3 atomX;
    float  atomCharge;

    const int blockIndex      = blockIdx.y * gridDim.x + blockIdx.x;
    const int atomIndexOffset = blockIndex * atomsPerBlock + kernelParams.pipelineAtomStart;

    /* Thread index w.r.t. block */
    const int threadLocalId =
            (threadIdx.z * (blockDim.x * blockDim.y)) + (threadIdx.y * blockDim.x) + threadIdx.x;
    /* Warp index w.r.t. block - could probably be obtained easier? */
    const int warpIndex = threadLocalId / parallelExecutionWidth;

    /* Atom index w.r.t. warp */
    const int atomWarpIndex = threadIdx.z % atomsPerWarp;
    /* Atom index w.r.t. block/shared memory */
    const int atomIndexLocal = warpIndex * atomsPerWarp + atomWarpIndex;
    /* Atom index w.r.t. global memory */
    const int atomIndexGlobal = atomIndexOffset + atomIndexLocal;

    /* Early return for fully empty blocks at the end
     * (should only happen for billions of input atoms)
     */
    if (atomIndexOffset >= kernelParams.atoms.nAtoms)
    {
        return;
    }
    /* Charges, required for both spline and spread */
    if constexpr (c_useAtomDataPrefetch)
    {
        pme_gpu_stage_atom_data<float, atomsPerBlock, 1>(
                sm_coefficients, &kernelParams.atoms.d_coefficients[0][kernelParams.pipelineAtomStart]);
        __syncthreads();
        atomCharge = sm_coefficients[atomIndexLocal];
    }
    else
    {
        atomCharge = kernelParams.atoms.d_coefficients[0][atomIndexGlobal];
    }

    if (computeSplines)
    {
        const float3* __restrict__ gm_coordinates = asFloat3(kernelParams.atoms.d_coordinates);
        if constexpr (c_useAtomDataPrefetch)
        {
            // Coordinates
            __shared__ float3 sm_coordinates[atomsPerBlock];

            /* Staging coordinates */
            pme_gpu_stage_atom_data<float3, atomsPerBlock, 1>(
                    sm_coordinates, &gm_coordinates[kernelParams.pipelineAtomStart]);
            __syncthreads();
            atomX = sm_coordinates[atomIndexLocal];
        }
        else
        {
            atomX = gm_coordinates[atomIndexGlobal];
        }
        calculate_splines<order, atomsPerBlock, atomsPerWarp, false, writeGlobal, numGrids, parallelExecutionWidth>(
                kernelParams, atomIndexOffset, atomX, atomCharge, sm_theta, &dtheta, sm_gridlineIndices);
        __builtin_amdgcn_wave_barrier();
    }
    else
    {
        /* Staging the data for spread
         * (the data is assumed to be in GPU global memory with proper layout already,
         * as in after running the spline kernel)
         */
        /* Spline data - only thetas (dthetas will only be needed in gather) */
        pme_gpu_stage_atom_data<float, atomsPerBlock, DIM * order>(sm_theta, kernelParams.atoms.d_theta);
        /* Gridline indices */
        pme_gpu_stage_atom_data<int, atomsPerBlock, DIM>(sm_gridlineIndices,
                                                         kernelParams.atoms.d_gridlineIndices);

        __syncthreads();
    }

    /* Spreading */
    if constexpr (spreadCharges)
    {
        const bool runSpreadChargeKernel =
                atomIndexGlobal < kernelParams.atoms.nAtoms
                && (!kernelParams.usePipeline || (atomIndexGlobal < kernelParams.pipelineAtomEnd));
        if (runSpreadChargeKernel)
        {
            spreadChargeKernel<order, wrapX, wrapY, 0, threadsPerAtom, parallelExecutionWidth>(
                    kernelParams, &atomCharge, sm_gridlineIndices, sm_theta);
        }
    }
    if (numGrids == 2)
    {
        __syncthreads();
        if constexpr (c_useAtomDataPrefetch)
        {
            pme_gpu_stage_atom_data<float, atomsPerBlock, 1>(sm_coefficients,
                                                             kernelParams.atoms.d_coefficients[1]);
            __syncthreads();
            atomCharge = sm_coefficients[atomIndexLocal];
        }
        else
        {
            atomCharge = kernelParams.atoms.d_coefficients[1][atomIndexGlobal];
        }
        if constexpr (spreadCharges)
        {
            const bool runSpreadChargeKernel =
                    atomIndexGlobal < kernelParams.atoms.nAtoms
                    && (!kernelParams.usePipeline || (atomIndexGlobal < kernelParams.pipelineAtomEnd));
            if (runSpreadChargeKernel)
            {
                spreadChargeKernel<order, wrapX, wrapY, 1, threadsPerAtom, parallelExecutionWidth>(
                        kernelParams, &atomCharge, sm_gridlineIndices, sm_theta);
            }
        }
    }
}

//! Kernel instantiations
/* Disable the "explicit template instantiation 'PmeSplineAndSpreadKernel<...>' will emit a vtable in every
 * translation unit [-Wweak-template-vtables]" warning.
 * It is only explicitly instantiated in this translation unit, so we should be safe.
 */
CLANG_DIAGNOSTIC_IGNORE("-Wweak-template-vtables")

#define INSTANTIATE_3(order, computeSplines, spreadCharges, numGrids, writeGlobal, threadsPerAtom, parallelExecutionWidth)                     \
    template __global__ void                                                                                                                   \
    pmeSplineAndSpreadKernel<order, computeSplines, spreadCharges, true, true, numGrids, writeGlobal, threadsPerAtom, parallelExecutionWidth>( \
            PmeGpuKernelParamsBase kernelParams);

#define INSTANTIATE_2(order, numGrids, threadsPerAtom, parallelExecutionWidth)                 \
    INSTANTIATE_3(order, true, true, numGrids, true, threadsPerAtom, parallelExecutionWidth);  \
    INSTANTIATE_3(order, true, false, numGrids, true, threadsPerAtom, parallelExecutionWidth); \
    INSTANTIATE_3(order, false, true, numGrids, true, threadsPerAtom, parallelExecutionWidth); \
    INSTANTIATE_3(order, true, true, numGrids, false, threadsPerAtom, parallelExecutionWidth);

#define INSTANTIATE(order, parallelExecutionWidth)                                 \
    INSTANTIATE_2(order, 1, ThreadsPerAtom::Order, parallelExecutionWidth);        \
    INSTANTIATE_2(order, 1, ThreadsPerAtom::OrderSquared, parallelExecutionWidth); \
    INSTANTIATE_2(order, 2, ThreadsPerAtom::Order, parallelExecutionWidth);        \
    INSTANTIATE_2(order, 2, ThreadsPerAtom::OrderSquared, parallelExecutionWidth);

INSTANTIATE(4, 32);
INSTANTIATE(4, 64);

CLANG_DIAGNOSTIC_RESET
