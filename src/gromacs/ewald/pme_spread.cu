/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013-2016,2017,2018,2019, by the GROMACS development team, led by
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
 *  \brief Implements PME GPU spline calculation and charge spreading in CUDA.
 *  TODO: consider always pre-sorting particles (as in DD case).
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 */

#include "gmxpre.h"

#include <cassert>

#include "gromacs/gpu_utils/cuda_kernel_utils.cuh"

#include "pme.cuh"
#include "pme_calculate_splines.cuh"
#include "pme_gpu_utils.h"
#include "pme_grid.h"

/*! \brief
 * Charge spreading onto the grid.
 * This corresponds to the CPU function spread_coefficients_bsplines_thread().
 * Optional second stage of the spline_and_spread_kernel.
 *
 * \tparam[in] order                PME interpolation order.
 * \tparam[in] wrapX                A boolean which tells if the grid overlap in dimension X should
 * be wrapped. \tparam[in] wrapY                A boolean which tells if the grid overlap in
 * dimension Y should be wrapped. \tparam[in] useOrderThreads      A boolean which Tells if we
 * should use order threads per atom (order*order used if false) \param[in]  kernelParams Input PME
 * CUDA data in constant memory. \param[in]  atomIndexOffset      Starting atom index for the
 * execution block w.r.t. global memory. \param[in]  atomCharge           Atom charge/coefficient of
 * atom processed by thread. \param[in]  sm_gridlineIndices   Atom gridline indices in the shared
 * memory. \param[in]  sm_theta             Atom spline values in the shared memory.
 */
template<const int order, const bool wrapX, const bool wrapY, const bool useOrderThreads>
__device__ __forceinline__ void spread_charges(const PmeGpuCudaKernelParams kernelParams,
                                               int                          atomIndexOffset,
                                               const float*                 atomCharge,
                                               const int* __restrict__ sm_gridlineIndices,
                                               const float* __restrict__ sm_theta)
{
    /* Global memory pointer to the output grid */
    float* __restrict__ gm_grid = kernelParams.grid.d_realGrid;


    const int atomsPerWarp = useOrderThreads ? c_pmeSpreadGatherAtomsPerWarp4ThPerAtom
                                             : c_pmeSpreadGatherAtomsPerWarp;

    const int nx  = kernelParams.grid.realGridSize[XX];
    const int ny  = kernelParams.grid.realGridSize[YY];
    const int nz  = kernelParams.grid.realGridSize[ZZ];
    const int pny = kernelParams.grid.realGridSizePadded[YY];
    const int pnz = kernelParams.grid.realGridSizePadded[ZZ];

    const int offx = 0, offy = 0, offz = 0; // unused for now

    const int atomIndexLocal  = threadIdx.z;
    const int atomIndexGlobal = atomIndexOffset + atomIndexLocal;

    const int globalCheck = pme_gpu_check_atom_data_index(atomIndexGlobal, kernelParams.atoms.nAtoms);
    const int chargeCheck = pme_gpu_check_atom_charge(*atomCharge);
    if (chargeCheck & globalCheck)
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
        const int ithyMin = useOrderThreads ? 0 : threadIdx.y;
        const int ithyMax = useOrderThreads ? order : threadIdx.y + 1;
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
 * \tparam[in] order                PME interpolation order.
 * \tparam[in] computeSplines       A boolean which tells if the spline parameter and
 *                                  gridline indices' computation should be performed.
 * \tparam[in] spreadCharges        A boolean which tells if the charge spreading should be performed.
 * \tparam[in] wrapX                A boolean which tells if the grid overlap in dimension X should be wrapped.
 * \tparam[in] wrapY                A boolean which tells if the grid overlap in dimension Y should be wrapped.
 * \tparam[in] writeGlobal          A boolean which tells if the theta values and gridlines should be written to global memory.
 * \tparam[in] useOrderThreads         A boolean which tells if we should use order threads per atom (order*order used if false).
 * \param[in]  kernelParams         Input PME CUDA data in constant memory.
 */
template<const int order, const bool computeSplines, const bool spreadCharges, const bool wrapX, const bool wrapY, const bool writeGlobal, const bool useOrderThreads>
__launch_bounds__(c_spreadMaxThreadsPerBlock) CLANG_DISABLE_OPTIMIZATION_ATTRIBUTE __global__
        void pme_spline_and_spread_kernel(const PmeGpuCudaKernelParams kernelParams)
{
    const int atomsPerBlock =
            useOrderThreads ? c_spreadMaxThreadsPerBlock / c_pmeSpreadGatherThreadsPerAtom4ThPerAtom
                            : c_spreadMaxThreadsPerBlock / c_pmeSpreadGatherThreadsPerAtom;
    // Gridline indices, ivec
    __shared__ int sm_gridlineIndices[atomsPerBlock * DIM];
    // Spline values
    __shared__ float sm_theta[atomsPerBlock * DIM * order];
    float            dtheta;

    const int atomsPerWarp = useOrderThreads ? c_pmeSpreadGatherAtomsPerWarp4ThPerAtom
                                             : c_pmeSpreadGatherAtomsPerWarp;

    float3 atomX;
    float  atomCharge;

    const int blockIndex      = blockIdx.y * gridDim.x + blockIdx.x;
    const int atomIndexOffset = blockIndex * atomsPerBlock;

    /* Thread index w.r.t. block */
    const int threadLocalId =
            (threadIdx.z * (blockDim.x * blockDim.y)) + (threadIdx.y * blockDim.x) + threadIdx.x;
    /* Warp index w.r.t. block - could probably be obtained easier? */
    const int warpIndex = threadLocalId / warp_size;

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
    if (c_useAtomDataPrefetch)
    {
        __shared__ float sm_coefficients[atomsPerBlock];
        pme_gpu_stage_atom_data<float, atomsPerBlock, 1>(kernelParams, sm_coefficients,
                                                         kernelParams.atoms.d_coefficients);
        __syncthreads();
        atomCharge = sm_coefficients[atomIndexLocal];
    }
    else
    {
        atomCharge = kernelParams.atoms.d_coefficients[atomIndexGlobal];
    }

    if (computeSplines)
    {
        if (c_useAtomDataPrefetch)
        {
            // Coordinates
            __shared__ float sm_coordinates[DIM * atomsPerBlock];

            /* Staging coordinates */
            pme_gpu_stage_atom_data<float, atomsPerBlock, DIM>(kernelParams, sm_coordinates,
                                                               kernelParams.atoms.d_coordinates);
            __syncthreads();
            atomX.x = sm_coordinates[atomIndexLocal * DIM + XX];
            atomX.y = sm_coordinates[atomIndexLocal * DIM + YY];
            atomX.z = sm_coordinates[atomIndexLocal * DIM + ZZ];
        }
        else
        {
            atomX.x = kernelParams.atoms.d_coordinates[atomIndexGlobal * DIM + XX];
            atomX.y = kernelParams.atoms.d_coordinates[atomIndexGlobal * DIM + YY];
            atomX.z = kernelParams.atoms.d_coordinates[atomIndexGlobal * DIM + ZZ];
        }
        calculate_splines<order, atomsPerBlock, atomsPerWarp, false, writeGlobal>(
                kernelParams, atomIndexOffset, atomX, atomCharge, sm_theta, &dtheta, sm_gridlineIndices);
        __syncwarp();
    }
    else
    {
        /* Staging the data for spread
         * (the data is assumed to be in GPU global memory with proper layout already,
         * as in after running the spline kernel)
         */
        /* Spline data - only thetas (dthetas will only be needed in gather) */
        pme_gpu_stage_atom_data<float, atomsPerBlock, DIM * order>(kernelParams, sm_theta,
                                                                   kernelParams.atoms.d_theta);
        /* Gridline indices */
        pme_gpu_stage_atom_data<int, atomsPerBlock, DIM>(kernelParams, sm_gridlineIndices,
                                                         kernelParams.atoms.d_gridlineIndices);

        __syncthreads();
    }

    /* Spreading */
    if (spreadCharges)
    {
        spread_charges<order, wrapX, wrapY, useOrderThreads>(
                kernelParams, atomIndexOffset, &atomCharge, sm_gridlineIndices, sm_theta);
    }
}

//! Kernel instantiations
template __global__ void pme_spline_and_spread_kernel<4, true, true, true, true, true, true>(const PmeGpuCudaKernelParams);
template __global__ void
pme_spline_and_spread_kernel<4, true, false, true, true, true, true>(const PmeGpuCudaKernelParams);
template __global__ void
pme_spline_and_spread_kernel<4, false, true, true, true, true, true>(const PmeGpuCudaKernelParams);

template __global__ void
pme_spline_and_spread_kernel<4, true, true, true, true, false, true>(const PmeGpuCudaKernelParams);

template __global__ void
pme_spline_and_spread_kernel<4, true, true, true, true, true, false>(const PmeGpuCudaKernelParams);
template __global__ void
pme_spline_and_spread_kernel<4, true, false, true, true, true, false>(const PmeGpuCudaKernelParams);
template __global__ void
pme_spline_and_spread_kernel<4, false, true, true, true, true, false>(const PmeGpuCudaKernelParams);

template __global__ void
pme_spline_and_spread_kernel<4, true, true, true, true, false, false>(const PmeGpuCudaKernelParams);
