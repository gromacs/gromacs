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
 *  \brief Implements PME force gathering in CUDA.
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
 * An inline CUDA function: unroll the dynamic index accesses to the constant grid sizes to avoid local memory operations.
 */
__device__ __forceinline__ float read_grid_size(const float* realGridSizeFP, const int dimIndex)
{
    switch (dimIndex)
    {
        case XX: return realGridSizeFP[XX];
        case YY: return realGridSizeFP[YY];
        case ZZ: return realGridSizeFP[ZZ];
    }
    assert(false);
    return 0.0f;
}

/*! \brief Reduce the partial force contributions.
 *
 * \tparam[in] order              The PME order (must be 4).
 * \tparam[in] atomDataSize       The number of partial force contributions for each atom (currently
 * order^2 == 16) \tparam[in] blockSize          The CUDA block size \param[out] sm_forces Shared
 * memory array with the output forces (number of elements is number of atoms per block) \param[in]
 * atomIndexLocal     Local atom index \param[in]  splineIndex        Spline index \param[in]
 * lineIndex          Line index (same as threadLocalId) \param[in]  realGridSizeFP     Local grid
 * size constant \param[in]  fx                 Input force partial component X \param[in]  fy Input
 * force partial component Y \param[in]  fz                 Input force partial component Z
 */
template<const int order, const int atomDataSize, const int blockSize>
__device__ __forceinline__ void reduce_atom_forces(float3* __restrict__ sm_forces,
                                                   const int    atomIndexLocal,
                                                   const int    splineIndex,
                                                   const int    lineIndex,
                                                   const float* realGridSizeFP,
                                                   float&       fx,
                                                   float&       fy,
                                                   float&       fz)
{
    if (!(order & (order - 1))) // Only for orders of power of 2
    {
        const unsigned int activeMask = c_fullWarpMask;

        // A tricky shuffle reduction inspired by reduce_force_j_warp_shfl
        // TODO: find out if this is the best in terms of transactions count
        static_assert(order == 4, "Only order of 4 is implemented");
        static_assert(atomDataSize <= warp_size,
                      "TODO: rework for atomDataSize > warp_size (order 8 or larger)");
        const int width = atomDataSize;

        fx += __shfl_down_sync(activeMask, fx, 1, width);
        fy += __shfl_up_sync(activeMask, fy, 1, width);
        fz += __shfl_down_sync(activeMask, fz, 1, width);

        if (splineIndex & 1)
        {
            fx = fy;
        }

        fx += __shfl_down_sync(activeMask, fx, 2, width);
        fz += __shfl_up_sync(activeMask, fz, 2, width);

        if (splineIndex & 2)
        {
            fx = fz;
        }

        // By now fx contains intermediate quad sums of all 3 components:
        // splineIndex    0            1            2 and 3      4            5            6 and 7 8...
        // sum of...      fx0 to fx3   fy0 to fy3   fz0 to fz3   fx4 to fx7   fy4 to fy7   fz4 to fz7 etc.

        // We have to just further reduce those groups of 4
        for (int delta = 4; delta < atomDataSize; delta <<= 1)
        {
            fx += __shfl_down_sync(activeMask, fx, delta, width);
        }

        const int dimIndex = splineIndex;
        if (dimIndex < DIM)
        {
            const float n = read_grid_size(realGridSizeFP, dimIndex);
            *((float*)(&sm_forces[atomIndexLocal]) + dimIndex) = fx * n;
        }
    }
    else
    {
        // We use blockSize shared memory elements to read fx, or fy, or fz, and then reduce them to
        // fit into smemPerDim elements which are stored separately (first 2 dimensions only)
        const int         smemPerDim   = warp_size;
        const int         smemReserved = (DIM)*smemPerDim;
        __shared__ float  sm_forceReduction[smemReserved + blockSize];
        __shared__ float* sm_forceTemp[DIM];

        const int numWarps = blockSize / smemPerDim;
        const int minStride =
                max(1, atomDataSize / numWarps); // order 4: 128 threads => 4, 256 threads => 2, etc

#pragma unroll
        for (int dimIndex = 0; dimIndex < DIM; dimIndex++)
        {
            int elementIndex = smemReserved + lineIndex;
            // Store input force contributions
            sm_forceReduction[elementIndex] = (dimIndex == XX) ? fx : (dimIndex == YY) ? fy : fz;
            // sync here because two warps write data that the first one consumes below
            __syncthreads();
            // Reduce to fit into smemPerDim (warp size)
#pragma unroll
            for (int redStride = atomDataSize / 2; redStride > minStride; redStride >>= 1)
            {
                if (splineIndex < redStride)
                {
                    sm_forceReduction[elementIndex] += sm_forceReduction[elementIndex + redStride];
                }
            }
            __syncthreads();
            // Last iteration - packing everything to be nearby, storing convenience pointer
            sm_forceTemp[dimIndex] = sm_forceReduction + dimIndex * smemPerDim;
            int redStride          = minStride;
            if (splineIndex < redStride)
            {
                const int packedIndex = atomIndexLocal * redStride + splineIndex;
                sm_forceTemp[dimIndex][packedIndex] =
                        sm_forceReduction[elementIndex] + sm_forceReduction[elementIndex + redStride];
            }
            __syncthreads();
        }

        assert((blockSize / warp_size) >= DIM);
        // assert (atomsPerBlock <= warp_size);

        const int warpIndex = lineIndex / warp_size;
        const int dimIndex  = warpIndex;

        // First 3 warps can now process 1 dimension each
        if (dimIndex < DIM)
        {
            int sourceIndex = lineIndex % warp_size;
#pragma unroll
            for (int redStride = minStride / 2; redStride > 1; redStride >>= 1)
            {
                if (!(splineIndex & redStride))
                {
                    sm_forceTemp[dimIndex][sourceIndex] += sm_forceTemp[dimIndex][sourceIndex + redStride];
                }
            }

            __syncwarp();

            const float n         = read_grid_size(realGridSizeFP, dimIndex);
            const int   atomIndex = sourceIndex / minStride;

            if (sourceIndex == minStride * atomIndex)
            {
                *((float*)(&sm_forces[atomIndex]) + dimIndex) =
                        (sm_forceTemp[dimIndex][sourceIndex] + sm_forceTemp[dimIndex][sourceIndex + 1]) * n;
            }
        }
    }
}

/*! \brief
 * A CUDA kernel which gathers the atom forces from the grid.
 * The grid is assumed to be wrapped in dimension Z.
 *
 * \tparam[in] order                The PME order (must be 4 currently).
 * \tparam[in] overwriteForces      True: the forces are written to the output buffer;
 *                                  False: the forces are added non-atomically to the output buffer (e.g. to the bonded forces).
 * \tparam[in] wrapX                Tells if the grid is wrapped in the X dimension.
 * \tparam[in] wrapY                Tells if the grid is wrapped in the Y dimension.
 * \tparam[in] readGlobal           Tells if we should read spline values from global memory
 * \tparam[in] useOrderThreads      Tells if we should use order threads per atom (order*order used if false)
 * \param[in]  kernelParams         All the PME GPU data.
 */
template<const int order, const bool overwriteForces, const bool wrapX, const bool wrapY, const bool readGlobal, const bool useOrderThreads>
__launch_bounds__(c_gatherMaxThreadsPerBlock, c_gatherMinBlocksPerMP) __global__
        void pme_gather_kernel(const PmeGpuCudaKernelParams kernelParams)
{
    /* Global memory pointers */
    const float* __restrict__ gm_coefficients = kernelParams.atoms.d_coefficients;
    const float* __restrict__ gm_grid         = kernelParams.grid.d_realGrid;
    float* __restrict__ gm_forces             = kernelParams.atoms.d_forces;

    /* Global memory pointers for readGlobal */
    const float* __restrict__ gm_theta         = kernelParams.atoms.d_theta;
    const float* __restrict__ gm_dtheta        = kernelParams.atoms.d_dtheta;
    const int* __restrict__ gm_gridlineIndices = kernelParams.atoms.d_gridlineIndices;

    float3 atomX;
    float  atomCharge;

    /* Some sizes */
    const int atomsPerBlock =
            useOrderThreads ? (c_gatherMaxThreadsPerBlock / c_pmeSpreadGatherThreadsPerAtom4ThPerAtom)
                            : (c_gatherMaxThreadsPerBlock / c_pmeSpreadGatherThreadsPerAtom);
    const int blockIndex = blockIdx.y * gridDim.x + blockIdx.x;

    /* Number of data components and threads for a single atom */
    const int atomDataSize = useOrderThreads ? c_pmeSpreadGatherThreadsPerAtom4ThPerAtom
                                             : c_pmeSpreadGatherThreadsPerAtom;
    const int atomsPerWarp = useOrderThreads ? c_pmeSpreadGatherAtomsPerWarp4ThPerAtom
                                             : c_pmeSpreadGatherAtomsPerWarp;

    const int blockSize = atomsPerBlock * atomDataSize;
    assert(blockSize == blockDim.x * blockDim.y * blockDim.z);

    /* These are the atom indices - for the shared and global memory */
    const int atomIndexLocal  = threadIdx.z;
    const int atomIndexOffset = blockIndex * atomsPerBlock;
    const int atomIndexGlobal = atomIndexOffset + atomIndexLocal;

    /* Early return for fully empty blocks at the end
     * (should only happen for billions of input atoms)
     */
    if (atomIndexOffset >= kernelParams.atoms.nAtoms)
    {
        return;
    }
    // 4 warps per block, 8 atoms per warp *3 *4
    const int        splineParamsSize    = atomsPerBlock * DIM * order;
    const int        gridlineIndicesSize = atomsPerBlock * DIM;
    __shared__ int   sm_gridlineIndices[gridlineIndicesSize];
    __shared__ float sm_theta[splineParamsSize];
    __shared__ float sm_dtheta[splineParamsSize];

    /* Spline Z coordinates */
    const int ithz = threadIdx.x;

    /* These are the spline contribution indices in shared memory */
    const int splineIndex = threadIdx.y * blockDim.x + threadIdx.x;
    const int lineIndex   = (threadIdx.z * (blockDim.x * blockDim.y))
                          + splineIndex; /* And to all the block's particles */

    const int threadLocalId =
            (threadIdx.z * (blockDim.x * blockDim.y)) + blockDim.x * threadIdx.y + threadIdx.x;
    const int threadLocalIdMax = blockDim.x * blockDim.y * blockDim.z;

    if (readGlobal)
    {
        /* Read splines */
        const int localGridlineIndicesIndex = threadLocalId;
        const int globalGridlineIndicesIndex = blockIndex * gridlineIndicesSize + localGridlineIndicesIndex;
        const int globalCheckIndices         = pme_gpu_check_atom_data_index(
                globalGridlineIndicesIndex, kernelParams.atoms.nAtoms * DIM);
        if ((localGridlineIndicesIndex < gridlineIndicesSize) & globalCheckIndices)
        {
            sm_gridlineIndices[localGridlineIndicesIndex] = gm_gridlineIndices[globalGridlineIndicesIndex];
            assert(sm_gridlineIndices[localGridlineIndicesIndex] >= 0);
        }
        /* The loop needed for order threads per atom to make sure we load all data values, as each thread must load multiple values
           with order*order threads per atom, it is only required for each thread to load one data value */

        const int iMin = 0;
        const int iMax = useOrderThreads ? 3 : 1;

        for (int i = iMin; i < iMax; i++)
        {
            int localSplineParamsIndex =
                    threadLocalId
                    + i * threadLocalIdMax; /* i will always be zero for order*order threads per atom */
            int globalSplineParamsIndex = blockIndex * splineParamsSize + localSplineParamsIndex;
            int globalCheckSplineParams = pme_gpu_check_atom_data_index(
                    globalSplineParamsIndex, kernelParams.atoms.nAtoms * DIM * order);
            if ((localSplineParamsIndex < splineParamsSize) && globalCheckSplineParams)
            {
                sm_theta[localSplineParamsIndex]  = gm_theta[globalSplineParamsIndex];
                sm_dtheta[localSplineParamsIndex] = gm_dtheta[globalSplineParamsIndex];
                assert(isfinite(sm_theta[localSplineParamsIndex]));
                assert(isfinite(sm_dtheta[localSplineParamsIndex]));
            }
        }
        __syncthreads();
    }
    else
    {
        /* Recaclulate  Splines  */
        if (c_useAtomDataPrefetch)
        {
            // charges
            __shared__ float sm_coefficients[atomsPerBlock];
            // Coordinates
            __shared__ float sm_coordinates[DIM * atomsPerBlock];
            /* Staging coefficients/charges */
            pme_gpu_stage_atom_data<float, atomsPerBlock, 1>(kernelParams, sm_coefficients,
                                                             kernelParams.atoms.d_coefficients);

            /* Staging coordinates */
            pme_gpu_stage_atom_data<float, atomsPerBlock, DIM>(kernelParams, sm_coordinates,
                                                               kernelParams.atoms.d_coordinates);
            __syncthreads();
            atomX.x    = sm_coordinates[atomIndexLocal * DIM + XX];
            atomX.y    = sm_coordinates[atomIndexLocal * DIM + YY];
            atomX.z    = sm_coordinates[atomIndexLocal * DIM + ZZ];
            atomCharge = sm_coefficients[atomIndexLocal];
        }
        else
        {
            atomCharge = gm_coefficients[atomIndexGlobal];
            atomX.x    = kernelParams.atoms.d_coordinates[atomIndexGlobal * DIM + XX];
            atomX.y    = kernelParams.atoms.d_coordinates[atomIndexGlobal * DIM + YY];
            atomX.z    = kernelParams.atoms.d_coordinates[atomIndexGlobal * DIM + ZZ];
        }
        calculate_splines<order, atomsPerBlock, atomsPerWarp, true, false>(
                kernelParams, atomIndexOffset, atomX, atomCharge, sm_theta, sm_dtheta, sm_gridlineIndices);
        __syncwarp();
    }
    float fx = 0.0f;
    float fy = 0.0f;
    float fz = 0.0f;

    const int globalCheck = pme_gpu_check_atom_data_index(atomIndexGlobal, kernelParams.atoms.nAtoms);
    const int chargeCheck = pme_gpu_check_atom_charge(gm_coefficients[atomIndexGlobal]);

    if (chargeCheck & globalCheck)
    {
        const int nx  = kernelParams.grid.realGridSize[XX];
        const int ny  = kernelParams.grid.realGridSize[YY];
        const int nz  = kernelParams.grid.realGridSize[ZZ];
        const int pny = kernelParams.grid.realGridSizePadded[YY];
        const int pnz = kernelParams.grid.realGridSizePadded[ZZ];

        const int atomWarpIndex = atomIndexLocal % atomsPerWarp;
        const int warpIndex     = atomIndexLocal / atomsPerWarp;

        const int splineIndexBase = getSplineParamIndexBase<order, atomsPerWarp>(warpIndex, atomWarpIndex);
        const int splineIndexZ = getSplineParamIndex<order, atomsPerWarp>(splineIndexBase, ZZ, ithz);
        const float2 tdz       = make_float2(sm_theta[splineIndexZ], sm_dtheta[splineIndexZ]);

        int       iz     = sm_gridlineIndices[atomIndexLocal * DIM + ZZ] + ithz;
        const int ixBase = sm_gridlineIndices[atomIndexLocal * DIM + XX];

        if (iz >= nz)
        {
            iz -= nz;
        }
        int constOffset, iy;

        const int ithyMin = useOrderThreads ? 0 : threadIdx.y;
        const int ithyMax = useOrderThreads ? order : threadIdx.y + 1;
        for (int ithy = ithyMin; ithy < ithyMax; ithy++)
        {
            const int splineIndexY = getSplineParamIndex<order, atomsPerWarp>(splineIndexBase, YY, ithy);
            const float2 tdy       = make_float2(sm_theta[splineIndexY], sm_dtheta[splineIndexY]);

            iy = sm_gridlineIndices[atomIndexLocal * DIM + YY] + ithy;
            if (wrapY & (iy >= ny))
            {
                iy -= ny;
            }
            constOffset = iy * pnz + iz;

#pragma unroll
            for (int ithx = 0; (ithx < order); ithx++)
            {
                int ix = ixBase + ithx;
                if (wrapX & (ix >= nx))
                {
                    ix -= nx;
                }
                const int gridIndexGlobal = ix * pny * pnz + constOffset;
                assert(gridIndexGlobal >= 0);
                const float gridValue = gm_grid[gridIndexGlobal];
                assert(isfinite(gridValue));
                const int splineIndexX =
                        getSplineParamIndex<order, atomsPerWarp>(splineIndexBase, XX, ithx);
                const float2 tdx  = make_float2(sm_theta[splineIndexX], sm_dtheta[splineIndexX]);
                const float  fxy1 = tdz.x * gridValue;
                const float  fz1  = tdz.y * gridValue;
                fx += tdx.y * tdy.x * fxy1;
                fy += tdx.x * tdy.y * fxy1;
                fz += tdx.x * tdy.x * fz1;
            }
        }
    }

    // Reduction of partial force contributions
    __shared__ float3 sm_forces[atomsPerBlock];
    reduce_atom_forces<order, atomDataSize, blockSize>(sm_forces, atomIndexLocal, splineIndex, lineIndex,
                                                       kernelParams.grid.realGridSizeFP, fx, fy, fz);
    __syncthreads();

    /* Calculating the final forces with no component branching, atomsPerBlock threads */
    const int forceIndexLocal  = threadLocalId;
    const int forceIndexGlobal = atomIndexOffset + forceIndexLocal;
    const int calcIndexCheck = pme_gpu_check_atom_data_index(forceIndexGlobal, kernelParams.atoms.nAtoms);
    if ((forceIndexLocal < atomsPerBlock) & calcIndexCheck)
    {
        const float3 atomForces     = sm_forces[forceIndexLocal];
        const float  negCoefficient = -gm_coefficients[forceIndexGlobal];
        float3       result;
        result.x = negCoefficient * kernelParams.current.recipBox[XX][XX] * atomForces.x;
        result.y = negCoefficient
                   * (kernelParams.current.recipBox[XX][YY] * atomForces.x
                      + kernelParams.current.recipBox[YY][YY] * atomForces.y);
        result.z = negCoefficient
                   * (kernelParams.current.recipBox[XX][ZZ] * atomForces.x
                      + kernelParams.current.recipBox[YY][ZZ] * atomForces.y
                      + kernelParams.current.recipBox[ZZ][ZZ] * atomForces.z);
        sm_forces[forceIndexLocal] = result;
    }

    __syncwarp();
    assert(atomsPerBlock <= warp_size);

    /* Writing or adding the final forces component-wise, single warp */
    const int blockForcesSize = atomsPerBlock * DIM;
    const int numIter         = (blockForcesSize + warp_size - 1) / warp_size;
    const int iterThreads     = blockForcesSize / numIter;
    if (threadLocalId < iterThreads)
    {
#pragma unroll
        for (int i = 0; i < numIter; i++)
        {
            int       outputIndexLocal  = i * iterThreads + threadLocalId;
            int       outputIndexGlobal = blockIndex * blockForcesSize + outputIndexLocal;
            const int globalOutputCheck =
                    pme_gpu_check_atom_data_index(outputIndexGlobal, kernelParams.atoms.nAtoms * DIM);
            if (globalOutputCheck)
            {
                const float outputForceComponent = ((float*)sm_forces)[outputIndexLocal];
                if (overwriteForces)
                {
                    gm_forces[outputIndexGlobal] = outputForceComponent;
                }
                else
                {
                    gm_forces[outputIndexGlobal] += outputForceComponent;
                }
            }
        }
    }
}

//! Kernel instantiations
template __global__ void pme_gather_kernel<4, true, true, true, true, true>(const PmeGpuCudaKernelParams);
template __global__ void pme_gather_kernel<4, true, true, true, true, false>(const PmeGpuCudaKernelParams);
template __global__ void pme_gather_kernel<4, false, true, true, true, true>(const PmeGpuCudaKernelParams);
template __global__ void pme_gather_kernel<4, false, true, true, true, false>(const PmeGpuCudaKernelParams);
template __global__ void pme_gather_kernel<4, true, true, true, false, true>(const PmeGpuCudaKernelParams);
template __global__ void pme_gather_kernel<4, true, true, true, false, false>(const PmeGpuCudaKernelParams);
template __global__ void pme_gather_kernel<4, false, true, true, false, true>(const PmeGpuCudaKernelParams);
template __global__ void pme_gather_kernel<4, false, true, true, false, false>(const PmeGpuCudaKernelParams);
