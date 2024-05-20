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
 *  \brief Implements PME force gathering in CUDA.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 */

#include "gmxpre.h"

#include <cassert>

#include "gromacs/gpu_utils/cuda_kernel_utils.cuh"
#include "gromacs/gpu_utils/typecasts_cuda_hip.h"

#include "pme.cuh"
#include "pme_gpu_calculate_splines.cuh"
#include "pme_grid.h"
#if GMX_NVSHMEM
#    include <nvshmemx.h>
#endif

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
    return 0.0F;
}

/*! \brief Reduce the partial force contributions.
 *
 * \tparam     order              The PME order (must be 4).
 * \tparam     atomDataSize       The number of partial force contributions for each atom (currently
 *                                order^2 == 16)
 * \tparam     blockSize          The CUDA block size
 *
 * \param[out] sm_forces          Shared memory array with the output forces (number of elements
 *                                is number of atoms per block)
 * \param[in]  atomIndexLocal     Local atom index
 * \param[in]  splineIndex        Spline index
 * \param[in]  lineIndex          Line index (same as threadLocalId)
 * \param[in]  realGridSizeFP     Local grid size constant
 * \param[in]  fx                 Input force partial component X
 * \param[in]  fy                 Input force partial component Y
 * \param[in]  fz                 Input force partial component Z
 */
template<int order, int atomDataSize, int blockSize>
__device__ __forceinline__ void reduce_atom_forces(float3* __restrict__ sm_forces,
                                                   const int    atomIndexLocal,
                                                   const int    splineIndex,
                                                   const int    lineIndex,
                                                   const float* realGridSizeFP,
                                                   float& fx, // NOLINT(google-runtime-references)
                                                   float& fy, // NOLINT(google-runtime-references)
                                                   float& fz) // NOLINT(google-runtime-references)
{
    if (gmx::isPowerOfTwo(order)) // Only for orders of power of 2
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
            float* __restrict__ sm_forcesAtomIndexOffset =
                    reinterpret_cast<float*>(&sm_forces[atomIndexLocal]);
            sm_forcesAtomIndexOffset[dimIndex] = fx * n;
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
                float* __restrict__ sm_forcesAtomIndexOffset =
                        reinterpret_cast<float*>(&sm_forces[atomIndex]);
                sm_forcesAtomIndexOffset[dimIndex] =
                        (sm_forceTemp[dimIndex][sourceIndex] + sm_forceTemp[dimIndex][sourceIndex + 1]) * n;
            }
        }
    }
}

/*! \brief Calculate the sum of the force partial components (in X, Y and Z)
 *
 * \tparam     order              The PME order (must be 4).
 * \tparam     atomsPerWarp       The number of atoms per GPU warp.
 * \tparam     wrapX              Tells if the grid is wrapped in the X dimension.
 * \tparam     wrapY              Tells if the grid is wrapped in the Y dimension.
 * \param[out] fx                 The force partial component in the X dimension.
 * \param[out] fy                 The force partial component in the Y dimension.
 * \param[out] fz                 The force partial component in the Z dimension.
 * \param[in] ithyMin             The thread minimum index in the Y dimension.
 * \param[in] ithyMax             The thread maximum index in the Y dimension.
 * \param[in] ixBase              The grid line index base value in the X dimension.
 * \param[in] iz                  The grid line index in the Z dimension.
 * \param[in] nx                  The grid real size in the X dimension.
 * \param[in] ny                  The grid real size in the Y dimension.
 * \param[in] pny                 The padded grid real size in the Y dimension.
 * \param[in] pnz                 The padded grid real size in the Z dimension.
 * \param[in] atomIndexLocal      The atom index for this thread.
 * \param[in] splineIndexBase     The base value of the spline parameter index.
 * \param[in] splineIndexZ        The theta and dtheta in the Z dimension.
 * \param[in] sm_gridlineIndices  Shared memory array of grid line indices.
 * \param[in] sm_theta            Shared memory array of atom theta values.
 * \param[in] sm_dtheta           Shared memory array of atom dtheta values.
 * \param[in] gm_grid             Global memory array of the grid to use.
 */
template<int order, int atomsPerWarp, bool wrapX, bool wrapY>
__device__ __forceinline__ void sumForceComponents(float* __restrict__ fx,
                                                   float* __restrict__ fy,
                                                   float* __restrict__ fz,
                                                   const int ithyMin,
                                                   const int ithyMax,
                                                   const int ixBase,
                                                   const int iz,
                                                   const int nx,
                                                   const int ny,
                                                   const int pny,
                                                   const int pnz,
                                                   const int atomIndexLocal,
                                                   const int splineIndexBase,
                                                   const int splineIndexZ,
                                                   const int* __restrict__ sm_gridlineIndices,
                                                   const float* __restrict__ sm_theta,
                                                   const float* __restrict__ sm_dtheta,
                                                   const float* __restrict__ gm_grid)
{
#pragma unroll
    for (int ithy = ithyMin; ithy < ithyMax; ithy++)
    {
        const int splineIndexY = getSplineParamIndex<order, atomsPerWarp>(splineIndexBase, YY, ithy);
        const float2 tdy       = make_float2(sm_theta[splineIndexY], sm_dtheta[splineIndexY]);

        int iy = sm_gridlineIndices[atomIndexLocal * DIM + YY] + ithy;
        if (wrapY & (iy >= ny))
        {
            iy -= ny;
        }
        const int constOffset = iy * pnz + iz;

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
            const int splineIndexX = getSplineParamIndex<order, atomsPerWarp>(splineIndexBase, XX, ithx);
            const float2 tdx       = make_float2(sm_theta[splineIndexX], sm_dtheta[splineIndexX]);
            const float  fxy1      = sm_theta[splineIndexZ] * gridValue;
            const float  fz1       = sm_dtheta[splineIndexZ] * gridValue;
            *fx += tdx.y * tdy.x * fxy1;
            *fy += tdx.x * tdy.y * fxy1;
            *fz += tdx.x * tdy.x * fz1;
        }
    }
}

/*! \brief Calculate the grid forces and store them in shared memory.
 *
 * \param[in,out] sm_forces       Shared memory array with the output forces.
 * \param[in] forceIndexLocal     The local (per thread) index in the sm_forces array.
 * \param[in] forceIndexGlobal    The index of the thread in the gm_coefficients array.
 * \param[in] recipBox            The reciprocal box.
 * \param[in] scale               The scale to use when calculating the forces. For gm_coefficientsB
 * (when using multiple coefficients on a single grid) the scale will be (1.0 - scale).
 * \param[in] gm_coefficients     Global memory array of the coefficients to use for an unperturbed
 * or FEP in state A if a single grid is used (\p multiCoefficientsSingleGrid == true).If two
 * separate grids are used this should be the coefficients of the grid in question.
 */
__device__ __forceinline__ void calculateAndStoreGridForces(float3* __restrict__ sm_forces,
                                                            const int   forceIndexLocal,
                                                            const int   forceIndexGlobal,
                                                            const float recipBox[DIM][DIM],
                                                            const float scale,
                                                            const float* __restrict__ gm_coefficients)
{
    const float3 atomForces     = sm_forces[forceIndexLocal];
    float        negCoefficient = -scale * gm_coefficients[forceIndexGlobal];
    float3       result;
    result.x = negCoefficient * recipBox[XX][XX] * atomForces.x;
    result.y = negCoefficient * (recipBox[XX][YY] * atomForces.x + recipBox[YY][YY] * atomForces.y);
    result.z = negCoefficient
               * (recipBox[XX][ZZ] * atomForces.x + recipBox[YY][ZZ] * atomForces.y
                  + recipBox[ZZ][ZZ] * atomForces.z);
    sm_forces[forceIndexLocal] = result;
}


template<int numIter, int iterThreads, int blockForcesSize>
__device__ void findPPRankAndPutForces(const int threadLocalId,
                                       const int ppRank,
                                       const int startAtoms,
                                       const int endAtoms,
                                       float* __restrict__ gm_forces, // NOLINT(readability-non-const-parameter)
                                       float* sm_forces) // NOLINT(readability-non-const-parameter)
{
#if GMX_NVSHMEM
    int       haveAnyWorkFromThread = 0;
    const int blockIndex            = blockIdx.y * gridDim.x + blockIdx.x;

#    pragma unroll
    for (int i = 0; i < numIter; i++)
    {
        const int outputIndexLocal  = i * iterThreads + threadLocalId;
        const int outputIndexGlobal = blockIndex * blockForcesSize + outputIndexLocal;
        haveAnyWorkFromThread += (((outputIndexGlobal < endAtoms)) && (outputIndexGlobal >= startAtoms));
    }

    const int haveAnyWorkInWarp = __any_sync(c_fullWarpMask, haveAnyWorkFromThread > 0);

    if (haveAnyWorkInWarp)
    {
        float* gm_forcesRemote = (float*)nvshmem_ptr(gm_forces, ppRank);
        // if gm_forcesRemote is non-null it implies the given ppRank is in the local node
        // so we can use normal ld/st over the remote ppRank memory with gm_forcesRemote
        const bool isPpRankIntraNode = (gm_forcesRemote != nullptr);
        const int  remoteOffset      = isPpRankIntraNode ? startAtoms : 0;
        float*     gm_forcesDest     = isPpRankIntraNode ? gm_forcesRemote : gm_forces;

#    pragma unroll
        for (int i = 0; i < numIter; i++)
        {
            const int outputIndexLocal  = i * iterThreads + threadLocalId;
            const int outputIndexGlobal = blockIndex * blockForcesSize + outputIndexLocal;
            if ((outputIndexGlobal < endAtoms) && (outputIndexGlobal >= startAtoms))
            {
                float outputForceComponent                      = sm_forces[outputIndexLocal];
                gm_forcesDest[outputIndexGlobal - remoteOffset] = outputForceComponent;
            }
        }

        if (!isPpRankIntraNode)
        {
            int totalAtomsWarpProcessed = haveAnyWorkFromThread;
            // reduce to find total atoms processed by this block for this PP rank.
            for (int mask = warp_size / 2; mask > 0; mask /= 2)
            {
                totalAtomsWarpProcessed += __shfl_xor_sync(c_fullWarpMask, totalAtomsWarpProcessed, mask);
            }

            __threadfence();

            // start offset for source buffer.
            const int startOffsetSrc = (startAtoms > (blockIndex * blockForcesSize))
                                               ? startAtoms
                                               : (blockIndex * blockForcesSize);
            // start offset for dest buffer.
            const int startOffsetDst = ((blockIndex * blockForcesSize) - startAtoms);

            nvshmemx_float_put_nbi_warp(gm_forces + startOffsetDst,     // dest buffer
                                        &gm_forcesDest[startOffsetSrc], // src buffer
                                        totalAtomsWarpProcessed,
                                        ppRank);
        }
    }
#else
    GMX_UNUSED_VALUE(threadLocalId);
    GMX_UNUSED_VALUE(ppRank);
    GMX_UNUSED_VALUE(startAtoms);
    GMX_UNUSED_VALUE(endAtoms);
    GMX_UNUSED_VALUE(gm_forces);
    GMX_UNUSED_VALUE(sm_forces);
#endif
}

// Add this prototype to fix clang warnings.
__global__ void nvshmemSignalKernel(PmeGpuCudaKernelParams kernelParams);

/*! \brief
 * A CUDA kernel which signals to the corresponding PP ranks of pme forces transfer completion.
 * This kernel should only be called after pme_gather kernel in the same stream.
 * \param[in]  kernelParams         All the PME GPU data.
 */
__global__ void nvshmemSignalKernel(const PmeGpuCudaKernelParams kernelParams)
{
#if GMX_NVSHMEM
    if (kernelParams.useNvshmem)
    {
        nvshmem_fence();
        for (int rankId = threadIdx.x; rankId < kernelParams.ppRanksInfoSize; rankId += blockDim.x)
        {
            auto ppRankInfo = kernelParams.ppRanksInfo[rankId];
            nvshmemx_signal_op(kernelParams.forcesReadyNvshmemFlags,
                               kernelParams.forcesReadyNvshmemFlagsCounter,
                               NVSHMEM_SIGNAL_SET,
                               ppRankInfo.rankId);
        }
    }
#else
    GMX_UNUSED_VALUE(kernelParams);
#endif
}

/*! \brief
 * A CUDA kernel which gathers the atom forces from the grid.
 * The grid is assumed to be wrapped in dimension Z.
 *
 * \tparam     order                The PME order (must be 4 currently).
 * \tparam     wrapX                Tells if the grid is wrapped in the X dimension.
 * \tparam     wrapY                Tells if the grid is wrapped in the Y dimension.
 * \tparam     numGrids             The number of grids to use in the kernel. Can be 1 or 2.
 * \tparam     readGlobal           Tells if we should read spline values from global memory
 * \tparam     threadsPerAtom       How many threads work on each atom
 *
 * \param[in]  kernelParams         All the PME GPU data.
 */
template<int order, bool wrapX, bool wrapY, int numGrids, bool readGlobal, ThreadsPerAtom threadsPerAtom>
__launch_bounds__(c_gatherMaxThreadsPerBlock, c_gatherMinBlocksPerMP) __global__
        void pme_gather_kernel(const PmeGpuCudaKernelParams kernelParams)
{
    assert(numGrids == 1 || numGrids == 2);

    /* Global memory pointers */
    const float* __restrict__ gm_coefficientsA = kernelParams.atoms.d_coefficients[0];
    const float* __restrict__ gm_coefficientsB = kernelParams.atoms.d_coefficients[1];
    const float* __restrict__ gm_gridA         = kernelParams.grid.d_realGrid[0];
    const float* __restrict__ gm_gridB         = kernelParams.grid.d_realGrid[1];
    static_assert(sizeof(*kernelParams.atoms.d_forces) == 3 * sizeof(float));
    float* __restrict__ gm_forces = reinterpret_cast<float*>(kernelParams.atoms.d_forces);

    /* Global memory pointers for readGlobal */
    const float* __restrict__ gm_theta         = kernelParams.atoms.d_theta;
    const float* __restrict__ gm_dtheta        = kernelParams.atoms.d_dtheta;
    const int* __restrict__ gm_gridlineIndices = kernelParams.atoms.d_gridlineIndices;

    float3 atomX;
    float  atomCharge;

    const int blockIndex = blockIdx.y * gridDim.x + blockIdx.x;

    /* Number of data components and threads for a single atom */
    const int threadsPerAtomValue = (threadsPerAtom == ThreadsPerAtom::Order) ? order : order * order;
    const int atomDataSize        = threadsPerAtomValue;
    const int atomsPerBlock       = c_gatherMaxThreadsPerBlock / atomDataSize;
    // Number of atoms processed by a single warp in spread and gather
    const int atomsPerWarp = warp_size / atomDataSize;

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
        if (localGridlineIndicesIndex < gridlineIndicesSize)
        {
            sm_gridlineIndices[localGridlineIndicesIndex] = gm_gridlineIndices[globalGridlineIndicesIndex];
            assert(sm_gridlineIndices[localGridlineIndicesIndex] >= 0);
        }
        /* The loop needed for order threads per atom to make sure we load all data values, as each thread must load multiple values
           with order*order threads per atom, it is only required for each thread to load one data value */

        const int iMin = 0;
        const int iMax = (threadsPerAtom == ThreadsPerAtom::Order) ? 3 : 1;

        for (int i = iMin; i < iMax; i++)
        {
            int localSplineParamsIndex =
                    threadLocalId
                    + i * threadLocalIdMax; /* i will always be zero for order*order threads per atom */
            int globalSplineParamsIndex = blockIndex * splineParamsSize + localSplineParamsIndex;
            if (localSplineParamsIndex < splineParamsSize)
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
        const float3* __restrict__ gm_coordinates = asFloat3(kernelParams.atoms.d_coordinates);
        /* Recalculate  Splines  */
        if (c_useAtomDataPrefetch)
        {
            // charges
            __shared__ float sm_coefficients[atomsPerBlock];
            // Coordinates
            __shared__ float3 sm_coordinates[atomsPerBlock];
            /* Staging coefficients/charges */
            pme_gpu_stage_atom_data<float, atomsPerBlock, 1>(sm_coefficients, gm_coefficientsA);

            /* Staging coordinates */
            pme_gpu_stage_atom_data<float3, atomsPerBlock, 1>(sm_coordinates, gm_coordinates);
            __syncthreads();
            atomX      = sm_coordinates[atomIndexLocal];
            atomCharge = sm_coefficients[atomIndexLocal];
        }
        else
        {
            atomX      = gm_coordinates[atomIndexGlobal];
            atomCharge = gm_coefficientsA[atomIndexGlobal];
        }
        calculate_splines<order, atomsPerBlock, atomsPerWarp, true, false, numGrids>(
                kernelParams, atomIndexOffset, atomX, atomCharge, sm_theta, sm_dtheta, sm_gridlineIndices);
        __syncwarp();
    }
    float fx = 0.0F;
    float fy = 0.0F;
    float fz = 0.0F;

    const int chargeCheck = pme_gpu_check_atom_charge(gm_coefficientsA[atomIndexGlobal]);

    const int nx  = kernelParams.grid.realGridSize[XX];
    const int ny  = kernelParams.grid.realGridSize[YY];
    const int nz  = kernelParams.grid.realGridSize[ZZ];
    const int pny = kernelParams.grid.realGridSizePadded[YY];
    const int pnz = kernelParams.grid.realGridSizePadded[ZZ];

    const int atomWarpIndex = atomIndexLocal % atomsPerWarp;
    const int warpIndex     = atomIndexLocal / atomsPerWarp;

    const int splineIndexBase = getSplineParamIndexBase<order, atomsPerWarp>(warpIndex, atomWarpIndex);
    const int splineIndexZ    = getSplineParamIndex<order, atomsPerWarp>(splineIndexBase, ZZ, ithz);

    int       iz     = sm_gridlineIndices[atomIndexLocal * DIM + ZZ] + ithz;
    const int ixBase = sm_gridlineIndices[atomIndexLocal * DIM + XX];

    if (iz >= nz)
    {
        iz -= nz;
    }

    const int ithyMin = (threadsPerAtom == ThreadsPerAtom::Order) ? 0 : threadIdx.y;
    const int ithyMax = (threadsPerAtom == ThreadsPerAtom::Order) ? order : threadIdx.y + 1;
    if (chargeCheck)
    {
        sumForceComponents<order, atomsPerWarp, wrapX, wrapY>(&fx,
                                                              &fy,
                                                              &fz,
                                                              ithyMin,
                                                              ithyMax,
                                                              ixBase,
                                                              iz,
                                                              nx,
                                                              ny,
                                                              pny,
                                                              pnz,
                                                              atomIndexLocal,
                                                              splineIndexBase,
                                                              splineIndexZ,
                                                              sm_gridlineIndices,
                                                              sm_theta,
                                                              sm_dtheta,
                                                              gm_gridA);
    }
    // Reduction of partial force contributions
    __shared__ float3 sm_forces[atomsPerBlock];
    reduce_atom_forces<order, atomDataSize, blockSize>(
            sm_forces, atomIndexLocal, splineIndex, lineIndex, kernelParams.grid.realGridSizeFP, fx, fy, fz);
    __syncthreads();

    /* Calculating the final forces with no component branching, atomsPerBlock threads */
    const int   forceIndexLocal  = threadLocalId;
    const int   forceIndexGlobal = atomIndexOffset + forceIndexLocal;
    const float scale            = kernelParams.current.scale;
    if (forceIndexLocal < atomsPerBlock)
    {
        calculateAndStoreGridForces(
                sm_forces, forceIndexLocal, forceIndexGlobal, kernelParams.current.recipBox, scale, gm_coefficientsA);
    }

    __syncwarp();
    assert(atomsPerBlock <= warp_size);

    /* Writing or adding the final forces component-wise, single warp */
    constexpr int blockForcesSize = atomsPerBlock * DIM;
    constexpr int numIter         = (blockForcesSize + warp_size - 1) / warp_size;
    constexpr int iterThreads     = blockForcesSize / numIter;

#if GMX_NVSHMEM
    if (!kernelParams.isVirialStep && kernelParams.useNvshmem)
    {
        if (threadLocalId < iterThreads)
        {
            for (int rankId = 0; rankId < kernelParams.ppRanksInfoSize; rankId++)
            {
                const auto ppRankInfo    = kernelParams.ppRanksInfo[rankId];
                const int  endAtomOffset = ppRankInfo.startAtomOffset + ppRankInfo.numAtoms;

                findPPRankAndPutForces<numIter, iterThreads, blockForcesSize>(
                        threadLocalId,
                        ppRankInfo.rankId,
                        ppRankInfo.startAtomOffset * DIM,
                        endAtomOffset * DIM,
                        gm_forces,
                        reinterpret_cast<float*>(sm_forces));
            }
        }
    }
    else
#endif
    {
        if (threadLocalId < iterThreads)
        {
#pragma unroll
            for (int i = 0; i < numIter; i++)
            {
                int   outputIndexLocal     = i * iterThreads + threadLocalId;
                int   outputIndexGlobal    = blockIndex * blockForcesSize + outputIndexLocal;
                float outputForceComponent = (reinterpret_cast<float*>(sm_forces)[outputIndexLocal]);
                gm_forces[outputIndexGlobal] = outputForceComponent;
            }
        }
    }

    if (numGrids == 2)
    {
        /* We must sync here since the same shared memory is used as above. */
        __syncthreads();
        fx                    = 0.0F;
        fy                    = 0.0F;
        fz                    = 0.0F;
        const int chargeCheck = pme_gpu_check_atom_charge(gm_coefficientsB[atomIndexGlobal]);
        if (chargeCheck)
        {
            sumForceComponents<order, atomsPerWarp, wrapX, wrapY>(&fx,
                                                                  &fy,
                                                                  &fz,
                                                                  ithyMin,
                                                                  ithyMax,
                                                                  ixBase,
                                                                  iz,
                                                                  nx,
                                                                  ny,
                                                                  pny,
                                                                  pnz,
                                                                  atomIndexLocal,
                                                                  splineIndexBase,
                                                                  splineIndexZ,
                                                                  sm_gridlineIndices,
                                                                  sm_theta,
                                                                  sm_dtheta,
                                                                  gm_gridB);
        }
        // Reduction of partial force contributions
        reduce_atom_forces<order, atomDataSize, blockSize>(
                sm_forces, atomIndexLocal, splineIndex, lineIndex, kernelParams.grid.realGridSizeFP, fx, fy, fz);
        __syncthreads();

        /* Calculating the final forces with no component branching, atomsPerBlock threads */
        if (forceIndexLocal < atomsPerBlock)
        {
            calculateAndStoreGridForces(sm_forces,
                                        forceIndexLocal,
                                        forceIndexGlobal,
                                        kernelParams.current.recipBox,
                                        1.0F - scale,
                                        gm_coefficientsB);
        }

        __syncwarp();

        /* Writing or adding the final forces component-wise, single warp */
        if (threadLocalId < iterThreads)
        {
#pragma unroll
            for (int i = 0; i < numIter; i++)
            {
                int   outputIndexLocal     = i * iterThreads + threadLocalId;
                int   outputIndexGlobal    = blockIndex * blockForcesSize + outputIndexLocal;
                float outputForceComponent = (reinterpret_cast<float*>(sm_forces)[outputIndexLocal]);
                gm_forces[outputIndexGlobal] += outputForceComponent;
            }
        }
    }
}

//! Kernel instantiations
// clang-format off
template __global__ void pme_gather_kernel<4, true, true, 1, true, ThreadsPerAtom::Order>        (const PmeGpuCudaKernelParams);
template __global__ void pme_gather_kernel<4, true, true, 1, true, ThreadsPerAtom::OrderSquared> (const PmeGpuCudaKernelParams);
template __global__ void pme_gather_kernel<4, true, true, 1, false, ThreadsPerAtom::Order>       (const PmeGpuCudaKernelParams);
template __global__ void pme_gather_kernel<4, true, true, 1, false, ThreadsPerAtom::OrderSquared>(const PmeGpuCudaKernelParams);
template __global__ void pme_gather_kernel<4, true, true, 2, true, ThreadsPerAtom::Order>        (const PmeGpuCudaKernelParams);
template __global__ void pme_gather_kernel<4, true, true, 2, true, ThreadsPerAtom::OrderSquared> (const PmeGpuCudaKernelParams);
template __global__ void pme_gather_kernel<4, true, true, 2, false, ThreadsPerAtom::Order>       (const PmeGpuCudaKernelParams);
template __global__ void pme_gather_kernel<4, true, true, 2, false, ThreadsPerAtom::OrderSquared>(const PmeGpuCudaKernelParams);
// clang-format on
