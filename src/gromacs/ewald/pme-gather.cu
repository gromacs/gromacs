/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018, by the GROMACS development team, led by
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
#include "pme-grid.h"
#include "pme-gpu-utils.h"

/*! \brief
 * An inline CUDA function: unroll the dynamic index accesses to the constant grid sizes to avoid local memory operations.
 */
__device__ __forceinline__ float read_grid_size(const float *realGridSizeFP,
                                                const int    dimIndex)
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

template<typename T,
         const int atomsPerBlock,
         const int dataCountPerAtom>
__device__  __forceinline__
void pme_gpu_stage_atom_data(const PmeGpuCudaKernelParams       kernelParams,
                             T * __restrict__                   sm_destination,
                             const T * __restrict__             gm_source)
{
    static_assert(c_usePadding, "With padding disabled, index checking should be fixed to account for spline theta/dtheta per-warp alignment");
    const int blockIndex       = blockIdx.y * gridDim.x + blockIdx.x;
    const int threadLocalIndex = ((threadIdx.z * blockDim.y + threadIdx.y) * blockDim.x) + threadIdx.x;
    const int localIndex       = threadLocalIndex;
    const int globalIndexBase  = blockIndex * atomsPerBlock * dataCountPerAtom;
    const int globalIndex      = globalIndexBase + localIndex;
    const int globalCheck      = pme_gpu_check_atom_data_index(globalIndex, kernelParams.atoms.nAtoms * dataCountPerAtom);
    if ((localIndex < atomsPerBlock * dataCountPerAtom) && globalCheck)
    {
        assert(isfinite(float(gm_source[globalIndex])));
        sm_destination[localIndex] = gm_source[globalIndex];
    }
}

#define PME_GPU_PARALLEL_SPLINE 0

template <const int order,
          const int atomsPerBlock>
__device__ __forceinline__ void calculate_splines_gather(const PmeGpuCudaKernelParams           kernelParams,
                                                  const int                              atomIndexOffset,
                                                  const float3 * __restrict__ atomX,
                                                  const float * __restrict__ atomCharge,
//                                                  const float3 * __restrict__            sm_coordinates,
//                                                  const float * __restrict__             sm_coefficients,
                                                  float2 * __restrict__                   sm_theta,
                                                  int * __restrict__                     sm_gridlineIndices)
{
    /* Global memory pointers for output */
 //   float * __restrict__ gm_theta           = kernelParams.atoms.d_theta;
 //   float * __restrict__ gm_dtheta          = kernelParams.atoms.d_dtheta;
    int * __restrict__   gm_gridlineIndices = kernelParams.atoms.d_gridlineIndices;

    const int            atomsPerWarp = c_pmeGatherAtomsPerWarp;

    /* Fractional coordinates */
    __shared__ float sm_fractCoords[atomsPerBlock * DIM];

    /* Thread index w.r.t. block */
    const int threadLocalId = (threadIdx.z * (blockDim.x * blockDim.y))
        + (threadIdx.y * blockDim.x) + threadIdx.x;
    /* Warp index w.r.t. block - could probably be obtained easier? */
    const int warpIndex = threadLocalId / warp_size;
    /* Thread index w.r.t. warp */
    const int threadWarpIndex = threadLocalId % warp_size;
    /* Atom index w.r.t. warp - alternating 0 1 0 1 .. */
//    const int atomWarpIndex = threadWarpIndex % atomsPerWarp;
    const int atomWarpIndex = threadIdx.z %atomsPerWarp;
    /* Atom index w.r.t. block/shared memory */
    const int atomIndexLocal = warpIndex * atomsPerWarp + atomWarpIndex;

    /* Atom index w.r.t. global memory */
    const int atomIndexGlobal = atomIndexOffset + atomIndexLocal;
    /* Spline contribution index in one dimension */
//    const int orderIndex = threadWarpIndex / (atomsPerWarp * DIM);
    const int threadLocalIdXY = (threadIdx.y * blockDim.x) + threadIdx.x;
    const int orderIndex = threadLocalIdXY / DIM ;
    /* Dimension index */
//    const int dimIndex = (threadWarpIndex - orderIndex * (atomsPerWarp * DIM)) / atomsPerWarp;
    const int dimIndex = threadLocalIdXY % DIM ;

    /* Multi-purpose index of rvec/ivec atom data */
    const int sharedMemoryIndex = atomIndexLocal * DIM + dimIndex;

    /* Spline parameters need a working buffer.
     * With PME_GPU_PARALLEL_SPLINE == 0 it is just a local array of (order) values for some of the threads, which is fine;
     * With PME_GPU_PARALLEL_SPLINE == 1 (order) times more threads are involved, so the shared memory is used to avoid overhead.
     * The buffer's size, striding and indexing are adjusted accordingly.
     * The buffer is accessed with SPLINE_DATA_PTR and SPLINE_DATA macros.
     */
#if PME_GPU_PARALLEL_SPLINE
    const int        splineDataStride  = atomsPerBlock * DIM;
   const int        splineDataIndex   = sharedMemoryIndex;
    __shared__ float sm_splineData[splineDataStride * order];
    float           *splineDataPtr = sm_splineData;
#else
    const int        splineDataStride = 1;
    const int        splineDataIndex  = 0;
    float            splineData[splineDataStride * order];
    float           *splineDataPtr = splineData;
#endif

#define SPLINE_DATA_PTR(i) (splineDataPtr + ((i) * splineDataStride + splineDataIndex))
#define SPLINE_DATA(i) (*SPLINE_DATA_PTR(i))

    const int localCheck  = (dimIndex < DIM) && (orderIndex < (PME_GPU_PARALLEL_SPLINE ? order : 1));
    const int globalCheck = pme_gpu_check_atom_data_index(atomIndexGlobal, kernelParams.atoms.nAtoms);

    if (localCheck && globalCheck)
    {
        /* Indices interpolation */

        if (orderIndex == 0)
        {
            int           tableIndex, tInt;
            float         n, t;
            assert(atomIndexLocal<DIM*atomsPerBlock);
//            const float3  x = sm_coordinates[atomIndexLocal];
             const float3  x = *atomX;
            /* Accessing fields in fshOffset/nXYZ/recipbox/... with dimIndex offset
             * puts them into local memory(!) instead of accessing the constant memory directly.
             * That's the reason for the switch, to unroll explicitly.
             * The commented parts correspond to the 0 components of the recipbox.
             */
            switch (dimIndex)
            {
                case XX:
                    tableIndex  = kernelParams.grid.tablesOffsets[XX];
                    n           = kernelParams.grid.realGridSizeFP[XX];
                    t           = x.x * kernelParams.current.recipBox[dimIndex][XX] + x.y * kernelParams.current.recipBox[dimIndex][YY] + x.z * kernelParams.current.recipBox[dimIndex][ZZ];
                    break;

                case YY:
                    tableIndex  = kernelParams.grid.tablesOffsets[YY];
                    n           = kernelParams.grid.realGridSizeFP[YY];
                    t           = /*x.x * kernelParams.current.recipbox[dimIndex][XX] + */ x.y * kernelParams.current.recipBox[dimIndex][YY] + x.z * kernelParams.current.recipBox[dimIndex][ZZ];
                    break;

                case ZZ:
                    tableIndex  = kernelParams.grid.tablesOffsets[ZZ];
                    n           = kernelParams.grid.realGridSizeFP[ZZ];
                    t           = /*x.x * kernelParams.current.recipbox[dimIndex][XX] + x.y * kernelParams.current.recipbox[dimIndex][YY] + */ x.z * kernelParams.current.recipBox[dimIndex][ZZ];
                    break;
            }
            const float shift = c_pmeMaxUnitcellShift;
            /* Fractional coordinates along box vectors, adding a positive shift to ensure t is positive for triclinic boxes */
            t    = (t + shift) * n;
            tInt = (int)t;
            assert(sharedMemoryIndex<atomsPerBlock*DIM);
            sm_fractCoords[sharedMemoryIndex] = t - tInt;
            tableIndex                       += tInt;
            assert(tInt >= 0);
            assert(tInt < c_pmeNeighborUnitcellCount * n);

            // TODO have shared table for both parameters to share the fetch, as index is always same?
            // TODO compare texture/LDG performance
            sm_fractCoords[sharedMemoryIndex] +=
                fetchFromParamLookupTable(kernelParams.grid.d_fractShiftsTable,
                                          kernelParams.fractShiftsTableTexture,
                                          tableIndex);
            sm_gridlineIndices[sharedMemoryIndex] =
                fetchFromParamLookupTable(kernelParams.grid.d_gridlineIndicesTable,
                                          kernelParams.gridlineIndicesTableTexture,
                                          tableIndex);
//            gm_gridlineIndices[atomIndexOffset * DIM + sharedMemoryIndex] = sm_gridlineIndices[sharedMemoryIndex];
        }

        /* B-spline calculation */

//        const int chargeCheck = pme_gpu_check_atom_charge(sm_coefficients[atomIndexLocal]);
        const int chargeCheck = pme_gpu_check_atom_charge(*atomCharge);
        if (chargeCheck)
        {
            float       div;
            int         o = orderIndex; // This is an index that is set once for PME_GPU_PARALLEL_SPLINE == 1

            const float dr = sm_fractCoords[sharedMemoryIndex];
            assert(isfinite(dr));

            /* dr is relative offset from lower cell limit */
            *SPLINE_DATA_PTR(order - 1) = 0.0f;
            *SPLINE_DATA_PTR(1)         = dr;
            *SPLINE_DATA_PTR(0)         = 1.0f - dr;

#pragma unroll
            for (int k = 3; k < order; k++)
            {
                div                     = 1.0f / (k - 1.0f);
                *SPLINE_DATA_PTR(k - 1) = div * dr * SPLINE_DATA(k - 2);
#pragma unroll
                for (int l = 1; l < (k - 1); l++)
                {
                    *SPLINE_DATA_PTR(k - l - 1) = div * ((dr + l) * SPLINE_DATA(k - l - 2) + (k - l - dr) * SPLINE_DATA(k - l - 1));
                }
                *SPLINE_DATA_PTR(0) = div * (1.0f - dr) * SPLINE_DATA(0);
            }

            const int thetaIndexBase        = getSplineParamIndexBase<order, atomsPerWarp>(warpIndex, atomWarpIndex);
            const int thetaGlobalOffsetBase = atomIndexOffset * DIM * order;

            /* Differentiation and storing the spline derivatives (dtheta) */
#if !PME_GPU_PARALLEL_SPLINE
            // With PME_GPU_PARALLEL_SPLINE == 1, o is already set to orderIndex;
            // With PME_GPU_PARALLEL_SPLINE == 0, we loop o over range(order).
#pragma unroll
            for (o = 0; o < order; o++)
#endif
            {
                const int   thetaIndex       = getSplineParamIndex<order, atomsPerWarp>(thetaIndexBase, dimIndex, o);
//                const int   thetaGlobalIndex = thetaGlobalOffsetBase + thetaIndex;

                const float dtheta = ((o > 0) ? SPLINE_DATA(o - 1) : 0.0f) - SPLINE_DATA(o);
                assert(isfinite(dtheta));
                assert(thetaIndex<order*DIM*atomsPerBlock);
                sm_theta[thetaIndex].y = dtheta;
            }

            div  = 1.0f / (order - 1.0f);
            *SPLINE_DATA_PTR(order - 1) = div * dr * SPLINE_DATA(order - 2);
#pragma unroll
            for (int k = 1; k < (order - 1); k++)
            {
                *SPLINE_DATA_PTR(order - k - 1) = div * ((dr + k) * SPLINE_DATA(order - k - 2) + (order - k - dr) * SPLINE_DATA(order - k - 1));
            }
            *SPLINE_DATA_PTR(0) = div * (1.0f - dr) * SPLINE_DATA(0);

            /* Storing the spline values (theta) */
#if !PME_GPU_PARALLEL_SPLINE
            // See comment for the loop above
#pragma unroll
            for (o = 0; o < order; o++)
#endif
            {
                const int thetaIndex       = getSplineParamIndex<order, atomsPerWarp>(thetaIndexBase, dimIndex, o);
                const int thetaGlobalIndex = thetaGlobalOffsetBase + thetaIndex;
                assert(thetaIndex<order*DIM*atomsPerBlock);
                sm_theta[thetaIndex].x       = SPLINE_DATA(o);
                
                assert(isfinite(sm_theta[thetaIndex].x));
//                gm_theta[thetaGlobalIndex] = SPLINE_DATA(o);
            }
        }
    }
}

/*! \brief Reduce the partial force contributions.
 *
 * \tparam[in] order              The PME order (must be 4).
 * \tparam[in] atomDataSize       The number of partial force contributions for each atom (currently order^2 == 16)
 * \tparam[in] blockSize          The CUDA block size
 * \param[out] sm_forces          Shared memory array with the output forces (number of elements is number of atoms per block)
 * \param[in]  atomIndexLocal     Local atom index
 * \param[in]  splineIndex        Spline index
 * \param[in]  lineIndex          Line index (same as threadLocalId)
 * \param[in]  realGridSizeFP     Local grid size constant
 * \param[in]  fx                 Input force partial component X
 * \param[in]  fy                 Input force partial component Y
 * \param[in]  fz                 Input force partial component Z
 */
template <
    const int order,
    const int atomDataSize,
    const int blockSize
    >
__device__ __forceinline__ void reduce_atom_forces(float3 * __restrict__ sm_forces,
                                                   const int             atomIndexLocal,
                                                   const int             splineIndex,
                                                   const int             lineIndex,
                                                   const float          *realGridSizeFP,
                                                   float                &fx,
                                                   float                &fy,
                                                   float                &fz)
{
    if (!(order & (order - 1))) // Only for orders of power of 2
    {
        const unsigned int activeMask = c_fullWarpMask;

        // A tricky shuffle reduction inspired by reduce_force_j_warp_shfl
        // TODO: find out if this is the best in terms of transactions count
        static_assert(order == 4, "Only order of 4 is implemented");
        static_assert(atomDataSize <= warp_size, "TODO: rework for atomDataSize > warp_size (order 8 or larger)");
        const int width = atomDataSize;

        fx += gmx_shfl_down_sync(activeMask, fx, 1, width);
        fy += gmx_shfl_up_sync  (activeMask, fy, 1, width);
        fz += gmx_shfl_down_sync(activeMask, fz, 1, width);

        if (splineIndex & 1)
        {
            fx = fy;
        }

        fx += gmx_shfl_down_sync(activeMask, fx, 2, width);
        fz += gmx_shfl_up_sync  (activeMask, fz, 2, width);

        if (splineIndex & 2)
        {
            fx = fz;
        }

        // By now fx contains intermediate quad sums of all 3 components:
        // splineIndex    0            1            2 and 3      4            5            6 and 7      8...
        // sum of...      fx0 to fx3   fy0 to fy3   fz0 to fz3   fx4 to fx7   fy4 to fy7   fz4 to fz7   etc.

        // We have to just further reduce those groups of 4
        for (int delta = 4; delta < atomDataSize; delta <<= 1)
        {
            fx += gmx_shfl_down_sync(activeMask, fx, delta, width);
        }

        const int dimIndex = splineIndex;
        if (dimIndex < DIM)
        {
            const float n = read_grid_size(realGridSizeFP, dimIndex);
            *((float *)(&sm_forces[atomIndexLocal]) + dimIndex) = fx * n;
        }
    }
    else
    {
        // We use blockSize shared memory elements to read fx, or fy, or fz, and then reduce them to fit into smemPerDim elements
        // which are stored separately (first 2 dimensions only)
        const int         smemPerDim   = warp_size;
        const int         smemReserved = (DIM - 1) * smemPerDim;
        __shared__ float  sm_forceReduction[smemReserved + blockSize];
        __shared__ float *sm_forceTemp[DIM];

        const int         numWarps  = blockSize / smemPerDim;
        const int         minStride = max(1, atomDataSize / numWarps); // order 4: 128 threads => 4, 256 threads => 2, etc

#pragma unroll
        for (int dimIndex = 0; dimIndex < DIM; dimIndex++)
        {
            int elementIndex = smemReserved + lineIndex;
            // Store input force contributions
            sm_forceReduction[elementIndex] = (dimIndex == XX) ? fx : (dimIndex == YY) ? fy : fz;
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
            int redStride = minStride;
            if (splineIndex < redStride)
            {
                const int packedIndex = atomIndexLocal * redStride + splineIndex;
                sm_forceTemp[dimIndex][packedIndex] = sm_forceReduction[elementIndex] + sm_forceReduction[elementIndex + redStride];
            }
        }

        __syncthreads();

        assert ((blockSize / warp_size) >= DIM);
        //assert (atomsPerBlock <= warp_size);

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

            const float n = read_grid_size(realGridSizeFP, dimIndex);

            const int   atomIndex = sourceIndex / minStride;
            if (sourceIndex == minStride * atomIndex)
            {
                *((float *)(&sm_forces[atomIndex]) + dimIndex) = (sm_forceTemp[dimIndex][sourceIndex] + sm_forceTemp[dimIndex][sourceIndex + 1]) * n;
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
 * \param[in]  kernelParams         All the PME GPU data.
 */
template <
    const int order,
    const bool overwriteForces,
    const bool wrapX,
    const bool wrapY
    >
__launch_bounds__(c_gatherMaxThreadsPerBlock, c_gatherMinBlocksPerMP)
__global__ void pme_gather_kernel(const PmeGpuCudaKernelParams    kernelParams)
{
    /* Global memory pointers */
    const float * __restrict__  gm_coefficients     = kernelParams.atoms.d_coefficients;
    const float * __restrict__  gm_grid             = kernelParams.grid.d_realGrid;
    float * __restrict__        gm_forces           = kernelParams.atoms.d_forces;

    const float * __restrict__  gm_theta            = kernelParams.atoms.d_theta;
    const float * __restrict__  gm_dtheta           = kernelParams.atoms.d_dtheta;
    const int * __restrict__    gm_gridlineIndices  = kernelParams.atoms.d_gridlineIndices;


    float3 atomX;
    float atomCharge;

    /* Some sizes */
    const int    atomsPerBlock  = (c_gatherMaxThreadsPerBlock / c_pmeGatherThreadsPerAtom);
    const int    blockIndex      = blockIdx.y * gridDim.x + blockIdx.x;
    const int    atomIndexOffset = blockIndex * atomsPerBlock;

    const int    atomDataSize   = c_pmeGatherThreadsPerAtom; /* Number of data components and threads for a single atom */
    const int    atomsPerWarp   = c_pmeGatherAtomsPerWarp;

    const int    blockSize      = atomsPerBlock * atomDataSize;
    assert(blockSize == blockDim.x * blockDim.y * blockDim.z);

    /* These are the atom indices - for the shared and global memory */
    const int         atomIndexLocal    = threadIdx.z;
    const int         atomIndexGlobal   = atomIndexOffset + atomIndexLocal;

    /* Early return for fully empty blocks at the end
     * (should only happen for billions of input atoms)
     */
    if (atomIndexOffset >= kernelParams.atoms.nAtoms)
    {
        return;
    }

    const int         splineParamsSize             = atomsPerBlock * DIM * order; // 96
    const int         gridlineIndicesSize          = atomsPerBlock * DIM; // 24
    __shared__ int    sm_gridlineIndices[gridlineIndicesSize];
    __shared__ float2 sm_splineParams[splineParamsSize]; /* Theta/dtheta pairs  as .x/.y */

//    __shared__ int    sm_gridlineIndices_tmp[gridlineIndicesSize];
//    __shared__ float2 sm_splineParams_tmp[splineParamsSize];

    /* Spline Y/Z coordinates */
//    const int ithy = threadIdx.y;
    const int ithz = threadIdx.x;

    int ithy;

    /* These are the spline contribution indices in shared memory */
//    const int splineIndex = threadIdx.y * blockDim.x + threadIdx.x;                  /* Relative to the current particle , 0..15 for order 4 */
     const int splineIndex = threadIdx.x;
    const int lineIndex   = (threadIdx.z * (blockDim.x * blockDim.y)) + splineIndex; /* And to all the block's particles */

    const int       threadLocalId = lineIndex;
//         (threadIdx.z * (blockDim.x * blockDim.y))
//        + (threadIdx.y * blockDim.x)
//        + threadIdx.x;

/*
    // * Staging the atom gridline indices, DIM * atomsPerBlock threads * /
    const int localGridlineIndicesIndex  = threadLocalId;
    const int globalGridlineIndicesIndex = blockIndex * gridlineIndicesSize + localGridlineIndicesIndex;
    const int globalCheckIndices         = pme_gpu_check_atom_data_index(globalGridlineIndicesIndex, kernelParams.atoms.nAtoms * DIM);
    if ((localGridlineIndicesIndex < gridlineIndicesSize) & globalCheckIndices)
    {
        sm_gridlineIndices[localGridlineIndicesIndex] = gm_gridlineIndices[globalGridlineIndicesIndex];
        assert(sm_gridlineIndices[localGridlineIndicesIndex] >= 0);
    }
    // * Staging the spline parameters, DIM * order * atomsPerBlock threads * /
    const int localSplineParamsIndex  = threadLocalId;
    const int globalSplineParamsIndex = blockIndex * splineParamsSize + localSplineParamsIndex;
    const int globalCheckSplineParams = pme_gpu_check_atom_data_index(globalSplineParamsIndex, kernelParams.atoms.nAtoms * DIM * order);
    if ((localSplineParamsIndex < splineParamsSize) && globalCheckSplineParams)
    {
        sm_splineParams[localSplineParamsIndex].x = gm_theta[globalSplineParamsIndex];
        sm_splineParams[localSplineParamsIndex].y = gm_dtheta[globalSplineParamsIndex];
        assert(isfinite(sm_splineParams[localSplineParamsIndex].x));
        assert(isfinite(sm_splineParams[localSplineParamsIndex].y));
    }
    __syncthreads();
*/

// RECALCULATE SPLINES

// GET CHARGE and positions
   atomCharge=kernelParams.atoms.d_coefficients[atomIndexGlobal];
   atomX.x =  kernelParams.atoms.d_coordinates[ atomIndexGlobal*DIM + XX];
   atomX.y =  kernelParams.atoms.d_coordinates[ atomIndexGlobal*DIM + YY];
   atomX.z =  kernelParams.atoms.d_coordinates[ atomIndexGlobal*DIM + ZZ];

// read local 

/*
const int blockIndex       = blockIdx.y * gridDim.x + blockIdx.x;
const int threadLocalIndex = ((threadIdx.z * blockDim.y + threadIdx.y) * blockDim.x) + threadIdx.x;
const int localIndex       = threadLocalIndex;
const int globalIndexBase  = blockIndex * atomsPerBlock * dataCountPerAtom;
const int globalIndex      = globalIndexBase + localIndex;

So we have (block number)*atomsPerBlock*dataCountPerAtom as the base
then each block is reading atomsPerBlock*dataCountPerAtom

    __syncthreads();
*/

/*
     __shared__ float sm_coefficients[atomsPerBlock];
     pme_gpu_stage_atom_data<float, atomsPerBlock, 1>(kernelParams, sm_coefficients, kernelParams.atoms.d_coefficients);
     __shared__ float sm_coordinates[DIM * atomsPerBlock];
     pme_gpu_stage_atom_data<float, atomsPerBlock, DIM>(kernelParams, sm_coordinates, kernelParams.atoms.d_coordinates);

     __syncthreads();
     calculate_splines_gather<order, atomsPerBlock>(kernelParams, atomIndexOffset, (const float3 *)sm_coordinates,
                                                sm_coefficients, sm_splineParams, sm_gridlineIndices);
 __syncthreads();
*/

calculate_splines_gather<order, atomsPerBlock>(kernelParams, atomIndexOffset, 
&atomX,&atomCharge,
//(const float3 *)sm_coordinates,sm_coefficients, 
sm_splineParams, sm_gridlineIndices);

/*
    if ((localGridlineIndicesIndex < gridlineIndicesSize) && globalCheckIndices)
    {
       assert(sm_gridlineIndices[localGridlineIndicesIndex] == sm_gridlineIndices_tmp[localGridlineIndicesIndex]);
       if(localGridlineIndicesIndex == (gridlineIndicesSize-1))
         if(blockIdx.x==25);
          printf("gridline - %d %d \n",sm_gridlineIndices[localGridlineIndicesIndex],
                sm_gridlineIndices_tmp[localGridlineIndicesIndex]); 
    }

    if ((localSplineParamsIndex < splineParamsSize) && globalCheckSplineParams)
    {
       assert( abs(sm_splineParams[localSplineParamsIndex].x - sm_splineParams_tmp[localSplineParamsIndex].x )  < 0.01);

    }
*/ 

    float           fx = 0.0f;
    float           fy = 0.0f;
    float           fz = 0.0f;

    const int       globalCheck = pme_gpu_check_atom_data_index(atomIndexGlobal, kernelParams.atoms.nAtoms);
    const int       chargeCheck = pme_gpu_check_atom_charge(gm_coefficients[atomIndexGlobal]);

    if (chargeCheck & globalCheck)
    {
        const int    nx        = kernelParams.grid.realGridSize[XX];
        const int    ny        = kernelParams.grid.realGridSize[YY];
        const int    nz        = kernelParams.grid.realGridSize[ZZ];
        const int    pny       = kernelParams.grid.realGridSizePadded[YY];
        const int    pnz       = kernelParams.grid.realGridSizePadded[ZZ];

        const int    atomWarpIndex = atomIndexLocal % atomsPerWarp;
        const int    warpIndex     = atomIndexLocal / atomsPerWarp;

        const int    splineIndexBase = getSplineParamIndexBase<order, atomsPerWarp>(warpIndex, atomWarpIndex);
        int    splineIndexY ;//   = getSplineParamIndex<order, atomsPerWarp>(splineIndexBase, YY, ithy);
        float2 tdy  ;//          = sm_splineParams[splineIndexY];
        const int    splineIndexZ    = getSplineParamIndex<order, atomsPerWarp>(splineIndexBase, ZZ, ithz);
        const float2 tdz             = sm_splineParams[splineIndexZ];

        int iz  = sm_gridlineIndices[atomIndexLocal * DIM + ZZ] + ithz;
        const int    ixBase         = sm_gridlineIndices[atomIndexLocal * DIM + XX];

        if (iz >= nz)
        {
            iz -= nz;
        }
        int constOffset, iy;       
        
        for (int ithy = 0; (ithy < order); ithy++)
        {
          splineIndexY    = getSplineParamIndex<order, atomsPerWarp>(splineIndexBase, YY, ithy);
          tdy = sm_splineParams[splineIndexY];

          iy             = sm_gridlineIndices[atomIndexLocal * DIM + YY] + ithy;
          if (wrapY & (iy >= ny))
          {
              iy -= ny;
          }
          constOffset    = iy * pnz + iz;

#pragma unroll
          for (int ithx = 0; (ithx < order); ithx++)
          {
            int ix = ixBase + ithx;
            if (wrapX & (ix >= nx))
            {
                ix -= nx;
            }
            const int     gridIndexGlobal  = ix * pny * pnz + constOffset;
            assert(gridIndexGlobal >= 0);
            const float   gridValue    = gm_grid[gridIndexGlobal];
            assert(isfinite(gridValue));
            const int     splineIndexX = getSplineParamIndex<order, atomsPerWarp>(splineIndexBase, XX, ithx);
            const float2  tdx          = sm_splineParams[splineIndexX];
            const float   fxy1         = tdz.x * gridValue;
            const float   fz1          = tdz.y * gridValue;
            fx += tdx.y * tdy.x * fxy1;
            fy += tdx.x * tdy.y * fxy1;
            fz += tdx.x * tdy.x * fz1;
          }
        }
    }

    // Reduction of partial force contributions
    __shared__ float3 sm_forces[atomsPerBlock];
    reduce_atom_forces<order, atomDataSize, blockSize>(sm_forces,
                                                       atomIndexLocal, splineIndex, lineIndex,
                                                       kernelParams.grid.realGridSizeFP,
                                                       fx, fy, fz);
    __syncthreads();

    /* Calculating the final forces with no component branching, atomsPerBlock threads */
    const int forceIndexLocal  = threadLocalId;
    const int forceIndexGlobal = atomIndexOffset + forceIndexLocal;
    const int calcIndexCheck   = pme_gpu_check_atom_data_index(forceIndexGlobal, kernelParams.atoms.nAtoms);
    if ((forceIndexLocal < atomsPerBlock) & calcIndexCheck)
    {
        const float3  atomForces               = sm_forces[forceIndexLocal];
        const float   negCoefficient           = -gm_coefficients[forceIndexGlobal];
        float3        result;
        result.x                   = negCoefficient * kernelParams.current.recipBox[XX][XX] * atomForces.x;
        result.y                   = negCoefficient * (kernelParams.current.recipBox[XX][YY] * atomForces.x + kernelParams.current.recipBox[YY][YY] * atomForces.y);
        result.z                   = negCoefficient * (kernelParams.current.recipBox[XX][ZZ] * atomForces.x + kernelParams.current.recipBox[YY][ZZ] * atomForces.y + kernelParams.current.recipBox[ZZ][ZZ] * atomForces.z);
        sm_forces[forceIndexLocal] = result;
    }

    gmx_syncwarp();
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
            int         outputIndexLocal  = i * iterThreads + threadLocalId;
            int         outputIndexGlobal = blockIndex * blockForcesSize + outputIndexLocal;
            const int   globalOutputCheck = pme_gpu_check_atom_data_index(outputIndexGlobal, kernelParams.atoms.nAtoms * DIM);
            if (globalOutputCheck)
            {
                const float outputForceComponent = ((float *)sm_forces)[outputIndexLocal];
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
template __global__ void pme_gather_kernel<4, true, true, true>(const PmeGpuCudaKernelParams);
template __global__ void pme_gather_kernel<4, false, true, true>(const PmeGpuCudaKernelParams);
