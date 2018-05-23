/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013-2016,2017,2018, by the GROMACS development team, led by
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
#include "pme-gpu-utils.h"
#include "pme-grid.h"

/*
 * This define affects the spline calculation behaviour in the kernel.
 * 0: a single GPU thread handles a single dimension of a single particle (calculating and storing (order) spline values and derivatives).
 * 1: (order) threads do redundant work on this same task, each one stores only a single theta and single dtheta into global arrays.
 * The only efficiency difference is less global store operations, countered by more redundant spline computation.
 *
 * TODO: estimate if this should be a boolean parameter (and add it to the unit test if so).
 */
#define PME_GPU_PARALLEL_SPLINE 0


/*! \brief
 * General purpose function for loading atom-related data from global to shared memory.
 *
 * \tparam[in] T                 Data type (float/int/...)
 * \tparam[in] atomsPerBlock     Number of atoms processed by a block - should be accounted for in the size of the shared memory array.
 * \tparam[in] dataCountPerAtom  Number of data elements per single atom (e.g. DIM for an rvec coordinates array).
 * \param[in]  kernelParams      Input PME CUDA data in constant memory.
 * \param[out] sm_destination    Shared memory array for output.
 * \param[in]  gm_source         Global memory array for input.
 */
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
    if ((localIndex < atomsPerBlock * dataCountPerAtom) & globalCheck)
    {
        assert(isfinite(float(gm_source[globalIndex])));
        sm_destination[localIndex] = gm_source[globalIndex];
    }
}

/*! \brief
 * PME GPU spline parameter and gridline indices calculation.
 * This corresponds to the CPU functions calc_interpolation_idx() and make_bsplines().
 * First stage of the whole kernel.
 *
 * \tparam[in] order                PME interpolation order.
 * \tparam[in] atomsPerBlock        Number of atoms processed by a block - should be accounted for
 *                                  in the sizes of the shared memory arrays.
 * \param[in]  kernelParams         Input PME CUDA data in constant memory.
 * \param[in]  atomIndexOffset      Starting atom index for the execution block w.r.t. global memory.
 * \param[in]  sm_coordinates       Atom coordinates in the shared memory.
 * \param[in]  sm_coefficients      Atom charges/coefficients in the shared memory.
 * \param[out] sm_theta             Atom spline values in the shared memory.
 * \param[out] sm_gridlineIndices   Atom gridline indices in the shared memory.
 */
template <const int order,
          const int atomsPerBlock>
__device__ __forceinline__ void calculate_splines(const PmeGpuCudaKernelParams           kernelParams,
                                                  const int                              atomIndexOffset,
                                                  const float3 * __restrict__            sm_coordinates,
                                                  const float * __restrict__             sm_coefficients,
                                                  float * __restrict__                   sm_theta,
                                                  int * __restrict__                     sm_gridlineIndices)
{
    /* Global memory pointers for output */
    float * __restrict__ gm_theta           = kernelParams.atoms.d_theta;
    float * __restrict__ gm_dtheta          = kernelParams.atoms.d_dtheta;
    int * __restrict__   gm_gridlineIndices = kernelParams.atoms.d_gridlineIndices;

    const int            atomsPerWarp = PME_SPREADGATHER_ATOMS_PER_WARP;

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
    const int atomWarpIndex = threadWarpIndex % atomsPerWarp;
    /* Atom index w.r.t. block/shared memory */
    const int atomIndexLocal = warpIndex * atomsPerWarp + atomWarpIndex;

    /* Atom index w.r.t. global memory */
    const int atomIndexGlobal = atomIndexOffset + atomIndexLocal;
    /* Spline contribution index in one dimension */
    const int orderIndex = threadWarpIndex / (atomsPerWarp * DIM);
    /* Dimension index */
    const int dimIndex = (threadWarpIndex - orderIndex * (atomsPerWarp * DIM)) / atomsPerWarp;

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
            const float3  x = sm_coordinates[atomIndexLocal];
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
            gm_gridlineIndices[atomIndexOffset * DIM + sharedMemoryIndex] = sm_gridlineIndices[sharedMemoryIndex];
        }

        /* B-spline calculation */

        const int chargeCheck = pme_gpu_check_atom_charge(sm_coefficients[atomIndexLocal]);
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
                const int   thetaGlobalIndex = thetaGlobalOffsetBase + thetaIndex;

                const float dtheta = ((o > 0) ? SPLINE_DATA(o - 1) : 0.0f) - SPLINE_DATA(o);
                assert(isfinite(dtheta));
                gm_dtheta[thetaGlobalIndex] = dtheta;
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

                sm_theta[thetaIndex]       = SPLINE_DATA(o);
                assert(isfinite(sm_theta[thetaIndex]));
                gm_theta[thetaGlobalIndex] = SPLINE_DATA(o);
            }
        }
    }
}

/*! \brief
 * Charge spreading onto the grid.
 * This corresponds to the CPU function spread_coefficients_bsplines_thread().
 * Second stage of the whole kernel.
 *
 * \tparam[in] order                PME interpolation order.
 * \tparam[in] wrapX                A boolean which tells if the grid overlap in dimension X should be wrapped.
 * \tparam[in] wrapY                A boolean which tells if the grid overlap in dimension Y should be wrapped.
 * \param[in]  kernelParams         Input PME CUDA data in constant memory.
 * \param[in]  atomIndexOffset      Starting atom index for the execution block w.r.t. global memory.
 * \param[in]  sm_coefficients      Atom coefficents/charges in the shared memory.
 * \param[in]  sm_gridlineIndices   Atom gridline indices in the shared memory.
 * \param[in]  sm_theta             Atom spline values in the shared memory.
 */
template <
    const int order, const bool wrapX, const bool wrapY>
__device__ __forceinline__ void spread_charges(const PmeGpuCudaKernelParams           kernelParams,
                                               int                                    atomIndexOffset,
                                               const float * __restrict__             sm_coefficients,
                                               const int * __restrict__               sm_gridlineIndices,
                                               const float * __restrict__             sm_theta)
{
    /* Global memory pointer to the output grid */
    float * __restrict__ gm_grid = kernelParams.grid.d_realGrid;


    const int atomsPerWarp = PME_SPREADGATHER_ATOMS_PER_WARP;

    const int nx  = kernelParams.grid.realGridSize[XX];
    const int ny  = kernelParams.grid.realGridSize[YY];
    const int nz  = kernelParams.grid.realGridSize[ZZ];
    const int pny = kernelParams.grid.realGridSizePadded[YY];
    const int pnz = kernelParams.grid.realGridSizePadded[ZZ];

    const int offx = 0, offy = 0, offz = 0; // unused for now

    const int atomIndexLocal  = threadIdx.z;
    const int atomIndexGlobal = atomIndexOffset + atomIndexLocal;

    const int globalCheck = pme_gpu_check_atom_data_index(atomIndexGlobal, kernelParams.atoms.nAtoms);
    const int chargeCheck = pme_gpu_check_atom_charge(sm_coefficients[atomIndexLocal]);
    if (chargeCheck & globalCheck)
    {
        // Spline Y/Z coordinates
        const int ithy   = threadIdx.y;
        const int ithz   = threadIdx.x;
        const int ixBase = sm_gridlineIndices[atomIndexLocal * DIM + XX] - offx;
        int       iy     = sm_gridlineIndices[atomIndexLocal * DIM + YY] - offy + ithy;
        if (wrapY & (iy >= ny))
        {
            iy -= ny;
        }
        int iz  = sm_gridlineIndices[atomIndexLocal * DIM + ZZ] - offz + ithz;
        if (iz >= nz)
        {
            iz -= nz;
        }

        /* Atom index w.r.t. warp - alternating 0 1 0 1 .. */
        const int    atomWarpIndex   = atomIndexLocal % atomsPerWarp;
        /* Warp index w.r.t. block - could probably be obtained easier? */
        const int    warpIndex       = atomIndexLocal / atomsPerWarp;

        const int    splineIndexBase = getSplineParamIndexBase<order, atomsPerWarp>(warpIndex, atomWarpIndex);
        const int    splineIndexZ    = getSplineParamIndex<order, atomsPerWarp>(splineIndexBase, ZZ, ithz);
        const float  thetaZ          = sm_theta[splineIndexZ];
        const int    splineIndexY    = getSplineParamIndex<order, atomsPerWarp>(splineIndexBase, YY, ithy);
        const float  thetaY          = sm_theta[splineIndexY];
        const float  constVal        = thetaZ * thetaY * sm_coefficients[atomIndexLocal];
        assert(isfinite(constVal));
        const int    constOffset     = iy * pnz + iz;

#pragma unroll
        for (int ithx = 0; (ithx < order); ithx++)
        {
            int ix = ixBase + ithx;
            if (wrapX & (ix >= nx))
            {
                ix -= nx;
            }
            const int   gridIndexGlobal = ix * pny * pnz + constOffset;
            const int   splineIndexX    = getSplineParamIndex<order, atomsPerWarp>(splineIndexBase, XX, ithx);
            const float thetaX          = sm_theta[splineIndexX];
            assert(isfinite(thetaX));
            assert(isfinite(gm_grid[gridIndexGlobal]));
            atomicAdd(gm_grid + gridIndexGlobal, thetaX * constVal);
        }
    }
}

/*! \brief
 * A spline computation and charge spreading kernel function.
 *
 * \tparam[in] order                PME interpolation order.
 * \tparam[in] computeSplines       A boolean which tells if the spline parameter and
 *                                  gridline indices' computation should be performed.
 * \tparam[in] spreadCharges        A boolean which tells if the charge spreading should be performed.
 * \tparam[in] wrapX                A boolean which tells if the grid overlap in dimension X should be wrapped.
 * \tparam[in] wrapY                A boolean which tells if the grid overlap in dimension Y should be wrapped.
 * \param[in]  kernelParams         Input PME CUDA data in constant memory.
 */
template <
    const int order,
    const bool computeSplines,
    const bool spreadCharges,
    const bool wrapX,
    const bool wrapY
    >
__launch_bounds__(c_spreadMaxThreadsPerBlock)
__global__ void pme_spline_and_spread_kernel(const PmeGpuCudaKernelParams kernelParams)
{
    const int        atomsPerBlock = c_spreadMaxThreadsPerBlock / PME_SPREADGATHER_THREADS_PER_ATOM;
    // Gridline indices, ivec
    __shared__ int   sm_gridlineIndices[atomsPerBlock * DIM];
    // Charges
    __shared__ float sm_coefficients[atomsPerBlock];
    // Spline values
    __shared__ float sm_theta[atomsPerBlock * DIM * order];

    const int        blockIndex      = blockIdx.y * gridDim.x + blockIdx.x;
    const int        atomIndexOffset = blockIndex * atomsPerBlock;

    /* Early return for fully empty blocks at the end
     * (should only happen on Fermi or billions of input atoms)
     */
    if (atomIndexOffset >= kernelParams.atoms.nAtoms)
    {
        return;
    }

    /* Staging coefficients/charges for both spline and spread */
    pme_gpu_stage_atom_data<float, atomsPerBlock, 1>(kernelParams, sm_coefficients, kernelParams.atoms.d_coefficients);

    if (computeSplines)
    {
        /* Staging coordinates */
        __shared__ float sm_coordinates[DIM * atomsPerBlock];
        pme_gpu_stage_atom_data<float, atomsPerBlock, DIM>(kernelParams, sm_coordinates, kernelParams.atoms.d_coordinates);

        __syncthreads();
        calculate_splines<order, atomsPerBlock>(kernelParams, atomIndexOffset, (const float3 *)sm_coordinates,
                                                sm_coefficients, sm_theta, sm_gridlineIndices);
        gmx_syncwarp();
    }
    else
    {
        /* Staging the data for spread
         * (the data is assumed to be in GPU global memory with proper layout already,
         * as in after running the spline kernel)
         */
        /* Spline data - only thetas (dthetas will only be needed in gather) */
        pme_gpu_stage_atom_data<float, atomsPerBlock, DIM * order>(kernelParams, sm_theta, kernelParams.atoms.d_theta);
        /* Gridline indices */
        pme_gpu_stage_atom_data<int, atomsPerBlock, DIM>(kernelParams, sm_gridlineIndices, kernelParams.atoms.d_gridlineIndices);

        __syncthreads();
    }

    /* Spreading */
    if (spreadCharges)
    {
        spread_charges<order, wrapX, wrapY>(kernelParams, atomIndexOffset, sm_coefficients, sm_gridlineIndices, sm_theta);
    }
}

//! Kernel instantiations
template __global__ void pme_spline_and_spread_kernel<4, true, true, true, true>(const PmeGpuCudaKernelParams);
template __global__ void pme_spline_and_spread_kernel<4, true, false, true, true>(const PmeGpuCudaKernelParams);
template __global__ void pme_spline_and_spread_kernel<4, false, true, true, true>(const PmeGpuCudaKernelParams);
