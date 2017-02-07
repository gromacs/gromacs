/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013-2016,2017, by the GROMACS development team, led by
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
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 */

#include "gmxpre.h"

#include <cassert>

#include "gromacs/ewald/pme.h"
#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

#include "pme.cuh"
#include "pme-grid.h"
#include "pme-timings.cuh"

/*
 * This define affects the spline calculation behaviour in the kernel.
 * 0: a single GPU thread handles a single dimension of a single particle (calculating and storing (order) spline values and derivatives).
 * 1: (order) threads do redundant work on this same task, each one stores only a single theta and single dtheta into global arrays.
 * The only efficiency difference is less global store operations, countered by more redundant spline computation.
 *
 * TODO: estimate if this should be a boolean parameter (and add it to the unit test if so).
 */
#define PME_GPU_PARALLEL_SPLINE 1

texture<int, 1, cudaReadModeElementType>   gridlineIndicesTableTextureRef;
texture<float, 1, cudaReadModeElementType> fractShiftsTableTextureRef;

void pme_gpu_make_fract_shifts_textures(pme_gpu_t *pmeGpu)
{
    pme_gpu_free_fract_shifts_textures(pmeGpu);

    if (pmeGpu->archSpecific->tableHandling != GpuTableHandling::NoTextures)
    {
        const int                     nx                  = pmeGpu->common->nk[XX];
        const int                     ny                  = pmeGpu->common->nk[YY];
        const int                     nz                  = pmeGpu->common->nk[ZZ];
        const int                     cellCount           = c_pmeNeighborUnitcellCount;
        const int                     newFractShiftsSize  = cellCount * (nx + ny + nz);

        pme_gpu_cuda_kernel_params_t *kernelParamsPtr = pmeGpu->kernelParams.get();

        float                        *fractShiftsArray      = kernelParamsPtr->grid.d_fractShiftsTable;
        int                          *gridlineIndicesArray  = kernelParamsPtr->grid.d_gridlineIndicesTable;

        cudaError_t                   stat;
        if (pmeGpu->archSpecific->tableHandling == GpuTableHandling::TextureObjects)
        {
            cudaResourceDesc rd;
            cudaTextureDesc  td;

            memset(&rd, 0, sizeof(rd));
            rd.resType                  = cudaResourceTypeLinear;
            rd.res.linear.devPtr        = fractShiftsArray;
            rd.res.linear.desc.f        = cudaChannelFormatKindFloat;
            rd.res.linear.desc.x        = 32;
            rd.res.linear.sizeInBytes   = newFractShiftsSize * sizeof(float);
            memset(&td, 0, sizeof(td));
            td.readMode                 = cudaReadModeElementType;
            stat = cudaCreateTextureObject(&kernelParamsPtr->fractShiftsTableTexture, &rd, &td, NULL);
            CU_RET_ERR(stat, "cudaCreateTextureObject on fractShiftsTableTexture failed");

            memset(&rd, 0, sizeof(rd));
            rd.resType                  = cudaResourceTypeLinear;
            rd.res.linear.devPtr        = gridlineIndicesArray;
            rd.res.linear.desc.f        = cudaChannelFormatKindSigned;
            rd.res.linear.desc.x        = 32;
            rd.res.linear.sizeInBytes   = newFractShiftsSize * sizeof(int);
            memset(&td, 0, sizeof(td));
            td.readMode                 = cudaReadModeElementType;
            stat = cudaCreateTextureObject(&kernelParamsPtr->gridlineIndicesTableTexture, &rd, &td, NULL);
            CU_RET_ERR(stat, "cudaCreateTextureObject on gridlineIndicesTableTexture failed");
        }
        else if (pmeGpu->archSpecific->tableHandling == GpuTableHandling::TextureReferences)
        {
            cudaChannelFormatDesc cd_fsh = cudaCreateChannelDesc<float>();
            stat = cudaBindTexture(NULL, &fractShiftsTableTextureRef, fractShiftsArray, &cd_fsh, newFractShiftsSize * sizeof(float));
            CU_RET_ERR(stat, "cudaBindTexture on fractShiftsTableTextureRef failed");

            cudaChannelFormatDesc cd_nn = cudaCreateChannelDesc<int>();
            stat = cudaBindTexture(NULL, &gridlineIndicesTableTextureRef, gridlineIndicesArray, &cd_nn, newFractShiftsSize * sizeof(int));
            CU_RET_ERR(stat, "cudaBindTexture on gridlineIndicesTableTextureRef failed");
        }
    }
}

void pme_gpu_free_fract_shifts_textures(const pme_gpu_t *pmeGpu)
{
    switch (pmeGpu->archSpecific->tableHandling)
    {
        case GpuTableHandling::TextureObjects:
            CU_RET_ERR(cudaDestroyTextureObject(pmeGpu->kernelParams->fractShiftsTableTexture),
                       "cudaDestroyTextureObject on fractShiftsTableTexture failed");
            CU_RET_ERR(cudaDestroyTextureObject(pmeGpu->kernelParams->gridlineIndicesTableTexture),
                       "cudaDestroyTextureObject on gridlineIndicesTableTexture failed");
            break;

        case GpuTableHandling::TextureReferences:
            CU_RET_ERR(cudaUnbindTexture(fractShiftsTableTextureRef),
                       "cudaUnbindTexture on fractShiftsTableTextureRef failed");
            CU_RET_ERR(cudaUnbindTexture(gridlineIndicesTableTextureRef),
                       "cudaUnbindTexture on gridlineIndicesTableTextureRef failed");
            break;

        default:
            break;
    }
}

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
void pme_gpu_stage_atom_data(const pme_gpu_cuda_kernel_params_t kernelParams,
                             T * __restrict__                   sm_destination,
                             const T * __restrict__             gm_source)
{
    // TODO: with padding disabled, this would ignore spline data alignment
    const int threadLocalIndex = ((threadIdx.z * blockDim.y + threadIdx.y) * blockDim.x) + threadIdx.x;
    const int localIndex       = threadLocalIndex;
    const int globalIndexBase  = blockIdx.x * atomsPerBlock * dataCountPerAtom;
    const int globalIndex      = globalIndexBase + localIndex;
    const int globalCheck      = pme_gpu_check_atom_data_index(globalIndex, kernelParams.atoms.nAtoms * dataCountPerAtom);
    if ((localIndex < atomsPerBlock * dataCountPerAtom) & globalCheck)
    {
        assert(!isnan(float(gm_source[globalIndex])));
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
__device__ __forceinline__ void calculate_splines(const pme_gpu_cuda_kernel_params_t     kernelParams,
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

    /* Fractional coordinates */
    __shared__ float sm_fractCoords[atomsPerBlock * DIM];

    /* Thread index w.r.t. block */
    const int                     threadLocalId = (threadIdx.z * (blockDim.x * blockDim.y))
        + (threadIdx.y * blockDim.x)
        + threadIdx.x;
    /* Warp index w.r.t. block - could probably be obtained easier? */
    const int warpIndex         = threadLocalId / warp_size;
    /* Thread index w.r.t. warp */
    const int threadWarpIndex   = threadLocalId % warp_size;
    /* Atom index w.r.t. warp - alternating 0 1 0 1 .. */
    const int atomWarpIndex     = threadWarpIndex % PME_SPREADGATHER_ATOMS_PER_WARP;
    /* Atom index w.r.t. block/shared memory */
    const int atomIndexLocal    = warpIndex * PME_SPREADGATHER_ATOMS_PER_WARP + atomWarpIndex;

    /* Atom index w.r.t. global memory */
    const int atomIndexGlobal   = atomIndexOffset + atomIndexLocal;
    /* Spline contribution index in one dimension */
    const int orderIndex = threadWarpIndex / (PME_SPREADGATHER_ATOMS_PER_WARP * DIM);
    /* Dimension index */
    const int dimIndex = (threadWarpIndex - orderIndex * (PME_SPREADGATHER_ATOMS_PER_WARP * DIM)) / PME_SPREADGATHER_ATOMS_PER_WARP;

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
                    t           = x.x * kernelParams.step.recipBox[dimIndex][XX] + x.y * kernelParams.step.recipBox[dimIndex][YY] + x.z * kernelParams.step.recipBox[dimIndex][ZZ];
                    break;

                case YY:
                    tableIndex  = kernelParams.grid.tablesOffsets[YY];
                    n           = kernelParams.grid.realGridSizeFP[YY];
                    t           = /*x.x * kernelParams.step.recipbox[dimIndex][XX] + */ x.y * kernelParams.step.recipBox[dimIndex][YY] + x.z * kernelParams.step.recipBox[dimIndex][ZZ];
                    break;

                case ZZ:
                    tableIndex  = kernelParams.grid.tablesOffsets[ZZ];
                    n           = kernelParams.grid.realGridSizeFP[ZZ];
                    t           = /*x.x * kernelParams.step.recipbox[dimIndex][XX] + x.y * kernelParams.step.recipbox[dimIndex][YY] + */ x.z * kernelParams.step.recipBox[dimIndex][ZZ];
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

#if PME_USE_TEXTURES
#if PME_USE_TEXOBJ
            sm_fractCoords[sharedMemoryIndex]    += tex1Dfetch<float>(kernelParams.fractShiftsTableTexture, tableIndex);
            sm_gridlineIndices[sharedMemoryIndex] = tex1Dfetch<int>(kernelParams.gridlineIndicesTableTexture, tableIndex);
#else
            sm_fractCoords[sharedMemoryIndex]    += tex1Dfetch(fractShiftsTableTextureRef, tableIndex);
            sm_gridlineIndices[sharedMemoryIndex] = tex1Dfetch(gridlineIndicesTableTextureRef, tableIndex);
#endif
#else
            const float * __restrict__  gm_fractShiftsTable      = kernelParams.grid.d_fractShiftsTable;
            const int * __restrict__    gm_gridlineIndicesTable  = kernelParams.grid.d_gridlineIndicesTable;
            sm_fractCoords[sharedMemoryIndex]    += gm_fractShiftsTable[tableIndex]; // TODO: try __ldg as well
            sm_gridlineIndices[sharedMemoryIndex] = gm_gridlineIndicesTable[tableIndex];
#endif
            gm_gridlineIndices[atomIndexOffset * DIM + sharedMemoryIndex] = sm_gridlineIndices[sharedMemoryIndex];
        }

        /* B-spline calculation */

        const int chargeCheck = pme_gpu_check_atom_charge(sm_coefficients[atomIndexLocal]);
        if (chargeCheck)
        {
            float       div;
            int         o = orderIndex; // This is an index that is set once for PME_GPU_PARALLEL_SPLINE == 1

            const float dr = sm_fractCoords[sharedMemoryIndex];
            assert(!isnan(dr));

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

            const int thetaGlobalOffsetBase = atomIndexOffset * DIM * order;

            /* Differentiation and storing the spline derivatives (dtheta) */
#if !PME_GPU_PARALLEL_SPLINE
            // With PME_GPU_PARALLEL_SPLINE == 1, o is already set to orderIndex;
            // With PME_GPU_PARALLEL_SPLINE == 0, we loop o over range(order).
#pragma unroll
            for (o = 0; o < order; o++)
#endif
            {
                const int   thetaIndex       = PME_SPLINE_THETA_STRIDE * (((o + order * warpIndex) * DIM + dimIndex) * PME_SPREADGATHER_ATOMS_PER_WARP + atomWarpIndex);
                const int   thetaGlobalIndex = thetaGlobalOffsetBase + thetaIndex;

                const float dtheta = ((o > 0) ? SPLINE_DATA(o - 1) : 0.0f) - SPLINE_DATA(o);
                assert(!isnan(dtheta));
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
                const int thetaIndex       = PME_SPLINE_THETA_STRIDE * (((o + order * warpIndex) * DIM + dimIndex) * PME_SPREADGATHER_ATOMS_PER_WARP + atomWarpIndex);
                const int thetaGlobalIndex = thetaGlobalOffsetBase + thetaIndex;

                sm_theta[thetaIndex]       = SPLINE_DATA(o);
                assert(!isnan(sm_theta[thetaIndex]));
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
__device__ __forceinline__ void spread_charges(const pme_gpu_cuda_kernel_params_t     kernelParams,
                                               int                                    atomIndexOffset,
                                               const float * __restrict__             sm_coefficients,
                                               const int * __restrict__               sm_gridlineIndices,
                                               const float * __restrict__             sm_theta)
{
    /* Global memory pointer to the output grid */
    float * __restrict__ gm_grid = kernelParams.grid.d_realGrid;

    const int            nx         = kernelParams.grid.realGridSize[XX];
    const int            ny         = kernelParams.grid.realGridSize[YY];
    const int            nz         = kernelParams.grid.realGridSize[ZZ];
    const int            pny        = kernelParams.grid.realGridSizePadded[YY];
    const int            pnz        = kernelParams.grid.realGridSizePadded[ZZ];

    const int            offx = 0, offy = 0, offz = 0; // unused for now

    const int            atomIndexLocal  = threadIdx.z;
    const int            atomIndexGlobal = atomIndexOffset + atomIndexLocal;

    const int            globalCheck = pme_gpu_check_atom_data_index(atomIndexGlobal, kernelParams.atoms.nAtoms);
    const int            chargeCheck = pme_gpu_check_atom_charge(sm_coefficients[atomIndexLocal]);
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
        const int    atomWarpIndex     = atomIndexLocal % PME_SPREADGATHER_ATOMS_PER_WARP;
        /* Warp index w.r.t. block - could probably be obtained easier? */
        const int    warpIndex         = atomIndexLocal / PME_SPREADGATHER_ATOMS_PER_WARP;
        const int    dimStride         = PME_SPLINE_THETA_STRIDE * PME_SPREADGATHER_ATOMS_PER_WARP;
        const int    orderStride       = dimStride * DIM;
        const int    thetaOffsetBase   = orderStride * order * warpIndex + atomWarpIndex;

        const float  thetaZ         = sm_theta[thetaOffsetBase + ithz * orderStride + ZZ * dimStride];
        const float  thetaY         = sm_theta[thetaOffsetBase + ithy * orderStride + YY * dimStride];
        const float  constVal       = thetaZ * thetaY * sm_coefficients[atomIndexLocal];
        assert(!isnan(constVal));
        const int    constOffset       = iy * pnz + iz;
        const float *sm_thetaX         = sm_theta + (thetaOffsetBase + XX * dimStride);

#pragma unroll
        for (int ithx = 0; (ithx < order); ithx++)
        {
            int ix = ixBase + ithx;
            if (wrapX & (ix >= nx))
            {
                ix -= nx;
            }
            const int gridIndexGlobal = ix * pny * pnz + constOffset;
            assert(!isnan(sm_thetaX[ithx * orderStride]));
            assert(!isnan(gm_grid[gridIndexGlobal]));
            atomicAdd(gm_grid + gridIndexGlobal, sm_thetaX[ithx * orderStride] * constVal);
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
__launch_bounds__(PME_SPREADGATHER_THREADS_PER_BLOCK, PME_MIN_BLOCKS_PER_MP)
__global__ void pme_spline_and_spread_kernel(const pme_gpu_cuda_kernel_params_t kernelParams)
{
    const int                     atomsPerBlock     = PME_SPREADGATHER_ATOMS_PER_BLOCK; // TODO: put it into a better function?
    // Gridline indices, ivec
    __shared__ int                sm_gridlineIndices[atomsPerBlock * DIM];
    // Charges
    __shared__ float              sm_coefficients[atomsPerBlock];
    // Spline values
    __shared__ float              sm_theta[atomsPerBlock * DIM * order];

    const int                     atomIndexOffset = blockIdx.x * atomsPerBlock;

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

void pme_gpu_spread(const pme_gpu_t *pmeGpu,
                    int gmx_unused   gridIndex,
                    real            *h_grid,
                    bool             computeSplines,
                    bool             spreadCharges)
{
    GMX_RELEASE_ASSERT(computeSplines || spreadCharges, "PME spline/spread kernel has invalid input (nothing to do)");

    cudaStream_t                        stream          = pmeGpu->archSpecific->pmeStream;
    const pme_gpu_cuda_kernel_params_t *kernelParamsPtr = pmeGpu->kernelParams.get();

    const int    order                   = pmeGpu->common->pme_order;
    const int    atomsPerBlock           = PME_SPREADGATHER_ATOMS_PER_BLOCK;

    GMX_RELEASE_ASSERT(kernelParamsPtr->atoms.nAtoms > 0, "No atom data in PME GPU spread");
    dim3 nBlocks(pmeGpu->nAtomsPadded / atomsPerBlock);
    dim3 dimBlock(order, order, atomsPerBlock);

    // These should later check for PME decomposition
    const bool wrapX = true;
    const bool wrapY = true;
    GMX_UNUSED_VALUE(wrapX);
    GMX_UNUSED_VALUE(wrapY);
    switch (order)
    {
        case 4:
        {
            // TODO: cleaner unroll with some template trick?
            if (computeSplines)
            {
                if (spreadCharges)
                {
                    pme_gpu_start_timing(pmeGpu, gtPME_SPLINEANDSPREAD);
                    pme_spline_and_spread_kernel<4, true, true, wrapX, wrapY> <<< nBlocks, dimBlock, 0, stream>>> (*kernelParamsPtr);
                    CU_LAUNCH_ERR("pme_spline_and_spread_kernel");
                    pme_gpu_stop_timing(pmeGpu, gtPME_SPLINEANDSPREAD);
                }
                else
                {
                    pme_gpu_start_timing(pmeGpu, gtPME_SPLINE);
                    pme_spline_and_spread_kernel<4, true, false, wrapX, wrapY> <<< nBlocks, dimBlock, 0, stream>>> (*kernelParamsPtr);
                    CU_LAUNCH_ERR("pme_spline_and_spread_kernel");
                    pme_gpu_stop_timing(pmeGpu, gtPME_SPLINE);
                }
            }
            else
            {
                pme_gpu_start_timing(pmeGpu, gtPME_SPREAD);
                pme_spline_and_spread_kernel<4, false, true, wrapX, wrapY> <<< nBlocks, dimBlock, 0, stream>>> (*kernelParamsPtr);
                CU_LAUNCH_ERR("pme_spline_and_spread_kernel");
                pme_gpu_stop_timing(pmeGpu, gtPME_SPREAD);
            }
        }
        break;

        default:
            GMX_THROW(gmx::NotImplementedError("The code for pme_order != 4 was not tested!"));
    }

    const bool copyBackAtomData = computeSplines && (pme_gpu_is_testing(pmeGpu) || !pme_gpu_performs_gather(pmeGpu));
    if (copyBackAtomData)
    {
        const int    alignment     = pme_gpu_get_atom_spline_data_alignment(pmeGpu);
        const size_t nAtomsPadded  = ((pmeGpu->nAtomsAlloc + alignment - 1) / alignment) * alignment;
        const size_t splinesSize   = DIM * nAtomsPadded * order * sizeof(float);
        cu_copy_D2H_async(pmeGpu->staging.h_dtheta, kernelParamsPtr->atoms.d_dtheta, splinesSize, stream);
        cu_copy_D2H_async(pmeGpu->staging.h_theta, kernelParamsPtr->atoms.d_theta, splinesSize, stream);
        cu_copy_D2H_async(pmeGpu->staging.h_gridlineIndices, kernelParamsPtr->atoms.d_gridlineIndices, kernelParamsPtr->atoms.nAtoms * DIM * sizeof(int), stream);
        cudaError_t stat = cudaEventRecord(pmeGpu->archSpecific->syncSplineAtomDataD2H, stream);
        CU_RET_ERR(stat, "PME spread atom data sync event record failure");
    }
    const bool copyBackGrid     = spreadCharges && (pme_gpu_is_testing(pmeGpu) || !pme_gpu_performs_FFT(pmeGpu));
    if (copyBackGrid)
    {
        const size_t gridSize = pmeGpu->archSpecific->realGridSize * sizeof(float);
        cu_copy_D2H_async(h_grid, kernelParamsPtr->grid.d_realGrid, gridSize, stream);
        cudaError_t  stat = cudaEventRecord(pmeGpu->archSpecific->syncSpreadGridD2H, stream);
        CU_RET_ERR(stat, "PME spread grid sync event record failure");
    }
}
