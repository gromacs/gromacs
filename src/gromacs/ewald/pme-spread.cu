/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013-2016, by the GROMACS development team, led by
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
#include "gromacs/utility/gmxassert.h"

#include "pme.cuh"
#include "pme-grid.h"
#include "pme-timings.cuh"

#define PME_GPU_PARALLEL_SPLINE 1
/* This define affects the spline calculation behaviour in the spreading kernel.
 * 0: a single GPU thread handles a single dimension of a single particle (calculating and storing (order) spline values and derivatives).
 * 1: (order) threads do redundant work on this same task, each one stores only a single theta and single dtheta into global arrays.
 * The only efficiency difference is less global store operations, countered by more redundant spline computation.
 */

#if PME_USE_TEXTURES
texture<int, 1, cudaReadModeElementType>   gridlineIndicesTableTextureRef;
texture<float, 1, cudaReadModeElementType> fractShiftsTableTextureRef;
#endif

/* This is the PME GPU spline calculation.
 * It corresponds to the CPU codepath functions calc_interpolation_idx and make_bsplines.
 */
template <const int order>
__device__ __forceinline__ void calculate_splines(const float3 * __restrict__            sm_coordinates,
                                                  float * __restrict__                   sm_coefficients,
                                                  float * __restrict__                   sm_theta,
                                                  int * __restrict__                     sm_gridlineIndices,
                                                  const pme_gpu_cuda_kernel_params_t     kernelParams,
                                                  const int                              globalIndexCalc,
                                                  const int                              localIndexCalc,
                                                  const int                              globalIndexBase,
                                                  const int                              dimIndex,
                                                  const int                              orderIndex)
{
    constexpr int        atomsPerBlock     = PME_SPREADGATHER_ATOMS_PER_BLOCK;
    /* Global memory pointers */
    float * __restrict__ gm_theta           = kernelParams.atoms.d_theta;
    float * __restrict__ gm_dtheta          = kernelParams.atoms.d_dtheta;
    int * __restrict__   gm_gridlineIndices = kernelParams.atoms.d_gridlineIndices;

    /* Fractional coordinates */
    __shared__ float sm_fractCoords[atomsPerBlock * DIM];

    const int        sharedMemoryIndex = localIndexCalc * DIM + dimIndex;

    const int        dataSize   = PME_GPU_PARALLEL_SPLINE ? (atomsPerBlock * DIM) : 1;
    const int        dataOffset = PME_GPU_PARALLEL_SPLINE ? sharedMemoryIndex : 0;
    /* Spline parameter storage, marked as shared for PME_GPU_PARALLEL_SPLINE == 1 to not overuse the local memory */
#if PME_GPU_PARALLEL_SPLINE
    __shared__
#endif
    float data[dataSize * order];

    const int localCheck  = (dimIndex < DIM) && (orderIndex < (PME_GPU_PARALLEL_SPLINE ? order : 1));
    const int globalCheck = pme_gpu_check_atom_data_index(globalIndexCalc, kernelParams.atoms.nAtoms);

    if (localCheck && globalCheck)
    {
        /* Indices interpolation */

        if (orderIndex == 0)
        {
            int           tableIndex, tInt;
            float         n, t;
            const float3  x = sm_coordinates[localIndexCalc];
            /* Accessing fields in fshOffset/nXYZ/recipbox/... with dimIndex offset
             * puts them into local memory(!) instead of accessing the constant memory directly.
             * That's the reason for the switch, to unroll explicitly.
             * The commented parts correspond to the 0 components of the recipbox.
             */
            switch (dimIndex)
            {
                case XX:
                    tableIndex  = kernelParams.grid.tablesOffsets[XX];
                    n           = kernelParams.grid.localGridSizeFP[XX];
                    t           = x.x * kernelParams.step.recipBox[dimIndex][XX] + x.y * kernelParams.step.recipBox[dimIndex][YY] + x.z * kernelParams.step.recipBox[dimIndex][ZZ];
                    break;

                case YY:
                    tableIndex  = kernelParams.grid.tablesOffsets[YY];
                    n           = kernelParams.grid.localGridSizeFP[YY];
                    t           = /*x.x * kernelParams.step.recipbox[dimIndex][XX] + */ x.y * kernelParams.step.recipBox[dimIndex][YY] + x.z * kernelParams.step.recipBox[dimIndex][ZZ];
                    break;

                case ZZ:
                    tableIndex  = kernelParams.grid.tablesOffsets[ZZ];
                    n           = kernelParams.grid.localGridSizeFP[ZZ];
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
            gm_gridlineIndices[globalIndexBase * DIM + sharedMemoryIndex] = sm_gridlineIndices[sharedMemoryIndex];
        }

        /* B-spline calculation */

        const int chargeCheck = pme_gpu_check_atom_charge(sm_coefficients[localIndexCalc]);
        if (chargeCheck)
        {
            float       div;
            int         k;

            const float dr = sm_fractCoords[sharedMemoryIndex];
            assert(!isnan(dr));

            /* dr is relative offset from lower cell limit */
            data[(order - 1) * dataSize + dataOffset] = 0.0f;
            data[1 * dataSize + dataOffset]           = dr;
            data[0 * dataSize + dataOffset]           = 1.0f - dr;

#pragma unroll
            for (int k = 3; k < order; k++)
            {
                div         = 1.0f / (k - 1.0f);
                data[(k - 1) * dataSize + dataOffset] = div * dr * data[(k - 2) * dataSize + dataOffset];
#pragma unroll
                for (int l = 1; l < (k - 1); l++)
                {
                    data[(k - l - 1) * dataSize + dataOffset] =
                        div * ((dr + l) * data[(k - l - 2) * dataSize + dataOffset] +
                               (k - l - dr) * data[(k - l - 1) * dataSize + dataOffset]);
                }
                data[0 * dataSize + dataOffset] = div * (1.0f - dr) * data[0 * dataSize + dataOffset];
            }

            const int particleWarpIndex = localIndexCalc % PME_SPREADGATHER_ATOMS_PER_WARP;
            const int warpIndex         = localIndexCalc / PME_SPREADGATHER_ATOMS_PER_WARP;

            const int thetaGlobalOffsetBase = globalIndexBase * DIM * order;

            /* Differentiation and storing the spline derivatives (dtheta) */
#if PME_GPU_PARALLEL_SPLINE
            k = orderIndex;
#else
#pragma unroll
            for (k = 0; k < order; k++)
#endif
            {
                const int   thetaIndex       = PME_SPLINE_THETA_STRIDE * (((k + order * warpIndex) * DIM + dimIndex) * PME_SPREADGATHER_ATOMS_PER_WARP + particleWarpIndex);
                const int   thetaGlobalIndex = thetaGlobalOffsetBase + thetaIndex;

                const float dtheta = ((k > 0) ? data[(k - 1) * dataSize + dataOffset] : 0.0f) - data[k * dataSize + dataOffset];
                assert(!isnan(dtheta));
                gm_dtheta[thetaGlobalIndex] = dtheta;
            }

            div             = 1.0f / (order - 1);
            data[(order - 1) * dataSize + dataOffset] = div * dr * data[(order - 2) * dataSize + dataOffset];
#pragma unroll
            for (int l = 1; l < (order - 1); l++)
            {
                data[(order - l - 1) * dataSize + dataOffset] = div * ((dr + l) * data[(order - l - 2) * dataSize + dataOffset] + (order - l - dr) * data[(order - l - 1) * dataSize + dataOffset]);
            }
            data[0 * dataSize + dataOffset] = div * (1.0f - dr) * data[0 * dataSize + dataOffset];

            /* Storing the spline values (theta) */
#if PME_GPU_PARALLEL_SPLINE
            k = orderIndex;
#else
#pragma unroll
            for (k = 0; k < order; k++)
#endif
            {
                const int thetaIndex       = PME_SPLINE_THETA_STRIDE * (((k + order * warpIndex) * DIM + dimIndex) * PME_SPREADGATHER_ATOMS_PER_WARP + particleWarpIndex);
                const int thetaGlobalIndex = thetaGlobalOffsetBase + thetaIndex;

                sm_theta[thetaIndex]       = data[k * dataSize + dataOffset];
                assert(!isnan(sm_theta[thetaIndex]));
                gm_theta[thetaGlobalIndex] = data[k * dataSize + dataOffset];
            }
        }
    }
}

template <
    const int order, const bool wrapX, const bool wrapY> // z is always wrapped
__device__ __forceinline__ void spread_charges(const float * __restrict__             sm_coefficient,
                                               const pme_gpu_cuda_kernel_params_t     kernelParams,
                                               const int                              globalIndex,
                                               const int                              localIndex,
                                               const int * __restrict__               sm_gridlineIndices,
                                               const float * __restrict__             sm_theta)
{
    /* Global memory pointer */
    float * __restrict__ gm_grid = kernelParams.grid.d_realGrid;

    const int            nx         = kernelParams.grid.localGridSize[XX];
    const int            ny         = kernelParams.grid.localGridSize[YY];
    const int            nz         = kernelParams.grid.localGridSize[ZZ];
    const int            pny        = kernelParams.grid.localGridSizePadded[YY];
    const int            pnz        = kernelParams.grid.localGridSizePadded[ZZ];

    const int            offx = 0, offy = 0, offz = 0;
    // unused for now

    const int globalCheck = pme_gpu_check_atom_data_index(globalIndex, kernelParams.atoms.nAtoms);
    const int chargeCheck = pme_gpu_check_atom_charge(sm_coefficient[localIndex]);
    if (chargeCheck & globalCheck)
    {
        // spline Y/Z coordinates
        const int    ithy      = threadIdx.y;
        const int    ithz      = threadIdx.x; //?
        const int    ixBase    = sm_gridlineIndices[localIndex * DIM + XX] - offx;
        const int    iy        = sm_gridlineIndices[localIndex * DIM + YY] - offy + ithy;
        const int    iyWrapped = (wrapY & (iy >= ny)) ? (iy - ny) : iy;
        const int    iz        = sm_gridlineIndices[localIndex * DIM + ZZ] - offz + ithz;
        const int    izWrapped = (iz >= nz) ? (iz - nz) : iz;                          //TODO: remove excess variables

        const int    particleWarpIndex = localIndex % PME_SPREADGATHER_ATOMS_PER_WARP; // index of particle w.r.t. the warp (so, 0 or 1)
        const int    warpIndex         = localIndex / PME_SPREADGATHER_ATOMS_PER_WARP; // should be just a normal warp index, actually!
        const int    dimStride         = PME_SPLINE_THETA_STRIDE * PME_SPREADGATHER_ATOMS_PER_WARP;
        const int    orderStride       = dimStride * DIM;
        const int    thetaOffsetBase   = orderStride * order * warpIndex + particleWarpIndex;

        const float  thz         = sm_theta[thetaOffsetBase + ithz * orderStride + ZZ * dimStride];
        const float  thy         = sm_theta[thetaOffsetBase + ithy * orderStride + YY * dimStride];
        const float  constVal    = thz * thy * sm_coefficient[localIndex];
        assert(!isnan(constVal));
        const int    constOffset    = iyWrapped * pnz + izWrapped;
        const float *sm_thx         = sm_theta + (thetaOffsetBase + XX * dimStride);

#pragma unroll
        for (int ithx = 0; (ithx < order); ithx++)
        {
            const int ix          = ixBase + ithx;
            const int ixWrapped   = (wrapX & (ix >= nx)) ? (ix - nx) : ix;
            const int globalIndex = ixWrapped * pny * pnz + constOffset;
            assert(!isnan(sm_thx[ithx * orderStride]));
            assert(!isnan(gm_grid[globalIndex]));
            atomicAdd(gm_grid + globalIndex, sm_thx[ithx * orderStride] * constVal);
        }
    }
}

template <const int atomsPerBlock>
__device__ __forceinline__ void stage_charges(const int                          threadLocalId,
                                              float * __restrict__               sm_coefficient,
                                              const pme_gpu_cuda_kernel_params_t kernelParams)
{
    /* Global memory pointer */
    const float * __restrict__   gm_coefficient = kernelParams.atoms.d_coefficients;

    const int                    globalIndexBase = blockIdx.x * atomsPerBlock;
    const int                    localIndex      = threadLocalId;
    const int                    globalIndex     = globalIndexBase + localIndex;
    const int                    globalCheck     = pme_gpu_check_atom_data_index(globalIndex, kernelParams.atoms.nAtoms);
    if ((localIndex < atomsPerBlock) & globalCheck)
    {
        assert(!isnan(gm_coefficient[globalIndex]));
        sm_coefficient[localIndex] = gm_coefficient[globalIndex];
    }
}

/*! \brief
 * Loads atom coordinates from global to shared memory.
 *
 * \tparam[in] atomsPerBlock   Number of atoms processed by a block - reflects the size of coordinates array.
 * \param[in]  threadLocalId   Index of thread with respect to the block.
 * \param[out] sm_coordinates  Shared memory array for atom coordinates - (atomsPerBlock * DIM) floats.
 * \param[in]  kernelParams    Input PME CUDA data in constant memory.
 */
template <const int atomsPerBlock>
__device__ __forceinline__ void stage_coordinates(int                                threadLocalId,
                                                  float * __restrict__               sm_coordinates,
                                                  const pme_gpu_cuda_kernel_params_t kernelParams)
{
    /* Global memory pointer */
    const float * __restrict__  gm_coordinates = kernelParams.atoms.d_coordinates;

    const int                   globalIndexBase = blockIdx.x * atomsPerBlock * DIM;
    const int                   localIndex      = threadLocalId - 1 * atomsPerBlock;
    const int                   globalIndex     = globalIndexBase + localIndex;
    const int                   globalCheck     = pme_gpu_check_atom_data_index(globalIndex, DIM * kernelParams.atoms.nAtoms); /* DIM floats per atom */
    if ((localIndex >= 0) && (localIndex < DIM * atomsPerBlock) && globalCheck)
    {
        sm_coordinates[localIndex] = gm_coordinates[globalIndex];
    }
}

template <
    const int order,
    const bool computeSplines,     // first part
    const bool spreadCharges,      // second part
    const bool wrapX,
    const bool wrapY
    >
__launch_bounds__(PME_SPREADGATHER_THREADS_PER_BLOCK, PME_MIN_BLOCKS_PER_MP)
__global__ void pme_spline_and_spread_kernel(const pme_gpu_cuda_kernel_params_t kernelParams)
{
    const int                     atomsPerBlock     = PME_SPREADGATHER_ATOMS_PER_BLOCK;
    // gridline indices, ivec
    __shared__ int                sm_gridlineIndices[atomsPerBlock * DIM];
    // charges
    __shared__ float              sm_coefficients[atomsPerBlock];
    // spline parameters
    __shared__ float              sm_theta[atomsPerBlock * DIM * order];

    const int                     threadLocalId = (threadIdx.z * (blockDim.x * blockDim.y))
        + (threadIdx.y * blockDim.x)
        + threadIdx.x;

    const int globalAtomIndexBase = blockIdx.x * atomsPerBlock;

    const int warpIndex         = threadLocalId / warp_size;
    const int threadWarpIndex   = threadLocalId % warp_size;
    const int particleWarpIndex = threadWarpIndex % PME_SPREADGATHER_ATOMS_PER_WARP;
    const int localCalcIndex    = warpIndex * PME_SPREADGATHER_ATOMS_PER_WARP + particleWarpIndex;
    const int globalCalcIndex   = globalAtomIndexBase + localCalcIndex;
    const int orderCalcIndex    = threadWarpIndex / (PME_SPREADGATHER_ATOMS_PER_WARP * DIM);
    const int dimCalcIndex      = (threadWarpIndex - orderCalcIndex * (PME_SPREADGATHER_ATOMS_PER_WARP * DIM)) / PME_SPREADGATHER_ATOMS_PER_WARP;

    /* Staging charges for both spline and spread */
    stage_charges<atomsPerBlock>(threadLocalId, sm_coefficients, kernelParams);

    if (computeSplines)
    {
        // coordinates
        __shared__ float sm_coordinates[DIM * atomsPerBlock];
        stage_coordinates<atomsPerBlock>(threadLocalId, sm_coordinates, kernelParams);
        __syncthreads();
        calculate_splines<order>((const float3 *)sm_coordinates, sm_coefficients,
                                 sm_theta, sm_gridlineIndices,
                                 kernelParams,
                                 globalCalcIndex,
                                 localCalcIndex,
                                 globalAtomIndexBase,
                                 dimCalcIndex,
                                 orderCalcIndex);
        if (!spreadCharges)
        {
            if (threadLocalId < atomsPerBlock)
            {
                sm_coefficients[threadLocalId] = -3253344453.0f;
            }
            if (threadLocalId < atomsPerBlock * DIM)
            {
                sm_gridlineIndices[threadLocalId] = -3254453;
            }
            if (threadLocalId < atomsPerBlock * DIM * order)
            {
                sm_theta[threadLocalId] = -6453254453.0f;
            }
            // fuck everything up!
        }
    }
    else
    {
        /* Staging all the data for spread (the data is assumed to be in GPU global memory already, as in after running the spline kernel) */
        /* TODO: maybe reuse these as functions in gather */

        /* Spline data */
        const float * __restrict__ gm_theta               = kernelParams.atoms.d_theta;
        const int                  splineDataCountPerAtom = DIM * order;
        const int                  localSplineIndex       = threadLocalId;
        const int                  globalSplineIndex      = globalAtomIndexBase * splineDataCountPerAtom + localSplineIndex;
        const int                  globalCheckSplines     = pme_gpu_check_atom_data_index(globalSplineIndex, kernelParams.atoms.nAtoms * splineDataCountPerAtom);
        if ((localSplineIndex < atomsPerBlock * splineDataCountPerAtom) && globalCheckSplines)
        {
            sm_theta[localSplineIndex] = gm_theta[globalSplineIndex];
        }

        /* Gridline indices */
        const int * __restrict__ gm_gridlineIndices      = kernelParams.atoms.d_gridlineIndices;
        const int                indicesDataCountPerAtom = DIM;
        const int                localIndicesIndex       = threadLocalId;
        const int                globalIndicesIndex      = globalAtomIndexBase * indicesDataCountPerAtom + localIndicesIndex;
        const int                globalCheckIndices      = pme_gpu_check_atom_data_index(globalIndicesIndex, kernelParams.atoms.nAtoms * indicesDataCountPerAtom);
        if ((localIndicesIndex < atomsPerBlock * indicesDataCountPerAtom) && globalCheckIndices)
        {
            sm_gridlineIndices[localIndicesIndex] = gm_gridlineIndices[globalIndicesIndex];
        }

        __syncthreads();
    }

    // SPREAD
    if (spreadCharges)
    {
        const int localSpreadIndex  = threadIdx.z;
        const int globalSpreadIndex = globalAtomIndexBase + localSpreadIndex;

        spread_charges<order, wrapX, wrapY>(sm_coefficients, kernelParams, globalSpreadIndex, localSpreadIndex,
                                            sm_gridlineIndices, sm_theta);
    }
}

void pme_gpu_make_fract_shifts_textures(pme_gpu_t *pmeGPU)
{
    pme_gpu_free_fract_shifts_textures(pmeGPU);

#if PME_USE_TEXTURES
    const int                     nx                  = pmeGPU->common->nk[XX];
    const int                     ny                  = pmeGPU->common->nk[YY];
    const int                     nz                  = pmeGPU->common->nk[ZZ];
    const int                     cellCount           = c_pmeNeighborUnitcellCount;
    const int                     newFractShiftsSize  = cellCount * (nx + ny + nz);

    pme_gpu_cuda_kernel_params_t *kernelParamsPtr = pmeGPU->kernelParams.get();

    float                        *fractShiftsArray      = kernelParamsPtr->grid.d_fractShiftsTable;
    int                          *gridlineIndicesArray  = kernelParamsPtr->grid.d_gridlineIndicesTable;

    cudaError_t                   stat;
    if (pmeGPU->archSpecific->useTextureObjects)
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
    else
    {
        cudaChannelFormatDesc cd_fsh = cudaCreateChannelDesc<float>();
        stat = cudaBindTexture(NULL, &fractShiftsTableTextureRef, fractShiftsArray, &cd_fsh, newFractShiftsSize * sizeof(float));
        CU_RET_ERR(stat, "cudaBindTexture on fractShiftsTableTextureRef failed");

        cudaChannelFormatDesc cd_nn = cudaCreateChannelDesc<int>();
        stat = cudaBindTexture(NULL, &gridlineIndicesTableTextureRef, gridlineIndicesArray, &cd_nn, newFractShiftsSize * sizeof(int));
        CU_RET_ERR(stat, "cudaBindTexture on gridlineIndicesTableTextureRef failed");
    }
#else // !PME_USE_TEXTURES
    GMX_UNUSED_VALUE(pmeGPU);
#endif
}

void pme_gpu_free_fract_shifts_textures(const pme_gpu_t *pmeGPU)
{
#if PME_USE_TEXTURES
    pme_gpu_cuda_kernel_params_t *kernelParamsPtr = pmeGPU->kernelParams.get();
    cudaError_t                   stat;
    if (pmeGPU->archSpecific->useTextureObjects)
    {
        stat = cudaDestroyTextureObject(kernelParamsPtr->fractShiftsTableTexture);
        CU_RET_ERR(stat, "cudaDestroyTextureObject on fractShiftsTableTexture failed");
        stat = cudaDestroyTextureObject(kernelParamsPtr->gridlineIndicesTableTexture);
        CU_RET_ERR(stat, "cudaDestroyTextureObject on gridlineIndicesTableTexture failed");
    }
    else
    {
        stat = cudaUnbindTexture(fractShiftsTableTextureRef);
        CU_RET_ERR(stat, "cudaUnbindTexture on fractShiftsTableTextureRef failed");
        stat = cudaUnbindTexture(gridlineIndicesTableTextureRef);
        CU_RET_ERR(stat, "cudaUnbindTexture on gridlineIndicesTableTextureRef failed");
    }
#else // !PME_USE_TEXTURES
    GMX_UNUSED_VALUE(pmeGPU);
#endif
}

void pme_gpu_spread(const pme_gpu_t *pmeGpu,
                    int gmx_unused   gridIndex,
                    real            *h_grid,
                    bool             computeSplines,
                    bool             spreadCharges)
{
    GMX_RELEASE_ASSERT(computeSplines || spreadCharges, "PME spread kernel has invalid input (nothing to do)");

    cudaStream_t                        s               = pmeGpu->archSpecific->pmeStream;
    const pme_gpu_cuda_kernel_params_t *kernelParamsPtr = pmeGpu->kernelParams.get();

    const int    order   = pmeGpu->common->pme_order;

    const size_t gridSize = pmeGpu->archSpecific->gridSize * sizeof(float);

    const int    atomsPerBlock           = PME_SPREADGATHER_ATOMS_PER_BLOCK;
    //const int    splineParticlesPerBlock = PME_SPREADGATHER_ATOMS_PER_BLOCK;
    //blockSize / DIM; - can be easily changed, just have to pass spread theta stride to the spline kernel!
    // duplicated below!

    GMX_RELEASE_ASSERT(kernelParamsPtr->atoms.nAtoms > 0, "No atom data in PME GPU spread");
    dim3 nBlocksSpread(pmeGpu->nAtomsPadded / atomsPerBlock);
    //dim3 nBlocksSpline((kernelParamsPtr->atoms.nAtoms + splineParticlesPerBlock - 1) / splineParticlesPerBlock); //???
    dim3 dimBlockSpread(order, order, atomsPerBlock); // used for spline_and_spread / separate spread
    //dim3 dimBlockSpline(splineParticlesPerBlock, DIM);                                                           // used for separate spline only

    // These should later check for PME decomposition
    const bool wrapX = true;
    const bool wrapY = true;
    GMX_UNUSED_VALUE(wrapX);
    GMX_UNUSED_VALUE(wrapY);
    switch (order)
    {
        case 4:
        {
            pme_gpu_start_timing(pmeGpu, gtPME_SPLINEANDSPREAD);
            if (computeSplines)
            {
                if (spreadCharges)
                {
                    pme_spline_and_spread_kernel<4, true, true, wrapX, wrapY> <<< nBlocksSpread, dimBlockSpread, 0, s>>> (*kernelParamsPtr);
                }
                else
                {
                    pme_spline_and_spread_kernel<4, true, false, wrapX, wrapY> <<< nBlocksSpread, dimBlockSpread, 0, s>>> (*kernelParamsPtr);

                }
            }
            else
            {
                pme_spline_and_spread_kernel<4, false, true, wrapX, wrapY> <<< nBlocksSpread, dimBlockSpread, 0, s>>> (*kernelParamsPtr);
            }
            CU_LAUNCH_ERR("pme_spline_and_spread_kernel");
            pme_gpu_stop_timing(pmeGpu, gtPME_SPLINEANDSPREAD);
        }
        break;

        default:
            gmx_fatal(FARGS, "the code for pme_order != 4 was not tested!");
    }

    const bool copyBackAtomData = computeSplines && (pme_gpu_is_testing(pmeGpu) || !pme_gpu_performs_gather(pmeGpu));
    if (copyBackAtomData)
    {
        const int    alignment     = pme_gpu_get_atom_spline_data_alignment(pmeGpu);
        const size_t nAtomsPadded  = ((pmeGpu->nAtomsAlloc + alignment - 1) / alignment) * alignment;
        const size_t splinesSize   = DIM * nAtomsPadded * order * sizeof(float);
        cu_copy_D2H_async(pmeGpu->staging.h_dtheta, kernelParamsPtr->atoms.d_dtheta, splinesSize, s);
        cu_copy_D2H_async(pmeGpu->staging.h_theta, kernelParamsPtr->atoms.d_theta, splinesSize, s);
        cu_copy_D2H_async(pmeGpu->staging.h_gridlineIndices, kernelParamsPtr->atoms.d_gridlineIndices, kernelParamsPtr->atoms.nAtoms * DIM * sizeof(int), s);
        cudaError_t stat = cudaEventRecord(pmeGpu->archSpecific->syncSplineAtomDataD2H, s);
        CU_RET_ERR(stat, "PME spread atom data sync fail");
    }
    const bool copyBackGrid     = spreadCharges && (pme_gpu_is_testing(pmeGpu) || !pme_gpu_performs_FFT(pmeGpu));
    if (copyBackGrid)
    {
        cu_copy_D2H_async(h_grid, kernelParamsPtr->grid.d_realGrid, gridSize, s);
        cudaError_t stat = cudaEventRecord(pmeGpu->archSpecific->syncSpreadGridD2H, s);
        CU_RET_ERR(stat, "PME spread grid sync fail");
    }
}
