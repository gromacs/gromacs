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
#include "gromacs/utility/fatalerror.h"

#include "pme.cuh"
#include "pme-grid.h"
#include "pme-internal.h"
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
    /* Spline parameter storage, shared for PME_GPU_PARALLEL_SPLINE==1 to not overuse the local memory */
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
            assert(tableIndex >= 0);

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

        const int       chargeCheck = pme_gpu_check_atom_charge(sm_coefficients[localIndexCalc]);
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

                sm_theta[thetaIndex]             = data[k * dataSize + dataOffset];
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
        const int ithy      = threadIdx.y;
        const int ithz      = threadIdx.x; //?
        const int ixBase    = sm_gridlineIndices[localIndex * DIM + XX] - offx;
        const int iy        = sm_gridlineIndices[localIndex * DIM + YY] - offy + ithy;
        const int iyWrapped = (wrapY & (iy >= ny)) ? (iy - ny) : iy;
        const int iz        = sm_gridlineIndices[localIndex * DIM + ZZ] - offz + ithz;
        const int izWrapped = (iz >= nz) ? (iz - nz) : iz;

        // copy
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

template <
    const int atomsPerBlock
    >
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

template <
    const int atomsPerBlock
    >
__device__ __forceinline__ void stage_coordinates(const int                          threadLocalId,
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
    const gmx_bool bCalcSplines,     // first part
    const gmx_bool bSpread,          // second part
    const bool wrapX,
    const bool wrapY
    >
__launch_bounds__(PME_SPREADGATHER_THREADS_PER_BLOCK, PME_MIN_BLOCKS_PER_MP)
__global__ void pme_spline_and_spread_kernel(const pme_gpu_cuda_kernel_params_t kernelParams)
{
    const int                     atomsPerBlock     = PME_SPREADGATHER_ATOMS_PER_BLOCK;
    // gridline indices, ivec
    __shared__ int                sm_gridlineIndices[atomsPerBlock * DIM]; // TODO clean tehse up?
    // charges
    __shared__ float              sm_coefficients[atomsPerBlock];
    // spline parameters
    __shared__ float              sm_theta[atomsPerBlock * DIM * order];

    const int                     threadLocalId = (threadIdx.z * (blockDim.x * blockDim.y))
        + (threadIdx.y * blockDim.x)
        + threadIdx.x;

    const int globalParticleIndexBase = blockIdx.x * atomsPerBlock;

    const int warpIndex         = threadLocalId / warp_size;
    const int threadWarpIndex   = threadLocalId % warp_size;
    const int particleWarpIndex = threadWarpIndex % PME_SPREADGATHER_ATOMS_PER_WARP;
    const int localCalcIndex    = warpIndex * PME_SPREADGATHER_ATOMS_PER_WARP + particleWarpIndex;
    const int globalCalcIndex   = globalParticleIndexBase + localCalcIndex;
    const int orderCalcIndex    = threadWarpIndex / (PME_SPREADGATHER_ATOMS_PER_WARP * DIM); // should be checked against order
    const int dimCalcIndex      = (threadWarpIndex - orderCalcIndex * (PME_SPREADGATHER_ATOMS_PER_WARP * DIM)) / PME_SPREADGATHER_ATOMS_PER_WARP;

    if (bCalcSplines)
    {
        // coordinates
        __shared__ float sm_coordinates[DIM * atomsPerBlock];

        stage_charges<atomsPerBlock>(threadLocalId, sm_coefficients, kernelParams);
        stage_coordinates<atomsPerBlock>(threadLocalId, sm_coordinates, kernelParams);
        __syncthreads();
        calculate_splines<order>((const float3 *)sm_coordinates, sm_coefficients,
                                 sm_theta, sm_gridlineIndices,
                                 kernelParams,
                                 globalCalcIndex,
                                 localCalcIndex,
                                 globalParticleIndexBase,
                                 dimCalcIndex,
                                 orderCalcIndex);
    }
    else if (bSpread) // staging for spread
    {
        //yupinov - unmaintained
        /*
           if ((globalIndexCalc < n) && (dimIndex < DIM) && (localIndexCalc < atomsPerBlock))
           {
           gridlineIndices[localIndexCalc * DIM + dimIndex] = gridlineIndicesGlobal[globalIndexCalc * DIM + dimIndex];

           const int thetaOffsetBase = localIndexCalc * DIM + dimIndex;
           const int thetaGlobalOffsetBase = globalIndexBase * DIM * order;
           #pragma unroll
           for (int k = 0; k < order; k++)
           {
            const int thetaIndex = thetaOffsetBase + k * thetaStride;
            theta[thetaIndex] = thetaGlobal[thetaGlobalOffsetBase + thetaIndex];
           }
           }
         */
        stage_charges<atomsPerBlock>(threadLocalId, sm_coefficients, kernelParams);
        __syncthreads();
    }

    // SPREAD
    if (bSpread)
    {
        const int localSpreadIndex  = threadIdx.z;
        const int globalSpreadIndex = globalParticleIndexBase + localSpreadIndex;
        spread_charges<order, wrapX, wrapY>(sm_coefficients, kernelParams, globalSpreadIndex, localSpreadIndex,
                                            sm_gridlineIndices, sm_theta);
    }
}


// pme_spline_and_spread split into pme_spline and pme_spread - as an experiment

template <
    const int order>
__launch_bounds__(PME_SPREADGATHER_THREADS_PER_BLOCK, PME_MIN_BLOCKS_PER_MP)
__global__ void pme_spline_kernel(const pme_gpu_cuda_kernel_params_t kernelParams)
{
    const int                    atomsPerBlock     = PME_SPREADGATHER_ATOMS_PER_BLOCK;
    // gridline indices
    __shared__ int               sm_gridlineIndices[DIM * atomsPerBlock];
    // charges
    __shared__ float             sm_coefficients[atomsPerBlock];
    // coordinates
    __shared__ float             sm_coordinates[DIM * atomsPerBlock];
    // spline parameters
    __shared__ float             sm_theta[DIM * atomsPerBlock * order];

    const int                    globalIndexBase = blockIdx.x * atomsPerBlock;

    const int                    threadLocalId = (threadIdx.z * (blockDim.x * blockDim.y))
        + (threadIdx.y * blockDim.x)
        + threadIdx.x;

    const int localIndexCalc  = threadIdx.x;
    const int orderIndex      = threadIdx.x; //yupinov - this is broken(?)
    const int dimIndex        = threadIdx.y;
    const int globalIndexCalc = globalIndexBase + localIndexCalc;

    stage_charges<atomsPerBlock>(threadLocalId, sm_coefficients, kernelParams);
    stage_coordinates<atomsPerBlock>(threadLocalId, sm_coordinates, kernelParams);
    __syncthreads();
    // TODO: clean up type casts

    calculate_splines<order>((const float3 *)sm_coordinates, sm_coefficients,
                             sm_theta, sm_gridlineIndices,
                             kernelParams,
                             globalIndexCalc,
                             localIndexCalc,
                             globalIndexBase,
                             dimIndex,
                             orderIndex);
}


template
<const int order, const bool wrapX, const bool wrapY>
__launch_bounds__(PME_SPREADGATHER_THREADS_PER_BLOCK, PME_MIN_BLOCKS_PER_MP)
__global__ void pme_spread_kernel(const pme_gpu_cuda_kernel_params_t kernelParams)
{
    constexpr size_t           atomsPerBlock     = PME_SPREADGATHER_ATOMS_PER_BLOCK;
    /* Global memory pointers */
    const int * __restrict__   gm_gridlineIndices = kernelParams.atoms.d_gridlineIndices;
    float * __restrict__       gm_theta           = kernelParams.atoms.d_theta;

    __shared__ int             sm_gridlineIndices[atomsPerBlock * DIM];
    __shared__ float           sm_coefficients[atomsPerBlock];
    __shared__ float           sm_theta[atomsPerBlock * DIM * order];

    const int                  localIndex              = threadIdx.x;
    const int                  globalParticleIndexBase = blockIdx.x * atomsPerBlock;
    const int                  globalIndex             = globalParticleIndexBase + localIndex;


    //yupinov - staging
    const int threadLocalId = (threadIdx.z * (blockDim.x * blockDim.y))
        + (threadIdx.y * blockDim.x)
        + threadIdx.x;


    stage_charges<atomsPerBlock>(threadLocalId, sm_coefficients, kernelParams);
    __syncthreads();

    const int localIndexCalc  = threadLocalId / DIM;
    const int dimIndex        = threadLocalId - localIndexCalc * DIM;
    const int globalIndexCalc = globalParticleIndexBase + localIndexCalc;
    const int globalCheck     = pme_gpu_check_atom_data_index(globalIndexCalc, kernelParams.atoms.nAtoms);

    if ((dimIndex < DIM) && (localIndexCalc < atomsPerBlock) && globalCheck)
    {
        sm_gridlineIndices[localIndexCalc * DIM + dimIndex] = gm_gridlineIndices[globalIndexCalc * DIM + dimIndex];

        //unmaintained...
        const int thetaOffsetBase       = localIndexCalc * DIM + dimIndex;
        const int thetaGlobalOffsetBase = globalParticleIndexBase * DIM * order;
#pragma unroll
        for (int k = 0; k < order; k++)
        {
            const int thetaIndex = thetaOffsetBase + k * PME_SPLINE_ORDER_STRIDE;
            sm_theta[thetaIndex] = gm_theta[thetaGlobalOffsetBase + thetaIndex];
        }
    }
    __syncthreads();

    // SPREAD
    spread_charges<order, wrapX, wrapY>(sm_coefficients, kernelParams, globalIndex, localIndex,
                                        sm_gridlineIndices, sm_theta);
}

void pme_gpu_make_fract_shifts_textures(pme_gpu_t *pmeGPU)
{
#if PME_USE_TEXTURES
    const int                     nx                  = pmeGPU->common->nk[XX];
    const int                     ny                  = pmeGPU->common->nk[YY];
    const int                     nz                  = pmeGPU->common->nk[ZZ];
    const int                     cellCount           = c_pmeNeighborUnitcellCount;
    const int                     newFractShiftsSize  = cellCount * (nx + ny + nz);

    pme_gpu_cuda_kernel_params_t *kernelParamsPtr = pmeGPU->kernelParams.get();

    float                        *fshArray = kernelParamsPtr->grid.d_fractShiftsTable;
    int                          *nnArray  = kernelParamsPtr->grid.d_gridlineIndicesTable;

    cudaError_t                   stat;
#if PME_USE_TEXOBJ
    //if (use_texobj(dev_info))
    // TODO: should check device info here for CC >= 3.0
    {
        cudaResourceDesc rd;
        cudaTextureDesc  td;

        memset(&rd, 0, sizeof(rd));
        rd.resType                  = cudaResourceTypeLinear;
        rd.res.linear.devPtr        = fshArray;
        rd.res.linear.desc.f        = cudaChannelFormatKindFloat;
        rd.res.linear.desc.x        = 32;
        rd.res.linear.sizeInBytes   = newFractShiftsSize * sizeof(float);
        memset(&td, 0, sizeof(td));
        td.readMode                 = cudaReadModeElementType;
        stat = cudaCreateTextureObject(&kernelParamsPtr->fractShiftsTableTexture, &rd, &td, NULL);
        CU_RET_ERR(stat, "cudaCreateTextureObject on fsh_d failed");


        memset(&rd, 0, sizeof(rd));
        rd.resType                  = cudaResourceTypeLinear;
        rd.res.linear.devPtr        = nnArray;
        rd.res.linear.desc.f        = cudaChannelFormatKindSigned;
        rd.res.linear.desc.x        = 32;
        rd.res.linear.sizeInBytes   = newFractShiftsSize * sizeof(int);
        memset(&td, 0, sizeof(td));
        td.readMode                 = cudaReadModeElementType;
        stat = cudaCreateTextureObject(&kernelParamsPtr->gridlineIndicesTableTexture, &rd, &td, NULL);
        CU_RET_ERR(stat, "cudaCreateTextureObject on nn_d failed");
    }
    //else
#else
    {
        cudaChannelFormatDesc cd_fsh = cudaCreateChannelDesc<float>();
        stat = cudaBindTexture(NULL, &fractShiftsTableTextureRef, fshArray, &cd_fsh, newFractsShiftSize * sizeof(float));
        CU_RET_ERR(stat, "cudaBindTexture on fsh failed");

        cudaChannelFormatDesc cd_nn = cudaCreateChannelDesc<int>();
        stat = cudaBindTexture(NULL, &gridlineIndicesTableTextureRef, nnArray, &cd_nn, newFractShiftsSize * sizeof(int));
        CU_RET_ERR(stat, "cudaBindTexture on nn failed");
    }
#endif
#else
    GMX_UNUSED_VALUE(pmeGPU);
#endif
}

void pme_gpu_free_fract_shifts_textures(const pme_gpu_t *pmeGPU)
{
    /* TODO: unbind textures here! */
    GMX_UNUSED_VALUE(pmeGPU);
}



void pme_gpu_spread(const gmx_pme_t *pme, pme_atomcomm_t gmx_unused *atc,
                    const int gmx_unused grid_index,
                    pmegrid_t *pmegrid,
                    const gmx_bool bCalcSplines,
                    const gmx_bool bSpread)
{
    const gmx_bool bSeparateKernels = FALSE;  // significantly slower if true
    if (!bCalcSplines && !bSpread)
    {
        gmx_fatal(FARGS, "No splining or spreading to be done?"); //yupinov use of gmx_fatal
    }
    const pme_gpu_t                    *pmeGPU = pme->gpu;

    cudaStream_t                        s               = pmeGPU->archSpecific->pmeStream;
    const pme_gpu_cuda_kernel_params_t *kernelParamsPtr = pmeGPU->kernelParams.get();

    const int    order   = pmeGPU->common->pme_order;

    const size_t gridSize = pmeGPU->archSpecific->gridSize * sizeof(float);

    const int    atomsPerBlock           = PME_SPREADGATHER_ATOMS_PER_BLOCK;
    const int    splineParticlesPerBlock = PME_SPREADGATHER_ATOMS_PER_BLOCK;
    //blockSize / DIM; - can be easily changed, just have to pass spread theta stride to the spline kernel!
    // duplicated below!

    GMX_RELEASE_ASSERT(kernelParamsPtr->atoms.nAtoms > 0, "No atom data in PME GPU spread");
    dim3 nBlocksSpread(pmeGPU->nAtomsPadded / atomsPerBlock);
    dim3 nBlocksSpline((kernelParamsPtr->atoms.nAtoms + splineParticlesPerBlock - 1) / splineParticlesPerBlock); //???
    dim3 dimBlockSpread(order, order, atomsPerBlock);                                                            // used for spline_and_spread / separate spread
    dim3 dimBlockSpline(splineParticlesPerBlock, DIM);                                                           // used for separate spline only
    const bool wrapX = true;                                                                                     //should check DD variables
    const bool wrapY = true;                                                                                     //should check DD variables
    switch (order)
    {
        case 4:
            if (bSeparateKernels)
            {
                if (bCalcSplines)
                {
                    pme_gpu_start_timing(pmeGPU, gtPME_SPLINE);
                    pme_spline_kernel<4> <<< nBlocksSpline, dimBlockSpline, 0, s>>> (*kernelParamsPtr);
                    CU_LAUNCH_ERR("pme_spline_kernel");
                    pme_gpu_stop_timing(pmeGPU, gtPME_SPLINE);
                }
                if (bSpread)
                {
                    pme_gpu_start_timing(pmeGPU, gtPME_SPREAD);
                    pme_spread_kernel<4, wrapX, wrapY> <<< nBlocksSpread, dimBlockSpread, 0, s>>> (*kernelParamsPtr);
                    CU_LAUNCH_ERR("pme_spread_kernel");
                    pme_gpu_stop_timing(pmeGPU, gtPME_SPREAD);
                }
            }
            else // a single monster kernel here
            {
                pme_gpu_start_timing(pmeGPU, gtPME_SPLINEANDSPREAD);
                if (bCalcSplines)
                {
                    if (bSpread)
                    {
                        pme_spline_and_spread_kernel<4, TRUE, TRUE, wrapX, wrapY> <<< nBlocksSpread, dimBlockSpread, 0, s>>> (*kernelParamsPtr);
                    }
                    else
                    {
                        gmx_fatal(FARGS, "the code for bSpread==false was not tested!");
                    }
                }
                else
                {
                    gmx_fatal(FARGS, "the code for bCalcSplines==false was not tested!");
                }
                CU_LAUNCH_ERR("pme_spline_and_spread_kernel");
                pme_gpu_stop_timing(pmeGPU, gtPME_SPLINEANDSPREAD);
            }
            break;

        default:
            gmx_fatal(FARGS, "the code for pme_order != 4 was not tested!");
    }

    const bool copyBackGrid     = pme_gpu_is_testing(pmeGPU) || !pme_gpu_performs_FFT(pmeGPU) && bSpread;
    const bool copyBackAtomData = pme_gpu_is_testing(pmeGPU) || !pme_gpu_performs_gather(pmeGPU);
    if (copyBackAtomData)
    {
        // FIXME: spline parameters layout is not the same on GPU => this would would fail with CPU gather.
        // Also no accounting for PME communication (bGPUSingle check?)
        const int    alignment     = pme_gpu_get_atom_spline_data_alignment(pmeGPU);
        const size_t nAtomsPadded  = ((pmeGPU->nAtomsAlloc + alignment - 1) / alignment) * alignment;
        const size_t splinesSize   = DIM * nAtomsPadded * order * sizeof(float);
        // TODO fractCoordinates as well?
        cu_copy_D2H_async(pmeGPU->staging.h_dtheta, kernelParamsPtr->atoms.d_dtheta, splinesSize, s);
        cu_copy_D2H_async(pmeGPU->staging.h_theta, kernelParamsPtr->atoms.d_theta, splinesSize, s); //hide it in a func???
        cu_copy_D2H_async(atc->idx, kernelParamsPtr->atoms.d_gridlineIndices, kernelParamsPtr->atoms.nAtoms * DIM * sizeof(int), s);
        cudaError_t stat = cudaEventRecord(pmeGPU->archSpecific->syncSpreadAtomDataD2H, s);
        CU_RET_ERR(stat, "PME spread atom data sync fail");
    }
    if (copyBackGrid)
    {
        cu_copy_D2H_async(pmegrid->grid, kernelParamsPtr->grid.d_realGrid, gridSize, s);
        cudaError_t stat = cudaEventRecord(pmeGPU->archSpecific->syncSpreadGridD2H, s);
        CU_RET_ERR(stat, "PME spread grid sync fail");
    }

}
