/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
 *  \brief Implements PME force gathering kernel in OpenCL.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 */

#include "pme-gpu-types.h"

#ifndef COMPILE_GATHER_HELPERS_ONCE
#define COMPILE_GATHER_HELPERS_ONCE

/*! \brief
 * Unrolls the dynamic index accesses to the constant grid sizes to avoid local memory operations.
 */
inline float read_grid_size(const float *realGridSizeFP,
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

/*! \brief Reduce the partial force contributions.
 *
 * \tparam[in] order              The PME order (must be 4).
 * \tparam[in] atomDataSize       The number of partial force contributions for each atom (currently order^2 == 16)
 * \tparam[in] blockSize          The kernel block size
 * \param[out] sm_forces          Local memory array with the output forces (rvec).
 * \param[in]  atomIndexLocal     Local atom index
 * \param[in]  splineIndex        Spline index
 * \param[in]  lineIndex          Line index (same as threadLocalId)
 * \param[in]  realGridSizeFP     Local grid size constant
 * \param[in]  fx                 Input force partial component X
 * \param[in]  fy                 Input force partial component Y
 * \param[in]  fz                 Input force partial component Z
 */
#if CAN_USE_TEMPLATES
template <
    const int order,
    const int atomDataSize,
    const int blockSize
    >
#endif
inline void reduce_atom_forces(__local float * __restrict__  sm_forces,
                               const int                     atomIndexLocal,
                               const int                     splineIndex,
                               const int                     lineIndex,
                               const float                  *realGridSizeFP,
                               float                         fx,
                               float                         fy,
                               float                         fz,
                               __local float * __restrict__  sm_forceReduction,
                               __local float ** __restrict__ sm_forceTemp
                               )

{
    // TODO: implement AMD intrinsics reduction, like with shuffles in CUDA version. #2514

    /* Number of data components and threads for a single atom */
#define atomDataSize PME_SPREADGATHER_THREADS_PER_ATOM
    // We use blockSize local memory elements to read fx, or fy, or fz, and then reduce them to fit into smemPerDim elements
    // All those guys are defines and not consts, because they go into the local memory array size.
#define blockSize (atomsPerBlock * atomDataSize)
#define smemPerDim warp_size
#define smemReserved  ((DIM - 1) * smemPerDim)

    const int numWarps  = blockSize / smemPerDim;
    const int minStride = max(1, atomDataSize / numWarps);

#pragma unroll DIM
    for (int dimIndex = 0; dimIndex < DIM; dimIndex++)
    {
        int elementIndex = smemReserved + lineIndex;
        // Store input force contributions
        sm_forceReduction[elementIndex] = (dimIndex == XX) ? fx : (dimIndex == YY) ? fy : fz;
        /* This barrier was not needed in CUDA. Different OpenCL compilers might have different ideas
         * about #pragma unroll, though. OpenCL 2 has _attribute__((opencl_unroll_hint)).
         * #2519
         */
        localMemoryBarrier();

        // Reduce to fit into smemPerDim (warp size)
#pragma unroll
        for (int redStride = atomDataSize / 2; redStride > minStride; redStride >>= 1)
        {
            if (splineIndex < redStride)
            {
                sm_forceReduction[elementIndex] += sm_forceReduction[elementIndex + redStride];
            }
        }
        localMemoryBarrier();
        // Last iteration - packing everything to be nearby, storing convenience pointer
        sm_forceTemp[dimIndex] = sm_forceReduction + dimIndex * smemPerDim;
        int redStride = minStride;
        if (splineIndex < redStride)
        {
            const int packedIndex = atomIndexLocal * redStride + splineIndex;
            sm_forceTemp[dimIndex][packedIndex] = sm_forceReduction[elementIndex] + sm_forceReduction[elementIndex + redStride];
        }
    }

    localMemoryBarrier();

    assert ((blockSize / warp_size) >= DIM);

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
            sm_forces[atomIndex * DIM + dimIndex] = (sm_forceTemp[dimIndex][sourceIndex] + sm_forceTemp[dimIndex][sourceIndex + 1]) * n;
        }
    }
}

#endif //COMPILE_GATHER_HELPERS_ONCE

/*! \brief
 * An OpenCL kernel which gathers the atom forces from the grid.
 * The grid is assumed to be wrapped in dimension Z.
 *
 * \tparam[in] order                The PME order (must be 4 currently).
 * \tparam[in] overwriteForces      True: the forces are written to the output buffer;
 *                                  False: the forces are added non-atomically to the output buffer (e.g. to the bonded forces).
 * \tparam[in] wrapX                Tells if the grid is wrapped in the X dimension.
 * \tparam[in] wrapY                Tells if the grid is wrapped in the Y dimension.
 * \param[in]  kernelParams         All the PME GPU data.
 */
#if CAN_USE_TEMPLATES
template <
    const int order,
    const bool overwriteForces,
    const bool wrapX,
    const bool wrapY
    >
#endif
//FIXME __launch_bounds__(c_gatherMaxThreadsPerBlock, c_gatherMinBlocksPerMP)
//__attribute__((reqd_work_group_size(order, order, atomsPerBlock))) //and the otehr work size hint maybe?????
__kernel void CUSTOMIZED_KERNEL_NAME(pme_gather_kernel)(const struct PmeOpenCLKernelParams       kernelParams
#if !                                                                                        CAN_USE_BUFFERS_IN_STRUCTS
                                                        ,
                                                        __global const float * __restrict__  gm_coefficients,
                                                        __global const float * __restrict__  gm_grid,
                                                        __global const float * __restrict__  gm_theta,
                                                        __global const float * __restrict__  gm_dtheta,
                                                        __global const int * __restrict__    gm_gridlineIndices,
                                                        __global float * __restrict__        gm_forces
#endif
                                                        )
{
#if CAN_USE_BUFFERS_IN_STRUCTS
    /* Global memory pointers */
    __global const float * __restrict__  gm_coefficients     = kernelParams.atoms.d_coefficients;
    __global const float * __restrict__  gm_grid             = kernelParams.grid.d_realGrid;
    __global const float * __restrict__  gm_theta            = kernelParams.atoms.d_theta;
    __global const float * __restrict__  gm_dtheta           = kernelParams.atoms.d_dtheta;
    __global const int * __restrict__    gm_gridlineIndices  = kernelParams.atoms.d_gridlineIndices;
    __global float * __restrict__        gm_forces           = kernelParams.atoms.d_forces;
#endif

    /* These are the atom indices - for the shared and global memory */
    const int         atomIndexLocal    = getThreadLocalIndex(ZZ);
    const int         atomIndexOffset   = getBlockIndex(XX) * atomsPerBlock;
    const int         atomIndexGlobal   = atomIndexOffset + atomIndexLocal;

#define gridlineIndicesSize (atomsPerBlock * DIM)
#define splineParamsSize (atomsPerBlock * DIM * order)

    __local int    sm_gridlineIndices[gridlineIndicesSize];
    __local float2 sm_splineParams[splineParamsSize]; /* Theta/dtheta pairs  as .x/.y */

    /* Spline Y/Z coordinates */
    const int ithy = getThreadLocalIndex(YY);
    const int ithz = getThreadLocalIndex(XX);

    const int threadLocalId = getThreadLocalIndex3d();

    /* These are the spline contribution indices in shared memory */
    const int splineIndex = getThreadLocalIndex2d(); /* Relative to the current particle , 0..15 for order 4 */
    const int lineIndex   = threadLocalId;           /* And to all the block's particles */

    /* Staging the atom gridline indices, DIM * atomsPerBlock threads */
    const int localGridlineIndicesIndex  = threadLocalId;
    const int globalGridlineIndicesIndex = getBlockIndex(XX) * gridlineIndicesSize + localGridlineIndicesIndex;
    const int globalCheckIndices         = pme_gpu_check_atom_data_index(globalGridlineIndicesIndex, kernelParams.atoms.nAtoms * DIM);
    if ((localGridlineIndicesIndex < gridlineIndicesSize) & globalCheckIndices)
    {
        sm_gridlineIndices[localGridlineIndicesIndex] = gm_gridlineIndices[globalGridlineIndicesIndex];
        assert(sm_gridlineIndices[localGridlineIndicesIndex] >= 0);
    }
    /* Staging the spline parameters, DIM * order * atomsPerBlock threads */
    const int localSplineParamsIndex  = threadLocalId;
    const int globalSplineParamsIndex = getBlockIndex(XX) * splineParamsSize + localSplineParamsIndex;
    const int globalCheckSplineParams = pme_gpu_check_atom_data_index(globalSplineParamsIndex, kernelParams.atoms.nAtoms * DIM * order);
    if ((localSplineParamsIndex < splineParamsSize) && globalCheckSplineParams)
    {
        sm_splineParams[localSplineParamsIndex].x = gm_theta[globalSplineParamsIndex];
        sm_splineParams[localSplineParamsIndex].y = gm_dtheta[globalSplineParamsIndex];
        assert(isfinite(sm_splineParams[localSplineParamsIndex].x));
        assert(isfinite(sm_splineParams[localSplineParamsIndex].y));
    }
    localMemoryBarrier();

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

        const int    particleWarpIndex = atomIndexLocal % PME_SPREADGATHER_ATOMS_PER_WARP;
        const int    warpIndex         = atomIndexLocal / PME_SPREADGATHER_ATOMS_PER_WARP;

        const int    thetaOffsetBase = order * warpIndex * DIM * PME_SPREADGATHER_ATOMS_PER_WARP + particleWarpIndex;
        const int    orderStride     = DIM * PME_SPREADGATHER_ATOMS_PER_WARP;
        const int    dimStride       = PME_SPREADGATHER_ATOMS_PER_WARP;

        const int    thetaOffsetY = thetaOffsetBase + ithy * orderStride + YY * dimStride;
        const float2 tdy          = sm_splineParams[thetaOffsetY];
        const int    thetaOffsetZ = thetaOffsetBase + ithz * orderStride + ZZ * dimStride;
        const float2 tdz          = sm_splineParams[thetaOffsetZ];

        const int    ixBase         = sm_gridlineIndices[atomIndexLocal * DIM + XX];
        int          iy             = sm_gridlineIndices[atomIndexLocal * DIM + YY] + ithy;
        if (wrapY & (iy >= ny))
        {
            iy -= ny;
        }
        int iz  = sm_gridlineIndices[atomIndexLocal * DIM + ZZ] + ithz;
        if (iz >= nz)
        {
            iz -= nz;
        }
        const int constOffset    = iy * pnz + iz;

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
            const int     thetaOffsetX = thetaOffsetBase + ithx * orderStride + XX * dimStride;
            const float2  tdx          = sm_splineParams[thetaOffsetX];
            const float   fxy1         = tdz.x * gridValue;
            const float   fz1          = tdz.y * gridValue;
            fx += tdx.y * tdy.x * fxy1;
            fy += tdx.x * tdy.y * fxy1;
            fz += tdx.x * tdy.x * fz1;
        }
    }

    // Reduction of partial force contributions
    __local float sm_forces[atomsPerBlock * DIM];


    /* Some sizes which are defines and not consts because they go into the array size */
    #define blockSize (atomsPerBlock * atomDataSize)
    #define smemPerDim warp_size
    #define smemReserved  ((DIM - 1) * smemPerDim)
    #define totalSharedMemory (smemReserved + blockSize)
    __local float  sm_forceReduction[totalSharedMemory];
    __local float *sm_forceTemp[DIM];

    reduce_atom_forces(sm_forces,
                       atomIndexLocal, splineIndex, lineIndex,
                       kernelParams.grid.realGridSizeFP,
                       fx, fy, fz,
                       sm_forceReduction,
                       sm_forceTemp);
    localMemoryBarrier();

    /* Calculating the final forces with no component branching, atomsPerBlock threads */
    const int forceIndexLocal  = threadLocalId;
    const int forceIndexGlobal = atomIndexOffset + forceIndexLocal;
    const int calcIndexCheck   = pme_gpu_check_atom_data_index(forceIndexGlobal, kernelParams.atoms.nAtoms);
    if ((forceIndexLocal < atomsPerBlock) & calcIndexCheck)
    {
        const float3  atomForces     = vload3(forceIndexLocal, sm_forces);
        const float   negCoefficient = -gm_coefficients[forceIndexGlobal];
        float3        result;
        result.x      = negCoefficient * kernelParams.current.recipBox[XX][XX] * atomForces.x;
        result.y      = negCoefficient * (kernelParams.current.recipBox[XX][YY] * atomForces.x +
                                          kernelParams.current.recipBox[YY][YY] * atomForces.y);
        result.z      = negCoefficient * (kernelParams.current.recipBox[XX][ZZ] * atomForces.x +
                                          kernelParams.current.recipBox[YY][ZZ] * atomForces.y +
                                          kernelParams.current.recipBox[ZZ][ZZ] * atomForces.z);
        vstore3(result, forceIndexLocal, sm_forces);
    }

    /* This is only here for execution of e.g. 32-sized warps on 16-wide hardware; this was gmx_syncwarp() in CUDA.
     * #2519
     */
    localMemoryBarrier();
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
            const int outputIndexLocal  = i * iterThreads + threadLocalId;
            const int outputIndexGlobal = getBlockIndex(XX) * blockForcesSize + outputIndexLocal;
            const int globalOutputCheck = pme_gpu_check_atom_data_index(outputIndexGlobal, kernelParams.atoms.nAtoms * DIM);
            if (globalOutputCheck)
            {
                const float outputForceComponent = sm_forces[outputIndexLocal];
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
