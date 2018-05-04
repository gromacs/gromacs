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

#include "../../ewald/pme-ocl-types-kernel.clh"

#ifndef COMPILE_GATHER_HELPERS_ONCE //FIXME
#define COMPILE_GATHER_HELPERS_ONCE

/*! \brief
 * An inline CUDA function: unroll the dynamic index accesses to the constant grid sizes to avoid local memory operations.
 */
DEVICE_INLINE float read_grid_size(const float *realGridSizeFP,
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
#if CAN_USE_TEMPLATES
template <
    const int order,
    const int atomDataSize,
    const int blockSize
    >
#endif
DEVICE_INLINE void reduce_atom_forces(SHARED float3 * __restrict__ sm_forces,
                                                   const int             atomIndexLocal,
                                                   const int             splineIndex,
                                                   const int             lineIndex,
                                                   const float          *realGridSizeFP,
                                                   float                fx,
                                                   float                fy,
                                                   float                fz,
                                                   //FIXME these 3 guys were references in CUDA, which is not alright in C99, only pointers would work.
						  SHARED float * __restrict__ sm_forceReduction,
						   SHARED float ** __restrict__ sm_forceTemp
						   )
                                                   
{
#if !CAN_USE_TEMPLATES
        //FIXME define or whatever
        #define atomDataSize PME_SPREADGATHER_THREADS_PER_ATOM
        /* Number of data components and threads for a single atom */
#endif


#if (GMX_PTX_ARCH >= 300) && defined(FIXME)
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
#endif
    {
        //const int    blockSize      = atomsPerBlock * atomDataSize;
        #define blockSize (atomsPerBlock * atomDataSize)

        #define smemPerDim warp_size
        #define smemReserved  ((DIM - 1) * smemPerDim)

        #define totalSharedMemory (smemReserved + blockSize)

        //FIXME

        // We use blockSize shared memory elements to read fx, or fy, or fz, and then reduce them to fit into smemPerDim elements
        // which are stored separately (first 2 dimensions only)
        /*
        const int         smemPerDim   = warp_size;
        const int         smemReserved = (DIM - 1) * smemPerDim;
        */

        /* //FIXME moved outside - can't be defined here on Intel
        SHARED float  sm_forceReduction[totalSharedMemory];
        SHARED float *sm_forceTemp[DIM];
        */

        const int         numWarps  = blockSize / smemPerDim;
        const int         minStride = max(1, atomDataSize / numWarps); // order 4: 128 threads => 4, 256 threads => 2, etc

#pragma unroll
//FIXME __attribute__((opencl_unroll_hint))
        for (int dimIndex = 0; dimIndex < DIM; dimIndex++)
        {
            int elementIndex = smemReserved + lineIndex;
            // Store input force contributions
            sm_forceReduction[elementIndex] = (dimIndex == XX) ? fx : (dimIndex == YY) ? fy : fz;
            sharedMemoryBarrier(); //FIXME why is this needed? unroll?

            // Reduce to fit into smemPerDim (warp size)
#pragma unroll
            for (int redStride = atomDataSize / 2; redStride > minStride; redStride >>= 1)
            {
                if (splineIndex < redStride)
                {
                    sm_forceReduction[elementIndex] += sm_forceReduction[elementIndex + redStride];
                }
            }
            sharedMemoryBarrier();
            // Last iteration - packing everything to be nearby, storing convenience pointer
            sm_forceTemp[dimIndex] = sm_forceReduction + dimIndex * smemPerDim;
            int redStride = minStride;
            if (splineIndex < redStride)
            {
                const int packedIndex = atomIndexLocal * redStride + splineIndex;
                sm_forceTemp[dimIndex][packedIndex] = sm_forceReduction[elementIndex] + sm_forceReduction[elementIndex + redStride];
            }
        }

        sharedMemoryBarrier();

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
                *((SHARED float *)(&sm_forces[atomIndex]) + dimIndex) = (sm_forceTemp[dimIndex][sourceIndex] + sm_forceTemp[dimIndex][sourceIndex + 1]) * n;
            }
        }
    }
}

#endif //COMPILE_GATHER_HELPERS_ONCE

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
KERNEL_FUNC void CUSTOMIZED_KERNEL_NAME(pme_gather_kernel)(const struct PmeGpuCudaKernelParams       kernelParams
#if !CAN_USE_BUFFERS_IN_STRUCTS
                                   ,
                                   GLOBAL const float * __restrict__  gm_coefficients,
                                   GLOBAL const float * __restrict__  gm_grid,
                                   GLOBAL const float * __restrict__  gm_theta,
                                   GLOBAL const float * __restrict__  gm_dtheta,
                                   GLOBAL const int * __restrict__    gm_gridlineIndices,
                                   GLOBAL float * __restrict__        gm_forces
#endif
)
{
#if CAN_USE_BUFFERS_IN_STRUCTS
    /* Global memory pointers */
    GLOBAL const float * __restrict__  gm_coefficients     = kernelParams.atoms.d_coefficients;
    GLOBAL const float * __restrict__  gm_grid             = kernelParams.grid.d_realGrid;
    GLOBAL const float * __restrict__  gm_theta            = kernelParams.atoms.d_theta;
    GLOBAL const float * __restrict__  gm_dtheta           = kernelParams.atoms.d_dtheta;
    GLOBAL const int * __restrict__    gm_gridlineIndices  = kernelParams.atoms.d_gridlineIndices;
    GLOBAL float * __restrict__        gm_forces           = kernelParams.atoms.d_forces;
#endif

    /* Some sizes */
#if CAN_USE_TEMPLATES
    const int    atomsPerBlock  = (c_gatherMaxThreadsPerBlock / PME_SPREADGATHER_THREADS_PER_ATOM);
#endif
    //const int    atomDataSize   = PME_SPREADGATHER_THREADS_PER_ATOM; /* Number of data components and threads for a single atom */
    //const int    blockSize      = atomsPerBlock * atomDataSize; //FIXME move into reduction

    /* These are the atom indices - for the shared and global memory */
    const int         atomIndexLocal    = getThreadLocalIndex(ZZ);
    const int         atomIndexOffset   = getBlockIndex(XX) * atomsPerBlock;
    const int         atomIndexGlobal   = atomIndexOffset + atomIndexLocal;

#if USE_C99_ONLY
#define gridlineIndicesSize (atomsPerBlock * DIM)
#define splineParamsSize (atomsPerBlock * DIM * order)
#else
    constexpr int splineParamsSize             = atomsPerBlock * DIM * order;
    constexpr int gridlineIndicesSize          = atomsPerBlock * DIM;
#endif
    SHARED int    sm_gridlineIndices[gridlineIndicesSize];
    SHARED float2 sm_splineParams[splineParamsSize]; /* Theta/dtheta pairs  as .x/.y */

    /* Spline Y/Z coordinates */
    const int ithy = getThreadLocalIndex(YY);
    const int ithz = getThreadLocalIndex(XX);

    const int threadLocalId = getThreadLocalIndex3d();

    /* These are the spline contribution indices in shared memory */
    const int splineIndex = getThreadLocalIndex2d();  /* Relative to the current particle , 0..15 for order 4 */
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
    sharedMemoryBarrier();

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

        const int    thetaOffsetBase = PME_SPLINE_THETA_STRIDE * order * warpIndex * DIM * PME_SPREADGATHER_ATOMS_PER_WARP + particleWarpIndex;
        const int    orderStride     = PME_SPLINE_THETA_STRIDE * DIM * PME_SPREADGATHER_ATOMS_PER_WARP;
        const int    dimStride       = PME_SPLINE_THETA_STRIDE * PME_SPREADGATHER_ATOMS_PER_WARP;

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
    SHARED float3 sm_forces[atomsPerBlock]; //FIXME float3?


    //FIXME moved from reduction function 
    #define blockSize (atomsPerBlock * atomDataSize)
    #define smemPerDim warp_size
    #define smemReserved  ((DIM - 1) * smemPerDim)
    #define totalSharedMemory (smemReserved + blockSize)
    SHARED float  sm_forceReduction[totalSharedMemory];
    SHARED float *sm_forceTemp[DIM];

    reduce_atom_forces TEMPLATE_PARAMETERS3(order, atomDataSize, blockSize) (sm_forces,
                                                       atomIndexLocal, splineIndex, lineIndex,
                                                       kernelParams.grid.realGridSizeFP,
                                                       fx, fy, fz,
						       sm_forceReduction,
						       sm_forceTemp);
    sharedMemoryBarrier();

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

    //TODO this is to make "warp_size" of 32 wrk on intel
    sharedMemoryBarrier();
    //gmx_syncwarp();
    assert(atomsPerBlock <= warp_size);

    /* Writing or adding the final forces component-wise, single warp */ //FIXME is this correct?
    const int blockForcesSize = atomsPerBlock * DIM;
    const int numIter         = (blockForcesSize + warp_size - 1) / warp_size;
    const int iterThreads     = blockForcesSize / numIter;
    if (threadLocalId < iterThreads)
    {
#pragma unroll
        for (int i = 0; i < numIter; i++)
        {
            int         outputIndexLocal  = i * iterThreads + threadLocalId;
            int         outputIndexGlobal = getBlockIndex(XX) * blockForcesSize + outputIndexLocal;
            const int   globalOutputCheck = pme_gpu_check_atom_data_index(outputIndexGlobal, kernelParams.atoms.nAtoms * DIM);
            if (globalOutputCheck)
            {
                const float outputForceComponent = ((SHARED float *)sm_forces)[outputIndexLocal * 4 / 3]; //FIXME
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

