/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 *  \brief Implements PME force gathering in SYCL.
 *
 *  \author Andrey Alekseenko <al42and@gmail.com>
 */

#include "gmxpre.h"

#include "pme_gather_sycl.h"

#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/gpu_utils/gputraits_sycl.h"
#include "gromacs/gpu_utils/sycl_kernel_utils.h"
#include "gromacs/gpu_utils/syclutils.h"
#include "gromacs/math/functions.h"
#include "gromacs/utility/basedefinitions.h"

#include "pme_gpu_calculate_splines_sycl.h"
#include "pme_gpu_constants.h"
#include "pme_gpu_types_host.h"
#include "pme_grid.h"

/*! \brief
 * Use loads from constant address space indexed by constant offsets rather
 * than dynamic index-based accesses to the grid size data to avoid
 * local memory operations and related large overhead.
 *
 * Drastically reduces register spills on AMD via hipSYCL, and improves performance 10x.
 *
 * \param[in]  realGridSizeFP     Local grid size constant
 * \param[in]  dimIndex           Dimension index (XX, YY, ZZ)
 *
 * \returns The grid size of the specified dimension.
 */
inline float readGridSize(const float* realGridSizeFP, const int dimIndex)
{
    switch (dimIndex)
    {
        case XX: return realGridSizeFP[XX];
        case YY: return realGridSizeFP[YY];
        case ZZ: return realGridSizeFP[ZZ];
    }
    SYCL_ASSERT(false);
    return 0.0F;
}


/*! \brief Reduce the partial force contributions.
 *
 * \tparam     order              The PME order (must be 4).
 * \tparam     atomDataSize       The number of partial force contributions for each atom (currently
 *                                order^2 == 16).
 * \tparam     workGroupSize      The size of a work-group.
 * \tparam     subGroupSize       The size of a sub-group.
 *
 * \param[in]  itemIdx            SYCL thread ID.
 * \param[out] sm_forces          Shared memory array with the output forces (number of elements
 *                                is number of atoms per block).
 * \param[in]  atomIndexLocal     Local atom index.
 * \param[in]  splineIndex        Spline index.
 * \param[in]  lineIndex          Line index (same as threadLocalId)
 * \param[in]  realGridSizeFP     Local grid size constant
 * \param[in]  fx                 Input force partial component X
 * \param[in]  fy                 Input force partial component Y
 * \param[in]  fz                 Input force partial component Z
 */
template<int order, int atomDataSize, int workGroupSize, int subGroupSize>
inline void reduceAtomForces(sycl::nd_item<3>        itemIdx,
                             sycl::local_ptr<Float3> sm_forces,
                             const int               atomIndexLocal,
                             const int               splineIndex,
                             const int gmx_unused    lineIndex,
                             const float             realGridSizeFP[3],
                             float&                  fx, // NOLINT(google-runtime-references)
                             float&                  fy, // NOLINT(google-runtime-references)
                             float&                  fz)                  // NOLINT(google-runtime-references)
{
    static_assert(gmx::isPowerOfTwo(order));
    // TODO: find out if this is the best in terms of transactions count
    static_assert(order == 4, "Only order of 4 is implemented");

    sycl::sub_group sg = itemIdx.get_sub_group();

    static_assert(atomDataSize <= subGroupSize,
                  "TODO: rework for atomDataSize > subGroupSize (order 8 or larger)");
    static_assert(gmx::isPowerOfTwo(atomDataSize));

    fx += sycl::shift_group_left(sg, fx, 1);
    fy += sycl::shift_group_right(sg, fy, 1);
    fz += sycl::shift_group_left(sg, fz, 1);
    if (splineIndex & 1)
    {
        fx = fy;
    }
    fx += sycl::shift_group_left(sg, fx, 2);
    fz += sycl::shift_group_right(sg, fz, 2);
    if (splineIndex & 2)
    {
        fx = fz;
    }
    static_assert(atomDataSize >= 4);
    // We have to just further reduce those groups of 4
    for (int delta = 4; delta < atomDataSize; delta *= 2)
    {
        fx += sycl::shift_group_left(sg, fx, delta);
    }
    const int dimIndex = splineIndex;
    if (dimIndex < DIM)
    {
        const float n                       = readGridSize(realGridSizeFP, dimIndex);
        sm_forces[atomIndexLocal][dimIndex] = fx * n;
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
 * \param[in] tdz                 The theta and dtheta in the Z dimension.
 * \param[in] sm_gridlineIndices  Shared memory array of grid line indices.
 * \param[in] sm_theta            Shared memory array of atom theta values.
 * \param[in] sm_dtheta           Shared memory array of atom dtheta values.
 * \param[in] gm_grid             Global memory array of the grid to use.
 */
template<int order, int atomsPerWarp, bool wrapX, bool wrapY>
inline void sumForceComponents(sycl::private_ptr<float>            fx,
                               sycl::private_ptr<float>            fy,
                               sycl::private_ptr<float>            fz,
                               const int                           ithyMin,
                               const int                           ithyMax,
                               const int                           ixBase,
                               const int                           iz,
                               const int                           nx,
                               const int                           ny,
                               const int                           pny,
                               const int                           pnz,
                               const int                           atomIndexLocal,
                               const int                           splineIndexBase,
                               const sycl::float2                  tdz,
                               const sycl::local_ptr<int>          sm_gridlineIndices,
                               const sycl::local_ptr<float>        sm_theta,
                               const sycl::local_ptr<float>        sm_dtheta,
                               const sycl::global_ptr<const float> gm_grid)
{
    for (int ithy = ithyMin; ithy < ithyMax; ithy++)
    {
        const int splineIndexY = getSplineParamIndex<order, atomsPerWarp>(splineIndexBase, YY, ithy);
        const sycl::float2 tdy{ sm_theta[splineIndexY], sm_dtheta[splineIndexY] };

        int iy = sm_gridlineIndices[atomIndexLocal * DIM + YY] + ithy;
        if (wrapY & (iy >= ny))
        {
            iy -= ny;
        }
        const int constOffset = iy * pnz + iz;

#pragma unroll
        for (int ithx = 0; ithx < order; ithx++)
        {
            int ix = ixBase + ithx;
            if (wrapX & (ix >= nx))
            {
                ix -= nx;
            }
            const int gridIndexGlobal = ix * pny * pnz + constOffset;
            SYCL_ASSERT(gridIndexGlobal >= 0);
            const float gridValue = gm_grid[gridIndexGlobal];
            assertIsFinite(gridValue);
            const int splineIndexX = getSplineParamIndex<order, atomsPerWarp>(splineIndexBase, XX, ithx);
            const sycl::float2 tdx{ sm_theta[splineIndexX], sm_dtheta[splineIndexX] };
            const float        fxy1 = tdz[XX] * gridValue;
            const float        fz1  = tdz[YY] * gridValue;
            *fx += tdx[YY] * tdy[XX] * fxy1;
            *fy += tdx[XX] * tdy[YY] * fxy1;
            *fz += tdx[XX] * tdy[XX] * fz1;
        }
    }
}


/*! \brief Calculate the grid forces and store them in shared memory.
 *
 * \param[in,out] sm_forces       Shared memory array with the output forces.
 * \param[in] forceIndexLocal     The local (per thread) index in the sm_forces array.
 * \param[in] forceIndexGlobal    The index of the thread in the gm_coefficients array.
 * \param[in] recipBox0           The reciprocal box (first vector).
 * \param[in] recipBox1           The reciprocal box (second vector).
 * \param[in] recipBox2           The reciprocal box (third vector).
 * \param[in] scale               The scale to use when calculating the forces. For gm_coefficientsB
 * (when using multiple coefficients on a single grid) the scale will be (1.0 - scale).
 * \param[in] gm_coefficients     Global memory array of the coefficients to use for an unperturbed
 * or FEP in state A if a single grid is used (\p multiCoefficientsSingleGrid == true).If two
 * separate grids are used this should be the coefficients of the grid in question.
 */
inline void calculateAndStoreGridForces(sycl::local_ptr<Float3>             sm_forces,
                                        const int                           forceIndexLocal,
                                        const int                           forceIndexGlobal,
                                        const Float3&                       recipBox0,
                                        const Float3&                       recipBox1,
                                        const Float3&                       recipBox2,
                                        const float                         scale,
                                        const sycl::global_ptr<const float> gm_coefficients)
{
    const Float3 atomForces     = sm_forces[forceIndexLocal];
    float        negCoefficient = -scale * gm_coefficients[forceIndexGlobal];
    Float3       result;
    result[XX] = negCoefficient * recipBox0[XX] * atomForces[XX];
    result[YY] = negCoefficient * (recipBox0[YY] * atomForces[XX] + recipBox1[YY] * atomForces[YY]);
    result[ZZ] = negCoefficient
                 * (recipBox0[ZZ] * atomForces[XX] + recipBox1[ZZ] * atomForces[YY]
                    + recipBox2[ZZ] * atomForces[ZZ]);
    sm_forces[forceIndexLocal] = result;
}

/*! \brief
 * A SYCL kernel which gathers the atom forces from the grid.
 * The grid is assumed to be wrapped in dimension Z.
 *
 * \tparam order          PME interpolation order.
 * \tparam wrapX          A boolean which tells if the grid overlap in dimension X should
 *                        be wrapped.
 * \tparam wrapY          A boolean which tells if the grid overlap in dimension Y should
 *                        be wrapped.
 * \tparam numGrids       The number of grids to use in the kernel. Can be 1 or 2.
 * \tparam readGlobal     Tells if we should read spline values from global memory.
 * \tparam threadsPerAtom How many threads work on each atom.
 * \tparam subGroupSize   Size of the sub-group.
 */
template<int order, bool wrapX, bool wrapY, int numGrids, bool readGlobal, ThreadsPerAtom threadsPerAtom, int subGroupSize>
auto pmeGatherKernel(sycl::handler& cgh,
                     const int      nAtoms,
                     const float* __restrict__ gm_gridA,
                     const float* __restrict__ gm_gridB,
                     const float* __restrict__ gm_coefficientsA,
                     const float* __restrict__ gm_coefficientsB,
                     const Float3* __restrict__ gm_coordinates,
                     Float3* __restrict__ gm_forces,
                     const float* __restrict__ gm_theta,
                     const float* __restrict__ gm_dtheta,
                     const int* __restrict__ gm_gridlineIndices,
                     const float* __restrict__ gm_fractShiftsTable,
                     const int* __restrict__ gm_gridlineIndicesTable,
                     const gmx::IVec tablesOffsets,
                     const gmx::IVec realGridSize,
                     const gmx::RVec realGridSizeFP,
                     const gmx::IVec realGridSizePadded,
                     const gmx::RVec currentRecipBox0,
                     const gmx::RVec currentRecipBox1,
                     const gmx::RVec currentRecipBox2,
                     const float     scale)
{
    static_assert(numGrids == 1 || numGrids == 2);

    constexpr int threadsPerAtomValue = (threadsPerAtom == ThreadsPerAtom::Order) ? order : order * order;
    constexpr int atomDataSize        = threadsPerAtomValue;
    constexpr int atomsPerBlock       = (c_gatherMaxWarpsPerBlock * subGroupSize) / atomDataSize;
    // Number of atoms processed by a single warp in spread and gather
    static_assert(subGroupSize >= atomDataSize);
    constexpr int atomsPerWarp        = subGroupSize / atomDataSize;
    constexpr int blockSize           = atomsPerBlock * atomDataSize;
    constexpr int splineParamsSize    = atomsPerBlock * DIM * order;
    constexpr int gridlineIndicesSize = atomsPerBlock * DIM;

    // Gridline indices, ivec
    sycl::local_accessor<int, 1> sm_gridlineIndices(sycl::range<1>(atomsPerBlock * DIM), cgh);
    // Spline values
    sycl::local_accessor<float, 1> sm_theta(sycl::range<1>(atomsPerBlock * DIM * order), cgh);
    // Spline derivatives
    sycl::local_accessor<float, 1> sm_dtheta(sycl::range<1>(atomsPerBlock * DIM * order), cgh);
    // Coefficients prefetch cache
    sycl::local_accessor<float, 1> sm_coefficients(sycl::range<1>(atomsPerBlock), cgh);
    // Coordinates prefetch cache
    sycl::local_accessor<Float3, 1> sm_coordinates(sycl::range<1>(atomsPerBlock), cgh);
    // Reduction of partial force contributions
    sycl::local_accessor<Float3, 1> sm_forces(sycl::range<1>(atomsPerBlock), cgh);

    auto sm_fractCoords = [&]() {
        if constexpr (!readGlobal)
        {
            return sycl::local_accessor<float, 1>(sycl::range<1>(atomsPerBlock * DIM), cgh);
        }
        else
        {
            return nullptr;
        }
    }();

    return [=](sycl::nd_item<3> itemIdx) [[intel::reqd_sub_group_size(subGroupSize)]]
    {
        if constexpr (skipKernelCompilation<subGroupSize>())
        {
            return;
        }
        SYCL_ASSERT(blockSize == itemIdx.get_local_range().size());
        /* These are the atom indices - for the shared and global memory */
        const int atomIndexLocal = itemIdx.get_local_id(XX);
        const int blockIndex =
                itemIdx.get_group(YY) * itemIdx.get_group_range(ZZ) + itemIdx.get_group(ZZ);
        // itemIdx.get_group_linear_id();

        const int atomIndexOffset = blockIndex * atomsPerBlock;
        const int atomIndexGlobal = atomIndexOffset + atomIndexLocal;
        /* Early return for fully empty blocks at the end
         * (should only happen for billions of input atoms)
         */
        if (atomIndexOffset >= nAtoms)
        {
            return;
        }
        /* Spline Z coordinates */
        const int ithz = itemIdx.get_local_id(ZZ);
        /* These are the spline contribution indices in shared memory */
        const int splineIndex =
                itemIdx.get_local_id(YY) * itemIdx.get_local_range(ZZ) + itemIdx.get_local_id(ZZ);

        const int threadLocalId    = itemIdx.get_local_linear_id();
        const int threadLocalIdMax = blockSize;
        SYCL_ASSERT(threadLocalId < threadLocalIdMax);

        const int lineIndex =
                (itemIdx.get_local_id(XX) * (itemIdx.get_local_range(ZZ) * itemIdx.get_local_range(YY)))
                + splineIndex; // And to all the block's particles
        SYCL_ASSERT(lineIndex == threadLocalId);

        if constexpr (readGlobal)
        {
            /* Read splines */
            const int localGridlineIndicesIndex = threadLocalId;
            const int globalGridlineIndicesIndex =
                    blockIndex * gridlineIndicesSize + localGridlineIndicesIndex;
            // itemIdx.get_group(ZZ) * gridlineIndicesSize + localGridlineIndicesIndex;
            if (localGridlineIndicesIndex < gridlineIndicesSize)
            {
                sm_gridlineIndices[localGridlineIndicesIndex] =
                        gm_gridlineIndices[globalGridlineIndicesIndex];
                SYCL_ASSERT(sm_gridlineIndices[localGridlineIndicesIndex] >= 0);
            }
            /* The loop needed for order threads per atom to make sure we load all data values, as each thread must load multiple values
               with order*order threads per atom, it is only required for each thread to load one data value */

            const int iMin = 0;
            const int iMax = (threadsPerAtom == ThreadsPerAtom::Order) ? 3 : 1;

            for (int i = iMin; i < iMax; i++)
            {
                // i will always be zero for order*order threads per atom
                const int localSplineParamsIndex = threadLocalId + i * threadLocalIdMax;
                const int globalSplineParamsIndex = blockIndex * splineParamsSize + localSplineParamsIndex;
                // const int globalSplineParamsIndex = itemIdx.get_group(ZZ) * splineParamsSize + localSplineParamsIndex;
                if (localSplineParamsIndex < splineParamsSize)
                {
                    sm_theta[localSplineParamsIndex]  = gm_theta[globalSplineParamsIndex];
                    sm_dtheta[localSplineParamsIndex] = gm_dtheta[globalSplineParamsIndex];
                    assertIsFinite(sm_theta[localSplineParamsIndex]);
                    assertIsFinite(sm_dtheta[localSplineParamsIndex]);
                }
            }

            itemIdx.barrier(fence_space::local_space);
        }
        else
        {
            /* Recalculate  Splines  */
            /* Staging coefficients/charges */
            pmeGpuStageAtomData<float, atomsPerBlock, 1>(
                    sm_coefficients.get_pointer(), gm_coefficientsA, itemIdx);
            /* Staging coordinates */
            pmeGpuStageAtomData<Float3, atomsPerBlock, 1>(
                    sm_coordinates.get_pointer(), gm_coordinates, itemIdx);
            itemIdx.barrier(fence_space::local_space);
            const Float3 atomX      = sm_coordinates[atomIndexLocal];
            const float  atomCharge = sm_coefficients[atomIndexLocal];

            calculateSplines<order, atomsPerBlock, atomsPerWarp, true, false, numGrids, subGroupSize>(
                    atomIndexOffset,
                    atomX,
                    atomCharge,
                    tablesOffsets,
                    realGridSizeFP,
                    currentRecipBox0,
                    currentRecipBox1,
                    currentRecipBox2,
                    nullptr,
                    nullptr,
                    nullptr,
                    gm_fractShiftsTable,
                    gm_gridlineIndicesTable,
                    sm_theta.get_pointer(),
                    sm_dtheta.get_pointer(),
                    sm_gridlineIndices.get_pointer(),
                    sm_fractCoords.get_pointer(),
                    itemIdx);
            subGroupBarrier(itemIdx);
        }
        float fx = 0.0F;
        float fy = 0.0F;
        float fz = 0.0F;

        const int chargeCheck = pmeGpuCheckAtomCharge(gm_coefficientsA[atomIndexGlobal]);

        const int nx  = realGridSize[XX];
        const int ny  = realGridSize[YY];
        const int nz  = realGridSize[ZZ];
        const int pny = realGridSizePadded[YY];
        const int pnz = realGridSizePadded[ZZ];

        const int atomWarpIndex = atomIndexLocal % atomsPerWarp;
        const int warpIndex     = atomIndexLocal / atomsPerWarp;

        const int splineIndexBase = getSplineParamIndexBase<order, atomsPerWarp>(warpIndex, atomWarpIndex);
        const int splineIndexZ = getSplineParamIndex<order, atomsPerWarp>(splineIndexBase, ZZ, ithz);
        const sycl::float2 tdz{ sm_theta[splineIndexZ], sm_dtheta[splineIndexZ] };

        int       iz     = sm_gridlineIndices[atomIndexLocal * DIM + ZZ] + ithz;
        const int ixBase = sm_gridlineIndices[atomIndexLocal * DIM + XX];

        if (iz >= nz)
        {
            iz -= nz;
        }

        const int ithyMin = (threadsPerAtom == ThreadsPerAtom::Order) ? 0 : itemIdx.get_local_id(YY);
        const int ithyMax =
                (threadsPerAtom == ThreadsPerAtom::Order) ? order : itemIdx.get_local_id(YY) + 1;
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
                                                                  tdz,
                                                                  sm_gridlineIndices.get_pointer(),
                                                                  sm_theta.get_pointer(),
                                                                  sm_dtheta.get_pointer(),
                                                                  gm_gridA);
        }
        reduceAtomForces<order, atomDataSize, blockSize, subGroupSize>(
                itemIdx, sm_forces.get_pointer(), atomIndexLocal, splineIndex, lineIndex, realGridSizeFP, fx, fy, fz);
        itemIdx.barrier(fence_space::local_space);

        /* Calculating the final forces with no component branching, atomsPerBlock threads */
        const int forceIndexLocal  = threadLocalId;
        const int forceIndexGlobal = atomIndexOffset + forceIndexLocal;
        if (forceIndexLocal < atomsPerBlock)
        {
            calculateAndStoreGridForces(sm_forces.get_pointer(),
                                        forceIndexLocal,
                                        forceIndexGlobal,
                                        currentRecipBox0,
                                        currentRecipBox1,
                                        currentRecipBox2,
                                        scale,
                                        gm_coefficientsA);
        }
        itemIdx.barrier(fence_space::local_space);

        static_assert(atomsPerBlock <= subGroupSize);

        /* Writing or adding the final forces component-wise, single warp */
        constexpr int blockForcesSize = atomsPerBlock * DIM;
        constexpr int numIter         = (blockForcesSize + subGroupSize - 1) / subGroupSize;
        constexpr int iterThreads     = blockForcesSize / numIter;
        if (threadLocalId < iterThreads)
        {
#pragma unroll
            for (int i = 0; i < numIter; i++)
            {
                const int floatIndexLocal  = i * iterThreads + threadLocalId;
                const int float3IndexLocal = floatIndexLocal / 3;
                const int dimLocal         = floatIndexLocal % 3;
                static_assert(blockForcesSize % DIM == 0); // Assures that dimGlobal == dimLocal
                const int float3IndexGlobal = blockIndex * atomsPerBlock + float3IndexLocal;
                // const int float3IndexGlobal = itemIdx.get_group(ZZ) * atomsPerBlock + float3IndexLocal;
                gm_forces[float3IndexGlobal][dimLocal] = sm_forces[float3IndexLocal][dimLocal];
            }
        }

        if constexpr (numGrids == 2)
        {
            /* We must sync here since the same shared memory is used as above. */
            itemIdx.barrier(fence_space::local_space);
            fx                     = 0.0F;
            fy                     = 0.0F;
            fz                     = 0.0F;
            const bool chargeCheck = pmeGpuCheckAtomCharge(gm_coefficientsB[atomIndexGlobal]);
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
                                                                      tdz,
                                                                      sm_gridlineIndices.get_pointer(),
                                                                      sm_theta.get_pointer(),
                                                                      sm_dtheta.get_pointer(),
                                                                      gm_gridB);
            }
            // Reduction of partial force contributions
            reduceAtomForces<order, atomDataSize, blockSize, subGroupSize>(
                    itemIdx, sm_forces.get_pointer(), atomIndexLocal, splineIndex, lineIndex, realGridSizeFP, fx, fy, fz);
            itemIdx.barrier(fence_space::local_space);

            /* Calculating the final forces with no component branching, atomsPerBlock threads */
            if (forceIndexLocal < atomsPerBlock)
            {
                calculateAndStoreGridForces(sm_forces.get_pointer(),
                                            forceIndexLocal,
                                            forceIndexGlobal,
                                            currentRecipBox0,
                                            currentRecipBox1,
                                            currentRecipBox2,
                                            1.0F - scale,
                                            gm_coefficientsB);
            }

            itemIdx.barrier(fence_space::local_space);

            /* Writing or adding the final forces component-wise, single warp */
            if (threadLocalId < iterThreads)
            {
#pragma unroll
                for (int i = 0; i < numIter; i++)
                {
                    const int floatIndexLocal  = i * iterThreads + threadLocalId;
                    const int float3IndexLocal = floatIndexLocal / 3;
                    const int dimLocal         = floatIndexLocal % 3;
                    static_assert(blockForcesSize % DIM == 0); // Assures that dimGlobal == dimLocal
                    const int float3IndexGlobal = blockIndex * atomsPerBlock + float3IndexLocal;
                    gm_forces[float3IndexGlobal][dimLocal] += sm_forces[float3IndexLocal][dimLocal];
                }
            }
        }
    };
}

template<int order, bool wrapX, bool wrapY, int numGrids, bool readGlobal, ThreadsPerAtom threadsPerAtom, int subGroupSize>
PmeGatherKernel<order, wrapX, wrapY, numGrids, readGlobal, threadsPerAtom, subGroupSize>::PmeGatherKernel()
{
    reset();
}

template<int order, bool wrapX, bool wrapY, int numGrids, bool readGlobal, ThreadsPerAtom threadsPerAtom, int subGroupSize>
void PmeGatherKernel<order, wrapX, wrapY, numGrids, readGlobal, threadsPerAtom, subGroupSize>::setArg(
        size_t argIndex,
        void*  arg)
{
    if (argIndex == 0)
    {
        auto* params   = reinterpret_cast<PmeGpuKernelParams*>(arg);
        gridParams_    = &params->grid;
        atomParams_    = &params->atoms;
        dynamicParams_ = &params->current;
    }
    else
    {
        GMX_RELEASE_ASSERT(argIndex == 0, "Trying to pass too many args to the kernel");
    }
}


template<int order, bool wrapX, bool wrapY, int numGrids, bool readGlobal, ThreadsPerAtom threadsPerAtom, int subGroupSize>
void PmeGatherKernel<order, wrapX, wrapY, numGrids, readGlobal, threadsPerAtom, subGroupSize>::launch(
        const KernelLaunchConfig& config,
        const DeviceStream&       deviceStream)
{
    GMX_RELEASE_ASSERT(gridParams_, "Can not launch the kernel before setting its args");
    GMX_RELEASE_ASSERT(atomParams_, "Can not launch the kernel before setting its args");
    GMX_RELEASE_ASSERT(dynamicParams_, "Can not launch the kernel before setting its args");

    using kernelNameType =
            PmeGatherKernel<order, wrapX, wrapY, numGrids, readGlobal, threadsPerAtom, subGroupSize>;

    // SYCL has different multidimensional layout than OpenCL/CUDA.
    const sycl::range<3> localSize{ config.blockSize[2], config.blockSize[1], config.blockSize[0] };
    const sycl::range<3> groupRange{ config.gridSize[2], config.gridSize[1], config.gridSize[0] };
    const sycl::nd_range<3> range{ groupRange * localSize, localSize };

    sycl::queue q = deviceStream.stream();

    q.submit(GMX_SYCL_DISCARD_EVENT[&](sycl::handler & cgh) {
        auto kernel = pmeGatherKernel<order, wrapX, wrapY, numGrids, readGlobal, threadsPerAtom, subGroupSize>(
                cgh,
                atomParams_->nAtoms,
                gridParams_->d_realGrid[0].get_pointer(),
                gridParams_->d_realGrid[1].get_pointer(),
                atomParams_->d_coefficients[0].get_pointer(),
                atomParams_->d_coefficients[1].get_pointer(),
                atomParams_->d_coordinates.get_pointer(),
                atomParams_->d_forces.get_pointer(),
                atomParams_->d_theta.get_pointer(),
                atomParams_->d_dtheta.get_pointer(),
                atomParams_->d_gridlineIndices.get_pointer(),
                gridParams_->d_fractShiftsTable.get_pointer(),
                gridParams_->d_gridlineIndicesTable.get_pointer(),
                gridParams_->tablesOffsets,
                gridParams_->realGridSize,
                gridParams_->realGridSizeFP,
                gridParams_->realGridSizePadded,
                dynamicParams_->recipBox[0],
                dynamicParams_->recipBox[1],
                dynamicParams_->recipBox[2],
                dynamicParams_->scale);
        cgh.parallel_for<kernelNameType>(range, kernel);
    });

    // Delete set args, so we don't forget to set them before the next launch.
    reset();
}


template<int order, bool wrapX, bool wrapY, int numGrids, bool readGlobal, ThreadsPerAtom threadsPerAtom, int subGroupSize>
void PmeGatherKernel<order, wrapX, wrapY, numGrids, readGlobal, threadsPerAtom, subGroupSize>::reset()
{
    gridParams_    = nullptr;
    atomParams_    = nullptr;
    dynamicParams_ = nullptr;
}


//! Kernel instantiations
/* Disable the "explicit template instantiation 'PmeSplineAndSpreadKernel<...>' will emit a vtable in every
 * translation unit [-Wweak-template-vtables]" warning.
 * It is only explicitly instantiated in this translation unit, so we should be safe.
 */
CLANG_DIAGNOSTIC_IGNORE("-Wweak-template-vtables")

#define INSTANTIATE_3(order, numGrids, readGlobal, threadsPerAtom, subGroupSize) \
    template class PmeGatherKernel<order, true, true, numGrids, readGlobal, threadsPerAtom, subGroupSize>;

#define INSTANTIATE_2(order, numGrids, threadsPerAtom, subGroupSize)    \
    INSTANTIATE_3(order, numGrids, true, threadsPerAtom, subGroupSize); \
    INSTANTIATE_3(order, numGrids, false, threadsPerAtom, subGroupSize);

#define INSTANTIATE(order, subGroupSize)                                 \
    INSTANTIATE_2(order, 1, ThreadsPerAtom::Order, subGroupSize);        \
    INSTANTIATE_2(order, 1, ThreadsPerAtom::OrderSquared, subGroupSize); \
    INSTANTIATE_2(order, 2, ThreadsPerAtom::Order, subGroupSize);        \
    INSTANTIATE_2(order, 2, ThreadsPerAtom::OrderSquared, subGroupSize);

#if GMX_SYCL_DPCPP
INSTANTIATE(4, 16); // TODO: Choose best value, Issue #4153.
#endif
INSTANTIATE(4, 32);
INSTANTIATE(4, 64);

CLANG_DIAGNOSTIC_RESET
