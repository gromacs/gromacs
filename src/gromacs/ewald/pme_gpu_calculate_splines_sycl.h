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
 *  \brief Implements helper routines for PME gather and spline routines.
 *
 *  \author Andrey Alekseenko <al42and@gmail.com>
 */

#include "gmxpre.h"

#include <cassert>

#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/gpu_utils/gputraits_sycl.h"
#include "gromacs/gpu_utils/sycl_kernel_utils.h"

#include "pme_gpu_constants.h"
#include "pme_gpu_types.h"
#include "pme_grid.h"

namespace
{

/*! \brief Asserts if the argument is finite.
 *
 *  The function works for any data type, that can be casted to float. Note that there is also
 *  a specialized implementation for float3 data type.
 *
 * \param[in] arg  Argument to check.
 */
template<typename T>
inline void assertIsFinite(T arg);

#if defined(NDEBUG) || GMX_SYCL_ACPP
// We have no sycl::isfinite in AdaptiveCpp yet
template<typename T>
inline void assertIsFinite(T /* arg */)
{
}
#else
template<>
inline void assertIsFinite(Float3 gmx_unused arg)
{
    SYCL_ASSERT(sycl::isfinite(arg[0]));
    SYCL_ASSERT(sycl::isfinite(arg[1]));
    SYCL_ASSERT(sycl::isfinite(arg[2]));
}

template<typename T>
inline void assertIsFinite(T gmx_unused arg)
{
    SYCL_ASSERT(sycl::isfinite(static_cast<float>(arg)));
}
#endif

} // namespace

using mode = sycl::access_mode;
using sycl::access::fence_space;

/*! \internal \brief
 * Gets a base of the unique index to an element in a spline parameter buffer (theta/dtheta),
 * which is laid out for GPU spread/gather kernels. The base only corresponds to the atom index within the execution block.
 * Feed the result into getSplineParamIndex() to get a full index.
 * TODO: it's likely that both parameters can be just replaced with a single atom index, as they are derived from it.
 * Do that, verifying that the generated code is not bloated, and/or revise the spline indexing scheme.
 * Removing warp dependency would also be nice (and would probably coincide with removing c_pmeSpreadGatherAtomsPerWarp).
 *
 * \tparam order                 PME order
 * \tparam atomsPerSubGroup      Number of atoms processed by a sub group
 * \param[in] subGroupIndex      Sub group index in the work group.
 * \param[in] atomSubGroupIndex  Atom index in the sub group (from 0 to atomsPerSubGroup - 1).
 *
 * \returns Index into theta or dtheta array using GPU layout.
 */
template<int order, int atomsPerSubGroup>
static inline int getSplineParamIndexBase(int subGroupIndex, int atomSubGroupIndex)
{
    SYCL_ASSERT((atomSubGroupIndex >= 0) && (atomSubGroupIndex < atomsPerSubGroup));
    constexpr int dimIndex    = 0;
    constexpr int splineIndex = 0;
    // The zeroes are here to preserve the full index formula for reference
    return (((splineIndex + order * subGroupIndex) * DIM + dimIndex) * atomsPerSubGroup + atomSubGroupIndex);
}

/*! \internal \brief
 * Gets a unique index to an element in a spline parameter buffer (theta/dtheta),
 * which is laid out for GPU spread/gather kernels. The index is wrt to the execution block,
 * in range(0, atomsPerBlock * order * DIM).
 * This function consumes result of getSplineParamIndexBase() and adjusts it for \p dimIndex and \p splineIndex.
 *
 * \tparam order               PME order
 * \tparam atomsPerSubGroup    Number of atoms processed by a sub group
 * \param[in] paramIndexBase   Must be result of getSplineParamIndexBase().
 * \param[in] dimIndex         Dimension index (from 0 to 2)
 * \param[in] splineIndex      Spline contribution index (from 0 to \p order - 1)
 *
 * \returns Index into theta or dtheta array using GPU layout.
 */
template<int order, int atomsPerSubGroup>
static inline int getSplineParamIndex(int paramIndexBase, int dimIndex, int splineIndex)
{
    SYCL_ASSERT((dimIndex >= XX) && (dimIndex < DIM));
    SYCL_ASSERT((splineIndex >= 0) && (splineIndex < order));
    return (paramIndexBase + (splineIndex * DIM + dimIndex) * atomsPerSubGroup);
}

/*! \internal \brief
 * An inline function for skipping the zero-charge atoms when we have \c c_skipNeutralAtoms set to \c true.
 *
 * \returns                   \c true if atom should be processed, \c false otherwise.
 * \param[in] charge          The atom charge.
 */
static inline bool pmeGpuCheckAtomCharge(const float charge)
{
    assertIsFinite(charge);
    return c_skipNeutralAtoms ? (charge != 0.0F) : true;
}

/*! \brief
 * General purpose function for loading atom-related data from global to shared memory.
 *
 * \tparam T Data type (float/int/...).
 * \tparam atomsPerWorkGroup Number of atoms processed by a block - should be accounted for
 *                           in the size of the shared memory array.
 * \tparam dataCountPerAtom Number of data elements
 *                          per single atom (e.g. \c DIM for an rvec coordinates array).
 * \param[out] sm_destination Shared memory array for output.
 * \param[in]  gm_source Global memory array for input.
 * \param[in]  itemIdx SYCL thread ID.
 */
template<typename T, int atomsPerWorkGroup, int dataCountPerAtom>
static inline void pmeGpuStageAtomData(sycl::local_ptr<T>              sm_destination,
                                       const sycl::global_ptr<const T> gm_source,
                                       sycl::nd_item<3>                itemIdx)
{
    const int blockIndex      = itemIdx.get_group_linear_id();
    const int localIndex      = itemIdx.get_local_linear_id();
    const int globalIndexBase = blockIndex * atomsPerWorkGroup * dataCountPerAtom;
    const int globalIndex     = globalIndexBase + localIndex;
    if (localIndex < atomsPerWorkGroup * dataCountPerAtom)
    {
        assertIsFinite(gm_source[globalIndex]);
        sm_destination[localIndex] = gm_source[globalIndex];
    }
}

/*! \brief
 * PME GPU spline parameter and gridline indices calculation.
 * This corresponds to the CPU functions calc_interpolation_idx() and make_bsplines().
 * First stage of the whole kernel.
 *
 * \tparam order                PME interpolation order.
 * \tparam atomsPerBlock        Number of atoms processed by a block - should be accounted for
 *                              in the sizes of the shared memory arrays.
 * \tparam atomsPerWarp         Number of atoms processed by a warp
 * \tparam writeSmDtheta        Bool controlling if the theta derivative should be written to
 *                              shared memory. Enables calculation of dtheta if set.
 * \tparam writeGlobal          A boolean which tells if the theta values and gridlines should
 *                              be written to global memory. Enables calculation of dtheta if set.
 * \tparam numGrids             The number of grids using the splines.
 * \tparam subGroupSize         The size of a sub-group (warp).
 * \param[in]  atomIndexOffset        Starting atom index for the execution block in the global
 *                                    memory.
 * \param[in]  atomX                  Coordinates of atom processed by thread.
 * \param[in]  atomCharge             Charge/coefficient of atom processed by thread.
 * \param[in]  tablesOffsets          Offsets for X/Y/Z components of \p gm_fractShiftsTable and
 *                                    \p gm_gridlineIndicesTable.
 * \param[in]  realGridSizeFP         Real-space grid dimensions, converted to floating point.
 * \param[in]  currentRecipBox0       Current reciprocal (inverted unit cell) box, vector 1.
 * \param[in]  currentRecipBox1       Current reciprocal (inverted unit cell) box, vector 2.
 * \param[in]  currentRecipBox2       Current reciprocal (inverted unit cell) box, vector 3.
 * \param[out] gm_theta               Atom spline values in the global memory.
 *                                    Used only if \p writeGlobal is \c true.
 * \param[out] gm_dtheta              Derivatives of atom spline values in the global memory.
 *                                    Used only if \p writeGlobal is \c true.
 * \param[out] gm_gridlineIndices     Atom gridline indices in the global memory.
 *                                    Used only if \p writeGlobal is \c true.
 * \param[in] gm_fractShiftsTable     Fractional shifts lookup table in the global memory.
 * \param[in] gm_gridlineIndicesTable Gridline indices lookup table in the global memory.
 * \param[out] sm_theta               Atom spline values in the local memory.
 * \param[out] sm_dtheta              Derivatives of atom spline values in the local memory.
 * \param[out] sm_gridlineIndices     Atom gridline indices in the local memory.
 * \param[out] sm_fractCoords         Fractional coordinates in the local memory.
 * \param[in]  itemIdx                SYCL thread ID.
 */

template<int order, int atomsPerBlock, int atomsPerWarp, bool writeSmDtheta, bool writeGlobal, int numGrids, int subGroupSize>
static inline void calculateSplines(const int                           atomIndexOffset,
                                    const Float3                        atomX,
                                    const float                         atomCharge,
                                    const gmx::IVec                     tablesOffsets,
                                    const gmx::RVec                     realGridSizeFP,
                                    const gmx::RVec                     currentRecipBox0,
                                    const gmx::RVec                     currentRecipBox1,
                                    const gmx::RVec                     currentRecipBox2,
                                    sycl::global_ptr<float>             gm_theta,
                                    sycl::global_ptr<float>             gm_dtheta,
                                    sycl::global_ptr<int>               gm_gridlineIndices,
                                    const sycl::global_ptr<const float> gm_fractShiftsTable,
                                    const sycl::global_ptr<const int>   gm_gridlineIndicesTable,
                                    sycl::local_ptr<float>              sm_theta,
                                    sycl::local_ptr<float>              sm_dtheta,
                                    sycl::local_ptr<int>                sm_gridlineIndices,
                                    sycl::local_ptr<float>              sm_fractCoords,
                                    sycl::nd_item<3>                    itemIdx)
{
    static_assert(numGrids == 1 || numGrids == 2);
    static_assert(numGrids == 1 || c_skipNeutralAtoms == false);

    /* Thread index w.r.t. block */
    const int threadLocalId = itemIdx.get_local_linear_id();
    /* Warp index w.r.t. block - could probably be obtained easier? */
    const int warpIndex = threadLocalId / subGroupSize;
    /* Atom index w.r.t. warp - alternating 0 1 0 1 ... */
    const int atomWarpIndex = itemIdx.get_local_id(0) % atomsPerWarp;
    /* Atom index w.r.t. block/shared memory */
    const int atomIndexLocal = warpIndex * atomsPerWarp + atomWarpIndex;

    /* Spline contribution index in one dimension */
    const int threadLocalIdXY =
            (itemIdx.get_local_id(1) * itemIdx.get_local_range(2)) + itemIdx.get_local_id(2);
    const int orderIndex = threadLocalIdXY / DIM;
    /* Dimension index */
    const int dimIndex = threadLocalIdXY % DIM;

    /* Multi-purpose index of rvec/ivec atom data */
    const int sharedMemoryIndex = atomIndexLocal * DIM + dimIndex;

    float splineData[order];

    const int localCheck = (dimIndex < DIM) && (orderIndex < 1);

    /* we have 4 threads per atom, but can only use 3 here for the dimensions */
    if (localCheck)
    {
        /* Indices interpolation */
        if (orderIndex == 0)
        {
            int   tableIndex, tInt;
            float n, t;
            SYCL_ASSERT(atomIndexLocal < DIM * atomsPerBlock);
            // Switch structure inherited from CUDA.
            // TODO: Issue #4153: Direct indexing with dimIndex can be better with SYCL
            switch (dimIndex)
            {
                case XX:
                    tableIndex = tablesOffsets[XX];
                    n          = realGridSizeFP[XX];
                    t          = atomX[XX] * currentRecipBox0[XX] + atomX[YY] * currentRecipBox0[YY]
                        + atomX[ZZ] * currentRecipBox0[ZZ];
                    break;

                case YY:
                    tableIndex = tablesOffsets[YY];
                    n          = realGridSizeFP[YY];
                    t = atomX[YY] * currentRecipBox1[YY] + atomX[ZZ] * currentRecipBox1[ZZ];
                    break;

                case ZZ:
                    tableIndex = tablesOffsets[ZZ];
                    n          = realGridSizeFP[ZZ];
                    t          = atomX[ZZ] * currentRecipBox2[ZZ];
                    break;
            }
            const float shift = c_pmeMaxUnitcellShift;
            /* Fractional coordinates along box vectors, adding a positive shift to ensure t is positive for triclinic boxes */
            t    = (t + shift) * n;
            tInt = static_cast<int>(t);
            SYCL_ASSERT(sharedMemoryIndex < atomsPerBlock * DIM);
            sm_fractCoords[sharedMemoryIndex] = t - tInt;
            tableIndex += tInt;
            SYCL_ASSERT(tInt >= 0);
            SYCL_ASSERT(tInt < c_pmeNeighborUnitcellCount * n);

            // TODO: Issue #4153: use shared table for both parameters to share the fetch, as index is always same.
            sm_fractCoords[sharedMemoryIndex] += gm_fractShiftsTable[tableIndex];
            sm_gridlineIndices[sharedMemoryIndex] = gm_gridlineIndicesTable[tableIndex];
            if constexpr (writeGlobal)
            {
                gm_gridlineIndices[atomIndexOffset * DIM + sharedMemoryIndex] =
                        sm_gridlineIndices[sharedMemoryIndex];
            }
        }

        /* B-spline calculation */
        const int chargeCheck = pmeGpuCheckAtomCharge(atomCharge);
        /* With FEP (numGrids == 2), we might have 0 charge in state A, but !=0 in state B, so we always calculate splines */
        if (numGrids == 2 || chargeCheck)
        {
            const float dr = sm_fractCoords[sharedMemoryIndex];
            assertIsFinite(dr);

            /* dr is relative offset from lower cell limit */
            splineData[order - 1] = 0.0F;
            splineData[1]         = dr;
            splineData[0]         = 1.0F - dr;

#pragma unroll
            for (int k = 3; k < order; k++)
            {
                const float div   = 1.0F / (k - 1.0F);
                splineData[k - 1] = div * dr * splineData[k - 2];
#pragma unroll
                for (int l = 1; l < (k - 1); l++)
                {
                    splineData[k - l - 1] =
                            div * ((dr + l) * splineData[k - l - 2] + (k - l - dr) * splineData[k - l - 1]);
                }
                splineData[0] = div * (1.0F - dr) * splineData[0];
            }

            const int thetaIndexBase =
                    getSplineParamIndexBase<order, atomsPerWarp>(warpIndex, atomWarpIndex);
            const int thetaGlobalOffsetBase = atomIndexOffset * DIM * order;
            /* only calculate dtheta if we are saving it to shared or global memory */
            if constexpr (writeSmDtheta || writeGlobal)
            {
                /* Differentiation and storing the spline derivatives (dtheta) */
#pragma unroll
                for (int o = 0; o < order; o++)
                {
                    const int thetaIndex =
                            getSplineParamIndex<order, atomsPerWarp>(thetaIndexBase, dimIndex, o);

                    const float dtheta = ((o > 0) ? splineData[o - 1] : 0.0F) - splineData[o];
                    assertIsFinite(dtheta);
                    SYCL_ASSERT(thetaIndex < order * DIM * atomsPerBlock);
                    if constexpr (writeSmDtheta)
                    {
                        sm_dtheta[thetaIndex] = dtheta;
                    }
                    if constexpr (writeGlobal)
                    {
                        const int thetaGlobalIndex  = thetaGlobalOffsetBase + thetaIndex;
                        gm_dtheta[thetaGlobalIndex] = dtheta;
                    }
                }
            }

            const float div       = 1.0F / (order - 1.0F);
            splineData[order - 1] = div * dr * splineData[order - 2];
#pragma unroll
            for (int k = 1; k < (order - 1); k++)
            {
                splineData[order - k - 1] = div
                                            * ((dr + k) * splineData[order - k - 2]
                                               + (order - k - dr) * splineData[order - k - 1]);
            }
            splineData[0] = div * (1.0F - dr) * splineData[0];

            /* Storing the spline values (theta) */
#pragma unroll
            for (int o = 0; o < order; o++)
            {
                const int thetaIndex =
                        getSplineParamIndex<order, atomsPerWarp>(thetaIndexBase, dimIndex, o);
                SYCL_ASSERT(thetaIndex < order * DIM * atomsPerBlock);
                sm_theta[thetaIndex] = splineData[o];
                assertIsFinite(sm_theta[thetaIndex]);
                if constexpr (writeGlobal)
                {
                    const int thetaGlobalIndex = thetaGlobalOffsetBase + thetaIndex;
                    gm_theta[thetaGlobalIndex] = splineData[o];
                }
            }
        }
    }
}
