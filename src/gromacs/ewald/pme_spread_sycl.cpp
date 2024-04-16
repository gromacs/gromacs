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
 *  \brief Implements PME GPU spline calculation and charge spreading in SYCL.
 *
 *  \author Andrey Alekseenko <al42and@gmail.com>
 */

#include "gmxpre.h"

#include "pme_spread_sycl.h"

#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/gpu_utils/gputraits_sycl.h"
#include "gromacs/gpu_utils/sycl_kernel_utils.h"
#include "gromacs/gpu_utils/syclutils.h"
#include "gromacs/utility/basedefinitions.h"

#include "pme_gpu_calculate_splines_sycl.h"
#include "pme_gpu_types_host.h"
#include "pme_grid.h"

/*! \brief
 * Charge spreading onto the grid.
 * This corresponds to the CPU function spread_coefficients_bsplines_thread().
 * Optional second stage of the spline_and_spread_kernel.
 *
 * \tparam     order                PME interpolation order.
 * \tparam     wrapX                Whether the grid overlap in dimension X should be wrapped.
 * \tparam     wrapY                Whether the grid overlap in dimension Y should be wrapped.
 * \tparam     threadsPerAtom       How many threads work on each atom.
 * \tparam     subGroupSize         Size of the sub-group.
 *
 * \param[in]  atomCharge           Atom charge/coefficient of atom processed by thread.
 * \param[in]  realGridSize         Size of the real grid.
 * \param[in]  realGridSizePadded   Padded of the real grid.
 * \param[in,out]  gm_grid          Device pointer to the real grid to which charges are added.
 * \param[in]  sm_gridlineIndices   Atom gridline indices in the local memory.
 * \param[in]  sm_theta             Atom spline values in the local memory.
 * \param[in]  itemIdx              Current thread ID.
 */
template<int order, bool wrapX, bool wrapY, ThreadsPerAtom threadsPerAtom, int subGroupSize>
inline void spread_charges(const float                  atomCharge,
                           const int                    realGridSize[DIM],
                           const int                    realGridSizePadded[DIM],
                           sycl::global_ptr<float>      gm_grid,
                           const sycl::local_ptr<int>   sm_gridlineIndices,
                           const sycl::local_ptr<float> sm_theta,
                           const sycl::nd_item<3>&      itemIdx)
{
    // Number of atoms processed by a single warp in spread and gather
    const int threadsPerAtomValue = (threadsPerAtom == ThreadsPerAtom::Order) ? order : order * order;
    const int atomsPerWarp        = subGroupSize / threadsPerAtomValue;

    const int nx  = realGridSize[XX];
    const int ny  = realGridSize[YY];
    const int nz  = realGridSize[ZZ];
    const int pny = realGridSizePadded[YY];
    const int pnz = realGridSizePadded[ZZ];

    const int atomIndexLocal = itemIdx.get_local_id(0);

    const int chargeCheck = pmeGpuCheckAtomCharge(atomCharge);

    if (chargeCheck)
    {
        // Spline Z coordinates
        const int ithz = itemIdx.get_local_id(2);

        const int ixBase = sm_gridlineIndices[atomIndexLocal * DIM + XX];
        const int iyBase = sm_gridlineIndices[atomIndexLocal * DIM + YY];
        int       iz     = sm_gridlineIndices[atomIndexLocal * DIM + ZZ] + ithz;
        if (iz >= nz)
        {
            iz -= nz;
        }
        /* Atom index w.r.t. warp - alternating 0 1 0 1 ... */
        const int atomWarpIndex = atomIndexLocal % atomsPerWarp;
        /* Warp index w.r.t. block - could probably be obtained easier? */
        const int warpIndex = atomIndexLocal / atomsPerWarp;

        const int splineIndexBase = getSplineParamIndexBase<order, atomsPerWarp>(warpIndex, atomWarpIndex);
        const int splineIndexZ = getSplineParamIndex<order, atomsPerWarp>(splineIndexBase, ZZ, ithz);
        const float thetaZ     = sm_theta[splineIndexZ];

        /* loop not used if order*order threads per atom */
        const int ithyMin = (threadsPerAtom == ThreadsPerAtom::Order) ? 0 : itemIdx.get_local_id(YY);
        const int ithyMax =
                (threadsPerAtom == ThreadsPerAtom::Order) ? order : itemIdx.get_local_id(YY) + 1;
        for (int ithy = ithyMin; ithy < ithyMax; ithy++)
        {
            int iy = iyBase + ithy;
            if (wrapY & (iy >= ny))
            {
                iy -= ny;
            }

            const int splineIndexY = getSplineParamIndex<order, atomsPerWarp>(splineIndexBase, YY, ithy);
            float       thetaY = sm_theta[splineIndexY];
            const float Val    = thetaZ * thetaY * (atomCharge);
            assertIsFinite(Val);
            const int offset = iy * pnz + iz;

#pragma unroll
            for (int ithx = 0; (ithx < order); ithx++)
            {
                int ix = ixBase + ithx;
                if (wrapX & (ix >= nx))
                {
                    ix -= nx;
                }
                const int gridIndexGlobal = ix * pny * pnz + offset;
                const int splineIndexX =
                        getSplineParamIndex<order, atomsPerWarp>(splineIndexBase, XX, ithx);
                const float thetaX = sm_theta[splineIndexX];
                assertIsFinite(thetaX);
                assertIsFinite(gm_grid[gridIndexGlobal]);
                atomicFetchAdd(gm_grid[gridIndexGlobal], thetaX * Val);
            }
        }
    }
}


/*! \brief
 * A spline computation and charge spreading kernel function.
 *
 * Two tuning parameters can be used for additional performance. For small systems and for debugging
 * writeGlobal should be used removing the need to recalculate the theta values in the gather kernel.
 * Similarly for large systems, with useOrderThreads, using order threads per atom gives higher
 * performance than order*order threads.
 *
 * \tparam order          PME interpolation order.
 * \tparam computeSplines A boolean which tells if the spline parameter and gridline indices'
 *                        computation should be performed.
 * \tparam spreadCharges  A boolean which tells if the charge spreading should be performed.
 * \tparam wrapX          A boolean which tells if the grid overlap in dimension X should be wrapped.
 * \tparam wrapY          A boolean which tells if the grid overlap in dimension Y should be wrapped.
 * \tparam numGrids       The number of grids to use in the kernel. Can be 1 or 2.
 * \tparam writeGlobal    A boolean which tells if the theta values and gridlines should be written
 *                        to global memory.
 * \tparam threadsPerAtom How many threads work on each atom.
 * \tparam subGroupSize   Size of the sub-group.
 */
template<int order, bool computeSplines, bool spreadCharges, bool wrapX, bool wrapY, int numGrids, bool writeGlobal, ThreadsPerAtom threadsPerAtom, int subGroupSize>
auto pmeSplineAndSpreadKernel(sycl::handler& cgh,
                              const int      nAtoms,
                              float* __restrict__ gm_realGrid_0,
                              float* __restrict__ gm_realGrid_1,
                              float* __restrict__ gm_theta,
                              float* __restrict__ gm_dtheta,
                              int* __restrict__ gm_gridlineIndices,
                              const float* __restrict__ gm_fractShiftsTable,
                              const int* __restrict__ gm_gridlineIndicesTable,
                              const float* __restrict__ gm_coefficients_0,
                              const float* __restrict__ gm_coefficients_1,
                              const Float3* __restrict__ gm_coordinates,
                              const gmx::IVec              tablesOffsets,
                              const gmx::IVec              realGridSize,
                              const gmx::RVec              realGridSizeFP,
                              const gmx::IVec              realGridSizePadded,
                              const gmx::RVec              currentRecipBox0,
                              const gmx::RVec              currentRecipBox1,
                              const gmx::RVec              currentRecipBox2,
                              const PmeGpuPipeliningParams pipeliningParams)
{
    constexpr int threadsPerAtomValue = (threadsPerAtom == ThreadsPerAtom::Order) ? order : order * order;
    constexpr int spreadMaxThreadsPerBlock = c_spreadMaxWarpsPerBlock * subGroupSize;
    constexpr int atomsPerBlock            = spreadMaxThreadsPerBlock / threadsPerAtomValue;
    // Number of atoms processed by a single warp in spread and gather
    static_assert(subGroupSize >= threadsPerAtomValue);
    constexpr int atomsPerWarp = subGroupSize / threadsPerAtomValue;

    // Gridline indices, ivec
    sycl::local_accessor<int, 1> sm_gridlineIndices(sycl::range<1>(atomsPerBlock * DIM), cgh);
    // Charges
    sycl::local_accessor<float, 1> sm_coefficients(sycl::range<1>(atomsPerBlock), cgh);
    // Spline values
    sycl::local_accessor<float, 1> sm_theta(sycl::range<1>(atomsPerBlock * DIM * order), cgh);
    auto                           sm_fractCoords = [&]() {
        if constexpr (computeSplines)
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
        const int blockIndex      = itemIdx.get_group_linear_id();
        const int atomIndexOffset = blockIndex * atomsPerBlock + pipeliningParams.pipelineAtomStart;

        /* Thread index w.r.t. block */
        const int threadLocalId = itemIdx.get_local_linear_id();
        /* Warp index w.r.t. block - could probably be obtained easier? */
        const int warpIndex = threadLocalId / subGroupSize;

        /* Atom index w.r.t. warp */
        const int atomWarpIndex = itemIdx.get_local_id(XX) % atomsPerWarp;
        /* Atom index w.r.t. block/shared memory */
        const int atomIndexLocal = warpIndex * atomsPerWarp + atomWarpIndex;
        /* Atom index w.r.t. global memory */
        const int atomIndexGlobal = atomIndexOffset + atomIndexLocal;

        /* Early return for fully empty blocks at the end
         * (should only happen for billions of input atoms) */
        if (atomIndexOffset >= nAtoms)
        {
            return;
        }

        /* Charges, required for both spline and spread */
        pmeGpuStageAtomData<float, atomsPerBlock, 1>(sm_coefficients.get_pointer(),
                                                     gm_coefficients_0 + pipeliningParams.pipelineAtomStart,
                                                     itemIdx);
        itemIdx.barrier(fence_space::local_space);
        const float atomCharge = sm_coefficients[atomIndexLocal];

        if constexpr (computeSplines)
        {
            // SYCL-TODO: Use prefetching? Issue #4153.
            const Float3 atomX = gm_coordinates[atomIndexGlobal];
            calculateSplines<order, atomsPerBlock, atomsPerWarp, false, writeGlobal, numGrids, subGroupSize>(
                    atomIndexOffset,
                    atomX,
                    atomCharge,
                    tablesOffsets,
                    realGridSizeFP,
                    currentRecipBox0,
                    currentRecipBox1,
                    currentRecipBox2,
                    gm_theta,
                    gm_dtheta,
                    gm_gridlineIndices,
                    gm_fractShiftsTable,
                    gm_gridlineIndicesTable,
                    sm_theta.get_pointer(),
                    nullptr,
                    sm_gridlineIndices.get_pointer(),
                    sm_fractCoords.get_pointer(),
                    itemIdx);
            subGroupBarrier(itemIdx);
        }
        else
        {
            /* Staging the data for spread
             * (the data is assumed to be in GPU global memory with proper layout already,
             * as in after running the spline kernel)
             */
            /* Spline data - only thetas (dthetas will only be needed in gather) */
            pmeGpuStageAtomData<float, atomsPerBlock, DIM * order>(sm_theta.get_pointer(), gm_theta, itemIdx);
            /* Gridline indices */
            pmeGpuStageAtomData<int, atomsPerBlock, DIM>(
                    sm_gridlineIndices.get_pointer(), gm_gridlineIndices, itemIdx);

            itemIdx.barrier(fence_space::local_space);
        }

        /* Spreading */
        if (spreadCharges && atomIndexGlobal < nAtoms)
        {
            if (!pipeliningParams.usePipeline || (atomIndexGlobal < pipeliningParams.pipelineAtomEnd))
            {
                spread_charges<order, wrapX, wrapY, threadsPerAtom, subGroupSize>(
                        atomCharge,
                        realGridSize,
                        realGridSizePadded,
                        gm_realGrid_0,
                        sm_gridlineIndices.get_pointer(),
                        sm_theta.get_pointer(),
                        itemIdx);
            }
        }
        if constexpr (numGrids == 2 && spreadCharges)
        {
            itemIdx.barrier(fence_space::local_space);
            pmeGpuStageAtomData<float, atomsPerBlock, 1>(sm_coefficients.get_pointer(),
                                                         gm_coefficients_1 + pipeliningParams.pipelineAtomStart,
                                                         itemIdx);
            itemIdx.barrier(fence_space::local_space);
            const float atomCharge = sm_coefficients[atomIndexLocal];
            if (atomIndexGlobal < nAtoms)
            {
                if (!pipeliningParams.usePipeline || (atomIndexGlobal < pipeliningParams.pipelineAtomEnd))
                {
                    spread_charges<order, wrapX, wrapY, threadsPerAtom, subGroupSize>(
                            atomCharge,
                            realGridSize,
                            realGridSizePadded,
                            gm_realGrid_1,
                            sm_gridlineIndices.get_pointer(),
                            sm_theta.get_pointer(),
                            itemIdx);
                }
            }
        }
    };
}

template<int order, bool computeSplines, bool spreadCharges, bool wrapX, bool wrapY, int numGrids, bool writeGlobal, ThreadsPerAtom threadsPerAtom, int subGroupSize>
PmeSplineAndSpreadKernel<order, computeSplines, spreadCharges, wrapX, wrapY, numGrids, writeGlobal, threadsPerAtom, subGroupSize>::PmeSplineAndSpreadKernel()
{
    reset();
}

template<int order, bool computeSplines, bool spreadCharges, bool wrapX, bool wrapY, int numGrids, bool writeGlobal, ThreadsPerAtom threadsPerAtom, int subGroupSize>
void PmeSplineAndSpreadKernel<order, computeSplines, spreadCharges, wrapX, wrapY, numGrids, writeGlobal, threadsPerAtom, subGroupSize>::setArg(
        size_t argIndex,
        void*  arg)
{
    if (argIndex == 0)
    {
        auto* params   = reinterpret_cast<PmeGpuKernelParams*>(arg);
        gridParams_    = &params->grid;
        atomParams_    = &params->atoms;
        dynamicParams_ = &params->current;
        pipeliningParams_ = { params->pipelineAtomStart, params->pipelineAtomEnd, params->usePipeline != 0 };
    }
    else
    {
        GMX_RELEASE_ASSERT(argIndex == 0, "Trying to pass too many args to the kernel");
    }
}


template<int order, bool computeSplines, bool spreadCharges, bool wrapX, bool wrapY, int numGrids, bool writeGlobal, ThreadsPerAtom threadsPerAtom, int subGroupSize>
void PmeSplineAndSpreadKernel<order, computeSplines, spreadCharges, wrapX, wrapY, numGrids, writeGlobal, threadsPerAtom, subGroupSize>::launch(
        const KernelLaunchConfig& config,
        const DeviceStream&       deviceStream)
{
    GMX_RELEASE_ASSERT(gridParams_, "Can not launch the kernel before setting its args");
    GMX_RELEASE_ASSERT(atomParams_, "Can not launch the kernel before setting its args");
    GMX_RELEASE_ASSERT(dynamicParams_, "Can not launch the kernel before setting its args");

    using kernelNameType =
            PmeSplineAndSpreadKernel<order, computeSplines, spreadCharges, wrapX, wrapY, numGrids, writeGlobal, threadsPerAtom, subGroupSize>;

    // SYCL has different multidimensional layout than OpenCL/CUDA.
    const sycl::range<3> localSize{ config.blockSize[2], config.blockSize[1], config.blockSize[0] };
    const sycl::range<3> groupRange{ config.gridSize[2], config.gridSize[1], config.gridSize[0] };
    const sycl::nd_range<3> range{ groupRange * localSize, localSize };

    sycl::queue q = deviceStream.stream();

    q.submit(GMX_SYCL_DISCARD_EVENT[&](sycl::handler & cgh) {
        auto kernel =
                pmeSplineAndSpreadKernel<order, computeSplines, spreadCharges, wrapX, wrapY, numGrids, writeGlobal, threadsPerAtom, subGroupSize>(
                        cgh,
                        atomParams_->nAtoms,
                        gridParams_->d_realGrid[0].get_pointer(),
                        gridParams_->d_realGrid[1].get_pointer(),
                        atomParams_->d_theta.get_pointer(),
                        atomParams_->d_dtheta.get_pointer(),
                        atomParams_->d_gridlineIndices.get_pointer(),
                        gridParams_->d_fractShiftsTable.get_pointer(),
                        gridParams_->d_gridlineIndicesTable.get_pointer(),
                        atomParams_->d_coefficients[0].get_pointer(),
                        atomParams_->d_coefficients[1].get_pointer(),
                        atomParams_->d_coordinates.get_pointer(),
                        gridParams_->tablesOffsets,
                        gridParams_->realGridSize,
                        gridParams_->realGridSizeFP,
                        gridParams_->realGridSizePadded,
                        dynamicParams_->recipBox[0],
                        dynamicParams_->recipBox[1],
                        dynamicParams_->recipBox[2],
                        pipeliningParams_);
        cgh.parallel_for<kernelNameType>(range, kernel);
    });

    // Delete set args, so we don't forget to set them before the next launch.
    reset();
}


template<int order, bool computeSplines, bool spreadCharges, bool wrapX, bool wrapY, int numGrids, bool writeGlobal, ThreadsPerAtom threadsPerAtom, int subGroupSize>
void PmeSplineAndSpreadKernel<order, computeSplines, spreadCharges, wrapX, wrapY, numGrids, writeGlobal, threadsPerAtom, subGroupSize>::reset()
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

#define INSTANTIATE_3(order, computeSplines, spreadCharges, numGrids, writeGlobal, threadsPerAtom, subGroupSize) \
    template class PmeSplineAndSpreadKernel<order, computeSplines, spreadCharges, true, true, numGrids, writeGlobal, threadsPerAtom, subGroupSize>;

#define INSTANTIATE_2(order, numGrids, threadsPerAtom, subGroupSize)                 \
    INSTANTIATE_3(order, true, true, numGrids, true, threadsPerAtom, subGroupSize);  \
    INSTANTIATE_3(order, true, false, numGrids, true, threadsPerAtom, subGroupSize); \
    INSTANTIATE_3(order, false, true, numGrids, true, threadsPerAtom, subGroupSize); \
    INSTANTIATE_3(order, true, true, numGrids, false, threadsPerAtom, subGroupSize);

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
