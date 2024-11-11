/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * \brief
 * Implements PmeGpuProgramImpl, which stores permanent PME GPU context-derived data,
 * such as (compiled) kernel handles.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \author Andrey Alekseenko <al42and@gmail.com>
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_ewald
 */
#include "gmxpre.h"

#include "gromacs/hardware/device_information.h"

#include "pme_gpu_constants.h"
#include "pme_gpu_internal.h" // for GridOrdering enum
#include "pme_gpu_internal.h"
#include "pme_gpu_program_impl.h"
#include "pme_gpu_types_host.h"

// PME interpolation order
constexpr int c_pmeOrder = 4;
// These hardcoded spread/gather parameters refer to not-implemented PME GPU 2D decomposition in X/Y
constexpr bool c_wrapX = true;
constexpr bool c_wrapY = true;

constexpr int c_stateA = 0;
constexpr int c_stateB = 1;

template<int parallelExecutionWidth>
constexpr int sc_spreadHipMaxWarpsPerBlock = (parallelExecutionWidth == 64) ? 8 : 4;

template<int parallelExecutionWidth>
static constexpr int sc_spreadMaxThreadsPerBlock =
        sc_spreadHipMaxWarpsPerBlock<parallelExecutionWidth> * parallelExecutionWidth;

template<int parallelExecutionWidth>
static constexpr int sc_solveMaxThreadsPerBlock = c_solveMaxWarpsPerBlock * parallelExecutionWidth;

template<int parallelExecutionWidth>
static constexpr int sc_gatherMaxThreadsPerBlock = c_gatherMaxWarpsPerBlock * parallelExecutionWidth;

static int deviceParallelExecutionSize(const DeviceInformation& deviceInfo)
{
    return deviceInfo.supportedSubGroupSizes[0];
}

//! PME HIP kernels forward declarations. Kernels are documented in their respective files.
template<int order, bool computeSplines, bool spreadCharges, bool wrapX, bool wrapY, int mode, bool writeGlobal, ThreadsPerAtom threadsPerAtom, int parallelExecutionWidth>
__global__ void pmeSplineAndSpreadKernel(PmeGpuKernelParamsBase kernelParams);

template<GridOrdering gridOrdering, bool computeEnergyAndVirial, const int gridIndex, int parallelExecutionWidth> /* It is significantly slower to pass gridIndex as a kernel parameter */
__global__ void pmeSolveKernel(PmeGpuKernelParamsBase kernelParams);

template<int order, bool wrapX, bool wrapY, int nGrids, bool readGlobal, ThreadsPerAtom threadsPerAtom, int parallelExecutionWidth>
__global__ void pmeGatherKernel(PmeGpuKernelParamsBase kernelParams);


#define INSTANTIATE_SPREAD_2(                                                                                                                  \
        order, computeSplines, spreadCharges, numGrids, writeGlobal, threadsPerAtom, parallelExecutionWidth)                                   \
    extern template __global__ void                                                                                                            \
    pmeSplineAndSpreadKernel<order, computeSplines, spreadCharges, true, true, numGrids, writeGlobal, threadsPerAtom, parallelExecutionWidth>( \
            PmeGpuKernelParamsBase kernelParams);

#define INSTANTIATE_SPREAD(order, numGrids, threadsPerAtom, parallelExecutionWidth)                   \
    INSTANTIATE_SPREAD_2(order, true, true, numGrids, true, threadsPerAtom, parallelExecutionWidth);  \
    INSTANTIATE_SPREAD_2(order, true, false, numGrids, true, threadsPerAtom, parallelExecutionWidth); \
    INSTANTIATE_SPREAD_2(order, false, true, numGrids, true, threadsPerAtom, parallelExecutionWidth); \
    INSTANTIATE_SPREAD_2(order, true, true, numGrids, false, threadsPerAtom, parallelExecutionWidth);

#define INSTANTIATE_GATHER_2(order, numGrids, readGlobal, threadsPerAtom, parallelExecutionWidth)                                     \
    extern template __global__ void pmeGatherKernel<order, true, true, numGrids, readGlobal, threadsPerAtom, parallelExecutionWidth>( \
            PmeGpuKernelParamsBase kernelParams);

#define INSTANTIATE_GATHER(order, numGrids, threadsPerAtom, parallelExecutionWidth)      \
    INSTANTIATE_GATHER_2(order, numGrids, true, threadsPerAtom, parallelExecutionWidth); \
    INSTANTIATE_GATHER_2(order, numGrids, false, threadsPerAtom, parallelExecutionWidth);

#define INSTANTIATE_X(x, order, parallelExecutionWidth)                              \
    INSTANTIATE_##x(order, 1, ThreadsPerAtom::Order, parallelExecutionWidth);        \
    INSTANTIATE_##x(order, 1, ThreadsPerAtom::OrderSquared, parallelExecutionWidth); \
    INSTANTIATE_##x(order, 2, ThreadsPerAtom::Order, parallelExecutionWidth);        \
    INSTANTIATE_##x(order, 2, ThreadsPerAtom::OrderSquared, parallelExecutionWidth);

#define INSTANTIATE_SOLVE(parallelExecutionWidth)                                                               \
    extern template __global__ void pmeSolveKernel<GridOrdering::XYZ, false, c_stateA, parallelExecutionWidth>( \
            PmeGpuKernelParamsBase kernelParams);                                                               \
    extern template __global__ void pmeSolveKernel<GridOrdering::XYZ, true, c_stateA, parallelExecutionWidth>(  \
            PmeGpuKernelParamsBase kernelParams);                                                               \
    extern template __global__ void pmeSolveKernel<GridOrdering::YZX, false, c_stateA, parallelExecutionWidth>( \
            PmeGpuKernelParamsBase kernelParams);                                                               \
    extern template __global__ void pmeSolveKernel<GridOrdering::YZX, true, c_stateA, parallelExecutionWidth>(  \
            PmeGpuKernelParamsBase kernelParams);                                                               \
    extern template __global__ void pmeSolveKernel<GridOrdering::XYZ, false, c_stateB, parallelExecutionWidth>( \
            PmeGpuKernelParamsBase kernelParams);                                                               \
    extern template __global__ void pmeSolveKernel<GridOrdering::XYZ, true, c_stateB, parallelExecutionWidth>(  \
            PmeGpuKernelParamsBase kernelParams);                                                               \
    extern template __global__ void pmeSolveKernel<GridOrdering::YZX, false, c_stateB, parallelExecutionWidth>( \
            PmeGpuKernelParamsBase kernelParams);                                                               \
    extern template __global__ void pmeSolveKernel<GridOrdering::YZX, true, c_stateB, parallelExecutionWidth>(  \
            PmeGpuKernelParamsBase kernelParams);

#define INSTANTIATE(order, parallelExecutionWidth)        \
    INSTANTIATE_X(SPREAD, order, parallelExecutionWidth); \
    INSTANTIATE_X(GATHER, order, parallelExecutionWidth); \
    INSTANTIATE_SOLVE(parallelExecutionWidth);

INSTANTIATE(4, 32);
INSTANTIATE(4, 64);

//! Helper function to set proper kernel functor pointers
template<int parallelExecutionWidth>
static void setKernelPointersAndParams(struct PmeGpuProgramImpl* pmeGpuProgram)
{
    /* Not all combinations of the splineAndSpread, spline and Spread kernels are required
     * If only the spline (without the spread) then it does not make sense not to write the data to global memory
     * Similarly the spread kernel (without the spline) implies that we should read the spline data from global memory
     */
    pmeGpuProgram->splineAndSpreadKernelSingle =
            pmeSplineAndSpreadKernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 1, false, ThreadsPerAtom::OrderSquared, parallelExecutionWidth>;
    pmeGpuProgram->splineAndSpreadKernelThPerAtom4Single =
            pmeSplineAndSpreadKernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 1, false, ThreadsPerAtom::Order, parallelExecutionWidth>;
    pmeGpuProgram->splineAndSpreadKernelWriteSplinesSingle =
            pmeSplineAndSpreadKernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::OrderSquared, parallelExecutionWidth>;
    pmeGpuProgram->splineAndSpreadKernelWriteSplinesThPerAtom4Single =
            pmeSplineAndSpreadKernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::Order, parallelExecutionWidth>;
    pmeGpuProgram->splineKernelSingle =
            pmeSplineAndSpreadKernel<c_pmeOrder, true, false, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::OrderSquared, parallelExecutionWidth>;
    pmeGpuProgram->splineKernelThPerAtom4Single =
            pmeSplineAndSpreadKernel<c_pmeOrder, true, false, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::Order, parallelExecutionWidth>;
    pmeGpuProgram->spreadKernelSingle =
            pmeSplineAndSpreadKernel<c_pmeOrder, false, true, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::OrderSquared, parallelExecutionWidth>;
    pmeGpuProgram->spreadKernelThPerAtom4Single =
            pmeSplineAndSpreadKernel<c_pmeOrder, false, true, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::Order, parallelExecutionWidth>;
    pmeGpuProgram->splineAndSpreadKernelDual =
            pmeSplineAndSpreadKernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 2, false, ThreadsPerAtom::OrderSquared, parallelExecutionWidth>;
    pmeGpuProgram->splineAndSpreadKernelThPerAtom4Dual =
            pmeSplineAndSpreadKernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 2, false, ThreadsPerAtom::Order, parallelExecutionWidth>;
    pmeGpuProgram->splineAndSpreadKernelWriteSplinesDual =
            pmeSplineAndSpreadKernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::OrderSquared, parallelExecutionWidth>;
    pmeGpuProgram->splineAndSpreadKernelWriteSplinesThPerAtom4Dual =
            pmeSplineAndSpreadKernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::Order, parallelExecutionWidth>;
    pmeGpuProgram->splineKernelDual =
            pmeSplineAndSpreadKernel<c_pmeOrder, true, false, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::OrderSquared, parallelExecutionWidth>;
    pmeGpuProgram->splineKernelThPerAtom4Dual =
            pmeSplineAndSpreadKernel<c_pmeOrder, true, false, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::Order, parallelExecutionWidth>;
    pmeGpuProgram->spreadKernelDual =
            pmeSplineAndSpreadKernel<c_pmeOrder, false, true, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::OrderSquared, parallelExecutionWidth>;
    pmeGpuProgram->spreadKernelThPerAtom4Dual =
            pmeSplineAndSpreadKernel<c_pmeOrder, false, true, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::Order, parallelExecutionWidth>;
    pmeGpuProgram->gatherKernelSingle =
            pmeGatherKernel<c_pmeOrder, c_wrapX, c_wrapY, 1, false, ThreadsPerAtom::OrderSquared, parallelExecutionWidth>;
    pmeGpuProgram->gatherKernelThPerAtom4Single =
            pmeGatherKernel<c_pmeOrder, c_wrapX, c_wrapY, 1, false, ThreadsPerAtom::Order, parallelExecutionWidth>;
    pmeGpuProgram->gatherKernelReadSplinesSingle =
            pmeGatherKernel<c_pmeOrder, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::OrderSquared, parallelExecutionWidth>;
    pmeGpuProgram->gatherKernelReadSplinesThPerAtom4Single =
            pmeGatherKernel<c_pmeOrder, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::Order, parallelExecutionWidth>;
    pmeGpuProgram->gatherKernelDual =
            pmeGatherKernel<c_pmeOrder, c_wrapX, c_wrapY, 2, false, ThreadsPerAtom::OrderSquared, parallelExecutionWidth>;
    pmeGpuProgram->gatherKernelThPerAtom4Dual =
            pmeGatherKernel<c_pmeOrder, c_wrapX, c_wrapY, 2, false, ThreadsPerAtom::Order, parallelExecutionWidth>;
    pmeGpuProgram->gatherKernelReadSplinesDual =
            pmeGatherKernel<c_pmeOrder, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::OrderSquared, parallelExecutionWidth>;
    pmeGpuProgram->gatherKernelReadSplinesThPerAtom4Dual =
            pmeGatherKernel<c_pmeOrder, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::Order, parallelExecutionWidth>;
    pmeGpuProgram->solveXYZKernelA =
            pmeSolveKernel<GridOrdering::XYZ, false, c_stateA, parallelExecutionWidth>;
    pmeGpuProgram->solveXYZEnergyKernelA =
            pmeSolveKernel<GridOrdering::XYZ, true, c_stateA, parallelExecutionWidth>;
    pmeGpuProgram->solveYZXKernelA =
            pmeSolveKernel<GridOrdering::YZX, false, c_stateA, parallelExecutionWidth>;
    pmeGpuProgram->solveYZXEnergyKernelA =
            pmeSolveKernel<GridOrdering::YZX, true, c_stateA, parallelExecutionWidth>;
    pmeGpuProgram->solveXYZKernelB =
            pmeSolveKernel<GridOrdering::XYZ, false, c_stateB, parallelExecutionWidth>;
    pmeGpuProgram->solveXYZEnergyKernelB =
            pmeSolveKernel<GridOrdering::XYZ, true, c_stateB, parallelExecutionWidth>;
    pmeGpuProgram->solveYZXKernelB =
            pmeSolveKernel<GridOrdering::YZX, false, c_stateB, parallelExecutionWidth>;
    pmeGpuProgram->solveYZXEnergyKernelB =
            pmeSolveKernel<GridOrdering::YZX, true, c_stateB, parallelExecutionWidth>;

    pmeGpuProgram->spreadWorkGroupSize   = sc_spreadMaxThreadsPerBlock<parallelExecutionWidth>;
    pmeGpuProgram->solveMaxWorkGroupSize = sc_solveMaxThreadsPerBlock<parallelExecutionWidth>;
    pmeGpuProgram->gatherWorkGroupSize   = sc_gatherMaxThreadsPerBlock<parallelExecutionWidth>;
}

PmeGpuProgramImpl::PmeGpuProgramImpl(const DeviceContext& deviceContext) :
    deviceContext_(deviceContext)
{
    // kernel parameters
    const DeviceInformation& deviceInfo = deviceContext.deviceInfo();
    warpSize_                           = deviceParallelExecutionSize(deviceInfo);

    switch (warpSize_)
    {
        case 32: setKernelPointersAndParams<32>(this); break;
        case 64: setKernelPointersAndParams<64>(this); break;
        default: GMX_RELEASE_ASSERT(false, "Invalid sub group size");
    }
}

PmeGpuProgramImpl::~PmeGpuProgramImpl() = default;
