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
 * \brief
 * Implements PmeGpuProgramImpl, which stores permanent PME GPU context-derived data,
 * such as (compiled) kernel handles.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \author Andrey Alekseenko <al42and@gmail.com>
 * \ingroup module_ewald
 */
#include "gmxpre.h"

#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/gpu_utils/syclutils.h"
#include "gromacs/hardware/device_information.h"

#include "pme_gather_sycl.h"
#include "pme_gpu_constants.h"
#include "pme_gpu_internal.h" // for GridOrdering enum
#include "pme_gpu_program_impl.h"
#include "pme_gpu_types_host.h"
#include "pme_solve_sycl.h"
#include "pme_spread_sycl.h"

// PME interpolation order
constexpr int c_pmeOrder = 4;
// These hardcoded spread/gather parameters refer to not-implemented PME GPU 2D decomposition in X/Y
constexpr bool c_wrapX = true;
constexpr bool c_wrapY = true;

constexpr int c_stateA = 0;
constexpr int c_stateB = 1;

static int chooseSubGroupSizeForDevice(const DeviceInformation& deviceInfo)
{
    if (deviceInfo.supportedSubGroupSizes.size() == 1)
    {
        return deviceInfo.supportedSubGroupSizes[0];
    }
    else if (deviceInfo.supportedSubGroupSizes.size() > 1)
    {
        switch (deviceInfo.deviceVendor)
        {
            case DeviceVendor::Intel: return 16; // TODO: Choose best value, Issue #4153.
            default:
                GMX_RELEASE_ASSERT(false, "Flexible sub-groups only supported for Intel GPUs");
                return 0;
        }
    }
    else
    {
        GMX_RELEASE_ASSERT(false, "Device has no known supported sub-group sizes");
        return 0;
    }
}

#define INSTANTIATE_SPREAD_2(                                                                      \
        order, computeSplines, spreadCharges, numGrids, writeGlobal, threadsPerAtom, subGroupSize) \
    extern template class PmeSplineAndSpreadKernel<order, computeSplines, spreadCharges, true, true, numGrids, writeGlobal, threadsPerAtom, subGroupSize>;

#define INSTANTIATE_SPREAD(order, numGrids, threadsPerAtom, subGroupSize)                   \
    INSTANTIATE_SPREAD_2(order, true, true, numGrids, true, threadsPerAtom, subGroupSize);  \
    INSTANTIATE_SPREAD_2(order, true, false, numGrids, true, threadsPerAtom, subGroupSize); \
    INSTANTIATE_SPREAD_2(order, false, true, numGrids, true, threadsPerAtom, subGroupSize); \
    INSTANTIATE_SPREAD_2(order, true, true, numGrids, false, threadsPerAtom, subGroupSize);

#define INSTANTIATE_GATHER_2(order, numGrids, readGlobal, threadsPerAtom, subGroupSize) \
    extern template class PmeGatherKernel<order, true, true, numGrids, readGlobal, threadsPerAtom, subGroupSize>;

#define INSTANTIATE_GATHER(order, numGrids, threadsPerAtom, subGroupSize)      \
    INSTANTIATE_GATHER_2(order, numGrids, true, threadsPerAtom, subGroupSize); \
    INSTANTIATE_GATHER_2(order, numGrids, false, threadsPerAtom, subGroupSize);

#define INSTANTIATE_X(x, order, subGroupSize)                              \
    INSTANTIATE_##x(order, 1, ThreadsPerAtom::Order, subGroupSize);        \
    INSTANTIATE_##x(order, 1, ThreadsPerAtom::OrderSquared, subGroupSize); \
    INSTANTIATE_##x(order, 2, ThreadsPerAtom::Order, subGroupSize);        \
    INSTANTIATE_##x(order, 2, ThreadsPerAtom::OrderSquared, subGroupSize);

#define INSTANTIATE_SOLVE(subGroupSize)                                                     \
    extern template class PmeSolveKernel<GridOrdering::XYZ, false, c_stateA, subGroupSize>; \
    extern template class PmeSolveKernel<GridOrdering::XYZ, true, c_stateA, subGroupSize>;  \
    extern template class PmeSolveKernel<GridOrdering::YZX, false, c_stateA, subGroupSize>; \
    extern template class PmeSolveKernel<GridOrdering::YZX, true, c_stateA, subGroupSize>;  \
    extern template class PmeSolveKernel<GridOrdering::XYZ, false, c_stateB, subGroupSize>; \
    extern template class PmeSolveKernel<GridOrdering::XYZ, true, c_stateB, subGroupSize>;  \
    extern template class PmeSolveKernel<GridOrdering::YZX, false, c_stateB, subGroupSize>; \
    extern template class PmeSolveKernel<GridOrdering::YZX, true, c_stateB, subGroupSize>;

#define INSTANTIATE(order, subGroupSize)        \
    INSTANTIATE_X(SPREAD, order, subGroupSize); \
    INSTANTIATE_X(GATHER, order, subGroupSize); \
    INSTANTIATE_SOLVE(subGroupSize);

#if GMX_SYCL_DPCPP
INSTANTIATE(4, 16);
#endif
INSTANTIATE(4, 32);
INSTANTIATE(4, 64);

//! Helper function to set proper kernel functor pointers
template<int subGroupSize>
static void setKernelPointers(struct PmeGpuProgramImpl* pmeGpuProgram)
{
    /* Not all combinations of the splineAndSpread, spline and Spread kernels are required
     * If only the spline (without the spread) then it does not make sense not to write the data to global memory
     * Similarly the spread kernel (without the spline) implies that we should read the spline data from global memory
     */
    pmeGpuProgram->splineAndSpreadKernelSingle =
            new PmeSplineAndSpreadKernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 1, false, ThreadsPerAtom::OrderSquared, subGroupSize>();
    pmeGpuProgram->splineAndSpreadKernelThPerAtom4Single =
            new PmeSplineAndSpreadKernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 1, false, ThreadsPerAtom::Order, subGroupSize>();
    pmeGpuProgram->splineAndSpreadKernelWriteSplinesSingle =
            new PmeSplineAndSpreadKernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::OrderSquared, subGroupSize>();
    pmeGpuProgram->splineAndSpreadKernelWriteSplinesThPerAtom4Single =
            new PmeSplineAndSpreadKernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::Order, subGroupSize>();
    pmeGpuProgram->splineKernelSingle =
            new PmeSplineAndSpreadKernel<c_pmeOrder, true, false, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::OrderSquared, subGroupSize>();
    pmeGpuProgram->splineKernelThPerAtom4Single =
            new PmeSplineAndSpreadKernel<c_pmeOrder, true, false, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::Order, subGroupSize>();
    pmeGpuProgram->spreadKernelSingle =
            new PmeSplineAndSpreadKernel<c_pmeOrder, false, true, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::OrderSquared, subGroupSize>();
    pmeGpuProgram->spreadKernelThPerAtom4Single =
            new PmeSplineAndSpreadKernel<c_pmeOrder, false, true, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::Order, subGroupSize>();
    pmeGpuProgram->splineAndSpreadKernelDual =
            new PmeSplineAndSpreadKernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 2, false, ThreadsPerAtom::OrderSquared, subGroupSize>();
    pmeGpuProgram->splineAndSpreadKernelThPerAtom4Dual =
            new PmeSplineAndSpreadKernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 2, false, ThreadsPerAtom::Order, subGroupSize>();
    pmeGpuProgram->splineAndSpreadKernelWriteSplinesDual =
            new PmeSplineAndSpreadKernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::OrderSquared, subGroupSize>();
    pmeGpuProgram->splineAndSpreadKernelWriteSplinesThPerAtom4Dual =
            new PmeSplineAndSpreadKernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::Order, subGroupSize>();
    pmeGpuProgram->splineKernelDual =
            new PmeSplineAndSpreadKernel<c_pmeOrder, true, false, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::OrderSquared, subGroupSize>();
    pmeGpuProgram->splineKernelThPerAtom4Dual =
            new PmeSplineAndSpreadKernel<c_pmeOrder, true, false, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::Order, subGroupSize>();
    pmeGpuProgram->spreadKernelDual =
            new PmeSplineAndSpreadKernel<c_pmeOrder, false, true, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::OrderSquared, subGroupSize>();
    pmeGpuProgram->spreadKernelThPerAtom4Dual =
            new PmeSplineAndSpreadKernel<c_pmeOrder, false, true, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::Order, subGroupSize>();
    pmeGpuProgram->gatherKernelSingle =
            new PmeGatherKernel<c_pmeOrder, c_wrapX, c_wrapY, 1, false, ThreadsPerAtom::OrderSquared, subGroupSize>();
    pmeGpuProgram->gatherKernelThPerAtom4Single =
            new PmeGatherKernel<c_pmeOrder, c_wrapX, c_wrapY, 1, false, ThreadsPerAtom::Order, subGroupSize>();
    pmeGpuProgram->gatherKernelReadSplinesSingle =
            new PmeGatherKernel<c_pmeOrder, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::OrderSquared, subGroupSize>();
    pmeGpuProgram->gatherKernelReadSplinesThPerAtom4Single =
            new PmeGatherKernel<c_pmeOrder, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::Order, subGroupSize>();
    pmeGpuProgram->gatherKernelDual =
            new PmeGatherKernel<c_pmeOrder, c_wrapX, c_wrapY, 2, false, ThreadsPerAtom::OrderSquared, subGroupSize>();
    pmeGpuProgram->gatherKernelThPerAtom4Dual =
            new PmeGatherKernel<c_pmeOrder, c_wrapX, c_wrapY, 2, false, ThreadsPerAtom::Order, subGroupSize>();
    pmeGpuProgram->gatherKernelReadSplinesDual =
            new PmeGatherKernel<c_pmeOrder, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::OrderSquared, subGroupSize>();
    pmeGpuProgram->gatherKernelReadSplinesThPerAtom4Dual =
            new PmeGatherKernel<c_pmeOrder, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::Order, subGroupSize>();
    pmeGpuProgram->solveXYZKernelA =
            new PmeSolveKernel<GridOrdering::XYZ, false, c_stateA, subGroupSize>();
    pmeGpuProgram->solveXYZEnergyKernelA =
            new PmeSolveKernel<GridOrdering::XYZ, true, c_stateA, subGroupSize>();
    pmeGpuProgram->solveYZXKernelA =
            new PmeSolveKernel<GridOrdering::YZX, false, c_stateA, subGroupSize>();
    pmeGpuProgram->solveYZXEnergyKernelA =
            new PmeSolveKernel<GridOrdering::YZX, true, c_stateA, subGroupSize>();
    pmeGpuProgram->solveXYZKernelB =
            new PmeSolveKernel<GridOrdering::XYZ, false, c_stateB, subGroupSize>();
    pmeGpuProgram->solveXYZEnergyKernelB =
            new PmeSolveKernel<GridOrdering::XYZ, true, c_stateB, subGroupSize>();
    pmeGpuProgram->solveYZXKernelB =
            new PmeSolveKernel<GridOrdering::YZX, false, c_stateB, subGroupSize>();
    pmeGpuProgram->solveYZXEnergyKernelB =
            new PmeSolveKernel<GridOrdering::YZX, true, c_stateB, subGroupSize>();
}

PmeGpuProgramImpl::PmeGpuProgramImpl(const DeviceContext& deviceContext) :
    deviceContext_(deviceContext)
{
    // kernel parameters
    const DeviceInformation& deviceInfo = deviceContext.deviceInfo();
    warpSize_                           = chooseSubGroupSizeForDevice(deviceInfo);
    GMX_RELEASE_ASSERT(std::find(deviceInfo.supportedSubGroupSizes.begin(),
                                 deviceInfo.supportedSubGroupSizes.end(),
                                 warpSize_)
                               != deviceInfo.supportedSubGroupSizes.end(),
                       "Device does not support selected sub-group size");
    spreadWorkGroupSize   = c_spreadMaxWarpsPerBlock * warpSize_;
    solveMaxWorkGroupSize = c_solveMaxWarpsPerBlock * warpSize_;
    gatherWorkGroupSize   = c_gatherMaxWarpsPerBlock * warpSize_;

    switch (warpSize_)
    {
#if GMX_SYCL_DPCPP
        case 16: setKernelPointers<16>(this); break;
#endif
        case 32: setKernelPointers<32>(this); break;
        case 64: setKernelPointers<64>(this); break;
        default: GMX_RELEASE_ASSERT(false, "Invalid sub group size");
    }
}

PmeGpuProgramImpl::~PmeGpuProgramImpl()
{
    delete splineKernelSingle;
    delete splineKernelThPerAtom4Single;
    delete spreadKernelSingle;
    delete spreadKernelThPerAtom4Single;
    delete splineAndSpreadKernelSingle;
    delete splineAndSpreadKernelThPerAtom4Single;
    delete splineAndSpreadKernelWriteSplinesSingle;
    delete splineAndSpreadKernelWriteSplinesThPerAtom4Single;
    delete splineKernelDual;
    delete splineKernelThPerAtom4Dual;
    delete spreadKernelDual;
    delete spreadKernelThPerAtom4Dual;
    delete splineAndSpreadKernelDual;
    delete splineAndSpreadKernelThPerAtom4Dual;
    delete splineAndSpreadKernelWriteSplinesDual;
    delete splineAndSpreadKernelWriteSplinesThPerAtom4Dual;
    delete gatherKernelSingle;
    delete gatherKernelThPerAtom4Single;
    delete gatherKernelReadSplinesSingle;
    delete gatherKernelReadSplinesThPerAtom4Single;
    delete gatherKernelDual;
    delete gatherKernelThPerAtom4Dual;
    delete gatherKernelReadSplinesDual;
    delete gatherKernelReadSplinesThPerAtom4Dual;
    delete solveYZXKernelA;
    delete solveXYZKernelA;
    delete solveYZXEnergyKernelA;
    delete solveXYZEnergyKernelA;
    delete solveYZXKernelB;
    delete solveXYZKernelB;
    delete solveYZXEnergyKernelB;
    delete solveXYZEnergyKernelB;
}
