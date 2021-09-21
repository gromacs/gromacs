/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2021, by the GROMACS development team, led by
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
 * \brief
 * Implements PmeGpuProgramImpl, which stores permanent PME GPU context-derived data,
 * such as (compiled) kernel handles.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \author Andrey Alekseenko <al42and@gmail.com>
 * \ingroup module_ewald
 */
#include "gmxpre.h"

#include "gromacs/hardware/device_information.h"
#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/gpu_utils/syclutils.h"

#include "pme_gpu_program_impl.h"
#include "pme_spread_sycl.h"

#include "pme_gpu_constants.h"
#include "pme_gpu_internal.h" // for GridOrdering enum
#include "pme_gpu_types_host.h"

// PME interpolation order
constexpr int c_pmeOrder = 4;
// These hardcoded spread/gather parameters refer to not-implemented PME GPU 2D decomposition in X/Y
constexpr bool c_wrapX = true;
constexpr bool c_wrapY = true;

static int subGroupSizeFromVendor(const DeviceInformation& deviceInfo)
{
    switch (deviceInfo.deviceVendor)
    {
        case DeviceVendor::Amd: return 64;   // Handle RDNA2 devices, Issue #3972.
        case DeviceVendor::Intel: return 16; // TODO: Choose best value, Issue #4153.
        case DeviceVendor::Nvidia: return 32;
        default: GMX_RELEASE_ASSERT(false, "Unknown device vendor"); return 0;
    }
}

#define INSTANTIATE_3(order, computeSplines, spreadCharges, numGrids, writeGlobal, threadsPerAtom, subGroupSize) \
    extern template class PmeSplineAndSpreadKernel<order, computeSplines, spreadCharges, true, true, numGrids, writeGlobal, threadsPerAtom, subGroupSize>;

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
INSTANTIATE(4, 16);
#elif GMX_SYCL_HIPSYCL
INSTANTIATE(4, 32);
INSTANTIATE(4, 64);
#endif


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
}

PmeGpuProgramImpl::PmeGpuProgramImpl(const DeviceContext& deviceContext) :
    deviceContext_(deviceContext)
{
    // kernel parameters
    warpSize_             = subGroupSizeFromVendor(deviceContext.deviceInfo());
    spreadWorkGroupSize   = c_spreadMaxWarpsPerBlock * warpSize_;
    solveMaxWorkGroupSize = c_solveMaxWarpsPerBlock * warpSize_;
    gatherWorkGroupSize   = c_gatherMaxWarpsPerBlock * warpSize_;

    switch (warpSize_)
    {
#if GMX_SYCL_DPCPP
        case 16: setKernelPointers<16>(this); break;
#elif GMX_SYCL_HIPSYCL
        case 32: setKernelPointers<32>(this); break;
        case 64: setKernelPointers<64>(this); break;
#endif
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
}
