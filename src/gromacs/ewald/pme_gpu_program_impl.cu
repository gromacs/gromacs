/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * \ingroup module_ewald
 */
#include "gmxpre.h"

#include "pme_gpu_program_impl.h"

#include "pme_gpu_constants.h"
#include "pme_gpu_internal.h" // for GridOrdering enum
#include "pme_gpu_types_host.h"

// PME interpolation order
constexpr int c_pmeOrder = 4;
// These hardcoded spread/gather parameters refer to not-implemented PME GPU 2D decomposition in X/Y
constexpr bool c_wrapX = true;
constexpr bool c_wrapY = true;

constexpr int c_stateA = 0;
constexpr int c_stateB = 1;

//! PME CUDA kernels forward declarations. Kernels are documented in their respective files.
template<int order, bool computeSplines, bool spreadCharges, bool wrapX, bool wrapY, int mode, bool writeGlobal, ThreadsPerAtom threadsPerAtom>
__global__ void pme_spline_and_spread_kernel(PmeGpuCudaKernelParams kernelParams);

// Add extern declarations to inform that there will be a definition
// provided in another translation unit.
// clang-format off
extern template void __global__
pme_spline_and_spread_kernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::Order>(const PmeGpuCudaKernelParams);
extern template void __global__
pme_spline_and_spread_kernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::OrderSquared>(const PmeGpuCudaKernelParams);
extern template void __global__
pme_spline_and_spread_kernel<c_pmeOrder, true, false, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::Order>(const PmeGpuCudaKernelParams);
extern template void __global__
pme_spline_and_spread_kernel<c_pmeOrder, true, false, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::OrderSquared>(const PmeGpuCudaKernelParams);
extern template __global__ void
pme_spline_and_spread_kernel<c_pmeOrder, false, true, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::Order>(const PmeGpuCudaKernelParams);
extern template __global__ void
pme_spline_and_spread_kernel<c_pmeOrder, false, true, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::OrderSquared>(const PmeGpuCudaKernelParams);
extern template __global__ void
pme_spline_and_spread_kernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 1, false, ThreadsPerAtom::Order>(const PmeGpuCudaKernelParams);
extern template __global__ void
pme_spline_and_spread_kernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 1, false, ThreadsPerAtom::OrderSquared>(const PmeGpuCudaKernelParams);
extern template __global__ void
pme_spline_and_spread_kernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::Order>(const PmeGpuCudaKernelParams);
extern template __global__ void
pme_spline_and_spread_kernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::OrderSquared>(const PmeGpuCudaKernelParams);
extern template __global__ void
pme_spline_and_spread_kernel<c_pmeOrder, true, false, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::Order>(const PmeGpuCudaKernelParams);
extern template __global__ void
pme_spline_and_spread_kernel<c_pmeOrder, true, false, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::OrderSquared>(const PmeGpuCudaKernelParams);
extern template __global__ void
pme_spline_and_spread_kernel<c_pmeOrder, false, true, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::Order>(const PmeGpuCudaKernelParams);
extern template __global__ void
pme_spline_and_spread_kernel<c_pmeOrder, false, true, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::OrderSquared>(const PmeGpuCudaKernelParams);
extern template __global__ void
pme_spline_and_spread_kernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 2, false, ThreadsPerAtom::Order>(const PmeGpuCudaKernelParams);
extern template __global__ void
pme_spline_and_spread_kernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 2, false, ThreadsPerAtom::OrderSquared>(const PmeGpuCudaKernelParams);

template<GridOrdering gridOrdering, bool computeEnergyAndVirial, const int gridIndex> /* It is significantly slower to pass gridIndex as a kernel parameter */
__global__ void pme_solve_kernel(PmeGpuCudaKernelParams kernelParams);

// Add extern declarations to inform that there will be a definition
// provided in another translation unit.
// clang-format off
extern template __global__ void pme_solve_kernel<GridOrdering::XYZ, false, c_stateA>(const PmeGpuCudaKernelParams);
extern template __global__ void pme_solve_kernel<GridOrdering::XYZ, true, c_stateA>(const PmeGpuCudaKernelParams);
extern template __global__ void pme_solve_kernel<GridOrdering::YZX, false, c_stateA>(const PmeGpuCudaKernelParams);
extern template __global__ void pme_solve_kernel<GridOrdering::YZX, true, c_stateA>(const PmeGpuCudaKernelParams);
extern template __global__ void pme_solve_kernel<GridOrdering::XYZ, false, c_stateB>(const PmeGpuCudaKernelParams);
extern template __global__ void pme_solve_kernel<GridOrdering::XYZ, true, c_stateB>(const PmeGpuCudaKernelParams);
extern template __global__ void pme_solve_kernel<GridOrdering::YZX, false, c_stateB>(const PmeGpuCudaKernelParams);
extern template __global__ void pme_solve_kernel<GridOrdering::YZX, true, c_stateB>(const PmeGpuCudaKernelParams);
// clang-format on

template<int order, bool wrapX, bool wrapY, int nGrids, bool readGlobal, ThreadsPerAtom threadsPerAtom>
__global__ void pme_gather_kernel(PmeGpuCudaKernelParams kernelParams);

__global__ void nvshmemSignalKernel(PmeGpuCudaKernelParams kernelParams);

// Add extern declarations to inform that there will be a definition
// provided in another translation unit.
// clang-format off
extern template __global__ void pme_gather_kernel<c_pmeOrder, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::Order>        (const PmeGpuCudaKernelParams);
extern template __global__ void pme_gather_kernel<c_pmeOrder, c_wrapX, c_wrapY, 1, false, ThreadsPerAtom::Order>       (const PmeGpuCudaKernelParams);
extern template __global__ void pme_gather_kernel<c_pmeOrder, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::OrderSquared> (const PmeGpuCudaKernelParams);
extern template __global__ void pme_gather_kernel<c_pmeOrder, c_wrapX, c_wrapY, 1, false, ThreadsPerAtom::OrderSquared>(const PmeGpuCudaKernelParams);
extern template __global__ void pme_gather_kernel<c_pmeOrder, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::Order>          (const PmeGpuCudaKernelParams);
extern template __global__ void pme_gather_kernel<c_pmeOrder, c_wrapX, c_wrapY, 2, false, ThreadsPerAtom::Order>         (const PmeGpuCudaKernelParams);
extern template __global__ void pme_gather_kernel<c_pmeOrder, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::OrderSquared>   (const PmeGpuCudaKernelParams);
extern template __global__ void pme_gather_kernel<c_pmeOrder, c_wrapX, c_wrapY, 2, false, ThreadsPerAtom::OrderSquared>  (const PmeGpuCudaKernelParams);
// clang-format on

PmeGpuProgramImpl::PmeGpuProgramImpl(const DeviceContext& deviceContext) :
    deviceContext_(deviceContext)
{
    // kernel parameters
    warpSize_             = warp_size;
    spreadWorkGroupSize   = c_spreadMaxThreadsPerBlock;
    solveMaxWorkGroupSize = c_solveMaxThreadsPerBlock;
    gatherWorkGroupSize   = c_gatherMaxThreadsPerBlock;

    /* Not all combinations of the splineAndSpread, spline and Spread kernels are required
     * If only the spline (without the spread) then it does not make sense not to write the data to global memory
     * Similarly the spread kernel (without the spline) implies that we should read the spline data from global memory
     */
    // clang-format off
    splineAndSpreadKernelSingle                       = pme_spline_and_spread_kernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 1, false, ThreadsPerAtom::OrderSquared>;
    splineAndSpreadKernelThPerAtom4Single             = pme_spline_and_spread_kernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 1, false, ThreadsPerAtom::Order>;
    splineAndSpreadKernelWriteSplinesSingle           = pme_spline_and_spread_kernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::OrderSquared>;
    splineAndSpreadKernelWriteSplinesThPerAtom4Single = pme_spline_and_spread_kernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::Order>;
    splineKernelSingle                                = pme_spline_and_spread_kernel<c_pmeOrder, true, false, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::OrderSquared>;
    splineKernelThPerAtom4Single                      = pme_spline_and_spread_kernel<c_pmeOrder, true, false, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::Order>;
    spreadKernelSingle                                = pme_spline_and_spread_kernel<c_pmeOrder, false, true, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::OrderSquared>;
    spreadKernelThPerAtom4Single                      = pme_spline_and_spread_kernel<c_pmeOrder, false, true, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::Order>;
    splineAndSpreadKernelDual                         = pme_spline_and_spread_kernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 2, false, ThreadsPerAtom::OrderSquared>;
    splineAndSpreadKernelThPerAtom4Dual               = pme_spline_and_spread_kernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 2, false, ThreadsPerAtom::Order>;
    splineAndSpreadKernelWriteSplinesDual             = pme_spline_and_spread_kernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::OrderSquared>;
    splineAndSpreadKernelWriteSplinesThPerAtom4Dual   = pme_spline_and_spread_kernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::Order>;
    splineKernelDual                                  = pme_spline_and_spread_kernel<c_pmeOrder, true, false, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::OrderSquared>;
    splineKernelThPerAtom4Dual                        = pme_spline_and_spread_kernel<c_pmeOrder, true, false, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::Order>;
    spreadKernelDual                                  = pme_spline_and_spread_kernel<c_pmeOrder, false, true, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::OrderSquared>;
    spreadKernelThPerAtom4Dual                        = pme_spline_and_spread_kernel<c_pmeOrder, false, true, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::Order>;
    gatherKernelSingle                                = pme_gather_kernel<c_pmeOrder, c_wrapX, c_wrapY, 1, false, ThreadsPerAtom::OrderSquared>;
    gatherKernelThPerAtom4Single                      = pme_gather_kernel<c_pmeOrder, c_wrapX, c_wrapY, 1, false, ThreadsPerAtom::Order>;
    gatherKernelReadSplinesSingle                     = pme_gather_kernel<c_pmeOrder, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::OrderSquared>;
    gatherKernelReadSplinesThPerAtom4Single           = pme_gather_kernel<c_pmeOrder, c_wrapX, c_wrapY, 1, true, ThreadsPerAtom::Order>;
    gatherKernelDual                                  = pme_gather_kernel<c_pmeOrder, c_wrapX, c_wrapY, 2, false, ThreadsPerAtom::OrderSquared>;
    gatherKernelThPerAtom4Dual                        = pme_gather_kernel<c_pmeOrder, c_wrapX, c_wrapY, 2, false, ThreadsPerAtom::Order>;
    gatherKernelReadSplinesDual                       = pme_gather_kernel<c_pmeOrder, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::OrderSquared>;
    gatherKernelReadSplinesThPerAtom4Dual             = pme_gather_kernel<c_pmeOrder, c_wrapX, c_wrapY, 2, true, ThreadsPerAtom::Order>;
    solveXYZKernelA                                   = pme_solve_kernel<GridOrdering::XYZ, false, c_stateA>;
    solveXYZEnergyKernelA                             = pme_solve_kernel<GridOrdering::XYZ, true, c_stateA>;
    solveYZXKernelA                                   = pme_solve_kernel<GridOrdering::YZX, false, c_stateA>;
    solveYZXEnergyKernelA                             = pme_solve_kernel<GridOrdering::YZX, true, c_stateA>;
    solveXYZKernelB                                   = pme_solve_kernel<GridOrdering::XYZ, false, c_stateB>;
    solveXYZEnergyKernelB                             = pme_solve_kernel<GridOrdering::XYZ, true, c_stateB>;
    solveYZXKernelB                                   = pme_solve_kernel<GridOrdering::YZX, false, c_stateB>;
    solveYZXEnergyKernelB                             = pme_solve_kernel<GridOrdering::YZX, true, c_stateB>;

    nvshmemSignalKern = nvshmemSignalKernel;

    // clang-format on
}

PmeGpuProgramImpl::~PmeGpuProgramImpl() = default;
