/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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

//! PME CUDA kernels forward declarations. Kernels are documented in their respective files.
template<const int order, const bool computeSplines, const bool spreadCharges, const bool wrapX, const bool wrapY, const bool writeGlobal, const bool orderThreads>
void pme_spline_and_spread_kernel(const PmeGpuCudaKernelParams kernelParams);

// Add extern declarations to inform that there will be a definition
// provided in another translation unit.
extern template void pme_spline_and_spread_kernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, true, true>(
        const PmeGpuCudaKernelParams);
extern template void pme_spline_and_spread_kernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, true, false>(
        const PmeGpuCudaKernelParams);
extern template void pme_spline_and_spread_kernel<c_pmeOrder, true, false, c_wrapX, c_wrapY, true, true>(
        const PmeGpuCudaKernelParams);
extern template void pme_spline_and_spread_kernel<c_pmeOrder, true, false, c_wrapX, c_wrapY, true, false>(
        const PmeGpuCudaKernelParams);
extern template void pme_spline_and_spread_kernel<c_pmeOrder, false, true, c_wrapX, c_wrapY, true, true>(
        const PmeGpuCudaKernelParams);
extern template void pme_spline_and_spread_kernel<c_pmeOrder, false, true, c_wrapX, c_wrapY, true, false>(
        const PmeGpuCudaKernelParams);
extern template void pme_spline_and_spread_kernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, false, true>(
        const PmeGpuCudaKernelParams);
extern template void pme_spline_and_spread_kernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, false, false>(
        const PmeGpuCudaKernelParams);

template<GridOrdering gridOrdering, bool computeEnergyAndVirial>
void pme_solve_kernel(const PmeGpuCudaKernelParams kernelParams);

// Add extern declarations to inform that there will be a definition
// provided in another translation unit.
extern template void pme_solve_kernel<GridOrdering::XYZ, false>(const PmeGpuCudaKernelParams);
extern template void pme_solve_kernel<GridOrdering::XYZ, true>(const PmeGpuCudaKernelParams);
extern template void pme_solve_kernel<GridOrdering::YZX, false>(const PmeGpuCudaKernelParams);
extern template void pme_solve_kernel<GridOrdering::YZX, true>(const PmeGpuCudaKernelParams);

template<const int order, const bool overwriteForces, const bool wrapX, const bool wrapY, const bool readGlobal, const bool orderThreads>
void pme_gather_kernel(const PmeGpuCudaKernelParams kernelParams);

// Add extern declarations to inform that there will be a definition
// provided in another translation unit.
extern template void pme_gather_kernel<c_pmeOrder, true, c_wrapX, c_wrapY, true, true>(const PmeGpuCudaKernelParams);
extern template void pme_gather_kernel<c_pmeOrder, true, c_wrapX, c_wrapY, false, true>(const PmeGpuCudaKernelParams);
extern template void pme_gather_kernel<c_pmeOrder, false, c_wrapX, c_wrapY, true, true>(const PmeGpuCudaKernelParams);
extern template void pme_gather_kernel<c_pmeOrder, false, c_wrapX, c_wrapY, false, true>(const PmeGpuCudaKernelParams);
extern template void pme_gather_kernel<c_pmeOrder, true, c_wrapX, c_wrapY, true, false>(const PmeGpuCudaKernelParams);
extern template void pme_gather_kernel<c_pmeOrder, true, c_wrapX, c_wrapY, false, false>(const PmeGpuCudaKernelParams);
extern template void pme_gather_kernel<c_pmeOrder, false, c_wrapX, c_wrapY, true, false>(const PmeGpuCudaKernelParams);
extern template void pme_gather_kernel<c_pmeOrder, false, c_wrapX, c_wrapY, false, false>(const PmeGpuCudaKernelParams);

PmeGpuProgramImpl::PmeGpuProgramImpl(const gmx_device_info_t*)
{
    // kernel parameters
    warpSize              = warp_size;
    spreadWorkGroupSize   = c_spreadMaxThreadsPerBlock;
    solveMaxWorkGroupSize = c_solveMaxThreadsPerBlock;
    gatherWorkGroupSize   = c_gatherMaxThreadsPerBlock;

    /*!
     * Not all combinations of the splineAndSpread, spline and Spread kernels are required
     * If only the spline (without the spread) then it does not make sense not to write the data to global memory
     * Similarly the spread kernel (without the spline) implies that we should read the spline data from global memory
     */
    splineAndSpreadKernel =
            pme_spline_and_spread_kernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, false, false>;
    splineAndSpreadKernelThPerAtom4 =
            pme_spline_and_spread_kernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, false, true>;
    splineAndSpreadKernelWriteSplines =
            pme_spline_and_spread_kernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, true, false>;
    splineAndSpreadKernelWriteSplinesThPerAtom4 =
            pme_spline_and_spread_kernel<c_pmeOrder, true, true, c_wrapX, c_wrapY, true, true>;
    splineKernel = pme_spline_and_spread_kernel<c_pmeOrder, true, false, c_wrapX, c_wrapY, true, false>;
    splineKernelThPerAtom4 =
            pme_spline_and_spread_kernel<c_pmeOrder, true, false, c_wrapX, c_wrapY, true, true>;
    spreadKernel = pme_spline_and_spread_kernel<c_pmeOrder, false, true, c_wrapX, c_wrapY, true, false>;
    spreadKernelThPerAtom4 =
            pme_spline_and_spread_kernel<c_pmeOrder, false, true, c_wrapX, c_wrapY, true, true>;
    gatherKernel            = pme_gather_kernel<c_pmeOrder, true, c_wrapX, c_wrapY, false, false>;
    gatherKernelThPerAtom4  = pme_gather_kernel<c_pmeOrder, true, c_wrapX, c_wrapY, false, true>;
    gatherKernelReadSplines = pme_gather_kernel<c_pmeOrder, true, c_wrapX, c_wrapY, true, false>;
    gatherKernelReadSplinesThPerAtom4 = pme_gather_kernel<c_pmeOrder, true, c_wrapX, c_wrapY, true, true>;
    gatherReduceWithInputKernel = pme_gather_kernel<c_pmeOrder, false, c_wrapX, c_wrapY, false, false>;
    gatherReduceWithInputKernelThPerAtom4 =
            pme_gather_kernel<c_pmeOrder, false, c_wrapX, c_wrapY, false, true>;
    gatherReduceWithInputKernelReadSplines =
            pme_gather_kernel<c_pmeOrder, false, c_wrapX, c_wrapY, true, false>;
    gatherReduceWithInputKernelReadSplinesThPerAtom4 =
            pme_gather_kernel<c_pmeOrder, false, c_wrapX, c_wrapY, true, true>;
    solveXYZKernel       = pme_solve_kernel<GridOrdering::XYZ, false>;
    solveXYZEnergyKernel = pme_solve_kernel<GridOrdering::XYZ, true>;
    solveYZXKernel       = pme_solve_kernel<GridOrdering::YZX, false>;
    solveYZXEnergyKernel = pme_solve_kernel<GridOrdering::YZX, true>;
}

PmeGpuProgramImpl::~PmeGpuProgramImpl() {}
