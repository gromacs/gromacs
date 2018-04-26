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
#include "pme-persistent-data.h"
#include "pme-gpu-types-host.h"
#include "pme-internal.h"
#include "pme-gpu-internal.h" // for GridOrdering enum

#if GMX_GPU == GMX_GPU_CUDA

//! PME CUDA kernels forward declarations. Kernels are documented in their respective files.
template <
    const int order,
    const bool computeSplines,
    const bool spreadCharges,
    const bool wrapX,
    const bool wrapY
    >
void pme_spline_and_spread_kernel(const PmeGpuCudaKernelParams kernelParams);
template<
    GridOrdering gridOrdering,
    bool computeEnergyAndVirial
    >
void pme_solve_kernel(const PmeGpuCudaKernelParams kernelParams);
template <
    const int order,
    const bool overwriteForces,
    const bool wrapX,
    const bool wrapY
    >
void pme_gather_kernel(const PmeGpuCudaKernelParams kernelParams);
#endif

PmeGpuPersistentData::PmeGpuPersistentData(gmx_pme_t *pme)
{
    //TODO in OpenCL this will create cl_context, compile kernels/program, initialize clFFT library

#if GMX_GPU == GMX_GPU_CUDA
    GMX_ASSERT(pme->pme_order == 4, "The only supported PME GPU order is 4");
    constexpr int  order = 4;
    // These hardcoded spread/gather parameters refer to not-implemented PME GPU 2D decomposition
    constexpr bool wrapX = true;
    constexpr bool wrapY = true;
    //GMX_UNUSED_VALUE(wrapX);
    //GMX_UNUSED_VALUE(wrapY);
    splineAndSpreadKernel       = pme_spline_and_spread_kernel<order, true, true, wrapX, wrapY>;
    splineKernel                = pme_spline_and_spread_kernel<order, true, false, wrapX, wrapY>;
    spreadKernel                = pme_spline_and_spread_kernel<order, false, true, wrapX, wrapY>;
    gatherKernel                = pme_gather_kernel<order, true, wrapX, wrapY>;
    gatherReduceWithInputKernel = pme_gather_kernel<order, false, wrapX, wrapY>;
    solveXYZKernel              = pme_solve_kernel<GridOrdering::XYZ, false>;
    solveXYZEnergyKernel        = pme_solve_kernel<GridOrdering::XYZ, true>;
    solveYZXKernel              = pme_solve_kernel<GridOrdering::YZX, false>;
    solveYZXEnergyKernel        = pme_solve_kernel<GridOrdering::YZX, true>;
#endif
}

PmeGpuPersistentData::~PmeGpuPersistentData()
{
    //TODO in OpenCL this will teardown clFFT, release cl_context, kernels, program (unless we use RAII for those)
}

PmePersistentDataHandle buildPersistentPmeData(gmx_device_info_t *deviceInfo)
{
#if GMX_GPU == GMX_GPU_NONE
    return nullptr;
#else
    // Currently a dummy structure will do to bypass proper initialization,
    // avoiding another clunky gmx_pme_init() call might get refactored anyway (?)
    gmx_pme_t dummy = {0};
    dummy.pme_order = 4;
    PmeGpu    dummyGpu = {0};
    dummy.gpu             = &dummyGpu;
    dummy.gpu->deviceInfo = deviceInfo;
    PmePersistentDataHandle result = std::make_shared<PmeGpuPersistentData>(&dummy);
    return result;
#endif
}
