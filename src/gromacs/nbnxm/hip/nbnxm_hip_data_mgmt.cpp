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
/*! \file
 *  \brief Define HIP implementation for GPU data transfer for NBNXM module
 *
 *  \author Szilard Pall <pall.szilard@gmail.com>
 *  \author Paul Bauer <paul.bauer.q@gmail.com>
 */

#include "gmxpre.h"

// The compiler generates the wrong code when calling rocprim for gfx1034 devices, so we need to make sure that it doesn't try to
// use unsupported dpp instructions. Tracked here, but not fixed even if the ticket says so: https://github.com/ROCm/rocPRIM/issues/452
#if __gfx1034__
#    define ROCPRIM_DISABLE_DPP
#    define ROCPRIM_DETAIL_USE_DPP false
#endif
#include <rocprim/rocprim.hpp>

// TODO We would like to move this down, but the way NbnxmGpu
//      is currently declared means this has to be before gpu_types.h
#include "nbnxm_hip_types.h"

// TODO Remove this comment when the above order issue is resolved
#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/pmalloc.h"
#include "gromacs/hardware/device_information.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/gpu_data_mgmt.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/nbnxm_gpu_data_mgmt.h"

namespace gmx
{

/* This is a heuristically determined parameter for the Kepler
 * and Maxwell architectures for the minimum size of ci lists by multiplying
 * this constant with the # of multiprocessors on the current device.
 * Since the maximum number of blocks per multiprocessor is 16, the ideal
 * count for small systems is 32 or 48 blocks per multiprocessor. Because
 * there is a bit of fluctuations in the generated block counts, we use
 * a target of 44 instead of the ideal value of 48.
 */
static const unsigned int gpu_min_ci_balanced_factor = 44;

void gpu_init_platform_specific(NbnxmGpu* /* nb */)
{
    // Nothing specific in HIP
}

void gpu_free_platform_specific(NbnxmGpu* /* nb */)
{
    // Nothing specific in HIP
}

int gpu_min_ci_balanced(NbnxmGpu* nb)
{
    return nb != nullptr ? gpu_min_ci_balanced_factor * nb->deviceContext_->deviceInfo().prop.multiProcessorCount
                         : 0;
}

namespace
{

size_t hipRocprimWrapper(size_t              temporaryBufferSize,
                         char*               temporaryBuffer,
                         gmx::GpuPairlist*   plist,
                         const DeviceStream& deviceStream)
{
    size_t size = temporaryBufferSize;

    hipError_t stat = rocprim::exclusive_scan(temporaryBuffer,
                                              size,
                                              plist->sorting.sciHistogram,
                                              plist->sorting.sciOffset,
                                              0,
                                              c_sciHistogramSize,
                                              rocprim::plus<int>(),
                                              deviceStream.stream());
    gmx::checkDeviceError(stat, "rocprim::exclusive_scan failed");

    return size;
}

} // namespace

size_t getExclusiveScanWorkingArraySize(GpuPairlist* plist, const DeviceStream& deviceStream)
{
    return hipRocprimWrapper(0, nullptr, plist, deviceStream);
}

void performExclusiveScan(size_t              temporaryBufferSize,
                          char*               temporaryBuffer,
                          GpuPairlist*        plist,
                          const DeviceStream& deviceStream)
{
    std::ignore = hipRocprimWrapper(temporaryBufferSize, temporaryBuffer, plist, deviceStream);
}

// reset our custom defines
#if __gfx1034__
#    undef ROCPRIM_DISABLE_DPP
#    undef ROCPRIM_DETAIL_USE_DPP
#endif

} // namespace gmx
