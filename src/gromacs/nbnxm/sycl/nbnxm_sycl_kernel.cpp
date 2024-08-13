/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 *  \brief
 *  NBNXM SYCL kernels
 *
 *  \ingroup module_nbnxm
 */
#include "gmxpre.h"

#include "nbnxm_sycl_kernel.h"

#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/utility/template_mp.h"

#include "nbnxm_sycl_types.h"

namespace gmx
{

static int getNbnxmSubGroupSize(const DeviceInformation& deviceInfo)
{
    if (deviceInfo.supportedSubGroupSizes.size() == 1)
    {
        return deviceInfo.supportedSubGroupSizes[0];
    }
    else if (deviceInfo.supportedSubGroupSizes.size() > 1)
    {
        switch (deviceInfo.deviceVendor)
        {
            /* For Intel, choose 8 for 4x4 clusters, and 32 for 8x8 clusters.
             * The optimal one depends on the hardware, but we cannot choose c_nbnxnGpuClusterSize
             * at runtime anyway yet. */
            case DeviceVendor::Intel:
                return c_nbnxnGpuClusterSize * c_nbnxnGpuClusterSize / c_nbnxnGpuClusterpairSplit;
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

template<int subGroupSize, bool doPruneNBL, bool doCalcEnergies>
void launchNbnxmKernelHelper(NbnxmGpu* nb, const gmx::StepWorkload& stepWork, const InteractionLocality iloc);

// clang-format off
#if SYCL_NBNXM_SUPPORTS_SUBGROUP_SIZE_8
extern template void launchNbnxmKernelHelper<8, false, false>(NbnxmGpu* nb,const gmx::StepWorkload&  stepWork, const InteractionLocality iloc);
extern template void launchNbnxmKernelHelper<8, false, true>(NbnxmGpu* nb, const gmx::StepWorkload&  stepWork, const InteractionLocality iloc);
extern template void launchNbnxmKernelHelper<8, true, false>(NbnxmGpu* nb, const gmx::StepWorkload&  stepWork, const InteractionLocality iloc);
extern template void launchNbnxmKernelHelper<8, true, true>(NbnxmGpu* nb, const gmx::StepWorkload&  stepWork, const InteractionLocality iloc);
#endif
#if SYCL_NBNXM_SUPPORTS_SUBGROUP_SIZE_32
extern template void launchNbnxmKernelHelper<32, false, false>(NbnxmGpu* nb, const gmx::StepWorkload&  stepWork, const InteractionLocality iloc);
extern template void launchNbnxmKernelHelper<32, false, true>(NbnxmGpu* nb, const gmx::StepWorkload&  stepWork, const InteractionLocality iloc);
extern template void launchNbnxmKernelHelper<32, true, false>(NbnxmGpu* nb, const gmx::StepWorkload&  stepWork, const InteractionLocality iloc);
extern template void launchNbnxmKernelHelper<32, true, true>(NbnxmGpu* nb, const gmx::StepWorkload&  stepWork, const InteractionLocality iloc);
#endif
#if SYCL_NBNXM_SUPPORTS_SUBGROUP_SIZE_64
extern template void launchNbnxmKernelHelper<64, false, false>(NbnxmGpu* nb, const gmx::StepWorkload&  stepWork, const InteractionLocality iloc);
extern template void launchNbnxmKernelHelper<64, false, true>(NbnxmGpu* nb, const gmx::StepWorkload&  stepWork, const InteractionLocality iloc);
extern template void launchNbnxmKernelHelper<64, true, false>(NbnxmGpu* nb, const gmx::StepWorkload&  stepWork, const InteractionLocality iloc);
extern template void launchNbnxmKernelHelper<64, true, true>(NbnxmGpu* nb, const gmx::StepWorkload&  stepWork, const InteractionLocality iloc);
#endif
// clang-format on

template<int subGroupSize>
void launchNbnxmKernel(NbnxmGpu* nb, const gmx::StepWorkload& stepWork, const InteractionLocality iloc, bool doPrune)
{
    const bool doCalcEnergies = stepWork.computeEnergy;

    gmx::dispatchTemplatedFunction(
            [&](auto doPruneNBL_, auto doCalcEnergies_) {
                launchNbnxmKernelHelper<subGroupSize, doPruneNBL_, doCalcEnergies_>(nb, stepWork, iloc);
            },
            doPrune,
            doCalcEnergies);
}

void launchNbnxmKernel(NbnxmGpu* nb, const gmx::StepWorkload& stepWork, const InteractionLocality iloc, bool doPrune)
{
    const int subGroupSize = getNbnxmSubGroupSize(nb->deviceContext_->deviceInfo());
    switch (subGroupSize)
    {
        // Ensure any changes are in sync with device_management_sycl.cpp, nbnxm_sycl_kernel_body.h, and the #if above
#if SYCL_NBNXM_SUPPORTS_SUBGROUP_SIZE_8
        case 8: launchNbnxmKernel<8>(nb, stepWork, iloc, doPrune); break;
#endif
#if SYCL_NBNXM_SUPPORTS_SUBGROUP_SIZE_32
        case 32: launchNbnxmKernel<32>(nb, stepWork, iloc, doPrune); break;
#endif
#if SYCL_NBNXM_SUPPORTS_SUBGROUP_SIZE_64
        case 64: launchNbnxmKernel<64>(nb, stepWork, iloc, doPrune); break;
#endif
        default: GMX_RELEASE_ASSERT(false, "Unsupported sub-group size");
    }
}

} // namespace gmx
