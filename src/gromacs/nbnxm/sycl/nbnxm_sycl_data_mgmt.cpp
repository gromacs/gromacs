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
 *  Stubs of functions that must be defined by nbnxm sycl implementation.
 *
 *  \ingroup module_nbnxm
 */
#include "gmxpre.h"

#include "gromacs/hardware/device_information.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/gpu_data_mgmt.h"
#include "gromacs/nbnxm/nbnxm_gpu.h"
#include "gromacs/nbnxm/nbnxm_gpu_data_mgmt.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/utility/exceptions.h"

#include "nbnxm_sycl_types.h"

namespace gmx
{

void gpu_init_platform_specific(NbnxmGpu* /* nb */)
{
    // Nothing specific in SYCL
}

void gpu_free_platform_specific(NbnxmGpu* /* nb */)
{
    // Nothing specific in SYCL
}

int gpu_min_ci_balanced(NbnxmGpu* nb)
{
    // SYCL-TODO: Logic and magic values taken from OpenCL
    static constexpr unsigned int balancedFactor = 50;
    if (nb == nullptr)
    {
        return 0;
    }
    const DeviceInformation& deviceInfo = nb->deviceContext_->deviceInfo();
    const sycl::device       device     = deviceInfo.syclDevice;
    const int numComputeUnits           = device.get_info<sycl::info::device::max_compute_units>();
    const int numComputeUnitsFactor     = getDeviceComputeUnitFactor(deviceInfo);
    return balancedFactor * numComputeUnits / numComputeUnitsFactor;
}

} // namespace gmx
