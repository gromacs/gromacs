/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2017 The GROMACS development team.
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
 *  \brief Defines the implementations of the device management that are common for CPU, CUDA and OpenCL.
 *
 *  \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_hardware
 */
#include "gmxpre.h"

#include <assert.h>

#include "gromacs/hardware/device_information.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/hardware/gpu_hw_info.h"
#include "gromacs/utility/smalloc.h"

bool canPerformGpuDetection()
{
    if (c_binarySupportsGpus && getenv("GMX_DISABLE_GPU_DETECTION") == nullptr)
    {
        return isGpuDetectionFunctional(nullptr);
    }
    else
    {
        return false;
    }
}

std::vector<int> getCompatibleGpus(const gmx_gpu_info_t& gpu_info)
{
    // Possible minor over-allocation here, but not important for anything
    std::vector<int> compatibleGpus;
    compatibleGpus.reserve(gpu_info.n_dev);
    for (int i = 0; i < gpu_info.n_dev; i++)
    {
        assert(gpu_info.deviceInfo);
        if (gpu_info_get_stat(gpu_info, i) == DeviceStatus::Compatible)
        {
            compatibleGpus.push_back(i);
        }
    }
    return compatibleGpus;
}

const char* getGpuCompatibilityDescription(const gmx_gpu_info_t& gpu_info, int index)
{
    return (index >= gpu_info.n_dev ? c_deviceStateString[DeviceStatus::Nonexistent]
                                    : c_deviceStateString[gpu_info_get_stat(gpu_info, index)]);
}

void free_gpu_info(const gmx_gpu_info_t* gpu_info)
{
    sfree(static_cast<void*>(gpu_info->deviceInfo)); // circumvent is_pod check in sfree
}
