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
#ifndef GMX_GPU_UTILS_DEVICE_CONTEXT_H
#define GMX_GPU_UTILS_DEVICE_CONTEXT_H

/*! \libinternal \file
 *
 * \brief Declarations for DeviceContext class.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_gpu_utils
 * \inlibraryapi
 */

#include "config.h"

#if GMX_GPU_OPENCL
#    include "gromacs/gpu_utils/gmxopencl.h"
#endif
#if GMX_GPU_SYCL
#    include "gromacs/gpu_utils/gmxsycl.h"
#endif

#include "gromacs/gpu_utils/pmalloc.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/utility/classhelpers.h"

struct DeviceInformation;

// Stub for device context
class DeviceContext
{
public:
    //! Constructs context and activates the device.
    DeviceContext(const DeviceInformation& deviceInfo);
    //! Destructor
    // NOLINTNEXTLINE(performance-trivially-destructible)
    ~DeviceContext();

    //! Get the associated device information
    const DeviceInformation& deviceInfo() const { return deviceInfo_; }

    //! Activate the device
    void activate() const
    {
        setActiveDevice(deviceInfo_);
        pmallocSetDefaultDeviceContext(this);
    }

private:
    //! A reference to the device information used upon context creation
    const DeviceInformation& deviceInfo_;

#if GMX_GPU_OPENCL
public:
    //! Getter
    cl_context context() const;

private:
    //! OpenCL context object
    cl_context context_ = nullptr;
#endif

#if GMX_GPU_SYCL
public:
    //! Const getter
    const sycl::context& context() const { return context_; }
    //! Getter
    sycl::context& context() { return context_; }

private:
    //! SYCL context object
    sycl::context context_;
#endif

    GMX_DISALLOW_COPY_MOVE_AND_ASSIGN(DeviceContext);
};

#endif // GMX_GPU_UTILS_DEVICE_CONTEXT_H
