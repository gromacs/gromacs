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
 *
 * \brief Implements the DeviceStream for OpenCL.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_gpu_utils
 */
#include "gmxpre.h"

#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/gputraits_ocl.h"
#include "gromacs/gpu_utils/oclutils.h"
#include "gromacs/hardware/device_information.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

DeviceStream::DeviceStream(const DeviceContext& deviceContext,
                           DeviceStreamPriority /* priority */,
                           const bool useTiming)
{
    const DeviceInformation&    deviceInfo      = deviceContext.deviceInfo();
    cl_command_queue_properties queueProperties = useTiming ? CL_QUEUE_PROFILING_ENABLE : 0;
    cl_device_id                deviceId        = deviceInfo.oclDeviceId;
    cl_int                      clError;
    stream_ = clCreateCommandQueue(deviceContext.context(), deviceId, queueProperties, &clError);
    if (clError != CL_SUCCESS)
    {
        GMX_THROW(gmx::InternalError(gmx::formatString(
                "Failed to create OpenCL command queue on GPU %s (OpenCL error ID %d).",
                deviceInfo.device_name,
                clError)));
    }
}

DeviceStream::~DeviceStream()
{
    if (isValid())
    {
        cl_int clError = clReleaseCommandQueue(stream_);
        GMX_RELEASE_ASSERT(
                clError == CL_SUCCESS,
                gmx::formatString("Failed to release OpenCL stream (OpenCL error ID %d).", clError).c_str());
        stream_ = nullptr;
    }
}

cl_command_queue DeviceStream::stream() const
{
    return stream_;
}

bool DeviceStream::isValid() const
{
    return (stream_ != nullptr);
}

void DeviceStream::synchronize() const
{
    cl_int clError = clFinish(stream_);
    GMX_RELEASE_ASSERT(
            CL_SUCCESS == clError,
            gmx::formatString("Error caught during clFinish (OpenCL error ID %d).", clError).c_str());
}

void issueClFlushInStream(const DeviceStream& deviceStream)
{
    cl_int cl_error = clFlush(deviceStream.stream());
    if (cl_error != CL_SUCCESS)
    {
        GMX_THROW(gmx::InternalError("clFlush failed: " + ocl_get_error_string(cl_error)));
    }
}
