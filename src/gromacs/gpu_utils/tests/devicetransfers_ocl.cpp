/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * \brief Defines helper functionality for device transfers for tests
 * for GPU host allocator.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */
#include "gmxpre.h"

#include "gromacs/gpu_utils/gmxopencl.h"
#include "gromacs/gpu_utils/oclutils.h"
#include "gromacs/hardware/gpu_hw_info.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "devicetransfers.h"

namespace gmx
{
namespace
{

/*! \brief Help give useful diagnostics about error \c status while doing \c message.
 *
 * \throws InternalError  If status indicates failure, supplying
 *                        descriptive text from \c message. */
static void throwUponFailure(cl_int status, const char *message)
{
    if (status != CL_SUCCESS)
    {
        GMX_THROW(InternalError(formatString("Failure while %s, error was %s", message, ocl_get_error_string(status).c_str())));
    }
}

}   // namespace

void doDeviceTransfers(const gmx_gpu_info_t &gpuInfo,
                       ArrayRef<const char>  input,
                       ArrayRef<char>        output)
{
    GMX_RELEASE_ASSERT(input.size() == output.size(), "Input and output must have matching size");
    if (gpuInfo.n_dev == 0)
    {
        std::copy(input.begin(), input.end(), output.begin());
        return;
    }
    cl_int                status;

    const auto           &device       = gpuInfo.gpu_dev[0];
    cl_context_properties properties[] = {
        CL_CONTEXT_PLATFORM,
        (cl_context_properties) device.ocl_gpu_id.ocl_platform_id,
        0
    };
    // Give uncrustify more space

    auto deviceId = device.ocl_gpu_id.ocl_device_id;
    auto context  = clCreateContext(properties, 1, &deviceId, NULL, NULL, &status);
    throwUponFailure(status, "creating context");
    auto commandQueue = clCreateCommandQueue(context, deviceId, 0, &status);
    throwUponFailure(status, "creating command queue");

    auto devicePointer = clCreateBuffer(context, CL_MEM_READ_WRITE, input.size(), nullptr, &status);
    throwUponFailure(status, "creating buffer");

    status = clEnqueueWriteBuffer(commandQueue, devicePointer, CL_TRUE, 0, input.size(), input.data(), 0, nullptr, nullptr);
    throwUponFailure(status, "transferring host to device");
    status = clEnqueueReadBuffer(commandQueue, devicePointer, CL_TRUE, 0, output.size(), output.data(), 0, nullptr, nullptr);
    throwUponFailure(status, "transferring device to host");

    status = clReleaseMemObject(devicePointer);
    throwUponFailure(status, "releasing buffer");
    status = clReleaseContext(context);
    throwUponFailure(status, "releasing context");
}

} // namespace gmx
