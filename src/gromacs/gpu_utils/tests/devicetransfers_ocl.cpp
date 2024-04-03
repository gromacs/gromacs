/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
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
 * \brief Defines helper functionality for device transfers for tests
 * for GPU host allocator.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */
#include "gmxpre.h"

#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/gmxopencl.h"
#include "gromacs/gpu_utils/oclutils.h"
#include "gromacs/hardware/device_information.h"
#include "gromacs/hardware/device_management.h"
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
void throwUponFailure(cl_int status, const char* message)
{
    if (status != CL_SUCCESS)
    {
        GMX_THROW(InternalError(formatString(
                "Failure while %s, error was %s", message, ocl_get_error_string(status).c_str())));
    }
}

} // namespace

void doDeviceTransfers(const DeviceContext& deviceContext, ArrayRef<const char> input, ArrayRef<char> output)
{
    GMX_RELEASE_ASSERT(input.size() == output.size(), "Input and output must have matching size");

    cl_int status;

    auto deviceId     = deviceContext.deviceInfo().oclDeviceId;
    auto context      = deviceContext.context();
    auto commandQueue = clCreateCommandQueue(context, deviceId, 0, &status);
    throwUponFailure(status, "creating command queue");

    auto devicePointer = clCreateBuffer(context, CL_MEM_READ_WRITE, input.size(), nullptr, &status);
    throwUponFailure(status, "creating buffer");

    status = clEnqueueWriteBuffer(
            commandQueue, devicePointer, CL_TRUE, 0, input.size(), input.data(), 0, nullptr, nullptr);
    throwUponFailure(status, "transferring host to device");
    status = clEnqueueReadBuffer(
            commandQueue, devicePointer, CL_TRUE, 0, output.size(), output.data(), 0, nullptr, nullptr);
    throwUponFailure(status, "transferring device to host");

    status = clReleaseMemObject(devicePointer);
    throwUponFailure(status, "releasing buffer");
    status = clReleaseCommandQueue(commandQueue);
    throwUponFailure(status, "releasing command queue");
}

} // namespace gmx
