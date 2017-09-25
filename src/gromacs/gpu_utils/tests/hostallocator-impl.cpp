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
 * \brief
 * Defines helper functionality for tests for GPU host allocator.
 *
 * Undefined symbols in Google Test, GROMACS use of -Wundef, and the
 * implementation of FindCUDA.cmake and/or nvcc mean that no
 * compilation unit should include a gtest header while being compiled
 * by nvcc. None of -isystem, -Wno-undef, nor pragma GCC diagnostic
 * work.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */
#include "gmxpre.h"

#include "hostallocator-impl.h"

#include "config.h"

#include "gromacs/hardware/gpu_hw_info.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#if GMX_GPU == GMX_GPU_CUDA

#include "gromacs/gpu_utils/cudautils.cuh"

#endif

#if GMX_GPU == GMX_GPU_OPENCL

#include "gromacs/gpu_utils/gmxopencl.h"
#include "gromacs/gpu_utils/oclutils.h"

#endif

namespace gmx
{
namespace
{

//! Convenience type
using StatusType =
#if GMX_GPU == GMX_GPU_CUDA
        cudaError_t
#elif GMX_GPU == GMX_GPU_OPENCL
        cl_int
#else
        void *
#endif
;

/*! \brief Help give useful diagnostics when throwing.
 *
 * \throws InternalError  If status indicates failure, supplying
 *                        descriptive text from \c message. */
static void throwUponFailure(StatusType status, const char *message)
{
#if GMX_GPU == GMX_GPU_CUDA
    if (status != cudaSuccess)
    {
        GMX_THROW(InternalError(formatString("Failure while %s", message)));;
    }
#elif GMX_GPU == GMX_GPU_OPENCL
    if (status != CL_SUCCESS)
    {
        GMX_THROW(InternalError(formatString("Failure while %s, error was %s", message, ocl_get_error_string(status).c_str())));
    }
#else
    GMX_UNUSED_VALUE(status);
    GMX_UNUSED_VALUE(message);
#endif
}

}   // namespace

void doDeviceTransfers(const gmx_gpu_info_t &gpuInfo,
                       ConstArrayRef<char>   input,
                       ArrayRef<char>       *output)
{
    GMX_RELEASE_ASSERT(input.size() == output->size(), "Input and output must have matching size");
    StatusType status;
    GMX_RELEASE_ASSERT(gpuInfo.n_dev > 0, "Must have a GPU device");

#if GMX_GPU == GMX_GPU_CUDA
    const auto &device = gpuInfo.gpu_dev[0];
    int         oldDeviceId;

    status = cudaGetDevice(&oldDeviceId);
    throwUponFailure(status, "getting old device id");
    status = cudaSetDevice(device.id);
    throwUponFailure(status, "setting device id to 0");

    void       *devicePointer;
    status = cudaMalloc(&devicePointer, input.size());
    throwUponFailure(status, "creating buffer");

    status = cudaMemcpy(devicePointer, input.data(), input.size(), cudaMemcpyHostToDevice);
    throwUponFailure(status, "transferring host to device");
    status = cudaMemcpy(output->data(), devicePointer, output->size(), cudaMemcpyDeviceToHost);
    throwUponFailure(status, "transferring device to host");

    status = cudaFree(devicePointer);
    throwUponFailure(status, "releasing buffer");

    status = cudaSetDevice(oldDeviceId);
    throwUponFailure(status, "setting old device id");

#elif GMX_GPU == GMX_GPU_OPENCL
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
    status = clEnqueueReadBuffer(commandQueue, devicePointer, CL_TRUE, 0, output->size(), output->data(), 0, nullptr, nullptr);
    throwUponFailure(status, "transferring device to host");

    status = clReleaseMemObject(devicePointer);
    throwUponFailure(status, "releasing buffer");
    status = clReleaseContext(context);
    throwUponFailure(status, "releasing context");
#else
    GMX_UNUSED_VALUE(gpuInfo);
    GMX_UNUSED_VALUE(input);
    GMX_UNUSED_VALUE(output);
    status = nullptr;
    // Avoid warning about unused throwUponFailure function
    throwUponFailure(status, "dummy");
#endif
}

} // namespace gmx
