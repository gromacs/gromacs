/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
#ifndef GMX_GPU_UTILS_DEVICEBUFFER_CUH
#define GMX_GPU_UTILS_DEVICEBUFFER_CUH

/*! \libinternal \file
 *  \brief Implements the DeviceBuffer type and routines for CUDA.
 *  Should only be included directly by the main DeviceBuffer file devicebuffer.h.
 *  TODO: the intent is for DeviceBuffer to become a class.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 *
 *  \inlibraryapi
 */

#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/gpu_utils/gpu_utils.h" //only for GpuApiCallBehavior
#include "gromacs/gpu_utils/gputraits.cuh"
#include "gromacs/utility/gmxassert.h"

/*! \brief
 * Allocates a device-side buffer.
 * It is currently a caller's responsibility to call it only on not-yet allocated buffers.
 *
 * \tparam        ValueType            Raw value type of the \p buffer.
 * \param[in,out] buffer               Pointer to the device-side buffer.
 * \param[in]     numValues            Number of values to accomodate.
 * \param[in]     deviceContext        The buffer's dummy device  context - not managed explicitly in CUDA RT.
 */
template<typename ValueType>
void allocateDeviceBuffer(DeviceBuffer<ValueType>* buffer, size_t numValues, DeviceContext /* deviceContext */)
{
    GMX_ASSERT(buffer, "needs a buffer pointer");
    cudaError_t stat = cudaMalloc((void**)buffer, numValues * sizeof(ValueType));
    GMX_RELEASE_ASSERT(stat == cudaSuccess, "cudaMalloc failure");
}

/*! \brief
 * Frees a device-side buffer.
 * This does not reset separately stored size/capacity integers,
 * as this is planned to be a destructor of DeviceBuffer as a proper class,
 * and no calls on \p buffer should be made afterwards.
 *
 * \param[in] buffer  Pointer to the buffer to free.
 */
template<typename DeviceBuffer>
void freeDeviceBuffer(DeviceBuffer* buffer)
{
    GMX_ASSERT(buffer, "needs a buffer pointer");
    if (*buffer)
    {
        GMX_RELEASE_ASSERT(cudaFree(*buffer) == cudaSuccess, "cudaFree failed");
    }
}

/*! \brief
 * Performs the host-to-device data copy, synchronous or asynchronously on request.
 *
 * TODO: This is meant to gradually replace cu/ocl_copy_h2d.
 *
 * \tparam        ValueType            Raw value type of the \p buffer.
 * \param[in,out] buffer               Pointer to the device-side buffer
 * \param[in]     hostBuffer           Pointer to the raw host-side memory, also typed \p ValueType
 * \param[in]     startingOffset       Offset (in values) at the device-side buffer to copy into.
 * \param[in]     numValues            Number of values to copy.
 * \param[in]     stream               GPU stream to perform asynchronous copy in.
 * \param[in]     transferKind         Copy type: synchronous or asynchronous.
 * \param[out]    timingEvent          A dummy pointer to the H2D copy timing event to be filled in.
 *                                     Not used in CUDA implementation.
 */
template<typename ValueType>
void copyToDeviceBuffer(DeviceBuffer<ValueType>* buffer,
                        const ValueType*         hostBuffer,
                        size_t                   startingOffset,
                        size_t                   numValues,
                        CommandStream            stream,
                        GpuApiCallBehavior       transferKind,
                        CommandEvent* /*timingEvent*/)
{
    if (numValues == 0)
    {
        return; // such calls are actually made with empty domains
    }
    GMX_ASSERT(buffer, "needs a buffer pointer");
    GMX_ASSERT(hostBuffer, "needs a host buffer pointer");
    cudaError_t  stat;
    const size_t bytes = numValues * sizeof(ValueType);

    switch (transferKind)
    {
        case GpuApiCallBehavior::Async:
            GMX_ASSERT(isHostMemoryPinned(hostBuffer),
                       "Source host buffer was not pinned for CUDA");
            stat = cudaMemcpyAsync(*((ValueType**)buffer) + startingOffset, hostBuffer, bytes,
                                   cudaMemcpyHostToDevice, stream);
            GMX_RELEASE_ASSERT(stat == cudaSuccess, "Asynchronous H2D copy failed");
            break;

        case GpuApiCallBehavior::Sync:
            stat = cudaMemcpy(*((ValueType**)buffer) + startingOffset, hostBuffer, bytes,
                              cudaMemcpyHostToDevice);
            GMX_RELEASE_ASSERT(stat == cudaSuccess, "Synchronous H2D copy failed");
            break;

        default: throw;
    }
}


/*! \brief
 * Performs the device-to-host data copy, synchronous or asynchronously on request.
 *
 * TODO: This is meant to gradually replace cu/ocl_copy_d2h.
 *
 * \tparam        ValueType            Raw value type of the \p buffer.
 * \param[in,out] hostBuffer           Pointer to the raw host-side memory, also typed \p ValueType
 * \param[in]     buffer               Pointer to the device-side buffer
 * \param[in]     startingOffset       Offset (in values) at the device-side buffer to copy from.
 * \param[in]     numValues            Number of values to copy.
 * \param[in]     stream               GPU stream to perform asynchronous copy in.
 * \param[in]     transferKind         Copy type: synchronous or asynchronous.
 * \param[out]    timingEvent          A dummy pointer to the H2D copy timing event to be filled in.
 *                                     Not used in CUDA implementation.
 */
template<typename ValueType>
void copyFromDeviceBuffer(ValueType*               hostBuffer,
                          DeviceBuffer<ValueType>* buffer,
                          size_t                   startingOffset,
                          size_t                   numValues,
                          CommandStream            stream,
                          GpuApiCallBehavior       transferKind,
                          CommandEvent* /*timingEvent*/)
{
    GMX_ASSERT(buffer, "needs a buffer pointer");
    GMX_ASSERT(hostBuffer, "needs a host buffer pointer");

    cudaError_t  stat;
    const size_t bytes = numValues * sizeof(ValueType);
    switch (transferKind)
    {
        case GpuApiCallBehavior::Async:
            GMX_ASSERT(isHostMemoryPinned(hostBuffer),
                       "Destination host buffer was not pinned for CUDA");
            stat = cudaMemcpyAsync(hostBuffer, *((ValueType**)buffer) + startingOffset, bytes,
                                   cudaMemcpyDeviceToHost, stream);
            GMX_RELEASE_ASSERT(stat == cudaSuccess, "Asynchronous D2H copy failed");
            break;

        case GpuApiCallBehavior::Sync:
            stat = cudaMemcpy(hostBuffer, *((ValueType**)buffer) + startingOffset, bytes,
                              cudaMemcpyDeviceToHost);
            GMX_RELEASE_ASSERT(stat == cudaSuccess, "Synchronous D2H copy failed");
            break;

        default: throw;
    }
}

/*! \brief
 * Clears the device buffer asynchronously.
 *
 * \tparam        ValueType        Raw value type of the \p buffer.
 * \param[in,out] buffer           Pointer to the device-side buffer
 * \param[in]     startingOffset   Offset (in values) at the device-side buffer to start clearing
 * at. \param[in]     numValues        Number of values to clear. \param[in]     stream GPU stream.
 */
template<typename ValueType>
void clearDeviceBufferAsync(DeviceBuffer<ValueType>* buffer, size_t startingOffset, size_t numValues, CommandStream stream)
{
    GMX_ASSERT(buffer, "needs a buffer pointer");
    const size_t bytes   = numValues * sizeof(ValueType);
    const char   pattern = 0;

    cudaError_t stat = cudaMemsetAsync(*((ValueType**)buffer) + startingOffset, pattern, bytes, stream);
    GMX_RELEASE_ASSERT(stat == cudaSuccess, "Couldn't clear the device buffer");
}

#endif
