/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
#ifndef DEVICEBUFFER_CUH
#define DEVICEBUFFER_CUH

#include "gromacs/gpu_utils/gputraits.cuh"
#include "gromacs/gpu_utils/gpu_utils.h" //only for GpuApiCallBehavior
#include "gromacs/utility/gmxassert.h"

//! \brief A device-side buffer of ValueTypes
template<typename ValueType>
using DeviceBuffer = ValueType *;

template <typename ValueType>
void allocateDeviceBuffer(DeviceBuffer<ValueType> *buffer,
                          size_t                   numValues,
                          Context                  context)
{
    GMX_ASSERT(buffer, "needs a buffer pointer");
    GMX_UNUSED_VALUE(context); // not used explicitly in CUDA RT
    cudaError_t stat = cudaMalloc((void **)buffer, numValues * sizeof(ValueType));
    GMX_RELEASE_ASSERT(stat == cudaSuccess, "cudaMalloc failure");
}

/*! \brief Free a device-side buffer.
 * This does not reset separately stored size/capacity integers,
 * as this is planned to be a destructor of DeviceBuffer as a proper class,
 * and no calls on \p buffer should be made afterwards.
 *
 * \param[in] buffer  Pointer to the buffer to free.
 */
template <typename DeviceBuffer>
void freeDeviceBuffer(DeviceBuffer *buffer)
{
    GMX_ASSERT(buffer, "needs a buffer pointer");
    if (*buffer)
    {
        GMX_RELEASE_ASSERT(cudaFree(*buffer) == cudaSuccess, "cudaFree failed");
    }
}

template <typename ValueType>
void copyToDeviceBuffer(DeviceBuffer<ValueType> *buffer,
                        const ValueType         *hostBuffer,
                        size_t                   startingValueIndex,
                        size_t                   numValues,
                        CommandStream            stream,
                        GpuApiCallBehavior       transferKind,
                        CommandEvent            *timingEvent)
{
    if (numValues == 0)
    {
        return; // such calls are actually made with empty domains
    }
    GMX_ASSERT(buffer, "needs a buffer pointer");
    GMX_ASSERT(hostBuffer, "needs a host buffer pointer");
    GMX_UNUSED_VALUE(timingEvent); // not applicable in CUDA
    cudaError_t  stat;
    const size_t bytes = numValues * sizeof(ValueType);

    switch (transferKind)
    {
        case GpuApiCallBehavior::Async:
            GMX_ASSERT(isHostMemoryPinned(hostBuffer), "Source host buffer was not pinned for CUDA");
            stat = cudaMemcpyAsync(*((ValueType **)buffer) + startingValueIndex, hostBuffer, bytes, cudaMemcpyHostToDevice, stream);
            GMX_RELEASE_ASSERT(stat == cudaSuccess, "Asynchronous H2D copy failed");
            break;

        case GpuApiCallBehavior::Sync:
            stat = cudaMemcpy(*((ValueType **)buffer) + startingValueIndex, hostBuffer, bytes, cudaMemcpyHostToDevice); //TODO is this type dependency the best?
            GMX_RELEASE_ASSERT(stat == cudaSuccess, "Synchronous H2D copy failed");
            break;

        default:
            throw;
    }
}

#include "devicebuffer.h"

#endif // DEVICEBUFFER_CUH
