/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
 *
 * \brief Implements the DeviceStream for CUDA.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_gpu_utils
 */
#include "gmxpre.h"

#include "device_stream.h"

#include "gromacs/gpu_utils/gputraits.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

DeviceStream::DeviceStream()
{
    stream_ = nullptr;
}

void DeviceStream::init(const DeviceContext& /* deviceContext */,
                        DeviceStreamPriority priority,
                        const bool /* useTiming */)
{
    cudaError_t stat;

    if (priority == DeviceStreamPriority::Normal)
    {
        stat = cudaStreamCreate(&stream_);
        if (stat != cudaSuccess)
        {
            GMX_THROW(gmx::InternalError(gmx::formatString(
                    "Could not create CUDA stream (CUDA error %d: %s).", stat, cudaGetErrorString(stat))));
        }
    }
    else if (priority == DeviceStreamPriority::High)
    {
        // Note that the device we're running on does not have to
        // support priorities, because we are querying the priority
        // range, which in that case will be a single value.
        int highestPriority;
        stat = cudaDeviceGetStreamPriorityRange(nullptr, &highestPriority);
        if (stat != cudaSuccess)
        {
            GMX_THROW(gmx::InternalError(gmx::formatString(
                    "Could not query CUDA stream priority range (CUDA error %d: %s).", stat,
                    cudaGetErrorString(stat))));
        }

        stat = cudaStreamCreateWithPriority(&stream_, cudaStreamDefault, highestPriority);
        if (stat != cudaSuccess)
        {
            GMX_THROW(gmx::InternalError(gmx::formatString(
                    "Could not create CUDA stream with high priority (CUDA error %d: %s).", stat,
                    cudaGetErrorString(stat))));
        }
    }
}

DeviceStream::~DeviceStream()
{
    if (isValid())
    {
        cudaError_t stat = cudaStreamDestroy(stream_);
        GMX_RELEASE_ASSERT(stat == cudaSuccess,
                           gmx::formatString("Failed to release CUDA stream (CUDA error %d: %s).",
                                             stat, cudaGetErrorString(stat))
                                   .c_str());
        stream_ = nullptr;
    }
}

cudaStream_t DeviceStream::stream() const
{
    return stream_;
}

bool DeviceStream::isValid() const
{
    return (stream_ != nullptr);
}

void DeviceStream::synchronize() const
{
    cudaError_t stat = cudaStreamSynchronize(stream_);
    GMX_RELEASE_ASSERT(stat == cudaSuccess,
                       gmx::formatString("cudaStreamSynchronize failed  (CUDA error %d: %s).", stat,
                                         cudaGetErrorString(stat))
                               .c_str());
}