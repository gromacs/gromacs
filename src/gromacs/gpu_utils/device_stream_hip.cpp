/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * \brief Implements the DeviceStream for HIP.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \author Julio Maia <julio.maia@amd.com>
 *
 * \ingroup module_gpu_utils
 */
#include "gmxpre.h"

#include "gromacs/gpu_utils/hiputils.h"

#include "device_stream.h"

DeviceStream::DeviceStream(const DeviceContext& /* deviceContext */,
                           DeviceStreamPriority priority,
                           const bool /* useTiming */)
{
    hipError_t stat;

    if (priority == DeviceStreamPriority::Normal)
    {
        stat = hipStreamCreate(&stream_);
        gmx::checkDeviceError(stat, "Could not create HIP stream.");
    }
    else if (priority == DeviceStreamPriority::High)
    {
        // Note that the device we're running on does not have to
        // support priorities, because we are querying the priority
        // range, which in that case will be a single value.
        int highestPriority;
        stat = hipDeviceGetStreamPriorityRange(nullptr, &highestPriority);
        gmx::checkDeviceError(stat, "Could not query HIP stream priority range.");

        stat = hipStreamCreateWithPriority(&stream_, hipStreamDefault, highestPriority);
        gmx::checkDeviceError(stat, "Could not create HIP stream with high priority.");
    }
}

DeviceStream::~DeviceStream()
{
    if (isValid())
    {
        hipError_t stat = hipStreamDestroy(stream_);
        if (stat != hipSuccess)
        {
            // Don't throw in the destructor, just print a warning
            std::fprintf(stderr,
                         "Failed to release HIP stream. %s\n",
                         gmx::getDeviceErrorString(stat).c_str());
        }
        stream_ = nullptr;
    }
}

hipStream_t DeviceStream::stream() const
{
    return stream_;
}

bool DeviceStream::isValid() const
{
    return (stream_ != nullptr);
}

void DeviceStream::synchronize() const
{
    hipError_t stat = hipStreamSynchronize(stream_);
    gmx::checkDeviceError(stat, "hipStreamSynchronize failed");
}

void issueClFlushInStream(const DeviceStream& /*deviceStream*/) {}
