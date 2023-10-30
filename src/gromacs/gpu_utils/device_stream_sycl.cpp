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
 * \brief Implements the DeviceStream for SYCL builds.
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \author Andrey Alekseenko <al42and@gmail.com>
 *
 * \ingroup module_gpu_utils
 */
#include "gmxpre.h"

#include <cstdio>

#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/utility/exceptions.h"


static sycl::property_list makeQueuePropertyList(bool enableProfiling, DeviceStreamPriority priority)
{
    // If hipSYCL priority extension is present, extract the priority range
    // and use the lowest or highest priority supported for DeviceStreamPriority::Low and
    // DeviceStreamPriority::High, respectively.
#ifdef HIPSYCL_EXT_QUEUE_PRIORITY
    // For simplicity, we assume 0 to be the default priority (guaranteed for CUDA, verified for HIP)
    int defaultPrioValue = 0;
    int highPrioValue    = 0;

#    if GMX_HIPSYCL_HAVE_CUDA_TARGET || GMX_HIPSYCL_HAVE_HIP_TARGET
#        if GMX_HIPSYCL_HAVE_CUDA_TARGET
    const auto status         = cudaDeviceGetStreamPriorityRange(nullptr, &highPrioValue);
    const auto apiSuccessCode = cudaSuccess;
#        elif GMX_HIPSYCL_HAVE_HIP_TARGET
    auto       status         = hipDeviceGetStreamPriorityRange(nullptr, &highPrioValue);
    const auto apiSuccessCode = hipSuccess;
#        endif
    if (status != apiSuccessCode)
    {
        GMX_THROW(gmx::InternalError("[cuda|hip]DeviceGetStreamPriorityRange failed"));
    }
#    endif

#    define HIPSYCL_PRIORITY_ATTRIBUTE_HIGH \
        sycl::property::queue::hipSYCL_priority { highPrioValue }
#    define HIPSYCL_PRIORITY_ATTRIBUTE_DEFAULT \
        sycl::property::queue::hipSYCL_priority { defaultPrioValue }
#elif defined(SYCL_EXT_ONEAPI_QUEUE_PRIORITY)
#    define HIPSYCL_PRIORITY_ATTRIBUTE_HIGH \
        sycl::ext::oneapi::property::queue::priority_high {}
#    define HIPSYCL_PRIORITY_ATTRIBUTE_DEFAULT \
        sycl::ext::oneapi::property::queue::priority_normal {}
#else
#    define HIPSYCL_PRIORITY_ATTRIBUTE_HIGH
#    define HIPSYCL_PRIORITY_ATTRIBUTE_DEFAULT
#endif

    if (enableProfiling)
    {
        if (priority == DeviceStreamPriority::High)
        {
            return { sycl::property::queue::in_order(),
                     sycl::property::queue::enable_profiling(),
                     HIPSYCL_PRIORITY_ATTRIBUTE_HIGH };
        }
        else
        {
            return { sycl::property::queue::in_order(),
                     sycl::property::queue::enable_profiling(),
                     HIPSYCL_PRIORITY_ATTRIBUTE_DEFAULT };
        }
    }
    else // !enableProfiling
    {
        if (priority == DeviceStreamPriority::High)
        {
            return { sycl::property::queue::in_order(), HIPSYCL_PRIORITY_ATTRIBUTE_HIGH };
        }
        else
        {
            return { sycl::property::queue::in_order(), HIPSYCL_PRIORITY_ATTRIBUTE_DEFAULT };
        }
    }
#undef HIPSYCL_PRIORITY_ATTRIBUTE_HIGH
#undef HIPSYCL_PRIORITY_ATTRIBUTE_DEFAULT
}

DeviceStream::DeviceStream(const DeviceContext& deviceContext, DeviceStreamPriority priority, const bool useTiming)
{
    const std::vector<sycl::device> devicesInContext = deviceContext.context().get_devices();
    // The context is constructed to have exactly one device
    const sycl::device device = devicesInContext[0];

    bool enableProfiling = false;
    if (useTiming)
    {
        const bool deviceSupportsTiming = device.has(sycl::aspect::queue_profiling);
        enableProfiling                 = deviceSupportsTiming;
    }
    stream_ = sycl::queue(
            deviceContext.context(), device, makeQueuePropertyList(enableProfiling, priority));
}

DeviceStream::~DeviceStream()
{
#if GMX_SYCL_HIPSYCL
    // Prevents use-after-free errors in hipSYCL's CUDA backend during unit tests
    try
    {
        synchronize();
    }
    catch (sycl::exception& e)
    {
        std::fprintf(stderr, "Error in SYCL queue: %s\n", e.what());
    }
#endif
}

// NOLINTNEXTLINE readability-convert-member-functions-to-static
bool DeviceStream::isValid() const
{
    return true;
}

void DeviceStream::synchronize()
{
    stream_.wait_and_throw();
};

void DeviceStream::synchronize() const
{
    /* sycl::queue::wait is a non-const function. However, a lot of code in GROMACS
     * assumes DeviceStream is const, yet wants to synchronize with it.
     * The chapter "4.3.2 Common reference semantics" of SYCL 1.2.1 specification says:
     * > Each of the following SYCL runtime classes: [...] queue, [...] must obey the following
     * > statements, where T is the runtime class type:
     * > - T must be copy constructible and copy assignable on the host application [...].
     * >   Any instance of T that is constructed as a copy of another instance, via either the
     * >   copy constructor or copy assignment operator, must behave as-if it were the original
     * >   instance and as-if any action performed on it were also performed on the original
     * >   instance [...].
     * Same in chapter "4.5.3" of provisional SYCL 2020 specification (June 30, 2020).
     * So, we can copy-construct a new queue and wait() on it.
     */
    sycl::queue(stream_).wait_and_throw();
}

void issueClFlushInStream(const DeviceStream& /*deviceStream*/) {}
