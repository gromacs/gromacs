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
/*! \libinternal \file
 *  \brief Implements a GpuEventSynchronizer class for SYCL.
 *
 *  There is no way to do this efficiently within the SYCL standard.
 *  With DPC++, this implementation relies on SYCL_EXT_ONEAPI_ENQUEUE_BARRIER extension,
 *  which directly enqueues a CUDA-like barrier event.
 *  With AdaptiveCpp, it relies on ACPP_EXT_ENQUEUE_CUSTOM_OPERATION extension to submit
 *  an empty operation (that can have dependencies and returns an event).
 *
 *  \author Erik Lindahl <erik.lindahl@gmail.com>
 *  \author Andrey Alekseenko <al42and@gmail.com>
 * \inlibraryapi
 */
#ifndef GMX_GPU_UTILS_DEVICE_EVENT_SYCL_H
#define GMX_GPU_UTILS_DEVICE_EVENT_SYCL_H

#include <algorithm>
#include <vector>

#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

#ifndef DOXYGEN

class DeviceEvent
{
public:
    //! A constructor.
    DeviceEvent() { events_.reserve(1); }
    //! A constructor from an existing event.
    DeviceEvent(const sycl::event& event) : events_{ event } {}
    //! A destructor.
    ~DeviceEvent() = default;
    // Disable copy, move, and assignment. They all can be allowed, but not needed yet.
    DeviceEvent& operator=(const DeviceEvent&) = delete;
    DeviceEvent(const DeviceEvent&)            = delete;
    DeviceEvent& operator=(DeviceEvent&&) = delete;
    DeviceEvent(DeviceEvent&&)            = delete;

    /*! \brief Marks the synchronization point in the \p deviceStream.
     * Should be called first and then followed by wait() or enqueueWait().
     */
    inline void mark(const DeviceStream& deviceStream)
    {
#    if defined(ACPP_EXT_ENQUEUE_CUSTOM_OPERATION) || defined(HIPSYCL_EXT_ENQUEUE_CUSTOM_OPERATION)
        // This will not launch any GPU operation, but it will mark an event which is returned
        events_ = { deviceStream.stream().hipSYCL_enqueue_custom_operation([=](sycl::interop_handle&) {}) };
#    elif defined(SYCL_EXT_ONEAPI_ENQUEUE_BARRIER)
        events_ = { deviceStream.stream().ext_oneapi_submit_barrier() };
#    else
        // We can do full stream synchronization here, but it's so inefficient we better bail out
#        error "Either ACPP_EXT_ENQUEUE_CUSTOM_OPERATION or SYCL_EXT_ONEAPI_ENQUEUE_BARRIER is needed"
#    endif
    }

    //! Synchronizes the host thread on the marked event.
    inline void wait()
    {
        // Note: this is not to prevent use-before-marking, but for checking the DPC++ vs hipSYCL consistency
        for (auto& event : events_)
        {
            event.wait_and_throw();
        }
    }

    inline void enqueueWait(const DeviceStream& deviceStream)
    {
#    if defined(ACPP_EXT_ENQUEUE_CUSTOM_OPERATION) || defined(HIPSYCL_EXT_ENQUEUE_CUSTOM_OPERATION)
        // Submit an empty operation that depends on all the events recorded.
        deviceStream.stream().submit(GMX_SYCL_DISCARD_EVENT[&](sycl::handler & cgh) {
            cgh.depends_on(events_);
            cgh.hipSYCL_enqueue_custom_operation([=](sycl::interop_handle&) {});
        });
#    elif defined(SYCL_EXT_ONEAPI_ENQUEUE_BARRIER)
        // Relies on sycl_ext_oneapi_enqueue_barrier extensions
        deviceStream.stream().ext_oneapi_submit_barrier(events_);
#    else
        // In over-synchronization case, we would do nothing, since we would have already synced
#        error "Either ACPP_EXT_ENQUEUE_CUSTOM_OPERATION or SYCL_EXT_ONEAPI_ENQUEUE_BARRIER is needed"
#    endif
    }

    //! Checks the completion of the underlying event.
    inline bool isReady()
    {
        bool allReady = std::all_of(events_.begin(), events_.end(), [](sycl::event& event) {
            auto info       = event.get_info<sycl::info::event::command_execution_status>();
            bool isComplete = (info == sycl::info::event_command_status::complete);
            return isComplete;
        });
        return allReady;
    }

    //! Checks whether this object encapsulates an underlying event.
    inline bool isMarked() const { return !events_.empty(); }

    //! Reset the event to unmarked state.
    inline void reset() { events_.clear(); }

private:
    /* We can only ever have 0 or 1 event, but the Intel extension uses std::vector as
     * a function argument, so we use it instead of an optional. */
    std::vector<sycl::event> events_;
};

#endif // DOXYGEN

#endif // GMX_GPU_UTILS_DEVICE_EVENT_SYCL_H
