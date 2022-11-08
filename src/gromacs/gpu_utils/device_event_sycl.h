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
 *  This implementation relies on SYCL_INTEL_enqueue_barrier proposal,
 *  https://github.com/intel/llvm/blob/sycl/sycl/doc/extensions/EnqueueBarrier/enqueue_barrier.asciidoc
 *
 *  Using event-based synchronization is not recommended for SYCL.
 *  SYCL queues are out-of-order and rely on data dependencies, allowing only to wait
 *  for a specific kernel (by capturing the \c event returned from \c queue.submit) or for all
 *  the tasks in the queue (\c queue.wait).
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
#    if GMX_SYCL_HIPSYCL
        // This will not launch any GPU operation, but it will mark an event which is returned;
        // it is functionally equivalent with ext_oneapi_submit_barrier().
        events_ = { deviceStream.stream().hipSYCL_enqueue_custom_operation([=](sycl::interop_handle&) {}) };
        isMarked_ = true;
#    else
        // Relies on SYCL_INTEL_enqueue_barrier
#        if __SYCL_COMPILER_VERSION >= 20211123
        events_ = { deviceStream.stream().ext_oneapi_submit_barrier() };
#        else
        events_ = { deviceStream.stream().submit_barrier() };
#        endif
#    endif
    }

    //! Synchronizes the host thread on the marked event.
    inline void wait()
    {
#    if GMX_SYCL_DPCPP
        // Note: this is not to prevent use-before-marking, but for checking the DPC++ vs hipSYCL consistency
        GMX_ASSERT(events_.size() <= 1, "One event expected in DPC++, but we have several!");
#    endif
        for (auto& event : events_)
        {
            event.wait_and_throw();
        }
    }

    inline void enqueueWait(const DeviceStream& deviceStream)
    {
#    if GMX_SYCL_HIPSYCL
        // Submit an empty operation that depends on all the events recorded.
        deviceStream.stream().submit(GMX_SYCL_DISCARD_EVENT[&](sycl::handler & cgh) {
            cgh.depends_on(events_);
            cgh.hipSYCL_enqueue_custom_operation([=](sycl::interop_handle&) {});
        });
#    else
        GMX_ASSERT(events_.size() <= 1, "One event expected in DPC++, but we have several!");
        // Relies on SYCL_INTEL_enqueue_barrier extensions
#        if __SYCL_COMPILER_VERSION >= 20211123
        deviceStream.stream().ext_oneapi_submit_barrier(events_);
#        else
        deviceStream.stream().submit_barrier(events_);
#        endif
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
    inline bool isMarked() const
    {
#    if GMX_SYCL_HIPSYCL
        return isMarked_;
#    else
        return !events_.empty();
#    endif
    }

    //! Reset the event to unmarked state.
    inline void reset()
    {
        events_.clear();
#    if GMX_SYCL_HIPSYCL
        isMarked_ = false;
#    endif
    }

private:
    std::vector<sycl::event> events_;
#    if GMX_SYCL_HIPSYCL
    /*! \brief Flag to track event marking in hipSYCL.
     *
     * In hipSYCL, we can have empty \ref events_ after marking if there were no pending tasks in
     * the queue. So, we use an explicit flag to check the event state. */
    bool isMarked_ = false;
#    endif
};

#endif // DOXYGEN

#endif // GMX_GPU_UTILS_DEVICE_EVENT_SYCL_H
