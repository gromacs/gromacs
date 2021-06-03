/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020,2021, by the GROMACS development team, led by
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
#ifndef GMX_GPU_UTILS_GPUEVENTSYNCHRONIZER_SYCL_H
#define GMX_GPU_UTILS_GPUEVENTSYNCHRONIZER_SYCL_H

#include <optional>

#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

#ifndef DOXYGEN
/*! \libinternal \brief
 * A class which allows for CPU thread to mark and wait for certain GPU stream execution point.
 * The event can be put into the stream with markEvent() and then later waited on with waitForEvent().
 * This can be repeated as necessary, but the current implementation does not allow waiting on
 * completed event more than once, expecting only exact pairs of markEvent(stream); waitForEvent().
 * The class generally attempts to track the correctness of its state transitions, but
 * please note that calling waitForEvent() right after the construction will fail with OpenCL
 * and SYCL but succeed with CUDA.
 *
 * Another possible mode of operation can be implemented if needed:
 * multiple calls to waitForEvent() after a single markEvent(). For this, event.reset() call
 * from waitForEvent() should instead happen conditionally at the beginning of markEvent(), replacing
 * the GMX_ASSERT(). This was tested to work both with CUDA, NVidia OpenCL, and Intel SYCL,
 * but not with AMD/Intel OpenCl.
 *
 *  \warning This class is offered for uniformity with other GPU implementations, but expect it to
 *  be deprecated in the future.
 *
 */
class GpuEventSynchronizer
{
public:
    //! A constructor.
    GpuEventSynchronizer()
    {
        doNotSynchronizeBetweenStreams_ = (std::getenv("GMX_GPU_SYCL_NO_SYNCHRONIZE") != nullptr);
        events_.reserve(1);
    }
    //! A constructor from an existing event.
    GpuEventSynchronizer(const cl::sycl::event& event) : events_{ event } {}
    //! A destructor.
    ~GpuEventSynchronizer() = default;
    //! No copying
    GpuEventSynchronizer(const GpuEventSynchronizer&) = delete;
    //! No assignment
    GpuEventSynchronizer& operator=(GpuEventSynchronizer&&) = delete;
    //! Moving is disabled but can be considered in the future if needed
    GpuEventSynchronizer(GpuEventSynchronizer&&) = delete;

    /*! \brief Marks the synchronization point in the \p deviceStream.
     * Should be called first and then followed by waitForEvent() or enqueueWaitEvent().
     */
    inline void markEvent(const DeviceStream& deviceStream)
    {
        GMX_ASSERT(!isMarked(), "Do not call markEvent more than once!");
#    if GMX_SYCL_HIPSYCL
        // Relies on HIPSYCL_EXT_QUEUE_WAIT_LIST extension
        events_ = deviceStream.stream().get_wait_list();
#    else
        // Relies on SYCL_INTEL_enqueue_barrier
        events_ = { deviceStream.stream().submit_barrier() };
#    endif
    }
    /*! \brief Synchronizes the host thread on the marked event.
     * As in the OpenCL implementation, the event is released.
     */
    inline void waitForEvent()
    {
        GMX_ASSERT(isMarked(), "Don't call waitForEvent before marking the event!");
#    if GMX_SYCL_DPCPP
        GMX_ASSERT(events_.size() == 1, "One event expected in DPCPP, but we have several!");
#    endif
        for (auto& event : events_)
        {
            event.wait_and_throw();
        }
        reset();
    }
    /*! \brief Checks the completion of the underlying event and resets the object if it was. */
    inline bool isReady()
    {
        bool allReady = std::all_of(events_.begin(), events_.end(), [](cl::sycl::event& event) {
            auto info       = event.get_info<cl::sycl::info::event::command_execution_status>();
            bool isComplete = (info == cl::sycl::info::event_command_status::complete);
            return isComplete;
        });
        if (allReady)
        {
            reset();
        }
        return allReady;
    }
    /*! \brief Enqueues a wait for the recorded event in stream \p deviceStream.
     * As in the OpenCL implementation, the event is released.
     */
    inline void enqueueWaitEvent(const DeviceStream& deviceStream)
    {
        if (!doNotSynchronizeBetweenStreams_)
        {
#    if GMX_SYCL_HIPSYCL
            // Submit an empty kernel that depends on all the events recorded.
            deviceStream.stream().single_task(events_, [=]() {});
#    else
            // Relies on SYCL_INTEL_enqueue_barrier extensions
            GMX_ASSERT(events_.size() == 1, "Only one event expected in DPCPP!");
            deviceStream.stream().submit_barrier(events_);
#    endif
        }
        reset();
    }
    //! Reset the event to unmarked state.
    inline void reset() { events_.clear(); }
    //! Check if the event is marked. Needed for some workarounds for #3988
    inline bool isMarked() const { return !events_.empty(); }

private:
    std::vector<cl::sycl::event> events_;
    /*! \brief Dev. setting to no-op enqueueWaitEvent
     *
     * In SYCL, dependencies between the GPU tasks are managed by the runtime, so manual
     * synchronization between GPU streams should be redundant, but we keep it on by default.
     *
     * Setting this to \c true via \c GMX_GPU_SYCL_NO_SYNCHRONIZE environment variable will
     * immediately return from \ref enqueueWaitEvent, without placing a barrier into the stream.
     */
    bool doNotSynchronizeBetweenStreams_;
};

#endif // !defined DOXYGEN

#endif // GMX_GPU_UTILS_GPUEVENTSYNCHRONIZER_SYCL_H
