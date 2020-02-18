/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
 *  \brief Implements a GpuEventSynchronizer class for OpenCL.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \inlibraryapi
 */
#ifndef GMX_GPU_UTILS_GPUEVENTSYNCHRONIZER_OCL_H
#define GMX_GPU_UTILS_GPUEVENTSYNCHRONIZER_OCL_H

#include "gromacs/gpu_utils/gputraits_ocl.h"
#include "gromacs/gpu_utils/oclutils.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

/*! \libinternal \brief
 * A class which allows for CPU thread to mark and wait for certain GPU stream execution point.
 * The event can be put into the stream with markEvent() and then later waited on with waitForEvent().
 * This can be repeated as necessary, but the current implementation does not allow waiting on
 * completed event more than once, expecting only exact pairs of markEvent(stream); waitForEvent().
 * The class generally attempts to track the correctness of its state transitions, but
 * please note that calling waitForEvent() right after the construction will fail with OpenCL but succeed with CUDA.
 *
 * Another possible mode of operation can be implemented if needed:
 * multiple calls to waitForEvent() after a single markEvent(). For this, clReleaseEvent() call
 * from waitForEvent() should instead happen conditionally at the beginning of markEvent(), replacing
 * the GMX_ASSERT(). This was tested to work both with CUDA and NVidia OpenCL, but not with AMD/Intel OpenCl.
 */
class GpuEventSynchronizer
{
public:
    //! A constructor
    GpuEventSynchronizer() : event_(nullptr) {}
    //! A destructor
    ~GpuEventSynchronizer()
    {
        // This additional code only prevents cl_event leak in an unlikely situation of destructor
        // being called after markEvent() but before waitForEvent() / enqueueWaitEvent().
        if (event_)
        {
            ensureReferenceCount(event_, 1);
            clReleaseEvent(event_);
        }
    }
    //! No copying
    GpuEventSynchronizer(const GpuEventSynchronizer&) = delete;
    //! No assignment
    GpuEventSynchronizer& operator=(GpuEventSynchronizer&&) = delete;
    //! Moving is disabled but can be considered in the future if needed
    GpuEventSynchronizer(GpuEventSynchronizer&&) = delete;

    /*! \brief Marks the synchronization point in the \p stream.
     * Should be called first and then followed by waitForEvent().
     */
    inline void markEvent(const DeviceStream& deviceStream)
    {
        GMX_ASSERT(nullptr == event_, "Do not call markEvent more than once!");
        cl_int clError = clEnqueueMarkerWithWaitList(deviceStream.stream(), 0, nullptr, &event_);
        if (CL_SUCCESS != clError)
        {
            GMX_THROW(gmx::InternalError("Failed to enqueue the GPU synchronization event: "
                                         + ocl_get_error_string(clError)));
        }
    }
    /*! \brief Synchronizes the host thread on the marked event. */
    inline void waitForEvent()
    {
        cl_int clError = clWaitForEvents(1, &event_);
        if (CL_SUCCESS != clError)
        {
            GMX_THROW(gmx::InternalError("Failed to synchronize on the GPU event: "
                                         + ocl_get_error_string(clError)));
        }

        releaseEvent();
    }
    /*! \brief Enqueues a wait for the recorded event in stream \p stream
     *
     *  After enqueue, the associated event is released, so this method should
     *  be only called once per markEvent() call.
     */
    inline void enqueueWaitEvent(const DeviceStream& deviceStream)
    {
        cl_int clError = clEnqueueBarrierWithWaitList(deviceStream.stream(), 1, &event_, nullptr);
        if (CL_SUCCESS != clError)
        {
            GMX_THROW(gmx::InternalError("Failed to enqueue device barrier for the GPU event: "
                                         + ocl_get_error_string(clError)));
        }

        releaseEvent();
    }

private:
    inline void releaseEvent()
    {
        // Reference count can't be checked after the event's released, it seems (segfault on NVIDIA).
        ensureReferenceCount(event_, 1);
        cl_int clError = clReleaseEvent(event_);
        if (CL_SUCCESS != clError)
        {
            GMX_THROW(gmx::InternalError("Failed to release the GPU event: "
                                         + ocl_get_error_string(clError)));
        }
        event_ = nullptr;
    }

    cl_event event_;
};

#endif
