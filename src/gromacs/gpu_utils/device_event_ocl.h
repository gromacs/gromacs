/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 *  \brief Implements a DeviceEvent class for OpenCL.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 *  \author Andrey Alekseenko <al42and@gmail.com>
 * \inlibraryapi
 */
#ifndef GMX_GPU_UTILS_DEVICE_EVENT_OCL_H
#define GMX_GPU_UTILS_DEVICE_EVENT_OCL_H

#include "gromacs/gpu_utils/gputraits_ocl.h"
#include "gromacs/gpu_utils/oclutils.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

#ifndef DOXYGEN

class DeviceEvent
{
public:
    //! A constructor
    DeviceEvent() : event_(sc_nullEvent) {}
    DeviceEvent(cl_event event) : event_(event) {}
    //! A destructor
    ~DeviceEvent()
    {
        if (isMarked())
        {
            // Can not throw in destructor, so not checking for any error
            clReleaseEvent(event_);
        }
    }
    // Disable copy, move, and assignment. Move can be allowed, but not needed yet.
    DeviceEvent& operator=(const DeviceEvent&) = delete;
    DeviceEvent(const DeviceEvent&)            = delete;
    DeviceEvent& operator=(DeviceEvent&&) = delete;
    DeviceEvent(DeviceEvent&&)            = delete;

    /*! \brief Marks the synchronization point in the \p stream.
     * Should be called first and then followed by wait().
     */
    inline void mark(const DeviceStream& deviceStream)
    {
        reset();
        cl_int clError = clEnqueueMarkerWithWaitList(deviceStream.stream(), 0, nullptr, &event_);
        if (CL_SUCCESS != clError)
        {
            GMX_THROW(gmx::InternalError("Failed to enqueue the GPU synchronization event: "
                                         + ocl_get_error_string(clError)));
        }
    }

    /*! \brief Synchronizes the host thread on the marked event. */
    inline void wait()
    {
        GMX_RELEASE_ASSERT(isMarked(), "Can not wait for an unmarked event");
        cl_int clError = clWaitForEvents(1, &event_);
        if (CL_SUCCESS != clError)
        {
            GMX_THROW(gmx::InternalError("Failed to synchronize on the GPU event: "
                                         + ocl_get_error_string(clError)));
        }
    }

    /*! \brief Enqueues a wait for the recorded event in stream \p stream. */
    inline void enqueueWait(const DeviceStream& deviceStream)
    {
        GMX_RELEASE_ASSERT(isMarked(), "Can not enqueue an unmarked event");
        cl_int clError = clEnqueueBarrierWithWaitList(deviceStream.stream(), 1, &event_, nullptr);
        if (CL_SUCCESS != clError)
        {
            GMX_THROW(gmx::InternalError("Failed to enqueue device barrier for the GPU event: "
                                         + ocl_get_error_string(clError)));
        }
    }

    //!  Checks the completion of the underlying event.
    inline bool isReady()
    {
        GMX_RELEASE_ASSERT(isMarked(), "Can not check the status of unmarked event");
        cl_int result;
        cl_int clError = clGetEventInfo(
                event_, CL_EVENT_COMMAND_EXECUTION_STATUS, sizeof(cl_int), &result, nullptr);
        if (CL_SUCCESS != clError)
        {
            GMX_THROW(gmx::InternalError("Failed to retrieve event info: " + ocl_get_error_string(clError)));
        }
        return (result == CL_COMPLETE);
    }

    //! Checks whether this object encapsulates an underlying event.
    inline bool isMarked() const { return event_ != sc_nullEvent; }

    //! Reset (release) the event to unmarked state.
    inline void reset()
    {
        if (isMarked())
        {
            cl_int clError = clReleaseEvent(event_);
            if (CL_SUCCESS != clError)
            {
                GMX_THROW(gmx::InternalError("Failed to release the GPU event: "
                                             + ocl_get_error_string(clError)));
            }
        }
        event_ = sc_nullEvent;
    }

private:
    cl_event event_;

    //! Magic value to indicate uninitialized state.
    static constexpr cl_event sc_nullEvent = nullptr;
};

#endif

#endif
