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
/*! \libinternal \file
 *  \brief Implements a DeviceEvent class for HIP.
 *
 *  \author Paul Bauer <paul.bauer.q@gmail.com>
 *  \author Julio Maia <julio.maia@amd.com>
 *  \inlibraryapi
 */
#ifndef GMX_GPU_UTILS_DEVICE_EVENT_HIP_H
#define GMX_GPU_UTILS_DEVICE_EVENT_HIP_H

#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/hiputils.h"

#ifndef DOXYGEN

class DeviceEvent
{
public:
    DeviceEvent()
    {
        hipError_t stat = hipEventCreateWithFlags(&event_, hipEventDisableTiming);
        if (stat != hipSuccess)
        {
            GMX_THROW(gmx::InternalError("hipEventCreate failed: " + gmx::getDeviceErrorString(stat)));
        }
    }
    ~DeviceEvent() { std::ignore = hipEventDestroy(event_); }
    // Disable copy, move, and assignment. Move can be allowed, but not needed yet.
    DeviceEvent& operator=(const DeviceEvent&) = delete;
    DeviceEvent(const DeviceEvent&)            = delete;
    DeviceEvent& operator=(DeviceEvent&&) = delete;
    DeviceEvent(DeviceEvent&&)            = delete;

    /*! \brief Marks the synchronization point in the \p stream.
     * Should be followed by waitForEvent().
     */
    inline void mark(const DeviceStream& deviceStream)
    {
        hipError_t stat = hipEventRecord(event_, deviceStream.stream());
        if (stat != hipSuccess)
        {
            GMX_THROW(gmx::InternalError("hipEventRecord failed: " + gmx::getDeviceErrorString(stat)));
        }
        isMarked_ = true;
    }
    //! Synchronizes the host thread on the marked event.
    inline void wait()
    {
        hipError_t gmx_used_in_debug stat = hipEventSynchronize(event_);
        if (stat != hipSuccess)
        {
            GMX_THROW(gmx::InternalError("hipEventSynchronize failed: " + gmx::getDeviceErrorString(stat)));
        }
    }
    //! Checks the completion of the underlying event.
    inline bool isReady()
    {
        hipError_t stat = hipEventQuery(event_);
        if (stat != hipSuccess && stat != hipErrorNotReady)
        {
            GMX_THROW(gmx::InternalError("hipEventQuery failed: " + gmx::getDeviceErrorString(stat)));
        }
        return (stat == hipSuccess);
    }
    //! Check if this event was marked
    inline bool isMarked() const { return isMarked_; }
    //! Enqueues a wait for the recorded event in stream \p stream
    inline void enqueueWait(const DeviceStream& deviceStream)
    {
        hipError_t stat = hipStreamWaitEvent(deviceStream.stream(), event_, 0);
        if (stat != hipSuccess)
        {
            GMX_THROW(gmx::InternalError("hipStreamWaitEvent failed: " + gmx::getDeviceErrorString(stat)));
        }
    }
    //! Reset the event
    inline void reset() { isMarked_ = false; }

private:
    hipEvent_t event_;
    bool       isMarked_ = false;
};

#endif

#endif
