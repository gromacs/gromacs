/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2021, by the GROMACS development team, led by
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
 *  \brief Implements a GpuEventSynchronizer class.
 *
 *  \author Andrey Alekseenko <al42and@gmail.com>
 *  \author Artem Zhmurov <zhmurov@gmail.com>
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \inlibraryapi
 */
#ifndef GMX_GPU_UTILS_GPUEVENTSYNCHRONIZER_H
#define GMX_GPU_UTILS_GPUEVENTSYNCHRONIZER_H

#include "config.h"

#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

#include "device_event.h"

/*! \libinternal \brief
 * A class which allows for CPU thread to mark and wait for certain GPU stream execution point.
 * The event can be put into the stream with \ref markEvent() and then later waited on with \ref waitForEvent().
 * This can be repeated as necessary, but the current implementation does not allow waiting on
 * completed event more than once, expecting only exact pairs of markEvent(stream); waitForEvent().
 * The class generally attempts to track the correctness of its state transitions, but
 * please note that calling waitForEvent() right after the construction will fail with OpenCL but succeed with CUDA.
 */
class GpuEventSynchronizer
{
public:
    //! A constructor
    GpuEventSynchronizer() = default;
    //! A destructor
    ~GpuEventSynchronizer() = default;
    //! Remove copy assignment, because we can not copy the underlying event object.
    GpuEventSynchronizer& operator=(const GpuEventSynchronizer&) = delete;
    //! Remove copy constructor, because we can not copy the underlying event object.
    GpuEventSynchronizer(const GpuEventSynchronizer&) = delete;
    //! Remove move assignment, because we don't allow moving the underlying event object.
    GpuEventSynchronizer& operator=(GpuEventSynchronizer&&) = delete;
    //! Remove move constructor, because we don't allow moving the underlying event object.
    GpuEventSynchronizer(GpuEventSynchronizer&&) = delete;

    /*! \brief Marks the synchronization point in the \p stream.
     * Should be called first and then followed by \ref waitForEvent().
     */
    inline void markEvent(const DeviceStream& deviceStream)
    {
#if !GMX_GPU_CUDA // For now, we have relaxed conditions for CUDA
        if (event_.isMarked())
        {
            GMX_THROW(gmx::InternalError("Trying to mark event before first consuming it"));
        }
#endif
        event_.mark(deviceStream);
    }
    /*! \brief Synchronizes the host thread on the marked event. */
    inline void waitForEvent()
    {
#if !GMX_GPU_CUDA // For now, we have relaxed conditions for CUDA
        if (!event_.isMarked())
        {
            GMX_THROW(gmx::InternalError(
                    "Trying to wait for event before marking it or after fully consuming it"));
        }
#endif
        event_.wait();
        reset();
    }
    /*! \brief Checks the completion of the underlying event and resets the object if it was. */
    inline bool isReady()
    {
#if !GMX_GPU_CUDA // For now, we have relaxed conditions for CUDA
        if (!event_.isMarked())
        {
            GMX_THROW(gmx::InternalError("Trying to check the status of event before marking it"));
        }
#endif
        bool isReady = event_.isReady();
        if (isReady)
        {
            reset();
        }
        return isReady;
    }
    /*! \brief Enqueues a wait for the recorded event in stream \p stream
     *
     *  After enqueue, the associated event is released, so this method should
     *  be only called once per \ref markEvent() call (not enforced in CUDA yet).
     */
    inline void enqueueWaitEvent(const DeviceStream& deviceStream)
    {
#if !GMX_GPU_CUDA // For now, we have relaxed conditions for CUDA
        if (!event_.isMarked())
        {
            GMX_THROW(
                    gmx::InternalError("Trying to enqueue wait for event before marking it or "
                                       "after fully consuming it"));
        }
#endif
        event_.enqueueWait(deviceStream);
        reset();
    }

    //! Resets the event to unmarked state, releasing the underlying event object if needed.
    inline void reset() { event_.reset(); }

private:
    DeviceEvent event_;
};

#endif
