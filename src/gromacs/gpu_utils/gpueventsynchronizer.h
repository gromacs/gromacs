/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 *
 * The event can be put into the stream with \ref markEvent and then later waited on with \ref
 * waitForEvent or \ref enqueueWaitEvent.
 *
 * Additionally, this class offers facilities for runtime checking of correctness by counting
 * how many times each marked event is used as a synchronization point.
 *
 * - When the class is constructed, a required minimal (\c minConsumptionCount) and maximal (\c maxConsumptionCount) number of
 * consumptions can be specified. By default, both are set to 1.
 * - The event is considered <em>fully consumed</em> if its current number of consumptions \c c equals
 * \c maxConsumptionCount.
 * - The event is considered <em>sufficiently consumed</em> if <tt>minConsumptionCount <= c <= maxConsumptionCount</tt>.
 * - The class is initialized in the <em>fully consumed</em> state, so it can not be consumed right away.
 * - Consuming the event is only possible if it is not <em>fully consumed</em> (<tt>c < maxConsumptionCount</tt>).
 * Consuming the event increments \c c by 1. Trying to consume <em>fully consumed</em> event
 * throws \ref gmx::InternalError.
 * - \ref reset returns object into the initial <em>fully consumed</em> state.
 * This function is intended to manually override the consumption limits.
 * - \ref consume \em consumes the event, without doing anything else.
 * This function is intended to manually override the consumption limits.
 * - \ref markEvent enqueues new event into the provided stream, and sets \c to 0. Marking is only
 * possible if the event is <em>sufficiently consumed</em>, otherwise \ref gmx::InternalError
 * is thrown.
 * - \ref waitForEvent \em consumes the event and blocks the host thread until the event
 * is ready (complete).
 * - \ref enqueueWaitEvent \em consumes the event and blocks the inserts a blocking barrier
 * into the provided stream which blocks the execution of all tasks later submitted to this
 * stream until the event is ready (completes).
 *
 * Default <tt>minConsumptionCount=maxConsumptionCount=1</tt> limits mean that each call to \ref markEvent must be followed
 * by exactly one \ref enqueueWaitEvent or \ref enqueueWaitEvent. This is the recommended pattern
 * for most use cases. By providing other constructor arguments, this requirement can be relaxed
 * as needed.
 */
class GpuEventSynchronizer
{
public:
    //! A constructor
    GpuEventSynchronizer(int minConsumptionCount, int maxConsumptionCount) :
        minConsumptionCount_(minConsumptionCount), maxConsumptionCount_(maxConsumptionCount)
    {
        reset();
    }
    GpuEventSynchronizer() : GpuEventSynchronizer(1, 1) {}
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

    /*! \brief Marks the synchronization point in the \p stream and reset the consumption counter.
     *
     * Should be called before implicitly consuming actions (\ref waitForEvent() or \ref enqueueWaitEvent()) are executed or explicit \ref consume() calls are made.
     *
     * If the event has been marked before and not fully consumed, throws \ref gmx::InternalError.
     */
    inline void markEvent(const DeviceStream& deviceStream)
    {
#if !GMX_GPU_CUDA // For now, we have relaxed conditions for CUDA
        if (consumptionCount_ < minConsumptionCount_)
        {
            GMX_THROW(gmx::InternalError("Trying to mark event before fully consuming it"));
        }
#endif
        event_.mark(deviceStream);
        consumptionCount_ = 0;
    }
    /*! \brief Synchronizes the host thread on the marked event.
     *
     * Consumes the event if able, otherwise throws \ref gmx::InternalError.
     */
    inline void waitForEvent()
    {
        consume();
        event_.wait();
        resetIfFullyConsumed();
    }
    //! Checks the completion of the underlying event and consumes the event if it is ready.
    inline bool isReady()
    {
        bool isReady = event_.isReady();
        if (isReady)
        {
            consume();
            resetIfFullyConsumed();
        }
        return isReady;
    }
    //! Checks whether the event was marked (and was not reset since then).
    inline bool isMarked() const { return event_.isMarked(); }
    /*! \brief Manually consume the event without waiting for it.
     *
     * If the event is already fully consumed, throws \ref gmx::InternalError.
     */
    inline void consume()
    {
#if !GMX_GPU_CUDA // For now, we have relaxed conditions for CUDA
        if (consumptionCount_ >= maxConsumptionCount_)
        {
            GMX_THROW(gmx::InternalError(
                    "Trying to consume an event before marking it or after fully consuming it"));
        }
#endif
        consumptionCount_++;
    }
    //! Helper function to reset the event when it is fully consumed.
    inline void resetIfFullyConsumed()
    {
        if (consumptionCount_ == maxConsumptionCount_)
        {
            event_.reset();
        }
    }
    /*! \brief Enqueues a wait for the recorded event in stream \p deviceStream.
     *
     * Consumes the event if able, otherwise throws \ref gmx::InternalError.
     */
    inline void enqueueWaitEvent(const DeviceStream& deviceStream)
    {
        consume();
        event_.enqueueWait(deviceStream);
        resetIfFullyConsumed();
    }

    //! Resets the event to unmarked state, releasing the underlying event object if needed.
    inline void reset()
    {
        // Set such that we can mark new event without triggering an exception, but can not consume.
        consumptionCount_ = maxConsumptionCount_;
        event_.reset();
    }

    //! Set the event consumption limits.
    inline void setConsumptionLimits(int minConsumptionCount, int maxConsumptionCount)
    {
        minConsumptionCount_ = minConsumptionCount;
        maxConsumptionCount_ = maxConsumptionCount;
    }

private:
    DeviceEvent event_;
    int         consumptionCount_;
#if defined(__clang__) && GMX_GPU_CUDA
    [[maybe_unused]]
#endif
    int minConsumptionCount_; // Unused in CUDA builds, yet
    int maxConsumptionCount_;
};

#endif
