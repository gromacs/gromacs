/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018,2019, by the GROMACS development team, led by
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
 *  \brief Defines the GPU region timer implementation/wrapper classes.
 *  The implementations live in gpuregiontimer.cuh for CUDA and gpuregiontimer_ocl.h for OpenCL.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 */

#ifndef GMX_GPU_UTILS_GPUREGIONTIMER_H
#define GMX_GPU_UTILS_GPUREGIONTIMER_H

#include <string>

#include "gromacs/utility/gmxassert.h"

//! Debug GPU timers in debug builds only
#if defined(NDEBUG)
static const bool c_debugTimerState = false;
#else
static const bool c_debugTimerState = true;
#endif

/*! \libinternal \brief
 * This is a GPU region timing wrapper class.
 * It allows for host-side tracking of the accumulated execution timespans in GPU code
 * (measuring kernel or transfers duration).
 * It also partially tracks the correctness of the timer state transitions,
 * as far as current implementation allows (see TODO in getLastRangeTime() for a disabled check).
 * Internally it uses GpuRegionTimerImpl for measuring regions.
 */
template<typename GpuRegionTimerImpl>
class GpuRegionTimerWrapper
{
    //! The timer state used for debug-only assertions
    enum class TimerState
    {
        Idle,
        Recording,
        Stopped
    } debugState_ = TimerState::Idle;

    //! The number of times the timespan has been measured
    unsigned int callCount_ = 0;
    //! The accumulated duration of the timespans measured (milliseconds)
    double totalMilliseconds_ = 0.0;
    //! The underlying region timer implementation
    GpuRegionTimerImpl impl_;

public:
    /*! \brief
     * To be called before the region start.
     *
     * \param[in] s   The GPU command stream where the event being measured takes place.
     */
    void openTimingRegion(CommandStream s)
    {
        if (c_debugTimerState)
        {
            std::string error = "GPU timer should be idle, but is "
                                + std::string((debugState_ == TimerState::Stopped) ? "stopped" : "recording")
                                + ".";
            GMX_ASSERT(debugState_ == TimerState::Idle, error.c_str());
            debugState_ = TimerState::Recording;
        }
        impl_.openTimingRegion(s);
    }
    /*! \brief
     * To be called after the region end.
     *
     * \param[in] s   The GPU command stream where the event being measured takes place.
     */
    void closeTimingRegion(CommandStream s)
    {
        if (c_debugTimerState)
        {
            std::string error = "GPU timer should be recording, but is "
                                + std::string((debugState_ == TimerState::Idle) ? "idle" : "stopped")
                                + ".";
            GMX_ASSERT(debugState_ == TimerState::Recording, error.c_str());
            debugState_ = TimerState::Stopped;
        }
        callCount_++;
        impl_.closeTimingRegion(s);
    }
    /*! \brief
     * Accumulates the last timespan of all the events used into the total duration,
     * and resets the internal timer state.
     * To be called after closeTimingRegion() and the command stream of the event having been
     * synchronized. \returns The last timespan (in milliseconds).
     */
    double getLastRangeTime()
    {
        if (c_debugTimerState)
        {
            /* The assertion below is commented because it is currently triggered on a special case:
             * the early return before the local non-bonded kernel launch if there is nothing to do.
             * This can be reproduced in OpenCL build by running
             * mdrun-mpi-test -ntmpi 2 --gtest_filter=*Empty*
             * Similarly, the GpuRegionTimerImpl suffers from non-nullptr
             * cl_event conditionals which ideally should only be asserts.
             * TODO: improve internal task scheduling, re-enable the assert, turn conditionals into asserts
             */
            /*
               std::string error = "GPU timer should be stopped, but is " + std::string((debugState_ == TimerState::Idle) ? "idle" : "recording") + ".";
               GMX_ASSERT(debugState_ == TimerState::Stopped, error.c_str());
             */
            debugState_ = TimerState::Idle;
        }
        double milliseconds = impl_.getLastRangeTime();
        totalMilliseconds_ += milliseconds;
        return milliseconds;
    }
    /*! \brief Resets the implementation and total time/call count to zeroes. */
    void reset()
    {
        if (c_debugTimerState)
        {
            debugState_ = TimerState::Idle;
        }
        totalMilliseconds_ = 0.0;
        callCount_         = 0;
        impl_.reset();
    }
    /*! \brief Gets total time recorded (in milliseconds). */
    double getTotalTime() const { return totalMilliseconds_; }
    /*! \brief Gets total call count recorded. */
    unsigned int getCallCount() const { return callCount_; }
    /*! \brief
     * Gets a pointer to a new timing event for passing into individual GPU API calls
     * within the region if they require it (e.g. on OpenCL).
     * \returns The pointer to the underlying single command timing event.
     */
    CommandEvent* fetchNextEvent()
    {
        if (c_debugTimerState)
        {
            std::string error = "GPU timer should be recording, but is "
                                + std::string((debugState_ == TimerState::Idle) ? "idle" : "stopped")
                                + ".";
            GMX_ASSERT(debugState_ == TimerState::Recording, error.c_str());
        }
        return impl_.fetchNextEvent();
    }
};

#endif
