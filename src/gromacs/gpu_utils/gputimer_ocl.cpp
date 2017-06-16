/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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

/*! \internal \file
 *  \brief Implements the GpuTimer class with a pair of CUDA events.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 */

#include "gmxpre.h"

#include "gputimer_ocl.h"

#include "gromacs/utility/gmxassert.h"

#if defined(NDEBUG)
static const bool c_debugTimerState = false;
#else
static const bool c_debugTimerState = true;
#endif

//FIXME copy
/*! \brief Returns the duration in milliseconds for the command associated with the event.
 *
 * It then releases the event and sets it to 0.
 * Before calling this function, make sure the command has finished either by
 * calling clFinish or clWaitForEvents.
 * The function returns 0.0 if the input event, *ocl_event, is 0.
 * Don't use this function when more than one wait will be issued for the event.
 */

GpuTimer::GpuTimer()
{
    reset();
}

GpuTimer::~GpuTimer()
{
}

void GpuTimer::startRecording(cl_command_queue)
{
    if (c_debugTimerState)
    {
        GMX_ASSERT(debugState_ == TimerState::Idle, "GPU timer should have been idle");
        debugState_ = TimerState::Recording;
    }
}

void GpuTimer::stopRecording(cl_command_queue)
{
    if (c_debugTimerState)
    {
        GMX_ASSERT(debugState_ == TimerState::Recording, "GPU timer should have been recording");
        debugState_ = TimerState::NeedsSynchronization;
    }
    callCount_++;
}

void GpuTimer::reset()
{
    totalMilliseconds_ = 0.0;
    callCount_         = 0;
    currentEvent_      = 0;
    if (c_debugTimerState)
    {
        debugState_ = TimerState::Idle;
    }
}

float GpuTimer::getLastTimeMilliseconds()
{
    if (c_debugTimerState)
    {
        GMX_ASSERT(debugState_ == TimerState::NeedsSynchronization, "GPU timer should have just stopped recording");
        debugState_ = TimerState::Idle;
    }
    float milliseconds = 0.0;
    if (callCount_ > 0) /* Only the touched events needed */
    {
        for (size_t i = 0; i < currentEvent_; i++)
        {
            cl_ulong start_ns, end_ns;
            cl_int   cl_error = clGetEventProfilingInfo(events_[i], CL_PROFILING_COMMAND_START,
                                                        sizeof(cl_ulong), &start_ns, nullptr);
            GMX_RELEASE_ASSERT(CL_SUCCESS == cl_error, "GPU timing update failure");
            cl_error = clGetEventProfilingInfo(events_[i], CL_PROFILING_COMMAND_END,
                                               sizeof(cl_ulong), &end_ns, nullptr);
            GMX_RELEASE_ASSERT(CL_SUCCESS == cl_error, "GPU timing update failure");
            clReleaseEvent(events_[i]);
            milliseconds += (end_ns - start_ns) / 1000000.0;
        }
        totalMilliseconds_ += milliseconds;
        currentEvent_       = 0;
    }
    return milliseconds;
}

float GpuTimer::getTotalTimeMilliseconds() const
{
    return totalMilliseconds_;
}

unsigned int GpuTimer::getCallCount() const
{
    return callCount_;
}

cl_event *GpuTimer::handle()
{
    GMX_ASSERT(currentEvent_ < c_maxEventNumber_, "Increase c_maxEventNumber_ if needed");
    auto result = &events_[currentEvent_];
    currentEvent_++;
    return result;
}
