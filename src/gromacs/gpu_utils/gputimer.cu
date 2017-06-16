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

#include "gputimer.cuh"

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/utility/gmxassert.h"

static const bool c_useCudaEventBlockingSync = false; /* makes the CPU thread block */

#if defined(NDEBUG)
static const bool c_debugTimerState = false;
#else
static const bool c_debugTimerState = true;
#endif

GpuTimer::GpuTimer()
{
    int eventFlags = (c_useCudaEventBlockingSync ? cudaEventBlockingSync : cudaEventDefault);
    CU_RET_ERR(cudaEventCreate(&eventStart_, eventFlags), "GPU timing creation failure");
    CU_RET_ERR(cudaEventCreate(&eventStop_, eventFlags), "GPU timing creation failure");
    reset();
}

GpuTimer::~GpuTimer()
{
    CU_RET_ERR(cudaEventDestroy(eventStart_), "GPU timing destruction failure");
    CU_RET_ERR(cudaEventDestroy(eventStop_), "GPU timing destruction failure");
}

void GpuTimer::startRecording(cudaStream_t s)
{
    if (c_debugTimerState)
    {
        GMX_ASSERT(debugState_ == TimerState::Idle, "GPU timer should have been idle");
        debugState_ = TimerState::Recording;
    }
    CU_RET_ERR(cudaEventRecord(eventStart_, s), "GPU timing recording failure");
}

void GpuTimer::stopRecording(cudaStream_t s)
{
    if (c_debugTimerState)
    {
        GMX_ASSERT(debugState_ == TimerState::Recording, "GPU timer should have been recording");
        debugState_ = TimerState::NeedsSynchronization;
    }
    CU_RET_ERR(cudaEventRecord(eventStop_, s), "GPU timing recording failure");
    callCount_++;
}

void GpuTimer::reset()
{
    totalMilliseconds_ = 0.0;
    callCount_         = 0;
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
        CU_RET_ERR(cudaEventElapsedTime(&milliseconds, eventStart_, eventStop_), "GPU timing update failure");
        totalMilliseconds_ += milliseconds;
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
