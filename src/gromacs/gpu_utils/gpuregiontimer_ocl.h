/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018, by the GROMACS development team, led by
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
 *  \brief Implements the GPU region timer for OpenCL.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 *
 *  \inlibraryapi
 */

#ifndef GMX_GPU_UTILS_GPUREGIONTIMER_OCL_H
#define GMX_GPU_UTILS_GPUREGIONTIMER_OCL_H

#include <array>

#include "gromacs/gpu_utils/oclutils.h"

#include "gpuregiontimer.h"

template <> struct GpuTraits<GpuFramework::OpenCL>
{
    using CommandStream = cl_command_queue;
    using CommandEvent  = cl_event;
};

//! Short-hand for external use
using GpuRegionTimer = GpuRegionTimerWrapper<GpuFramework::OpenCL>;

// cppcheck-suppress noConstructor
/*! \libinternal \brief
 * The OpenCL implementation of the GPU code region timing.
 * With OpenCL, one has to use cl_event handle for each API call that has to be timed, and
 * accumulate the timing afterwards. As we would like to avoid overhead on API calls,
 * we only query and accumulate cl_event timing at the end of time steps, not after the API calls.
 * Thus, this implementation does not reuse a single cl_event for multiple calls, but instead
 * maintains an array of cl_events to be used within any single code region.
 * The array size is fixed at a small but sufficiently large value for the number of cl_events
 * that might contribute to a timer region, currently 10.
 */
template <> class GpuRegionTimerImpl<GpuFramework::OpenCL>
{
    //! Short-hands
    using CommandStream = typename GpuTraits<GpuFramework::OpenCL>::CommandStream;
    using CommandEvent  = typename GpuTraits<GpuFramework::OpenCL>::CommandEvent;

    /*! \brief The underlying individual timing events array.
     * The maximum size is chosen arbitrarily to work with current code, and can be changed.
     * There is simply no need for run-time resizing, and it's unlikely we'll ever need more than 10.
     */
    std::array<cl_event, 10> events_ = {{nullptr}};
    //! Index of the active event
    size_t                   currentEvent_ = 0;

    public:

        /*! \brief Should be called before the region start. */
        inline void openTimingRegion(CommandStream){}
        /*! \brief Should be called after the region end. */
        inline void closeTimingRegion(CommandStream){}
        /*! \brief Returns the last measured region timespan (in milliseconds) and calls reset(). */
        inline double getLastRangeTime()
        {
            double milliseconds = 0.0;
            for (size_t i = 0; i < currentEvent_; i++)
            {
                if (events_[i]) // This conditional is ugly, but is required to make some tests (e.g. empty domain) pass
                {
                    cl_ulong          start_ns, end_ns;
                    cl_int gmx_unused cl_error;

                    cl_error = clGetEventProfilingInfo(events_[i], CL_PROFILING_COMMAND_START,
                                                       sizeof(cl_ulong), &start_ns, nullptr);
                    GMX_ASSERT(CL_SUCCESS == cl_error, "GPU timing update failure");
                    cl_error = clGetEventProfilingInfo(events_[i], CL_PROFILING_COMMAND_END,
                                                       sizeof(cl_ulong), &end_ns, nullptr);
                    GMX_ASSERT(CL_SUCCESS == cl_error, "GPU timing update failure");
                    milliseconds += (end_ns - start_ns) / 1000000.0;
                }
            }
            reset();
            return milliseconds;
        }
        /*! \brief Resets the internal state, releasing the used cl_events. */
        inline void reset()
        {
            for (size_t i = 0; i < currentEvent_; i++)
            {
                if (events_[i]) // This conditional is ugly, but is required to make some tests (e.g. empty domain) pass
                {
                    cl_int gmx_unused cl_error = clReleaseEvent(events_[i]);
                    GMX_ASSERT(CL_SUCCESS == cl_error, "OpenCL event release failure");
                }
            }
            currentEvent_ = 0;
            // As long as we're doing nullptr checks, we might want to be extra cautious.
            events_.fill(nullptr);
        }
        /*! \brief Provides next unused cl_event for OpenCL API consumption. */
        inline CommandEvent *fetchNextEvent()
        {
            GMX_ASSERT(currentEvent_ < events_.size(), "Increase c_maxEventNumber_ if needed");
            cl_event *result = &events_[currentEvent_];
            currentEvent_++;
            return result;
        }
};

#endif
