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

#include "gromacs/gpu_utils/gputraits_ocl.h"
#include "gromacs/gpu_utils/oclutils.h"

#include "gpuregiontimer.h"

/*! \libinternal \brief
 * This is a GPU region timing implementation for OpenCL.
 * It provides methods for measuring the last timespan.
 * Copying/assignment is disabled since the underlying timing events are owned by this.
 */
// cppcheck-suppress noConstructor
class GpuRegionTimerImpl
{
    /*! \brief The underlying individual timing events array.
     * The maximum size is chosen arbitrarily to work with current code, and can be changed.
     * There is simply no need for run-time resizing, and it's unlikely we'll ever need more than 10.
     */
    std::array<cl_event, 10> events_ = {{nullptr}};
    //! Index of the active event
    size_t                   currentEvent_ = 0;

    public:

        GpuRegionTimerImpl()  = default;
        ~GpuRegionTimerImpl() = default;
        //! No copying
        GpuRegionTimerImpl(const GpuRegionTimerImpl &)       = delete;
        //! No assignment
        GpuRegionTimerImpl &operator=(GpuRegionTimerImpl &&) = delete;
        //! Moving is disabled but can be considered in the future if needed
        GpuRegionTimerImpl(GpuRegionTimerImpl &&)            = delete;

        /*! \brief Will be called before the region start. */
        inline void openTimingRegion(CommandStream){}

        /*! \brief Will be called after the region end. */
        inline void closeTimingRegion(CommandStream){}

        /*! \brief Returns the last measured region timespan (in milliseconds) and calls reset() */
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

        /*! \brief Resets internal state */
        inline void reset()
        {
            for (size_t i = 0; i < currentEvent_; i++)
            {
                if (events_[i]) // This conditional is ugly, but is required to make some tests (e.g. empty domain) pass
                {
                    GMX_ASSERT(CL_SUCCESS == clReleaseEvent(events_[i]), "OpenCL event release failure");
                }
            }
            currentEvent_ = 0;
            // As long as we're doing nullptr checks, we might want to be extra cautious.
            events_.fill(nullptr);
        }

        /*! \brief Returns a new raw timing event
         * for passing into individual GPU API calls
         * within the region if the API requires it (e.g. on OpenCL).
         */
        inline CommandEvent *fetchNextEvent()
        {
            GMX_ASSERT(currentEvent_ < events_.size(), "Increase c_maxEventNumber_ if needed");
            cl_event *result = &events_[currentEvent_];
            currentEvent_++;
            return result;
        }
};

//! Short-hand for external use
using GpuRegionTimer = GpuRegionTimerWrapper<GpuRegionTimerImpl>;

#endif
