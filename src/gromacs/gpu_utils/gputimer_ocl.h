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

/*! \libinternal \file
 *  \brief Defines the GPU timer class in CUDA.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 */

#ifndef GMX_GPU_UTILS_GPUTIMER_OCL_H
#define GMX_GPU_UTILS_GPUTIMER_OCL_H

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/gpu_utils/oclutils.h"

/*! \libinternal \brief
 * This is a GPU timing class, based on OpenCL events. FIXME
 * It allows for host-side tracking of the accumulated execution timespans in the CUDA code
 * (measuring kernel or transfers duration).
 * Note however, that the data reported by CUDA events is not correct with multiple CUDA streams (e.g. PME and NB).
 * The reason is that currently the underlying CUDA events use the host-side queue.
 */
class GpuTimer
{
    //! The timer state used for debug-only assertions
    enum class TimerState
    {
        Idle,
        Recording,
        NeedsSynchronization
    } debugState_;

    //! brief Maximum number of cl_events in the buffer
    static constexpr size_t c_maxEventNumber_ = 5;

    //! The underlying OpenCL timing events
    cl_event      events_[c_maxEventNumber_];
    //! Index of the active event
    size_t        currentEvent_;
    //! The number of times the timespan has been measured
    unsigned int  callCount_;
    //! The accumulated duration of the timespans measured (milliseconds)
    double        totalMilliseconds_;

    public:
        GpuTimer();
        ~GpuTimer();

        /*! \brief
         * To be called before the kernel/transfer launch.
         *
         * \param[in] s   The OpenCL command stream where the event being measured takes place.
         */
        void startRecording(cl_command_queue s);
        /*! \brief
         * To be called after the kernel/transfer launch.
         *
         * \param[in] s   The OpenCL command stream where the event being measured takes place.
         */
        void stopRecording(cl_command_queue s);
        /*! \brief
         * Accumulates the last timespan into the the total duration.
         * To be called after stopRecording() and the CUDA stream of the event having been synchronized.
         * \returns The last timespan.
         */
        float getLastTimeMilliseconds();
        /*! \brief Resets the timing event total time/call count to zero values. */
        void reset();
        /*! \brief Gets total time recorded. */
        float getTotalTimeMilliseconds() const;
        /*! \brief Gets total call count recorded. */
        unsigned int getCallCount() const;

        /*! \brief
         * Gets the underlying timing event for passing into OpenCL functions.
         * \returns The pointer to the underlying OpenCL timing event.
         */
        cl_event *handle();
};

#endif
