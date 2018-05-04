/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
 *  \brief Implements a GpuHostSynchronizer class for CUDA.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 *  \inlibraryapi
 */
#ifndef GMX_GPU_UTILS_GPUHOSTSYNCHRONIZER_CUH
#define GMX_GPU_UTILS_GPUHOSTSYNCHRONIZER_CUH

#include "gromacs/gpu_utils/gputraits.cuh"
#include "gromacs/utility/gmxassert.h"

/*! \libinternal \brief
 * A class which allows for CPU thread to mark and wait for certain GPU stream execution point.
 */
class GpuHostSynchronizer
{
    public:
        GpuHostSynchronizer()
        {
            cudaError_t stat = cudaEventCreateWithFlags(&event_, cudaEventDisableTiming);
            GMX_RELEASE_ASSERT(stat == cudaSuccess, "cudaEventCreate failed");
        }
        ~GpuHostSynchronizer()
        {
            stat = cudaEventDestroy(&event_);
            GMX_ASSERT(stat == cudaSuccess, "cudaEventDestroy failed");
        }
        //! No copying
        GpuHostSynchronizer(const GpuHostSynchronizer &)       = delete;
        //! No assignment
        GpuHostSynchronizer &operator = (GpuHostSynchronizer &&) = delete;
        //! Moving is disabled but can be considered in the future if needed
        GpuHostSynchronizer(GpuHostSynchronizer &&)            = delete;

        /*! \brief Marks the synchronization point in the \p stream.
         * Should be followed by waitForSyncPoint().
         */
        inline void markSyncPoint(CommandStream stream)
        {
            cudaError_t stat = cudaEventRecord(event_, stream);
            GMX_ASSERT(stat == cudaSuccess, "cudaEventRecord failed");
        }
        /*! \brief Synchronizes the host thread on the marked event. */
        inline void waitForSyncPoint()
        {
            cudaError_t stat = cudaEventSynchronize(event_);
            GMX_ASSERT(stat == cudaSuccess, "cudaEventSynchronize failed");
        }

    private:
        cudaEvent_t event_;
};

#endif
