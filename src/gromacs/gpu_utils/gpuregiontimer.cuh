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
 *  \brief Implements the GPU region timer for CUDA.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 *
 *  \inlibraryapi
 */

#ifndef GMX_GPU_UTILS_GPUREGIONTIMER_CUH
#define GMX_GPU_UTILS_GPUREGIONTIMER_CUH

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/gputraits.cuh"

#include "gpuregiontimer.h"

/*! \libinternal \brief
 * This is a GPU region timing implementation for CUDA.
 * It provides methods for measuring the last timespan.
 * Copying/assignment is disabled since the underlying timing events are owned by this.
 */
class GpuRegionTimerImpl
{
    //! The underlying timing event pair - the beginning and the end of the timespan
    cudaEvent_t eventStart_, eventStop_;

public:
    GpuRegionTimerImpl()
    {
        const int eventFlags = cudaEventDefault;
        CU_RET_ERR(cudaEventCreate(&eventStart_, eventFlags), "GPU timing creation failure");
        CU_RET_ERR(cudaEventCreate(&eventStop_, eventFlags), "GPU timing creation failure");
    }
    ~GpuRegionTimerImpl()
    {
        CU_RET_ERR(cudaEventDestroy(eventStart_), "GPU timing destruction failure");
        CU_RET_ERR(cudaEventDestroy(eventStop_), "GPU timing destruction failure");
    }
    //! No copying
    GpuRegionTimerImpl(const GpuRegionTimerImpl&) = delete;
    //! No assignment
    GpuRegionTimerImpl& operator=(GpuRegionTimerImpl&&) = delete;
    //! Moving is disabled but can be considered in the future if needed
    GpuRegionTimerImpl(GpuRegionTimerImpl&&) = delete;

    /*! \brief Will be called before the region start. */
    inline void openTimingRegion(CommandStream s)
    {
        CU_RET_ERR(cudaEventRecord(eventStart_, s), "GPU timing recording failure");
    }

    /*! \brief Will be called after the region end. */
    inline void closeTimingRegion(CommandStream s)
    {
        CU_RET_ERR(cudaEventRecord(eventStop_, s), "GPU timing recording failure");
    }

    /*! \brief Returns the last measured region timespan (in milliseconds) and calls reset() */
    inline double getLastRangeTime()
    {
        float milliseconds = 0.0;
        CU_RET_ERR(cudaEventElapsedTime(&milliseconds, eventStart_, eventStop_),
                   "GPU timing update failure");
        reset();
        return milliseconds;
    }

    /*! \brief Resets internal state */
    inline void reset() {}

    /*! \brief Returns a new raw timing event
     * for passing into individual GPU API calls.
     * This is just a dummy in CUDA.
     */
    inline CommandEvent* fetchNextEvent() { return nullptr; }
};

//! Short-hand for external use
using GpuRegionTimer = GpuRegionTimerWrapper<GpuRegionTimerImpl>;

#endif
