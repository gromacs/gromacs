/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 *  \brief Implements the GPU region timer for SYCL.
 *
 *  \author Andrey Alekseenko <al42and@gmail.com>
 *
 *  \inlibraryapi
 */

#ifndef GMX_GPU_UTILS_GPUREGIONTIMER_SYCL_H
#define GMX_GPU_UTILS_GPUREGIONTIMER_SYCL_H

#include "gromacs/gpu_utils/devicebuffer_sycl.h"
#include "gromacs/gpu_utils/gputraits.h"
#include "gromacs/utility/stringutil.h"

#include "gpuregiontimer.h"

// Disabling Doxygen to avoid it having GpuRegionTimerImpl from both OpenCL and SYCL
#ifndef DOXYGEN

/*! \libinternal \brief
 * The stub of SYCL implementation of the GPU code region timing.
 *
 * Does not do anything.
 *
 * \todo Implement
 */
class GpuRegionTimerImpl
{
public:
    GpuRegionTimerImpl()  = default;
    ~GpuRegionTimerImpl() = default;
    //! No copying
    GpuRegionTimerImpl(const GpuRegionTimerImpl&) = delete;
    //! No assignment
    GpuRegionTimerImpl& operator=(GpuRegionTimerImpl&&) = delete;
    //! Moving is disabled but can be considered in the future if needed
    GpuRegionTimerImpl(GpuRegionTimerImpl&&) = delete;

    /*! \brief Should be called before the region start. */
    inline void openTimingRegion(const DeviceStream& /*unused*/) {}
    /*! \brief Should be called after the region end. */
    inline void closeTimingRegion(const DeviceStream& /*unused*/) {}
    /*! \brief Returns the last measured region timespan (in milliseconds) and calls \c reset(). */
    // NOLINTNEXTLINE readability-convert-member-functions-to-static
    inline double getLastRangeTime() { return 0; }
    /*! \brief Resets the internal state, releasing the used handles, if any. */
    inline void reset() {}
    /*! \brief Returns a new raw timing event
     * for passing into individual GPU API calls
     * within the region if the API requires it (e.g. on OpenCL).
     */
    inline CommandEvent* fetchNextEvent() { return nullptr; }
};

//! Short-hand for external use
using GpuRegionTimer = GpuRegionTimerWrapper<GpuRegionTimerImpl>;

#endif // !DOXYGEN

#endif
