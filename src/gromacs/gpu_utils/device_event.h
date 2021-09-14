/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020,2021, by the GROMACS development team, led by
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
 *  \brief Declares DeviceEvent for all build configuraitons
 *
 *  This header may be included from any build configuration and
 *  defers valid GPU declarations to headers valid only in such
 *  build configurations.
 *
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 */
#ifndef GMX_GPU_UTILS_DEVICE_EVENT_H
#define GMX_GPU_UTILS_DEVICE_EVENT_H

#include "config.h"

#include "gromacs/utility/exceptions.h"

#if !GMX_GPU || defined(DOXYGEN)

class DeviceStream;

// [[noreturn]] attributes must be added to the methods, so it's
// easier to silence the warning here and avoid them appearing in
// the Doxygen
#    ifdef __clang__
#        pragma clang diagnostic push
#        pragma clang diagnostic ignored "-Wmissing-noreturn"
#    endif

class DeviceEvent
{
public:
    DeviceEvent() = default;
    // Disable copy, move, and assignment. Move can be allowed, but not needed yet.
    DeviceEvent& operator=(const DeviceEvent&) = delete;
    DeviceEvent(const DeviceEvent&)            = delete;
    DeviceEvent& operator=(DeviceEvent&&) = delete;
    DeviceEvent(DeviceEvent&&)            = delete;

    /*! \brief Marks the synchronization point in the \p stream.
     * Should be followed by waitForEvent().
     */
    inline void mark(const DeviceStream& /*deviceStream*/) // NOLINT readability-convert-member-functions-to-static
    {
        GMX_THROW(gmx::NotImplementedError("Not implemented for non-GPU build"));
    }
    //! Synchronizes the host thread on the marked event.
    inline void wait() // NOLINT readability-convert-member-functions-to-static
    {
        GMX_THROW(gmx::NotImplementedError("Not implemented for non-GPU build"));
    }
    //! Checks the completion of the underlying event.
    inline bool isReady() // NOLINT readability-convert-member-functions-to-static
    {
        GMX_THROW(gmx::NotImplementedError("Not implemented for non-GPU build"));
    }
    //! Enqueues a wait for the recorded event in stream \p stream
    // NOLINTNEXTLINE readability-convert-member-functions-to-static
    inline void enqueueWait(const DeviceStream& /*deviceStream*/)
    {
        GMX_THROW(gmx::NotImplementedError("Not implemented for non-GPU build"));
    }
    //! Checks whether the underlying event was marked.
    inline bool isMarked() const // NOLINT readability-convert-member-functions-to-static
    {
        GMX_THROW(gmx::NotImplementedError("Not implemented for non-GPU build"));
    }

    //! Reset the event (not needed in CUDA)
    // NOLINTNEXTLINE readability-convert-member-functions-to-static
    inline void reset() // NOLINT readability-convert-member-functions-to-static
    {
        GMX_THROW(gmx::NotImplementedError("Not implemented for non-GPU build"));
    }
};

#    ifdef __clang__
#        pragma clang diagnostic pop
#    endif

#elif GMX_GPU_CUDA
#    include "device_event.cuh"
#elif GMX_GPU_OPENCL
#    include "device_event_ocl.h"
#elif GMX_GPU_SYCL
#    include "device_event_sycl.h"
#endif

#endif // GMX_GPU_UTILS_DEVICE_EVENT_H
