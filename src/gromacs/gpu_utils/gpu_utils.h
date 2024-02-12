/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
 *  \brief Declare functions for detection and initialization for GPU devices.
 *
 *  \author Szilard Pall <pall.szilard@gmail.com>
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 *  \inlibraryapi
 */
#ifndef GMX_GPU_UTILS_GPU_UTILS_H
#define GMX_GPU_UTILS_GPU_UTILS_H

#include "gromacs/utility/arrayref.h"

namespace gmx
{
class MDLogger;
}

//! Enum which is only used to describe transfer calls at the moment
enum class GpuApiCallBehavior : int
{
    //! Synchronous
    Sync,
    //! Asynchronous
    Async,
    //! Size of the enumeration
    Count
};

//! String corresponding to GPU API call behavior
const char* enumValueToString(GpuApiCallBehavior enumValue);

//! Types of actions associated to waiting or checking the completion of GPU tasks
enum class GpuTaskCompletion
{
    Wait, /*<< Issue a blocking wait for the task */
    Check /*<< Only check whether the task has completed */
};

/*! \brief Starts the GPU profiler if mdrun is being profiled.
 *
 *  When a profiler run is in progress (based on the presence of the NVPROF_ID
 *  env. var.), the profiler is started to begin collecting data during the
 *  rest of the run (or until stopGpuProfiler is called).
 *
 *  Note that this is implemented only for the CUDA API.
 */
void startGpuProfiler();


/*! \brief Resets the GPU profiler if mdrun is being profiled.
 *
 * When a profiler run is in progress (based on the presence of the NVPROF_ID
 * env. var.), the profiler data is restet in order to eliminate the data collected
 * from the preceding part fo the run.
 *
 * This function should typically be called at the mdrun counter reset time.
 *
 * Note that this is implemented only for the CUDA API.
 */
void resetGpuProfiler();


/*! \brief Stops the CUDA profiler if mdrun is being profiled.
 *
 *  This function can be called at cleanup when skipping recording
 *  recording subsequent API calls from being traces/profiled is desired,
 *  e.g. before uninitialization.
 *
 *  Note that this is implemented only for the CUDA API.
 */
void stopGpuProfiler();

//! Tells whether the host buffer was pinned for non-blocking transfers. Only implemented for CUDA.
bool isHostMemoryPinned(const void* h_ptr);

/*! \brief Enable peer access between GPUs where supported
 * \param[in] gpuIdsToUse   List of GPU IDs in use
 * \param[in] mdlog         Logger object
 */
void setupGpuDevicePeerAccess(gmx::ArrayRef<const int> gpuIdsToUse, const gmx::MDLogger& mdlog);

/*! \brief Check the platform-defaults and environment variable to decide whether GPU timings
 * should be enabled.
 *
 * Currently, timings are enabled for OpenCL, but disabled for CUDA and SYCL. This can be overridden
 * by \c GMX_ENABLE_GPU_TIMING and \c GMX_DISABLE_GPU_TIMING environment variables.
 */
bool decideGpuTimingsUsage();

/*! \brief Check for API errors to avoid propagating these across e.g. MD steps.
 */
void checkPendingDeviceErrorBetweenSteps();

#endif
