/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2010- The GROMACS Authors
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
/*! \file
 *  \brief Define functions for detection and initialization for CUDA devices.
 *
 *  \author Szilard Pall <pall.szilard@gmail.com>
 */

#include "gmxpre.h"

#include "gpu_utils.h"

#include "config.h"

#include <cuda_profiler_api.h>

#include <cassert>
#include <cstdio>
#include <cstdlib>

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/hardware/device_information.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/snprintf.h"
#include "gromacs/utility/stringutil.h"

// Set profiler run flag true if nvprof, Nsight Systems, or Nsight Compute tools have been detected.
static const bool cudaProfilerRun =
        ((std::getenv("NVPROF_ID") != nullptr) || (std::getenv("NSYS_PROFILING_SESSION_ID") != nullptr)
         || (std::getenv("NV_NSIGHT_INJECTION_TRANSPORT_TYPE") != nullptr));

bool isHostMemoryPinned(const void* h_ptr)
{
    cudaPointerAttributes memoryAttributes;
    cudaError_t           stat = cudaPointerGetAttributes(&memoryAttributes, h_ptr);

    bool isPinned = false;
    switch (stat)
    {
        case cudaSuccess:
            // In CUDA 11.0, the field called memoryType in
            // cudaPointerAttributes was replaced by a field called
            // type, along with a documented change of behavior when the
            // pointer passed to cudaPointerGetAttributes is to
            // non-registered host memory. That change means that this
            // code needs conditional compilation and different
            // execution paths to function with all supported versions.
#if CUDART_VERSION < 11 * 1000
            isPinned = true;
#else
            isPinned = (memoryAttributes.type == cudaMemoryTypeHost);
#endif
            break;

        case cudaErrorInvalidValue:
            // If the buffer was not pinned, then it will not be recognized by CUDA at all
            isPinned = false;
            // Reset the last error status
            cudaGetLastError();
            break;

        default: CU_RET_ERR(stat, "Unexpected CUDA error");
    }
    return isPinned;
}

void startGpuProfiler()
{
    /* An environment variable set by the CUDA tool in use indicates that
       mdrun is executed in the CUDA profiler.
       With the profiler runtime flag "--profile-from-start off", the profiler will
       be started here. This way we can avoid tracing the CUDA events from the
       first part of the run. Starting the profiler again does nothing.
     */
    if (cudaProfilerRun)
    {
        cudaError_t stat;
        stat = cudaProfilerStart();
        CU_RET_ERR(stat, "cudaProfilerStart failed");
    }
}

void stopGpuProfiler()
{
    /* Stopping the nvidia here allows us to eliminate the subsequent
       API calls from the trace, e.g. uninitialization and cleanup. */
    if (cudaProfilerRun)
    {
        cudaError_t stat;
        stat = cudaProfilerStop();
        CU_RET_ERR(stat, "cudaProfilerStop failed");
    }
}

void resetGpuProfiler()
{
    /* With CUDA <=7.5 the profiler can't be properly reset; we can only start
     *  the profiling here (can't stop it) which will achieve the desired effect if
     *  the run was started with the profiling disabled.
     *
     * TODO: add a stop (or replace it with reset) when this will work correctly in CUDA.
     * stopGpuProfiler();
     */
    if (cudaProfilerRun)
    {
        startGpuProfiler();
    }
}

/*! \brief Check and act on status returned from peer access CUDA call
 *
 * If status is "cudaSuccess", we continue. If
 * "cudaErrorPeerAccessAlreadyEnabled", then peer access has already
 * been enabled so we ignore. If "cudaErrorInvalidDevice" then the
 * run is trying to access an invalid GPU, so we throw an error. If
 * "cudaErrorInvalidValue" then there is a problem with the arguments
 * to the CUDA call, and we throw an error. These cover all expected
 * statuses, but if any other is returned we issue a warning and
 * continue.
 *
 * \param[in] stat           CUDA call return status
 * \param[in] gpuA           ID for GPU initiating peer access call
 * \param[in] gpuB           ID for remote GPU
 * \param[in] mdlog          Logger object
 * \param[in] cudaCallName   name of CUDA peer access call
 */
static void peerAccessCheckStat(const cudaError_t    stat,
                                const int            gpuA,
                                const int            gpuB,
                                const gmx::MDLogger& mdlog,
                                const char*          cudaCallName)
{

    if (stat == cudaErrorPeerAccessAlreadyEnabled)
    {
        // Since peer access has already been enabled, this error can safely be ignored.
        // Now clear the error internally within CUDA:
        cudaGetLastError();
        return;
    }
    if ((stat == cudaErrorInvalidDevice) || (stat == cudaErrorInvalidValue))
    {
        std::string errorString =
                gmx::formatString("%s from GPU %d to GPU %d failed", cudaCallName, gpuA, gpuB);
        CU_RET_ERR(stat, errorString);
    }
    if (stat != cudaSuccess)
    {
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendTextFormatted(
                        "GPU peer access not enabled between GPUs %d and %d due to unexpected "
                        "return value from %s. %s",
                        gpuA,
                        gpuB,
                        cudaCallName,
                        gmx::getDeviceErrorString(stat).c_str());
        // Clear the error internally within CUDA
        cudaGetLastError();
    }
}

void setupGpuDevicePeerAccess(gmx::ArrayRef<const int> gpuIdsToUse, const gmx::MDLogger& mdlog)
{
    cudaError_t stat;

    // take a note of currently-set GPU
    int currentGpu;
    stat = cudaGetDevice(&currentGpu);
    CU_RET_ERR(stat, "cudaGetDevice in setupGpuDevicePeerAccess failed");

    std::string message = gmx::formatString(
            "Note: Peer access enabled between the following GPU pairs in the node:\n ");
    bool peerAccessEnabled = false;

    for (unsigned int i = 0; i < gpuIdsToUse.size(); i++)
    {
        int gpuA = gpuIdsToUse[i];
        stat     = cudaSetDevice(gpuA);
        if (stat != cudaSuccess)
        {
            GMX_LOG(mdlog.warning)
                    .asParagraph()
                    .appendTextFormatted(
                            "GPU peer access not enabled due to unexpected return value from "
                            "cudaSetDevice(%d). %s",
                            gpuA,
                            gmx::getDeviceErrorString(stat).c_str());
            return;
        }
        for (unsigned int j = 0; j < gpuIdsToUse.size(); j++)
        {
            if (j != i)
            {
                int gpuB          = gpuIdsToUse[j];
                int canAccessPeer = 0;
                stat              = cudaDeviceCanAccessPeer(&canAccessPeer, gpuA, gpuB);
                peerAccessCheckStat(stat, gpuA, gpuB, mdlog, "cudaDeviceCanAccessPeer");

                if (canAccessPeer)
                {
                    stat = cudaDeviceEnablePeerAccess(gpuB, 0);
                    peerAccessCheckStat(stat, gpuA, gpuB, mdlog, "cudaDeviceEnablePeerAccess");

                    message           = gmx::formatString("%s%d->%d ", message.c_str(), gpuA, gpuB);
                    peerAccessEnabled = true;
                }
            }
        }
    }

    // re-set GPU to that originally set
    stat = cudaSetDevice(currentGpu);
    if (stat != cudaSuccess)
    {
        CU_RET_ERR(stat, "cudaSetDevice in setupGpuDevicePeerAccess failed");
        return;
    }

    if (peerAccessEnabled)
    {
        GMX_LOG(mdlog.info).asParagraph().appendTextFormatted("%s", message.c_str());
    }
}

void checkPendingDeviceErrorBetweenSteps()
{
    std::string errorPrefix =
            "An unhandled error from a CUDA operation during the current MD step was detected:";
    gmx::checkDeviceError(cudaGetLastError(), errorPrefix);
}
