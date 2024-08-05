/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 *  \brief Define functions for detection and initialization for HIP devices.
 *
 *  \author Paul Bauer <paul.bauer.q@gmail.com>
 *  \author Julio Maia <julio.maia@amd.com>
 */

#include "gmxpre.h"

#include <cstdlib>

#include <hip/hip_profile.h>

#include "gromacs/gpu_utils/hiputils.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/stringutil.h"

#include "gpu_utils.h"

bool isHostMemoryPinned(const void* h_ptr)
{
    hipPointerAttribute_t memoryAttributes;
    hipError_t            stat = hipPointerGetAttributes(&memoryAttributes, h_ptr);

    bool isPinned = false;
    switch (stat)
    {
        case hipSuccess:
#if HIP_VERSION_MAJOR > 5
            isPinned = (memoryAttributes.type == hipMemoryTypeHost);
#else
            isPinned = (memoryAttributes.memoryType == hipMemoryTypeHost);
#endif
            break;

        case hipErrorInvalidValue:
            // If the buffer was not allocated by ROCm, then it will not be
            // recognized at all, and an error flag set.
            isPinned = false;
            // Reset error
            std::ignore = hipGetLastError();
            break;

        default: gmx::checkDeviceError(stat, "Unexpected HIP error");
    }
    return isPinned;
}

// Profiling is handled in timing instead of here
void startGpuProfiler() {}

void stopGpuProfiler() {}

void resetGpuProfiler() {}

/*! \brief Check and act on status returned from peer access HIP call
 *
 * If status is "hipSuccess", we continue. If
 * "hipErrorPeerAccessAlreadyEnabled", then peer access has already
 * been enabled so we ignore. If "hipErrorInvalidDevice" then the
 * run is trying to access an invalid GPU, so we throw an error. If
 * "hipErrorInvalidValue" then there is a problem with the arguments
 * to the HIP call, and we throw an error. These cover all expected
 * statuses, but if any other is returned we issue a warning and
 * continue.
 *
 * \param[in] stat           HIP call return status
 * \param[in] gpuA           ID for GPU initiating peer access call
 * \param[in] gpuB           ID for remote GPU
 * \param[in] mdlog          Logger object
 * \param[in] hipCallName    name of HIP peer access call
 */
static void peerAccessCheckStat(const hipError_t     stat,
                                const int            gpuA,
                                const int            gpuB,
                                const gmx::MDLogger& mdlog,
                                const char*          hipCallName)
{

    if (stat == hipErrorPeerAccessAlreadyEnabled)
    {
        // Since peer access has already been enabled, this error can safely be ignored.
        // Now clear the error internally within HIP:
        std::ignore = hipGetLastError();
        return;
    }
    if ((stat == hipErrorInvalidDevice) || (stat == hipErrorInvalidValue))
    {
        std::string errorString =
                gmx::formatString("%s from GPU %d to GPU %d failed", hipCallName, gpuA, gpuB);
        gmx::checkDeviceError(stat, errorString);
    }
    if (stat != hipSuccess)
    {
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendTextFormatted(
                        "GPU peer access not enabled between GPUs %d and %d due to unexpected "
                        "return value from %s. %s",
                        gpuA,
                        gpuB,
                        hipCallName,
                        gmx::getDeviceErrorString(stat).c_str());
        // Clear the error internally within HIP
        std::ignore = hipGetLastError();
    }
}

void setupGpuDevicePeerAccess(gmx::ArrayRef<const int> gpuIdsToUse, const gmx::MDLogger& mdlog)
{
    hipError_t stat;

    // take a note of currently-set GPU
    int currentGpu;
    stat = hipGetDevice(&currentGpu);
    gmx::checkDeviceError(stat, "hipGetDevice in setupGpuDevicePeerAccess failed");

    std::string message = gmx::formatString(
            "Note: Peer access enabled between the following GPU pairs in the node:\n ");
    bool peerAccessEnabled = false;

    for (unsigned int i = 0; i < gpuIdsToUse.size(); i++)
    {
        int gpuA = gpuIdsToUse[i];
        stat     = hipSetDevice(gpuA);
        if (stat != hipSuccess)
        {
            GMX_LOG(mdlog.warning)
                    .asParagraph()
                    .appendTextFormatted(
                            "GPU peer access not enabled due to unexpected return value from "
                            "hipSetDevice(%d). %s",
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
                stat              = hipDeviceCanAccessPeer(&canAccessPeer, gpuA, gpuB);
                peerAccessCheckStat(stat, gpuA, gpuB, mdlog, "hipDeviceCanAccessPeer");

                if (canAccessPeer)
                {
                    stat = hipDeviceEnablePeerAccess(gpuB, 0);
                    peerAccessCheckStat(stat, gpuA, gpuB, mdlog, "hipDeviceEnablePeerAccess");

                    message           = gmx::formatString("%s%d->%d ", message.c_str(), gpuA, gpuB);
                    peerAccessEnabled = true;
                }
            }
        }
    }

    // re-set GPU to that originally set
    stat = hipSetDevice(currentGpu);
    gmx::checkDeviceError(stat, "hipSetDevice in setupGpuDevicePeerAccess failed");

    if (peerAccessEnabled)
    {
        GMX_LOG(mdlog.info).asParagraph().appendText(message);
    }
}

void checkPendingDeviceErrorBetweenSteps()
{
    std::string errorPrefix =
            "An unhandled error from a HIP operation during the current MD step was detected:";
    gmx::checkDeviceError(hipGetLastError(), errorPrefix);
}
