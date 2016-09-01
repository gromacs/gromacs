/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014,2015,2016, by the GROMACS development team, led by
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
/*! \file
 *  \brief Define functions for managing NVIDIA GPU application clocks via NVML.
 *
 *  \author Szilard Pall <pall.szilard@gmail.com>
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 */

#include "gmxpre.h"

#include "nvml.cuh"

#include "config.h"

#include <string>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/stringutil.h"

#if HAVE_NVML
#include <nvml.h>
#endif

#if HAVE_NVML
#  define HAVE_NVML_APPLICATION_CLOCKS (NVML_API_VERSION >= 6)
#else
#  define HAVE_NVML_APPLICATION_CLOCKS 0
#endif

namespace gmx
{

namespace
{

//! Controls whether any error checking takes place (should perhaps mirror the behaviour of CHECK_CUDA_ERRORS).
const bool gmx_unused checkNvmlErrors = true;
//! Whether NVML support was compiled in.
const bool            haveNvmlSupport = HAVE_NVML;
//! Whether NVML application clock support was compiled in.
const bool            haveApplicationClockSupport = HAVE_NVML_APPLICATION_CLOCKS;

#if HAVE_NVML || defined DOXYGEN
/*! Check for NVML error on the return status of a NVML API call.
 *
 * \param[in]  status          Status returned by NVML API
 * \param[in]  message         Descriptive message to issue if status is not NVML_SUCCESS
 * \throws     NvmlException   If status was checked and found to be != NVML_SUCCESS
 *             std::bad_alloc  If memory allocation fails.
 */
void handleNvmlError(nvmlReturn_t status, const char *message)
{
    if (checkNvmlErrors && haveApplicationClockSupport && status != NVML_SUCCESS)
    {
        throw (NvmlException(gmx::formatString("%s: %s", message, nvmlErrorString(status))));
    }
}

/*! \brief Determines the NVML device ID corresponding to the passed \prop.
 *
 * This is done by matching PCI-E information from \prop with
 * the available NVML devices.
 *
 * \param[in] prop           CUDA device information to enrich with NVML device info.
 * \param[in] name           CUDA device name for decorating error messages
 * \returns                  The NVML device ID.
 *
 * \throws    NvmlException  If an NVML API call failed, or a suitable device doesn't have an NVML device ID.
 */
nvmlDevice_t findNvmlDeviceId(const cudaDeviceProp &prop, const char *name)
{
    // We need conditional compilation around the actual API usages
#if HAVE_NVML_APPLICATION_CLOCKS
    unsigned int deviceCount = 0;
    handleNvmlError(nvmlDeviceGetCount(&deviceCount),
                    "nvmlDeviceGetCount failed");

    for (unsigned int i = 0; i < deviceCount; ++i)
    {
        nvmlDevice_t deviceId;
        handleNvmlError(nvmlDeviceGetHandleByIndex(i, &deviceId),
                        "nvmlDeviceGetHandleByIndex failed");

        nvmlPciInfo_t pci_info;
        handleNvmlError(nvmlDeviceGetPciInfo(deviceId, &pci_info),
                        "nvmlDeviceGetPciInfo failed");

        if (static_cast<unsigned int>(prop.pciBusID) == pci_info.bus &&
            static_cast<unsigned int>(prop.pciDeviceID) == pci_info.device &&
            static_cast<unsigned int>(prop.pciDomainID) == pci_info.domain)
        {
            return deviceId;
        }
    }
    throw (NvmlException(gmx::formatString("NVML device ID not found for %s", name)));
#else  /* HAVE_NVML_APPLICATION_CLOCKS */
    GMX_UNUSED_VALUE(prop);
#endif /* HAVE_NVML_APPLICATION_CLOCKS */
}
#endif /* HAVE_NVML */

}      // namespace

void NvmlManager::setup(const cudaDeviceProp &prop, std::string *message)
{
    name_ = prop.name;
    int  cudaVersionNumber    = prop.major * 10 + prop.minor;
    // In practice, CUDA versions always exceed these numbers
    bool clocksMightBeChanged =
        ((0 == gmx_wcmatch("*Tesla*", name_) && cudaVersionNumber >= 35 ) ||
         (0 == gmx_wcmatch("*Quadro*", name_) && cudaVersionNumber >= 52 ));
    message->clear();
    if (!clocksMightBeChanged)
    {
        // This GPU (or possibly CUDA) doesn't support changing
        // clocks, so there's nothing to tell the user.
        return;
    }

    if (!haveNvmlSupport)
    {
        message->assign(gmx::formatString
                            ("NOTE: GROMACS was configured without NVML support hence it can not exploit\n"
                            "      the capabilities of the detected %s GPU to improve performance.\n"
                            "      Recompile with an NVML library compatible with the driver used, or set the GPU clocks manually.",
                            name_));
        return;
    }

#if HAVE_NVML
    if (!haveApplicationClockSupport)
    {
        message->assign(gmx::formatString
                            ("NOTE: GROMACS was compiled with an old NVML library which does not support\n"
                            "      managing the locks of the detected %s GPU to improve performance.\n"
                            "      Consider upgrading NVML (and/or GPU driver) and recompile GROMACS, or set\n"
                            "      the GPU clocks manually.",
                            name_));
        return;
    }

    /* We've compiled with NVML application clocks support, and have a GPU that can use it */
    char        *env;
    //TODO: GMX_GPU_APPLICATION_CLOCKS is currently only used to enable/disable setting of application clocks
    //      this variable can be later used to give a user more fine grained control.
    env = getenv("GMX_GPU_APPLICATION_CLOCKS");
    if (env != NULL && ( strcmp( env, "0") == 0 ||
                         gmx_strcasecmp( env, "OFF") == 0 ||
                         gmx_strcasecmp( env, "DISABLE") == 0 ))
    {
        return;
    }

    // We need conditional compilation around the actual API usages
#if HAVE_NVML_APPLICATION_CLOCKS
    handleNvmlError(nvmlInit(), "nvmlInit failed");

    deviceId_ = findNvmlDeviceId(prop, name_);

    nvmlEnableState_t applicationClocksApiRestrictionStatus;
    handleNvmlError(nvmlDeviceGetAPIRestriction(deviceId_, NVML_RESTRICTED_API_SET_APPLICATION_CLOCKS, &applicationClocksApiRestrictionStatus),
                    "nvmlDeviceGetAPIRestriction failed");

    /* Per NVML docs, applicationClocksApiRestrictionStatus is
       NVML_FEATURE_ENABLED / NVML_FEATURE_DISABLED when the API is
       accessible to only root / all users. (Yes, that looks
       strange, but that's how it is.) */
    clocksCanBeChanged_ = (applicationClocksApiRestrictionStatus == NVML_FEATURE_DISABLED);
#endif /* HAVE_NVML_APPLICATION_CLOCKS */
#endif /* HAVE_NVML */
}

void NvmlManager::changeClocks(std::string *message)
{
    message->clear();
    if (!haveApplicationClockSupport || !clocksCanBeChanged_)
    {
        return;
    }

    // We need conditional compilation around the actual API usages
#if HAVE_NVML_APPLICATION_CLOCKS
    handleNvmlError(nvmlDeviceGetApplicationsClock(deviceId_, NVML_CLOCK_SM, &originalSmClock_),
                    "nvmlDeviceGetApplicationsClock failed");

    handleNvmlError(nvmlDeviceGetApplicationsClock(deviceId_, NVML_CLOCK_MEM, &originalMemoryClock_),
                    "nvmlDeviceGetApplicationsClock failed");

    // Get max application clocks
    unsigned int maxSmClock, maxMemoryClock;
    handleNvmlError(nvmlDeviceGetMaxClockInfo(deviceId_, NVML_CLOCK_SM, &maxSmClock),
                    "nvmlDeviceGetMaxClockInfo failed");

    handleNvmlError(nvmlDeviceGetMaxClockInfo(deviceId_, NVML_CLOCK_MEM, &maxMemoryClock),
                    "nvmlDeviceGetMaxClockInfo failed");

    if (originalSmClock_ >= maxSmClock)
    {
        // There's no reason to change the clock settings, they're already fast

        //TODO: This should probably be integrated into the GPU Properties table.
        message->assign(gmx::formatString
                            ("GPU clocks for %s are (%d,%d)",
                            name_, originalMemoryClock_, originalSmClock_));
        return;
    }

    /* Note: Distinguishing between different types of GPUs here might be necessary in the future,
       e.g. if max application clocks should not be used for certain GPUs. */

    handleNvmlError(nvmlDeviceSetApplicationsClocks(deviceId_, maxMemoryClock, maxSmClock),
                    "nvmlDeviceGetApplicationsClock failed");

    message->assign(gmx::formatString
                        ("Changing GPU clocks for %s to (%d,%d)",
                        name_, maxMemoryClock, maxSmClock));
    clocksWereChanged_ = true;
    ourSmClock_        = maxSmClock;
    ourMemoryClock_    = maxMemoryClock;

    return;
#endif  /* HAVE_NVML_APPLICATION_CLOCKS */
}

void NvmlManager::resetClocks()
{
    if (!haveApplicationClockSupport || !clocksCanBeChanged_ || !clocksWereChanged_)
    {
        return;
    }

    // We need conditional compilation around the actual API usages
#if HAVE_NVML_APPLICATION_CLOCKS
    /* Check if the clocks are still what we set them to.
     * If so, set them back to the state we originally found them in.
     * If not, don't touch them, because something else set them later. */
    unsigned int currentSmClock, currentMemoryClock;
    handleNvmlError(nvmlDeviceGetApplicationsClock(deviceId_, NVML_CLOCK_SM, &currentSmClock),
                    "nvmlDeviceGetApplicationsClock failed");

    handleNvmlError(nvmlDeviceGetApplicationsClock(deviceId_, NVML_CLOCK_MEM, &currentMemoryClock),
                    "nvmlDeviceGetApplicationsClock failed");

    if (currentSmClock == ourSmClock_ &&
        currentMemoryClock == ourMemoryClock_)
    {
        handleNvmlError(nvmlDeviceSetApplicationsClocks(deviceId_, originalMemoryClock_, originalSmClock_),
                        "nvmlDeviceGetApplicationsClock failed");
    }

    handleNvmlError(nvmlShutdown(), "nvmlShutdown failed");
#endif  /* HAVE_NVML_APPLICATION_CLOCKS */
}

} // namespace
