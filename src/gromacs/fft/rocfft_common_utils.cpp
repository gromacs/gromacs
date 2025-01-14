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

/*! \internal \file
 *  \brief Implements GPU 3D FFT routines for hipSYCL via rocFFT.
 *
 *  \author Andrey Alekseenko <al42and@gmail.com>
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * For hipSYCL, in order to call FFT APIs from the respective vendors
 * using the same DeviceStream as other operations, a vendor extension
 * called "custom operations" is used (see hipSYCL
 * doc/enqueue-custom-operation.md). That effectively enqueues an
 * asynchronous host-side lambda into the same queue. The body of the
 * lambda unpacks the runtime data structures to get the native
 * handles and calls the native FFT APIs.
 *
 * hipSYCL queues operate at a higher level of abstraction than hip
 * streams, with the runtime distributing work to the latter to
 * balance load. It is possible to set the HIP stream in
 * rocfft_execution_info, but then there is no guarantee that a
 * subsequent queue item will run using the same stream. So we
 * currently do not attempt to set the stream.
 *
 *  \ingroup module_fft
 */

#include "gmxpre.h"

#include "rocfft_common_utils.h"

#include "config.h"

#include <vector>

#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

namespace
{

//! Strings that match enum rocfft_status_e in rocfft.h
constexpr std::array<const char*, rocfft_status_invalid_work_buffer + 1> c_rocfftErrorStrings = {
    "success",
    "failure",
    "invalid argument value",
    "invalid dimensions",
    "invalid array type",
    "invalid strides",
    "invalid distance",
    "invalid offset",
    "invalid work buffer"
};

#if (!(GMX_GPU_SYCL && GMX_GPU_FFT_ROCFFT)) && (!(GMX_GPU_HIP && GMX_GPU_FFT_ROCFFT))
#    error Including rocfft common utils in unsupported build config
#endif

} // namespace

//! Helper for consistent error handling
void handleRocFftError(rocfft_status result, const std::string& msg)
{
    if (result != rocfft_status_success)
    {
        if (result <= rocfft_status_invalid_work_buffer)
        {
            GMX_THROW(gmx::InternalError(gmx::formatString(
                    "%s: (error code %d - %s)\n", msg.c_str(), result, c_rocfftErrorStrings[result])));
        }
        else
        {
            GMX_THROW(gmx::InternalError(gmx::formatString("%s: (error code %d)\n", msg.c_str(), result)));
        }
    }
}

//! Helper for consistent error handling
void handleRocFftError(rocfft_status result, const std::string& direction, const std::string& msg)
{
    if (result != rocfft_status_success)
    {
        handleRocFftError(result, msg + " doing " + direction);
    }
}

RocfftInitializer::RocfftInitializer()
{
    rocfft_status result;
    result = rocfft_setup();
    handleRocFftError(result, "rocfft_setup failure");
}

RocfftInitializer::~RocfftInitializer()
{
    // No need to handle any errors in a destructor, and
    // anyway one cannot throw.
    rocfft_cleanup();
}

//! Destructor
RocfftPlan::~RocfftPlan()
{
    // No need to handle any errors in a destructor,
    // and anyway one cannot throw.
    if (plan)
    {
        rocfft_plan_destroy(plan);
        plan = nullptr;
    }
    if (info)
    {
        rocfft_execution_info_destroy(info);
        info = nullptr;
    }
    if (workBuffer)
    {
        GMX_UNUSED_VALUE(hipFree(workBuffer));
        workBuffer = nullptr;
    }
}

} // namespace gmx
