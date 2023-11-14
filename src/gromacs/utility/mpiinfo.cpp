/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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

#include "gmxpre.h"

#include "gromacs/utility/mpiinfo.h"

#include <cstdlib>

#include <regex>
#include <string>
#include <string_view>

// need to include gmxapi.h here as mpi.h needs to be included before mpi-ext.h
#include "gromacs/utility/gmxmpi.h"

#if HAVE_MPI_EXT
#    include <mpi-ext.h>
#endif

namespace gmx
{

namespace
{

//! Return a copy of the version string from the MPI library or an empty string
std::string findMpiLibraryVersionString()
{
    // All MPI libraries should define
    // MPI_MAX_LIBRARY_VERSION_STRING because it is part of the
    // standard, but we may as well check for it. This also makes
    // this code safe for thread-MPI, which does not define this
    // version string.
#ifdef MPI_MAX_LIBRARY_VERSION_STRING
    char versionString[MPI_MAX_LIBRARY_VERSION_STRING];
    int  length = -1;
    // Note this returns a zero-terminated string when called from
    // C/C++.
    MPI_Get_library_version(versionString, &length);
    return versionString;
#else
    return "";
#endif
}

} // namespace

std::string_view mpiLibraryVersionString()
{
    // Avoid calling into the MPI library and then making a new
    // std::string every time this function is called. C++ standard
    // guarantees that this is thread-safe.
    static const std::string cachedVersionString = findMpiLibraryVersionString();

    return cachedVersionString;
}

bool usingIntelMpi()
{
    // Sample output from Intel MPI:
    // 'Intel(R) MPI Library 2021.10 for Linux* OS
    // '
    // ie. it includes a newline!
    return mpiLibraryVersionString().find("Intel(R) MPI Library") != std::string::npos;
}

GpuAwareMpiStatus checkMpiCudaAwareSupport()
{
#if MPI_SUPPORTS_CUDA_AWARE_DETECTION
    // With OMPI version <=4.x, this function doesn't check if UCX PML is built with CUDA-support
    // or if CUDA is disabled at runtime.
    // Expect this function to work only if OMPI uses OB1 PML
    // This is a known issue (https://github.com/open-mpi/ompi/issues/7963) and fix for this is
    // expected soon (written March 2021)
    GpuAwareMpiStatus status = (MPIX_Query_cuda_support() == 1) ? GpuAwareMpiStatus::Supported
                                                                : GpuAwareMpiStatus::NotSupported;
#else
    GpuAwareMpiStatus status = GpuAwareMpiStatus::NotSupported;
#endif

    if (status != GpuAwareMpiStatus::Supported && getenv("GMX_FORCE_GPU_AWARE_MPI") != nullptr)
    {
        status = GpuAwareMpiStatus::Forced;
    }
    return status;
}

GpuAwareMpiStatus checkMpiHipAwareSupport()
{
#if MPI_SUPPORTS_HIP_AWARE_DETECTION
    GpuAwareMpiStatus status = (MPIX_Query_hip_support() == 1) ? GpuAwareMpiStatus::Supported
                                                               : GpuAwareMpiStatus::NotSupported;
#elif MPI_SUPPORTS_ROCM_AWARE_DETECTION
    GpuAwareMpiStatus status = (MPIX_Query_rocm_support() == 1) ? GpuAwareMpiStatus::Supported
                                                                : GpuAwareMpiStatus::NotSupported;
#else
    GpuAwareMpiStatus status = GpuAwareMpiStatus::NotSupported;
#endif

    if (status != GpuAwareMpiStatus::Supported && getenv("GMX_FORCE_GPU_AWARE_MPI") != nullptr)
    {
        status = GpuAwareMpiStatus::Forced;
    }
    return status;
}


GpuAwareMpiStatus checkMpiZEAwareSupport()
{
    GpuAwareMpiStatus status = GpuAwareMpiStatus::NotSupported;
#if MPI_SUPPORTS_ZE_AWARE_DETECTION
    if (MPIX_Query_ze_support() == 1)
    {
        status = GpuAwareMpiStatus::Supported;
    }
#else
    // Find if we are using Intel MPI
    if (usingIntelMpi())
    {
        // If so, then we can decide whether it supports GPU-aware
        // MPI. First, find the library major and minor version
        // numbers.
        //
        // Sample output from Intel MPI:
        // 'Intel(R) MPI Library 2021.10 for Linux* OS
        // '
        // ie. it includes a newline!
        std::string versionString{ mpiLibraryVersionString() };
        std::regex  re("Intel\\(R\\) MPI Library (.*)\\.(.*) for");
        std::smatch match;
        if (std::regex_search(versionString, match, re) && match.size() == 3)
        {
            // These will be zero in case of error, which will work
            // correctly in the following logic.
            int majorVersion = std::atoi(match.str(1).c_str());
            int minorVersion = std::atoi(match.str(2).c_str());

            // Technically 2021.8 had support, but it's so old and
            // slow that GROMACS doesn't consider that to be support.
            if ((majorVersion > 2021) || (majorVersion == 2021 && minorVersion >= 9))
            {
                // Now check whether the user may have asked Intel MPI to do
                // GPU-aware MPI.
                if (const char* environmentValueAsCString = std::getenv("I_MPI_OFFLOAD"))
                {
                    // Now check whether they did ask for GPU-aware MPI.  Note
                    // std::atoi returns zero on error, so that works OK here
                    // because I_MPI_OFFLOAD=0 means no GPU-aware support.
                    if (std::atoi(environmentValueAsCString))
                    {
                        status = GpuAwareMpiStatus::Supported;
                    }
                }
            }
        }
    }
#endif

    if (status != GpuAwareMpiStatus::Supported && getenv("GMX_FORCE_GPU_AWARE_MPI") != nullptr)
    {
        status = GpuAwareMpiStatus::Forced;
    }
    return status;
}

} // namespace gmx
