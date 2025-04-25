/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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
 * \brief Implements functionality for printing information about the
 * GPU support in the currently running binary
 *
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/gpu_utils/binary_information.h"

#include "config.h"

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "buildinfo.h"
#include "gpuinfo.h"

#if GMX_GPU_SYCL
#    include "gromacs/gpu_utils/gmxsycl.h"
#endif

#if GMX_GPU_HIP
#    include "gromacs/gpu_utils/hiputils.h"
#endif

namespace
{

#if GMX_GPU_CUDA
std::string getCudaDriverVersionString()
{
    int cuda_driver = 0;
    if (cudaDriverGetVersion(&cuda_driver) != cudaSuccess)
    {
        return "N/A";
    }
    return gmx::formatString("%d.%d", cuda_driver / 1000, cuda_driver % 100);
}

std::string getCudaRuntimeVersionString()
{
    int cuda_runtime = 0;
    if (cudaRuntimeGetVersion(&cuda_runtime) != cudaSuccess)
    {
        return "N/A";
    }
    return gmx::formatString("%d.%d", cuda_runtime / 1000, cuda_runtime % 100);
}
#endif

#if GMX_GPU_SYCL
std::string getSyclVersion()
{
#    if GMX_SYCL_DPCPP
#        ifdef __LIBSYCL_MAJOR_VERSION
    return gmx::formatString("%d (libsycl %d.%d.%d)",
                             __SYCL_COMPILER_VERSION,
                             __LIBSYCL_MAJOR_VERSION,
                             __LIBSYCL_MINOR_VERSION,
                             __LIBSYCL_PATCH_VERSION);
#        else
    return gmx::formatString("%d", __SYCL_COMPILER_VERSION);
#        endif
#    elif GMX_SYCL_ACPP
    return hipsycl::sycl::detail::version_string();
#    else
    GMX_THROW(gmx::InternalError("Not implemented for this SYCL build"));
#    endif
}

std::vector<std::string> getSyclOptionalFeatures()
{
    std::vector<std::string> optionalFeatures;
#    if GMX_SYCL_ENABLE_EXPERIMENTAL_SUBMIT_API
    optionalFeatures.push_back("experimental_submit_api");
#    endif
#    if GMX_HAVE_GPU_GRAPH_SUPPORT
    optionalFeatures.push_back("graphs");
#    endif
    return optionalFeatures;
}

std::string getSyclCompilerVersion()
{
    std::string                    versionStr       = getSyclVersion();
    const std::vector<std::string> optionalFeatures = getSyclOptionalFeatures();
    if (optionalFeatures.empty())
    {
        return versionStr;
    }
    else
    {
        return versionStr + " with " + gmx::joinStrings(optionalFeatures, ",");
    }
}

#endif // GMX_SYCL_DPCPP || GMX_SYCL_ACPP

#if GMX_GPU_HIP
std::string getHipDriverAndRuntimeVersionString()
{
    int hipDriver = 0;
    if (hipDriverGetVersion(&hipDriver) != hipSuccess)
    {
        std::ignore = hipGetLastError();
        return "N/A";
    }
    return gmx::formatString("%d.%d.%d", hipDriver / 10000000, hipDriver / 100000 % 100, hipDriver % 100000);
}
#endif

const char* getGpuImplementationString()
{
    // Some flavors of clang complain about unreachable returns.
    CLANG_DIAGNOSTIC_IGNORE("-Wunreachable-code-return")
    if (GMX_GPU)
    {
        if (GMX_GPU_CUDA)
        {
            return "CUDA";
        }
        else if (GMX_GPU_OPENCL)
        {
            return "OpenCL";
        }
        else if (GMX_GPU_HIP)
        {
            return "HIP";
        }
        else if (GMX_GPU_SYCL)
        {
            if (GMX_SYCL_DPCPP)
            {
                return "SYCL (oneAPI DPC++)";
            }
            else if (GMX_SYCL_ACPP)
            {
                return "SYCL (AdaptiveCpp)";
            }
            else
            {
                return "SYCL (unknown)";
            }
        }
        else
        {
            GMX_RELEASE_ASSERT(false, "Unknown GPU configuration");
            return "impossible";
        }
    }
    else
    {
        return "disabled";
    }
    CLANG_DIAGNOSTIC_RESET
}

} // namespace

namespace gmx
{

std::unordered_map<std::string, std::string> gpuDescriptions()
{
    std::unordered_map<std::string, std::string> descriptions;
    // Note that these string keys must be kept in sync with
    // those in mdrun/binary_information.cpp
    descriptions["GPU support"] = getGpuImplementationString();
#if GMX_GPU_OPENCL
    descriptions["OpenCL include dir"] = OPENCL_INCLUDE_DIR;
    descriptions["OpenCL library"]     = OPENCL_LIBRARY;
    descriptions["OpenCL version"]     = OPENCL_VERSION_STRING;
#elif GMX_GPU_CUDA
    descriptions["CUDA compiler"] = CUDA_COMPILER_INFO;
    descriptions["CUDA compiler flags"] =
            std::string(CUDA_COMPILER_FLAGS) + " " + CMAKE_BUILD_CONFIGURATION_CXX_FLAGS;
    descriptions["CUDA driver"]  = getCudaDriverVersionString();
    descriptions["CUDA runtime"] = getCudaRuntimeVersionString();
#elif GMX_SYCL_DPCPP
    descriptions["SYCL version"]        = "oneAPI DPC++ " + getSyclCompilerVersion();
    descriptions["SYCL compiler flags"] = SYCL_DPCPP_COMPILER_FLAGS;
    descriptions["SYCL linker flags"]   = SYCL_DPCPP_LINKER_FLAGS;
#elif GMX_SYCL_ACPP
    descriptions["SYCL version"]        = getSyclCompilerVersion();
    descriptions["SYCL compiler"]       = SYCL_ACPP_COMPILER_LAUNCHER;
    descriptions["SYCL compiler flags"] = SYCL_ACPP_COMPILER_FLAGS;
    descriptions["SYCL GPU flags"]      = SYCL_ACPP_DEVICE_COMPILER_FLAGS;
    descriptions["SYCL targets"]        = SYCL_ACPP_TARGETS;
#elif GMX_GPU_HIP
    descriptions["HIP compiler"] = HIP_COMPILER_INFO;
    descriptions["HIP compiler flags"] =
            std::string(HIP_COMPILER_FLAGS) + " " + CMAKE_BUILD_CONFIGURATION_CXX_FLAGS;
    descriptions["HIP driver/runtime"] = getHipDriverAndRuntimeVersionString();
#endif
    return descriptions;
};

} // namespace gmx
