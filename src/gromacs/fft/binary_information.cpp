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
 * FFT support in the currently running binary
 *
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/fft/binary_information.h"

#include "config.h"

#if GMX_FFT_FFTW3 || GMX_FFT_ARMPL_FFTW3
// Needed for construction of the FFT library description string
#    include <fftw3.h>
#endif

#if GMX_FFT_MKL
#    include <mkl.h>
#endif

#if GMX_GPU_FFT_ONEMATH
#    include <oneapi/math/dft.hpp>
#endif

#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{
namespace
{

// This function is duplicated in linearalgebra, keep both in sync.
std::string describeMkl()
{
#if GMX_FFT_MKL
    MKLVersion mklVersion;
    mkl_get_version(&mklVersion);
    auto description = gmx::formatString("Intel MKL version %d.%d.%d Build %s",
                                         mklVersion.MajorVersion,
                                         mklVersion.MinorVersion,
                                         mklVersion.UpdateVersion,
                                         mklVersion.Build);
    if (mklVersion.ProductStatus != std::string("Product"))
    {
        description += " ";
        description += mklVersion.ProductStatus;
    }
    return description;
#else
    return "Intel MKL";
#endif
}

std::string describeOneMath()
{
#if GMX_GPU_FFT_ONEMATH
    std::string description = "oneMath interface library (backends:";
#    ifdef ONEMATH_ENABLE_CUFFT_BACKEND
    description += " cuFFT";
#    endif
#    ifdef ONEMATH_ENABLE_MKLGPU_BACKEND
    description += " MKLGPU";
#    endif
#    ifdef ONEMATH_ENABLE_ROCFFT_BACKEND
    description += " rocFFT";
#    endif
#    ifdef ONEMATH_ENABLE_PORTFFT_BACKEND
    description += " portFFT";
#    endif
    description += ")";
    return description;
#else
    GMX_RELEASE_ASSERT(false, "describeOneMath called in a build without oneMath");
    return "";
#endif
}

} // namespace

//! Construct a string that describes the library that provides CPU FFT support to this build
std::string cpuFftDescription()
{
// Define the FFT description string
#if GMX_FFT_FFTW3 || GMX_FFT_ARMPL_FFTW3
#    if GMX_NATIVE_WINDOWS
    // Don't buy trouble
    return "fftw3";
#    else
    // Use the version string provided by libfftw3
#        if GMX_DOUBLE
    return fftw_version;
#        else
    return fftwf_version;
#        endif
#    endif
#endif
#if GMX_FFT_MKL
    return describeMkl();
#endif
#if GMX_FFT_FFTPACK
    return "fftpack (built-in)";
#endif
}

//! Construct a string that describes the library that provides GPU FFT support to this build
std::string gpuFftDescription()
{
    if (GMX_GPU)
    {
        if (GMX_GPU_FFT_CUFFT)
        {
            return "cuFFT";
        }
        else if (GMX_GPU_FFT_CLFFT)
        {
            return "clFFT";
        }
        else if (GMX_GPU_FFT_VKFFT)
        {
            return std::string("VkFFT ") + vkfft_VERSION;
        }
        else if (GMX_GPU_FFT_MKL)
        {
            return describeMkl();
        }
        else if (GMX_GPU_FFT_ONEMATH)
        {
            return describeOneMath();
        }
        else if (GMX_GPU_FFT_ROCFFT)
        {
            return std::string("rocFFT ") + rocfft_VERSION;
        }
        else if (GMX_GPU_FFT_HIPFFT)
        {
            return std::string("hipFFT ") + hipfft_VERSION;
        }
        else if (GMX_GPU_FFT_BBFFT)
        {
            return std::string("Double-Batched FFT Library ") + bbfft_VERSION;
        }
        else
        {
            /* Some SYCL builds have no support for GPU FFT,
             * but that's a corner case not intended for general users */
            GMX_RELEASE_ASSERT(GMX_GPU_SYCL,
                               "Only the SYCL build can function without a GPU FFT library");
            return "none / unknown";
        }
    }
    else
    {
        return "none";
    }
}

/*! \brief Construct a string that describes the library (if any)
 * that provides multi-GPU FFT support to this build */
std::string multiGpuFftDescription()
{
    if (GMX_USE_Heffte)
    {
        if (GMX_GPU_FFT_CUFFT)
        {
            // This could be either in a CUDA or SYCL build, but the
            // distinction does not matter here.
            return gmx::formatString("HeFFTe %s with cuFFT backend", Heffte_VERSION);
        }
        else if (GMX_GPU_HIP && GMX_GPU_FFT_HIPFFT)
        {
            return gmx::formatString("HeFFTe %s with hipFFT backend", Heffte_VERSION);
        }
        else if (GMX_GPU_SYCL && GMX_GPU_FFT_MKL)
        {
            return gmx::formatString("HeFFTe %s with oneMKL backend", Heffte_VERSION);
        }
        else if ((GMX_GPU_SYCL || GMX_GPU_HIP) && GMX_GPU_FFT_ROCFFT)
        {
            return gmx::formatString("HeFFTe %s with rocFFT backend", Heffte_VERSION);
        }
        else
        {
            return gmx::formatString("HeFFTe %s with unknown backend", Heffte_VERSION);
        }
    }
    else if (GMX_USE_cuFFTMp)
    {
        return "cuFFTMp";
    }
    else
    {
        return "none";
    }
}

} // namespace gmx
