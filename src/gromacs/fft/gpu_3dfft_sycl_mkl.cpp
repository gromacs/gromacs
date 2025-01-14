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

/*! \internal \file
 *  \brief Implements GPU 3D FFT routines for SYCL.
 *
 *  \author Andrey Alekseenko <al42and@gmail.com>
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 *  \ingroup module_fft
 *
 *  In DPC++, we use Intel oneMKL to perform the FFT.
 */

#include "gmxpre.h"

#include "gpu_3dfft_sycl_mkl.h"

#include "config.h"

#include <vector>

#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/devicebuffer_sycl.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

class DeviceContext;

#if (!GMX_SYCL_DPCPP)
#    error This file can only be compiled with Intel DPC++ compiler
#endif

#if (!GMX_FFT_MKL && !GMX_GPU_FFT_ONEMATH)
#    error Must use MKL library CPU FFT when using MKL for GPU FFT and Intel DPC++ compiler
#endif

#include <cstddef>
#pragma clang diagnostic ignored "-Wsuggest-override" // can be removed when support for 2022.0 is dropped
#pragma clang diagnostic ignored "-Wundefined-func-template"

#if GMX_GPU_FFT_MKL
#    include <mkl_version.h>
#    if INTEL_MKL_VERSION < 20250000
#        include <oneapi/mkl/dfti.hpp>
#        define PLACEMENT_INPLACE DFTI_INPLACE
#        define PLACEMENT_NOT_INPLACE DFTI_NOT_INPLACE
#        define COMPLEX_COMPLEX_STORAGE DFTI_COMPLEX_COMPLEX
#    else
#        include <oneapi/mkl/dft.hpp>
#        define PLACEMENT_INPLACE mkl_dft::config_value::INPLACE
#        define PLACEMENT_NOT_INPLACE mkl_dft::config_value::NOT_INPLACE
#        define COMPLEX_COMPLEX_STORAGE mkl_dft::config_value::COMPLEX_COMPLEX
#    endif
#    include <oneapi/mkl/exceptions.hpp>
using mkl_exception = oneapi::mkl::exception;
#elif GMX_GPU_FFT_ONEMATH
// Using oneMath interface library. (GMX_GPU_FFT_ONEMATH)
#    include <oneapi/math/dft.hpp>
#    include <oneapi/math/exceptions.hpp>
#    define PLACEMENT_INPLACE mkl_dft::config_value::INPLACE
#    define PLACEMENT_NOT_INPLACE mkl_dft::config_value::NOT_INPLACE
#    define COMPLEX_COMPLEX_STORAGE mkl_dft::config_value::COMPLEX_COMPLEX
using mkl_exception = oneapi::math::exception;
#endif // GMX_GPU_FFT_MKL


namespace gmx
{

Gpu3dFft::ImplSyclMkl::Descriptor Gpu3dFft::ImplSyclMkl::initDescriptor(const ivec realGridSize)
{
    try
    {
        const std::vector<std::int64_t> realGridDimensions{ realGridSize[XX],
                                                            realGridSize[YY],
                                                            realGridSize[ZZ] };
        return { realGridDimensions };
    }
    catch (mkl_exception& exc)
    {
        GMX_THROW(InternalError(formatString("MKL failure while constructing descriptor: %s", exc.what())));
    }
}

Gpu3dFft::ImplSyclMkl::ImplSyclMkl(bool allocateRealGrid,
                                   MPI_Comm /*comm*/,
                                   ArrayRef<const int> gridSizesInXForEachRank,
                                   ArrayRef<const int> gridSizesInYForEachRank,
                                   int /*nz*/,
                                   const bool           performOutOfPlaceFFT,
                                   const DeviceContext& context,
                                   const DeviceStream&  pmeStream,
                                   ivec                 realGridSize,
                                   ivec                 realGridSizePadded,
                                   ivec                 complexGridSizePadded,
                                   DeviceBuffer<float>* realGrid,
                                   DeviceBuffer<float>* complexGrid) :
    Gpu3dFft::Impl::Impl(performOutOfPlaceFFT),
    realGrid_(*realGrid->buffer_),
    queue_(pmeStream.stream()),
    r2cDescriptor_(initDescriptor(realGridSize)),
    c2rDescriptor_(initDescriptor(realGridSize))
{
    GMX_RELEASE_ASSERT(!allocateRealGrid, "Grids needs to be pre-allocated");
    GMX_RELEASE_ASSERT(gridSizesInXForEachRank.size() == 1 && gridSizesInYForEachRank.size() == 1,
                       "Multi-rank FFT decomposition not implemented with the SYCL MKL backend");

    GMX_ASSERT(checkDeviceBuffer(*realGrid,
                                 realGridSizePadded[XX] * realGridSizePadded[YY] * realGridSizePadded[ZZ]),
               "Real grid buffer is too small for the declared padded size");

    allocateComplexGrid(complexGridSizePadded, realGrid, complexGrid, context);

    GMX_ASSERT(checkDeviceBuffer(*complexGrid,
                                 complexGridSizePadded[XX] * complexGridSizePadded[YY]
                                         * complexGridSizePadded[ZZ] * 2),
               "Complex grid buffer is too small for the declared padded size");

    // GROMACS doesn't use ILP64 for non-GPU interfaces (BLAS, FFT). The GPU/oneAPI interface assumes it.
    // With LP64 on Windows MKL_LONG is 32-bit. Therefore we need to use int64_t and not MKL_LONG.

    // MKL expects row-major
    const std::vector<int64_t> realGridStrides = {
        0, static_cast<int64_t>(realGridSizePadded[YY] * realGridSizePadded[ZZ]), realGridSizePadded[ZZ], 1
    };
    const std::vector<int64_t> complexGridStrides = { 0,
                                                      static_cast<int64_t>(complexGridSizePadded[YY]
                                                                           * complexGridSizePadded[ZZ]),
                                                      complexGridSizePadded[ZZ],
                                                      1 };

    const auto placement = performOutOfPlaceFFT ? PLACEMENT_NOT_INPLACE : PLACEMENT_INPLACE;

#if defined(INTEL_MKL_VERSION) && (INTEL_MKL_VERSION < 20240001) // Intel oneMKL < 2024.1
    try
    {
        using mkl_dft::config_param;
        r2cDescriptor_.set_value(config_param::INPUT_STRIDES, realGridStrides.data());
        r2cDescriptor_.set_value(config_param::OUTPUT_STRIDES, complexGridStrides.data());
        r2cDescriptor_.set_value(config_param::CONJUGATE_EVEN_STORAGE, COMPLEX_COMPLEX_STORAGE);
        r2cDescriptor_.set_value(config_param::PLACEMENT, placement);
        r2cDescriptor_.commit(queue_);
    }
    catch (mkl_exception& exc)
    {
        GMX_THROW(InternalError(
                formatString("MKL failure while configuring R2C descriptor: %s", exc.what())));
    }

    try
    {
        using mkl_dft::config_param;
        c2rDescriptor_.set_value(config_param::INPUT_STRIDES, complexGridStrides.data());
        c2rDescriptor_.set_value(config_param::OUTPUT_STRIDES, realGridStrides.data());
        c2rDescriptor_.set_value(config_param::CONJUGATE_EVEN_STORAGE, COMPLEX_COMPLEX_STORAGE);
        c2rDescriptor_.set_value(config_param::PLACEMENT, placement);
        c2rDescriptor_.commit(queue_);
    }
    catch (mkl_exception& exc)
    {
        GMX_THROW(InternalError(
                formatString("MKL failure while configuring C2R descriptor: %s", exc.what())));
    }
#else // Open-source oneMath or Intel oneMKL >= 2024.1
    // Unify the two descriptors once we get rid of (IN|OUT)PUT_STRIDES API above
    for (auto* descriptor : { &r2cDescriptor_, &c2rDescriptor_ })
    {
        try
        {
            using mkl_dft::config_param;
#    if GMX_GPU_FFT_MKL && INTEL_MKL_VERSION >= 20250000
            descriptor->set_value(config_param::FWD_STRIDES, realGridStrides);
            descriptor->set_value(config_param::BWD_STRIDES, complexGridStrides);
            // config_param::CONJUGATE_EVEN_STORAGE is deprecated, it's complex-complex by default
#    else
            descriptor->set_value(config_param::FWD_STRIDES, realGridStrides.data());
            descriptor->set_value(config_param::BWD_STRIDES, complexGridStrides.data());
            descriptor->set_value(config_param::CONJUGATE_EVEN_STORAGE, COMPLEX_COMPLEX_STORAGE);
#    endif
            descriptor->set_value(config_param::PLACEMENT, placement);
            descriptor->commit(queue_);
        }
        catch (mkl_exception& exc)
        {
            GMX_THROW(InternalError(
                    formatString("MKL failure while configuring FFT descriptor: %s", exc.what())));
        }
    }
#endif
}

Gpu3dFft::ImplSyclMkl::~ImplSyclMkl()
{
    deallocateComplexGrid();
}

void Gpu3dFft::ImplSyclMkl::perform3dFft(gmx_fft_direction dir, CommandEvent* /*timingEvent*/)
{
    float* complexGrid = *complexGrid_.buffer_;
    switch (dir)
    {
        case GMX_FFT_REAL_TO_COMPLEX:
            try
            {
                mkl_dft::compute_forward<Descriptor, float, float>(r2cDescriptor_, realGrid_, complexGrid);
            }
            catch (mkl_exception& exc)
            {
                GMX_THROW(InternalError(
                        formatString("MKL failure while executing R2C transform: %s", exc.what())));
            }
            break;
        case GMX_FFT_COMPLEX_TO_REAL:
            try
            {
                mkl_dft::compute_backward<Descriptor, float, float>(c2rDescriptor_, complexGrid, realGrid_);
            }
            catch (mkl_exception& exc)
            {
                GMX_THROW(InternalError(
                        formatString("MKL failure while executing C2R transform: %s", exc.what())));
            }
            break;
        default:
            GMX_THROW(NotImplementedError("The chosen 3D-FFT case is not implemented on GPUs"));
    }
}

} // namespace gmx
