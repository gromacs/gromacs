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
 *  In DPC++, we use Intel oneMKL to perform the FFT. Currently, we only support binary version
 *  of oneMKL, see #4744.
 */

#include "gmxpre.h"

#include "gpu_3dfft_sycl_mkl.h"

#include "config.h"

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

#if (!GMX_FFT_MKL)
#    error Must use MKL library for FFT when compiling with Intel DPC++ compiler
#endif

#include <cstddef>
#pragma clang diagnostic ignored "-Wsuggest-override" // can be removed when support for 2022.0 is dropped
#pragma clang diagnostic ignored "-Wundefined-func-template"

#include <oneapi/mkl/dfti.hpp>
#include <oneapi/mkl/exceptions.hpp>

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
    catch (oneapi::mkl::exception& exc)
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
    const std::array<int64_t, 4> realGridStrides = {
        0, static_cast<int64_t>(realGridSizePadded[YY] * realGridSizePadded[ZZ]), realGridSizePadded[ZZ], 1
    };
    const std::array<int64_t, 4> complexGridStrides = {
        0,
        static_cast<int64_t>(complexGridSizePadded[YY] * complexGridSizePadded[ZZ]),
        complexGridSizePadded[ZZ],
        1
    };

    const auto placement = performOutOfPlaceFFT ? DFTI_NOT_INPLACE : DFTI_INPLACE;

    try
    {
        using oneapi::mkl::dft::config_param;
        r2cDescriptor_.set_value(config_param::INPUT_STRIDES, realGridStrides.data());
        r2cDescriptor_.set_value(config_param::OUTPUT_STRIDES, complexGridStrides.data());
        r2cDescriptor_.set_value(config_param::CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
        r2cDescriptor_.set_value(config_param::PLACEMENT, placement);
        r2cDescriptor_.commit(queue_);
    }
    catch (oneapi::mkl::exception& exc)
    {
        GMX_THROW(InternalError(
                formatString("MKL failure while configuring R2C descriptor: %s", exc.what())));
    }

    try
    {
        using oneapi::mkl::dft::config_param;
        c2rDescriptor_.set_value(config_param::INPUT_STRIDES, complexGridStrides.data());
        c2rDescriptor_.set_value(config_param::OUTPUT_STRIDES, realGridStrides.data());
        c2rDescriptor_.set_value(config_param::CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
        c2rDescriptor_.set_value(config_param::PLACEMENT, placement);
        c2rDescriptor_.commit(queue_);
    }
    catch (oneapi::mkl::exception& exc)
    {
        GMX_THROW(InternalError(
                formatString("MKL failure while configuring C2R descriptor: %s", exc.what())));
    }
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
                oneapi::mkl::dft::compute_forward(r2cDescriptor_, realGrid_, complexGrid);
            }
            catch (oneapi::mkl::exception& exc)
            {
                GMX_THROW(InternalError(
                        formatString("MKL failure while executing R2C transform: %s", exc.what())));
            }
            break;
        case GMX_FFT_COMPLEX_TO_REAL:
            try
            {
                oneapi::mkl::dft::compute_backward(c2rDescriptor_, complexGrid, realGrid_);
            }
            catch (oneapi::mkl::exception& exc)
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
