/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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
 *  using the double-batched FFT library found at
 *  https://github.com/intel/double-batched-fft-library.
 *  MKL is expected to perform a bit better than bbfft
 *  except for extremely large simulations.
 *
 *  \author Carsten Uphoff <carsten.uphoff@intel.com>
 *  \ingroup module_fft
 */

#include "gmxpre.h"

#include "gpu_3dfft_sycl_bbfft.h"

#include "config.h"

#include <bbfft/bad_configuration.hpp>
#include <bbfft/configuration.hpp>
#include <bbfft/sycl/make_plan.hpp>

#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/devicebuffer_sycl.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

class DeviceContext;

#if (!GMX_SYCL_DPCPP)
#    error This file is only supported with Intel DPC++ compiler
#endif

#include <cstddef>

namespace gmx
{

Gpu3dFft::ImplSyclBbfft::ImplSyclBbfft(bool allocateRealGrid,
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
    Gpu3dFft::Impl::Impl(performOutOfPlaceFFT), realGrid_(*realGrid->buffer_), queue_(pmeStream.stream())
{
    GMX_RELEASE_ASSERT(!allocateRealGrid, "Grids needs to be pre-allocated");
    GMX_RELEASE_ASSERT(gridSizesInXForEachRank.size() == 1 && gridSizesInYForEachRank.size() == 1,
                       "Multi-rank FFT decomposition not implemented with the SYCL backend");

    GMX_ASSERT(checkDeviceBuffer(*realGrid,
                                 realGridSizePadded[XX] * realGridSizePadded[YY] * realGridSizePadded[ZZ]),
               "Real grid buffer is too small for the declared padded size");

    allocateComplexGrid(complexGridSizePadded, realGrid, complexGrid, context);

    GMX_ASSERT(checkDeviceBuffer(*complexGrid,
                                 complexGridSizePadded[XX] * complexGridSizePadded[YY]
                                         * complexGridSizePadded[ZZ] * 2),
               "Complex grid buffer is too small for the declared padded size");

    std::array<size_t, bbfft::max_tensor_dim> shape   = { 1,
                                                        static_cast<size_t>(realGridSize[ZZ]),
                                                        static_cast<size_t>(realGridSize[YY]),
                                                        static_cast<size_t>(realGridSize[XX]),
                                                        1 };
    std::array<size_t, bbfft::max_tensor_dim> rstride = {
        1,
        1,
        static_cast<size_t>(realGridSizePadded[ZZ]),
        static_cast<size_t>(realGridSizePadded[ZZ] * realGridSizePadded[YY]),
        static_cast<size_t>(realGridSizePadded[ZZ] * realGridSizePadded[YY] * realGridSizePadded[XX])
    };
    std::array<size_t, bbfft::max_tensor_dim> cstride = {
        1,
        1,
        static_cast<size_t>(complexGridSizePadded[ZZ]),
        static_cast<size_t>(complexGridSizePadded[ZZ] * complexGridSizePadded[YY]),
        static_cast<size_t>(complexGridSizePadded[ZZ] * complexGridSizePadded[YY] * complexGridSizePadded[XX])
    };

    try
    {
        bbfft::configuration cfg = { 3,
                                     shape,
                                     bbfft::precision::f32,
                                     bbfft::direction::forward,
                                     bbfft::transform_type::r2c,
                                     rstride,
                                     cstride };
        r2cDescriptor_           = bbfft::make_plan(cfg, queue_);
    }
    catch (bbfft::bad_configuration& exc)
    {
        GMX_THROW(InternalError(
                formatString("bbfft failure while configuring R2C descriptor: %s", exc.what())));
    }

    try
    {
        bbfft::configuration cfg = { 3,
                                     shape,
                                     bbfft::precision::f32,
                                     bbfft::direction::backward,
                                     bbfft::transform_type::c2r,
                                     cstride,
                                     rstride };
        c2rDescriptor_           = bbfft::make_plan(cfg, queue_);
    }
    catch (bbfft::bad_configuration& exc)
    {
        GMX_THROW(InternalError(
                formatString("bbfft failure while configuring C2R descriptor: %s", exc.what())));
    }
}

Gpu3dFft::ImplSyclBbfft::~ImplSyclBbfft()
{
    deallocateComplexGrid();
}

void Gpu3dFft::ImplSyclBbfft::perform3dFft(gmx_fft_direction dir, CommandEvent* /*timingEvent*/)
{
    float* complexGrid = *complexGrid_.buffer_;
    switch (dir)
    {
        case GMX_FFT_REAL_TO_COMPLEX: r2cDescriptor_.execute(realGrid_, complexGrid); break;
        case GMX_FFT_COMPLEX_TO_REAL: c2rDescriptor_.execute(complexGrid, realGrid_); break;
        default:
            GMX_THROW(NotImplementedError("The chosen 3D-FFT case is not implemented on GPUs"));
    }
}

} // namespace gmx
