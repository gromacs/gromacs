/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
 *  \brief Implements GPU 3D FFT routines for OpenCL.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 *  \ingroup module_fft
 */

#include "gmxpre.h"

#include "gpu_3dfft_ocl.h"

#include <clFFT.h>

#include <array>
#include <vector>

#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gmxopencl.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{
//! Throws the exception on clFFT error
static void handleClfftError(clfftStatus status, const char* msg)
{
    // Supposedly it's just a superset of standard OpenCL errors
    if (status != CLFFT_SUCCESS)
    {
        GMX_THROW(InternalError(formatString("%s: %d", msg, status)));
    }
}

Gpu3dFft::ImplOcl::ImplOcl(bool allocateRealGrid,
                           MPI_Comm /*comm*/,
                           ArrayRef<const int> gridSizesInXForEachRank,
                           ArrayRef<const int> gridSizesInYForEachRank,
                           const int /*nz*/,
                           bool                 performOutOfPlaceFFT,
                           const DeviceContext& context,
                           const DeviceStream&  pmeStream,
                           ivec                 realGridSize,
                           ivec                 realGridSizePadded,
                           ivec                 complexGridSizePadded,
                           DeviceBuffer<float>* realGrid,
                           DeviceBuffer<float>* complexGrid) :
    Gpu3dFft::Impl::Impl(performOutOfPlaceFFT), realGrid_(*realGrid)
{
    GMX_RELEASE_ASSERT(allocateRealGrid == false, "Grids needs to be pre-allocated");
    GMX_RELEASE_ASSERT(gridSizesInXForEachRank.size() == 1 && gridSizesInYForEachRank.size() == 1,
                       "FFT decomposition not implemented with OpenCL backend");

    allocateComplexGrid(complexGridSizePadded, realGrid, complexGrid, context);

    cl_context clContext = context.context();
    commandStreams_.push_back(pmeStream.stream());

    // clFFT expects row-major, so dimensions/strides are reversed (ZYX instead of XYZ)
    std::array<size_t, DIM> realGridDimensions = { size_t(realGridSize[ZZ]),
                                                   size_t(realGridSize[YY]),
                                                   size_t(realGridSize[XX]) };
    std::array<size_t, DIM> realGridStrides    = {
        1, size_t(realGridSizePadded[ZZ]), size_t(realGridSizePadded[YY] * realGridSizePadded[ZZ])
    };
    std::array<size_t, DIM> complexGridStrides = {
        1, size_t(complexGridSizePadded[ZZ]), size_t(complexGridSizePadded[YY] * complexGridSizePadded[ZZ])
    };

    constexpr clfftDim dims = CLFFT_3D;
    handleClfftError(clfftCreateDefaultPlan(&planR2C_, clContext, dims, realGridDimensions.data()),
                     "clFFT planning failure");
    handleClfftError(clfftSetResultLocation(planR2C_, performOutOfPlaceFFT ? CLFFT_OUTOFPLACE : CLFFT_INPLACE),
                     "clFFT planning failure");
    handleClfftError(clfftSetPlanPrecision(planR2C_, CLFFT_SINGLE), "clFFT planning failure");
    constexpr cl_float scale = 1.0;
    handleClfftError(clfftSetPlanScale(planR2C_, CLFFT_FORWARD, scale),
                     "clFFT coefficient setup failure");
    handleClfftError(clfftSetPlanScale(planR2C_, CLFFT_BACKWARD, scale),
                     "clFFT coefficient setup failure");

    // The only difference between 2 plans is direction
    handleClfftError(clfftCopyPlan(&planC2R_, clContext, planR2C_), "clFFT plan copying failure");

    handleClfftError(clfftSetLayout(planR2C_, CLFFT_REAL, CLFFT_HERMITIAN_INTERLEAVED),
                     "clFFT R2C layout failure");
    handleClfftError(clfftSetLayout(planC2R_, CLFFT_HERMITIAN_INTERLEAVED, CLFFT_REAL),
                     "clFFT C2R layout failure");

    handleClfftError(clfftSetPlanInStride(planR2C_, dims, realGridStrides.data()),
                     "clFFT stride setting failure");
    handleClfftError(clfftSetPlanOutStride(planR2C_, dims, complexGridStrides.data()),
                     "clFFT stride setting failure");

    handleClfftError(clfftSetPlanInStride(planC2R_, dims, complexGridStrides.data()),
                     "clFFT stride setting failure");
    handleClfftError(clfftSetPlanOutStride(planC2R_, dims, realGridStrides.data()),
                     "clFFT stride setting failure");

    handleClfftError(clfftBakePlan(planR2C_, commandStreams_.size(), commandStreams_.data(), nullptr, nullptr),
                     "clFFT precompiling failure");
    handleClfftError(clfftBakePlan(planC2R_, commandStreams_.size(), commandStreams_.data(), nullptr, nullptr),
                     "clFFT precompiling failure");

    // TODO: implement solve kernel as R2C FFT callback
    // TODO: disable last transpose (clfftSetPlanTransposeResult)
}

Gpu3dFft::ImplOcl::~ImplOcl()
{
    deallocateComplexGrid();

    clfftDestroyPlan(&planR2C_);
    clfftDestroyPlan(&planC2R_);
}

void Gpu3dFft::ImplOcl::perform3dFft(gmx_fft_direction dir, CommandEvent* timingEvent)
{
    cl_mem                            tempBuffer = nullptr;
    constexpr std::array<cl_event, 0> waitEvents{ {} };

    clfftPlanHandle plan;
    clfftDirection  direction;
    cl_mem *        inputGrids, *outputGrids;
    cl_mem          complexGrid = complexGrid_;

    switch (dir)
    {
        case GMX_FFT_REAL_TO_COMPLEX:
            plan        = planR2C_;
            direction   = CLFFT_FORWARD;
            inputGrids  = &realGrid_;
            outputGrids = &complexGrid;
            break;
        case GMX_FFT_COMPLEX_TO_REAL:
            plan        = planC2R_;
            direction   = CLFFT_BACKWARD;
            inputGrids  = &complexGrid;
            outputGrids = &realGrid_;
            break;
        default:
            GMX_THROW(NotImplementedError("The chosen 3D-FFT case is not implemented on GPUs"));
    }
    handleClfftError(clfftEnqueueTransform(plan,
                                           direction,
                                           commandStreams_.size(),
                                           commandStreams_.data(),
                                           waitEvents.size(),
                                           waitEvents.data(),
                                           timingEvent,
                                           inputGrids,
                                           outputGrids,
                                           tempBuffer),
                     "clFFT execution failure");
}

} // namespace gmx
