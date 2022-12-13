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

/*
 * OpenCL VkFFT support to GROMACS is based heavily on the SYCL VkFFT backend.
 * View the respective source files in this directory for author attribution.
 */

/*! \internal \file
 *  \brief Implements GPU 3D FFT routines for OpenCL.
 *
 *  \author Philip Turner <philipturner.ar@gmail.com>
 *  \ingroup module_fft
 */

#include "gmxpre.h"

#include "gpu_3dfft_ocl_vkfft.h"

#include "config.h"

#include <vkFFT.h>

#include <string>

#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gmxopencl.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

namespace
{

//! Helper for consistent error handling
void handleFftError(VkFFTResult result, const std::string& msg)
{
    if (result == VKFFT_SUCCESS)
    {
        return;
    }
    GMX_THROW(gmx::InternalError(gmx::formatString(
            "%s: (error code %d - %s)\n", msg.c_str(), result, getVkFFTErrorString(result))));
}
} // namespace

//! Impl class
class Gpu3dFft::ImplOclVkfft::Impl
{
public:
    //! \copydoc Gpu3dFft::Impl::Impl
    Impl(bool allocateGrids,
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
         DeviceBuffer<float>* complexGrid);

    ~Impl();


    VkFFTConfiguration configuration_;
    VkFFTApplication   application_;
    VkFFTLaunchParams  launchParams_;
    uint64_t           bufferSize_;
    uint64_t           inputBufferSize_;
    cl_context         context_;
    cl_device_id       device_;
    cl_command_queue   commandQueue_;

    DeviceBuffer<float> realGrid_;
};


Gpu3dFft::ImplOclVkfft::Impl::Impl(bool allocateGrids,
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
    realGrid_(*realGrid)
{
    GMX_RELEASE_ASSERT(allocateGrids == false, "Grids needs to be pre-allocated");
    GMX_RELEASE_ASSERT(gridSizesInXForEachRank.size() == 1 && gridSizesInYForEachRank.size() == 1,
                       "FFT decomposition not implemented with the OpenCL VkFFT backend");
    GMX_RELEASE_ASSERT(realGrid, "Bad (null) input real-space grid");
    GMX_RELEASE_ASSERT(complexGrid, "Bad (null) input complex grid");

    context_      = context.context();
    commandQueue_ = pmeStream.stream();

    // Fetch OpenCL device from the command queue; this backend assumes only one GPU is present.
    cl_device_id devices[1];
    cl_int       status =
            clGetContextInfo(context_, CL_CONTEXT_DEVICES, sizeof(cl_device_id*), &devices, nullptr);
    GMX_RELEASE_ASSERT(status == CL_SUCCESS, "Could not fetch OpenCL device for VkFFT");
    device_ = devices[0];


    configuration_         = {};
    launchParams_          = {};
    application_           = {};
    configuration_.FFTdim  = 3;
    configuration_.size[0] = realGridSize[ZZ];
    configuration_.size[1] = realGridSize[YY];
    configuration_.size[2] = realGridSize[XX];

    configuration_.performR2C   = 1;
    configuration_.commandQueue = &commandQueue_;
    configuration_.device       = &device_;
    configuration_.context      = &context_;

    bufferSize_ = complexGridSizePadded[XX] * complexGridSizePadded[YY] * complexGridSizePadded[ZZ]
                  * sizeof(float) * 2;
    configuration_.bufferSize      = &bufferSize_;
    configuration_.bufferStride[0] = complexGridSizePadded[ZZ];
    configuration_.bufferStride[1] = complexGridSizePadded[ZZ] * complexGridSizePadded[YY];
    configuration_.bufferStride[2] =
            complexGridSizePadded[ZZ] * complexGridSizePadded[YY] * complexGridSizePadded[XX];

    configuration_.isInputFormatted           = 1;
    configuration_.inverseReturnToInputBuffer = 1;
    inputBufferSize_ =
            realGridSizePadded[XX] * realGridSizePadded[YY] * realGridSizePadded[ZZ] * sizeof(float);
    configuration_.inputBufferSize      = &inputBufferSize_;
    configuration_.inputBufferStride[0] = realGridSizePadded[ZZ];
    configuration_.inputBufferStride[1] = realGridSizePadded[ZZ] * realGridSizePadded[YY];
    configuration_.inputBufferStride[2] =
            realGridSizePadded[ZZ] * realGridSizePadded[YY] * realGridSizePadded[XX];

    VkFFTResult result = initializeVkFFT(&application_, configuration_);
    handleFftError(result, "Initializing VkFFT");
}

Gpu3dFft::ImplOclVkfft::Impl::~Impl()
{
    deleteVkFFT(&application_);
}

Gpu3dFft::ImplOclVkfft::ImplOclVkfft(bool                 allocateGrids,
                                     MPI_Comm             comm,
                                     ArrayRef<const int>  gridSizesInXForEachRank,
                                     ArrayRef<const int>  gridSizesInYForEachRank,
                                     const int            nz,
                                     bool                 performOutOfPlaceFFT,
                                     const DeviceContext& context,
                                     const DeviceStream&  pmeStream,
                                     ivec                 realGridSize,
                                     ivec                 realGridSizePadded,
                                     ivec                 complexGridSizePadded,
                                     DeviceBuffer<float>* realGrid,
                                     DeviceBuffer<float>* complexGrid) :
    Gpu3dFft::Impl::Impl(performOutOfPlaceFFT),
    impl_(std::make_unique<Impl>(allocateGrids,
                                 comm,
                                 gridSizesInXForEachRank,
                                 gridSizesInYForEachRank,
                                 nz,
                                 performOutOfPlaceFFT,
                                 context,
                                 pmeStream,
                                 realGridSize,
                                 realGridSizePadded,
                                 complexGridSizePadded,
                                 realGrid,
                                 complexGrid))
{
    allocateComplexGrid(complexGridSizePadded, realGrid, complexGrid, context);
}

Gpu3dFft::ImplOclVkfft::~ImplOclVkfft()
{
    deallocateComplexGrid();
}

void Gpu3dFft::ImplOclVkfft::perform3dFft(gmx_fft_direction dir, CommandEvent* /*timingEvent*/)
{
    impl_->launchParams_.commandQueue = &impl_->commandQueue_;
    // In the OpenCL backend, `DeviceBuffer` is a wrapper for `cl_mem`.
    impl_->launchParams_.inputBuffer = reinterpret_cast<cl_mem*>(&impl_->realGrid_);
    impl_->launchParams_.buffer      = reinterpret_cast<cl_mem*>(&complexGrid_);
    VkFFTResult result               = VKFFT_SUCCESS;
    if (dir == GMX_FFT_REAL_TO_COMPLEX)
    {
        result = VkFFTAppend(&impl_->application_, -1, &impl_->launchParams_);
        handleFftError(result, "VkFFT: Real to complex");
    }
    else
    {
        result = VkFFTAppend(&impl_->application_, 1, &impl_->launchParams_);
        handleFftError(result, "VkFFT: Complex to real");
    }
}

} // namespace gmx
