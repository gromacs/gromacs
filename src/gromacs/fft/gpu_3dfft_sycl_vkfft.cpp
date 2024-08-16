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
 * VkFFT hipSYCL support to GROMACS was contributed by Advanced Micro Devices, Inc.
 * Copyright (c) 2022, Advanced Micro Devices, Inc.  All rights reserved.
 */

/*! \internal \file
 *  \brief Implements GPU 3D FFT routines for hipSYCL using vkFFT.
 *
 *  \author Bálint Soproni <balint@streamhpc.com>
 *  \author Anton Gorenko <anton@streamhpc.com>
 *  \ingroup module_fft
 */

#include "gmxpre.h"

#include "gpu_3dfft_sycl_vkfft.h"

#include "config.h"

#include <vkFFT.h>

#include <string>

#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

#if VKFFT_BACKEND == 1 // CUDA
using NativeDevice = int;
#    if GMX_SYCL_ACPP
static constexpr auto sc_syclBackend = sycl::backend::cuda;
#    elif GMX_SYCL_DPCPP
static constexpr auto sc_syclBackend = sycl::backend::ext_oneapi_cuda;
#    endif
#elif VKFFT_BACKEND == 2 // HIP
using NativeDevice                   = hipDevice_t;
#    if GMX_SYCL_ACPP
static constexpr auto sc_syclBackend = sycl::backend::hip;
#    elif GMX_SYCL_DPCPP
#        include <sycl/ext/oneapi/backend/hip.hpp>
static constexpr auto sc_syclBackend = sycl::backend::ext_oneapi_hip;
#    endif
#else
#    error "VkFFT launcher in GROMACS does not support the selected VkFFT backend"
#endif

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
class Gpu3dFft::ImplSyclVkfft::Impl
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
    NativeDevice       queue_device_;

    DeviceBuffer<float> realGrid_;

    /*! \brief Copy of PME stream
     *
     * This copy is guaranteed by the SYCL standard to work as if
     * it was the original. */
    sycl::queue queue_;
};


Gpu3dFft::ImplSyclVkfft::Impl::Impl(bool allocateGrids,
                                    MPI_Comm /*comm*/,
                                    ArrayRef<const int> gridSizesInXForEachRank,
                                    ArrayRef<const int> gridSizesInYForEachRank,
                                    const int /*nz*/,
                                    bool performOutOfPlaceFFT,
                                    const DeviceContext& /*context*/,
                                    const DeviceStream&  pmeStream,
                                    ivec                 realGridSize,
                                    ivec                 realGridSizePadded,
                                    ivec                 complexGridSizePadded,
                                    DeviceBuffer<float>* realGrid,
                                    DeviceBuffer<float>* complexGrid) :
    realGrid_(*realGrid), queue_(pmeStream.stream())
{
    GMX_RELEASE_ASSERT(allocateGrids == false, "Grids needs to be pre-allocated");
    GMX_RELEASE_ASSERT(gridSizesInXForEachRank.size() == 1 && gridSizesInYForEachRank.size() == 1,
                       "FFT decomposition not implemented with the SYCL VkFFT backend");
    GMX_RELEASE_ASSERT(performOutOfPlaceFFT, "Only out-of-place FFT is implemented in hipSYCL");
    GMX_RELEASE_ASSERT(realGrid, "Bad (null) input real-space grid");
    GMX_RELEASE_ASSERT(complexGrid, "Bad (null) input complex grid");


    configuration_         = {};
    launchParams_          = {};
    application_           = {};
    configuration_.FFTdim  = 3;
    configuration_.size[0] = realGridSize[ZZ];
    configuration_.size[1] = realGridSize[YY];
    configuration_.size[2] = realGridSize[XX];

    configuration_.performR2C  = 1;
    queue_device_              = sycl::get_native<sc_syclBackend>(queue_.get_device());
    configuration_.device      = &queue_device_;
    configuration_.num_streams = 1;

    bufferSize_ = complexGridSizePadded[XX] * complexGridSizePadded[YY] * complexGridSizePadded[ZZ]
                  * sizeof(float) * 2;
    configuration_.bufferSize      = &bufferSize_;
    configuration_.aimThreads      = 64; // Tuned for AMD GCN architecture
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
#if GMX_SYCL_DPCPP
    VkFFTResult result = initializeVkFFT(&application_, configuration_);
    handleFftError(result, "Initializing VkFFT");
#elif GMX_SYCL_ACPP
    queue_.submit([&](sycl::handler& cgh) {
              cgh.hipSYCL_enqueue_custom_operation([=](sycl::interop_handle& gmx_unused h) {
                  VkFFTResult result = initializeVkFFT(&application_, configuration_);
                  handleFftError(result, "Initializing VkFFT");
              });
          }).wait();
#endif
}


Gpu3dFft::ImplSyclVkfft::Impl::~Impl()
{
    deleteVkFFT(&application_);
}

Gpu3dFft::ImplSyclVkfft::ImplSyclVkfft(bool                 allocateGrids,
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

Gpu3dFft::ImplSyclVkfft::~ImplSyclVkfft()
{
    deallocateComplexGrid();
}

template<typename NativeQueue>
static void launchVkFft(const DeviceBuffer<float>& realGrid,
                        const DeviceBuffer<float>& complexGrid,
                        NativeQueue                queue,
                        gmx_fft_direction          fftDirection,
                        VkFFTApplication*          application,
                        VkFFTLaunchParams*         launchParams)
{
    void* d_complexGrid = reinterpret_cast<void*>(complexGrid.buffer_->ptr_);
    void* d_realGrid    = reinterpret_cast<void*>(realGrid.buffer_->ptr_);
    // based on: https://github.com/DTolm/VkFFT/issues/78
    application->configuration.stream = &queue;
    launchParams->inputBuffer         = &d_realGrid;
    launchParams->buffer              = &d_complexGrid;
    VkFFTResult result                = VKFFT_SUCCESS;
    if (fftDirection == GMX_FFT_REAL_TO_COMPLEX)
    {
        result = VkFFTAppend(application, -1, launchParams);
        handleFftError(result, "VkFFT: Real to complex");
    }
    else
    {
        result = VkFFTAppend(application, 1, launchParams);
        handleFftError(result, "VkFFT: Complex to real");
    }
}

void Gpu3dFft::ImplSyclVkfft::perform3dFft(gmx_fft_direction dir, CommandEvent* /*timingEvent*/)
{
#if GMX_SYCL_ACPP // use hipSYCL_enqueue_custom_operation
    impl_->queue_.submit(GMX_SYCL_DISCARD_EVENT[&](sycl::handler & cgh) {
        cgh.hipSYCL_enqueue_custom_operation([=](sycl::interop_handle& h) {
            launchVkFft(impl_->realGrid_,
                        complexGrid_,
                        h.get_native_queue<sc_syclBackend>(),
                        dir,
                        &impl_->application_,
                        &impl_->launchParams_);
        });
    });
#elif GMX_SYCL_DPCPP // submit directly
    launchVkFft(impl_->realGrid_,
                complexGrid_,
                sycl::get_native<sc_syclBackend>(impl_->queue_),
                dir,
                &impl_->application_,
                &impl_->launchParams_);
#endif
}

} // namespace gmx
