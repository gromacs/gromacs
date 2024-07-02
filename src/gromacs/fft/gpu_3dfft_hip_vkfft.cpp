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

/*
 * VkFFT support to GROMACS was contributed by Advanced Micro Devices, Inc.
 * Copyright (c) 2022, Advanced Micro Devices, Inc.  All rights reserved.
 */

/*! \internal \file
 *  \brief Implements GPU 3D FFT routines for VkFFT with HIP.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 *  \author Paul Bauer <paul.bauer.q@gmail.com>
 *  \ingroup module_fft
 */

#include "gmxpre.h"

#include "gpu_3dfft_hip_vkfft.h"

#include <vkFFT.h>

#include <string>

#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/hardware/device_information.h"
#include "gromacs/utility/arrayref.h"
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
class Gpu3dFft::ImplHipVkFft::Impl
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
    uint64_t           bufferSize_;
    uint64_t           inputBufferSize_;
    hipStream_t        deviceStream_;

    DeviceBuffer<float> realGrid_;
};


Gpu3dFft::ImplHipVkFft::Impl::Impl(bool allocateRealGrid,
                                   MPI_Comm /*comm*/,
                                   ArrayRef<const int> gridSizesInXForEachRank,
                                   ArrayRef<const int> gridSizesInYForEachRank,
                                   const int /*nz*/,
                                   bool performOutOfPlaceFFT,
                                   const DeviceContext& /* context */,
                                   const DeviceStream&  pmeStream,
                                   ivec                 realGridSize,
                                   ivec                 realGridSizePadded,
                                   ivec                 complexGridSizePadded,
                                   DeviceBuffer<float>* realGrid,
                                   DeviceBuffer<float>* complexGrid) :
    realGrid_(*realGrid)
{
    GMX_RELEASE_ASSERT(allocateRealGrid == false, "Grids needs to be pre-allocated");
    GMX_RELEASE_ASSERT(gridSizesInXForEachRank.size() == 1 && gridSizesInYForEachRank.size() == 1,
                       "FFT decomposition not implemented with hipFFT backend");
    GMX_RELEASE_ASSERT(performOutOfPlaceFFT, "Only out-of-place FFT is implemented in HIP");
    GMX_RELEASE_ASSERT(realGrid, "Bad (null) input real-space grid");
    GMX_RELEASE_ASSERT(complexGrid, "Bad (null) input complex grid");

    GMX_RELEASE_ASSERT(realGrid, "Bad (null) input real-space grid");
    GMX_RELEASE_ASSERT(complexGrid, "Bad (null) input complex grid");

    configuration_         = {};
    application_           = {};
    configuration_.FFTdim  = 3;
    configuration_.size[0] = realGridSize[ZZ];
    configuration_.size[1] = realGridSize[YY];
    configuration_.size[2] = realGridSize[XX];

    configuration_.performR2C = 1;
    configuration_.device     = new hipDevice_t;
    gmx::checkDeviceError(hipGetDevice(configuration_.device),
                          "Trying to get hip device for VkFFT failed");
    deviceStream_              = pmeStream.stream();
    configuration_.stream      = &deviceStream_;
    configuration_.num_streams = 1;

    bufferSize_ = complexGridSizePadded[XX] * complexGridSizePadded[YY] * complexGridSizePadded[ZZ]
                  * 2 * sizeof(float);
    configuration_.bufferSize      = &bufferSize_;
    configuration_.aimThreads      = 64;
    configuration_.bufferStride[0] = complexGridSizePadded[ZZ];
    configuration_.bufferStride[1] = complexGridSizePadded[ZZ] * complexGridSizePadded[YY];
    configuration_.bufferStride[2] =
            complexGridSizePadded[ZZ] * complexGridSizePadded[YY] * complexGridSizePadded[XX];
    configuration_.buffer                     = reinterpret_cast<void**>(complexGrid);
    configuration_.isInputFormatted           = 1;
    configuration_.inverseReturnToInputBuffer = 1;
    inputBufferSize_ =
            realGridSizePadded[XX] * realGridSizePadded[YY] * realGridSizePadded[ZZ] * sizeof(float);
    configuration_.inputBufferSize      = &inputBufferSize_;
    configuration_.inputBufferStride[0] = realGridSizePadded[ZZ];
    configuration_.inputBufferStride[1] = realGridSizePadded[ZZ] * realGridSizePadded[YY];
    configuration_.inputBufferStride[2] =
            realGridSizePadded[ZZ] * realGridSizePadded[YY] * realGridSizePadded[XX];
    configuration_.inputBuffer = reinterpret_cast<void**>(realGrid);
    VkFFTResult result         = initializeVkFFT(&application_, configuration_);
    handleFftError(result, "Initializing VkFFT");
}

Gpu3dFft::ImplHipVkFft::Impl::~Impl()
{
    delete configuration_.device;
    deleteVkFFT(&application_);
}

Gpu3dFft::ImplHipVkFft::ImplHipVkFft(bool                 allocateGrids,
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

Gpu3dFft::ImplHipVkFft::~ImplHipVkFft()
{
    deallocateComplexGrid();
}

void Gpu3dFft::ImplHipVkFft::perform3dFft(gmx_fft_direction dir, CommandEvent* /*timingEvent*/)
{
    VkFFTResult result = VKFFT_SUCCESS;
    if (dir == GMX_FFT_REAL_TO_COMPLEX)
    {
        result = VkFFTAppend(&impl_->application_, -1, NULL);
        handleFftError(result, "VkFFT: Real to complex");
    }
    else
    {
        result = VkFFTAppend(&impl_->application_, 1, NULL);
        handleFftError(result, "VdkFFT: Complex to real");
    }
}

} // namespace gmx
