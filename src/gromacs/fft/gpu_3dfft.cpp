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
 *  \brief Implements stub GPU 3D FFT routines for CPU-only builds
 *
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 *  \author Gaurav Garg <gaugarg@nvidia.com>
 *  \ingroup module_fft
 */

#include "gmxpre.h"

#include "gpu_3dfft.h"

#include "gpu_3dfft_impl.h"

#if GMX_GPU_CUDA
#    include "gpu_3dfft_cufft.h"
#elif GMX_GPU_HIP
#    if GMX_GPU_FFT_VKFFT
#        include "gpu_3dfft_hip_vkfft.h"
#    else
#        include "gpu_3dfft_hipfft.h"
#    endif
#elif GMX_GPU_OPENCL
#    if GMX_GPU_FFT_VKFFT
#        include "gpu_3dfft_ocl_vkfft.h"
#    else
#        include "gpu_3dfft_ocl.h"
#    endif
#elif GMX_GPU_SYCL
#    include "gpu_3dfft_sycl.h"
#    if GMX_GPU_FFT_MKL || GMX_GPU_FFT_ONEMKL
#        include "gpu_3dfft_sycl_mkl.h"
#    elif GMX_GPU_FFT_BBFFT
#        include "gpu_3dfft_sycl_bbfft.h"
#    elif GMX_GPU_FFT_ROCFFT
#        include "gpu_3dfft_sycl_rocfft.h"
#    elif GMX_GPU_FFT_VKFFT
#        include "gpu_3dfft_sycl_vkfft.h"
#    endif
#endif

#if GMX_USE_Heffte
#    include "gpu_3dfft_heffte.h"
#endif

#if GMX_USE_cuFFTMp
#    include "gpu_3dfft_cufftmp.h"
#endif

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"

namespace gmx
{

// [[noreturn]] attributes must be added in the common headers, so it's easier to silence the warning here
CLANG_DIAGNOSTIC_IGNORE("-Wmissing-noreturn")

Gpu3dFft::Gpu3dFft(FftBackend           backend,
                   bool                 allocateRealGrid,
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
                   DeviceBuffer<float>* complexGrid)
{
#if GMX_GPU_CUDA
    switch (backend)
    {
        case FftBackend::Cufft:
            impl_ = std::make_unique<Gpu3dFft::ImplCuFft>(allocateRealGrid,
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
                                                          complexGrid);
            break;
#    if GMX_USE_cuFFTMp
        case FftBackend::CuFFTMp:
            impl_ = std::make_unique<Gpu3dFft::ImplCuFftMp>(allocateRealGrid,
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
                                                            complexGrid);
            break;
#    endif
        default:
            GMX_RELEASE_ASSERT(backend == FftBackend::HeFFTe_CUDA,
                               "Unsupported FFT backend requested");
    }
#elif GMX_GPU_HIP
    switch (backend)
    {
#    if GMX_GPU_FFT_HIPFFT
        case FftBackend::Hipfft:
            impl_ = std::make_unique<Gpu3dFft::ImplHipFft>(allocateRealGrid,
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
                                                           complexGrid);
            break;
#    elif GMX_GPU_FFT_VKFFT
        case FftBackend::HipVkfft:
            impl_ = std::make_unique<Gpu3dFft::ImplHipVkFft>(allocateRealGrid,
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
                                                             complexGrid);
            break;
#    endif
        default: GMX_THROW(InternalError("Unsupported FFT backend requested"));
    }
#elif GMX_GPU_OPENCL
    switch (backend)
    {
#    if GMX_GPU_FFT_VKFFT
        case FftBackend::OclVkfft:
            impl_ = std::make_unique<Gpu3dFft::ImplOclVkfft>(allocateRealGrid,
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
                                                             complexGrid);
            break;
#    elif GMX_GPU_FFT_CLFFT
        case FftBackend::Ocl:
            impl_ = std::make_unique<Gpu3dFft::ImplOcl>(allocateRealGrid,
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
                                                        complexGrid);
            break;
#    endif
        default: GMX_THROW(InternalError("Unsupported FFT backend requested"));
    }
#elif GMX_GPU_SYCL
    switch (backend)
    {
#    if GMX_GPU_FFT_MKL || GMX_GPU_FFT_ONEMKL
        case FftBackend::SyclMkl: [[fallthrough]];
        case FftBackend::SyclOneMkl:
            impl_ = std::make_unique<Gpu3dFft::ImplSyclMkl>(allocateRealGrid,
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
                                                            complexGrid);
            break;
#    elif GMX_GPU_FFT_BBFFT
        case FftBackend::SyclBbfft:
            impl_ = std::make_unique<Gpu3dFft::ImplSyclBbfft>(allocateRealGrid,
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
                                                              complexGrid);
            break;
#    elif GMX_GPU_FFT_VKFFT
        case FftBackend::SyclVkfft:
            impl_ = std::make_unique<Gpu3dFft::ImplSyclVkfft>(allocateRealGrid,
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
                                                              complexGrid);
            break;
#    elif GMX_GPU_FFT_ROCFFT
        case FftBackend::SyclRocfft:
            impl_ = std::make_unique<Gpu3dFft::ImplSyclRocfft>(allocateRealGrid,
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
                                                               complexGrid);
            break;
#    endif
        case FftBackend::Sycl:
            impl_ = std::make_unique<Gpu3dFft::ImplSycl>(allocateRealGrid,
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
                                                         complexGrid);
            break;
        default:
            if (backend != FftBackend::HeFFTe_Sycl_OneMkl && backend != FftBackend::HeFFTe_Sycl_Rocfft
                && backend != FftBackend::HeFFTe_Sycl_cuFFT)
            {
                GMX_THROW(NotImplementedError("Unsupported FFT backend requested"));
            }
    }
#endif

#if GMX_USE_Heffte
    switch (backend)
    {
        case FftBackend::HeFFTe_CUDA:
#    if GMX_GPU_CUDA
            GMX_RELEASE_ASSERT(heffte::backend::is_enabled<heffte::backend::cufft>::value,
                               "HeFFTe not compiled with CUDA support");
            impl_ = std::make_unique<Gpu3dFft::ImplHeFfte<heffte::backend::cufft>>(
                    allocateRealGrid,
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
                    complexGrid);
#    else
            GMX_RELEASE_ASSERT(
                    false,
                    "HeFFTe_CUDA FFT backend is supported only with GROMACS compiled with CUDA");
#    endif
            break;
        case FftBackend::HeFFTe_Sycl_OneMkl:
#    if GMX_GPU_SYCL && GMX_GPU_FFT_MKL
            GMX_RELEASE_ASSERT(heffte::backend::is_enabled<heffte::backend::onemkl>::value,
                               "HeFFTe was not compiled with oneMKL support");
            impl_ = std::make_unique<Gpu3dFft::ImplHeFfte<heffte::backend::onemkl>>(
                    allocateRealGrid,
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
                    complexGrid);
#    else
            GMX_RELEASE_ASSERT(false,
                               "HeFFTe multi-GPU FFT backend is supported in GROMACS SYCL "
                               "build configurations only with oneMKL, rocFFT, or cuFFT");
#    endif
            break;
        case FftBackend::HeFFTe_Sycl_Rocfft:
#    if GMX_GPU_SYCL && GMX_GPU_FFT_ROCFFT
            GMX_RELEASE_ASSERT(heffte::backend::is_enabled<heffte::backend::rocfft>::value,
                               "HeFFTe was not compiled with rocFFT support");
            impl_ = std::make_unique<Gpu3dFft::ImplHeFfte<heffte::backend::rocfft>>(
                    allocateRealGrid,
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
                    complexGrid);
#    else
            GMX_RELEASE_ASSERT(false,
                               "HeFFTe multi-GPU FFT backend is supported in GROMACS SYCL "
                               "build configurations only with oneMKL, rocFFT, or cuFFT");
#    endif
            break;
        case FftBackend::HeFFTe_Sycl_cuFFT:
#    if GMX_GPU_SYCL && GMX_GPU_FFT_CUFFT
            GMX_RELEASE_ASSERT(heffte::backend::is_enabled<heffte::backend::cufft>::value,
                               "HeFFTe was not compiled with cuFFT support");
            impl_ = std::make_unique<Gpu3dFft::ImplHeFfte<heffte::backend::cufft>>(
                    allocateRealGrid,
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
                    complexGrid);
#    else
            GMX_RELEASE_ASSERT(false,
                               "HeFFTe multi-GPU FFT backend is supported in GROMACS SYCL "
                               "build configurations only with oneMKL, rocFFT, or cuFFT");
#    endif
            break;

        default: GMX_RELEASE_ASSERT(impl_ != nullptr, "Unsupported FFT backend requested");
    }
#endif
}

Gpu3dFft::~Gpu3dFft() = default;

void Gpu3dFft::perform3dFft(gmx_fft_direction dir, CommandEvent* timingEvent)
{
    GMX_RELEASE_ASSERT(impl_ != nullptr, "Cannot run GPU routines in a CPU-only configuration");
    impl_->perform3dFft(dir, timingEvent);
}

CLANG_DIAGNOSTIC_RESET

} // namespace gmx
