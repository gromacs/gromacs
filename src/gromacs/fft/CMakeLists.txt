#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2013- The GROMACS Authors
# and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
# Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# https://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at https://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out https://www.gromacs.org.

# Set up the module library
add_library(fft INTERFACE)

gmx_add_libgromacs_sources(
     calcgrid.cpp
     fft.cpp
     fft5d.cpp
     parallel_3dfft.cpp
     )

if (GMX_FFT_FFTPACK)
    gmx_add_libgromacs_sources(
        fft_fftpack.cpp
        ${CMAKE_SOURCE_DIR}/src/external/fftpack/fftpack.cpp)
endif()
if (GMX_FFT_FFTW3 OR GMX_FFT_ARMPL_FFTW3)
    gmx_add_libgromacs_sources(fft_fftw3.cpp)
endif()
if (GMX_FFT_MKL)
    gmx_add_libgromacs_sources(fft_mkl.cpp)
endif()
if(GMX_USE_Heffte)
    gmx_add_libgromacs_sources(gpu_3dfft_heffte.cpp)
    if (GMX_GPU_SYCL)
        _gmx_add_files_to_property(SYCL_SOURCES gpu_3dfft_heffte.cpp)
    endif()
endif()
if (GMX_GPU_CUDA)
    if(GMX_USE_cuFFTMp)
        gmx_add_libgromacs_sources(gpu_3dfft_cufftmp.cpp)
    endif()
    gmx_add_libgromacs_sources(
        # CUDA-specific sources
        gpu_3dfft_cufft.cu
        gpu_3dfft_impl.cpp
        gpu_3dfft.cpp
        )
    _gmx_add_files_to_property(CUDA_SOURCES
        # Must add these files so they can include cuda_runtime.h
        gpu_3dfft.cpp
        gpu_3dfft_impl.cpp
        )
    if(GMX_USE_cuFFTMp)
        _gmx_add_files_to_property(CUDA_SOURCES
            # Must add this files as it requires path set in cuda_include_directories
            gpu_3dfft_cufftmp.cpp
        )
    endif()
elseif (GMX_GPU_HIP)
    if (GMX_GPU_FFT_VKFFT)
        gmx_add_libgromacs_sources(
            # HIP-specific sources
            gpu_3dfft_hip_vkfft.cpp
            gpu_3dfft_impl.cpp
            gpu_3dfft.cpp
            )
    elseif(GMX_GPU_FFT_HIPFFT)
        gmx_add_libgromacs_sources(
            # HIP-specific sources
            gpu_3dfft_hipfft.cpp
            gpu_3dfft_impl.cpp
            gpu_3dfft.cpp
            )
    endif()
    _gmx_add_files_to_property(HIP_SOURCES
        # Must add these files so they can include hip_runtime.h
        gpu_3dfft.cpp
        gpu_3dfft_impl.cpp
        gpu_3dfft_hip_vkfft.cpp
        gpu_3dfft_hipfft.cpp
        )
elseif (GMX_GPU_OPENCL)
    if (GMX_GPU_FFT_VKFFT)
        gmx_add_libgromacs_sources(
            gpu_3dfft_ocl_vkfft.cpp
            gpu_3dfft_impl.cpp
            gpu_3dfft.cpp
            )
    elseif(GMX_GPU_FFT_CLFFT)
        gmx_add_libgromacs_sources(
            gpu_3dfft_ocl.cpp
            gpu_3dfft_impl.cpp
            gpu_3dfft.cpp
            )
    endif()
elseif (GMX_GPU_SYCL)
    set(_sycl_fft_sources gpu_3dfft.cpp gpu_3dfft_impl.cpp gpu_3dfft_sycl.cpp)
    if (GMX_GPU_FFT_MKL OR GMX_GPU_FFT_ONEMKL)
        list(APPEND _sycl_fft_sources gpu_3dfft_sycl_mkl.cpp)
    elseif(GMX_GPU_FFT_BBFFT)
        list(APPEND _sycl_fft_sources gpu_3dfft_sycl_bbfft.cpp)
    elseif(GMX_GPU_FFT_ROCFFT)
        list(APPEND _sycl_fft_sources gpu_3dfft_sycl_rocfft.cpp)
    elseif(GMX_GPU_FFT_VKFFT)
        list(APPEND _sycl_fft_sources gpu_3dfft_sycl_vkfft.cpp)
    endif()
    gmx_add_libgromacs_sources(${_sycl_fft_sources})
    _gmx_add_files_to_property(SYCL_SOURCES ${_sycl_fft_sources})
endif()

# Source files have the following private module dependencies.
target_link_libraries(fft PRIVATE
                      #                      gmxlib
                      #                      math
                      #                      mdtypes
                      #                      tng_io
                      )

# Public interface for modules, including dependencies and interfaces
#target_include_directories(fft PUBLIC
#                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
#target_link_libraries(fft PUBLIC
target_link_libraries(fft INTERFACE
                      legacy_api
                      gpu_utils
                      math
                      timing
                      utility
                      )

if (GMX_GPU_FFT_BBFFT)
    target_link_libraries(fft INTERFACE bbfft::bbfft-sycl)
endif()
if (GMX_GPU_FFT_VKFFT)
    target_link_libraries(fft INTERFACE VkFFT)
endif()

# TODO: when fft is an OBJECT target
#target_link_libraries(fft PUBLIC legacy_api)
#target_link_libraries(fft PRIVATE common)

# Source files have the following private module dependencies.
#target_link_libraries(fft PRIVATE tng_io)
# TODO: Explicitly link specific modules.
#target_link_libraries(fft PRIVATE legacy_modules)

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()
