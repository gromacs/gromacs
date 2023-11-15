#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2022- The GROMACS Authors
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

# Manage VkFFT, GPU FFT library used with OpenCL and SYCL.

function(gmx_manage_vkfft BACKEND_NAME)
    set(GMX_EXTERNAL_VKFFT FALSE CACHE BOOL "Use VkFFT library that is external to GROMACS (ON), or the bundled one (OFF)")
    mark_as_advanced(GMX_EXTERNAL_VKFFT)

    if (NOT GMX_EXTERNAL_VKFFT)
        set(vkfft_DIR ${PROJECT_SOURCE_DIR}/src/external/vkfft)
        set(vkfft_VERSION "internal (1.3.1) with ${BACKEND_NAME} backend" PARENT_SCOPE)
    else()
        find_path(vkfft_DIR
            NAMES vkFFT.h
            HINTS "${VKFFT_INCLUDE_DIR}"
            DOC "vkFFT directory"
        )
        if(NOT vkfft_DIR)
            message(FATAL_ERROR "External VkFFT requested, but could not be found. Please set VKFFT_INCLUDE_DIR to the directory containing vkFFT.h")
        endif()
        set(vkfft_VERSION "external (from ${vkfft_DIR}) with ${BACKEND_NAME} backend" PARENT_SCOPE)
    endif()

    add_library(VkFFT INTERFACE)
    target_include_directories(VkFFT INTERFACE ${vkfft_DIR})

    # The "-Wcast-qual" warning appears when compiling VkFFT for OpenCL, but not for HIP. It cannot be suppressed.
    gmx_target_interface_warning_suppression(VkFFT "-Wno-unused-parameter" HAS_WARNING_NO_UNUSED_PARAMETER)
    gmx_target_interface_warning_suppression(VkFFT "-Wno-unused-variable" HAS_WARNING_NO_UNUSED_VARIABLE)
    gmx_target_interface_warning_suppression(VkFFT "-Wno-newline-eof" HAS_WARNING_NO_NEWLINE_EOF)
    gmx_target_interface_warning_suppression(VkFFT "-Wno-old-style-cast" HAS_WARNING_NO_OLD_STYLE_CAST)
    gmx_target_interface_warning_suppression(VkFFT "-Wno-zero-as-null-pointer-constant" HAS_WARNING_NO_ZERO_AS_NULL_POINTER_CONSTANT)
    gmx_target_interface_warning_suppression(VkFFT "-Wno-unused-but-set-variable" HAS_WARNING_NO_UNUSED_BUT_SET_VARIABLE)
    gmx_target_interface_warning_suppression(VkFFT "-Wno-sign-compare" HAS_WARNING_NO_SIGN_COMPARE)

    if (APPLE)
        # macOS Ventura because `sprintf` was deprecated in favor of `snprintf`.
        gmx_target_interface_warning_suppression(VkFFT "-Wno-deprecated-declarations" HAS_WARNING_NO_DEPRECATED_DECLARATIONS)
    endif()

    # Backend-specific settings and workarounds
    if (BACKEND_NAME STREQUAL "CUDA")
        target_compile_definitions(VkFFT INTERFACE VKFFT_BACKEND=1)
        # This is not ideal, because it uses some random version of CUDA. See #4621.
        find_package(CUDAToolkit REQUIRED)
        target_link_libraries(VkFFT INTERFACE CUDA::cuda_driver CUDA::nvrtc)
        list(APPEND GMX_PUBLIC_LIBRARIES CUDA::cuda_driver) # Workaround for #4902, #4922
        set(GMX_PUBLIC_LIBRARIES ${GMX_PUBLIC_LIBRARIES} PARENT_SCOPE)
        if (GMX_SYCL_DPCPP)
            if(NOT DEFINED ENV{GITLAB_CI}) # Don't warn in CI builds
                message(WARNING "The use of VkFFT with CUDA backend is experimental and not intended for production use")
            endif()
            target_link_libraries(VkFFT INTERFACE CUDA::cudart) # Needed only with DPC++
        endif()
    elseif(BACKEND_NAME STREQUAL "HIP")
        target_compile_definitions(VkFFT INTERFACE VKFFT_BACKEND=2)
        if (GMX_SYCL_DPCPP)
            # HIP does not include hiprtc CMake config prior to version 5.6
            # https://github.com/ROCm-Developer-Tools/HIP/issues/3131
            # Using find_package(HIP) pulls in too many dependencies, in particular clang_rt.
            # Once we require ROCm 5.6 or newer, we can simply do
            # find_package(hiprtc REQUIRED)
            # target_link_libraries(VkFFT INTERFACE hiprtc::hiprtc)
            # But for now, we use our custom cmake/FindHip.cmake module:
            find_package(Hip REQUIRED COMPONENTS hiprtc)
            target_link_libraries(VkFFT INTERFACE Hip::amdhip Hip::hiprtc)
        endif()
        # hipFree is marked `nodiscard` but VkFFT ignores it
        gmx_target_interface_warning_suppression(VkFFT "-Wno-unused-result" HAS_WARNING_NO_UNUSED_RESULT)
    elseif(BACKEND_NAME STREQUAL "OpenCL")
        target_compile_definitions(VkFFT INTERFACE VKFFT_BACKEND=3)
        # The "-Wcast-qual" warning appears when compiling VkFFT for OpenCL, but not for HIP.
        gmx_target_interface_warning_suppression(VkFFT "-Wno-cast-qual" HAS_WARNING_NO_CAST_QUAL)
    elseif(BACKEND_NAME STREQUAL "LevelZero")
        target_compile_definitions(VkFFT INTERFACE VKFFT_BACKEND=4)
    else()
        message(FATAL_ERROR "Unknown VkFFT backend name ${BACKEND_NAME}")
    endif()

endfunction()

