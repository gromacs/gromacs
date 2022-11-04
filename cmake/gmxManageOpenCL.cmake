#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2012- The GROMACS Authors
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

# OpenCL required version: 1.2 or newer
set(REQUIRED_OPENCL_MIN_VERSION_MAJOR 1)
set(REQUIRED_OPENCL_MIN_VERSION_MINOR 2)
set(REQUIRED_OPENCL_MIN_VERSION ${REQUIRED_OPENCL_MIN_VERSION_MAJOR}.${REQUIRED_OPENCL_MIN_VERSION_MINOR})

set(GMX_GPU_OPENCL ON)

if(GMX_DOUBLE)
    message(FATAL_ERROR "OpenCL acceleration is not available in double precision")
endif()

# for some reason FindOpenCL checks CUDA_PATH but not CUDA_HOME
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH};$ENV{CUDA_HOME})
find_package(OpenCL)

if (OpenCL_FOUND)
    if (OpenCL_VERSION_STRING VERSION_LESS REQUIRED_OPENCL_MIN_VERSION)
        message(FATAL_ERROR "OpenCL " "${OpenCL_VERSION_STRING}" " is not supported. OpenCL version " "${REQUIRED_OPENCL_MIN_VERSION}" " or newer is required.")
        return ()
    endif()
else ()
    message(FATAL_ERROR "OpenCL not found.")
    return()
endif()

# Tell compiler to hide warnings for comments caused by cl_gl_ext.h on Linux
if (UNIX)
    set(OpenCL_DEFINITIONS ${OpenCL_DEFINITIONS} " -Wno-comment")
endif()

# Yes Virginia, Darwin kernel version 14.4 corresponds to OS X 10.4.
if (APPLE AND ${CMAKE_SYSTEM_VERSION} VERSION_LESS "14.4")
        message(WARNING "OS X prior to version 10.10.4 produces incorrect AMD OpenCL code at runtime. You will not be able to use AMD GPUs on this host unless you upgrade your operating system.")
endif()

if (APPLE)
    # VkFFT is required for OpenCL acceleration; clFFT makes the shader compiler silently fail at runtime.
    set(GMX_GPU_FFT_VKFFT ON)
endif()

if (GMX_GPU_FFT_VKFFT)
    # Use VkFFT with OpenCL back end as header-only library
    set(vkfft_VERSION "1.2.26-b15cb0ca3e884bdb6c901a12d87aa8aadf7637d8")
    add_library(VkFFT INTERFACE)
    set(_backend 3)
    target_compile_definitions(VkFFT INTERFACE VKFFT_BACKEND=${_backend})
    target_include_directories(VkFFT INTERFACE ${CMAKE_PROJECT_ROOT}/src/external/VkFFT)
    
    # The "-Wcast-qual" warning appears when compiling VkFFT for OpenCL, but not for HIP. It cannot be suppressed.
    gmx_target_interface_warning_suppression(VkFFT "-Wno-unused-parameter" HAS_WARNING_NO_UNUSED_PARAMETER)
    gmx_target_interface_warning_suppression(VkFFT "-Wno-unused-variable" HAS_WARNING_NO_UNUSED_VARIABLE)
    gmx_target_interface_warning_suppression(VkFFT "-Wno-newline-eof" HAS_WARNING_NO_NEWLINE_EOF)
    gmx_target_interface_warning_suppression(VkFFT "-Wno-old-style-cast" HAS_WARNING_NO_OLD_STYLE_CAST)
    gmx_target_interface_warning_suppression(VkFFT "-Wno-zero-as-null-pointer-constant" HAS_WARNING_NO_ZERO_AS_NULL_POINTER_CONSTANT)
    gmx_target_interface_warning_suppression(VkFFT "-Wno-unused-but-set-variable" HAS_WARNING_NO_UNUSED_BUT_SET_VARIABLE)
    gmx_target_interface_warning_suppression(VkFFT "-Wno-sign-compare" HAS_WARNING_NO_SIGN_COMPARE)
endif()

add_definitions(${OpenCL_DEFINITIONS})

include_directories(SYSTEM ${OpenCL_INCLUDE_DIRS})

# Ensure the OpenCL implementation is 64-bit, because we only support that;
# Note that this can only be revised if the cl_mem size assumptions made
# (originally in pme-gpu-types.h) are relieved.
if (NOT CMAKE_SIZEOF_VOID_P EQUAL 8)
    message(FATAL_ERROR "The OpenCL implementation is only supported on 64-bit platforms.")
endif()

set(GMX_INSTALL_OCLDIR       ${GMX_INSTALL_GMXDATADIR}/opencl)
