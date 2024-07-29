#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2015- The GROMACS Authors
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

# Note, source files should be listed in alphabetical order, to reduce
# the incidence of textual clashes when adding/moving files that
# otherwise make the end of the list a hotspot.

add_library(gpu_utils INTERFACE)

gmx_add_libgromacs_sources(
        clfftinitializer.cpp
        device_stream_manager.cpp
        gpu_utils.cpp
        gpueventsynchronizer_helpers.cpp
        hostallocator.cpp
        nvshmem_utils.cpp
        )
if(GMX_GPU_OPENCL)
    gmx_add_libgromacs_sources(
        device_context_ocl.cpp
        device_stream_ocl.cpp
        gpu_utils_impl.cpp
        pmalloc.cpp
        ocl_compiler.cpp
        ocl_caching.cpp
        oclutils.cpp
        )
elseif(GMX_GPU_CUDA)
    gmx_add_libgromacs_sources(
        device_context.cpp
        device_stream.cu
        gpu_utils.cu
        pmalloc.cu
        )
    _gmx_add_files_to_property(CUDA_SOURCES
        device_stream_manager.cpp
        )
elseif(GMX_GPU_HIP)
    gmx_add_libgromacs_sources(
        device_context.cpp
        device_stream_hip.cpp
        gpu_utils_hip.cpp
        pmalloc_hip.cpp
        )
    _gmx_add_files_to_property(HIP_SOURCES
        device_stream_hip.cpp
        gpu_utils_hip.cpp
        pmalloc_hip.cpp
        )
elseif(GMX_GPU_SYCL)
    gmx_add_libgromacs_sources(
        devicebuffer_sycl.cpp
        device_context_sycl.cpp
        device_stream_sycl.cpp
        gpu_utils_impl.cpp
        pmalloc_sycl.cpp
        )
    _gmx_add_files_to_property(SYCL_SOURCES
        devicebuffer_sycl.cpp
        device_context_manager.cpp
        device_context_sycl.cpp
        device_stream_manager.cpp
        device_stream_sycl.cpp
        pmalloc_sycl.cpp
        )
else()
    gmx_add_libgromacs_sources(
        device_context.cpp
        device_stream.cpp
        gpu_utils_impl.cpp
        pmalloc.cpp
    )
endif()

# Source files have the following private module dependencies.
target_link_libraries(gpu_utils PRIVATE
#                      gmxlib
#                      math
#                      mdtypes
#                      tng_io
                      )

# Public interface for modules, including dependencies and interfaces
#target_include_directories(gpu_utils PUBLIC
#                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
#target_link_libraries(gpu_utils PUBLIC
target_link_libraries(gpu_utils INTERFACE
                      legacy_api
                      math
                      )

# TODO: when gpu_utils is an OBJECT target
#target_link_libraries(gpu_utils PUBLIC legacy_api)
#target_link_libraries(gpu_utils PRIVATE common)

# Source files have the following private module dependencies.
#target_link_libraries(gpu_utils PRIVATE tng_io)
# TODO: Explicitly link specific modules.
#target_link_libraries(gpu_utils PRIVATE legacy_modules)

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()
