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

add_library(hardware INTERFACE)
gmx_add_libgromacs_sources(
    cpuinfo.cpp
    detecthardware.cpp
    device_management_common.cpp
    hardwaretopology.cpp
    printhardware.cpp
    identifyavx512fmaunits.cpp
    simd_support.cpp
    )

if(GMX_GPU_OPENCL)
    gmx_add_libgromacs_sources(
        device_management_ocl.cpp
        )
elseif(GMX_GPU_CUDA)
    gmx_add_libgromacs_sources(
        device_management.cu
        )
    _gmx_add_files_to_property(CUDA_SOURCES
        # Must add these files so they can include device_information.h
        device_management_common.cpp
        detecthardware.cpp
	    )
elseif(GMX_GPU_HIP)
    gmx_add_libgromacs_sources(
        device_management_hip.cpp
        )
    _gmx_add_files_to_property(HIP_SOURCES
        # Must add these files so they can include device_information.h
        device_management_common.cpp
        device_management_hip.cpp
        detecthardware.cpp
	    )
elseif(GMX_GPU_SYCL)
    gmx_add_libgromacs_sources(
        device_management_sycl.cpp
        device_management_sycl_intel_device_ids.cpp
        )
    _gmx_add_files_to_property(SYCL_SOURCES
        device_management_sycl.cpp
        # Must add these files so they can include device_information.h
        device_management.cpp
        device_management_common.cpp
        detecthardware.cpp
        )
else()
    gmx_add_libgromacs_sources(
        device_management.cpp
        )
endif()

# Source files have the following private module dependencies.
target_link_libraries(hardware PRIVATE
                      )

# Public interface for modules, including dependencies and interfaces
#target_include_directories(hardware PUBLIC
#                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
#target_link_libraries(hardware PUBLIC
target_link_libraries(hardware INTERFACE
                      legacy_api
                      )
# The interface conveys transitive dependence on other module(s)
target_link_libraries(hardware INTERFACE utility)

# TODO: when hardware is an OBJECT target
#target_link_libraries(hardware PUBLIC legacy_api)
#target_link_libraries(hardware PRIVATE common)

# Source files have the following private module dependencies.
#target_link_libraries(hardware PRIVATE tng_io)
# TODO: Explicitly link specific modules.
#target_link_libraries(hardware PRIVATE legacy_modules)

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()
