#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2019- The GROMACS Authors
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
add_library(nbnxm INTERFACE)

add_subdirectory(kernels_simd_4xm)
add_subdirectory(kernels_simd_2xmm)

file(GLOB NBNXM_SOURCES
    # Source files
    atomdata.cpp
    exclusionchecker.cpp
    freeenergydispatch.cpp
    grid.cpp
    gridset.cpp
    kernel_common.cpp
    kerneldispatch.cpp
    nbnxm.cpp
    nbnxm_geometry.cpp
    nbnxm_setup.cpp
    pairlist.cpp
    pairlistparams.cpp
    pairlistset.cpp
    pairlistwork.cpp
    pairlist_simd_kernel.cpp
    pairlist_tuning.cpp
    pairsearch.cpp
    prunekerneldispatch.cpp
    # Reference kernel source files
    kernels_reference/kernel_gpu_ref.cpp
    kernels_reference/kernel_ref.cpp
    kernels_reference/kernel_ref_prune.cpp
    simd_energy_accumulator.cpp
    simd_prune_kernel.cpp
    # Benchmark source files
    # TODO these should not be in libgromacs
    benchmark/bench_setup.cpp
    benchmark/bench_system.cpp
    )

if(GMX_GPU_CUDA)
    add_subdirectory(cuda)
    gmx_add_libgromacs_sources(nbnxm_gpu_data_mgmt.cpp nbnxm_gpu_buffer_ops.cpp)
    _gmx_add_files_to_property(CUDA_SOURCES
        nbnxm_gpu_data_mgmt.cpp
        nbnxm_gpu_buffer_ops.cpp
        )
endif()

if(GMX_GPU_OPENCL)
    add_subdirectory(opencl)
    set(NBNXM_OPENCL_KERNELS ${NBNXM_OPENCL_KERNELS} PARENT_SCOPE)
    gmx_add_libgromacs_sources(nbnxm_gpu_data_mgmt.cpp nbnxm_gpu_buffer_ops.cpp)
endif()

if(GMX_GPU_SYCL)
    add_subdirectory(sycl)
    gmx_add_libgromacs_sources(nbnxm_gpu_data_mgmt.cpp nbnxm_gpu_buffer_ops.cpp)
    _gmx_add_files_to_property(SYCL_SOURCES nbnxm_gpu_data_mgmt.cpp nbnxm_gpu_buffer_ops.cpp)
endif()

set(LIBGROMACS_SOURCES ${LIBGROMACS_SOURCES} ${NBNXM_SOURCES} PARENT_SCOPE)


# Source files have the following private module dependencies.
target_link_libraries(nbnxm PRIVATE
        )

# Public interface for modules, including dependencies and interfaces
#target_include_directories(nbnxm PUBLIC
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
#target_link_libraries(nbnxm PUBLIC
target_link_libraries(nbnxm INTERFACE
        legacy_api
        gpu_utils
        )

# TODO: when nbnxm is an OBJECT target
#target_link_libraries(nbnxm PUBLIC legacy_api)
#target_link_libraries(nbnxm PRIVATE common)


if (BUILD_TESTING)
    add_subdirectory(tests)
endif()
