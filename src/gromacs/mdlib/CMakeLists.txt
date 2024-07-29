#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2010- The GROMACS Authors
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

add_library(mdlib INTERFACE)

file(GLOB MDLIB_SOURCES *.cpp)
# To avoid listing all the necessary files manually, we will remove SYCL-specific files here:
list(REMOVE_ITEM MDLIB_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/gpuforcereduction_impl.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gpuforcereduction_impl_internal_sycl.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/leapfrog_gpu.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/leapfrog_gpu_internal_sycl.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/lincs_gpu.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/lincs_gpu_internal_sycl.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/settle_gpu.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/settle_gpu_internal_sycl.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/update_constrain_gpu_impl.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/update_constrain_gpu_internal_sycl.cpp)

set(MDLIB_SOURCES ${MDLIB_SOURCES} PARENT_SCOPE)
if(GMX_GPU_CUDA)
    gmx_add_libgromacs_sources(
       gpuforcereduction_impl.cpp
       gpuforcereduction_impl_internal.cu
       leapfrog_gpu.cpp
       leapfrog_gpu_internal.cu
       lincs_gpu.cpp
       lincs_gpu_internal.cu
       settle_gpu.cpp
       settle_gpu_internal.cu
       update_constrain_gpu_impl.cpp
       update_constrain_gpu_internal.cu
       mdgraph_gpu_impl.cu
       )
    _gmx_add_files_to_property(CUDA_SOURCES
       gpuforcereduction_impl.cpp
       leapfrog_gpu.cpp
       lincs_gpu.cpp
       settle_gpu.cpp
       update_constrain_gpu_impl.cpp
       )
endif()

if(GMX_GPU_SYCL)
    gmx_add_libgromacs_sources(
        gpuforcereduction_impl.cpp
        gpuforcereduction_impl_internal_sycl.cpp
        leapfrog_gpu.cpp
        leapfrog_gpu_internal_sycl.cpp
        lincs_gpu.cpp
        lincs_gpu_internal_sycl.cpp
        settle_gpu.cpp
        settle_gpu_internal_sycl.cpp
        update_constrain_gpu_impl.cpp
        update_constrain_gpu_internal_sycl.cpp
    )

    _gmx_add_files_to_property(SYCL_SOURCES
        gpuforcereduction_impl.cpp
        gpuforcereduction_impl_internal_sycl.cpp
        leapfrog_gpu.cpp
        leapfrog_gpu_internal_sycl.cpp
        lincs_gpu.cpp
        lincs_gpu_internal_sycl.cpp
        settle_gpu.cpp
        settle_gpu_internal_sycl.cpp
        update_constrain_gpu_impl.cpp
        update_constrain_gpu_internal_sycl.cpp
    )
endif()

# Source files have the following private module dependencies.
target_link_libraries(mdlib PRIVATE
#                      gmxlib
#                      math
#                      mdtypes
#                      tng_io
                      )

# Public interface for modules, including dependencies and interfaces
#target_include_directories(mdlib PUBLIC
#                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
#target_link_libraries(mdlib PUBLIC
target_link_libraries(mdlib INTERFACE
                      legacy_api
                      pbcutil
                      pulling
                      topology
                      timing
                      )

# TODO: when mdlib is an OBJECT target
#target_link_libraries(mdlib PUBLIC legacy_api)
#target_link_libraries(mdlib PRIVATE common)

# Source files have the following private module dependencies.
#target_link_libraries(mdlib PRIVATE tng_io)
# TODO: Explicitly link specific modules.
#target_link_libraries(mdlib PRIVATE legacy_modules)

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()
