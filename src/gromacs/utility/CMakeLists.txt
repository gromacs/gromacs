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

add_library(utility INTERFACE)

file(GLOB UTILITY_SOURCES *.cpp)
if (GMX_GPU_CUDA)
    set_source_files_properties(cuda_version_information.cu PROPERTIES LANGUAGE CUDA)
    set_property(GLOBAL APPEND PROPERTY GMX_LIBGROMACS_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/cuda_version_information.cu)
endif ()
if (GMX_GPU_HIP)
    set_source_files_properties(hip_version_information_hip.cpp PROPERTIES LANGUAGE HIP)
endif ()
if (GMX_GPU_SYCL)
    _gmx_add_files_to_property(SYCL_SOURCES
                               sycl_version_information.cpp
                               )
endif ()
set(LIBGROMACS_SOURCES ${LIBGROMACS_SOURCES} ${UTILITY_SOURCES} PARENT_SCOPE)

# Public interface for modules, including dependencies and interfaces
#target_include_directories(utility PUBLIC
target_include_directories(utility INTERFACE
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
#target_link_libraries(utility PUBLIC
target_link_libraries(utility INTERFACE
                      legacy_api
                      )
if (GMX_LIB_MPI)
    target_link_libraries(utility INTERFACE MPI::MPI_CXX)
else ()
    target_include_directories(
        utility INTERFACE
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/external/thread_mpi/include>
    )
    target_compile_definitions(utility INTERFACE TMPI_USE_VISIBILITY)
endif ()

# TODO(#3288): when utility is an OBJECT target, explicitly link source dependencies.
# Source files have the following private module dependencies.
#target_link_libraries(utility PRIVATE common)
# TODO(#3288): Explicitly link specific modules as legacy_modules is removed.
#target_link_libraries(utility PRIVATE legacy_modules)

if (BUILD_TESTING)
    add_subdirectory(tests)
endif ()
