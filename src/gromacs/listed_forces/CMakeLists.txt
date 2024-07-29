#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2014- The GROMACS Authors
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

add_library(listed_forces INTERFACE)
gmx_add_libgromacs_sources(
    bonded.cpp
    disre.cpp
    listed_forces_gpu_impl.cpp
    listed_forces.cpp
    listed_internal.cpp
    manage_threading.cpp
    orires.cpp
    pairs.cpp
    position_restraints.cpp
    restcbt.cpp
    )

if(GMX_GPU_CUDA)
    gmx_add_libgromacs_sources(
       listed_forces_gpu_impl_gpu.cpp
       listed_forces_gpu_internal.cu
       )
    _gmx_add_files_to_property(CUDA_SOURCES
      listed_forces_gpu_impl_gpu.cpp
      )
endif()

if (GMX_GPU_SYCL)
    gmx_add_libgromacs_sources(
        listed_forces_gpu_impl_gpu.cpp
        listed_forces_gpu_internal_sycl.cpp
        )
    _gmx_add_files_to_property(SYCL_SOURCES
        listed_forces_gpu_impl_gpu.cpp
        listed_forces_gpu_internal_sycl.cpp
        )
endif()

# Source files have the following private module dependencies.
target_link_libraries(listed_forces PRIVATE
#                      gmxlib
#                      math
#                      mdtypes
                      )

# Public interface for modules, including dependencies and interfaces
#target_include_directories(listed_forces PUBLIC
#                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
#target_link_libraries(listed_forces PUBLIC
target_link_libraries(listed_forces INTERFACE
                      legacy_api
                      )

# TODO: when listed_forces is an OBJECT target
#target_link_libraries(listed_forces PUBLIC legacy_api)
#target_link_libraries(listed_forces PRIVATE common)

# Source files have the following private module dependencies.
#target_link_libraries(listed_forces PRIVATE tng_io)
# TODO: Explicitly link specific modules.
#target_link_libraries(listed_forces PRIVATE legacy_modules)

if (BUILD_TESTING)
     add_subdirectory(tests)
endif()
