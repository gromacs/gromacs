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

# Set up the module library
add_library(mdtypes INTERFACE)

file(GLOB MDTYPES_SOURCES
    awh_history.cpp
    checkpointdata.cpp
    commrec.cpp
    df_history.cpp
    energyhistory.cpp
    forcebuffers.cpp
    group.cpp
    iforceprovider.cpp
    inputrec.cpp
    interaction_const.cpp
    md_enums.cpp
    multipletimestepping.cpp
    observableshistory.cpp
    observablesreducer.cpp
    state.cpp
    threaded_force_buffer.cpp)

if(GMX_GPU)
    gmx_add_libgromacs_sources(
       state_propagator_data_gpu_impl_gpu.cpp
       )
   if(GMX_GPU_CUDA)
       _gmx_add_files_to_property(CUDA_SOURCES
           state_propagator_data_gpu_impl_gpu.cpp
           )
   endif()
   if(GMX_GPU_SYCL)
       _gmx_add_files_to_property(SYCL_SOURCES
           state_propagator_data_gpu_impl_gpu.cpp
           )
   endif()
   if(GMX_GPU_HIP)
	  gmx_add_libgromacs_sources(
		   state_propagator_data_gpu_impl.cpp
		   )
	endif()
else()
    gmx_add_libgromacs_sources(
      state_propagator_data_gpu_impl.cpp
      )
endif()

# Source files have the following private module dependencies.
target_link_libraries(mdtypes PRIVATE
                      #                      gmxlib
                      #                      math
                      #                      mdtypes
                      #                      tng_io
                      )

# Public interface for modules, including dependencies and interfaces
#target_include_directories(mdtypes PUBLIC
#                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
#target_link_libraries(mdtypes PUBLIC
target_link_libraries(mdtypes INTERFACE
                      legacy_api
                      gpu_utils
                      math
                      )
# Interface conveys transitive dependence on other interface(s)
target_link_libraries(mdtypes INTERFACE utility)

# TODO: when fileio is an OBJECT target
#target_link_libraries(mdtypes PUBLIC legacy_api)
#target_link_libraries(mdtypes PRIVATE common)

# Source files have the following private module dependencies.
#target_link_libraries(mdtypes PRIVATE tng_io)

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()

set(LIBGROMACS_SOURCES ${LIBGROMACS_SOURCES} ${MDTYPES_SOURCES} PARENT_SCOPE)
