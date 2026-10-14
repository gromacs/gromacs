#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2023- The GROMACS Authors
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

# Always create an object-library target, so that libgromacs code
# doesn't need to know whether Colvars support is being compiled.
add_library(colvars_objlib OBJECT
    colvarsoptions.cpp
    )
# Ensure that colvars_objlib can be used in both STATIC and SHARED libraries.
if (BUILD_SHARED_LIBS)
    set_target_properties(colvars_objlib PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

gmx_target_compile_options(colvars_objlib)
target_compile_definitions(colvars_objlib PRIVATE HAVE_CONFIG_H)

if (GMX_HAVE_COLVARS AND GMX_USE_COLVARS STREQUAL "INTERNAL")

    # Create an object library for the Colvars sources. This
    # would be a real library if we were depending on
    # an installed Colvars
    set(COLVARS_DIR "${CMAKE_SOURCE_DIR}/src/external/colvars")
    file(GLOB COLVARS_SOURCES ${COLVARS_DIR}/colvars/*.cpp)
    add_library(colvars_library OBJECT ${COLVARS_SOURCES})

    if(GMX_LIB_MPI)
        target_compile_definitions(colvars_library PUBLIC -DCOLVARS_MPI)
        target_link_libraries(colvars_library PUBLIC MPI::MPI_CXX)
    endif()
    if(GMX_OPENMP)
        target_link_libraries(colvars_library PUBLIC OpenMP::OpenMP_CXX)
    endif()
    if(GMX_TORCH)
        target_compile_options(colvars_library PUBLIC ${TORCH_CXX_FLAGS})
        target_include_directories(colvars_library PUBLIC ${TORCH_INCLUDE_DIRS})
        target_compile_definitions(colvars_library PUBLIC -DCOLVARS_TORCH)
    endif()

    # Set correctly the value of __cplusplus, which MSVC doesn't do by default
    target_compile_options(colvars_library PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/Zc:__cplusplus>)
    # Ensure the target can be used in both STATIC and SHARED libraries.
    set_target_properties(colvars_library PROPERTIES POSITION_INDEPENDENT_CODE ON)

    # Ensure we don't get warnings or linting suggestions when
    # including these headers in GROMACS source code.
    target_include_directories(colvars_library SYSTEM INTERFACE $<BUILD_INTERFACE:${COLVARS_DIR}>)

    target_sources(colvars_objlib PRIVATE
        colvarsMDModule.cpp
        colvarspreprocessor.cpp
        colvarssimulationsparameters.cpp
        colvarsforceprovider.cpp
        colvarproxygromacs.cpp
        )
    target_link_libraries(colvars_objlib PRIVATE
        colvars_library
        common
        domdec
        legacy_modules
        mdtypes
        topology
        utility
        )
else()
    target_sources(colvars_objlib PRIVATE
        colvarsMDModule_stub.cpp
        )
    target_link_libraries(colvars_objlib PRIVATE
        common
        domdec
        legacy_modules
        mdtypes
        topology
        utility
        )
endif()

# Add a target that can pass on transitive requirements like include
# headers or compile definitions (should that be needed)
add_library(colvars_module INTERFACE)
target_link_libraries(colvars_module INTERFACE
    $<BUILD_INTERFACE:colvars_objlib>
    $<TARGET_OBJECTS:colvars_objlib>
    )
if (TARGET colvars_library)
    target_link_libraries(colvars_module INTERFACE
        $<BUILD_INTERFACE:colvars_library>
        $<TARGET_OBJECTS:colvars_library>
        )
endif()

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()
