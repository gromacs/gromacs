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

# Build Colvars library as bundled in a GROMACS worktree; not supporting external linkage yet
gmx_option_multichoice(GMX_USE_COLVARS
    "Build the collective variables (Colvars) library interfaced with GROMACS"
    INTERNAL
    INTERNAL NONE)
mark_as_advanced(GMX_USE_COLVARS)

function(gmx_manage_colvars)
    if(GMX_USE_COLVARS STREQUAL "INTERNAL")
        # Create an object library for the colvars sources
        set(COLVARS_DIR "${CMAKE_SOURCE_DIR}/src/external/colvars")
        file(GLOB COLVARS_SOURCES ${COLVARS_DIR}/*.cpp)
        add_library(colvars_objlib OBJECT ${COLVARS_SOURCES})
        # Set correctly the value of __cplusplus, which MSVC doesn't do by default
        target_compile_options(colvars_objlib PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/Zc:__cplusplus>)
        # Ensure that colvars_objlib can be used in both STATIC and SHARED libraries.
        set_target_properties(colvars_objlib PROPERTIES POSITION_INDEPENDENT_CODE ON)

        # Create an INTERFACE library for colvars with the object library as a dependency
        add_library(colvars INTERFACE)
        target_sources(colvars INTERFACE $<TARGET_OBJECTS:colvars_objlib>)
        target_include_directories(colvars SYSTEM INTERFACE $<BUILD_INTERFACE:${COLVARS_DIR}>)

        if(GMX_OPENMP)
            target_compile_options(colvars_objlib PRIVATE ${OpenMP_CXX_FLAGS})
            target_link_libraries(colvars_objlib PRIVATE OpenMP::OpenMP_CXX)
        endif()

    else()
        # Create a dummy link target so the calling code doesn't need to know
        # whether colvars support is being compiled.
        add_library(colvars INTERFACE)
    endif()
endfunction()
