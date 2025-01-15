#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2024- The GROMACS Authors
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







# Builds the interface to plumed and add the linkage to libPlumed

gmx_option_multichoice(GMX_USE_PLUMED
    "Build the PLUMED wrapper with GROMACS"
    AUTO
    AUTO ON OFF)
mark_as_advanced(GMX_USE_PLUMED)

function(gmx_manage_plumed)
    # Create a link target, leave it empty if the plumed option is not active
    add_library(plumedgmx INTERFACE)
    set (GMX_PLUMED_ACTIVE OFF CACHE INTERNAL "Cache entry for PLUMED activation")
    if(WIN32)
        if(GMX_USE_PLUMED STREQUAL "ON")
                message(FATAL_ERROR "PLUMED is not supported on Windows. Reconfigure with -DGMX_USE_PLUMED=OFF.")
        endif()
    elseif(NOT GMX_USE_PLUMED STREQUAL "OFF")
        include(CMakePushCheckState)
        cmake_push_check_state(RESET)
        set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_DL_LIBS})
        include(CheckFunctionExists)
        check_function_exists(dlopen HAVE_DLOPEN)
        cmake_pop_check_state()
        if(HAVE_DLOPEN)
            # Plumed.h  compiled in c++ mode creates a fully inlined interface
            # So we just need to activate the directory in applied_forces
            set(PLUMED_DIR "${CMAKE_SOURCE_DIR}/src/external/plumed")
            target_link_libraries( plumedgmx INTERFACE ${CMAKE_DL_LIBS} )
            # The plumedgmx already exists, now we set it up:
            target_include_directories(plumedgmx SYSTEM INTERFACE $<BUILD_INTERFACE:${PLUMED_DIR}>)
            set (GMX_PLUMED_ACTIVE ON CACHE INTERNAL "Cache entry for PLUMED activation")
        else()
            if(GMX_USE_PLUMED STREQUAL "ON")
                message(FATAL_ERROR "PLUMED needs dlopen or anything equivalent. Reconfigure with -DGMX_USE_PLUMED=OFF.")
            else() # "AUTO"
                message(STATUS "PLUMED needs dlopen or anything equivalent. Disabling support.")
            endif()

        endif()

    endif()

endfunction()
