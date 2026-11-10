#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2026- The GROMACS Authors
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

set(GMX_SZ3_MINIMUM_REQUIRED_VERSION "3.4.0")

# Make a four-state enumeration (auto = use external library if available, else internal)
gmx_option_multichoice(GMX_USE_SZ3
    "How to handle the SZ3 dependency of GROMACS"
    AUTO
    AUTO INTERNAL EXTERNAL NONE)
mark_as_advanced(GMX_USE_SZ3)

# Macro which parses GMX_USE_SZ3 to configure GROMACS to build
# with internal or external SZ3 library
macro(gmx_manage_sz3)
    if(GMX_USE_HDF5)
        if (GMX_USE_SZ3 STREQUAL "AUTO" OR GMX_USE_SZ3 STREQUAL "EXTERNAL")
            if(SZ3_ALREADY_SEARCHED)
                set(SZ3_FIND_QUIETLY ON)
            endif()

            # Once we require CMake 4.x QUIET can be replaced with the more
            # appropriate OPTIONAL
            find_package(SZ3 ${GMX_SZ3_MINIMUM_REQUIRED_VERSION} QUIET COMPONENTS hdf5sz3)

            # For "auto" we configure for an internal build below if
            # the library is not found, for "external" we report an error
            if(GMX_USE_SZ3 STREQUAL "EXTERNAL" AND NOT SZ3_FOUND)
                message(FATAL_ERROR "Building with an external SZ3 library was selected, but a library with its HDF5 filter could not be found. Use GMX_USE_SZ3=internal to build with a GROMACS-managed library.")
            endif()

            set(SZ3_ALREADY_SEARCHED TRUE CACHE INTERNAL
                "True if a search for SZ3 has already been done")
        endif()

        if (SZ3_FOUND)
            # At the time of writing, find_package does not print an
            # affirmative message if SZ3 is found, so we print it ourselves
            if(NOT SZ3_ALREADY_PRINTED_EXTERNAL_STATUS_MESSAGE)
                message(STATUS "Found external SZ3 library (found version ${SZ3_VERSION})")
                set(SZ3_ALREADY_PRINTED_EXTERNAL_STATUS_MESSAGE TRUE CACHE INTERNAL
                    "True if a status message has already been printed")
            endif()
            set(GMX_HAVE_SZ3 1 CACHE INTERNAL "Are we building with SZ3?")

            # Reset status messages for internal library configuration
            set(SZ3_ALREADY_PRINTED_INTERNAL_STATUS_MESSAGE FALSE CACHE INTERNAL "")
        elseif (GMX_USE_SZ3 STREQUAL "AUTO" OR GMX_USE_SZ3 STREQUAL "INTERNAL")
            if(NOT SZ3_ALREADY_PRINTED_INTERNAL_STATUS_MESSAGE)
                message(STATUS "Using internal SZ3 library")
                set(SZ3_ALREADY_PRINTED_INTERNAL_STATUS_MESSAGE TRUE CACHE INTERNAL
                    "True if a status message has already been printed")
            endif()

            # Reset status messages for external library configuration
            set(SZ3_ALREADY_SEARCHED FALSE)
            set(SZ3_ALREADY_PRINTED_EXTERNAL_STATUS_MESSAGE FALSE CACHE INTERNAL "")

            gmx_setup_internal_sz3()
            set(GMX_HAVE_SZ3 1 CACHE INTERNAL "Are we building with SZ3?")
        else()
            set(GMX_HAVE_SZ3 0 CACHE INTERNAL "Are we building with SZ3?")
        endif()
    endif()

    # Mark SZ3-related CMake options as "advanced"
    get_cmake_property(vars VARIABLES)
    foreach (varname ${vars})
        if (varname MATCHES "^SZ3_.*" OR varname MATCHES "^FETCHCONTENT_")
            mark_as_advanced(${varname})
        endif()
    endforeach()
endmacro()

function(gmx_setup_internal_sz3)
    # CMake options for building the SZ3 library
    #
    # These must not modify variables set elsewhere during the CMake invocation,
    # which is handled here by running the setup in function scope.
    set(BUILD_H5Z_FILTER ON)
    set(BUILD_MDZ OFF)
    set(BUILD_PARAVIEW_PLUGIN OFF)
    set(BUILD_SZ3_BINARY OFF)
    set(BUILD_TESTING OFF)
    set(SZ3_INSTALL OFF)
    set(SZ3_DEBUG_TIMINGS OFF)
    # We only link against SZ3core and hdf5sz3, neither of which require OpenMP
    set(CMAKE_DISABLE_FIND_PACKAGE_OpenMP ON)
    set(CMAKE_MESSAGE_LOG_LEVEL WARNING)

    include(FetchContent)
    set(FETCHCONTENT_QUIET ON)

    FetchContent_Declare(SZ3 SOURCE_DIR ${CMAKE_SOURCE_DIR}/src/external/sz3 EXCLUDE_FROM_ALL SYSTEM)
    FetchContent_GetProperties(SZ3)
    if (NOT ${sz3_POPULATED})
        FetchContent_MakeAvailable(SZ3)
    endif()

    if (BUILD_SHARED_LIBS)
        install(TARGETS hdf5sz3
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libraries
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT libraries
        )
    endif()
endfunction()
