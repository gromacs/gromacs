#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2021- The GROMACS Authors
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

set(GMX_MUPARSER_REQUIRED_VERSION "2.3.4")

include(gmxOptionUtilities)

# Make a three-state enumeration, defaulting to 
gmx_option_multichoice(GMX_USE_MUPARSER
    "How to handle the muparser dependency of GROMACS"
    INTERNAL
    INTERNAL EXTERNAL NONE)
mark_as_advanced(GMX_USE_MUPARSER)

# Make a fully functional muparser library target that libgromacs can
# depend on regardless of how the user directed muparser support and/or
# linking to work.
function(gmx_manage_muparser)
    if(GMX_USE_MUPARSER STREQUAL "INTERNAL")
        # Use cmake's FetchContent to organize the build, even though
        # the content is already present in src/external. In
        # particular, it sets up an easy call to add_subdirectory().
        include(FetchContent)
        FetchContent_Declare(muparser SOURCE_DIR ${CMAKE_SOURCE_DIR}/src/external/muparser)
        if (NOT muparser_POPULATED)
            if (OpenMP_CXX_FLAGS)
                set(OpenMP_FIND_QUIETLY ON)
            endif()
            FetchContent_Populate(muparser)
        endif()
        set(ENABLE_SAMPLES OFF)
        set(ENABLE_OPENMP ${GMX_OPENMP})
        set(ENABLE_WIDE_CHAR OFF)
        set(BUILD_TESTING OFF)
        add_subdirectory(${muparser_SOURCE_DIR} ${muparser_BINARY_DIR} EXCLUDE_FROM_ALL)
        if (BUILD_SHARED_LIBS)
            # Ensure muparser is in the export set called libgromacs,
            # so that it gets installed along with libgromacs.
            install(TARGETS muparser EXPORT libgromacs)
        endif()
        if (WIN32)
            gmx_target_warning_suppression(muparser /wd4310 HAS_NO_MSVC_CAST_TRUNCATES_CONSTANT_VALUE)
        endif()

        # Hide FETCHCONTENT* CMake variables
        get_cmake_property(_VARS VARIABLES)
        foreach (_VARNAME ${_VARS})
            if (_VARNAME MATCHES "^FETCHCONTENT_")
                mark_as_advanced(${_VARNAME})
            endif()
        endforeach()

        set(HAVE_MUPARSER 1 CACHE INTERNAL "Is muparser found?")
    elseif(GMX_USE_MUPARSER STREQUAL "EXTERNAL")
        # Find an external muparser library.
        find_package(muparser ${GMX_MUPARSER_REQUIRED_VERSION} REQUIRED)
        set(HAVE_MUPARSER 1 CACHE INTERNAL "Is muparser found?")
        get_target_property(_muparser_compile_defs muparser::muparser COMPILE_DEFINITIONS)
        if(_muparser_compile_defs MATCHES ".*_UNICODE.*")
            message(FATAL_ERROR "External muParser library with wide char enabled not supported.")
        endif()
    else()
        # Create a dummy link target so the calling code doesn't need to know
        # whether muparser support is being compiled.
        add_library(muparser::muparser INTERFACE)
        # Add the muparser interface library to the libgromacs Export name, even though
        # we will not be installing any content.
        install(TARGETS muparser EXPORT libgromacs)

        set(HAVE_MUPARSER 0 CACHE INTERNAL "Is muparser found?")
    endif()
    mark_as_advanced(HAVE_MUPARSER)
endfunction()
