#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2015,2016,2017, by the GROMACS development team, led by
# Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
# and including many others, as listed in the AUTHORS file in the
# top-level source directory and at http://www.gromacs.org.
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
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.

# FindHwloc
#
# - Locate headers and libraries for the Portable Hardware Locality (hwloc)
#   library.
#
# Usage: find_package(Hwloc)
#
# FindHwloc defines the following variables:
#
# HWLOC_FOUND          - True if both headers and libraries were located
# HWLOC_INCLUDE_DIRS   - Where to find hwloc headers
# HWLOC_LIBRARIES      - Libraries to link with to enable hwloc usage
# HWLOC_VERSION        - Version (string) of the hwloc library
#

# Function for converting hex version numbers from HWLOC_API_VERSION if necessary
function(HEX2DEC str res)
    string(LENGTH "${str}" len)
    if("${len}" EQUAL 1)
        if("${str}" MATCHES "[0-9]")
            set(${res} "${str}" PARENT_SCOPE)
        elseif( "${str}" MATCHES "[aA]")
            set(${res} 10 PARENT_SCOPE)
        elseif( "${str}" MATCHES "[bB]")
            set(${res} 11 PARENT_SCOPE)
        elseif( "${str}" MATCHES "[cC]")
            set(${res} 12 PARENT_SCOPE)
        elseif( "${str}" MATCHES "[dD]")
            set(${res} 13 PARENT_SCOPE)
        elseif( "${str}" MATCHES "[eE]")
            set(${res} 14 PARENT_SCOPE)
        elseif( "${str}" MATCHES "[fF]")
            set(${res} 15 PARENT_SCOPE)
        else()
            return()
        endif()
    else()
        string(SUBSTRING "${str}" 0 1 str1)
        string(SUBSTRING "${str}" 1 -1 str2)
        hex2dec(${str1} res1)
        hex2dec(${str2} res2)
        math(EXPR val "16 * ${res1} + ${res2}")
        set(${res} "${val}" PARENT_SCOPE)
    endif()
endfunction()

find_path(HWLOC_INCLUDE_DIRS "hwloc.h")
find_library(HWLOC_LIBRARIES "hwloc")

if(HWLOC_INCLUDE_DIRS)
    # If we are not cross-compiling we try to use the hwloc-info program
    if(NOT CMAKE_CROSSCOMPILING)
        find_program(HWLOC_INFO "hwloc-info")

        if(HWLOC_INFO)
            execute_process(COMMAND ${HWLOC_INFO} "--version"
                            RESULT_VARIABLE HWLOC_INFO_RES
                            OUTPUT_VARIABLE HWLOC_INFO_OUT
                            ERROR_VARIABLE  HWLOC_INFO_ERR)

            if(HWLOC_INFO_ERR)
	        message(STATUS "Error executing hwloc-info: ${HWLOC_INFO_ERR}")
            endif()
            string(REGEX MATCH "[0-9]+.*[0-9]+" HWLOC_INFO_OUT "${HWLOC_INFO_OUT}")
            set(HWLOC_VERSION ${HWLOC_INFO_OUT} CACHE STRING "Hwloc library version")
        endif()
    endif()

    if (NOT Hwloc_FIND_QUIETLY)
        message(STATUS "hwloc version: ${HWLOC_VERSION}")
    endif()

    # Parse header if cross-compiling, or if hwloc-info was not found
    if(NOT HWLOC_VERSION)
        # Hwloc is never installed as a framework on OS X, so this should always work.
        file(READ "${HWLOC_INCLUDE_DIRS}/hwloc.h"
             HEADER_CONTENTS LIMIT 16384)
        string(REGEX REPLACE ".*#define HWLOC_API_VERSION (0[xX][0-9a-fA-F]+).*" "\\1"
               HWLOC_API_VERSION "${HEADER_CONTENTS}")
        string(SUBSTRING "${HWLOC_API_VERSION}" 4 2 HEX_MAJOR)
        string(SUBSTRING "${HWLOC_API_VERSION}" 6 2 HEX_MINOR)
        string(SUBSTRING "${HWLOC_API_VERSION}" 8 2 HEX_PATCH)
        hex2dec(${HEX_MAJOR} DEC_MAJOR)
        hex2dec(${HEX_MINOR} DEC_MINOR)
        hex2dec(${HEX_PATCH} DEC_PATCH)
        set(HWLOC_VERSION "${DEC_MAJOR}.${DEC_MINOR}.${DEC_PATCH}" CACHE STRING "Hwloc library version")
    endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Hwloc
                                  REQUIRED_VARS HWLOC_LIBRARIES HWLOC_INCLUDE_DIRS
                                  VERSION_VAR HWLOC_VERSION)

mark_as_advanced(HWLOC_INCLUDE_DIRS HWLOC_LIBRARIES HWLOC_VERSION)

