#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2016,2017, by the GROMACS development team, led by
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

# Helper macros to encapsulate some usage of CTest
#
# This file is intended to contain CTest workarounds and such.
include(CMakeParseArguments)

macro (gmx_ctest_init)
    # Set a default valgrind suppression file.
    # This unfortunately needs to duplicate information from CTest to work as
    # expected...
    #set(MEMORYCHECK_SUPPRESSIONS_FILE
    #    "${CMAKE_SOURCE_DIR}/cmake/legacy_and_external.supp"
    #    CACHE FILEPATH
    #    "File that contains suppressions for the memory checker")
    string(TOUPPER "${CMAKE_BUILD_TYPE}" _cmake_build_type)
    if (_cmake_build_type STREQUAL "ASAN")
        set(MEMORYCHECK_TYPE "AddressSanitizer")
    endif()
    include(CTest)
    if(_cmake_build_type STREQUAL "ASAN")
        # AddressSanitizer support in CTest
        # does not work without this...
        set(_ctest_config_file "${PROJECT_BINARY_DIR}/DartConfiguration.tcl")
        file(STRINGS ${_ctest_config_file} _existing REGEX "^CMakeCommand: ")
        if (NOT _existing)
            file(APPEND ${_ctest_config_file} "\nCMakeCommand: ${CMAKE_COMMAND}\n")
        endif()
    endif()
endmacro()

function (gmx_get_test_prefix_cmd VAR)
    set(_options IGNORE_LEAKS)
    cmake_parse_arguments(ARG "${_options}" "" "" ${ARGN})
    set(_opts "")
    if (ARG_IGNORE_LEAKS OR APPLE)
        list(APPEND _opts "detect_leaks=0")
    endif()
    set(_cmd "")
    if (MEMORYCHECK_TYPE STREQUAL "AddressSanitizer")
        string(REPLACE ";" " " _opts "${_opts}")
        set(_cmd ${PROJECT_SOURCE_DIR}/cmake/with_asan_opts.sh ${_opts} --)
    endif()
    set(${VAR} "${_cmd}" PARENT_SCOPE)
endfunction()
