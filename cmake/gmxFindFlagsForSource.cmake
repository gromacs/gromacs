#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2013,2014,2016,2017, by the GROMACS development team, led by
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

# Helper function to call the correct language version of the CMake
# check_*_compiler_flag function.
function(gmx_check_compiler_flag FLAGS LANGUAGE RESULT_VARIABLE)
    if (LANGUAGE STREQUAL "C")
        check_c_compiler_flag("${FLAGS}" ${RESULT_VARIABLE})
    elseif (LANGUAGE STREQUAL "CXX")
        check_cxx_compiler_flag("${FLAGS}" ${RESULT_VARIABLE})
    else()
        message(FATAL_ERROR "Language '${LANGUAGE}' is not supported by gmx_check_compiler_flag")
    endif()
endfunction()

# Helper function to call the correct language version of the CMake
# check_*_source_compiles function.
function(gmx_check_source_compiles_with_flags SOURCE FLAGS LANGUAGE RESULT_VARIABLE)
   set(CMAKE_REQUIRED_FLAGS "${FLAGS}")
   if (LANGUAGE STREQUAL "C")
       check_c_source_compiles("${SOURCE}" ${RESULT_VARIABLE})
   elseif (LANGUAGE STREQUAL "CXX")
       check_cxx_source_compiles("${SOURCE}" ${RESULT_VARIABLE})
   else()
       message(FATAL_ERROR "Language '${LANGUAGE}' is not supported by gmx_check_source_compiles_with_flags")
   endif()
endfunction()

# Helper routine to find flag (from a list) to compile a specific source (in C or C++).
# RESULT_VARIABLE           Name of variable to set in the parent scope to true if
#                           we have found a flag that works (which could be "")
# SOURCE                    Source code to test
# LANGUAGE                  Specifies the language as "C" or "CXX"
# TOOLCHAIN_FLAGS_VARIABLE  Name of a variable that contains any flags already known
#                           to be needed by the toolchain (unchanged)
# NEW_FLAGS_VARIABLE        The first working flag will be set to this variable.
# Args 6 through N          Multiple strings with compiler flags to test
#
# If gmx_check_compiler_flag() finds a working compiler flag, but the project in
# gmx_check_source_compiles_with_flags() fails to build source code that needs,
# the flag, this function sets SUGGEST_BINUTILS_UPDATE in the parent scope to
# suggest that the calling code tell the user about this issue if needed.
FUNCTION(GMX_FIND_FLAG_FOR_SOURCE RESULT_VARIABLE SOURCE LANGUAGE TOOLCHAIN_FLAGS_VARIABLE NEW_FLAGS_VARIABLE)
    # Insert a blank element last in the list (ie. try without any flags too)
    # This must come last, since some compilers (Intel) might try to emulate
    # emulate AVX instructions with SSE4.1 otherwise.
    foreach(_testflag ${ARGN} "")
        # make valid variable names from the flag string: replace all non-alphanumerical chars
        string(REGEX REPLACE "[^a-zA-Z0-9]+" "_" FLAG_ACCEPTED_VARIABLE "${LANGUAGE}_${_testflag}_FLAG_ACCEPTED")
        string(REGEX REPLACE "[^a-zA-Z0-9]+" "_" COMPILE_WORKS_VARIABLE "${LANGUAGE}_${_testflag}_COMPILE_WORKS")

        # Check that the flag itself is fine, and that it does not generate warnings either, but don't
        # confuse the user by saying an empty flag was "accepted"
        if("${_testflag}" STREQUAL "")
            set(${FLAG_ACCEPTED_VARIABLE} 1)
            # Give the empty flag variable a better name since we cannot control the CMake status messages
            set(COMPILE_WORKS_VARIABLE "${LANGUAGE}_COMPILE_WORKS_WITHOUT_SPECIAL_FLAGS")
        else()
            gmx_check_compiler_flag("${${TOOLCHAIN_FLAGS_VARIABLE}} ${_testflag}" "${LANGUAGE}" ${FLAG_ACCEPTED_VARIABLE})
        endif()

        if(${FLAG_ACCEPTED_VARIABLE})
            IF(DEFINED ${COMPILE_WORKS_VARIABLE})
                # This is a subsequent call to CMake, don't spam the status line
                set(RUN_QUIETLY TRUE)
            endif()
            # Check that we can really compile source with the full
            # toolchain (this does not catch compiler warnings).
            gmx_check_source_compiles_with_flags("${SOURCE}" "${${TOOLCHAIN_FLAGS_VARIABLE}} ${_testflag}" "${LANGUAGE}" ${COMPILE_WORKS_VARIABLE})

            if(NOT ${COMPILE_WORKS_VARIABLE})
                if (NOT RUN_QUIETLY)
                    if("${_testflag}" STREQUAL "")
                        message(STATUS "Could not find any flag to build test source (this could be due to either the compiler or binutils)")
                    else()
                        message(STATUS "Flag was accepted, but it did not build test source (this could be due to either the compiler or binutils)")
                        set(SUGGEST_BINUTILS_UPDATE 1 PARENT_SCOPE)
                    endif()
                endif()
            endif()
            if (${FLAG_ACCEPTED_VARIABLE} AND ${COMPILE_WORKS_VARIABLE})
                set(${RESULT_VARIABLE} 1 PARENT_SCOPE)
                set(${NEW_FLAGS_VARIABLE} "${_testflag}" PARENT_SCOPE)
                break()
            endif()
        endif()
    endforeach()
    # If no flag has been found, then leaving ${RESULT_VARIABLE} unset
    # will be interpreted by CMake as false.
ENDFUNCTION()

# Helper routine to find a flag (from a list) that will compile a specific source (in both C and C++).
# C_RESULT_VARIABLE             Names a variable that will be set true if a way
#                               to compile the source as C was found
# CXX_RESULT_VARIABLE           Names a variable that will be set true if a way
#                               to compile the source as C++ was found
# SOURCE                        Source code to test
# TOOLCHAIN_C_FLAGS_VARIABLE    As input, names a variable that contains flags needed
#                               by the C toolchain.
# TOOLCHAIN_CXX_FLAGS_VARIABLE  As input, names a variable that contains flags needed
#                               by the C++ toolchain.
# NEW_C_FLAGS_VARIABLE          The first working C flag will be set to this variable
# NEW_CXX_FLAGS_VARIABLE        The first working C++ flag will be set to this variable
# Args 8 through N              Multiple strings with compiler flags to test
#
# If a compile flag is found, but the project in check_c/cxx_source_compiles
# fails to build, sets SUGGEST_BINUTILS_UPDATE in parent scope to suggest
# that the calling code tell the user about this issue if needed.
macro(gmx_find_flags C_RESULT_VARIABLE CXX_RESULT_VARIABLE SOURCE TOOLCHAIN_C_FLAGS_VARIABLE TOOLCHAIN_CXX_FLAGS_VARIABLE C_FLAGS_VARIABLE CXX_FLAGS_VARIABLE)
    gmx_find_flag_for_source(${C_RESULT_VARIABLE} "${SOURCE}" "C" ${TOOLCHAIN_C_FLAGS_VARIABLE} ${C_FLAGS_VARIABLE} ${ARGN})
    gmx_find_flag_for_source(${CXX_RESULT_VARIABLE} "${SOURCE}" "CXX" ${TOOLCHAIN_CXX_FLAGS_VARIABLE} ${CXX_FLAGS_VARIABLE} ${ARGN})
endmacro()
