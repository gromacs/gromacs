#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2013,2014,2016, by the GROMACS development team, led by
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
        message(FATAL_ERROR "Language '${LANGUAGE}' is not supported by check_compiler_flag")
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
       message(FATAL_ERROR "Language '${LANGUAGE}' is not supported by check_source_compiles_with_flags")
   endif()
endfunction()

# Helper routine to find flag (from a list) to compile a specific source (in C or C++).
# RESULT_VARIABLE     This will be set when we have found a flag that works
# SOURCE              Source code to test
# LANGUAGE            Specifies the language as "C" or "CXX"
# TOOLCHAIN_FLAGS     Flags needed by the toolchain.
# Args 5 through N    Multiple strings with optimization flags to test
#
# If a compile flag is found, but the project in check_c/cxx_source_compiles
# fails to build, sets SUGGEST_BINUTILS_UPDATE in parent scope to suggest
# that the calling code tell the user about this issue if needed.
FUNCTION(GMX_FIND_FLAG_FOR_SOURCE RESULT_VARIABLE SOURCE LANGUAGE TOOLCHAIN_FLAGS)
    # Insert a blank element last in the list (try without any flags too)
    # This must come last, since some compilers (Intel) might try to emulate
    # emulate AVX instructions with SSE4.1 otherwise.
    foreach(_testflag ${ARGN} "")
        # make valid variable names from the flag string: replace all non-alphanumerical chars
        string(REGEX REPLACE "[^a-zA-Z0-9]+" "_" SIMD_FLAG_WORKS_VARIABLE "${LANGUAGE}_SIMD_FLAG_WORKS_${_testflag}")
        string(REGEX REPLACE "[^a-zA-Z0-9]+" "_" SIMD_COMPILE_WORKS_VARIABLE "${LANGUAGE}_SIMD_COMPILE_WORKS_${_testflag}")

        # Check that the flag itself is fine, and that is does not generate warnings either
        gmx_check_compiler_flag("${_testflag}" "${LANGUAGE}" ${SIMD_FLAG_WORKS_VARIABLE})

        if(${SIMD_FLAG_WORKS_VARIABLE})
            IF(DEFINED ${SIMD_COMPILE_WORKS_VARIABLE})
                # This is a subsequent call to CMake, don't spam the status line
                set(RUN_QUIETLY TRUE)
            endif()
            # Check that we can really compile SIMD source with the
            # full toolchain (this does not catch compiler warnings).
            gmx_check_source_compiles_with_flags("${SOURCE}" "${TOOLCHAIN_FLAGS} ${_testflag}" "${LANGUAGE}" ${SIMD_COMPILE_WORKS_VARIABLE})

            if(NOT ${SIMD_COMPILE_WORKS_VARIABLE})
                if (NOT RUN_QUIETLY)
                    message(STATUS "Compiler flag was valid, but executable did not build - perhaps update the binutils package")
                endif()
                set(SUGGEST_BINUTILS_UPDATE 1 PARENT_SCOPE)
            endif()
        endif()

        if(${SIMD_FLAG_WORKS_VARIABLE} AND ${SIMD_COMPILE_WORKS_VARIABLE})
            set(${RESULT_VARIABLE} "${_testflag}" PARENT_SCOPE)
            break()
        endif()
    endforeach()
endfunction()

# Helper routine to find flag (from a list) to compile a specific source (in both C and C++).
# SOURCE                   Source code to test
# TOOLCHAIN_C_FLAGS        Flags needed by the C toolchain.
# TOOLCHAIN_CXX_FLAGS      Flags needed by the C++ toolchain.
# SIMD_C_FLAGS_VARIABLE    Variable into which to write the C compiler flag needed to compile the source
# SIMD_CXX_FLAGS_VARIABLE  Variable into which to write the C++ compiler flag needed to compile the source
# Args 6 through N         Multiple strings with optimization flags to test
#
# If a compile flag is found, but the project in check_c/cxx_source_compiles
# fails to build, sets SUGGEST_BINUTILS_UPDATE in parent scope to suggest
# that the calling code tell the user about this issue if needed.
function(gmx_find_simd_flags SOURCE TOOLCHAIN_C_FLAGS TOOLCHAIN_CXX_FLAGS SIMD_C_FLAGS_VARIABLE SIMD_CXX_FLAGS_VARIABLE)
    gmx_find_flag_for_source(${SIMD_C_FLAGS_VARIABLE} "${SOURCE}" "C" "${TOOLCHAIN_C_FLAGS}" ${ARGN})
    set(${SIMD_C_FLAGS_VARIABLE} ${${SIMD_C_FLAGS_VARIABLE}} PARENT_SCOPE)
    gmx_find_flag_for_source(${SIMD_CXX_FLAGS_VARIABLE} "${SOURCE}" "CXX" "${TOOLCHAIN_CXX_FLAGS}" ${ARGN})
    set(${SIMD_CXX_FLAGS_VARIABLE} ${${SIMD_CXX_FLAGS_VARIABLE}} PARENT_SCOPE)
    set(SUGGEST_BINUTILS_UPDATE ${SUGGEST_BINUTILS_UPDATE} PARENT_SCOPE)
endfunction()
