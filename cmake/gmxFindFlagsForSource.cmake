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

# Helper routine to find flag (from a list) to compile a specific C source.
# VARIABLE            Name of variable to set in the parent scope to true if
#                     we have found a flag that works (and set
#                     ${VARIABLE}_FLAG likewise to contain that flag).
# DESCRIPTION         Text string describing what flag we are trying to find
# SOURCE              Source code to test
#                     The compiler is chosen based on the extension of this file
# CFLAGS              Any flags already known to be needed by the toolchain.
# Args 5 through N    Multiple strings with compiler flags to test
#
# If a compile flag is found, but the project in check_c_source_compiles
# fails to build, sets SUGGEST_BINUTILS_UPDATE in parent scope to suggest
# that the calling code tell the user about this issue if needed.
#
# Note that if you update this function, you probably want to make
# matching changes to the C++ version below.
FUNCTION(GMX_FIND_CFLAG_FOR_SOURCE VARIABLE DESCRIPTION SOURCE CFLAGS)
    if(DEFINED ${VARIABLE})
        # This is a subsequent call to CMake, don't spam the status line
        set(RUN_QUIETLY TRUE)
    endif()
    # TODO This is just here to add a level of indentation, remove in master branch
    if(TRUE)
        # Insert a blank element last in the list (try without any flags too)
        # This must come last, since some compilers (Intel) might try to emulate
        # emulate AVX instructions with SSE4.1 otherwise.
        foreach(_testflag ${ARGN} "")
            if (NOT RUN_QUIETLY)
                message(STATUS "Try ${DESCRIPTION} = [${_testflag}]")
            endif()
            # TODO This sets up implicit context for both the compilation
            # checks below, but the implementation could be improved to
            # make this more explicit, or cache such state.
            set(CMAKE_REQUIRED_FLAGS "${CFLAGS} ${_testflag}")
            # make valid variable names from the flag string: replace all non-alphanumerical chars
            string(REGEX REPLACE "[^a-zA-Z0-9]+" "_" COMPILE_FLAG_VARIABLE "C_FLAG_${_testflag}")
            string(REGEX REPLACE "[^a-zA-Z0-9]+" "_" COMPILE_SOURCE_VARIABLE "C_SOURCE_COMPILES_FLAG_${_testflag}")

            # Check that the flag itself is fine, and that is does not generate warnings either
            check_c_compiler_flag("${_testflag}" ${COMPILE_FLAG_VARIABLE})

            if(${COMPILE_FLAG_VARIABLE})
                # Check that we can compile source (this does not catch compiler warnings)
                check_c_source_compiles("${SOURCE}" ${COMPILE_SOURCE_VARIABLE})
                if(NOT ${COMPILE_SOURCE_VARIABLE})
                    if (NOT RUN_QUIETLY)
                        message(STATUS "Compiler flag was valid, but executable did not build - perhaps update the binutils package")
                    endif()
                    set(SUGGEST_BINUTILS_UPDATE 1 PARENT_SCOPE)
                endif()
            endif()

            if(${COMPILE_FLAG_VARIABLE} AND ${COMPILE_SOURCE_VARIABLE})
                set(${VARIABLE}_FLAG "${_testflag}" PARENT_SCOPE)
                set(${VARIABLE} 1 PARENT_SCOPE)
                break()
            endif()
        endforeach()
        # If no flag has been found, then leaving ${VARIABLE} and
        # ${VARIABLE}_FLAG unset will be interpreted by CMake as
        # false.
    endif()
ENDFUNCTION()


# Helper routine to find flag (from list) to compile a specific C++ source.
# VARIABLE            Name of variable to set in the parent scope to true if
#                     we have found a flag that works (and set
#                     ${VARIABLE}_FLAG likweise to contain that flag).
# DESCRIPTION         Text string describing what flag we are trying to find
# SOURCE              Source code to test
#                     The compiler is chosen based on the extension of this file
# CXXFLAGS            Any flags already known to be needed by the toolchain.
# Args 5 through N    Multiple strings with compiler flags to test
#
# If a compile flag is found, but the project in check_cxx_source_compiles
# fails to build, sets SUGGEST_BINUTILS_UPDATE in parent scope to suggest
# that the calling code tell the user about this issue if needed.
#
# Note that if you update this function, you probably want to make
# matching changes to the C version above.
FUNCTION(GMX_FIND_CXXFLAG_FOR_SOURCE VARIABLE DESCRIPTION SOURCE CXXFLAGS)
    IF(DEFINED ${VARIABLE})
        # This is a subsequent call to CMake, don't spam the status line
        set(RUN_QUIETLY TRUE)
    endif()
    # TODO This is just here to add a level of indentation, remove in master branch
    if(TRUE)
        # Insert a blank element last in the list (try without any flags too)
        # This must come last, since some compilers (Intel) might try to
        # emulate AVX instructions with SSE4.1 otherwise.
        foreach(_testflag ${ARGN} "")
            if (NOT RUN_QUIETLY)
                message(STATUS "Try ${DESCRIPTION} = [${_testflag}]")
            endif()
            # TODO This sets up implicit context for both the compilation
            # checks below, but the implementation could be improved to
            # make this more explicit, or cache such state.
            set(CMAKE_REQUIRED_FLAGS "${CXXFLAGS} ${_testflag}")
            # make valid variable names from the flag string: replace all non-alphanumerical chars
            string(REGEX REPLACE "[^a-zA-Z0-9]+" "_" COMPILE_FLAG_VARIABLE "CXX_FLAG_${_testflag}")
            string(REGEX REPLACE "[^a-zA-Z0-9]+" "_" COMPILE_SOURCE_VARIABLE "CXX_SOURCE_COMPILES_FLAG_${_testflag}")
            
            # Check that the flag itself is fine, and that is does not generate warnings either
            check_cxx_compiler_flag("${_testflag}" ${COMPILE_FLAG_VARIABLE})

            if(${COMPILE_FLAG_VARIABLE})
                # Check that we can compile source (this does not catch compiler warnings)
                check_cxx_source_compiles("${SOURCE}" ${COMPILE_SOURCE_VARIABLE})
                if(NOT ${COMPILE_SOURCE_VARIABLE})
                    if (NOT RUN_QUIETLY)
                        message(STATUS "Compiler flag was valid, but executable did not build - perhaps update the binutils package")
                    endif()
                    set(SUGGEST_BINUTILS_UPDATE 1 PARENT_SCOPE)
                endif()
            endif()

            if(${COMPILE_FLAG_VARIABLE} AND ${COMPILE_SOURCE_VARIABLE})
                set(${VARIABLE}_FLAG "${_testflag}" PARENT_SCOPE)
                set(${VARIABLE} 1 PARENT_SCOPE)
                break()
            endif()
        endforeach()
        # If no flag has been found, then leaving ${VARIABLE} and
        # ${VARIABLE}_FLAG unset will be interpreted by CMake as
        # false.
    endif()
ENDFUNCTION()

