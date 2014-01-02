#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2013,2014, by the GROMACS development team, led by
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
# VARIABLE            This will be set when we have found a flag that works
# DESCRIPTION         Text string describing what flag we are trying to find
# SOURCE              Source code to test
#                     The compiler is chosen based on the extension of this file
# FLAGSVAR            Variable (string) to which we should add the correct flag
# Args 5 through N    Multiple strings with acceleration flags to test
FUNCTION(GMX_FIND_CFLAG_FOR_SOURCE VARIABLE DESCRIPTION SOURCE CFLAGSVAR)
    IF(NOT DEFINED ${VARIABLE})
        # Insert a blank element last in the list (try without any flags too)
        # This must come last, since some compilers (Intel) might try to emulate
        # emulate AVX instructions with SSE4.1 otherwise.
        foreach(_testflag ${ARGN} "")
            message(STATUS "Try ${DESCRIPTION} = [${_testflag}]")
            set(CMAKE_REQUIRED_FLAGS "${${CFLAGSVAR}} ${_testflag}")
            # make a valid variable name from the flag string: replace all non-alphanumerical chars
            string(REGEX REPLACE "[^a-zA-Z0-9]+" "_" COMPILE_VARIABLE "C_FLAG_${_testflag}")
            check_c_source_compiles("${SOURCE}" ${COMPILE_VARIABLE})
            if(${${COMPILE_VARIABLE}})
                set(${VARIABLE}_FLAG "${_testflag}" CACHE INTERNAL "${DESCRIPTION}")
                set(${VARIABLE} 1 CACHE INTERNAL "Result of test for ${DESCRIPTION}" FORCE)
                break()
            else(${${COMPILE_VARIABLE}})
                set(${VARIABLE} 0 CACHE INTERNAL "Result of test for ${DESCRIPTION}" FORCE)
            endif(${${COMPILE_VARIABLE}})
        endforeach()
    ENDIF(NOT DEFINED ${VARIABLE})
    IF (${VARIABLE})
        SET (${CFLAGSVAR} "${${CFLAGSVAR}} ${${VARIABLE}_FLAG}" PARENT_SCOPE)
    ENDIF (${VARIABLE})
ENDFUNCTION(GMX_FIND_CFLAG_FOR_SOURCE VARIABLE DESCRIPTION SOURCE CFLAGSVAR)


# Helper routine to find flag (from list) to compile a specific C++ source.
# VARIABLE            This will be set when we have found a flag that works
# DESCRIPTION         Text string describing what flag we are trying to find
# SOURCE              Source code to test
#                     The compiler is chosen based on the extension of this file
# FLAGSVAR            Variable (string) to which we should add the correct flag
# Args 5 through N    Multiple strings with acceleration flags to test
FUNCTION(GMX_FIND_CXXFLAG_FOR_SOURCE VARIABLE DESCRIPTION SOURCE CXXFLAGSVAR)
    IF(NOT DEFINED ${VARIABLE})
        # Insert a blank element last in the list (try without any flags too)
        # This must come last, since some compilers (Intel) might try to
        # emulate AVX instructions with SSE4.1 otherwise.
        foreach(_testflag ${ARGN} "")
            message(STATUS "Try ${DESCRIPTION} = [${_testflag}]")
            set(CMAKE_REQUIRED_FLAGS "${${CXXFLAGSVAR}} ${_testflag}")
            # make a valid variable name from the flag string: replace all non-alphanumerical chars
            string(REGEX REPLACE "[^a-zA-Z0-9]+" "_" COMPILE_VARIABLE "CXX_FLAG_${_testflag}")
            check_cxx_source_compiles("${SOURCE}" ${COMPILE_VARIABLE})
            if(${${COMPILE_VARIABLE}})
                set(${VARIABLE}_FLAG "${_testflag}" CACHE INTERNAL "${DESCRIPTION}")
                set(${VARIABLE} 1 CACHE INTERNAL "Result of test for ${DESCRIPTION}" FORCE)
                break()
            else(${${COMPILE_VARIABLE}})
                set(${VARIABLE} 0 CACHE INTERNAL "Result of test for ${DESCRIPTION}" FORCE)
            endif(${${COMPILE_VARIABLE}})
        endforeach()
    ENDIF(NOT DEFINED ${VARIABLE})
    IF (${VARIABLE})
        SET (${CXXFLAGSVAR} "${${CXXFLAGSVAR}} ${${VARIABLE}_FLAG}" PARENT_SCOPE)
    ENDIF (${VARIABLE})
ENDFUNCTION(GMX_FIND_CXXFLAG_FOR_SOURCE VARIABLE DESCRIPTION SOURCE CXXFLAGSVAR)

