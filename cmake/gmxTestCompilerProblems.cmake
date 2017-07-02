#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

# Macro that runs through a number of tests for buggy compiler
# versions, or other potential problems.
macro(gmx_test_compiler_problems)

    # Warn if C and C++ compilers do not match
    if(NOT CMAKE_C_COMPILER_ID STREQUAL CMAKE_CXX_COMPILER_ID)
        message(WARNING "The ids of the C and C++ compilers do not match (${CMAKE_C_COMPILER_ID} and ${CMAKE_CXX_COMPILER_ID}, respectively). Mixing different C/C++ compilers can cause problems.")
    endif()
    if(NOT CMAKE_C_COMPILER_VERSION STREQUAL CMAKE_CXX_COMPILER_VERSION)
        message(WARNING "The versions of the C and C++ compilers do not match (${CMAKE_C_COMPILER_VERSION} and ${CMAKE_CXX_COMPILER_VERSION}, respectively). Mixing different C/C++ compilers can cause problems.")
    endif()

    # Note that we've already tested that the compiler works with C++11
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        if(WIN32 AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.5.0")
            message(WARNING "Using Clang on Windows requires Clang 3.5.0")
        endif()
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
        message(WARNING "Currently tested PGI compiler versions (up to 15.7) generate binaries that do not pass all regression test, and the generated binaries are significantly slower than with GCC, ICC or Clang. For now we do not recommend PGI beyond development testing - make sure to run the regressiontests.")
    endif()

endmacro(gmx_test_compiler_problems)
