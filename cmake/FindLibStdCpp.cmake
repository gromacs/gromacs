#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2019, by the GROMACS development team, led by
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

if (NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND NOT CMAKE_CXX_COMPILER_ID MATCHES "Intel") # Compilers supported
    return()
endif()

# Check whether we are using libstdc++ (compiler could be using libc++ or MSVC)
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("#include <new>
int main() { return __GLIBCXX__; }" USING_LIBSTDCXX)

if (NOT USING_LIBSTDCXX)
    return()
endif()

find_program(GCC_PATH g++)
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
   set(EXTRA_MESSAGE "Clang supports to use libc++ (CXXFLAGS=--stdlib=libc++). With libc++ GCC isn't required. ")
endif()
if (NOT GCC_PATH)
   message(FATAL_ERROR "Couldn't find g++. Please set GCC_PATH, PATH or CMAKE_PREFIX_PATH for cmake to find it."
        "${EXTRA_MESSAGE}")
endif()
execute_process(COMMAND ${GCC_PATH} --version OUTPUT_VARIABLE GCC_VERSION)
if (NOT "${GCC_VERSION}" MATCHES "\\) ([0-9]+)\\.([0-9]+\\.[0-9]+)") #Should never happen
   message(FATAL_ERROR "Couldn't detect g++ version for ${GCC_PATH}. Version output: ${GCC_VERSION}. "
           "Please report to developers. ${EXTRA_MESSAGE}")
endif()
if (${CMAKE_MATCH_1} LESS 5)
   set(GCC_PATH ${GCC_PATH})
   unset(GCC_PATH CACHE)
   message(FATAL_ERROR "Found g++ at ${GCC_PATH}. Its version is ${CMAKE_MATCH_1}.${CMAKE_MATCH_2}. "
           "GROMACS requires at least version 5.1. "
           "Please specify a different g++ using GCC_PATH, PATH or CMAKE_PREFIX_PATH. ${EXTRA_MESSAGE}")
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
   get_filename_component(GCC_PATH "${GCC_PATH}" REALPATH)
   get_filename_component(GCC_PATH "${GCC_PATH}" DIRECTORY) #strip g++
   get_filename_component(GCC_PATH "${GCC_PATH}" DIRECTORY) #strip bin
   if (NOT EXISTS "${GCC_PATH}/include/c++")
       message(FATAL_ERROR "${GCC_PATH}/include/c++ doesn't exist even though it should. Please report to developers.")
   endif()
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --gcc-toolchain=${GCC_PATH}")
else() #Intel
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -gcc-name=${GCC_PATH}")
endif()

check_cxx_source_compiles("#include <iterator>
int main() { int a[2]; std::cbegin(a); }" CXX14_COMPILES)

if (NOT CXX14_COMPILES)
   message(FATAL_ERROR "C++14 test failed to compile. Shouldn't happend. Please report to developers.")
endif()
