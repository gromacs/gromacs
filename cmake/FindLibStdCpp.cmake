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

# For compilers which might require libstdc++ (Clang and Intel), find it and set CMAKE_CXX_FLAGS.
#
# Does nothing if compiler includes std-library (e.g. GCC) or compiler uses different std-library
# (either because of different defaults (e.g. on MacOS) or user flags (e.g. -stdlib=libc++)).
# The heuristic by the compiler of how to find libstdc++ is ignored. So are any user provided flags
# CXXFLAGS for the location of libstdc++. The user can choose the libstdc++ by setting
# GMX_GPLUSPLUS_PATH, PATH or CMAKE_PREFIX_PATH to make sure the correct the g++ is found.
# Gives error if no g++ is found or the g++ found isn't new enough (5.1 is required).
# The location of g++ is cached as GMX_GPLUSPLUS_PATH making sure that the same libstdc++ for different
# builds using the same cache file.

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

find_program(GMX_GPLUSPLUS_PATH g++)
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
   set(EXTRA_MESSAGE
       "Clang supports using libc++ (CXXFLAGS=--stdlib=libc++), and if so there will be no need to find g++.")
endif()
if (NOT GMX_GPLUSPLUS_PATH)
   message(FATAL_ERROR "Couldn't find g++. Please set GMX_GPLUSPLUS_PATH, PATH or CMAKE_PREFIX_PATH "
           "for cmake to find it. ${EXTRA_MESSAGE}")
endif()
execute_process(COMMAND ${GMX_GPLUSPLUS_PATH} -dumpfullversion -dumpversion OUTPUT_VARIABLE GCC_VERSION
                OUTPUT_STRIP_TRAILING_WHITESPACE)
if (NOT "${GCC_VERSION}" MATCHES "^[0-9]+\\.[0-9]+\\.[0-9]+$") #Should never happen
   message(FATAL_ERROR "Couldn't detect g++ version for ${GMX_GPLUSPLUS_PATH}. Version output: ${GCC_VERSION}. "
           "Please report to developers. ${EXTRA_MESSAGE}")
endif()
if (${GCC_VERSION} VERSION_LESS 5.1)
   set(GMX_GPLUSPLUS_PATH ${GMX_GPLUSPLUS_PATH})
   unset(GMX_GPLUSPLUS_PATH CACHE)
   message(FATAL_ERROR "Found g++ at ${GMX_GPLUSPLUS_PATH}. Its version is ${GCC_VERSION}. "
           "GROMACS requires at least version 5.1. "
           "Please specify a different g++ using GMX_GPLUSPLUS_PATH, PATH or CMAKE_PREFIX_PATH. ${EXTRA_MESSAGE}")
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
   get_filename_component(GMX_GPLUSPLUS_PATH "${GMX_GPLUSPLUS_PATH}" REALPATH)
   get_filename_component(GMX_GPLUSPLUS_PATH "${GMX_GPLUSPLUS_PATH}" DIRECTORY) #strip g++
   get_filename_component(GMX_GPLUSPLUS_PATH "${GMX_GPLUSPLUS_PATH}" DIRECTORY) #strip bin
   if (NOT EXISTS "${GMX_GPLUSPLUS_PATH}/include/c++")
       message(FATAL_ERROR "${GMX_GPLUSPLUS_PATH}/include/c++ doesn't exist even though it should. "
               "Please report to developers.")
   endif()
   #Appending flags overwrite potential user provided flags
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --gcc-toolchain=${GMX_GPLUSPLUS_PATH}")
else() #Intel
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -gcc-name=${GMX_GPLUSPLUS_PATH}")
endif()

# Test a feature which was added in libstdc++ 5
check_cxx_source_compiles("#include <iterator>
int main() { int a[2]; std::cbegin(a); }" CXX14_COMPILES)

if (NOT CXX14_COMPILES)
   message(FATAL_ERROR "C++14 test failed to compile. Shouldn't happen. Please report to developers.")
endif()
