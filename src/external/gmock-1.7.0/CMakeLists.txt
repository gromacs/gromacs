#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2011,2012,2013,2014, by the GROMACS development team, led by
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

# As stated in README.Gromacs, this file is not part of GMock, but is written
# specifically for the GROMACS build system from scratch.

include(gmxGetGmockTupleWorkaround)
get_gmock_tuple_workaround(GMOCK_COMPILE_DEFINITIONS)
set(GMOCK_COMPILE_DEFINITIONS "_GNU_SOURCE=1;${GMOCK_COMPILE_DEFINITIONS}")
set(GMOCK_COMPILE_DEFINITIONS ${GMOCK_COMPILE_DEFINITIONS} PARENT_SCOPE)

# GTest/GMock suggest linking with pthreads when available for thread safety
set(CMAKE_THREAD_PREFER_PTHREAD 1)
find_package(Threads)
set(PTHREADS_LIBRARIES)
if (CMAKE_USE_PTHREADS_INIT)
    set(PTHREADS_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
endif()

set(GMOCK_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(GTEST_DIR ${GMOCK_DIR}/gtest)
set(GTEST_SOURCES ${GTEST_DIR}/src/gtest-all.cc)
set(GMOCK_SOURCES ${GMOCK_DIR}/src/gmock-all.cc)

set(GTEST_INCLUDE_DIRS ${GTEST_DIR}/include)
set(GMOCK_INCLUDE_DIRS ${GMOCK_DIR}/include ${GTEST_INCLUDE_DIRS})

include_directories(BEFORE ${GTEST_INCLUDE_DIRS})
include_directories(BEFORE ${GTEST_DIR})
include_directories(BEFORE ${GMOCK_INCLUDE_DIRS})
include_directories(BEFORE ${GMOCK_DIR})

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-Wno-unused-variable HAS_NO_UNUSED_VARIABLE)
if (HAS_NO_UNUSED_VARIABLE)
    set_source_files_properties(${GTEST_SOURCES} PROPERTIES COMPILE_FLAGS "-Wno-unused-variable")
endif()

add_library(gmock STATIC ${UNITTEST_TARGET_OPTIONS} ${GMOCK_SOURCES} ${GTEST_SOURCES})
set_property(TARGET gmock APPEND PROPERTY COMPILE_DEFINITIONS "${GMOCK_COMPILE_DEFINITIONS};GTEST_CAN_STREAM_RESULTS=0")

set(GMOCK_LIBRARIES gmock ${PTHREADS_LIBRARIES} PARENT_SCOPE)
set(GTEST_LIBRARIES ${GMOCK_LIBRARIES} PARENT_SCOPE)
set(GMOCK_INCLUDE_DIRS ${GMOCK_INCLUDE_DIRS} PARENT_SCOPE)
set(GTEST_INCLUDE_DIRS ${GTEST_INCLUDE_DIRS} PARENT_SCOPE)
