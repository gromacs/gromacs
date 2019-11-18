#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2011,2012,2013,2014,2015,2016,2019, by the GROMACS development team, led by
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

if(BUILD_TESTING AND GMX_BUILD_UNITTESTS)
    # GTest/GMock suggest linking with pthreads when available for thread safety
    set(CMAKE_THREAD_PREFER_PTHREAD 1)
    find_package(Threads)
    set(PTHREADS_LIBRARIES)
    if (CMAKE_USE_PTHREADS_INIT)
        set(PTHREADS_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
        # Input to googletest CMake scripts
        set(GTEST_HAS_PTHREAD 1)
        # Output for consumption in GROMACS build tree
        set(GTEST_IS_THREADSAFE 1)
    else()
        # Input to googletest CMake scripts
        set(GTEST_HAS_PTHREAD 0)
        # Output for consumption in GROMACS build tree
        set(GTEST_IS_THREADSAFE 0)
    endif()

    add_subdirectory(googletest EXCLUDE_FROM_ALL)
    target_compile_definitions(gmock PRIVATE
                               _GNU_SOURCE=1
                               GTEST_CAN_STREAM_RESULTS=0)
    if (WIN32)
        target_compile_definitions(gmock PRIVATE _SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING)
        target_compile_definitions(gtest PRIVATE _SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING)
        target_compile_definitions(gtest_main PRIVATE _SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING)
    endif()
    gmx_target_warning_suppression(gmock -Wno-deprecated-copy HAVE_NO_DEPRECATED_COPY)

    set(GTEST_IS_THREADSAFE "${GTEST_IS_THREADSAFE}" PARENT_SCOPE)
endif()
