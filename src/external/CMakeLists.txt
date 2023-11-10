#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2011- The GROMACS Authors
# and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
# Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
# https://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at https://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out https://www.gromacs.org.

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

    option(BUILD_GMOCK "googletest: Builds the googlemock subproject. GROMACS uses googlemock." ON)
    option(INSTALL_GTEST "googletest: Enable installation of googletest. GROMACS uses it internally; installation not needed." OFF)
    add_subdirectory(googletest EXCLUDE_FROM_ALL)
    mark_as_advanced(BUILD_GMOCK)
    mark_as_advanced(INSTALL_GTEST)
    target_compile_definitions(gmock
        PRIVATE
            _GNU_SOURCE=1
            GTEST_CAN_STREAM_RESULTS=0
        PUBLIC
            GTEST_LANG_CXX11
            )
    target_compile_definitions(gtest PUBLIC GTEST_LANG_CXX11)
    if(BUILD_SHARED_LIBS)
        foreach(_lib gmock gtest)
            target_compile_definitions(${_lib}
                PRIVATE
                    GTEST_CREATE_SHARED_LIBRARY=1
                INTERFACE
                    GTEST_LINKED_AS_SHARED_LIBRARY=1
                )
        endforeach()
    endif()
    if (CYGWIN)
        # Ensure GoogleTest and things that use it can find
        # POSIX things needed by GoogleTest
        target_compile_definitions(gtest PUBLIC _POSIX_C_SOURCE=200809L)
    endif()
    gmx_target_warning_suppression(gmock -Wno-deprecated-copy HAVE_NO_DEPRECATED_COPY)
    if (CMAKE_BUILD_TYPE MATCHES "Debug")
        gmx_target_warning_suppression(gtest -Wno-debug-disables-optimization HAVE_NO_DEBUG_DISABLES_OPTIMIZATION)
        gmx_target_warning_suppression(gmock -Wno-debug-disables-optimization HAVE_NO_DEBUG_DISABLES_OPTIMIZATION)
    endif()
    if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        # GoogleTest disables this warning for IntelLLVM, but old CMake versions identify it as plain Clang
        # Remove when CMake 3.20+ is required
        # https://github.com/google/googletest/issues/4315
        gmx_target_warning_suppression(gtest -Wno-implicit-float-size-conversion HAVE_NO_IMPLICIT_FLOAT_SIZE_CONVERSION)
        gmx_target_warning_suppression(gtest -Wno-implicit-int-float-conversion HAVE_NO_IMPLICIT_INT_FLOAT_CONVERSION)
    endif()

    set(GTEST_IS_THREADSAFE "${GTEST_IS_THREADSAFE}" PARENT_SCOPE)
endif()
add_subdirectory(boost)
