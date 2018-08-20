#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2018, by the GROMACS development team, led by
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

# For an executable target already defined, CMake tests are added with the
# following syntax:
#
# add_test(NAME <name> [CONFIGURATIONS [Debug|Release|...]]
#            [WORKING_DIRECTORY dir]
#            COMMAND <command> [arg1 [arg2 ...]])
#
# The Gromacs convention to get the appropriate compilation environment for
# unit tests is the gmx_add_unit_test macro, defined in
# src/testutils/TestMacros.cmake, which does not allow arbitrary argument
# passing to add_test. Instead, use set_tests_properties().

# Copy test files to test directory
configure_file(data/topol.tpr topol.tpr COPYONLY)
configure_file(data/water.tpr water.tpr COPYONLY)

configure_file(testingconfiguration.in.h testingconfiguration.h)

#
# Test public interface.
#

gmx_add_gtest_executable(
    gmxapi-test
    test_restraint.cpp
    test_runner.cpp
    test_status.cpp
    test_system.cpp
    test_version_helpers.cpp
    )

target_include_directories(gmxapi-test PRIVATE
                           ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(gmxapi-test Gromacs::gmxapi)

# Unit tests for implementation details may require access to private headers.
#target_link_libraries(gmxapi-impl-test Gromacs::gmxapi Gromacs::gmxapiImpl)

gmx_register_gtest_test(GmxapiExternalInterfaceTests gmxapi-test INTEGRATION_TEST)

set_tests_properties(GmxapiExternalInterfaceTests PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

#
# Test public interface with MPI.
#

gmx_add_gtest_executable(
    gmxapi-mpi-test MPI
    test_restraint.cpp
    test_runner.cpp
    test_status.cpp
    test_system.cpp
    test_version_helpers.cpp
    )

target_include_directories(gmxapi-mpi-test PRIVATE
                           ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(gmxapi-mpi-test Gromacs::gmxapi)

# Tests that specify MPI_RANKS will be passed a ``-ntmpi`` command-line option
# when GMX_THREAD_MPI is true, but the test executables in this directory do not
# yet understand that command-line option.
if(NOT GMX_THREAD_MPI)
    gmx_register_gtest_test(GmxapiMpiTests gmxapi-mpi-test MPI_RANKS 2 INTEGRATION_TEST)

    if(BUILD_TESTING)
        # Apparently gmx_register_gtest_test has some conditional like the following
        # when MPI_RANKS is set or when gmx_add_gtest_executable has MPI set.
        if(GMX_MPI OR GMX_THREAD_MPI)
            set_tests_properties(GmxapiMpiTests PROPERTIES
                                 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
        endif(GMX_MPI OR GMX_THREAD_MPI)
    endif()
endif()

#
# Unit testing for implementation details
#
gmx_add_gtest_executable(
    gmxapi-impl-test
    test_workflow.cpp
)

target_include_directories(gmxapi-impl-test PRIVATE
                           ${CMAKE_CURRENT_BINARY_DIR})

# Unit tests for implementation details require access to private headers.
target_link_libraries(gmxapi-impl-test Gromacs::gmxapi Gromacs::gmxapiImpl)

gmx_register_gtest_test(GmxapiImplTests gmxapi-impl-test INTEGRATION_TEST)

set_tests_properties(GmxapiImplTests PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

#
# Unit testing for implementation details with (t)MPI parallelism
#

gmx_add_gtest_executable(
    gmxapi-mpi-impl-test MPI
    test_workflow.cpp
)

target_include_directories(gmxapi-mpi-impl-test PRIVATE
                           ${CMAKE_CURRENT_BINARY_DIR})

# Unit tests for implementation details require access to private headers.
target_link_libraries(gmxapi-mpi-impl-test Gromacs::gmxapi Gromacs::gmxapiImpl)

# Tests that specify MPI_RANKS will be passed a ``-ntmpi`` command-line option
# when GMX_THREAD_MPI is true, but the test executables in this directory do not
# yet understand that command-line option.
if(NOT GMX_THREAD_MPI)
    gmx_register_gtest_test(GmxapiMpiImplTests gmxapi-mpi-impl-test MPI_RANKS 2 INTEGRATION_TEST)

    if(BUILD_TESTING)
        # Apparently gmx_register_gtest_test has some conditional like the following
        # when MPI_RANKS is set or when gmx_add_gtest_executable has MPI set.
        if(GMX_MPI OR GMX_THREAD_MPI)
            set_tests_properties(GmxapiMpiImplTests PROPERTIES
                                 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
        endif(GMX_MPI OR GMX_THREAD_MPI)
    endif()
endif()

# If we want to override the API version in the version.h for testing purposes, we'll do it another way...
#target_compile_definitions(gmxapi-test PUBLIC
#    GMXAPI_MAJOR=${GMXAPI_MAJOR}
#    GMXAPI_MINOR=${GMXAPI_MINOR}
#    GMXAPI_PATCH=${GMXAPI_PATCH})
