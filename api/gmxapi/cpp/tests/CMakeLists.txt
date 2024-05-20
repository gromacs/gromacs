#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2018- The GROMACS Authors
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

# For an executable target already defined, CMake tests are added with the
# following syntax:
#
# add_test(NAME <name> [CONFIGURATIONS [Debug|Release|...]]
#            [WORKING_DIRECTORY dir]
#            COMMAND <command> [arg1 [arg2 ...]])
#
# The GROMACS convention to get the appropriate compilation environment for
# unit tests is the gmx_add_unit_test macro, defined in
# src/testutils/TestMacros.cmake, which does not allow arbitrary argument
# passing to add_test. Instead, use set_tests_properties().

#
# Test public interface.
#

add_library(gmxapi-testsupport INTERFACE)
target_include_directories(gmxapi-testsupport INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

# TODO: Test tMPI build in a MPI-enabled client context.
gmx_add_gtest_executable(gmxapi-test
    CPP_SOURCE_FILES
        restraint.cpp
        runner.cpp
        status.cpp
        stopsignaler.cpp
        system.cpp
        version.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
        )

# Link against the gmxapi libraries and get access to its public (installed) headers.
target_link_libraries(gmxapi-test PRIVATE Gromacs::gmxapi gmxapi-testsupport mdrun_test_infrastructure)

gmx_register_gtest_test(GmxapiExternalInterfaceTests gmxapi-test OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS QUICK_GPU_TEST)

if (GMX_CAN_RUN_MPI_TESTS)
    set_tests_properties(GmxapiExternalInterfaceTests PROPERTIES
                         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()

#
# Test public interface with MPI.
#
# gmxapi is expected to work whether or not GROMACS was built with MPI (or tMPI)
# but we don't assume that MPI is available to the build tree or testing tools
# unless CMake was configured with MPI.
#
# For MPI-enabled testing of gmxapi clients using non-MPI-enabled GROMACS, we
# defer testing to the Python gmxapi client package tests.
#
if (GMX_MPI)
    gmx_add_gtest_executable(gmxapi-mpi-test MPI
                             CPP_SOURCE_FILES
                             context.cpp
                             restraint.cpp
                             runner.cpp
                             status.cpp
                             stopsignaler.cpp
                             system.cpp
                             version.cpp
                             # pseudo-library for code for mdrun
                             $<TARGET_OBJECTS:mdrun_objlib>
                             )

    target_include_directories(gmxapi-mpi-test PRIVATE
                               ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(gmxapi-mpi-test PRIVATE Gromacs::gmxapi gmxapi-testsupport mdrun_test_infrastructure)

    gmx_register_gtest_test(GmxapiMpiTests gmxapi-mpi-test MPI_RANKS 2 OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS)

    if (TEST GmxapiMpiTests AND GMX_CAN_RUN_MPI_TESTS)
        set_tests_properties(GmxapiMpiTests PROPERTIES
                             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    endif ()
    if (TEST GmxapiExternalInterfaceTests AND TEST GmxapiMpiTests)
        # Because the gmxapi tests depend on writing files whose names are not
        # unique across test cases, prevent CTest from running these tests
        # concurrently. See #4654
        set_tests_properties(GmxapiMpiTests PROPERTIES DEPENDS GmxapiExternalInterfaceTests)
    endif()
endif ()
