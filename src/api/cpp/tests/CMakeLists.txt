#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
# The GROMACS convention to get the appropriate compilation environment for
# unit tests is the gmx_add_unit_test macro, defined in
# src/testutils/TestMacros.cmake, which does not allow arbitrary argument
# passing to add_test. Instead, use set_tests_properties().

#
# Test public interface.
#

gmx_add_gtest_executable(gmxapi-test
    CPP_SOURCE_FILES
        restraint.cpp
        status.cpp
        system.cpp
        version.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
        )
if(NOT GMX_USE_OPENCL)
    # GPU resources may not be properly reinitialized between simulations in
    # the same process.
    # TODO: include this with the other test sources once the issue is resolved
    # Ref https://gitlab.com/gromacs/gromacs/-/issues/2689
    target_sources(gmxapi-test PRIVATE runner.cpp stopsignaler.cpp)
endif()

target_include_directories(gmxapi-test PRIVATE
                           ${CMAKE_CURRENT_SOURCE_DIR})

# Link against the gmxapi libraries and get access to its public (installed) headers.
target_link_libraries(gmxapi-test PRIVATE Gromacs::gmxapi mdrun_test_infrastructure)

gmx_register_gtest_test(GmxapiExternalInterfaceTests gmxapi-test OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS)

set_tests_properties(GmxapiExternalInterfaceTests PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

#
# Test public interface with MPI.
#

gmx_add_gtest_executable(gmxapi-mpi-test MPI
    CPP_SOURCE_FILES
        restraint.cpp
        status.cpp
        system.cpp
        version.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
        )
if(NOT GMX_USE_OPENCL)
    # GPU resources may not be properly reinitialized between simulations in
    # the same process.
    # TODO: include this with the other test sources once the issue is resolved
    # Ref https://gitlab.com/gromacs/gromacs/-/issues/2689
    target_sources(gmxapi-mpi-test PRIVATE runner.cpp stopsignaler.cpp)
endif()

target_include_directories(gmxapi-mpi-test PRIVATE
                           ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(gmxapi-mpi-test PRIVATE Gromacs::gmxapi mdrun_test_infrastructure)

gmx_register_gtest_test(GmxapiMpiTests gmxapi-mpi-test MPI_RANKS 2 OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS)
