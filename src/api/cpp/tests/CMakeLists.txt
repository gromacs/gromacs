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
# The GROMACS convention to get the appropriate compilation environment for
# unit tests is the gmx_add_unit_test macro, defined in
# src/testutils/TestMacros.cmake, which does not allow arbitrary argument
# passing to add_test. Instead, use set_tests_properties().

# Those tests can only be run if we are not just building the mdrun binary.
if (NOT GMX_BUILD_MDRUN_ONLY)

# Build test input files to test directory.
# Note that CMake interprets `gmx` in `COMMAND gmx` as the executable produced
# by the `add_executable` command that defined the `gmx` target.
# \todo replace command line wrapper with API utilities or test framework facilities.
add_custom_command(
    OUTPUT topol.tpr
    COMMAND gmx grompp -o topol.tpr
        -f ${CMAKE_CURRENT_SOURCE_DIR}/data/spc-and-methanol.mdp
        -c ${CMAKE_CURRENT_SOURCE_DIR}/data/spc-and-methanol.gro
        -n ${CMAKE_CURRENT_SOURCE_DIR}/data/spc-and-methanol.ndx
        -p ${CMAKE_CURRENT_SOURCE_DIR}/data/spc-and-methanol.top
    DEPENDS gmx
            data/spc-and-methanol.gro
            data/spc-and-methanol.mdp
            data/spc-and-methanol.top
            data/spc-and-methanol.ndx
)

add_custom_target(spc_and_methanol_tpr
                  DEPENDS topol.tpr
                  COMMENT "Preparing test reference data."
                  )

configure_file(testingconfiguration.in.cpp testingconfiguration.cpp)

#
# Test public interface.
#

gmx_add_gtest_executable(
    gmxapi-test
    status.cpp
    testingconfiguration.cpp
    version.cpp
    )
add_dependencies(gmxapi-test spc_and_methanol_tpr)

target_include_directories(gmxapi-test PRIVATE
                           ${CMAKE_CURRENT_SOURCE_DIR})

# Link against the gmxapi libraries and get access to its public (installed) headers.
target_link_libraries(gmxapi-test Gromacs::gmxapi)

gmx_register_gtest_test(GmxapiExternalInterfaceTests gmxapi-test INTEGRATION_TEST)

set_tests_properties(GmxapiExternalInterfaceTests PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

#
# Test public interface with MPI.
#

gmx_add_gtest_executable(
    gmxapi-mpi-test MPI
    status.cpp
    testingconfiguration.cpp
    version.cpp
    )
add_dependencies(gmxapi-mpi-test spc_and_methanol_tpr)

target_include_directories(gmxapi-mpi-test PRIVATE
                           ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(gmxapi-mpi-test Gromacs::gmxapi)

# Tests that specify MPI_RANKS will be passed a ``-ntmpi`` command-line option
# when GMX_THREAD_MPI is true, but the test executables in this directory do not
# yet understand that command-line option.
if(NOT GMX_THREAD_MPI)
    gmx_register_gtest_test(GmxapiMpiTests gmxapi-mpi-test MPI_RANKS 2 INTEGRATION_TEST)
endif()

endif()
