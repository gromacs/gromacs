#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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

# Make an "object library" for test infrastructure code that we re-use
# in multiple test executables.
add_library(mdrun_test_objlib OBJECT
    energyreader.cpp
    energycomparison.cpp
    moduletest.cpp
    simulationdatabase.cpp
    terminationhelper.cpp
    trajectorycomparison.cpp
    trajectoryreader.cpp
    )

set(testname "MdrunTests")
set(exename "mdrun-test")

gmx_add_gtest_executable(
    ${exename}
    # files with code for tests
    compressed_x_output.cpp
    grompp.cpp
    initialconstraints.cpp
    interactiveMD.cpp
    pmetest.cpp
    swapcoords.cpp
    tabulated_bonded_interactions.cpp
    termination.cpp
    trajectory_writing.cpp
    mimic.cpp
    # pseudo-library for code for testing mdrun
    $<TARGET_OBJECTS:mdrun_test_objlib>
    # pseudo-library for code for mdrun
    $<TARGET_OBJECTS:mdrun_objlib>
    )
gmx_register_gtest_test(${testname} ${exename} OPENMP_THREADS 2 INTEGRATION_TEST)

# To avoid running into test timeouts, some end-to-end tests of mdrun
# functionality are split off. This can be rearranged in future as we
# see fit.
set(testname "MdrunNonIntegratorTests")
set(exename "mdrun-non-integrator-test")

gmx_add_gtest_executable(
    ${exename}
    # files with code for tests
    minimize.cpp
    rerun.cpp
    # pseudo-library for code for testing mdrun
    $<TARGET_OBJECTS:mdrun_test_objlib>
    # pseudo-library for code for mdrun
    $<TARGET_OBJECTS:mdrun_objlib>
    )
gmx_register_gtest_test(${testname} ${exename} OPENMP_THREADS 2 INTEGRATION_TEST)

# To permit other end-to-end tests to run with multiple OpenMP
# threads, test cases that still use the group scheme are separated.
set(testname "LegacyGroupSchemeMdrunTests")
set(exename "legacy-mdrun-test")

gmx_add_gtest_executable(
    ${exename}
    # When TPI supports the Verlet scheme, move this into
    # MdrunNonIntegratorTest
    tpitest.cpp
    # pseudo-library for code for testing mdrun
    $<TARGET_OBJECTS:mdrun_test_objlib>
    # pseudo-library for code for mdrun
    $<TARGET_OBJECTS:mdrun_objlib>
    )
gmx_register_gtest_test(${testname} ${exename} INTEGRATION_TEST)

# Tests that only make sense to run with multiple ranks and/or real
# MPI are implemented here.
set(testname "MdrunMpiTests")
set(exename "mdrun-mpi-test")

gmx_add_gtest_executable(
    ${exename} MPI
    # files with code for tests
    domain_decomposition.cpp
    minimize.cpp
    mimic.cpp
    multisim.cpp
    multisimtest.cpp
    pmetest.cpp
    replicaexchange.cpp
    # pseudo-library for code for testing mdrun
    $<TARGET_OBJECTS:mdrun_test_objlib>
    # pseudo-library for code for mdrun
    $<TARGET_OBJECTS:mdrun_objlib>
    )
gmx_register_gtest_test(${testname} ${exename} OPENMP_THREADS 2 INTEGRATION_TEST)
