#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2013,2014,2015,2016,2017 The GROMACS development team.
# Copyright (c) 2018,2019,2020,2021, by the GROMACS development team, led by
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

# Make a static library for test infrastructure code that we re-use
# in multiple test executables across the repository.
gmx_add_unit_test_library(mdrun_test_infrastructure
    energyreader.cpp
    energycomparison.cpp
    moduletest.cpp
    simulatorcomparison.cpp
    terminationhelper.cpp
    trajectorycomparison.cpp
    trajectoryreader.cpp
    # pseudo-library for code for mdrun
    $<TARGET_OBJECTS:mdrun_objlib>
    )
target_include_directories(mdrun_test_infrastructure SYSTEM PRIVATE ${PROJECT_SOURCE_DIR}/src/external)
target_link_libraries(mdrun_test_infrastructure PUBLIC legacy_api utility)

# To avoid running into test timeouts, some end-to-end tests of mdrun
# functionality are split off. This can be rearranged in future as we
# see fit.
set(testname "MdrunOutputTests")
set(exename "mdrun-output-test")

gmx_add_gtest_executable(${exename}
    CPP_SOURCE_FILES
        compressed_x_output.cpp
        helpwriting.cpp
        outputfiles.cpp
        trajectory_writing.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
    )
target_link_libraries(${exename} PRIVATE mdrun_test_infrastructure)
gmx_register_gtest_test(${testname} ${exename} OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS)

set(testname "MdrunModulesTests")
set(exename "mdrun-modules-test")

gmx_add_gtest_executable(${exename}
    CPP_SOURCE_FILES
        densityfittingmodule.cpp
        interactiveMD.cpp
        mimic.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
    )
target_link_libraries(${exename} PRIVATE mdrun_test_infrastructure)
gmx_register_gtest_test(${testname} ${exename} OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS)

set(testname "MdrunIOTests")
set(exename "mdrun-io-test")

gmx_add_gtest_executable(${exename}
    CPP_SOURCE_FILES
        checkpoint.cpp
        exactcontinuation.cpp
        grompp.cpp
        initialconstraints.cpp
        termination.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
    )
target_link_libraries(${exename} PRIVATE mdrun_test_infrastructure)
gmx_register_gtest_test(${testname} ${exename} OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS)

# To avoid running into test timeouts, some end-to-end tests of mdrun
# functionality are split off. This can be rearranged in future as we
# see fit.
set(testname "MdrunTests")
set(exename "mdrun-test")

gmx_add_gtest_executable(${exename}
    CPP_SOURCE_FILES
        ewaldsurfaceterm.cpp
        multiple_time_stepping.cpp
        orires.cpp
        swapcoords.cpp
        tabulated_bonded_interactions.cpp
        freezegroups.cpp
	constantacceleration.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
    )
target_link_libraries(${exename} PRIVATE mdrun_test_infrastructure)
gmx_register_gtest_test(${testname} ${exename} OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS)


set(testname "MdrunPmeTests")
set(exename "mdrun-pme-test")

gmx_add_gtest_executable(${exename}
    CPP_SOURCE_FILES
        pmetest.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
    )
target_link_libraries(${exename} PRIVATE mdrun_test_infrastructure)
gmx_register_gtest_test(${testname} ${exename} OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS)


set(testname "MdrunNonIntegratorTests")
set(exename "mdrun-non-integrator-test")

gmx_add_gtest_executable(${exename}
    CPP_SOURCE_FILES
        # files with code for tests
        minimize.cpp
        nonbonded_bench.cpp
        normalmodes.cpp
        rerun.cpp
        simple_mdrun.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
        )
target_link_libraries(${exename} PRIVATE mdrun_test_infrastructure)
gmx_register_gtest_test(${testname} ${exename} OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS)

# TPI does not support OpenMP, so we need a separate test binary
set(testname "MdrunTpiTests")
set(exename "mdrun-tpi-test")

gmx_add_gtest_executable(${exename}
    CPP_SOURCE_FILES
        # files with code for tests
        tpitest.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
        )
target_link_libraries(${exename} PRIVATE mdrun_test_infrastructure)
gmx_register_gtest_test(${testname} ${exename} INTEGRATION_TEST IGNORE_LEAKS)

# Tests that only make sense to run with multiple ranks and/or real
# MPI are implemented here.
set(testname "MdrunMpiTests")
set(exename "mdrun-mpi-test")

gmx_add_gtest_executable(${exename} MPI
    CPP_SOURCE_FILES
        # files with code for tests
        domain_decomposition.cpp
        minimize.cpp
        mimic.cpp
        pmetest.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
        )
target_link_libraries(${exename} PRIVATE mdrun_test_infrastructure)
gmx_register_gtest_test(${testname} ${exename} MPI_RANKS 2 OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS)

# Multi sim only makes sense with real MPI, and ideally at least 4 ranks,
# to allow for multiple simulations (>= 2 sims) each using DD (>= 2 ranks per sim)
set(testname "MdrunMultiSimTests")
set(exename "mdrun-multisim-test")

gmx_add_gtest_executable(${exename} MPI
    CPP_SOURCE_FILES
        # files with code for tests
        multisim.cpp
        multisimtest.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
        )
target_link_libraries(${exename} PRIVATE mdrun_test_infrastructure)
gmx_register_gtest_test(${testname} ${exename} MPI_RANKS 4 INTEGRATION_TEST IGNORE_LEAKS)

set(testname "MdrunMultiSimReplexTests")
set(exename "mdrun-multisim-replex-test")

gmx_add_gtest_executable(${exename} MPI
    CPP_SOURCE_FILES
        # files with code for tests
	multisimtest.cpp
        replicaexchange.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
        )
target_link_libraries(${exename} PRIVATE mdrun_test_infrastructure)
gmx_register_gtest_test(${testname} ${exename} MPI_RANKS 4 INTEGRATION_TEST IGNORE_LEAKS)

set(testname "MdrunMultiSimReplexEquivalenceTests")
set(exename "mdrun-multisim-replex-equivalence-test")

gmx_add_gtest_executable(${exename} MPI
    CPP_SOURCE_FILES
        # files with code for tests
	multisimtest.cpp
        replicaexchange_equivalence.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
        )
target_link_libraries(${exename} PRIVATE mdrun_test_infrastructure)
gmx_register_gtest_test(${testname} ${exename} MPI_RANKS 4 INTEGRATION_TEST IGNORE_LEAKS)

# Tests that only make sense to run with multiple ranks and/or real
# MPI are implemented here. Special case for slow PME tests
set(testname "MdrunMpiPmeTests")
set(exename "mdrun-mpi-pme-test")

gmx_add_gtest_executable(${exename} MPI
    CPP_SOURCE_FILES
        # files with code for tests
        pmetest.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
        )
target_link_libraries(${exename} PRIVATE mdrun_test_infrastructure)
gmx_register_gtest_test(${testname} ${exename} MPI_RANKS 2 OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS)

# Slow-running tests that target testing multiple-rank coordination behaviors
# These tests are extremely slow without optimization or OpenMP, so only run them for
# build types like Release or RelWithDebInfo and if the build has been configured
# with OpenMP enabled
set(exename "mdrun-mpi-coordination-basic-test")
gmx_add_gtest_executable(${exename} MPI
    CPP_SOURCE_FILES
        # files with code for tests
        periodicactions.cpp
	periodicactions_basic.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
        )
target_link_libraries(${exename} PRIVATE mdrun_test_infrastructure)

if (CMAKE_BUILD_TYPE MATCHES "Rel" AND GMX_OPENMP)
    set(testname "MdrunMpiCoordinationBasicTestsOneRank")
    gmx_register_gtest_test(${testname} ${exename} MPI_RANKS 1 SLOW_TEST IGNORE_LEAKS)
    set(testname "MdrunMpiCoordinationBasicTestsTwoRanks")
    gmx_register_gtest_test(${testname} ${exename} MPI_RANKS 2 SLOW_TEST IGNORE_LEAKS)
endif()

set(exename "mdrun-mpi-coordination-coupling-test")
gmx_add_gtest_executable(${exename} MPI
    CPP_SOURCE_FILES
        # files with code for tests
        periodicactions.cpp
	periodicactions_coupling.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
        )
target_link_libraries(${exename} PRIVATE mdrun_test_infrastructure)

if (CMAKE_BUILD_TYPE MATCHES "Rel" AND GMX_OPENMP)
    set(testname "MdrunMpiCoordinationCouplingTestsOneRank")
    gmx_register_gtest_test(${testname} ${exename} MPI_RANKS 1 SLOW_TEST IGNORE_LEAKS)
    set(testname "MdrunMpiCoordinationCouplingTestsTwoRanks")
    gmx_register_gtest_test(${testname} ${exename} MPI_RANKS 2 SLOW_TEST IGNORE_LEAKS)
endif()

set(exename "mdrun-mpi-coordination-constraints-test")
gmx_add_gtest_executable(${exename} MPI
    CPP_SOURCE_FILES
        # files with code for tests
        periodicactions.cpp
	periodicactions_constraints.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
        )
target_link_libraries(${exename} PRIVATE mdrun_test_infrastructure)

if (CMAKE_BUILD_TYPE MATCHES "Rel" AND GMX_OPENMP)
    set(testname "MdrunMpiCoordinationConstraintsTestsOneRank")
    gmx_register_gtest_test(${testname} ${exename} MPI_RANKS 1 SLOW_TEST IGNORE_LEAKS)
    set(testname "MdrunMpiCoordinationConstraintsTestsTwoRanks")
    gmx_register_gtest_test(${testname} ${exename} MPI_RANKS 2 SLOW_TEST IGNORE_LEAKS)
endif()

# Keeping the FEP tests separate for now to be able to judge runtime more easily
# Can be included in mdrun tests later
set(testname "MdrunFEPTests")
set(exename "mdrun-fep-test")

gmx_add_gtest_executable(${exename}
        CPP_SOURCE_FILES
        # files with code for tests
        freeenergy.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
)
target_link_libraries(${exename} PRIVATE mdrun_test_infrastructure)
# TODO: Link specific modules: topology

# End-to-end test for pull functionality
target_link_libraries(${exename} PRIVATE legacy_modules)
gmx_register_gtest_test(${testname} ${exename} OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS)
set(testname "MdrunPullTests")
set(exename "mdrun-pull-test")

gmx_add_gtest_executable(${exename}
        CPP_SOURCE_FILES
        # files with code for tests
        pull.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
)
target_link_libraries(${exename} PRIVATE mdrun_test_infrastructure)
target_link_libraries(${exename} PRIVATE legacy_modules)
gmx_register_gtest_test(${testname} ${exename} OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS)

# End-to-end tests comparing different simulator code paths
set(testname "MdrunSimulatorComparison")
set(exename "mdrun-simulator-comparison-test")

gmx_add_gtest_executable(${exename}
        CPP_SOURCE_FILES
        # files with code for tests
        simulator.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
        )
target_link_libraries(${exename} PRIVATE mdrun_test_infrastructure)
gmx_register_gtest_test(${testname} ${exename} OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS)

# Tests checking virtual sites - we want to run them with multiple ranks
set(testname "MdrunVirtualSiteTests")
set(exename "mdrun-vsites-test")

gmx_add_gtest_executable(${exename} MPI
        CPP_SOURCE_FILES
        # files with code for tests
        virtualsites.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
        )
target_link_libraries(${exename} PRIVATE mdrun_test_infrastructure)
gmx_register_gtest_test(${testname} ${exename} MPI_RANKS 2 INTEGRATION_TEST IGNORE_LEAKS)
