#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2013- The GROMACS Authors
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

# Make a static library for test infrastructure code that we re-use
# in multiple test executables across the repository.
gmx_add_unit_test_library(mdrun_test_infrastructure
    energyreader.cpp
    energycomparison.cpp
    moduletest.cpp
    simulatorcomparison.cpp
    terminationhelper.cpp
    trajectorycomparison.cpp
    # pseudo-library for code for mdrun
    $<TARGET_OBJECTS:mdrun_objlib>
    )
target_include_directories(mdrun_test_infrastructure SYSTEM PRIVATE ${PROJECT_SOURCE_DIR}/src/external)
target_link_libraries(mdrun_test_infrastructure PRIVATE legacy_api utility)

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
gmx_register_gtest_test(${testname} ${exename} OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS QUICK_GPU_TEST)

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
gmx_register_gtest_test(${testname} ${exename} OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS QUICK_GPU_TEST)

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
gmx_register_gtest_test(${testname} ${exename} OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS SLOW_GPU_TEST)

# To avoid running into test timeouts, some end-to-end tests of mdrun
# functionality are split off. This can be rearranged in future as we
# see fit.
set(exename "mdrun-test")

gmx_add_gtest_executable(${exename} MPI
    CPP_SOURCE_FILES
        multiple_time_stepping.cpp
        swapcoords.cpp
        tabulated_bonded_interactions.cpp
        freezegroups.cpp
        constantacceleration.cpp
        boxdeformation.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
    )
target_link_libraries(${exename} PRIVATE mdrun_test_infrastructure)
gmx_register_gtest_test(MdrunTestsOneRank ${exename} MPI_RANKS 1 OPENMP_THREADS 2 INTEGRATION_TEST SLOW_TEST IGNORE_LEAKS SLOW_GPU_TEST)
gmx_register_gtest_test(MdrunTestsTwoRanks ${exename} MPI_RANKS 2 OPENMP_THREADS 2 INTEGRATION_TEST SLOW_TEST IGNORE_LEAKS SLOW_GPU_TEST)
# Because the mdrun tests depend on writing files whose names are not
# unique across test cases, prevent CTest from running these tests
# concurrently. See #4654
if (TEST MdrunTestsTwoRanks AND TEST MdrunTestsOneRank)
    set_tests_properties(MdrunTestsTwoRanks PROPERTIES DEPENDS MdrunTestsOneRank)
endif()

# The orires test is separate, as it supports only a single MPI rank
set(testname "MdrunSingleRankAlgorithmsTests")
set(exename "mdrun-single-rank-algorithms-test")

gmx_add_gtest_executable(${exename}
    CPP_SOURCE_FILES
        dispersion_correction.cpp
        ewaldsurfaceterm.cpp
        orires.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
    )
target_link_libraries(${exename} PRIVATE mdrun_test_infrastructure)
gmx_register_gtest_test(${testname} ${exename} OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS SLOW_GPU_TEST)

set(exename "minimize-test")

gmx_add_gtest_executable(${exename} MPI
    CPP_SOURCE_FILES
        # files with code for tests
        minimize.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
        )
target_link_libraries(${exename} PRIVATE mdrun_test_infrastructure)
gmx_register_gtest_test(Minimize1RankTests ${exename} MPI_RANKS 1 OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS QUICK_GPU_TEST)
gmx_register_gtest_test(Minimize2RankTests ${exename} MPI_RANKS 2 OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS QUICK_GPU_TEST)
# Because the minimizer tests depend on writing files whose names are not
# unique across test cases, prevent CTest from running these tests
# concurrently. See #4654
if (TEST Minimize2RankTests AND TEST Minimize1RankTests)
    set_tests_properties(Minimize2RankTests PROPERTIES DEPENDS Minimize1RankTests)
endif()

set(testname "MdrunNonIntegratorTests")
set(exename "mdrun-non-integrator-test")

gmx_add_gtest_executable(${exename}
    CPP_SOURCE_FILES
        # files with code for tests
        nonbonded_bench.cpp
        normalmodes.cpp
        rerun.cpp
        simple_mdrun.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
        )
target_link_libraries(${exename} PRIVATE mdrun_test_infrastructure)
gmx_register_gtest_test(${testname} ${exename} OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS SLOW_GPU_TEST)

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
gmx_register_gtest_test(${testname} ${exename} INTEGRATION_TEST IGNORE_LEAKS QUICK_GPU_TEST)

# Tests that only make sense to run with multiple ranks and/or real
# MPI are implemented here.
set(testname "MdrunMpiTests")
set(exename "mdrun-mpi-test")

gmx_add_gtest_executable(${exename} MPI
    CPP_SOURCE_FILES
        # files with code for tests
        domain_decomposition.cpp
        mimic.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
        )
target_link_libraries(${exename} PRIVATE mdrun_test_infrastructure)
gmx_register_gtest_test(${testname} ${exename} MPI_RANKS 2 OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS QUICK_GPU_TEST)

# Because the mdrun tests depend on writing files whose names are not
# unique across test cases, prevent CTest from running these tests
# concurrently. See #4654. These test binaries both test aspects of MiMiC.
if (TEST MdrunMpiTests AND TEST MdrunModulesTests)
    set_tests_properties(MdrunMpiTests PROPERTIES DEPENDS MdrunModulesTests)
endif()

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
gmx_register_gtest_test(${testname} ${exename} MPI_RANKS 4 INTEGRATION_TEST IGNORE_LEAKS QUICK_GPU_TEST)

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
gmx_register_gtest_test(${testname} ${exename} MPI_RANKS 4 INTEGRATION_TEST IGNORE_LEAKS QUICK_GPU_TEST)

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
gmx_register_gtest_test(${testname} ${exename} MPI_RANKS 4 INTEGRATION_TEST IGNORE_LEAKS SLOW_GPU_TEST)

# Tests that only make sense to run with multiple ranks and/or real
# MPI are implemented here. Special case for slow PME tests

set(exename "mdrun-mpi-pme-test")
gmx_add_gtest_executable(${exename} MPI
    CPP_SOURCE_FILES
        # files with code for tests
        pmetest.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
        )
target_link_libraries(${exename} PRIVATE mdrun_test_infrastructure)
gmx_register_gtest_test(MdrunMpi1RankPmeTests ${exename} MPI_RANKS 1 OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS SLOW_GPU_TEST)
gmx_register_gtest_test(MdrunMpi2RankPmeTests ${exename} MPI_RANKS 2 OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS SLOW_GPU_TEST)
# Because the mdrun tests depend on writing files whose names are not
# unique across test cases, prevent CTest from running these tests
# concurrently. See #4654
if (TEST MdrunMpi2RankPmeTests AND TEST MdrunMpi1RankPmeTests)
    set_tests_properties(MdrunMpi2RankPmeTests PROPERTIES DEPENDS MdrunMpi1RankPmeTests)
endif()

# Slow-running tests that target testing multiple-rank coordination behaviors
# These tests are extremely slow without optimization or OpenMP, so only run them for
# build types like Release or RelWithDebInfo and if the build has been configured
# with OpenMP enabled
set(exename "mdrun-coordination-basic-test")
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
    set(testname "MdrunCoordinationBasicTests1Rank")
    gmx_register_gtest_test(${testname} ${exename} MPI_RANKS 1 SLOW_TEST IGNORE_LEAKS QUICK_GPU_TEST)
    set(testname "MdrunCoordinationBasicTests2Ranks")
    gmx_register_gtest_test(${testname} ${exename} MPI_RANKS 2 SLOW_TEST IGNORE_LEAKS QUICK_GPU_TEST)
    # Because the mdrun tests depend on writing files whose names are not
    # unique across test cases, prevent CTest from running these tests
    # concurrently. See #4654
    if (TEST MdrunCoordinationBasicTests2Ranks AND TEST MdrunCoordinationBasicTests1Rank)
        set_tests_properties(MdrunCoordinationBasicTests2Ranks PROPERTIES DEPENDS MdrunCoordinationBasicTests1Rank)
    endif()
endif()

set(exename "mdrun-coordination-coupling-test")
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
    set(testname "MdrunCoordinationCouplingTests1Rank")
    gmx_register_gtest_test(${testname} ${exename} MPI_RANKS 1 SLOW_TEST IGNORE_LEAKS SLOW_GPU_TEST)
    set(testname "MdrunCoordinationCouplingTests2Ranks")
    gmx_register_gtest_test(${testname} ${exename} MPI_RANKS 2 SLOW_TEST IGNORE_LEAKS SLOW_GPU_TEST)
    # Because the mdrun tests depend on writing files whose names are not
    # unique across test cases, prevent CTest from running these tests
    # concurrently. See #4654
    if (TEST MdrunCoordinationCouplingTests2Ranks AND TEST MdrunCoordinationCouplingTests1Rank)
        set_tests_properties(MdrunCoordinationCouplingTests2Ranks PROPERTIES DEPENDS MdrunCoordinationCouplingTests1Rank)
    endif()
endif()

set(exename "mdrun-coordination-constraints-test")
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
    set(testname "MdrunCoordinationConstraintsTests1Rank")
    gmx_register_gtest_test(${testname} ${exename} MPI_RANKS 1 SLOW_TEST IGNORE_LEAKS SLOW_GPU_TEST)
    set(testname "MdrunCoordinationConstraintsTests2Ranks")
    gmx_register_gtest_test(${testname} ${exename} MPI_RANKS 2 SLOW_TEST IGNORE_LEAKS SLOW_GPU_TEST)
    # Because the mdrun tests depend on writing files whose names are not
    # unique across test cases, prevent CTest from running these tests
    # concurrently. See #4654
    if (TEST MdrunCoordinationConstraintsTests2Ranks AND TEST MdrunCoordinationConstraintsTests1Rank)
        set_tests_properties(MdrunCoordinationConstraintsTests2Ranks PROPERTIES DEPENDS MdrunCoordinationConstraintsTests1Rank)
    endif()
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
gmx_register_gtest_test(${testname} ${exename} OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS SLOW_GPU_TEST)

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
gmx_register_gtest_test(${testname} ${exename} OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS QUICK_GPU_TEST)

# End-to-end test for enforced rotation ("rotational pulling") functionality
set(testname "MdrunRotationTests")
set(exename "mdrun-rotation-test")

gmx_add_gtest_executable(${exename}
        CPP_SOURCE_FILES
        # files with code for tests
        pull_rotation.cpp
        # pseudo-library for code for mdrun
        $<TARGET_OBJECTS:mdrun_objlib>
)
target_include_directories(${exename} PRIVATE ${PROJECT_SOURCE_DIR}/src/gromacs/math/include)
target_link_libraries(${exename} PRIVATE mdrun_test_infrastructure)
target_link_libraries(${exename} PRIVATE legacy_modules)
gmx_register_gtest_test(${testname} ${exename} OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS QUICK_GPU_TEST)

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
gmx_register_gtest_test(${testname} ${exename} OPENMP_THREADS 2 INTEGRATION_TEST IGNORE_LEAKS SLOW_GPU_TEST)

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
gmx_register_gtest_test(${testname} ${exename} MPI_RANKS 2 INTEGRATION_TEST IGNORE_LEAKS QUICK_GPU_TEST)
