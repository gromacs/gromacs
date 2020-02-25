#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2020, by the GROMACS development team, led by
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
#
# \author Victor Holanda <victor.holanda@cscs.ch>
# \author Joe Jordan <ejjordan@kth.se>
# \author Prashanth Kanduri <kanduri@cscs.ch>
# \author Sebastian Keller <keller@cscs.ch>
#

# Make a static library for test infrastructure code that we re-use
# in multiple test executables across the repository.
gmx_add_unit_test_library(nblib_test_infrastructure
  testhelpers.cpp
  testsystems.cpp
    )
target_include_directories(nblib_test_infrastructure SYSTEM PRIVATE ${PROJECT_SOURCE_DIR}/src/external)

set(testname "NbLibSetupTests")
set(exename "nblib-setup-test")

gmx_add_gtest_executable(
    ${exename}
    # files with code for tests
    particletype.cpp
    molecules.cpp
    box.cpp
    util.cpp
    topology.cpp
    # pseudo-library for code for mdrun
    #$<TARGET_OBJECTS:nblib_objlib>
    )
target_link_libraries(${exename} PRIVATE nblib_test_infrastructure)
gmx_register_gtest_test(${testname} ${exename} OPENMP_THREADS 1 INTEGRATION_TEST)


set(testname "NbLibTests")
set(exename "nblib-test")

gmx_add_gtest_executable(
    ${exename}
    # files with code for tests
    simstate.cpp
    nbkernelsystem.cpp
    nbnxnsetup.cpp
    # pseudo-library for code for mdrun
    #$<TARGET_OBJECTS:nblib_objlib>
    )
target_link_libraries(${exename} PRIVATE nblib_test_infrastructure)
gmx_register_gtest_test(${testname} ${exename} OPENMP_THREADS 1 INTEGRATION_TEST)
