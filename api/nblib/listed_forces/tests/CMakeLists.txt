#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2020- The GROMACS Authors
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

# \author Victor Holanda <victor.holanda@cscs.ch>
# \author Joe Jordan <ejjordan@kth.se>
# \author Prashanth Kanduri <kanduri@cscs.ch>
# \author Sebastian Keller <keller@cscs.ch>

# Make a static library for test infrastructure code that we re-use
# in multiple test executables across the repository.

set(testname "NbLibListedForcesTests")
set(exename "nblib-listed-forces-test")

gmx_add_gtest_executable(
        ${exename}
        CPP_SOURCE_FILES
        # files with code for tests
        aggregate_transformations.cpp
        bondtypes.cpp
        free_energy.cpp
        gmxcalculator.cpp
        helpers.cpp
        listedtesthelpers.cpp
        kernels.cpp
        positionrestraints.cpp
        posresdx.cpp
        typetests.cpp
        calculator.cpp
        conversions.cpp
        shiftforces.cpp
        transformations.cpp
)
target_link_libraries(${exename} PRIVATE nblib_test_infrastructure nblib
                      listed_forces
                      mdtypes
                      pbcutil
                      timing)
gmx_register_gtest_test(${testname} ${exename} INTEGRATION_TEST)
add_dependencies(check-nblib ${exename})

unset(testname)
unset(exename)

# separate compile target, makes it easy to choose compiler optimization and flags
# and avoid dependencies on testing frameworks
add_executable(listedperformance ../../topologyhelpers.cpp gmxcalculator.cpp performance.cpp)
target_link_libraries(listedperformance PRIVATE nblib
                                                common
                                                libgromacs
                                                math
                                                pbcutil
                                                timing
                                                topology
                                                utility)
target_include_directories(listedperformance SYSTEM PRIVATE ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/src/external/thread_mpi/include)

if(CMAKE_CUDA_COMPILER AND GMX_GPU_CUDA AND NOT ${CMAKE_VERSION} VERSION_LESS "3.17.0" AND NOT ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "10.0.0")

    set(exename "nblib-listed-gpu-test")
    add_executable(${exename} listedgpu.cu)

    target_link_libraries(${exename} PRIVATE nblib
                                             common
                                             libgromacs
                                             math
                                             pbcutil
                                             timing
                                             topology
                                             utility
                                             gtest_main)
    target_include_directories(${exename} PRIVATE
            ${PROJECT_SOURCE_DIR}
            ${PROJECT_SOURCE_DIR}/api/legacy/include
            ${PROJECT_SOURCE_DIR}/src
            ${PROJECT_SOURCE_DIR}/src/external/thread_mpi/include)

    set_property(TARGET ${exename} PROPERTY CUDA_ARCHITECTURES "35;60")
    set_property(TARGET ${exename} PROPERTY CUDA_STANDARD 17)
    target_compile_options(${exename} PRIVATE -Wno-deprecated-gpu-targets)

    unset(exename)

    set(exename gpu-listedperformance)
    add_executable(${exename} ../gpu_interface.cu gpu_performance.cu)

    target_link_libraries(${exename} PRIVATE nblib
                                             common
                                             libgromacs
                                             math
                                             pbcutil
                                             timing
                                             topology
                                             utility)
    target_include_directories(${exename} PRIVATE
            ${PROJECT_SOURCE_DIR}
            ${PROJECT_SOURCE_DIR}/api/legacy/include
            ${PROJECT_SOURCE_DIR}/src
            ${PROJECT_SOURCE_DIR}/src/external/thread_mpi/include)

    set_property(TARGET ${exename} PROPERTY CUDA_ARCHITECTURES "35;60")
    set_property(TARGET ${exename} PROPERTY CUDA_STANDARD 17)
    target_compile_options(${exename} PRIVATE -Wno-deprecated-gpu-targets)

    unset(exename)
endif()
