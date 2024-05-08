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
#

add_custom_target(nblib-tests
        COMMENT "Target to build all nblib tests including samples"
        )
# Ensure that "make tests" builds all nblib tests so the top-level
# "make check" will work.
if (BUILD_TESTING)
        add_dependencies(tests nblib-tests)

        # this allows all nblib tests to be run with "make check-nblib"
        add_custom_target(check-nblib
                COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -R NbLib
                COMMENT "Running nblib tests"
                USES_TERMINAL VERBATIM)
        add_dependencies(check-nblib nblib-tests)
endif()

set(NBLIB_MAJOR 0)
set(NBLIB_MINOR 1)
set(NBLIB_RELEASE ${NBLIB_MAJOR}.${NBLIB_MINOR}.0)

add_library(nblib)
set_target_properties(nblib PROPERTIES
        VERSION_MAJOR ${NBLIB_MAJOR}
        VERSION_MINOR ${NBLIB_MINOR}
        SOVERSION ${NBLIB_MAJOR}
        RELEASE ${NBLIB_RELEASE}
        VERSION ${NBLIB_RELEASE}
        LINKER_LANGUAGE CXX
        OUTPUT_NAME "nblib_gmx")

if (GMX_GPU_CUDA)
    target_link_libraries(nblib PRIVATE CUDA::cudart_static)
    if(GMX_USE_cuFFTMp OR GMX_NVSHMEM)
        # We need to PUBLIC link to the stub libraries nvml/cuda to WAR an issue
        # with `bfd` linker which expects you always have dependencies on the link line.
        # Also leaking path via the -rpath does not work with a stub library as the SONAME
        # present inside these stub library (used by runtime) is *.so.1 instead of *.so as these
        # stub library are meant only for compile time and the runtime version is a driver provided
        # version which is not present on the system where driver is not installed.
        # At some future release cmake will provide a workaround for this under the hood at which
        # point we can use PRIVATE linking for these libs until then we stick with PUBLIC linking.
        target_link_libraries(nblib PUBLIC CUDA::nvml CUDA::cuda_driver)
    endif()
endif()

# Workaround for #4922
target_link_libraries(nblib PUBLIC ${GMX_PUBLIC_LIBRARIES})

target_sources(nblib
        PRIVATE
        box.cpp
        gmxcalculatorcpu.cpp
        integrator.cpp
        interactions.cpp
        molecules.cpp
        nbnxmsetuphelpers.cpp
        particlesequencer.cpp
        particletype.cpp
        simulationstate.cpp
        topologyhelpers.cpp
        topology.cpp
        tpr.cpp
        virials.cpp
        )

gmx_target_compile_options(nblib)

target_link_libraries(nblib PRIVATE libgromacs)
target_include_directories(nblib PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:include>)
target_link_libraries(nblib PRIVATE common)
target_link_libraries(nblib PRIVATE
                      fileio
                      gmxlib
                      listed_forces
                      math
                      mdtypes
                      mdlib
                      pbcutil
                      simd
                      topology
                      utility)

# There are transitive interface dependencies on the legacy GROMACS headers.
# TODO(#3288): Explicitly link specific modules for public API dependencies as legacy_api is removed.
target_link_libraries(nblib PUBLIC legacy_api)
# TODO: Explicitly link specific modules.
target_link_libraries(nblib PRIVATE legacy_modules)
if (GMX_OPENMP)
    target_link_libraries(nblib PUBLIC OpenMP::OpenMP_CXX)
endif()

install(TARGETS nblib
        EXPORT nblib
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        INCLUDES DESTINATION include
        )

if(GMX_INSTALL_NBLIB_API)
    install(DIRECTORY include/nblib DESTINATION include)
endif()

add_subdirectory(listed_forces)
add_subdirectory(samples)
add_subdirectory(util)

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
