#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2012- The GROMACS Authors
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

# Custom build type "Reference", to be used for creating new
# reference values in the GROMACS regression tests.
set( CMAKE_CXX_FLAGS_REFERENCE "-O0 -g" CACHE STRING "C++ flags for regressiontests reference runs." FORCE)
set( CMAKE_C_FLAGS_REFERENCE "-O0 -g" CACHE STRING "C flags for regressiontests reference runs." FORCE)
set( CMAKE_EXE_LINKER_FLAGS_REFERENCE "" CACHE STRING "Linker flags for regressiontests reference runs.")
mark_as_advanced( CMAKE_CXX_FLAGS_REFERENCE CMAKE_C_FLAGS_REFERENCE CMAKE_EXE_LINKER_FLAGS_REFERENCE)

# turn off all fancy options for the regressiontests reference build
if("${CMAKE_BUILD_TYPE}" STREQUAL "Reference")
    set(GMX_GPU OFF CACHE BOOL "Disabled for regressiontests reference builds" FORCE)
    set(GMX_OPENMP OFF CACHE BOOL "Disabled for regressiontests reference builds" FORCE)
    set(GMX_SIMD "None" CACHE STRING "Disabled for regressiontests reference builds" FORCE)
    set(GMX_FFT_LIBRARY "fftpack" CACHE STRING "Use fftpack for regressiontests reference builds" FORCE)
    set(GMX_THREAD_MPI OFF CACHE BOOL "Disabled for regressiontests reference builds" FORCE)

    math(EXPR _major_version_too_high "${GMX_GCC_MINIMUM_REQUIRED_VERSION}+1")
    if(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_CXX_COMPILER_VERSION}" VERSION_GREATER_EQUAL ${_major_version_too_high})
        message(WARNING "Reference values for regressiontests should use GROMACS compiled with "
            "gcc ${GMX_GCC_MINIMUM_REQUIRED_VERSION}.x, but your configuration is using ${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}.")
    endif()
endif()
