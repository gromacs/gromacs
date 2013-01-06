#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012, by the GROMACS development team, led by
# David van der Spoel, Berk Hess, Erik Lindahl, and including many
# others, as listed in the AUTHORS file in the top-level source
# directory and at http://www.gromacs.org.
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
# the name of the target operating system
set(CMAKE_SYSTEM_NAME Linux CACHE STRING "Cross-compiling for Fujitsu Sparc64")

set_property(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS FALSE)

# set the compiler
set(CMAKE_C_COMPILER mpifccpx)
set(CMAKE_CXX_COMPILER mpiFCCpx)
set(CMAKE_C_COMPILER_ID "Fujitsu" CACHE STRING "Prevent CMake from adding GNU-specific linker flags (-rdynamic)" FORCE)

set(CMAKE_C_FLAGS "-Kopenmp -Kfast,reduction,swp,simd=2,uxsimd -x500 -Xg -DGMX_RELAXED_DOUBLE_PRECISION -w" CACHE STRING "Fujitsu Sparc64 C Flags" FORCE)
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS}" CACHE STRING "Fujitsu Sparc64 C++ Flags" FORCE)
set(GMX_SOFTWARE_INVSQRT OFF CACHE BOOL "Use native 1.0/sqrt(x) on Fujitsu Sparc64" FORCE)

set(GMX_THREAD_MPI OFF CACHE BOOL "Use real MPI instead" FORCE)
set(GMX_MPI ON CACHE BOOL "Use MPI library" FORCE)
set(GMX_DOUBLE ON CACHE BOOL "Use double by default on Fujitsu Sparc64 (due to HPC-ACE)" FORCE)
set(GMX_GPU OFF CACHE BOOL "Cannot do GPU acceleration on Fujitsu Sparc64" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Use static linking by default on Fujitsu Sparc64" FORCE)

set(GMX_CPU_ACCELERATION "Sparc64_HPC_ACE" CACHE STRING "Enabling Sparc64 HPC-ACE acceleration when using Fujitsu Sparc64 toolchain")
