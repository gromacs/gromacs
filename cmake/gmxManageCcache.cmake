#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2018- The GROMACS Authors
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

# Permit the use of ccache (when available in the system path or
# CMAKE_PREFIX_PATH), which wraps the CMAKE_C_COMPILER and
# CMAKE_CXX_COMPILER to speed up build times. Reference
# https://ccache.samba.org and
# https://crascit.com/2016/04/09/using-ccache-with-cmake/

option(GMX_ENABLE_CCACHE "Allow CMake to use Ccache compiler wrappers if available." OFF)

if(NOT GMX_ENABLE_CCACHE)
    return()
endif()

if(GMX_CLANG_TIDY)
    message(FATAL_ERROR "Ccache does not work with the wrapper script used for "
        "clang-tidy builds. Use -DGMX_ENABLE_CCACHE=off.")
endif()
if(GMX_CLANG_ANALYZER)
    message(FATAL_ERROR "Ccache does not work with the wrapper script used for "
        "clang-analyzer builds. Use -DGMX_ENABLE_CCACHE=off.")
endif()


# Find ccache and set it as the compiler launcher if available
find_program(CCACHE_PROGRAM ccache)

if(CCACHE_PROGRAM)
    # Set ccache as the launcher for both C and C++ compilers
    set(CMAKE_C_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
    set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
else()
    message(FATAL_ERROR "Ccache is not installed or could not be found. Please install Ccache to proceed, or configure the project to build without Ccache support.")
endif()
