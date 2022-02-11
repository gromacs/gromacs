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

# - Find clFFT, AMD's OpenCL FFT library
#
# This script does define cache variables
#   CLFFT_INCLUDE_DIR    - Location of clFFT's include directory.
#   CLFFT_LIBRARY        - Location of clFFT's libraries
# however the preferred use is simply to link to the imported
# CMake target clFFT, which is constructed to behave as if clFFT
# was built within the parent project, and correctly populates
# the required include directories and linker flags.
#
# If your clFFT installation is not in a standard installation
# directory, you may provide a hint to where it may be found. Simply
# set the value clFFT_ROOT to the directory containing
# 'include/clFFT.h" prior to calling find_package().

if(clFFT_INCLUDE_DIR)
  # Already in cache, be silent
  set (clFFT_FIND_QUIETLY TRUE)
endif()

find_package(PkgConfig)
pkg_check_modules(PC_clFFT QUIET clFFT)

find_path(clFFT_ROOT_DIR
    NAMES include/clFFT.h
    HINTS ${clFFT_ROOT}
    DOC "clFFT root directory.")

find_path(clFFT_INCLUDE_DIR
    NAMES clFFT.h
    HINTS
    ${clFFT_ROOT_DIR}/include
    ${PC_clFFT_INCLUDE_DIRS}
    DOC "clFFT Include directory")

find_library(clFFT_LIBRARY
    NAMES clFFT
    HINTS
    ${clFFT_ROOT_DIR}/lib64
    ${clFFT_ROOT_DIR}/lib
    ${PC_clFFT_LIBRARY_DIRS}
    )

# handle the QUIETLY and REQUIRED arguments and set clFFT_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(clFFT
    REQUIRED_VARS clFFT_INCLUDE_DIR clFFT_LIBRARY
)

mark_as_advanced(clFFT_ROOT_DIR clFFT_LIBRARY clFFT_INCLUDE_DIR)

# Prepare a faux link target that can be used as if the clFFT
# that was found was actually built in this project.
if(clFFT_FOUND)
    add_library(clFFT INTERFACE IMPORTED)
    target_link_libraries(clFFT INTERFACE "${clFFT_LIBRARY}" "${CMAKE_DL_LIBS}")
    target_include_directories(clFFT SYSTEM BEFORE INTERFACE "${clFFT_INCLUDE_DIR}")
endif()

