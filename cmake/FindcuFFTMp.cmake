#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2022- The GROMACS Authors
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

# - Find cuFFTMp, NVIDIA's FFT library for multi-GPU multi-node
#
# This script does define cache variables
#   CUFFTMP_INCLUDE_DIR    - Location of cuFFTMP's include directory.
#   CUFFTMP_LIBRARY        - Location of cuFFTMP's libraries
#
# If your cuFFTMp installation is not in a standard installation
# directory, you may provide a hint to where it may be found. Simply
# set the value cuFFTMp_ROOT to the directory containing
# 'include/cufftMp.h" prior to calling find_package().

if(cuFFTMp_INCLUDE_DIR)
  # Already in cache, be silent
  set (cuFFTMp_FIND_QUIETLY TRUE)
endif()

find_package(PkgConfig)
pkg_check_modules(PC_cuFFTMp QUIET cuFFTMp)

find_path(cuFFTMp_ROOT_DIR
    NAMES 
    include/cufftmp_ea/cufftMp.h
    include/cufftmp/cufftMp.h
    HINTS ${cuFFTMp_ROOT}
    DOC "cuFFTMp root directory.")

find_path(cuFFTMp_INCLUDE_DIR
    NAMES cufftMp.h
    HINTS
    ${cuFFTMp_ROOT_DIR}/include/cufftmp_ea
    ${cuFFTMp_ROOT_DIR}/include/cufftmp
    ${PC_cuFFTMp_INCLUDE_DIRS}
    DOC "cuFFTMp Include directory")

find_library(cuFFTMp_LIBRARY
    NAMES cufftMp
    HINTS
    ${cuFFTMp_ROOT_DIR}/lib64
    ${cuFFTMp_ROOT_DIR}/lib
    ${PC_cuFFTMp_LIBRARY_DIRS}
    )

# handle the QUIETLY and REQUIRED arguments and set cuFFTMp_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(cuFFTMp
    REQUIRED_VARS cuFFTMp_INCLUDE_DIR cuFFTMp_LIBRARY
)

mark_as_advanced(cuFFTMp_ROOT_DIR cuFFTMp_LIBRARY cuFFTMp_INCLUDE_DIR)

