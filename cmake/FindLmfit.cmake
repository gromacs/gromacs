#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2016, by the GROMACS development team, led by
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

# This package tries to find an external lmfit library. It is intended
# to work with pkg-config, because that is the mechanism supported in
# lmfit. Upon exit, the following variables may be set:
#
# LMFIT_FOUND       - lmfit was found
# LMFIT_INCLUDE_DIR - lmfit include directory
# LMFIT_LIBRARIES   - lmfit libraries
# LMFIT_LINKS_OK    - lmfit libraries link correctly
# LMFIT_VERSION     - lmfit version string as "major.minor"
#
# If you cannot use pkg-config for some reason, then setting
# LMFIT_INCLUDE_DIRS and LMFIT_LIBRARY_DIRS on the cmake command line
# to suitable values will work.

include(CMakePushCheckState)
cmake_push_check_state()

find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
    if(LMFIT_FIND_VERSION)
        # lmfit doesn't support CMake-based find_package version
        # checking in 6.1, so this code does nothing.
        if(LMFIT_FIND_VERSION_EXACT)
            pkg_check_modules(PC_LMFIT lmfit=${LMFIT_FIND_VERSION})
        else()
            pkg_check_modules(PC_LMFIT lmfit>=${LMFIT_FIND_VERSION})
        endif()
    else()
        pkg_check_modules(PC_LMFIT lmfit)
        if (PC_LMFIT_VERSION)
            string(REGEX REPLACE "^([0-9]+):([0-9]+)" "\\1.\\2" LMFIT_VERSION "${PC_LMFIT_VERSION}")
        endif()
    endif()
endif()

# Try to find lmfit, perhaps with help from pkg-config
find_path(LMFIT_INCLUDE_DIRS lmcurve.h HINTS "${PC_LMFIT_INCLUDE_DIRS}" PATH_SUFFIXES include)
find_library(LMFIT_LIBRARY_DIRS NAMES lmfit HINTS "${PC_LMFIT_LIBRARY_DIRS}" PATH_SUFFIXES lib64 lib)

# Make sure we can also link, so that cross-compilation is properly supported
if (LMFIT_INCLUDE_DIRS AND LMFIT_LIBRARY_DIRS)
    include(CheckCXXSourceCompiles)
    set(CMAKE_REQUIRED_INCLUDES ${LMFIT_INCLUDE_DIRS})
    set(CMAKE_REQUIRED_LIBRARIES ${LMFIT_LIBRARY_DIRS})
    check_cxx_source_compiles("#include <lmcurve.h>\nint main(){lmcurve(0,0,0,0,0,0,0,0);}" LMFIT_LINKS_OK)
endif()

if (LMFIT_LINKS_OK)
    set(LMFIT_INCLUDE_DIR ${LMFIT_INCLUDE_DIRS})
    set(LMFIT_LIBRARIES ${LMFIT_LIBRARY_DIRS})
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(lmfit
    FOUND_VAR
    LMFIT_FOUND
    REQUIRED_VARS
    LMFIT_INCLUDE_DIR
    LMFIT_LIBRARIES
    LMFIT_LINKS_OK
    VERSION_VAR
    LMFIT_VERSION)

mark_as_advanced(LMFIT_INCLUDE_DIRS LMFIT_LIBRARY_DIRS)

cmake_pop_check_state()
