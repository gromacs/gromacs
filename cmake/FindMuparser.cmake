#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2021, by the GROMACS development team, led by
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

# This package tries to find an external muparser library, version
# 2.3.
#
# MUPARSER_FOUND       - muparser was found
# MUPARSER_INCLUDE_DIR - muparser include directory
# MUPARSER_LIBRARY     - muparser library
# MUPARSER_LINKS_OK    - muparser libraries link correctly
# MUPARSER_VERSION     - muparser version string as "major.minor"
#
# CMake will search the CMAKE_PREFIX_PATH in the usual way, but if you
# need more control then setting MUPARSER_INCLUDE_DIR and MUPARSER_LIBRARY
# on the cmake command line to suitable values will work.

include(CMakePushCheckState)
cmake_push_check_state()

find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
    if(MUPARSER_FIND_VERSION)
        if(MUPARSER_FIND_VERSION_EXACT)
            pkg_check_modules(PC_MUPARSER QUIET muparser=${MUPARSER_FIND_VERSION})
        else()
            pkg_check_modules(PC_MUPARSER QUIET muparser>=${MUPARSER_FIND_VERSION})
        endif()
    else()
        pkg_check_modules(PC_MUPARSER QUIET muparser)
        if (PC_MUPARSER_VERSION)
            string(REGEX REPLACE "^([0-9]+):([0-9]+)" "\\1.\\2" MUPARSER_VERSION "${PC_MUPARSER_VERSION}")
        endif()
    endif()
endif()

# Try to find muparser, perhaps with help from pkg-config
find_path(MUPARSER_INCLUDE_DIR muParser.h HINTS "${PC_MUPARSER_INCLUDE_DIRS}" PATH_SUFFIXES include)
find_library(MUPARSER_LIBRARY NAMES muparser HINTS "${PC_MUPARSER_LIBRARY_DIRS}" PATH_SUFFIXES lib64 lib)

# Make sure we can also link, so that cross-compilation is properly supported
if (MUPARSER_INCLUDE_DIR AND MUPARSER_LIBRARY)
    include(CheckCXXSourceCompiles)
    set(CMAKE_REQUIRED_INCLUDES ${MUPARSER_INCLUDE_DIR})
    set(CMAKE_REQUIRED_LIBRARIES ${MUPARSER_LIBRARY})
    check_cxx_source_compiles("#include <muParser.h>\nint main(){mu::Parser parser;}" MUPARSER_LINKS_OK)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(muparser
    FOUND_VAR
    MUPARSER_FOUND
    REQUIRED_VARS
    MUPARSER_INCLUDE_DIR
    MUPARSER_LIBRARY
    MUPARSER_LINKS_OK
    VERSION_VAR
    MUPARSER_VERSION)

mark_as_advanced(MUPARSER_INCLUDE_DIR MUPARSER_LIBRARY)

# Make a target that other targets can depend on just like this was a
# library built in the main project.
if (MUPARSER_FOUND)
    add_library(muparser INTERFACE IMPORTED)
    set_target_properties(muparser PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${MUPARSER_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "${MUPARSER_LIBRARY}"
        )
endif()

cmake_pop_check_state()
