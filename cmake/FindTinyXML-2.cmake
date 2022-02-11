#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2016- The GROMACS Authors
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

# This package tries to find the TinyXML-2 library. It will consider
# the directory ${TinyXML2_DIR} in its search. Upon exit, the
# following variables may be set:
#
# TinyXML2_FOUND       - TinyXML-2 was found
# TinyXML2_INCLUDE_DIR - TinyXML-2 include directory
# TinyXML2_LIBRARIES   - TinyXML-2 libraries
# TinyXML2_LINKS_OK    - TinyXML-2 libraries link correctly
# TinyXML2_VERSION     - TinyXML-2 version string as "major.minor"

include(CMakePushCheckState)
cmake_push_check_state()

find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
    if(TinyXML2_FIND_VERSION)
        if(TinyXML2_FIND_VERSION_EXACT)
            pkg_check_modules(PC_TINYXML2 QUIET tinyxml2=${TinyXML2_FIND_VERSION})
        else()
            pkg_check_modules(PC_TINYXML2 QUIET tinyxml2>=${TinyXML2_FIND_VERSION})
        endif()
    else()
        pkg_check_modules(PC_TINYXML2 QUIET tinyxml2)
    endif()
endif()

# Try to find tinyxml2, perhaps with help from pkg-config
find_path(TinyXML2_INCLUDE_DIR tinyxml2.h HINTS "${TinyXML2_DIR}" "${PC_TINYXML2_INCLUDE_DIRS}" PATH_SUFFIXES include)
find_library(TinyXML2_LIBRARIES NAMES tinyxml2 HINTS "${TinyXML2_DIR}" "${PC_TINYXML2_LIBRARIES}" PATH_SUFFIXES lib64 lib)

if(TinyXML2_INCLUDE_DIR AND EXISTS "${TinyXML2_INCLUDE_DIR}/tinyxml2.h" AND TinyXML2_LIBRARIES)
    file(STRINGS "${TinyXML2_INCLUDE_DIR}/tinyxml2.h" _TinyXML2_H_MAJOR REGEX "TIXML2_MAJOR_VERSION = [0-9]+;")
    file(STRINGS "${TinyXML2_INCLUDE_DIR}/tinyxml2.h" _TinyXML2_H_MINOR REGEX "TIXML2_MINOR_VERSION = [0-9]+;")
    string(REGEX REPLACE "TIXML2_MAJOR_VERSION = ([0-9]+);" "\\1" _TinyXML2_MAJOR_VERSION "${_TinyXML2_H_MAJOR}")
    string(REGEX REPLACE "TIXML2_MINOR_VERSION = ([0-9]+);" "\\1" _TinyXML2_MINOR_VERSION "${_TinyXML2_H_MINOR}")
    set(TinyXML2_VERSION "${_TinyXML2_MAJOR_VERSION}.${_TinyXML2_MINOR_VERSION}")
endif()

# Make sure we can also link, so that cross-compilation is properly supported
if (TinyXML2_INCLUDE_DIR AND TinyXML2_LIBRARIES)
    include(CheckCXXSourceCompiles)
    set(CMAKE_REQUIRED_INCLUDES ${TinyXML2_INCLUDE_DIR})
    set(CMAKE_REQUIRED_LIBRARIES ${TinyXML2_LIBRARIES})
    check_cxx_source_compiles("#include <tinyxml2.h>\nint main(){tinyxml2::XMLDocument doc;}" TinyXML2_LINKS_OK)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TinyXML2
    FOUND_VAR
    TinyXML2_FOUND
    REQUIRED_VARS
    TinyXML2_INCLUDE_DIR
    TinyXML2_LIBRARIES
    TinyXML2_LINKS_OK
    VERSION_VAR
    TinyXML2_VERSION)

mark_as_advanced(TinyXML2_INCLUDE_DIR TinyXML2_LIBRARIES)

cmake_pop_check_state()
