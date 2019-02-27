#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2019, by the GROMACS development team, led by
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

cmake_minimum_required(VERSION 3.9.6)

# This needs to be set before project() in order to pick up toolchain files
#list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake)

# OS X deployment target should be >=10.9 for modern C++ compatibility.
set(CMAKE_OSX_DEPLOYMENT_TARGET 10.9 CACHE STRING
    "OS X deployment target below 10.9 does not use C++11 standard library"
    FORCE)
set(CMAKE_OSX_ARCHITECTURES x86_64 CACHE STRING
    "Python package should only be built in 64-bit architecture"
    FORCE)

project(gmxapi VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Only interpret if() arguments as variables or keywords when unquoted.
cmake_policy(SET CMP0054 NEW)
# honor the language standard settings for try_compile()
cmake_policy(SET CMP0067 NEW)

# It may help to disable extensions to encourage consistent standard library implementations,
# but I don't know whether this would disable necessary attributes for symbol visibility
# or unused variable annotation, etc. TBD. Also, may be incompatible with pybind config (also TBD).
# See https://crascit.com/2015/03/28/enabling-cxx11-in-cmake/
#set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory(external/pybind)

pybind11_add_module(_gmxapi gmxapi/_gmxapi.cpp)

install(TARGETS _gmxapi LIBRARY DESTINATION gmxapi)

if(SKBUILD)
    install(TARGETS _gmxapi LIBRARY DESTINATION gmxapi)
else()
    message(WARNING "CMake should be run by scikit-build through setup.py to get pure Py sources.")
endif()
