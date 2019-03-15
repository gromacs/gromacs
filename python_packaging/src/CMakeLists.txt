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

# This CMakeLists.txt allows source distributions of the gmxapi Python package
# to rely on scikit-build for support of various Python packaging systems. The
# simplest use case is to allow the `setup.py` file to invoke skbuild to
# configure and run CMake. CMake could be invoked directly by the user or a
# parent package, but the Python distribution would not be packaged automatically.
# Reference https://redmine.gromacs.org/issues/2896 for additional discussion.
cmake_minimum_required(VERSION 3.9.6)

# This needs to be set before project() in order to pick up toolchain files
#list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake)

# OS X deployment target should be >=10.9 for modern C++ compatibility.
# Reference https://scikit-build.readthedocs.io/en/latest/generators.html#macosx
# and https://github.com/MacPython/wiki/wiki/Spinning-wheels
set(CMAKE_OSX_DEPLOYMENT_TARGET 10.9 CACHE STRING
    "OS X deployment target below 10.9 does not use modern standard library"
    FORCE)
set(CMAKE_OSX_ARCHITECTURES x86_64 CACHE STRING
    "OS X should build Python package for 64-bit architecture"
    FORCE)

project(gmxapi VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Only interpret if() arguments as variables or keywords when unquoted.
cmake_policy(SET CMP0054 NEW)
# honor the language standard settings for try_compile()
cmake_policy(SET CMP0067 NEW)

find_package(gmxapi 0.0.8 REQUIRED
             HINTS "$ENV{GROMACS_DIR}"
             )
if(gmxapi_FOUND)
    set(_suffix "")
    # GROMACS master branch and development branches may have divergent
    # pre-release APIs. This check allows us to distinguish them and behave
    # differently if needed. github.com/kassonlab/gromacs-gmxapi devel branch
    # sets gmxapi_EXPERIMENTAL=TRUE. Upstream GROMACS master branch does not.
    # Ref: https://github.com/kassonlab/gmxapi/issues/166
    if(gmxapi_EXPERIMENTAL)
        set(_suffix " (unofficial)")
    endif()
    message(STATUS "Found gmxapi version ${gmxapi_VERSION}${_suffix}")
endif()

# TODO: Account for externally-provided pybind headers where recommended by packaging system.
# Reference https://redmine.gromacs.org/issues/2896
add_subdirectory(external/pybind)

set(GMXAPI_PYTHON_EXTENSION_SOURCES
    gmxapi/module.cpp
    )

pybind11_add_module(_gmxapi
                    ${GMXAPI_PYTHON_EXTENSION_SOURCES}
                    )

target_include_directories(_gmxapi PRIVATE
                           ${CMAKE_CURRENT_SOURCE_DIR}/gmxapi
                           ${CMAKE_CURRENT_BINARY_DIR}/gmxapi
                           )

# RPATH management: make sure build artifacts can find GROMACS library.
set_target_properties(_gmxapi PROPERTIES SKIP_BUILD_RPATH FALSE)
set_target_properties(_gmxapi PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)
set_target_properties(_gmxapi PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

target_link_libraries(_gmxapi PRIVATE Gromacs::gmxapi)

install(TARGETS _gmxapi LIBRARY DESTINATION gmxapi)

if(NOT SKBUILD)
    # TODO: Determine packaging and installation cases and implementation.
    # Reference https://redmine.gromacs.org/issues/2896 for additional discussion.
    message(WARNING "CMake should be run by scikit-build through setup.py for proper Python packaging.")
    # Try to do something sensible. Collect Python sources as setup.py would.
    file(GLOB_RECURSE _py_sources CONFIGURE_DEPENDS $CMAKE_CURRENT_SOURCE_DIR/gmxapi/*.py)
    install(FILES ${_py_sources} DESTINATION gmxapi)
endif()
