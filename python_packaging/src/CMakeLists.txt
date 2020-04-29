#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2019,2020, by the GROMACS development team, led by
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
# Reference https://gitlab.com/gromacs/gromacs/-/issues/2896 for additional discussion.
cmake_minimum_required(VERSION 3.13.0)

# This needs to be set before project() in order to pick up toolchain files
#list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake)

# OS X deployment target should be >=10.14 for modern C++ compatibility.
# Reference https://scikit-build.readthedocs.io/en/latest/generators.html#macosx
# and https://github.com/MacPython/wiki/wiki/Spinning-wheels
set(CMAKE_OSX_DEPLOYMENT_TARGET 10.14 CACHE STRING
    "OS X deployment target below 10.14 does not use modern standard library"
    FORCE)
set(CMAKE_OSX_ARCHITECTURES x86_64 CACHE STRING
    "OS X should build Python package for 64-bit architecture"
    FORCE)

project(gmxapi VERSION 0.1.0)

# Check if Python package is being built directly or via add_subdirectory
set(GMXAPI_MASTER_PROJECT OFF)
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(GMXAPI_MASTER_PROJECT ON)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Only interpret if() arguments as variables or keywords when unquoted.
cmake_policy(SET CMP0054 NEW)
# honor the language standard settings for try_compile()
cmake_policy(SET CMP0067 NEW)
if(POLICY CMP0074) #3.12
    # Allow gmxapi_ROOT hint.
    cmake_policy(SET CMP0074 NEW)
endif()

if(GMXAPI_MASTER_PROJECT)
    find_package(gmxapi 0.0.8 REQUIRED
                 HINTS "$ENV{GROMACS_DIR}"
                 )
endif()
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

# The Gromacs::gmxapi target could be imported from an existing installation or
# provided as an alias target within the GROMACS build tree.
if (NOT TARGET Gromacs::gmxapi)
    message(FATAL_ERROR "Cannot build Python package without GROMACS gmxapi support.")
endif ()

# TODO: Provide user hints for mpi4py installation.
# Note that neither the Python package nor the Gromacs::gmxapi CMake target are
# built with MPI in any case, but they _should_ be built with a C++ compiler
# that is compatible with the available MPI compiler wrappers, and technically
# _that_ is what we want to help the user identify when installing mpi4py, even
# if libgromacs is not built with MPI support either.
# For convenience, it is fine if libgmxapi and _gmxapi are built with the mpi
# compiler wrapper.

option(GMXAPI_USE_BUNDLED_PYBIND
       "Use pybind11 headers bundled with this repository. If OFF, CMake does `find_package(pybind11)`."
       ON)
if(GMXAPI_USE_BUNDLED_PYBIND)
    add_subdirectory(external/pybind)
else()
    # Reference https://gitlab.com/gromacs/gromacs/-/issues/2896
    find_package(pybind11 2.2 REQUIRED)
endif()

set(GMXAPI_PYTHON_EXTENSION_SOURCES
    gmxapi/module.cpp
    gmxapi/export_context.cpp
    gmxapi/export_exceptions.cpp
    gmxapi/export_system.cpp
    gmxapi/export_tprfile.cpp
    gmxapi/pycontext.cpp
    gmxapi/pysystem.cpp
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

# Python sources (*.py) will be packaged by scikit-build and setuptools.
# Note that library targets are built in CMAKE_LIBRARY_OUTPUT_DIRECTORY if not otherwise specified.
# This may be an unexpected location, whether inherited from the GROMACS build tree
# or the SKBUILD framework. Note, also, that when scikit-build is invoked with setup.py,
# the CMake build takes place in a subdirectory of ./_skbuild/.
set(GMXAPI_PYTHON_STAGING_DIR ${CMAKE_CURRENT_BINARY_DIR}/gmxapi_staging)
set_target_properties(_gmxapi PROPERTIES
                      LIBRARY_OUTPUT_DIRECTORY ${GMXAPI_PYTHON_STAGING_DIR}/gmxapi)

# scikit-build sets SKBUILD when running Python packaging tools through setup.py
# (e.g. with pip)
if(SKBUILD)
    # The Python module is being built for a GROMACS installation.
    set_target_properties(_gmxapi PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)
    set_target_properties(_gmxapi PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
    target_link_libraries(_gmxapi PRIVATE Gromacs::gmxapi)
    # By default, scikit-build expects the library to be installed into a directory
    # named for the Python package as in setup.py.
    install(TARGETS _gmxapi LIBRARY DESTINATION gmxapi)
else()
    # The Python module is being built against GROMACS in its build tree.
    # Note: we do not have plans to install the staged package when SKBUILD != TRUE
    set_target_properties(_gmxapi PROPERTIES BUILD_WITH_INSTALL_RPATH FALSE)
    target_link_libraries(_gmxapi PRIVATE Gromacs::gmxapi)

    # TODO: Determine packaging and installation cases and implementation.
    # Reference https://gitlab.com/gromacs/gromacs/-/issues/2896 for additional discussion.
    # Currently, CMake should be run by scikit-build through setup.py for proper Python packaging.
    # We don't want to install by default in the outer scope of the GROMACS
    # CMake procedure because we could end up trying to install to a system directory
    # the user did not intend, or a user might install a Python-installation-specific
    # package into an overly generic GROMACS path.

    # Instead, we should probably build a source package and alert the user of its location.
    # We can use CMake to call the Python packaging tools to create an 'sdist'
    # source distribution archive to be installed in the GROMACS installation
    # destination. We can use the build directory as the working directory for
    # easier clean-up, as well.
    # TODO: (ref Issue #2896) Build and install 'sdist' with GROMACS.

    # However, we can still produce an importable package for documentation builds and
    # basic testing in ${CMAKE_CURRENT_BINARY_DIR}/gmxapi_staging
    if(CMAKE_VERSION VERSION_LESS 3.12)
        file(GLOB_RECURSE _py_sources
             ${CMAKE_CURRENT_SOURCE_DIR}/gmxapi/*.py)
    else()
        # CONFIGURE_DEPENDS appears in CMake 3.12 and can help to more robustly detect
        # the need to update anything depending on the staged package.
        file(GLOB_RECURSE _py_sources
             CONFIGURE_DEPENDS
             ${CMAKE_CURRENT_SOURCE_DIR}/gmxapi/*.py)
    endif()
    foreach(_package_file IN LISTS _py_sources)
        get_filename_component(_absolute_dir ${_package_file} DIRECTORY)
        file(RELATIVE_PATH _relative_dir ${CMAKE_CURRENT_SOURCE_DIR} ${_absolute_dir})
        file(COPY ${_package_file} DESTINATION gmxapi_staging/${_relative_dir})
    endforeach()
    file(COPY setup.py CMakeLists.txt DESTINATION ${GMXAPI_PYTHON_STAGING_DIR})
    # Set CMake variable pybind11_DIR to ${CMAKE_CURRENT_SOURCE_DIR}/external/pybind/tools
    # if re-invoking CMake (including via Python setuptools) for the files in gmxapi_staging.

    # Unit test and build docs using PYTHONPATH=$CMAKE_CURRENT_BINARY_DIR/gmxapi_staging
    set_target_properties(_gmxapi PROPERTIES staging_dir ${GMXAPI_PYTHON_STAGING_DIR})
    # Note: Integration testing for multiple Python versions and/or CMake-driven
    # sdist preparation could be performed with CMake custom_commands and custom_targets.
endif()

# When building as part of GROMACS umbrella project, add a testing target
# to the `check` target. Normal usage is to first install the Python package,
# then run `pytest` on the `tests` directory. Refer to gmxapi package documentation.
if(NOT GMXAPI_MASTER_PROJECT)
    add_subdirectory(test)
endif()
