#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2019- The GROMACS Authors
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

# This CMakeLists.txt is not intended to be used directly, but either through
# setup.py or as an inclusion of the full GROMACS project.
# See https://manual.gromacs.org/current/gmxapi/userguide/install.html for more.
cmake_minimum_required(VERSION 3.18.4)

# Note that this is the gmxapi._gmxapi Python bindings package version,
# not the C++ API version. It is not essential that it match the pure Python
# package version, but is likely to do so.
project(gmxapi)

# Check if Python package is being built directly or via add_subdirectory
set(GMXAPI_MASTER_PROJECT OFF)
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(GMXAPI_MASTER_PROJECT ON)
    if (NOT Python3_FIND_STRATEGY)
        # If the user provides a hint for the Python installation with Python3_ROOT_DIR,
        # prevent FindPython3 from overriding the choice with a newer Python version
        # when CMP0094 is set to OLD.
        set(Python3_FIND_STRATEGY LOCATION)
    endif ()
    if(NOT Python3_FIND_VIRTUALENV)
        # We advocate using Python venvs to manage package availability, so by default
        # we want to preferentially discover user-space software.
        set(Python3_FIND_VIRTUALENV FIRST)
    endif()
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

find_package(Python3 3.7 COMPONENTS Interpreter Development)
find_package(pybind11 2.6 QUIET CONFIG)
# If we are not running through setup.py, we may need to look for the pybind11 headers.
if (NOT pybind11_FOUND)
    execute_process(
        COMMAND
        "${Python3_EXECUTABLE}" -c
        "import pybind11; print(pybind11.get_cmake_dir())"
        OUTPUT_VARIABLE _tmp_dir
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    list(APPEND CMAKE_PREFIX_PATH "${_tmp_dir}")
    # The following should only run once, and if it runs more often than that, we would want to
    # know about it. So we don't use `QUIET` here.
    find_package(pybind11 2.6 CONFIG)
endif ()
if (NOT pybind11_FOUND)
    message(FATAL_ERROR "Python package build dependencies not found with interpreter ${Python3_EXECUTABLE}. "
            "See https://manual.gromacs.org/current/gmxapi/userguide/install.html")
endif ()

# Workaround for issue #4563 for GROMACS releases that won't be patched:
# Find GROMACS package early.
if(GMXAPI_MASTER_PROJECT)
    find_package(GROMACS 2021 REQUIRED
                 NAMES gromacs gromacs_mpi gromacs_d gromacs_mpi_d
                 HINTS "$ENV{GROMACS_DIR}"
                 )
    if (NOT DEFINED GROMACS_IS_DOUBLE)
        message(AUTHOR_WARNING "GROMACS_IS_DOUBLE undefined.")
    endif ()
else()
    if (NOT DEFINED GMX_DOUBLE)
        message(AUTHOR_WARNING "GMX_DOUBLE undefined.")
    endif ()
    set(GROMACS_IS_DOUBLE ${GMX_DOUBLE})
endif()

if(GMXAPI_MASTER_PROJECT)
    find_package(gmxapi REQUIRED)
    if (gmxapi_VERSION VERSION_LESS 0.2.1)
        message(
            FATAL_ERROR
            "Your GROMACS installation is too old. This package requires GROMACS 2021.3 or higher.")
    endif ()
else()
    # Building as part of the GROMACS master project. GROMACS CMake logic should
    # not be processing this unless Python3 was appropriately detected.
    if (NOT Python3_FOUND)
        message(FATAL_ERROR "Error in CMake script. Please report GROMACS bug.")
    endif ()

    get_target_property(gmxapi_VERSION gmxapi VERSION)
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
endif()

message(STATUS "Configuring Python package for gmxapi version ${gmxapi_VERSION}${_suffix}")

# The Gromacs::gmxapi target could be imported from an existing installation or
# provided as an alias target within the GROMACS build tree.
if (NOT TARGET Gromacs::gmxapi)
    message(FATAL_ERROR "Cannot build Python package without GROMACS gmxapi support.")
endif ()

# TODO(#3279): Provide user hints for mpi4py installation.
# Note that neither the Python package nor the Gromacs::gmxapi CMake target are
# built with MPI in any case, but they _should_ be built with a C++ compiler
# that is compatible with the available MPI compiler wrappers, and technically
# _that_ is what we want to help the user identify when installing mpi4py, even
# if libgromacs is not built with MPI support either.
# For convenience, it is fine if libgmxapi and _gmxapi are built with the mpi
# compiler wrapper.

pybind11_add_module(_gmxapi
                    src/cpp/launch_0_2_1.cpp
                    src/cpp/module.cpp
                    src/cpp/export_context.cpp
                    src/cpp/export_system.cpp
                    src/cpp/export_tprfile.cpp
                    src/cpp/gmxpy_exceptions.cpp
                    src/cpp/pycontext.cpp
                    src/cpp/pysystem.cpp
                    )

if (NOT DEFINED gmxapi_VERSION)
    get_target_property(gmxapi_VERSION Gromacs::gmxapi VERSION)
endif ()
if (NOT DEFINED gmxapi_VERSION)
    message(FATAL_ERROR "gmxapi library version not detected.")
endif ()
if (gmxapi_VERSION VERSION_LESS 0.3.1)
    message(WARNING "Found an old gmxapi library version. Please consider updating your GROMACS installation.")
    target_sources(_gmxapi PRIVATE src/cpp/wrapped_exceptions_0_1_0.cpp)
else()
    target_sources(_gmxapi PRIVATE src/cpp/wrapped_exceptions_0_3_1.cpp)
endif()

target_include_directories(_gmxapi PRIVATE
                           ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp
                           ${CMAKE_CURRENT_BINARY_DIR}/src/cpp
                           )

# RPATH management: make sure build artifacts can find GROMACS library.
set_target_properties(_gmxapi PROPERTIES SKIP_BUILD_RPATH FALSE)

# Get details of GROMACS installation needed by the Python package at run time.

# Get the MPI capability.
get_target_property(_gmx_mpi Gromacs::gmxapi MPI)
if (${_gmx_mpi} STREQUAL "library")
    set(_gmx_mpi_type "\"library\"")
elseif(${_gmx_mpi} STREQUAL "tmpi")
    set(_gmx_mpi_type "\"tmpi\"")
elseif(${_gmx_mpi} STREQUAL "none")
    set(_gmx_mpi_type "null")
else()
    message(FATAL_ERROR "Unrecognized gmxapi MPI value: ${_gmx_mpi}")
endif ()
unset(_gmx_mpi)
# Get the path of the command line entry point and binary install directory.
if (NOT TARGET Gromacs::gmx)
    message(FATAL_ERROR "GROMACS command line tool not found.")
endif ()
get_target_property(_gmx_executable_imported Gromacs::gmx IMPORTED)
if (_gmx_executable_imported)
    get_target_property(_gmx_executable Gromacs::gmx LOCATION)
    get_filename_component(_gmx_bindir ${_gmx_executable} DIRECTORY)
    message(STATUS "Imported ${_gmx_bindir} executable.")
    unset(_gmx_executable_imported)
else()
    get_target_property(_gmx_bindir Gromacs::gmx RUNTIME_OUTPUT_DIRECTORY)
    get_target_property(_gmx_executable Gromacs::gmx OUTPUT_NAME)
    set(_gmx_executable "${_gmx_bindir}/${_gmx_executable}")
    message(STATUS "Using ${_gmx_executable} from build tree.")
endif ()
if (NOT _gmx_bindir OR NOT _gmx_executable)
    message(FATAL_ERROR "Could not get path for gmx wrapper binary.")
endif ()
if (GROMACS_IS_DOUBLE)
    set(_gmx_double "true")
else()
    set(_gmx_double "false")
endif()
configure_file(src/gmxapi/gmxconfig.json.in src/gmxapi/gmxconfig.json)
unset(_gmx_executable)
unset(_gmx_bindir)
unset(_gmx_mpi_type)
unset(_gmx_double)

if (GMXAPI_MASTER_PROJECT)
    set_target_properties(_gmxapi PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)
    set_target_properties(_gmxapi PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
    target_link_libraries(_gmxapi PRIVATE Gromacs::gmxapi)
    # The Python setup.py sets CMAKE_LIBRARY_OUTPUT_DIRECTORY and will be looking for generated files there.
    file(COPY ${CMAKE_CURRENT_BINARY_DIR}/src/gmxapi/gmxconfig.json
         DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
else()
    # The rest of the logic in this conditional is to support the GMX_PYTHON_PACKAGE option
    # for testing the gmxapi Python packages within a full GROMACS project build_command and.
    # for building full GROMACS project documentation.
    #
    # Note that file copying occurs during the CMake configure phase. During development,
    # edits may not trigger a rerun of the CMake configure phase for targets like
    # `gmxapi-pytest` (enabled by -DGMX_PYTHON_PACKAGE=ON). Manual re-run of `cmake`
    # may be necessary between edits and testing.

    set(GMXAPI_PYTHON_STAGING_DIR ${CMAKE_CURRENT_BINARY_DIR}/gmxapi_staging)
    # Instead, we should probably build a source package and alert the user of its location.
    # We can use CMake to call the Python packaging tools to create an 'sdist'
    # source distribution archive to be installed in the GROMACS installation
    # destination. We can use the build directory as the working directory for
    # easier clean-up, as well.
    # TODO: (ref Issue #2896) Build and install 'sdist' with GROMACS.

    # The Python module is being built against GROMACS in its build tree, so we will not install.
    set_target_properties(_gmxapi PROPERTIES BUILD_WITH_INSTALL_RPATH FALSE)
    target_link_libraries(_gmxapi PRIVATE Gromacs::gmxapi)
    # However, we can still produce an importable package for documentation builds and
    # basic testing in ${CMAKE_CURRENT_BINARY_DIR}/gmxapi_staging
    set_target_properties(_gmxapi PROPERTIES
                          LIBRARY_OUTPUT_DIRECTORY ${GMXAPI_PYTHON_STAGING_DIR}/gmxapi)
    file(GLOB_RECURSE _py_sources
         CONFIGURE_DEPENDS
         ${CMAKE_CURRENT_SOURCE_DIR}/src/gmxapi/*.py)
    foreach(_package_file IN LISTS _py_sources)
        get_filename_component(_absolute_dir ${_package_file} DIRECTORY)
        file(RELATIVE_PATH _relative_dir ${CMAKE_CURRENT_SOURCE_DIR}/src ${_absolute_dir})
        file(COPY ${_package_file} DESTINATION ${GMXAPI_PYTHON_STAGING_DIR}/${_relative_dir})
    endforeach()
    file(COPY setup.py CMakeLists.txt DESTINATION ${GMXAPI_PYTHON_STAGING_DIR})
    file(COPY ${CMAKE_CURRENT_BINARY_DIR}/src/gmxapi/gmxconfig.json DESTINATION ${GMXAPI_PYTHON_STAGING_DIR}/gmxapi)

    # Unit test and build docs using PYTHONPATH=$CMAKE_CURRENT_BINARY_DIR/gmxapi_staging
    set_target_properties(_gmxapi PROPERTIES staging_dir ${GMXAPI_PYTHON_STAGING_DIR})
    # Note: Integration testing for multiple Python versions and/or CMake-driven
    # sdist preparation could be performed with CMake custom_commands and custom_targets.
endif()

# When building as part of GROMACS umbrella project, add a testing target
# to the `check` target. Normal usage is to first install the Python package,
# then run `pytest` on the `tests` directory. Refer to gmxapi package documentation.
if(NOT GMXAPI_MASTER_PROJECT)
	if (BUILD_TESTING)
		add_subdirectory(test)
	endif()
endif()
