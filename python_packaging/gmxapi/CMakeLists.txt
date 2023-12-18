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
set(GMXAPI_MAIN_PROJECT OFF)
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(GMXAPI_MAIN_PROJECT ON)
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
if (POLICY CMP0139) #3.24
    cmake_policy(SET CMP0139 NEW)
endif ()

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

if(NOT gmxapi_ROOT AND ENV{gmxapi_ROOT})
    set(gmxapi_ROOT $ENV{gmxapi_ROOT})
endif()

# Workaround for issue #4563 for GROMACS releases (<2022) that won't be patched:
# Find GROMACS package early.
if(GMXAPI_MAIN_PROJECT)
    find_package(GROMACS 2021 REQUIRED
                 NAMES gromacs${GROMACS_SUFFIX} gromacs gromacs_mpi gromacs_d gromacs_mpi_d
                 HINTS "$ENV{GROMACS_DIR}" ${gmxapi_ROOT}
                 )
    if (NOT DEFINED GROMACS_IS_DOUBLE)
        message(AUTHOR_WARNING "GROMACS_IS_DOUBLE undefined.")
    endif ()
    if (NOT DEFINED GROMACS_SUFFIX)
        message(AUTHOR_WARNING "GROMACS_SUFFIX undefined.")
    endif ()
else()
    if (NOT DEFINED GMX_DOUBLE)
        message(AUTHOR_WARNING "GMX_DOUBLE undefined.")
    endif ()
    set(GROMACS_IS_DOUBLE ${GMX_DOUBLE})
    if (NOT DEFINED GMX_LIBS_SUFFIX)
        message(AUTHOR_WARNING "GMX_LIBS_SUFFIX undefined.")
    endif ()
    set(GROMACS_SUFFIX ${GMX_LIBS_SUFFIX})
endif()

if(GMXAPI_MAIN_PROJECT)
    find_package(gmxapi REQUIRED
                 HINTS "$ENV{GROMACS_DIR}" ${GROMACS_ROOT}
                 )
    if (gmxapi_VERSION VERSION_LESS 0.2.1)
        message(
            FATAL_ERROR
            "Your GROMACS installation is too old. This package requires GROMACS 2021.3 or higher.")
    endif ()
else()
    # Building as part of the GROMACS base project. GROMACS CMake logic should
    # not be processing this unless Python3 was appropriately detected.
    if (NOT Python3_FOUND)
        message(FATAL_ERROR "Error in CMake script. Please report GROMACS bug.")
    endif ()

    get_target_property(gmxapi_VERSION gmxapi VERSION)
endif()

# Note: this isn't exactly what we want, and we also need a GROMACS_DIR hint.
# TODO: Consider the necessary hints and documentation.
if(NOT gmxapi_ROOT AND gmxapi_DIR)
    set(gmxapi_ROOT "${gmxapi_DIR}")
endif()

if(gmxapi_FOUND)
    set(_suffix "")
    # GROMACS main branch and development branches may have divergent
    # pre-release APIs. This check allows us to distinguish them and behave
    # differently if needed. github.com/kassonlab/gromacs-gmxapi devel branch
    # sets gmxapi_EXPERIMENTAL=TRUE. Upstream GROMACS main branch does not.
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

#
# Get details of GROMACS installation needed by the Python package at run time.
#

# Get the MPI capability.
get_target_property(_gmx_mpi Gromacs::gmxapi MPI)
if (${_gmx_mpi} STREQUAL "library")
    set(_gmx_mpi_type "\"library\"")
    set(HAS_GROMACS_MPI TRUE)
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
unset(_gmxapi_level)

# If libgromacs is not using an MPI library, then no MPI CMake targets
# are defined. But gmxapi can optionally be built against MPI to
# direct work done by multiple instances of such libgromacs, and if so
# we need to find MPI and have it agree with that used for mpi4py.
if (NOT TARGET MPI::MPI_CXX)
    if (HAS_GROMACS_MPI)
        message(WARNING "Inconsistency found! GROMACS does not have MPI support but "
            "find_package(MPI) has been called and created the MPI::MPI_CXX target")
    endif()
    find_package(MPI COMPONENTS CXX)
endif()

# Detect mpi4py. Try to confirm compatibility of linked libraries.
# For MPI-enabled GROMACS, enable mpi4py.h and bindings for
# a `createContext(const ResourceAssignment& resources)` version.
execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c
        "import mpi4py; print(mpi4py.get_include());"
    RESULT_VARIABLE _mpi4py_result
    OUTPUT_VARIABLE _mpi4py_include_path
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT _mpi4py_result EQUAL 0)
    # If mpi4py headers are not available when the gmxapi bindings are built,
    # it is confusing and hard to document robustly how to rebuild and
    # reinstall the package later. To ensure that all gmxapi installations are
    # compatible with mpi4py and libgromacs, we require mpi4py at build time.
    string(
        CONCAT _msg
        "gmxapi core features includes MPI Python bindings and mpi4py compatibility. "
        "Please install mpi4py with a compatible toolchain, then try again to install this package.\n"
    )
    if (MPI_CXX_COMPILER)
        string(
            CONCAT _msg
            ${_msg}
            "The GROMACS C++ compiler appears to be compatible with ${MPI_CXX_COMPILER} "
            "so you might try:\n"
            "  MPICC=${MPI_C_COMPILER} pip install --upgrade --no-cache-dir mpi4py"
        )
    endif ()
    message(FATAL_ERROR ${_msg})
else()
    # NOTE: Even if GROMACS build does not support MPI, we can still accept offers
    # of communicators. In a more advanced case, we could inspect the cgroups or
    # hwloc details to determine how many cores to make available for multithreading,
    # even though we won't be using the MPI facilities of the offered communicator.
    # However, we don't want to refuse to compile, completely, if mpi4py was built
    # with a toolchain that is incompatible with the GROMACS toolchain.
    execute_process(
        COMMAND "${Python3_EXECUTABLE}" -c
            "import mpi4py; print(mpi4py.get_config()['mpicc'])"
        RESULT_VARIABLE _mpi4py_result
        OUTPUT_VARIABLE _mpi4py_mpicc
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    execute_process(
        COMMAND "${Python3_EXECUTABLE}" -c
        "import mpi4py; print(mpi4py.get_config()['mpicxx'])"
        RESULT_VARIABLE _mpi4py_result
        OUTPUT_VARIABLE _mpi4py_mpicxx
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if ((NOT _mpi4py_mpicc OR NOT EXISTS "${_mpi4py_mpicc}") AND (NOT _mpi4py_mpicxx OR NOT EXISTS "${_mpi4py_mpicxx}"))
        execute_process(
            COMMAND "${Python3_EXECUTABLE}" -c
            "import mpi4py; print(mpi4py.get_config()['mpi_dir'])"
            RESULT_VARIABLE _mpi4py_result
            OUTPUT_VARIABLE _mpi4py_mpi_dir
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    endif ()
    # NOTE: The reported compiler wrappers may be falsely reported or unavailable
    # if mpi4py was built in a different environment. We will rely on checks in
    # FindMPI.cmake to confirm tool chain compatibility, later.

    if (GMXAPI_MAIN_PROJECT)
        if (NOT MPI_C_COMPILER AND NOT MPI_CXX_COMPILER)
            if (_mpi4py_mpicc)
                set(MPI_C_COMPILER ${_mpi4py_mpicc})
            endif ()

            if (_mpi4py_mpicxx)
                set(MPI_CXX_COMPILER ${_mpi4py_mpicxx})
            endif ()

            if (_mpi4py_mpi_dir)
                if (NOT MPI_HOME)
                    set(MPI_HOME ${_mpi4py_mpi_dir})
                endif ()
                if (NOT MPI_ROOT)
                    set(MPI_ROOT ${_mpi4py_mpi_dir})
                endif ()
            endif ()
        endif ()
        find_package(MPI COMPONENTS CXX)
    else()
        # Building as part of GROMACS build tree.
        if (_mpi4py_mpicxx AND MPI_CXX_COMPILER AND NOT _mpi4py_compiler_warned)
            if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.19")
                file(REAL_PATH ${_mpi4py_mpicxx} _mpi4py_mpicxx)
                file(REAL_PATH ${MPI_CXX_COMPILER} MPI_CXX_COMPILER)
            endif ()

            set(_warning "mpi4py reports building with ${_mpi4py_mpicxx}, but MPI_CXX_COMPILER=${MPI_CXX_COMPILER}")
            if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24")
                if (NOT ("${_mpi4py_mpicxx}" PATH_EQUAL "${MPI_CXX_COMPILER}"))
                    message(WARNING "${_warning}")
                    set(_mpi4py_compiler_warned TRUE CACHE BOOL "Suppress warning on subsequent runs.")
                endif ()
            else()
                if (NOT ("${_mpi4py_mpicxx}" STREQUAL "${MPI_CXX_COMPILER}"))
                    message(WARNING "${_warning}")
                    set(_mpi4py_compiler_warned TRUE CACHE BOOL "Suppress warning on subsequent runs.")
                endif ()
            endif ()
            unset(_warning)
        endif ()
    endif ()

    # TODO: Try to confirm `mpi4py.MPI.Get_library_version()` and `mpi4py.MPI.Get_library_version()`
    #       is consistent with the version in use by GROMACS. Sample output of Get_library_version():
    #       'Open MPI v4.1.4, package: Open MPI brew@Monterey Distribution, ident: 4.1.4, repo rev: v4.1.4, May 26, 2022\x00'
    # Note that different implementations have different C API calls to support this query.

    set(HAS_MPI4PY TRUE)
endif()

target_link_libraries(_gmxapi PRIVATE Gromacs::gmxapi)

# Enable MPI bindings, to the extent allowed by the Python and GROMACS environments.
if (HAS_MPI4PY AND TARGET MPI::MPI_CXX)
    target_include_directories(_gmxapi PRIVATE ${_mpi4py_include_path})
    target_sources(_gmxapi PRIVATE src/cpp/mpi_bindings.cpp)
    target_sources(_gmxapi PRIVATE src/cpp/pycontext_create.cpp)
    # Conditionally activate the offer_comm() API.
    if (HAS_GROMACS_MPI)
        target_sources(_gmxapi PRIVATE src/cpp/mpi_gromacs_support.cpp)
    else()
        target_sources(_gmxapi PRIVATE src/cpp/mpi_no_gromacs_support.cpp)
    endif ()
    target_link_libraries(_gmxapi PRIVATE MPI::MPI_CXX)
else()
    target_sources(_gmxapi PRIVATE src/cpp/pycontext_create_no_mpi.cpp)
endif()

if (GMXAPI_MAIN_PROJECT)
    set_target_properties(_gmxapi PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)
    set_target_properties(_gmxapi PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
    # Install phase does not affect sdist. Only wheel.
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/src/gmxapi/gmxconfig.json DESTINATION ${SKBUILD_PROJECT_NAME})
    install(TARGETS _gmxapi DESTINATION ${SKBUILD_PLATLIB_DIR}/${SKBUILD_PROJECT_NAME})
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/gmxapi/ DESTINATION ${SKBUILD_PLATLIB_DIR}/${SKBUILD_PROJECT_NAME}
            FILES_MATCHING PATTERN "*.py")
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
    file(COPY ${CMAKE_CURRENT_BINARY_DIR}/src/gmxapi/gmxconfig.json DESTINATION ${GMXAPI_PYTHON_STAGING_DIR}/gmxapi)

    # Unit test and build docs using PYTHONPATH=$CMAKE_CURRENT_BINARY_DIR/gmxapi_staging
    set_target_properties(_gmxapi PROPERTIES staging_dir ${GMXAPI_PYTHON_STAGING_DIR})
    # Note: Integration testing for multiple Python versions and/or CMake-driven
    # sdist preparation could be performed with CMake custom_commands and custom_targets.
endif()

# When building as part of GROMACS umbrella project, add a testing target
# to the `check` target. Normal usage is to first install the Python package,
# then run `pytest` on the `test` directory. Refer to gmxapi package documentation.
if(NOT GMXAPI_MAIN_PROJECT)
    if (BUILD_TESTING)
        add_subdirectory(test)
    endif()
endif()
