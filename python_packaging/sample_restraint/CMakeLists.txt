cmake_minimum_required(VERSION 3.16.3)
# If you are using this repository as a template, you should probably change the
# project name and adopt your own versioning scheme.
project(sample_restraint)

# This project requires a GROMACS supporting gmxapi 0.0.8 or higher. It should
# be sufficient to source the GMXRC, but you can also set the GROMACS_DIR or
# gmxapi_DIR environment variable to help CMake find the GROMACS installation.

# Note that the code will need to be built separately for different versions of Python and for substantially different
# versions of GROMACS. If building from the command line, you can specify a Python executable with the PYTHON_EXECUTABLE
# variable. For instance, to make sure you are building for your default Python, cmake -DPYTHON_EXECUTABLE=`which python`.

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)

# CMake modules are in a subdirectory to keep this file cleaner
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Check if Python package is being built directly or via add_subdirectory.
# I.e. is this being built as a standalone project or as part of the GROMACS
# build tree (for testing)?
set(GMXAPI_EXTENSION_MASTER_PROJECT OFF)
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(GMXAPI_EXTENSION_MASTER_PROJECT ON)
endif()

option(GMXAPI_EXTENSION_DOWNLOAD_PYBIND ON)
if(GMXAPI_EXTENSION_DOWNLOAD_PYBIND)
    configure_file(CMakeLists.pybind.in pybind-download/CMakeLists.txt)
    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
                    RESULT_VARIABLE result
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/pybind-download)
    if(result)
        message(FATAL_ERROR "CMake step for pybind download failed: ${result}")
    endif()
    execute_process(COMMAND ${CMAKE_COMMAND} --build .
                    RESULT_VARIABLE result
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/pybind-download)
    if(result)
        message(FATAL_ERROR "Build step for pybind failed: ${result}")
    endif()
    add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/pybind-src)
else()
    find_package(pybind11 2.6 REQUIRED)
endif()
# Assuming GROMACS is in our path or that we have set either the gmxapi_DIR or GROMACS_DIR environment variables,
# this will find the CMake configuration for the GROMACS libraries we need and define the CMake library objects
# Gromacs::gmxapi
if(GMXAPI_EXTENSION_MASTER_PROJECT)
    find_package(gmxapi
                 0.2.1 REQUIRED CONFIG
                 PATHS "$ENV{GROMACS_DIR}"
                 )

    message(STATUS "Found gmxapi version ${gmxapi_VERSION_MAJOR}.${gmxapi_VERSION_MINOR}.${gmxapi_VERSION_PATCH}")
endif()

if(GMXAPI_EXTENSION_MASTER_PROJECT)
    ######################################################
    # The following is boiler-plate recommended by GROMACS
    ######################################################
    # In principle, this could be deduced from GROMACS_IS_DOUBLE returned by
    # find_package(GROMACS) based on the suffix alone, but it is clearer that the
    # user explicitly sets what they want to get, and then need to provide a suffix
    # to match.
    option(GMX_DOUBLE "Use double precision" OFF)
    set(GMX_SUFFIX "" CACHE STRING "Suffix for the GROMACS installation to use (empty for default)")

    # This does not allow for a non-suffixed double-precision libgromacs, but
    # that should be rare enough for demonstration purposes.
    if (GMX_DOUBLE AND NOT GMX_SUFFIX)
        set(GROMACS_SUFFIX "_d")
    else()
        set(GROMACS_SUFFIX ${GMX_SUFFIX})
    endif()

    find_package(GROMACS REQUIRED)
    gromacs_check_double(GMX_DOUBLE)
    gromacs_check_compiler(CXX)
    include_directories(${GROMACS_INCLUDE_DIRS})
    add_definitions(${GROMACS_DEFINITIONS})

    # Use static linking on MSVC
    if (CMAKE_GENERATOR MATCHES "Visual Studio")
        string(REPLACE /MD /MT CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE})
        set(CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE} CACHE STRING "" FORCE)
        string(REPLACE /MD /MT CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG})
        set(CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG} CACHE STRING "" FORCE)
    endif()

    # Activate the CTest module for this directory.
    include(CTest)
endif()

########################################################

# Stuff for our plugin:
#
# If the user is not in a virtual environment and is not a privileged user and has not specified an install location
# for the Python module (GMXPLUGIN_INSTALL_PATH), this option causes the automatic install location to query the user
# site-packages directory instead of using the default site-packages directory for the interpreter.
option(GMXPLUGIN_USER_INSTALL
       "Override the default site-packages directory with the user-specific Python packages directory. \
       (Do not use with virtual environments.) \
       Has no effect if GMXPLUGIN_INSTALL_PATH is defined or cached. \
       Use -UGMXPLUGIN_INSTALL_PATH to force recalculation."
       OFF)

# Since a user may have multiple virtual environments with different Python interpreters, it is generally confusing to
# have a package for a virtual environment installed in the user's default user site-packages directory. If virtual
# environments are in use at all, we recommend you do _not_ perform a "user" install in or out of a virtual env. If you do
# not use any Python virtual environments, we recommend you _do_ perform "user" installs exclusively. Overall, we
# we recommend you use Python virtual environments and activate one before performing a regular (non-"user") install.
if (PYTHON_EXECUTABLE)
    message(STATUS "Found Python interpreter: ${PYTHON_EXECUTABLE}")
    if (PYTHON_LIBRARIES OR PythonLibs_FOUND OR pybind11_FOUND)
        if (GMXPLUGIN_USER_INSTALL)
            execute_process(COMMAND ${PYTHON_EXECUTABLE} "-m" "site" "--user-site"
                            OUTPUT_VARIABLE GMXPLUGIN_DEFAULT_SITE_PACKAGES
                            OUTPUT_STRIP_TRAILING_WHITESPACE)
            message(STATUS "Python user site-packages directory is ${GMXPLUGIN_DEFAULT_SITE_PACKAGES}")
        else()
            execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
                                "import sys; import os; \
                                print(os.path.abspath(os.path.join(sys.prefix, \
                                    'lib', \
                                    'python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}', \
                                    'site-packages')))"
                            OUTPUT_VARIABLE GMXPLUGIN_DEFAULT_SITE_PACKAGES
                            OUTPUT_STRIP_TRAILING_WHITESPACE)
            message(STATUS "Python site-packages directory is ${GMXPLUGIN_DEFAULT_SITE_PACKAGES}")
        endif()
    else()
        message(FATAL_ERROR
            "Found Python interpreter ${PYTHON_EXECUTABLE} but this Python installation does not have developer tools."
            "Set PYTHON_EXECUTABLE to the Python interpreter that was installed with a working Python.h header file.")
    endif()
else()
    message(FATAL_ERROR "Could not find Python interpreter. Set CMake flag -DPYTHON_EXECUTABLE=/path/to/python to hint.")
endif()

# At some point this may be part of a CMake package with several components for which a single CMAKE_INSTALL_PREFIX does
# not make sense, so let's manage the install path separately.
set(GMXPLUGIN_INSTALL_PATH ${GMXPLUGIN_DEFAULT_SITE_PACKAGES} CACHE PATH
    "Path to Python module install location (site-packages). For an automatically determined install location based on \
    the Python installation, leave undefined or explicitly undefined with -UGMXPLUGIN_INSTALL_PATH and, optionally, set \
    GMXPLUGIN_USER_INSTALL on or off to specify the installation's site-packages directory or the 'user' site-packages \
    directory.")

if(GMXAPI_EXTENSION_MASTER_PROJECT)
    message(STATUS "Python module will be installed to GMXPLUGIN_INSTALL_PATH cache value ${GMXPLUGIN_INSTALL_PATH}")
endif()

# TODO: (Issue #3027) Handle Googletest sources both for forked projects and
#  for GROMACS-project-internal testing.
# Projects based on this subtree should bundle googletest sources.
# The GROMACS project already bundles googletest sources for internal use, but
# they will only be available to this project with some additional management by
# the parent project CMake configuration.
option(DOWNLOAD_GOOGLETEST "Download the latest master branch of googletest." OFF)
mark_as_advanced(DOWNLOAD_GOOGLETEST)

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)


# Now move on to building the custom code.
add_subdirectory(src)

# Set up documentation build targets (work in progress).
add_subdirectory(docs)

# Process CMake configuration for Python and C++ tests.
include(CTest)
if(BUILD_TESTING)
    include(GoogleTest)
    add_subdirectory(tests)
endif()
