cmake_minimum_required(VERSION 3.13.0)
# If you are using this repository as a template, you should probably change the
# project name and adopt your own versioning scheme.
project(sample_restraint VERSION 0.2.1)

find_package(PythonInterp)

# Allow gmxapi_ROOT hint (available with CMake 3.12).
cmake_policy(SET CMP0074 NEW)

# Check if Python package is being built directly or via add_subdirectory.
# I.e. is this being built as a standalone project or as part of the GROMACS
# build tree (for testing)?
set(GMXAPI_EXTENSION_MASTER_PROJECT OFF)
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(GMXAPI_EXTENSION_MASTER_PROJECT ON)
endif()

# TODO: (Issue #3027) Handle pybind sources both for forked projects and
#  for GROMACS-project-internal testing.
# Projects based on this subtree should bundle pybind11 sources.
# The GROMACS project already bundles one copy of pybind sources for the gmxapi
# Python package, and should package them with source distribution archives of
# the paackage that should be installed with GROMACS, but the pybind sources
# would only be available to this project with some additional management by
# the parent project CMake configuration.
# Also reference https://gitlab.com/gromacs/gromacs/-/issues/2896
if(GMXAPI_EXTENSION_MASTER_PROJECT)
    # TODO: (Issue #3027) Handle locally available sources.
    set(GMXAPI_EXTENSION_USE_BUNDLED_PYBIND OFF CACHE BOOL
        "Use pybind11 headers bundled with this repository. If OFF, CMake does `find_package(pybind11)`.")
    if(GMXAPI_EXTENSION_USE_BUNDLED_PYBIND)
        add_subdirectory(external/pybind)
    else()
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
            find_package(pybind11 2.2 REQUIRED)
        endif()
    endif()
else()
    # If configuring as part of the GROMACS project, pybind11 2.2 has already been found.
endif()


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

# Assuming GROMACS is in our path or that we have set either the gmxapi_DIR or GROMACS_DIR environment variables,
# this will find the CMake configuration for the GROMACS libraries we need and define the CMake library objects
# Gromacs::gmxapi
if(GMXAPI_EXTENSION_MASTER_PROJECT)
    find_package(gmxapi
                 0.0.8 REQUIRED CONFIG
                 PATHS "$ENV{GROMACS_DIR}"
                 )

    message(STATUS "Found gmxapi version ${gmxapi_VERSION_MAJOR}.${gmxapi_VERSION_MINOR}.${gmxapi_VERSION_PATCH}")
endif()

if(GMXAPI_EXTENSION_MASTER_PROJECT)
    ######################################################
    # The following is boiler-plate recommended by GROMACS
    ######################################################

    find_package(GROMACS REQUIRED
                 NAMES gromacs gromacs_mpi gromacs_d
                 HINTS "$ENV{GROMACS_DIR}"
                 )

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

if (PYTHONINTERP_FOUND)
    message(STATUS "Found Python interpreter: ${PYTHON_EXECUTABLE}")
    if (PYTHON_LIBRARIES)
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

# Projects based on this subtree should bundle googletest sources.
# The GROMACS project already bundles googletest sources for internal use, but
# they will only be available to this project with some additional management by
# the parent project CMake configuration.
option(DOWNLOAD_GOOGLETEST "Download the version of googletest used for GROMACS testing." OFF)
mark_as_advanced(DOWNLOAD_GOOGLETEST)


# Now move on to building the custom code.
add_subdirectory(src)

# Set up documentation build targets (work in progress).
add_subdirectory(docs)

if(DOWNLOAD_GOOGLETEST)
    # Download and unpack googletest at configure time. For simplicity, we pin to the same
    # googletest ref as in
    #    https://gitlab.com/gromacs/gromacs/-/blob/master/src/external/googletest/README.Gromacs
    # Reference https://google.github.io/googletest/quickstart-cmake.html#set-up-a-project
    include(FetchContent)
    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/release-1.8.1.zip
    )
    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    # Note: as of CMake 3.16.3, setting policy CMP0077 to NEW does not seem to allow
    # `set(INSTALL_GTEST OFF)` to work as desired, so we have to override the option definition.
    option(INSTALL_GTEST "Override default behavior" OFF)
    FetchContent_MakeAvailable(googletest)

    add_library(GTest::GTest ALIAS gtest)
    add_library(GTest::Main ALIAS gtest_main)
else()
    if(GMXAPI_EXTENSION_MASTER_PROJECT)
        # Allow project maintainers to bundle or provide googletest sources in
        # a subdirectory `external`.
        if (DEFINED ENV{CI})
            # In a CI environment, assume that the supporting software can be provided.
            if ("$ENV{CI_PROJECT_NAME}" STREQUAL "gromacs")
                # We are running in a GitLab CI environment for a gromacs repository.
                # See also admin/gitlab-ci/sample_restraint.gitlab-ci.yml
                FILE(COPY
                     $ENV{CI_PROJECT_DIR}/src/external/googletest
                     DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/external/
                     )
            endif ()
            # Note: as of CMake 3.16.3, setting policy CMP0077 to NEW does not seem to allow
            # `set(INSTALL_GTEST OFF)` to work as desired, so we have to override the option definition.
            option(INSTALL_GTEST "Override default behavior" OFF)
            add_subdirectory(external/googletest)
            add_library(GTest::GTest ALIAS gtest)
            add_library(GTest::Main ALIAS gtest_main)
        else()
            # Logic could be added here for other CI environments. Otherwise, the user is
            # responsible for helping to find GTest or copying sources to external/googletest
            find_package(GTest REQUIRED)
        endif ()
    else() # building as part of a GROMACS build...
        if(TARGET gtest)
            add_library(GTest::GTest ALIAS gtest)
            add_library(GTest::Main ALIAS gtest_main)
        else()
            find_package(GTest)
        endif()
        if(NOT TARGET GTest::GTest)
            message(FATAL_ERROR "GTest::GTest target should have been supplied by the parent project.")
        endif()
        if(NOT TARGET GTest::Main)
            message(FATAL_ERROR "GTest::Main target should have been supplied by the parent project.")
        endif()
    endif()
endif() # DOWNLOAD_GOOGLETEST

# Process CMake configuration for Python and C++ tests.
include(CTest)
if(BUILD_TESTING)
    include(GoogleTest)
    add_subdirectory(tests)
endif()
