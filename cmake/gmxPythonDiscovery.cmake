#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2020,2021, by the GROMACS development team, led by
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

# Perform Python installation discovery early and in one place, for consistency.

if($<VERSION_GREATER_EQUAL:${CMAKE_VERSION},"3.18">)
    # With embedded packages in this repository, the module-scoped default behavior of
    # FindPython3 version/component requirements and artifact specification is likely
    # confusing (at best) or leading to bad logic (at worst).
    # https://cmake.org/cmake/help/latest/module/FindPython3.html#artifacts-specification
    option(Python3_ARTIFACTS_INTERACTIVE TRUE
           "Make artifacts specification global and cached.")
elseif(NOT FIND_PACKAGE_MESSAGE_DETAILS_Python3)
    # On first run, check whether the user has triggered the behavior described above.
    if (Python3_EXECUTABLE
            OR Python3_COMPILER
            OR Python3_DOTNET_LAUNCHER
            OR Python3_LIBRARY
            OR Python3_INCLUDE_DIR
            OR Python3_NumPy_INCLUDE_DIR)
        message(WARNING
                "Specifying FindPython3 artifacts may complicate consistent Python detection "
                "in projects bundled with GROMACS. See "
                "https://cmake.org/cmake/help/latest/module/FindPython3.html#artifacts-specification"
                )
    endif ()
endif()

# Note: If necessary, the Python location can be hinted with Python3_ROOT_DIR
# For additional parameters affecting Python installation discovery, see
# https://cmake.org/cmake/help/latest/module/FindPython3.html#hints
if(FIND_PACKAGE_MESSAGE_DETAILS_Python3)
    # Keep quiet on subsequent runs of cmake
    set(Python3_FIND_QUIETLY ON)
    set(PythonInterp_FIND_QUIETLY ON)
endif()
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
find_package(Python3 3.7 COMPONENTS Interpreter Development)
if (GMX_PYTHON_PACKAGE AND (NOT Python3_FOUND OR NOT Python3_Development_FOUND))
    message(FATAL_ERROR "Could not locate Python development requirements. \
            Provide appropriate CMake hints or set GMX_PYTHON_PACKAGE=OFF")
endif ()

# Provide hints for other Python detection that may occur later.
#
# Other components, such as pybind and googletest, may expect the
# PYTHON_EXECUTABLE variable from pre-3.12 FindPythonInterp.cmake.
if (Python3_Interpreter_FOUND)
    set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE} CACHE FILEPATH "Location hint for Python interpreter.")
endif ()
# We've already generated all of the output we need, even though other subcomponents
# may call find_package(PythonInterp) later on.
set(Python3_FIND_QUIETLY ON)
set(PythonInterp_FIND_QUIETLY ON)
# Older versions of FindPythonLibs and FindPythonInterp might not search for
# Python newer than 3.7 by default. (as of CMake 3.9)
# Note that this hint is not used by the newer FindPython module we rely on above.
set(Python_ADDITIONAL_VERSIONS 3.8 3.9 3.10)
