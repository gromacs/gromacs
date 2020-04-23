#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2020, by the GROMACS development team, led by
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
#
# Note: If necessary, the Python location can be hinted with Python3_ROOT_DIR
# For additional parameters affecting Python installation discovery, see
# https://cmake.org/cmake/help/latest/module/FindPython3.html#hints
if(FIND_PACKAGE_MESSAGE_DETAILS_Python3)
    # Keep quiet on subsequent runs of cmake
    set(Python3_FIND_QUIETLY ON)
    set(PythonInterp_FIND_QUIETLY ON)
endif()
# Older CMake versions might not search for Python newer than 3.7.
set(Python_ADDITIONAL_VERSIONS 3.8)
# We advocate using Python venvs to manage package availability, so by default
# we want to preferentially discover user-space software.
set(Python3_FIND_REGISTRY LAST)
# Make package discovery consistent with Unix behavior and our documented
# suggestions for installing dependencies.
set(CMAKE_FIND_FRAMEWORK LAST)
if(GMX_PYTHON_PACKAGE)
    find_package(Python3 3.6 COMPONENTS Interpreter Development)
    if (NOT Python3_FOUND OR NOT Python3_Development_FOUND)
        message(FATAL_ERROR "Could not locate Python development requirements. \
                Provide appropriate CMake hints or set GMX_PYTHON_PACKAGE=OFF")
    endif ()
else()
    find_package(Python3 3.6 COMPONENTS Interpreter)
endif()
# Other components, such as pybind and googletest, may expect the
# PYTHON_EXECUTABLE variable from pre-3.12 FindPythonInterp.cmake.
if (Python3_Interpreter_FOUND)
    set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE} CACHE FILEPATH "Location hint for Python interpreter.")
endif ()
# We've already generated all of the output we need, even though other subcomponents
# may call find_package(PythonInterp) later on.
set(Python3_FIND_QUIETLY ON)
set(PythonInterp_FIND_QUIETLY ON)
