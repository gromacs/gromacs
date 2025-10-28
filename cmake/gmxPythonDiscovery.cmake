#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2020- The GROMACS Authors
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

# Perform Python installation discovery early and in one place, for consistency.

# With embedded packages in this repository, the module-scoped default behavior of
# FindPython3 version/component requirements and artifact specification is likely
# confusing (at best) or leading to bad logic (at worst).
# https://cmake.org/cmake/help/latest/module/FindPython3.html#artifacts-specification
option(Python3_ARTIFACTS_INTERACTIVE TRUE
       "Make artifacts specification global and cached.")

# Note: If necessary, the Python location can be hinted with Python3_ROOT_DIR
# For additional parameters affecting Python installation discovery, see
# https://cmake.org/cmake/help/latest/module/FindPython3.html#hints
if (Python3_FIND_QUIETLY_AFTER_FIRST_RUN)
    # Keep quiet on subsequent runs of cmake
    set (Python3_FIND_QUIETLY TRUE)
else()
    set (Python3_FIND_QUIETLY FALSE)
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
find_package(Python3 3.9 COMPONENTS Interpreter Development)
if (GMX_PYTHON_PACKAGE AND (NOT Python3_FOUND OR NOT Python3_Development_FOUND))
    message(FATAL_ERROR "Could not locate Python development requirements. \
            Provide appropriate CMake hints or set GMX_PYTHON_PACKAGE=OFF")
endif ()
set(Python3_FIND_QUIETLY_AFTER_FIRST_RUN TRUE CACHE INTERNAL "Be quiet during future attempts to find Python3")

# Provide hints for other Python detection that may occur later.
#
# Other components, such as pybind, googletest, and our python_packaging
# expect the PYTHON_EXECUTABLE variable from pre-3.12 FindPythonInterp.cmake.
if (Python3_Interpreter_FOUND)
    set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE} CACHE FILEPATH "Location hint for Python interpreter.")
endif ()
# We've already generated all of the output we need, even though other subcomponents
# may call find_package(PythonInterp) later on.
set(Python3_FIND_QUIETLY ON)
set(PythonInterp_FIND_QUIETLY ON)

# The standard FindPython3.cmake is not fully reproducible when run a second time in the
# same build tree, so we work around that here
unset(_Python3_Interpreter_REASON_FAILURE CACHE)
