#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2017, by the GROMACS development team, led by
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

set(GMX_API_MAJOR 0)
set(GMX_API_MINOR 0)
set(GMX_API_RELEASE 0)

###############################
# New public C++ API components
###############################
add_subdirectory(cpp)

###############
# Python module
###############
option(GMX_PYTHON_MODULE "Build Gromacs Python module" ON)
if (GMX_PYTHON_MODULE)
    find_package(PythonInterp REQUIRED)
    find_package(PythonLibs REQUIRED)
    # Double-check we found Python.h...
    MESSAGE(STATUS "Checking for headers required for Python C++ extensions.")
    set(CMAKE_REQUIRED_INCLUDES "${PYTHON_INCLUDE_DIRS}")
    MESSAGE(STATUS "Looking in ${PYTHON_INCLUDE_DIRS}")
    CHECK_INCLUDE_FILE(Python.h FOUND_PYTHON_H)

    check_cxx_source_compiles(
      "#include <Python.h>
    int main(void) {
    return 0;
    }" PYTHON_H)

    set(CMAKE_REQUIRED_INCLUDES)
    IF(${PYTHON_H})
    ELSE()
      MESSAGE(FATAL_ERROR "No Python.h?")
    ENDIF()

    add_subdirectory(python)
endif()

##############################################
# Sphinx documentation for external interfaces
##############################################
# TODO: enforce sphinx associated with consistent Python library version.
# E.g. SPHINX_EXECUTABLE
unset(Sphinx_FOUND)
find_package(Sphinx 1.3 COMPONENTS pygments)
if(Sphinx_FOUND)
    add_subdirectory(docs)
endif()
