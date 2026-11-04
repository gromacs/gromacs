#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2015- The GROMACS Authors
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

# Check the Python installation for sphinx-build, even if not on PATH.
if(Python3_EXECUTABLE)
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c "import sys; print(sys.exec_prefix)"
        OUTPUT_VARIABLE _python_exec_prefix
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    list(APPEND CMAKE_PROGRAM_PATH ${_python_exec_prefix})
    unset(_python_exec_prefix)
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -m site --user-base
        OUTPUT_VARIABLE _python_user_base
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if (_python_user_base)
        list(APPEND CMAKE_PROGRAM_PATH ${_python_user_base})
    endif ()
    unset(_python_user_base)
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c "import sys; print(sys.base_exec_prefix)"
        OUTPUT_VARIABLE _python_base_exec_prefix
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if (_python_base_exec_prefix)
        list(APPEND CMAKE_PROGRAM_PATH ${_python_base_exec_prefix})
    endif ()
    unset(_python_base_exec_prefix)
endif()

find_program(SPHINX_EXECUTABLE NAMES sphinx-build sphinx-build2
             HINTS ENV SPHINX_DIR
             PATH_SUFFIXES bin
             DOC "Sphinx documentation generator"
             )
mark_as_advanced(SPHINX_EXECUTABLE)

# Detect Sphinx version
if (SPHINX_EXECUTABLE AND NOT DEFINED SPHINX_EXECUTABLE_VERSION)
    execute_process(
        COMMAND ${SPHINX_EXECUTABLE} --version
        ERROR_VARIABLE  SPHINX_VERSION_OUTPUT_VARIABLE
        OUTPUT_VARIABLE SPHINX_VERSION_OUTPUT_VARIABLE
        RESULT_VARIABLE SPHINX_VERSION_RESULT_VARIABLE
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )

    # Detect the sphinx version. First try to match the error message
    # from old versions that didn't even support --version, then try
    # to detect more modern sphinx versions. If nothing is found, then
    # the cache variable is set to an empty value.
    set(_version "")
    if(SPHINX_VERSION_OUTPUT_VARIABLE MATCHES "Sphinx v([0-9\.]+)\n.*")
        set(_version ${CMAKE_MATCH_1})
    elseif (SPHINX_VERSION_OUTPUT_VARIABLE MATCHES ".*build[ )]*(.*)")
        set(_version ${CMAKE_MATCH_1})
    endif()

    set(SPHINX_EXECUTABLE_VERSION ${_version} CACHE INTERNAL "Version of ${SPHINX_EXECUTABLE}")
endif()


# Check that the Python modules requested as components (Sphinx extensions
# and themes, given by their import names) can be imported.
include(FindPythonModule)
foreach(_module ${Sphinx_FIND_COMPONENTS})
    find_python_module(${_module} QUIET)
    string(TOUPPER ${_module} _module_upper)
    if (PYTHONMODULE_${_module_upper})
        set(Sphinx_${_module}_FOUND TRUE)
    endif()
endforeach()
unset(_module_upper)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Sphinx
    REQUIRED_VARS SPHINX_EXECUTABLE
    VERSION_VAR SPHINX_EXECUTABLE_VERSION
    HANDLE_COMPONENTS
    )
