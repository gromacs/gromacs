#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2015,2018,2019, by the GROMACS development team, led by
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

find_program(SPHINX_EXECUTABLE NAMES sphinx-build sphinx-build2
    HINTS
    $ENV{SPHINX_DIR}
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

set(_find_deps_options)
if (Sphinx_FIND_QUIETLY)
    set(_find_deps_options ERROR_QUIET)
endif()

if (NOT Sphinx_pygments_FOUND)
    # Check if pygments module is available via the Unix error code (ie. 0
    # for success)
    execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
        "import pygments"
        RESULT_VARIABLE _pygments_status
        ${_find_deps_options}
        )
    if (_pygments_status EQUAL 0)
        set(Sphinx_pygments_FOUND TRUE CACHE BOOL "Whether pygments module is available for Sphinx")
    endif()
endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Sphinx
    REQUIRED_VARS SPHINX_EXECUTABLE
    VERSION_VAR SPHINX_EXECUTABLE_VERSION
    HANDLE_COMPONENTS
    )
