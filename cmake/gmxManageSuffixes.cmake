#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2013,2014, by the GROMACS development team, led by
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

function(gmx_create_suffix_options)
    set(_defaults "_d for double, _mpi for MPI, _mpi_d for both")
    option(
        GMX_DEFAULT_SUFFIX
        "Use default suffixes for GROMACS binaries and libs (${_defaults}; rerun cmake after changing to see relevant options)"
        ON)
    gmx_dependent_cache_variable(
        GMX_BINARY_SUFFIX
        "Suffix for GROMACS binaries (default: ${_defaults})."
        STRING ""
        "NOT GMX_DEFAULT_SUFFIX")
    gmx_dependent_cache_variable(
        GMX_LIBS_SUFFIX
        "Suffix for GROMACS libraries (default: ${_defaults})."
        STRING ""
        "NOT GMX_DEFAULT_SUFFIX")
endfunction()

gmx_create_suffix_options()
if (GMX_DEFAULT_SUFFIX)
    gmx_check_if_changed(SUFFIXES_CHANGED GMX_DEFAULT_SUFFIX)
    set(GMX_BINARY_SUFFIX "")
    set(GMX_LIBS_SUFFIX "")
    if (GMX_LIB_MPI)
        set(GMX_BINARY_SUFFIX "_mpi")
        set(GMX_LIBS_SUFFIX "_mpi")
    endif()
    if (GMX_DOUBLE)
        set(GMX_BINARY_SUFFIX "${GMX_BINARY_SUFFIX}_d")
        set(GMX_LIBS_SUFFIX "${GMX_LIBS_SUFFIX}_d")
    endif()
    if (SUFFIXES_CHANGED)
        message(STATUS "Using default binary suffix: \"${GMX_BINARY_SUFFIX}\"")
        message(STATUS "Using default library suffix: \"${GMX_LIBS_SUFFIX}\"")
    endif()
else()
    if ("${GMX_LIBS_SUFFIX}" MATCHES "^_mdrun")
        message(FATAL_ERROR "Library suffix (GMX_LIBS_SUFFIX) cannot start with '_mdrun', as that is reserved for internal use")
    endif()
    gmx_check_if_changed(SUFFIXES_CHANGED
                         GMX_DEFAULT_SUFFIX GMX_BINARY_SUFFIX GMX_LIBS_SUFFIX)
    if (SUFFIXES_CHANGED)
        message(STATUS "Using manually set binary suffix: \"${GMX_BINARY_SUFFIX}\"")
        message(STATUS "Using manually set library suffix: \"${GMX_LIBS_SUFFIX}\"")
    endif()
endif()
unset(SUFFIXES_CHANGED)

if (GMX_BUILD_MDRUN_ONLY)
    set(GMX_LIBS_SUFFIX "_mdrun${GMX_LIBS_SUFFIX}")
endif()
