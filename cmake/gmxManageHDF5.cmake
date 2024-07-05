#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2023- The GROMACS Authors
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

set(GMX_HDF5_REQUIRED_VERSION "1.10.1")

macro(gmx_manage_hdf5)
    # Find an external hdf5 library.
    find_package(HDF5 ${GMX_HDF5_REQUIRED_VERSION})
    set(GMX_USE_HDF5 ${HDF5_FOUND} CACHE BOOL "Use HDF5 (is it available?)")
    if(GMX_USE_HDF5)
        if(NOT HDF5_FOUND OR HDF5_VERSION VERSION_LESS GMX_HDF5_REQUIRED_VERSION)
            message(FATAL_ERROR "Cannot find HDF5 (version required ${GMX_HDF5_REQUIRED_VERSION}). Disable HDF5 by setting the option GMX_USE_HDF5 to OFF.")
        endif()
    endif()
endmacro()

macro(gmx_manage_sz3)
    # FIXME: H5Z-SZ3 cannot be built using gcc on Mac OS X
    if(GMX_USE_HDF5)
        include(FetchContent)
        set(FETCHCONTENT_QUIET OFF)
        set(BUILD_H5Z_FILTER ON)
        FetchContent_Declare(SZ3 SOURCE_DIR ${CMAKE_SOURCE_DIR}/src/external/SZ3-bio)
        FetchContent_GetProperties(SZ3)
        if (NOT ${sz3_POPULATED})
            FetchContent_Populate(SZ3)
            add_subdirectory(${sz3_SOURCE_DIR} ${sz3_BINARY_DIR} EXCLUDE_FROM_ALL)
        endif()
        if (BUILD_SHARED_LIBS)
            install(TARGETS SZ3 hdf5sz3 EXPORT fileio)
        endif()
        install(FILES ${CMAKE_SOURCE_DIR}/src/gromacs/fileio/sz3.config DESTINATION share/SZ3)
    endif()
endmacro()
