#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2021- The GROMACS Authors
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

if(GMX_CP2K)
    
    # CMake flags for CP2K linking 
    set(CP2K_DIR "" CACHE STRING "Path to the directory with libcp2k.a library")
    set(CP2K_LINKER_FLAGS "" CACHE STRING "List of flags and libraries required for linking libcp2k. Typically this should be combination of LDFLAGS and LIBS variables from ARCH file used to compile CP2K")

    # Check is CP2K_DIR present (this flags is required)
    if (NOT CP2K_DIR)
        message(FATAL_ERROR "To build GROMACS with CP2K Interface CP2K_DIR should be defined")
    endif()

    # if CP2K_LINKER_FLAGS defined then it should be used for linking instead pkg-config
    if (CP2K_LINKER_FLAGS)
        message(STATUS "CP2K_LINKER_FLAGS will be used to link libcp2k")

        # Add directory with libcp2k.h into system include directories
        include_directories(SYSTEM "${CP2K_DIR}/../../../src/start")

        # Add libcp2k and DBCSR for linking 
        list(APPEND GMX_COMMON_LIBRARIES "-Wl,--allow-multiple-definition -L${CP2K_DIR} -lcp2k -L${CP2K_DIR}/exts/dbcsr -ldbcsr")

        # Add User provided libraries
        list(APPEND GMX_COMMON_LIBRARIES ${CP2K_LINKER_FLAGS})
    else()
        # In other case pkg-config should be used to search for libcp2k
        find_package(PkgConfig QUIET)

        # if pkg-config not found then ask user to provide CP2K_LINKER_FLAGS
        if(NOT PKG_CONFIG_FOUND)
            message(FATAL_ERROR "pkg-config not found, define CP2K_LINKER_FLAGS for custom linking of libcp2k")
        endif()

        # Append PKG_CONFIG_PATH_PATH with ${CP2K_DIR}/pkgconfig which should contain libcp2k.pc
        set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${CP2K_DIR}/pkgconfig")
 
        # Search for libcp2k
        pkg_check_modules(LIBCP2K QUIET libcp2k)

        # if libcp2k not found then ask user to provide CP2K_LINKER_FLAGS
        if(NOT LIBCP2K_FOUND)
            message(FATAL_ERROR "pkg-config could not find libcp2k, define CP2K_LINKER_FLAGS for custom linking of libcp2k")
        endif()

        # libcp2k found: add libraries and paths to the respecting GROMACS variables
        message(STATUS "Found libcp2k in ${LIBCP2K_LIBDIR}")
        include_directories(SYSTEM "${LIBCP2K_INCLUDE_DIRS}")
        link_directories(${LIBCP2K_LIBRARY_DIRS})
        list(APPEND GMX_COMMON_LIBRARIES ${LIBCP2K_LIBRARIES})
    endif()
 
    # If we use GNU compilers then also libgfortran should be linked
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        list(APPEND GMX_COMMON_LIBRARIES "gfortran")
    endif()

    # If we use external MPI then Fortran MPI library should be linked
    if (GMX_LIB_MPI)
        # If Fortran MPI library is found then add it to GMX_COMMON_LIBRARIES
        if (MPI_Fortran_FOUND)
            list(APPEND GMX_COMMON_LIBRARIES ${MPI_Fortran_LIBRARIES})
        else()
            message(FATAL_ERROR "Could not find Fortran MPI library which is required for CP2K")
        endif()
    endif()

endif()
