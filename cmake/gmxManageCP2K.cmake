#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2021, by the GROMACS development team, led by
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

if(GMX_CP2K)
    
    # Check that all necessary CMake flags are present 
    set(CP2K_DIR "" CACHE STRING "Path to the directory with libcp2k.a library")
    set(CP2K_LINKER_FLAGS "" CACHE STRING "List of flags and libraries required for linking libcp2k. Typically this should be combination of LDFLAGS and LIBS variables from ARCH file used to compile CP2K")
    if ((CP2K_DIR STREQUAL "") OR (CP2K_LINKER_FLAGS STREQUAL ""))
        message(FATAL_ERROR "To build GROMACS with CP2K Interface both CP2K_DIR and CP2K_LINKER_FLAGS should be defined")
    endif()

    # Add directory with libcp2k.h into system include directories
    include_directories(SYSTEM "${CP2K_DIR}/../../../src/start")

    # Add libcp2k and DBCSR for linking 
    set(CMAKE_CXX_STANDARD_LIBRARIES "${CMAKE_CXX_STANDARD_LIBRARIES} -Wl,--allow-multiple-definition -L${CP2K_DIR} -lcp2k -L${CP2K_DIR}/exts/dbcsr -ldbcsr")

    # Add User provided libraries
    set(CMAKE_CXX_STANDARD_LIBRARIES "${CMAKE_CXX_STANDARD_LIBRARIES} ${CP2K_LINKER_FLAGS}")

    # If we use GNU compilers then also libgfortran should be linked
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(CMAKE_CXX_STANDARD_LIBRARIES "${CMAKE_CXX_STANDARD_LIBRARIES} -lgfortran")
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

