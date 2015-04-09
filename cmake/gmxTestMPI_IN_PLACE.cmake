#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2009,2011,2012,2014,2015, by the GROMACS development team, led by
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

# - Define macro to check if MPI_IN_PLACE exists
#
#  GMX_TEST_MPI_IN_PLACE(VARIABLE)
#
#  VARIABLE will be set to true if MPI_IN_PLACE exists
#

include(CheckCSourceCompiles)
MACRO(GMX_TEST_MPI_IN_PLACE VARIABLE)
  if(NOT DEFINED MPI_IN_PLACE_COMPILE_OK)
    MESSAGE(STATUS "Checking for MPI_IN_PLACE")

    set(CMAKE_REQUIRED_FLAGS ${MPI_COMPILE_FLAGS})
    set(CMAKE_REQUIRED_INCLUDES ${MPI_INCLUDE_PATH})
    set(CMAKE_REQUIRED_LIBRARIES ${MPI_LIBRARIES})
    check_c_source_compiles(
      "#include <mpi.h>
int main(void) {
  void* buf;
  MPI_Allreduce(MPI_IN_PLACE, buf, 10, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
}" MPI_IN_PLACE_COMPILE_OK)

    if(MPI_IN_PLACE_COMPILE_OK)
        MESSAGE(STATUS "Checking for MPI_IN_PLACE - yes")
    else()
        MESSAGE(STATUS "Checking for MPI_IN_PLACE - no")
    endif()
    set(MPI_IN_PLACE_COMPILE_OK "${MPI_IN_PLACE_COMPILE_OK}" CACHE INTERNAL "Result of mpi_in_place check")
    set(CMAKE_REQUIRED_FLAGS)
    set(CMAKE_REQUIRED_INCLUDES)
    set(CMAKE_REQUIRED_LIBRARIES)
  endif()
  if (MPI_IN_PLACE_COMPILE_OK)
    set(${VARIABLE} ${MPI_IN_PLACE_COMPILE_OK}
      "Result of test for MPI_IN_PLACE")
  endif()
ENDMACRO(GMX_TEST_MPI_IN_PLACE VARIABLE)



