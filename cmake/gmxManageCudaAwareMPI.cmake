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

# - Define function to check if underlying MPI is CUDA-aware
#
#  GMX_TEST_CUDA_AWARE_MPI()
#
#  GMX_TEST_CUDA_AWARE_MPI puts HAVE_MPI_EXT, MPI_SUPPORTS_CUDA_AWARE_DETECTION variables in cache
#
include(CheckCXXSourceCompiles)
function(GMX_TEST_CUDA_AWARE_MPI)
  if (NOT DEFINED HAVE_MPI_EXT OR NOT DEFINED MPI_SUPPORTS_CUDA_AWARE_DETECTION)
    MESSAGE(STATUS "Checking for CUDA_AWARE_MPI")
    list(JOIN MPI_COMPILE_FLAGS " " CMAKE_REQUIRED_FLAGS)
    set(CMAKE_REQUIRED_INCLUDES ${MPI_INCLUDE_PATH})
    set(CMAKE_REQUIRED_LIBRARIES ${MPI_LIBRARIES})
    # cannot use check_include_file here as mpi.h needs to be included
    # before mpi-ext.h for compilation, check_include_file doesn't support 
    # this use-case
    check_cxx_source_compiles(
      "#include <mpi.h>
      #include <mpi-ext.h>
      int main(void) 
      {
        return 0;
      }" HAVE_MPI_EXT)

    if(NOT HAVE_MPI_EXT)
      set(HAVE_MPI_EXT 0)
    endif()

    list(APPEND CMAKE_REQUIRED_DEFINITIONS -DHAVE_MPI_EXT=${HAVE_MPI_EXT})

    check_cxx_source_compiles(
    "#include <mpi.h>
    #if HAVE_MPI_EXT
    #include <mpi-ext.h>
    #endif
    int main(void) 
    {
      return MPIX_Query_cuda_support();
    }" MPI_SUPPORTS_CUDA_AWARE_DETECTION)

    list(REMOVE_ITEM CMAKE_REQUIRED_DEFINITIONS -DHAVE_MPI_EXT)

    if(MPI_SUPPORTS_CUDA_AWARE_DETECTION)
      MESSAGE(STATUS "Checking for MPI_SUPPORTS_CUDA_AWARE_DETECTION - yes")
    else()
      MESSAGE(STATUS "Checking for MPI_SUPPORTS_CUDA_AWARE_DETECTION - no")
      MESSAGE(WARNING "GROMACS cannot determine if underlying MPI is CUDA-aware, " 
      "for better multi-GPU performance consider using a more recent CUDA-aware MPI.")
    endif()
  endif()
endfunction()

# Test if CUDA-aware MPI is supported
gmx_test_cuda_aware_mpi()




