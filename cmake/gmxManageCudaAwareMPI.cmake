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

# - Define function to check if underlying MPI is CUDA-aware
#
#  GMX_TEST_CUDA_AWARE_MPI()
#
#  GMX_TEST_CUDA_AWARE_MPI puts HAVE_CUDA_AWARE_MPI variable in cache
#
include(CheckCXXSourceCompiles)
function(GMX_TEST_CUDA_AWARE_MPI)
  if (NOT DEFINED HAVE_CUDA_AWARE_MPI)
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
      #if defined(MPIX_CUDA_AWARE_SUPPORT) && (MPIX_CUDA_AWARE_SUPPORT==1)
        return 0;
      #else
      #error MPI implementation is not CUDA-aware
      #endif
      }" HAVE_CUDA_AWARE_MPI)

    if(HAVE_CUDA_AWARE_MPI)
      MESSAGE(STATUS "Checking for CUDA_AWARE_MPI - yes")
    else()
      MESSAGE(STATUS "Checking for CUDA_AWARE_MPI - no")
      MESSAGE(WARNING "GROMACS cannot determine if underlying MPI is CUDA-aware, " 
      "for better multi-GPU performance consider using a more recent CUDA-aware MPI.")
    endif()
  endif()
endfunction()

# Test if CUDA-aware MPI is supported
gmx_test_cuda_aware_mpi()




