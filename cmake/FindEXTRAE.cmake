#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2014, by the GROMACS development team, led by
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

# - Try to find EXTRAE
# Once done this will define
#  EXTRAE_FOUND - System has EXTRAE
#  EXTRAE_INCLUDE_DIRS - The EXTRAE include directories
#  EXTRAE_LIBRARIES - The libraries needed to use EXTRAE

find_path(EXTRAE_INCLUDE_DIR extrae_user_events.h)

# EXTRAE libraries have different names depending on the supported features,
# and in extrae-2.5.0 the following combinations are supported:

# seqtrace:        seq code
# mpitrace:        MPI
# omptrace:        OpenMP
# ompitrace:       MPI + OpenMP
# pttrace:         pthreads
# ptmpitrace       pthreads + MPI (unsupported combination in Gromacs)
# cudatrace:       CUDA
# cudampitrace:    CUDA + MPI
# cudaompitrace:   CUDA+OPENMP+MPI (in the dev version)

# TODO: Add support for the following combinations when available in a future release:

# cudaomptrace:    CUDA + OPENMP
# cudapttrace:        CUDA + pthreads
# cudaptmpitrace:    CUDA + pthreads + MPI (unsupported combination in Gromacs)

set (extraelib "trace")

# libs with MPI support
if (GMX_MPI)
  if (GMX_OPENMP)
    set (extraelib "ompi${extraelib}")
  else()
    set (extraelib "mpi${extraelib}")
  endif()
  if (GMX_GPU)
    set (extraelib "cuda${extraelib}")
  endif()

# other libs with OpenMP support
elseif (GMX_OPENMP)
  set (extraelib "omp${extraelib}")
    if (GMX_GPU)
      set (extraelib "cuda${extraelib}")
    endif()

# library with CUDA only support
elseif (GMX_GPU)
    set (extraelib "cuda${extraelib}")

# library with PThreads support
elseif (GMX_THREAD_MPI)
    set (extraelib "pt${extraelib}")

else()
  set (extraelib "seq${extraelib}")
endif()

find_library(EXTRAE_LIBRARY NAMES  ${extraelib})

set(EXTRAE_LIBRARIES ${EXTRAE_LIBRARY} )
set(EXTRAE_INCLUDE_DIRS ${EXTRAE_INCLUDE_DIR} )

# handle the QUIETLY and REQUIRED arguments and set EXTRAE_FOUND to TRUE
# if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(EXTRAE  DEFAULT_MSG
                                  EXTRAE_LIBRARY EXTRAE_INCLUDE_DIR)

mark_as_advanced(EXTRAE_INCLUDE_DIR EXTRAE_LIBRARY )

