#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2014,2015, by the GROMACS development team, led by
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
# and in extrae-3.1.0 the following combinations are supported:

# seqtrace:        seq code
# mpitrace:        MPI
# ptmpitrace       pthreads + MPI (unsupported combination in Gromacs)
# ompitrace:       OPENMP + MPI
# oclmpitrace:     OpenCL + MPI
# cudampitrace:    CUDA + MPI
# cudaompitrace:   CUDA+OPENMP+MPI
# omptrace:        OpenMP
# pttrace:         pthreads
# cudatrace:       CUDA
# ocltrace:        OpenCL

option(EXTRAE_LIB_MANUAL "Custom Extrae library to use for profiling" OFF)
mark_as_advanced(EXTRAE_LIB_MANUAL)

if (EXTRAE_LIB_MANUAL)
    set (extraelib "${EXTRAE_LIB_MANUAL}")
else()
    set (extraelib "trace")

    # Checks for what features we have enabled in Gromacs
    # The order of the prefix seems to be: [ocl|cuda]pt[omp|mpi|ompi]
    # default sequential code is "seq"

    # first check for MPI variants
    if (GMX_MPI)
        if (GMX_OPENMP)
            if(GMX_GPU)
                if(GMX_USE_OPENCL)
                    message(FATAL_ERROR "Extrae doesn't have currently support for MPI + OPENMP + OPENCL code")
                else()
                    set (extraelib "cudaompi${extraelib}")
                endif()
            else()
                set (extraelib "ompi${extraelib}")
            endif()
        elseif(GMX_GPU)
            if (GMX_USE_OPENCL)
                set (extraelib "oclmpi${extraelib}")
            else()
                set (extraelib "cudampi${extraelib}")
            endif()
        else()
            set (extraelib "mpi${extraelib}")
        endif()

    # check for OpenMP only support
    elseif (GMX_OPENMP AND NOT GMX_GPU)
        if (GMX_THREAD_MPI)
            message(FATAL_ERROR "Extrae doesn't have currently support for OPENMP + PT code")
        else()
            set (extraelib "omp${extraelib}")
        endif()

    # check for pthreads only support
    elseif (GMX_THREAD_MPI)
        if (GMX_OPENMP)
            message(FATAL_ERROR "Extrae doesn't have currently support for PT + OPENMP code")
        else()
            set (extraelib "pt${extraelib}")
        endif()

    # check for cuda/opencl support
    # libcudatrace and libocltrace have the same name whether or not they have OpenMP support
    elseif (GMX_GPU)
        if (GMX_THREAD_MPI)
            message(FATAL_ERROR "Extrae doesn't have currently support for CUDA/OPENCL + PT code")
        elseif(GMX_USE_OPENCL)
            set (extraelib "ocl${extraelib}")
        else ()
            set (extraelib "cuda${extraelib}")
        endif()

    # default sequential code
    else()
      set (extraelib "seq${extraelib}")
    endif()

endif()

gmx_check_if_changed (EXTRAE_LIBNAME_CHANGED extraelib)

if(EXTRAE_LIBNAME_CHANGED)
    message (STATUS, "Extrae support: will use lib${extraelib} for tracing")
    unset (EXTRAE_LIBRARY)
endif()

find_library(EXTRAE_LIBRARY NAMES ${extraelib})

set(EXTRAE_LIBRARIES ${EXTRAE_LIBRARY} )
set(EXTRAE_INCLUDE_DIRS ${EXTRAE_INCLUDE_DIR} )

# handle the QUIETLY and REQUIRED arguments and set EXTRAE_FOUND to TRUE
# if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(EXTRAE  DEFAULT_MSG
                                  EXTRAE_LIBRARY EXTRAE_INCLUDE_DIR)

mark_as_advanced(EXTRAE_INCLUDE_DIR EXTRAE_LIBRARY )
