#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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

# Manage the OpenMP setup. This wrapper file checks for some known-bad compiler
# versions before trying to detect OpenMP with the standard find-package-module,
# and then does some additional tests for flags afterwards.

# First exclude compilers known to not work with OpenMP although claim to support it:
# gcc 4.2.1 and gcc-llvm 4.2.1 (also claims to be 4.2.1) on Mac OS X
# This fixes redmine 900 and needs to run before OpenMP flags are set below.
if(GMX_OPENMP)
    if (CMAKE_SYSTEM_NAME STREQUAL "Darwin" AND
        (CMAKE_COMPILER_IS_GNUCC AND CMAKE_C_COMPILER_VERSION AND CMAKE_C_COMPILER_VERSION VERSION_LESS 4.3))
        message(STATUS "OpenMP multithreading not supported with gcc/llvm-gcc 4.2 on Mac OS X, disabled")
        set(GMX_OPENMP OFF CACHE BOOL
            "OpenMP multithreading not not supported with gcc/llvm-gcc 4.2 on Mac OS X, disabled!" FORCE)
    elseif(CMAKE_C_COMPILER_ID MATCHES "Cray" AND CMAKE_VERSION VERSION_LESS 3)
        message(STATUS "OpenMP multithreading is not detected correctly for the Cray compiler with CMake before version 3.0 (see http://public.kitware.com/Bug/view.php?id=14567)")
        set(GMX_OPENMP OFF CACHE BOOL
            "OpenMP multithreading is not detected correctly for the Cray compiler with CMake before version 3.0 (see http://public.kitware.com/Bug/view.php?id=14567)" FORCE)
    else()
        # We should do OpenMP detection if we get here
        # OpenMP check must come before other CFLAGS!
        find_package(OpenMP)
        if(OPENMP_FOUND)
            # CMake on Windows doesn't support linker flags passed to target_link_libraries
            # (i.e. it treats /openmp as \openmp library file). Also, no OpenMP linker flags are needed.
            if(NOT (WIN32 AND NOT MINGW))
                if(CMAKE_COMPILER_IS_GNUCC AND GMX_PREFER_STATIC_OPENMP AND NOT APPLE)
                    set(OpenMP_LINKER_FLAGS "-Wl,-static -lgomp -lrt -Wl,-Bdynamic -lpthread")
                    set(OpenMP_SHARED_LINKER_FLAGS "")
                else()
                    # Only set a linker flag if the user didn't set them manually
                    if(NOT DEFINED OpenMP_LINKER_FLAGS)
                        set(OpenMP_LINKER_FLAGS "${OpenMP_C_FLAGS}")
                    endif()
                    if(NOT DEFINED OpenMP_SHARED_LINKER_FLAGS)
                        set(OpenMP_SHARED_LINKER_FLAGS "${OpenMP_C_FLAGS}")
                    endif()
                endif()
            endif()
            if(MINGW)
                #GCC Bug 48659
                set(OpenMP_C_FLAGS "${OpenMP_C_FLAGS} -mstackrealign")
            endif()
        else()
            message(WARNING
                    "The compiler you are using does not support OpenMP parallelism. This might hurt your performance a lot, in particular with GPUs. Try using a more recent version, or a different compiler. For now, we are proceeding by turning off OpenMP.")
            set(GMX_OPENMP OFF CACHE STRING "Whether GROMACS will use OpenMP parallelism." FORCE)
        endif()
    endif()
endif()
gmx_dependent_cache_variable(GMX_OPENMP_MAX_THREADS
    "Maximum number of OpenMP Threads supported. Has to be 32 or a multiple of 64."
    STRING 32 GMX_OPENMP)
mark_as_advanced(GMX_OPENMP_MAX_THREADS)
math(EXPR MAX_THREAD_MOD "${GMX_OPENMP_MAX_THREADS} % 64")
if (NOT GMX_OPENMP_MAX_THREADS EQUAL 32 AND NOT ${MAX_THREAD_MOD} EQUAL 0)
    message(FATAL_ERROR "Only 32 or multiples of 64 supported for GMX_OPENMP_MAX_THREADS.")
endif()
