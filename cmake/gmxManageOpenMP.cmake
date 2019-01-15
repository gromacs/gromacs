#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015,2016,2017,2019, by the GROMACS development team, led by
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

if(GMX_OPENMP)
    # We should do OpenMP detection if we get here
    # OpenMP check must come before other CFLAGS!
    find_package(OpenMP)
    if(NOT OPENMP_FOUND)
        message(WARNING
                "The compiler you are using does not support OpenMP parallelism. This might hurt your performance a lot, in particular with GPUs. Try using a more recent version, or a different compiler. For now, we are proceeding by turning off OpenMP.")
        set(GMX_OPENMP OFF CACHE STRING "Whether GROMACS will use OpenMP parallelism." FORCE)
    endif()
endif()
gmx_dependent_cache_variable(GMX_OPENMP_MAX_THREADS
    "Maximum number of OpenMP Threads supported. Has to be 32 or a multiple of 64."
    STRING 64 GMX_OPENMP)
mark_as_advanced(GMX_OPENMP_MAX_THREADS)
math(EXPR MAX_THREAD_MOD "${GMX_OPENMP_MAX_THREADS} % 64")
if (NOT GMX_OPENMP_MAX_THREADS EQUAL 32 AND NOT ${MAX_THREAD_MOD} EQUAL 0)
    message(FATAL_ERROR "Only 32 or multiples of 64 supported for GMX_OPENMP_MAX_THREADS.")
endif()
