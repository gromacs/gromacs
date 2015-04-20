/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \libinternal \file
 * \brief
 * Wraps <mpi.h> usage in Gromacs.
 *
 * This header wraps the MPI header <mpi.h>, and should be included instead of
 * that one.  It abstracts away the case that depending on compilation
 * settings, MPI routines may be provided by <mpi.h> or by thread-MPI.
 * In the absence of MPI, this header still declares some types for
 * convenience.  It also disables MPI C++ bindings that can cause compilation
 * issues.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_GMXMPI_H
#define GMX_UTILITY_GMXMPI_H

#include "config.h"

/*! \cond */
#ifdef GMX_LIB_MPI
/* MPI C++ binding is deprecated and can cause name conflicts (e.g. stdio/mpi seek) */
#define MPICH_SKIP_MPICXX 1
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
/* Starting with 2.2 MPI_INT64_T is required. Earlier version still might have it.
   In theory MPI_Datatype doesn't have to be a #define, but current available MPI
   implementations (OpenMPI + MPICH (+derivates)) use #define and future versions
   should support 2.2. */
#if (MPI_VERSION == 1 || (MPI_VERSION == 2 && MPI_SUBVERSION < 2)) && !defined MPI_INT64_T
#include <limits.h>
#if LONG_MAX == 9223372036854775807L
#define MPI_INT64_T MPI_LONG
#elif LONG_LONG_MAX == 9223372036854775807L
#define MPI_INT64_T MPI_LONG_LONG
#else
#error No MPI_INT64_T and no 64 bit integer found.
#endif
#endif /*MPI_INT64_T*/
#else
#ifdef GMX_THREAD_MPI
#include "thread_mpi/mpi_bindings.h" /* IWYU pragma: export */
#include "thread_mpi/tmpi.h"         /* IWYU pragma: export */
#else
typedef void* MPI_Comm;
typedef void* MPI_Request;
typedef void* MPI_Group;
#define MPI_COMM_NULL NULL
#endif
#endif
//! \endcond

#endif
