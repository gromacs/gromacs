/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
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

#ifndef GMX_TIMING_WALLTIME_ACCOUNTING_H
#define GMX_TIMING_WALLTIME_ACCOUNTING_H

#include "gromacs/utility/basedefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

#if 0
}
#endif

/*! \brief
 * Contains per-process and per-thread data about elapsed wall-clock
 * times and integration steps performed. */
typedef struct gmx_walltime_accounting *gmx_walltime_accounting_t;

//! Constructor
gmx_walltime_accounting_t
walltime_accounting_init(int numOpenMPThreads);

//! Destructor
void
walltime_accounting_destroy(gmx_walltime_accounting_t walltime_accounting);

/*! \brief
 * Record initial time stamps, e.g. at run end or counter re-initalization time
 */
void
walltime_accounting_start(gmx_walltime_accounting_t walltime_accounting);

/*! \brief
 * Measure and cache the elapsed wall-clock time since
 * walltime_accounting_start() */
void
walltime_accounting_end(gmx_walltime_accounting_t walltime_accounting);

/*! \brief
 * Measure and return the elapsed wall-clock time since
 * walltime_accounting_start() */
double
walltime_accounting_get_current_elapsed_time(gmx_walltime_accounting_t walltime_accounting);

//! Get the cached wall-clock time for this node
double
walltime_accounting_get_elapsed_time(gmx_walltime_accounting_t walltime_accounting);

//! Get the cached wall-clock time, multiplied by the number of OpenMP threads
double
walltime_accounting_get_elapsed_time_over_all_threads(gmx_walltime_accounting_t walltime_accounting);

//! Get the cached initial time stamp for this node
double
walltime_accounting_get_start_time_stamp(gmx_walltime_accounting_t walltime_accounting);

//! Get the number of integration steps done
gmx_int64_t
walltime_accounting_get_nsteps_done(gmx_walltime_accounting_t walltime_accounting);

/*! \brief Set the number of integration steps done
 *
 * TODO consider whether this should get done in walltime_accounting_end */
void
walltime_accounting_set_nsteps_done(gmx_walltime_accounting_t   walltime_accounting,
                                    gmx_int64_t                 nsteps_done);

/*! \brief
 * Calls system timing routines (e.g. clock_gettime) to get the (fractional)
 * number of seconds elapsed since the epoch.
 *
 * Resolution is implementation-dependent, but typically nanoseconds
 * or microseconds. */
double gmx_gettime();

#ifdef __cplusplus
}
#endif

#endif  /* GMX_TIMING_WALLTIME_ACCOUNTING_H */
