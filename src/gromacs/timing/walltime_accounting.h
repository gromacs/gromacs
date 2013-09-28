/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 *
 * And Hey:
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef GMX_TIMING_WALLTIME_ACCOUNTING_H
#define GMX_TIMING_WALLTIME_ACCOUNTING_H

#include "gromacs/legacyheaders/types/simple.h"

#ifdef __cplusplus
extern "C" {
#endif

#if 0
}
#endif

/*! Contains per-process and per-thread data about elapsed wall-clock
 *  times and integration steps performed.
 */
typedef struct gmx_walltime_accounting *gmx_walltime_accounting_t;

//! Constructor
gmx_walltime_accounting_t
walltime_accounting_init(int numOpenMPThreads);

//! Destructor
void
walltime_accounting_destroy(gmx_walltime_accounting_t walltime_accounting);

//! Record initial time stamps, e.g. at run end or counter
//! re-initalization time
void
walltime_accounting_start(gmx_walltime_accounting_t walltime_accounting);

//! Measure and cache the elapsed wall-clock time since
//! walltime_accounting_start
void
walltime_accounting_end(gmx_walltime_accounting_t walltime_accounting);

//! Measure and return the elapsed wall-clock time since
//! walltime_accounting_start
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
double
walltime_accounting_get_nsteps_done(gmx_walltime_accounting_t walltime_accounting);

//! Set the number of integration steps done
//
// TODO consider whether this should get done in walltime_accounting_end
void
walltime_accounting_set_nsteps_done(gmx_walltime_accounting_t   walltime_accounting,
                                    gmx_large_int_t             nsteps_done);

//! Calls system timing routines (e.g. clock_gettime) to get the
// (fractional) number of seconds elapsed since the epoch.
//
// Resolution is implementation-dependent, but typically nanoseconds
// or microseconds.
double gmx_gettime();

#ifdef __cplusplus
}
#endif

#endif  /* GMX_TIMING_WALLTIME_ACCOUNTING_H */
