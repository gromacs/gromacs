/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
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
 * GROwing Monsters And Cloning Shrimps
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "gromacs/legacyheaders/runtime.h"
#include "gromacs/legacyheaders/sim_util.h"
#include "gromacs/legacyheaders/types/simple.h"

/*! /brief Manages measuring wall clock times for simulations */
typedef struct gmx_runtime {
    double          start_time_stamp; //!< Seconds since the epoch recorded at the start of the simulation
    double          elapsed_run_time; //!< Total seconds elapsed over the simulation
    gmx_large_int_t nsteps_done;      //!< Used by integrators to report the amount of work they did
} t_gmx_runtime;

// TODO In principle, all this should get protected by checks that
// runtime is not null. In practice, that seems impossible, and future
// refactoring will likely enforce it.

void
runtime_start(gmx_runtime_t runtime)
{
    runtime->start_time_stamp = gmx_gettime();
    runtime->elapsed_run_time = 0;
    runtime->nsteps_done      = 0;
}

void
runtime_end(gmx_runtime_t runtime)
{
    double now;

    now = gmx_gettime();

    runtime->elapsed_run_time = now - runtime->start_time_stamp;
    runtime->start_time_stamp = now;
}

double
runtime_get_elapsed_run_time(gmx_runtime_t runtime)
{
    return gmx_gettime() - runtime->start_time_stamp;
}

double
runtime_get_start_time_stamp(gmx_runtime_t runtime)
{
    return runtime->start_time_stamp;
}

double
runtime_get_nsteps_done(gmx_runtime_t runtime)
{
    return runtime->nsteps_done;
}

void
runtime_set_nsteps_done(gmx_runtime_t   runtime,
                        gmx_large_int_t _nsteps_done)
{
    runtime->nsteps_done = _nsteps_done;
}

double
gmx_gettime()
{
#ifdef HAVE_GETTIMEOFDAY
    struct timeval t;
    double         seconds;

    // TODO later: gettimeofday() is deprecated by POSIX. We could use
    // clock_gettime in POSIX (which also offers nanosecond resolution
    // if the hardware supports it), but that requires linking with
    // -lrt. Maybe a better option will come along before we have to
    // really change from gettimeofday().
    gettimeofday(&t, NULL);

    seconds = (double) t.tv_sec + 1e-6*(double)t.tv_usec;

    return seconds;
#else
    double  seconds;

    seconds = time(NULL);

    return seconds;
#endif
}
