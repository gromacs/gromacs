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
#include "gromacs/legacyheaders/smalloc.h"
#include "gromacs/legacyheaders/sim_util.h"
#include "gromacs/legacyheaders/types/simple.h"

#include <time.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#elif defined HAVE_SYS_TIME_H
#include <sys/time.h>
#endif

/*! /brief Manages measuring wall clock times for simulations */
typedef struct gmx_runtime {
    double          start_time_stamp;            //!< Seconds since the epoch recorded at the start of the simulation
    double          start_time_stamp_per_thread; //!< Seconds since the epoch recorded at the start of the simulation for this thread
    double          elapsed_run_time;            //!< Total seconds elapsed over the simulation
    double          elapsed_run_time_per_thread; //!< Total seconds elapsed over the simulation running this thread
    gmx_large_int_t nsteps_done;                 //!< Used by integrators to report the amount of work they did
} t_gmx_runtime;

static double gmx_gettime_per_thread();

// TODO In principle, all this should get protected by checks that
// runtime is not null. In practice, that NULL condition does not
// happen, and future refactoring will likely enforce it by having
// the gmx_runtime_t object be owned by the runner object.

gmx_runtime_t
runtime_init()
{
    gmx_runtime_t runtime;

    snew(runtime, 1);
    runtime->start_time_stamp            = 0;
    runtime->start_time_stamp_per_thread = 0;
    runtime->elapsed_run_time            = 0;
    runtime->nsteps_done                 = 0;

    return runtime;
}

void
runtime_destroy(gmx_runtime_t runtime)
{
    sfree(runtime);
}

void
runtime_start(gmx_runtime_t runtime)
{
    runtime->start_time_stamp            = gmx_gettime();
    runtime->start_time_stamp_per_thread = gmx_gettime_per_thread();
    runtime->elapsed_run_time            = 0;
    runtime->nsteps_done                 = 0;
}

void
runtime_end(gmx_runtime_t runtime)
{
    double now, now_per_thread;

    now            = gmx_gettime();
    now_per_thread = gmx_gettime_per_thread();

    runtime->elapsed_run_time            = now - runtime->start_time_stamp;
    runtime->elapsed_run_time_per_thread = now_per_thread - runtime->start_time_stamp_per_thread;
}

double
runtime_get_current_elapsed_run_time(gmx_runtime_t runtime)
{
    return gmx_gettime() - runtime->start_time_stamp;
}

double
runtime_get_elapsed_run_time(gmx_runtime_t runtime)
{
    return runtime->elapsed_run_time;
}

double
runtime_get_elapsed_run_time_per_thread(gmx_runtime_t runtime)
{
    return runtime->elapsed_run_time_per_thread;
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
#ifdef _POSIX_TIMERS
    struct timespec t;
    double          seconds;

    clock_gettime(CLOCK_REALTIME, &t);
    seconds = (double) t.tv_sec + 1e-9*(double)t.tv_nsec;

    return seconds;
#elif defined HAVE_GETTIMEOFDAY
    // Note that gettimeofday() is deprecated by POSIX.
    struct timeval t;
    double         seconds;

    gettimeofday(&t, NULL);
    seconds = (double) t.tv_sec + 1e-6*(double)t.tv_usec;

    return seconds;
#else
    double  seconds;

    seconds = time(NULL);

    return seconds;
#endif
}

static double
gmx_gettime_per_thread()
{
#ifdef _POSIX_THREAD_CPUTIME
    struct timespec t;
    double          seconds;

    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t);
    seconds = (double) t.tv_sec + 1e-9*(double)t.tv_nsec;

    return seconds;
#else
    return gmx_gettime();
#endif
}
