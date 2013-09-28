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
 * Copyright (c) 2013, The GROMACS development team,
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

#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/legacyheaders/smalloc.h"
#include "gromacs/legacyheaders/types/simple.h"

#include <time.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif

/* TODO in future: convert gmx_walltime_accounting to a class,
 * resolve who should have responsibility for recording the number of
 * steps done, consider whether parts of finish_time, print_perf,
 * wallcycle_print belong in this module.
 *
 * If/when any kind of task parallelism is implemented (even OpenMP
 * regions simultaneously assigned to different tasks), consider
 * whether this data structure (and/or cycle counters) should be
 * maintained on a per-OpenMP-thread basis. */

/*! Manages caching wall-clock time measurements for simulations */
typedef struct gmx_walltime_accounting {
    //! Seconds since the epoch recorded at the start of the simulation
    double          start_time_stamp;
    //! Seconds since the epoch recorded at the start of the simulation for this thread
    double          start_time_stamp_per_thread;
    //! Total seconds elapsed over the simulation
    double          elapsed_time;
    //! Total seconds elapsed over the simulation running this thread
    double          elapsed_time_over_all_threads;
    //! Number of OpenMP threads that will be launched by this MPI
    // rank. This is used to scale elapsed_time_over_all_threads so that any
    // combination of real MPI, thread MPI and OpenMP (even mdrun
    // -ntomp_pme) processes/threads would (when run at maximum
    // efficiency) return values such that the sum of
    // elapsed_time_over_all_threads over all threads was constant with
    // respect to parallelism implementation.
    int             numOpenMPThreads;
    //! Set by integrators to report the amount of work they did
    gmx_large_int_t nsteps_done;
} t_gmx_walltime_accounting;

//! Calls system timing routines (e.g. clock_gettime) to get the
// (fractional) number of seconds elapsed since the epoch when this
// thread was executing (if supported by the POSIX-compliant OS). This
// can be used to measure system load. This can be unreliable if
// threads migrate between sockets. If thread-specific timers are not
// supported by the OS, is implemented by gmx_gettime.
static double gmx_gettime_per_thread();

// TODO In principle, all this should get protected by checks that
// walltime_accounting is not null. In practice, that NULL condition
// does not happen, and future refactoring will likely enforce it by
// having the gmx_walltime_accounting_t object be owned by the runner
// object. When these become member functions, existence will be
// guaranteed.

gmx_walltime_accounting_t
walltime_accounting_init(int numOpenMPThreads)
{
    gmx_walltime_accounting_t walltime_accounting;

    snew(walltime_accounting, 1);
    walltime_accounting->start_time_stamp            = 0;
    walltime_accounting->start_time_stamp_per_thread = 0;
    walltime_accounting->elapsed_time                = 0;
    walltime_accounting->nsteps_done                 = 0;
    walltime_accounting->numOpenMPThreads            = numOpenMPThreads;

    return walltime_accounting;
}

void
walltime_accounting_destroy(gmx_walltime_accounting_t walltime_accounting)
{
    sfree(walltime_accounting);
}

void
walltime_accounting_start(gmx_walltime_accounting_t walltime_accounting)
{
    walltime_accounting->start_time_stamp            = gmx_gettime();
    walltime_accounting->start_time_stamp_per_thread = gmx_gettime_per_thread();
    walltime_accounting->elapsed_time                = 0;
    walltime_accounting->nsteps_done                 = 0;
}

void
walltime_accounting_end(gmx_walltime_accounting_t walltime_accounting)
{
    double now, now_per_thread;

    now            = gmx_gettime();
    now_per_thread = gmx_gettime_per_thread();

    walltime_accounting->elapsed_time                  = now - walltime_accounting->start_time_stamp;
    walltime_accounting->elapsed_time_over_all_threads = now_per_thread - walltime_accounting->start_time_stamp_per_thread;
    /* For thread-MPI, the per-thread CPU timer makes this just
     * work. For OpenMP threads, the per-thread CPU timer measurement
     * needs to be multiplied by the number of OpenMP threads used,
     * under the current assumption that all regions ever opened
     * within a process are of the same size, and each thread should
     * keep one core busy.
     */
    walltime_accounting->elapsed_time_over_all_threads *= walltime_accounting->numOpenMPThreads;
}

double
walltime_accounting_get_current_elapsed_time(gmx_walltime_accounting_t walltime_accounting)
{
    return gmx_gettime() - walltime_accounting->start_time_stamp;
}

double
walltime_accounting_get_elapsed_time(gmx_walltime_accounting_t walltime_accounting)
{
    return walltime_accounting->elapsed_time;
}

double
walltime_accounting_get_elapsed_time_over_all_threads(gmx_walltime_accounting_t walltime_accounting)
{
    return walltime_accounting->elapsed_time_over_all_threads;
}

double
walltime_accounting_get_start_time_stamp(gmx_walltime_accounting_t walltime_accounting)
{
    return walltime_accounting->start_time_stamp;
}

double
walltime_accounting_get_nsteps_done(gmx_walltime_accounting_t walltime_accounting)
{
    return walltime_accounting->nsteps_done;
}

void
walltime_accounting_set_nsteps_done(gmx_walltime_accounting_t   walltime_accounting,
                                    gmx_large_int_t             nsteps_done)
{
    walltime_accounting->nsteps_done = nsteps_done;
}

double
gmx_gettime()
{
#if defined HAVE_CLOCK_GETTIME && _POSIX_TIMERS >= 0
    /* Mac and Windows do not support this. For added fun, Windows
     * defines _POSIX_TIMERS without actually providing the
     * implementation. */
    struct timespec t;
    double          seconds;

    clock_gettime(CLOCK_REALTIME, &t);
    seconds = (double) t.tv_sec + 1e-9*(double)t.tv_nsec;

    return seconds;
#elif defined HAVE_GETTIMEOFDAY
    // Note that gettimeofday() is deprecated by POSIX, but since Mac
    // and Windows do not yet support POSIX, we are still stuck.
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
#if defined HAVE_CLOCK_GETTIME && _POSIX_THREAD_CPUTIME >= 0
    struct timespec t;
    double          seconds;

    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t);
    seconds = (double) t.tv_sec + 1e-9*(double)t.tv_nsec;

    return seconds;
#else
    return gmx_gettime();
#endif
}
