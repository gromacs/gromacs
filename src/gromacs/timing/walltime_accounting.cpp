/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "walltime_accounting.h"

#include "config.h"

#include <ctime>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/smalloc.h"

/* TODO in future: convert gmx_walltime_accounting to a class,
 * resolve who should have responsibility for recording the number of
 * steps done, consider whether parts of finish_time, print_perf,
 * wallcycle_print belong in this module.
 *
 * If/when any kind of task parallelism is implemented (even OpenMP
 * regions simultaneously assigned to different tasks), consider
 * whether this data structure (and/or cycle counters) should be
 * maintained on a per-OpenMP-thread basis.
 *
 * Consider also replacing this with std::chrono. */

/*! \brief Manages caching wall-clock time measurements for
 * simulations */
typedef struct gmx_walltime_accounting {
    //! Seconds since the epoch recorded at the start of the simulation
    double          start_time_stamp;
    //! Seconds since the epoch recorded at the start of the simulation for this thread
    double          start_time_stamp_per_thread;
    //! Total seconds elapsed over the simulation
    double          elapsed_time;
    //! Total seconds elapsed over the simulation running this thread
    double          elapsed_time_over_all_threads;
    /*! \brief Number of OpenMP threads that will be launched by this
     * MPI rank.
     *
     * This is used to scale elapsed_time_over_all_threads so
     * that any combination of real MPI, thread MPI and OpenMP (even
     * mdrun -ntomp_pme) processes/threads would (when run at maximum
     * efficiency) return values such that the sum of
     * elapsed_time_over_all_threads over all threads was constant
     * with respect to parallelism implementation. */
    int             numOpenMPThreads;
    //! Set by integrators to report the amount of work they did
    gmx_int64_t     nsteps_done;
    //! Whether the simulation has finished in a way valid for walltime reporting.
    bool            isValidFinish;
} t_gmx_walltime_accounting;

/*! \brief Calls system timing routines (e.g. clock_gettime) to get
 * the (fractional) number of seconds elapsed since the epoch when
 * this thread was executing.
 *
 * This can be used to measure system load. This can be unreliable if
 * threads migrate between sockets. If thread-specific timers are not
 * supported by the OS (e.g. if the OS is not POSIX-compliant), this
 * function is implemented by gmx_gettime. */
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
    walltime_accounting->isValidFinish               = false;

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

gmx_int64_t
walltime_accounting_get_nsteps_done(gmx_walltime_accounting_t walltime_accounting)
{
    return walltime_accounting->nsteps_done;
}

void
walltime_accounting_set_nsteps_done(gmx_walltime_accounting_t   walltime_accounting,
                                    gmx_int64_t                 nsteps_done)
{
    walltime_accounting->nsteps_done = nsteps_done;
}

double
gmx_gettime()
{
    /* Use clock_gettime only if we know linking the C run-time
       library will work (which is not trivial on e.g. Crays), and its
       headers claim sufficient support for POSIX (ie not Mac and
       Windows), and it isn't BG/Q (whose compute node kernel only
       supports gettimeofday, and bgclang doesn't provide a fully
       functional implementation clock_gettime). */
#if HAVE_CLOCK_GETTIME && defined(_POSIX_TIMERS) && _POSIX_TIMERS > 0 && !(defined __bgq__ && defined __clang__)
    struct timespec t;
    double          seconds;

    clock_gettime(CLOCK_REALTIME, &t);
    seconds = static_cast<double>(t.tv_sec) + 1e-9*t.tv_nsec;

    return seconds;
#elif HAVE_GETTIMEOFDAY
    // Note that gettimeofday() is deprecated by POSIX, but since Mac
    // and Windows do not yet support POSIX, we are still stuck.
    // Also, this is the only supported API call on Bluegene/Q.
    struct timeval t;
    double         seconds;

    gettimeofday(&t, NULL);
    seconds = static_cast<double>(t.tv_sec) + 1e-6*t.tv_usec;

    return seconds;
#else
    double  seconds;

    seconds = time(NULL);

    return seconds;
#endif
}

void
walltime_accounting_set_valid_finish(gmx_walltime_accounting_t walltime_accounting)
{
    walltime_accounting->isValidFinish = true;
}

//! Return whether the simulation finished in a way valid for reporting walltime.
bool
walltime_accounting_get_valid_finish(const gmx_walltime_accounting_t walltime_accounting)
{
    return walltime_accounting->isValidFinish;
}

static double
gmx_gettime_per_thread()
{
    /* Use clock_gettime only if we know linking the C run-time
       library will work (which is not trivial on e.g. Crays), and its
       headers claim sufficient support for POSIX (ie not Mac and
       Windows), and it isn't BG/Q (whose compute node kernel only
       supports gettimeofday, and bgclang doesn't provide a fully
       functional implementation clock_gettime). */
#if HAVE_CLOCK_GETTIME && defined(_POSIX_THREAD_CPUTIME) && _POSIX_THREAD_CPUTIME > 0 && !(defined __bgq__ && defined __clang__)
    struct timespec t;
    double          seconds;

    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t);
    seconds = static_cast<double>(t.tv_sec) + 1e-9*t.tv_nsec;

    return seconds;
#else
    return gmx_gettime();
#endif
}
