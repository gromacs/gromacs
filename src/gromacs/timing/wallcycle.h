/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017,2018 by the GROMACS development team.
 * Copyright (c) 2019,2020,2021, by the GROMACS development team, led by
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
#ifndef GMX_TIMING_WALLCYCLE_H
#define GMX_TIMING_WALLCYCLE_H

/* NOTE: None of the routines here are safe to call within an OpenMP
 * region */

#include <stdio.h>

#include "gromacs/timing/cyclecounter.h"
#include "gromacs/utility/basedefinitions.h"

typedef struct gmx_wallcycle* gmx_wallcycle_t;
struct t_commrec;
static constexpr gmx_wallcycle* nullWallcycle = nullptr;

#ifndef DEBUG_WCYCLE
/*! \brief Enables consistency checking for the counters.
 *
 * If the macro is set to 1, code checks if you stop a counter different from the last
 * one that was opened and if you do nest too deep.
 */
#    define DEBUG_WCYCLE 0
#endif

enum
{
    ewcRUN,
    ewcSTEP,
    ewcPPDURINGPME,
    ewcDOMDEC,
    ewcDDCOMMLOAD,
    ewcDDCOMMBOUND,
    ewcVSITECONSTR,
    ewcPP_PMESENDX,
    ewcNS,
    ewcLAUNCH_GPU,
    ewcMOVEX,
    ewcFORCE,
    ewcMOVEF,
    ewcPMEMESH,
    ewcPME_REDISTXF,
    ewcPME_SPREAD,
    ewcPME_GATHER,
    ewcPME_FFT,
    ewcPME_FFTCOMM,
    ewcLJPME,
    ewcPME_SOLVE,
    ewcPMEWAITCOMM,
    ewcPP_PMEWAITRECVF,
    ewcWAIT_GPU_PME_SPREAD,
    ewcPME_FFT_MIXED_MODE,
    ewcPME_SOLVE_MIXED_MODE,
    ewcWAIT_GPU_PME_GATHER,
    ewcWAIT_GPU_BONDED,
    ewcPME_GPU_F_REDUCTION,
    ewcWAIT_GPU_NB_NL,
    ewcWAIT_GPU_NB_L,
    ewcWAIT_GPU_STATE_PROPAGATOR_DATA,
    ewcNB_XF_BUF_OPS,
    ewcVSITESPREAD,
    ewcPULLPOT,
    ewcAWH,
    ewcTRAJ,
    ewcUPDATE,
    ewcCONSTR,
    ewcMoveE,
    ewcROT,
    ewcROTadd,
    ewcSWAP,
    ewcIMD,
    ewcTEST,
    ewcNR
};

enum
{
    ewcsDD_REDIST,
    ewcsDD_GRID,
    ewcsDD_SETUPCOMM,
    ewcsDD_MAKETOP,
    ewcsDD_MAKECONSTR,
    ewcsDD_TOPOTHER,
    ewcsDD_GPU,
    ewcsNBS_GRID_LOCAL,
    ewcsNBS_GRID_NONLOCAL,
    ewcsNBS_SEARCH_LOCAL,
    ewcsNBS_SEARCH_NONLOCAL,
    ewcsLISTED,
    ewcsLISTED_FEP,
    ewcsRESTRAINTS,
    ewcsLISTED_BUF_OPS,
    ewcsNONBONDED_PRUNING,
    ewcsNONBONDED_KERNEL,
    ewcsNONBONDED_CLEAR,
    ewcsNONBONDED_FEP,
    ewcsLAUNCH_GPU_NONBONDED,
    ewcsLAUNCH_GPU_BONDED,
    ewcsLAUNCH_GPU_PME,
    ewcsLAUNCH_STATE_PROPAGATOR_DATA,
    ewcsEWALD_CORRECTION,
    ewcsNB_X_BUF_OPS,
    ewcsNB_F_BUF_OPS,
    ewcsCLEAR_FORCE_BUFFER,
    ewcsLAUNCH_GPU_NB_X_BUF_OPS,
    ewcsLAUNCH_GPU_NB_F_BUF_OPS,
    ewcsLAUNCH_GPU_MOVEX,
    ewcsLAUNCH_GPU_MOVEF,
    ewcsLAUNCH_GPU_UPDATE_CONSTRAIN,
    ewcsTEST,
    ewcsNR
};

static constexpr const bool sc_useCycleSubcounters = GMX_CYCLE_SUBCOUNTERS;

struct wallcc_t
{
    int          n;
    gmx_cycles_t c;
    gmx_cycles_t start;
};

#if DEBUG_WCYCLE
static constexpr int c_MaxWallCycleDepth = 6;
#endif


struct gmx_wallcycle
{
    wallcc_t* wcc;
    /* did we detect one or more invalid cycle counts */
    bool haveInvalidCount;
    /* variables for testing/debugging */
    bool      wc_barrier;
    wallcc_t* wcc_all;
    int       wc_depth;
#if DEBUG_WCYCLE
    int* counterlist;
    int  count_depth;
    bool isMasterRank;
#endif
    int              ewc_prev;
    gmx_cycles_t     cycle_prev;
    int64_t          reset_counters;
    const t_commrec* cr;
    wallcc_t*        wcsc;
};

//! Returns whether cycle counting is supported.
bool wallcycle_have_counter();

/*! \brief
 * Returns a wallcycle datastructure.
 *
 * If cycle counting is not supported, returns nullptr instead.
 */
gmx_wallcycle_t wallcycle_init(FILE* fplog, int resetstep, const t_commrec* cr);

//! Cleans up wallcycle structure.
void wallcycle_destroy(gmx_wallcycle_t wc);

//! Adds custom barrier for wallcycle counting.
void wallcycleBarrier(gmx_wallcycle* wc);

inline void wallcycle_all_start(gmx_wallcycle* wc, int ewc, gmx_cycles_t cycle)
{
    wc->ewc_prev   = ewc;
    wc->cycle_prev = cycle;
}

inline void wallcycle_all_stop(gmx_wallcycle* wc, int ewc, gmx_cycles_t cycle)
{
    const int prev    = wc->ewc_prev;
    const int current = ewc;
    wc->wcc_all[prev * ewcNR + current].n += 1;
    wc->wcc_all[prev * ewcNR + current].c += cycle - wc->cycle_prev;
}

//! Starts the cycle counter for \c ewc (and increases the call count).
inline void wallcycle_start(gmx_wallcycle_t wc, int ewc)
{
    if (wc == nullptr)
    {
        return;
    }

    wallcycleBarrier(wc);

#if DEBUG_WCYCLE
    debug_start_check(wc, ewc);
#endif
    gmx_cycles_t cycle = gmx_cycles_read();
    wc->wcc[ewc].start = cycle;
    if (wc->wcc_all)
    {
        wc->wc_depth++;
        if (ewc == ewcRUN)
        {
            wallcycle_all_start(wc, ewc, cycle);
        }
        else if (wc->wc_depth == 3)
        {
            wallcycle_all_stop(wc, ewc, cycle);
        }
    }
}

//! Starts the cycle counter for \c ewc without increasing the call count.
inline void wallcycle_start_nocount(gmx_wallcycle_t wc, int ewc)
{
    if (wc == nullptr)
    {
        return;
    }
    wc->wcc[ewc].n++;
}

//! Stop the cycle count for \c ewc, returns the last cycle count.
inline double wallcycle_stop(gmx_wallcycle_t wc, int ewc)
{
    gmx_cycles_t cycle, last;

    if (wc == nullptr)
    {
        return 0;
    }

    wallcycleBarrier(wc);

#if DEBUG_WCYCLE
    debug_stop_check(wc, ewc);
#endif

    /* When processes or threads migrate between cores, the cycle counting
     * can get messed up if the cycle counter on different cores are not
     * synchronized. When this happens we expect both large negative and
     * positive cycle differences. We can detect negative cycle differences.
     * Detecting too large positive counts if difficult, since count can be
     * large, especially for ewcRUN. If we detect a negative count,
     * we will not print the cycle accounting table.
     */
    cycle = gmx_cycles_read();
    if (cycle >= wc->wcc[ewc].start)
    {
        last = cycle - wc->wcc[ewc].start;
    }
    else
    {
        last                 = 0;
        wc->haveInvalidCount = true;
    }
    wc->wcc[ewc].c += last;
    wc->wcc[ewc].n++;
    if (wc->wcc_all)
    {
        wc->wc_depth--;
        if (ewc == ewcRUN)
        {
            wallcycle_all_stop(wc, ewc, cycle);
        }
        else if (wc->wc_depth == 2)
        {
            wallcycle_all_start(wc, ewc, cycle);
        }
    }

    return last;
}

//! Only increment call count for \c ewc by one.
inline void wallcycle_increment_event_count(gmx_wallcycle_t wc, int ewc)
{
    if (wc == nullptr)
    {
        return;
    }
    wc->wcc[ewc].n++;
}

//! Returns the cumulative count and cycle count for \c ewc.
void wallcycle_get(gmx_wallcycle_t wc, int ewc, int* n, double* c);

//! Returns the cumulative count and sub cycle count for \c ewcs.
void wallcycle_sub_get(gmx_wallcycle_t wc, int ewcs, int* n, double* c);

//! Resets all cycle counters to zero.
void wallcycle_reset_all(gmx_wallcycle_t wc);

//! Scale the cycle counts to reflect how many threads run for that number of cycles.
void wallcycle_scale_by_num_threads(gmx_wallcycle_t wc, bool isPmeRank, int nthreads_pp, int nthreads_pme);

//! Return reset_counters.
int64_t wcycle_get_reset_counters(gmx_wallcycle_t wc);

//! Set reset_counters.
void wcycle_set_reset_counters(gmx_wallcycle_t wc, int64_t reset_counters);

//! Set the start sub cycle count for \c ewcs.
inline void wallcycle_sub_start(gmx_wallcycle_t wc, int ewcs)
{
    if (sc_useCycleSubcounters && wc != nullptr)
    {
        wc->wcsc[ewcs].start = gmx_cycles_read();
    }
}

//! Set the start sub cycle count for \c ewcs without increasing the call count.
inline void wallcycle_sub_start_nocount(gmx_wallcycle_t wc, int ewcs)
{
    if (sc_useCycleSubcounters && wc != nullptr)
    {
        wc->wcsc[ewcs].start = gmx_cycles_read();
    }
}

//! Stop the sub cycle count for \c ewcs.
inline void wallcycle_sub_stop(gmx_wallcycle_t wc, int ewcs)
{
    if (sc_useCycleSubcounters && wc != nullptr)
    {
        wc->wcsc[ewcs].c += gmx_cycles_read() - wc->wcsc[ewcs].start;
        wc->wcsc[ewcs].n++;
    }
}

#endif
