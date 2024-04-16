/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#ifndef GMX_TIMING_WALLCYCLE_H
#define GMX_TIMING_WALLCYCLE_H

/* NOTE: None of the routines here are safe to call within an OpenMP
 * region */

#include <cstdio>

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/timing/cyclecounter.h"
#include "gromacs/timing/instrumentation.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/enumerationhelpers.h"

#ifndef DEBUG_WCYCLE
/*! \brief Enables consistency checking for the counters.
 *
 * If the macro is set to 1, code checks if you stop a counter different from the last
 * one that was opened and if you do nest too deep.
 */
#    define DEBUG_WCYCLE 0
#endif

struct t_commrec;

enum class WallCycleCounter : int
{
    Run,
    Step,
    PpDuringPme,
    Domdec,
    DDCommLoad,
    DDCommBound,
    VsiteConstr,
    PpPmeSendX,
    NS,
    LaunchGpuPp,
    MoveX,
    Force,
    MoveF,
    PmeMesh,
    PmeGpuMesh, /* PmeGpuMesh is used for GPU code and similar to PmeMesh on CPU. It includes WaitGpuPmePPRecvX cycles too. */
    PmeRedistXF,
    PmeSpread,
    PmeGather,
    PmeFft,
    PmeFftComm,
    LJPme,
    PmeSolve,
    WaitGpuPmeGridD2hCopy, /* Time for PME grid D2H transfer. Used in mixed mode. */
    PmeFftMixedMode,
    PmeSolveMixedMode,
    WaitGpuPmeGather,
    PmeGpuFReduction,
    LaunchGpuPme, /* Time for launching PME specific GPU operations*/
    WaitGpuPmePPRecvX, /* Time for waiting on receiving PP X on PME rank. Used only when GPU direct comm is active.*/
    WaitGpuPmeSpread, /* Time taken to finish PME spread on GPU. Used only when PME halo-exchange is active with PME decomposition*/
    WaitGpuFftToPmeGrid, /* Time taken to convert to PME grid after FFTs are complete. Used only when PME halo-exchange is active with PME decomposition*/
    PmeHaloExchangeComm, /* Time taken in PME halo-exchange, active with PME decomposition*/
    PmeWaitComm, /* PmeWaitComm = Run - PmeMesh. Without GPU direct comm, this includes time spent in waiting for coord and force comm.
                 With GPU direct comm, waiting for coord comm is part of PME mesh and is measured with WaitGpuPmePPRecvX sub-counter*/
    PpPmeWaitRecvF,
    WaitGpuBonded,
    WaitGpuNbNL,
    WaitGpuNbL,
    WaitGpuStatePropagatorData,
    NbXFBufOps,
    VsiteSpread,
    PullPot,
    Awh,
    Traj,
    Update,
    Constr,
    MoveE,
    Rot,
    RotAdd,
    Swap,
    Imd,
    MdGpuGraph,
    Test,
    Count
};

enum class WallCycleSubCounter : int
{
    DDRedist,
    DDGrid,
    DDSetupComm,
    DDMakeTop,
    DDMakeConstr,
    DDTopOther,
    DDGpu,
    NBSGridLocal,
    NBSGridNonLocal,
    NBSSearchLocal,
    NBSSearchNonLocal,
    GpuBondedListUpdate,
    Listed,
    ListedFep,
    Restraints,
    ListedBufOps,
    NonbondedPruning,
    NonbondedKernel,
    NonbondedClear,
    NonbondedFep,
    NonbondedFepReduction,
    LaunchGpuNonBonded,
    LaunchGpuBonded,
    LaunchStatePropagatorData,
    EwaldCorrection,
    NBXBufOps,
    NBFBufOps,
    ClearForceBuffer,
    LaunchGpuNBXBufOps,
    LaunchGpuNBFBufOps,
    LaunchGpuMoveX,
    LaunchGpuMoveF,
    LaunchGpuUpdateConstrain,
    LaunchGpuPmeFft, /* Time for launching FFT operations on GPU*/
    MdGpuGraphWaitBeforeCapture,
    MdGpuGraphCapture,
    MdGpuGraphInstantiateOrUpdate,
    MdGpuGraphWaitBeforeLaunch,
    MdGpuGraphLaunch,
    ConstrComm,
    Test,
    Count
};

template<int maxLength, typename Container>
static constexpr bool checkStringsLengths(const Container& strings)
{
    // NOLINTNEXTLINE(readability-use-anyofallof) // std::all_of is constexpr only since C++20
    for (const char* str : strings)
    {
        if (std::char_traits<char>::length(str) > maxLength)
        {
            return false;
        }
    }
    return true;
}

/* Each name should not exceed 22 printing characters
   (ie. terminating null can be twentieth) */
static const char* enumValuetoString(WallCycleCounter enumValue)
{
    constexpr gmx::EnumerationArray<WallCycleCounter, const char*> wallCycleCounterNames = {
        "Run",
        "Step",
        "PP during PME",
        "Domain decomp.",
        "DD comm. load",
        "DD comm. bounds",
        "Vsite constr.",
        "Send X to PME",
        "Neighbor search",
        "Launch PP GPU ops.",
        "Comm. coord.",
        "Force",
        "Wait + Comm. F",
        "PME mesh",
        "PME GPU mesh",
        "PME redist. X/F",
        "PME spread",
        "PME gather",
        "PME 3D-FFT",
        "PME 3D-FFT Comm.",
        "PME solve LJ",
        "PME solve Elec",
        "Wait PME GPU D2H",
        "PME 3D-FFT",
        "PME solve",
        "Wait PME GPU gather",
        "Reduce GPU PME F",
        "Launch PME GPU ops.",
        "Wait PME Recv. PP X",
        "Wait PME GPU spread",
        "Wait GPU FFT to PME",
        "PME Halo exch comm",
        "PME wait for PP",
        "Wait + Recv. PME F",
        "Wait Bonded GPU",
        "Wait GPU NB nonloc.",
        "Wait GPU NB local",
        "Wait GPU state copy",
        "NB X/F buffer ops.",
        "Vsite spread",
        "COM pull force",
        "AWH",
        "Write traj.",
        "Update",
        "Constraints",
        "Comm. energies",
        "Enforced rotation",
        "Add rot. forces",
        "Position swapping",
        "IMD",
        "MD Graph",
        "Test"
    };
    static_assert(checkStringsLengths<22>(wallCycleCounterNames));
    return wallCycleCounterNames[enumValue];
}

CLANG_DIAGNOSTIC_IGNORE("-Wunneeded-internal-declaration")
static const char* enumValuetoString(WallCycleSubCounter enumValue)
{
    constexpr gmx::EnumerationArray<WallCycleSubCounter, const char*> wallCycleSubCounterNames = {
        "DD redist.",
        "DD NS grid + sort",
        "DD setup comm.",
        "DD make top.",
        "DD make constr.",
        "DD top. other",
        "DD GPU ops.",
        "NS grid local",
        "NS grid non-local",
        "NS search local",
        "NS search non-local",
        "GPU Bonded list update",
        "Bonded F",
        "Bonded-FEP F",
        "Restraints F",
        "Listed buffer ops.",
        "NB pruning",
        "NB F kernel",
        "NB F clear",
        "NB FEP",
        "NB FEP reduction",
        "Launch GPU NB tasks",
        "Launch GPU Bonded",
        "Launch state copy",
        "Ewald F correction",
        "NB X buffer ops.",
        "NB F buffer ops.",
        "Clear force buffer",
        "Launch GPU NB X ops.",
        "Launch GPU NB F ops.",
        "Launch GPU Comm. X",
        "Launch GPU Comm. F",
        "Launch GPU update",
        "Launch PME GPU FFT",
        "Graph wait pre-capture",
        "Graph capture",
        "Graph instantiate/upd.",
        "Graph wait pre-launch",
        "Graph launch",
        "Constraints Comm.", // constraints communication time, note that this counter will contain load imbalance
        "Test subcounter"
    };
    static_assert(checkStringsLengths<22>(wallCycleSubCounterNames));
    return wallCycleSubCounterNames[enumValue];
}
CLANG_DIAGNOSTIC_RESET


//! Number of all main counters.
static constexpr int sc_numWallCycleCounters = static_cast<int>(WallCycleCounter::Count);
//! Number of all subcyclecounters.
static constexpr int sc_numWallCycleSubCounters = static_cast<int>(WallCycleSubCounter::Count);
//! Scare of all counters for keeping track.
static constexpr int sc_numWallCycleCountersSquared = sc_numWallCycleCounters * sc_numWallCycleCounters;
//! Do we keep track of sub cycle counters.
static constexpr bool sc_useCycleSubcounters = GMX_CYCLE_SUBCOUNTERS;
//! Whether wallcycle debugging is enabled.
constexpr bool sc_enableWallcycleDebug = (DEBUG_WCYCLE != 0);

//! Counters for an individual wallcycle timing region
struct wallcc_t
{
    //! Counter for number of times this timing region has been opened
    int n = 0;
    //! Counter for total number of cycles in this timing region
    gmx_cycles_t c = 0;
    //! Start time (in cycles) for the last time this timing region was opened
    gmx_cycles_t start;
};

struct gmx_wallcycle
{
    /*! \brief Methods used when debugging wallcycle counting
     *
     *  \todo Make these private when the functions they are called
     *  from become class methods. */
    //! \{
    //! Check the start preconditions
    void checkStart(WallCycleCounter ewc);
    //! Check the stop preconditions
    void checkStop(WallCycleCounter ewc);
    //! \}

public:
    //! Storage for wallcycle counters
    gmx::EnumerationArray<WallCycleCounter, wallcc_t> wcc;
    //! The step count at which counter reset will happen
    int64_t reset_counters;
    //! Storage for wallcycle subcounters
    gmx::EnumerationArray<WallCycleSubCounter, wallcc_t> wcsc;

    // The remaining fields are only used in special cases or in
    // printing summary output.

    //! Commrec for communicator used when using wallcycle barriers
    const t_commrec* cr;

    //! Used when doing "all" wallcycle counting
    //! \{
    //! All counters
    std::vector<wallcc_t> wcc_all;
    //! Counter depth
    int wc_depth = 0;
    //! Previous counter index
    WallCycleCounter ewc_prev = WallCycleCounter::Count;
    //! Previous cycle count value
    gmx_cycles_t cycle_prev;
    //! \}

    //! Did we detect one or more invalid cycle counts?
    bool haveInvalidCount = false;
    //! Add extra barriers during testing
    bool wc_barrier = false;

    //! Used when debugging wallcycle counting
    //! \{
    //! Maximum depth of counters
    static constexpr int sc_maxDepth = sc_enableWallcycleDebug ? 6 : 0;
    //! Counters
    std::array<WallCycleCounter, sc_maxDepth> counterlist;
    //! Counter depth
    int count_depth = 0;
    //! Whether this rank is the main rank of the simulation
    bool isMainRank = false;
    //! \}
};

//! Returns if cycle counting is supported
bool wallcycle_have_counter();

//! Returns the wall cycle structure.
std::unique_ptr<gmx_wallcycle> wallcycle_init(FILE* fplog, int resetstep, const t_commrec* cr);

//! Adds custom barrier for wallcycle counting.
void wallcycleBarrier(gmx_wallcycle* wc);

void wallcycle_sub_get(gmx_wallcycle* wc, WallCycleSubCounter ewcs, int* n, double* c);
/* Returns the cumulative count and sub cycle count for ewcs */

inline void wallcycle_all_start(gmx_wallcycle* wc, WallCycleCounter ewc, gmx_cycles_t cycle)
{
    wc->ewc_prev   = ewc;
    wc->cycle_prev = cycle;
}

inline void wallcycle_all_stop(gmx_wallcycle* wc, WallCycleCounter ewc, gmx_cycles_t cycle)
{
    const int prev    = static_cast<int>(wc->ewc_prev);
    const int current = static_cast<int>(ewc);
    wc->wcc_all[prev * sc_numWallCycleCounters + current].n += 1;
    wc->wcc_all[prev * sc_numWallCycleCounters + current].c += cycle - wc->cycle_prev;
}

//! Starts the cycle counter (and increases the call count)
inline void wallcycle_start(gmx_wallcycle* wc, WallCycleCounter ewc)
{
    if (ewc >= WallCycleCounter::Step)
    {
        traceRangeStart(enumValuetoString(ewc), static_cast<int>(ewc));
    }

    if (wc == nullptr)
    {
        return;
    }

    wallcycleBarrier(wc);

    if constexpr (sc_enableWallcycleDebug)
    {
        wc->checkStart(ewc);
    }
    gmx_cycles_t cycle = gmx_cycles_read();
    wc->wcc[ewc].start = cycle;
    if (!wc->wcc_all.empty())
    {
        wc->wc_depth++;
        if (ewc == WallCycleCounter::Run)
        {
            wallcycle_all_start(wc, ewc, cycle);
        }
        else if (wc->wc_depth == 3)
        {
            wallcycle_all_stop(wc, ewc, cycle);
        }
    }
}

//! Starts the cycle counter without increasing the call count
inline void wallcycle_start_nocount(gmx_wallcycle* wc, WallCycleCounter ewc)
{
    if (wc == nullptr)
    {
        return;
    }
    wallcycle_start(wc, ewc);
    wc->wcc[ewc].n--;
}

//! Stop the cycle count for ewc , returns the last cycle count
inline double wallcycle_stop(gmx_wallcycle* wc, WallCycleCounter ewc)
{
    if (ewc >= WallCycleCounter::Step)
    {
        traceRangeEnd();
    }

    gmx_cycles_t cycle, last;

    if (wc == nullptr)
    {
        return 0;
    }

    wallcycleBarrier(wc);

    if constexpr (sc_enableWallcycleDebug)
    {
        wc->checkStop(ewc);
    }

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
    if (!wc->wcc_all.empty())
    {
        wc->wc_depth--;
        if (ewc == WallCycleCounter::Run)
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

//! Only increment call count for ewc by one
inline void wallcycle_increment_event_count(gmx_wallcycle* wc, WallCycleCounter ewc)
{
    if (wc == nullptr)
    {
        return;
    }
    wc->wcc[ewc].n++;
}

//! Returns the cumulative count and cycle count for ewc
void wallcycle_get(gmx_wallcycle* wc, WallCycleCounter ewc, int* n, double* c);

//! Resets all cycle counters to zero
void wallcycle_reset_all(gmx_wallcycle* wc);

//! Scale the cycle counts to reflect how many threads run for that number of cycles
void wallcycle_scale_by_num_threads(gmx_wallcycle* wc, bool isPmeRank, int nthreads_pp, int nthreads_pme);

//! Return reset_counters from wc struct
int64_t wcycle_get_reset_counters(gmx_wallcycle* wc);

//! Set reset_counters
void wcycle_set_reset_counters(gmx_wallcycle* wc, int64_t reset_counters);

//! Set the start sub cycle count for ewcs
inline void wallcycle_sub_start(gmx_wallcycle* wc, WallCycleSubCounter ewcs)
{
    if constexpr (sc_useCycleSubcounters)
    {
        traceSubRangeStart(enumValuetoString(ewcs), static_cast<int>(ewcs));

        if (wc != nullptr)
        {
            wc->wcsc[ewcs].start = gmx_cycles_read();
        }
    }
}

//! Set the start sub cycle count for ewcs without increasing the call count
inline void wallcycle_sub_start_nocount(gmx_wallcycle* wc, WallCycleSubCounter ewcs)
{
    if constexpr (sc_useCycleSubcounters)
    {
        if (wc != nullptr)
        {
            wallcycle_sub_start(wc, ewcs);
            wc->wcsc[ewcs].n--;
        }
    }
}

//! Stop the sub cycle count for ewcs
inline void wallcycle_sub_stop(gmx_wallcycle* wc, WallCycleSubCounter ewcs)
{
    if constexpr (sc_useCycleSubcounters)
    {
        traceSubRangeEnd();

        if (wc != nullptr)
        {
            wc->wcsc[ewcs].c += gmx_cycles_read() - wc->wcsc[ewcs].start;
            wc->wcsc[ewcs].n++;
        }
    }
}

#endif
