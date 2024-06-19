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
#include "gmxpre.h"

#include "gromacs/timing/wallcycle.h"

#include "config.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include <array>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/math/functions.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/timing/cyclecounter.h"
#include "gromacs/timing/gpu_timing.h"
#include "gromacs/timing/wallcyclereporting.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/snprintf.h"
#include "gromacs/utility/stringutil.h"

//! True if only the main rank should print debugging output
static constexpr bool sc_onlyMainDebugPrints = true;
//! True if cycle counter nesting depth debugging prints are enabled
static constexpr bool sc_debugPrintDepth = false;


#if GMX_USE_ITT
CLANG_DIAGNOSTIC_IGNORE("-Wzero-as-null-pointer-constant")
//! Heler for Intel tracing tools instrumentation
const __itt_domain*  g_ittDomain = __itt_domain_create("GMX");
__itt_string_handle* g_ittCounterHandles[static_cast<int>(WallCycleCounter::Count)];
__itt_string_handle* g_ittSubCounterHandles[static_cast<int>(WallCycleSubCounter::Count)];
#endif


/* PME GPU timing events' names - correspond to the enum in the gpu_timing.h */
static const char* enumValuetoString(PmeStage enumValue)
{
    constexpr gmx::EnumerationArray<PmeStage, const char*> pmeStageNames = {
        "PME spline", "PME spread",     "PME spline + spread", "PME 3D-FFT r2c",
        "PME solve",  "PME 3D-FFT c2r", "PME gather"
    };
    return pmeStageNames[enumValue];
}

bool wallcycle_have_counter()
{
    return gmx_cycles_have_counter();
}

std::unique_ptr<gmx_wallcycle> wallcycle_init(FILE* fplog, int resetstep, const t_commrec* cr)
{
    std::unique_ptr<gmx_wallcycle> wc;


    if (!wallcycle_have_counter())
    {
        return wc;
    }

    wc = std::make_unique<gmx_wallcycle>();

    wc->reset_counters = resetstep;
    wc->cr             = cr;


#if GMX_MPI
    if (cr != nullptr && PAR(cr) && getenv("GMX_CYCLE_BARRIER") != nullptr)
    {
        if (fplog)
        {
            fprintf(fplog, "\nWill call MPI_Barrier before each cycle start/stop call\n\n");
        }
        wc->wc_barrier = true;
    }
#endif

    if (getenv("GMX_CYCLE_ALL") != nullptr)
    {
        if (fplog)
        {
            fprintf(fplog, "\nWill time all the code during the run\n\n");
        }
        wc->wcc_all.resize(sc_numWallCycleCountersSquared);
    }

    // NOLINTNEXTLINE(readability-misleading-indentation)
    if constexpr (sc_enableWallcycleDebug)
    {
        wc->isMainRank = (cr == nullptr) || MAIN(cr);
    }

#if GMX_USE_ITT
    for (auto wcc : gmx::EnumerationWrapper<WallCycleCounter>{})
    {
        g_ittCounterHandles[static_cast<int>(wcc)] = __itt_string_handle_create(enumValuetoString(wcc));
    }
    for (auto wcsc : gmx::EnumerationWrapper<WallCycleSubCounter>{})
    {
        g_ittSubCounterHandles[static_cast<int>(wcsc)] =
                __itt_string_handle_create(enumValuetoString(wcsc));
    }
#endif

    return wc;
}

#if GMX_USE_ITT
CLANG_DIAGNOSTIC_RESET
#endif

void gmx_wallcycle::checkStart(WallCycleCounter ewc)
{
    // NOLINTNEXTLINE(readability-misleading-indentation)
    if constexpr (sc_enableWallcycleDebug)
    {
        // NOLINTNEXTLINE(misc-redundant-expression)
        if (count_depth < 0 || count_depth >= sc_maxDepth)
        {
            gmx_fatal(FARGS, "wallcycle counter depth out of range: %d", count_depth + 1);
        }
        counterlist[count_depth] = ewc;
        count_depth++;

        if (sc_debugPrintDepth && (!sc_onlyMainDebugPrints || isMainRank))
        {
            std::string indentStr(4 * count_depth, ' ');
            fprintf(stderr, "%swcycle_start depth %d, %s\n", indentStr.c_str(), count_depth, enumValuetoString(ewc));
        }
    }
}

void gmx_wallcycle::checkStop(WallCycleCounter ewc)
{
    // NOLINTNEXTLINE(readability-misleading-indentation)
    if constexpr (sc_enableWallcycleDebug)
    {
        if (sc_debugPrintDepth && (!sc_onlyMainDebugPrints || isMainRank))
        {
            std::string indentStr(4 * count_depth, ' ');
            fprintf(stderr, "%swcycle_stop  depth %d, %s\n", indentStr.c_str(), count_depth, enumValuetoString(ewc));
        }

        count_depth--;

        if (count_depth < 0)
        {
            gmx_fatal(FARGS,
                      "wallcycle counter depth out of range when stopping %s: %d",
                      enumValuetoString(ewc),
                      count_depth);
        }
        if (counterlist[count_depth] != ewc)
        {
            gmx_fatal(FARGS,
                      "wallcycle mismatch at stop, start %s, stop %s",
                      enumValuetoString(counterlist[count_depth]),
                      enumValuetoString(ewc));
        }
    }
}

void wallcycle_get(gmx_wallcycle* wc, WallCycleCounter ewc, int* n, double* c)
{
    *n = wc->wcc[ewc].n;
    *c = static_cast<double>(wc->wcc[ewc].c);
}

void wallcycle_sub_get(gmx_wallcycle* wc, WallCycleSubCounter ewcs, int* n, double* c)
{
    // NOLINTNEXTLINE(readability-misleading-indentation)
    if constexpr (sc_useCycleSubcounters)
    {
        if (wc != nullptr)
        {
            *n = wc->wcsc[ewcs].n;
            *c = static_cast<double>(wc->wcsc[ewcs].c);
        }
    }
}

void wallcycle_reset_all(gmx_wallcycle* wc)
{
    if (wc == nullptr)
    {
        return;
    }

    for (auto& counter : wc->wcc)
    {
        counter.n = 0;
        counter.c = 0;
    }
    wc->haveInvalidCount = false;

    if (!wc->wcc_all.empty())
    {
        for (int i = 0; i < sc_numWallCycleCountersSquared; i++)
        {
            wc->wcc_all[i].n = 0;
            wc->wcc_all[i].c = 0;
        }
    }

    // NOLINTNEXTLINE(readability-misleading-indentation)
    if constexpr (sc_useCycleSubcounters)
    {
        for (auto& counter : wc->wcsc)
        {
            counter.n = 0;
            counter.c = 0;
        }
    }
}

static bool is_pme_counter(WallCycleCounter ewc)
{
    return (ewc >= WallCycleCounter::PmeMesh && ewc <= WallCycleCounter::PmeWaitComm);
}

static bool is_pme_subcounter(WallCycleCounter ewc)
{
    return (ewc >= WallCycleCounter::PmeRedistXF && ewc < WallCycleCounter::PmeWaitComm);
}

void wallcycleBarrier(gmx_wallcycle* wc)
{
#if GMX_MPI
    if (wc->wc_barrier)
    {
        MPI_Barrier(wc->cr->mpi_comm_mygroup);
    }
#else
    GMX_UNUSED_VALUE(wc);
#endif
}

/* Subtract counter ewc_sub timed inside a timing block for ewc_main */
// NOLINTNEXTLINE(google-runtime-references)
static void subtract_cycles(gmx::EnumerationArray<WallCycleCounter, wallcc_t>& wcc,
                            WallCycleCounter                                   ewc_main,
                            WallCycleCounter                                   ewc_sub)
{
    if (wcc[ewc_sub].n > 0)
    {
        if (wcc[ewc_main].c >= wcc[ewc_sub].c)
        {
            wcc[ewc_main].c -= wcc[ewc_sub].c;
        }
        else
        {
            /* Something is wrong with the cycle counting */
            wcc[ewc_main].c = 0;
        }
    }
}

void wallcycle_scale_by_num_threads(gmx_wallcycle* wc, bool isPmeRank, int nthreads_pp, int nthreads_pme)
{
    if (wc == nullptr)
    {
        return;
    }

    for (auto key : keysOf(wc->wcc))
    {
        if (is_pme_counter(key) || (key == WallCycleCounter::Run && isPmeRank))
        {
            wc->wcc[key].c *= nthreads_pme;

            if (!wc->wcc_all.empty())
            {
                const int current = static_cast<int>(key);
                for (int j = 0; j < sc_numWallCycleCounters; j++)
                {
                    wc->wcc_all[current * sc_numWallCycleCounters + j].c *= nthreads_pme;
                }
            }
        }
        else
        {
            wc->wcc[key].c *= nthreads_pp;

            if (!wc->wcc_all.empty())
            {
                const int current = static_cast<int>(key);
                for (int j = 0; j < sc_numWallCycleCounters; j++)
                {
                    wc->wcc_all[current * sc_numWallCycleCounters + j].c *= nthreads_pp;
                }
            }
        }
    }
    // NOLINTNEXTLINE(readability-misleading-indentation)
    if constexpr (sc_useCycleSubcounters)
    {
        if (!isPmeRank)
        {
            for (auto& counter : wc->wcsc)
            {
                counter.c *= nthreads_pp;
            }
        }
    }
}

/* TODO Make an object for this function to return, containing some
 * vectors of something like wallcc_t for the summed wcc, wcc_all and
 * wcsc, AND the original wcc for rank 0.
 *
 * The GPU timing is reported only for rank 0, so we want to preserve
 * the original wcycle on that rank. Rank 0 also reports the global
 * counts before that, so needs something to contain the global data
 * without over-writing the rank-0 data. The current implementation
 * uses cycles_sum to manage this, which works OK now because wcsc and
 * wcc_all are unused by the GPU reporting, but it is not satisfactory
 * for the future. Also, there's no need for MPI_Allreduce, since
 * only MAINRANK uses any of the results. */
WallcycleCounts wallcycle_sum(const t_commrec* cr, gmx_wallcycle* wc)
{
    WallcycleCounts                                    cycles_sum;
    gmx::EnumerationArray<WallCycleCounter, double>    cyclesMain;
    gmx::EnumerationArray<WallCycleSubCounter, double> cyclesSub;
#if GMX_MPI
    gmx::EnumerationArray<WallCycleCounter, double>    cyclesMainOnNode;
    gmx::EnumerationArray<WallCycleSubCounter, double> cyclesSubOnNode;
#endif

    if (wc == nullptr)
    {
        /* Default construction of std::array of non-class T can leave
           the values indeterminate, just like a C array */
        cycles_sum.fill(0);
        return cycles_sum;
    }

    auto& wcc = wc->wcc;

    subtract_cycles(wcc, WallCycleCounter::Domdec, WallCycleCounter::DDCommLoad);
    subtract_cycles(wcc, WallCycleCounter::Domdec, WallCycleCounter::DDCommBound);

    subtract_cycles(wcc, WallCycleCounter::PmeFft, WallCycleCounter::PmeFftComm);

    if (cr->npmenodes == 0)
    {
        /* All nodes do PME (or no PME at all) */
        subtract_cycles(wcc, WallCycleCounter::Force, WallCycleCounter::PmeMesh);
    }
    else
    {
        /* The are PME-only nodes */
        if (wcc[WallCycleCounter::PmeMesh].n > 0)
        {
            GMX_ASSERT(wcc[WallCycleCounter::PmeGpuMesh].c == 0,
                       "PME mesh GPU ticks should be 0 when PME mesh is running on CPU");
            /* This must be a PME only node, calculate the Wait + Comm. time */
            GMX_ASSERT(wcc[WallCycleCounter::Run].c >= wcc[WallCycleCounter::PmeMesh].c,
                       "Total run ticks must be greater than PME-only ticks");
            wcc[WallCycleCounter::PmeWaitComm].c =
                    wcc[WallCycleCounter::Run].c - wcc[WallCycleCounter::PmeMesh].c;
        }

        if (wcc[WallCycleCounter::PmeGpuMesh].n > 0)
        {
            GMX_ASSERT(wcc[WallCycleCounter::PmeMesh].c == 0,
                       "PME mesh CPU ticks should be 0 when PME mesh is running on GPU");

            /* This must be a PME only node, calculate the Wait + Comm. time */
            GMX_ASSERT(wcc[WallCycleCounter::Run].c >= wcc[WallCycleCounter::PmeGpuMesh].c,
                       "Total run ticks must be greater than PME-only ticks");
            wcc[WallCycleCounter::PmeWaitComm].c =
                    wcc[WallCycleCounter::Run].c - wcc[WallCycleCounter::PmeGpuMesh].c;
        }
    }

    /* Store the cycles in a double buffer for summing */
    for (auto key : keysOf(wcc))
    {
#if GMX_MPI
        cyclesMainOnNode[key] = static_cast<double>(wcc[key].n);
#endif
        cyclesMain[key] = static_cast<double>(wcc[key].c);
    }
    // NOLINTNEXTLINE(readability-misleading-indentation)
    if constexpr (sc_useCycleSubcounters)
    {
        for (auto key : keysOf(wc->wcsc))
        {
#if GMX_MPI
            cyclesSubOnNode[key] = static_cast<double>(wc->wcsc[key].n);
#endif
            cyclesSub[key] = static_cast<double>(wc->wcsc[key].c);
        }
    }

#if GMX_MPI
    if (cr->nnodes > 1)
    {
        gmx::EnumerationArray<WallCycleCounter, double>    bufMain;
        gmx::EnumerationArray<WallCycleSubCounter, double> bufSub;

        // TODO this code is used only at the end of the run, so we
        // can just do a simple reduce of haveInvalidCount in
        // wallcycle_print, and avoid bugs
        double haveInvalidCount = (wc->haveInvalidCount ? 1 : 0);
        // TODO Use MPI_Reduce
        MPI_Allreduce(cyclesMainOnNode.data(), bufMain.data(), bufMain.size(), MPI_DOUBLE, MPI_MAX, cr->mpi_comm_mysim);
        // NOLINTNEXTLINE(readability-misleading-indentation)
        if constexpr (sc_useCycleSubcounters)
        {
            MPI_Allreduce(cyclesSubOnNode.data(), bufSub.data(), bufSub.size(), MPI_DOUBLE, MPI_MAX, cr->mpi_comm_mysim);
        }
        MPI_Allreduce(MPI_IN_PLACE, &haveInvalidCount, 1, MPI_DOUBLE, MPI_MAX, cr->mpi_comm_mysim);
        for (auto key : keysOf(wcc))
        {
            wcc[key].n = gmx::roundToInt(bufMain[key]);
        }
        wc->haveInvalidCount = (haveInvalidCount > 0);
        // NOLINTNEXTLINE(readability-misleading-indentation)
        if constexpr (sc_useCycleSubcounters)
        {
            for (auto key : keysOf(wc->wcsc))
            {
                wc->wcsc[key].n = gmx::roundToInt(bufSub[key]);
            }
        }

        // TODO Use MPI_Reduce
        MPI_Allreduce(cyclesMain.data(), cycles_sum.data(), cyclesMain.size(), MPI_DOUBLE, MPI_SUM, cr->mpi_comm_mysim);
        // NOLINTNEXTLINE(readability-misleading-indentation)
        if constexpr (sc_useCycleSubcounters)
        {
            MPI_Allreduce(cyclesSub.data(),
                          cycles_sum.data() + sc_numWallCycleCounters,
                          cyclesSub.size(),
                          MPI_DOUBLE,
                          MPI_SUM,
                          cr->mpi_comm_mysim);
        }

        if (!wc->wcc_all.empty())
        {
            std::array<double, sc_numWallCycleCountersSquared> cyc_all;
            std::array<double, sc_numWallCycleCountersSquared> buf_all;

            for (int i = 0; i < sc_numWallCycleCountersSquared; i++)
            {
                cyc_all[i] = wc->wcc_all[i].c;
            }
            // TODO Use MPI_Reduce
            MPI_Allreduce(cyc_all.data(),
                          buf_all.data(),
                          sc_numWallCycleCountersSquared,
                          MPI_DOUBLE,
                          MPI_SUM,
                          cr->mpi_comm_mysim);
            for (int i = 0; i < sc_numWallCycleCountersSquared; i++)
            {
                wc->wcc_all[i].c = static_cast<gmx_cycles_t>(buf_all[i]);
            }
        }
    }
    else
#endif
    {
        for (auto key : keysOf(cyclesMain))
        {
            cycles_sum[static_cast<int>(key)] = cyclesMain[key];
        }
        // NOLINTNEXTLINE(readability-misleading-indentation)
        if constexpr (sc_useCycleSubcounters)
        {
            for (auto key : keysOf(cyclesSub))
            {
                const int offset   = static_cast<int>(key) + sc_numWallCycleCounters;
                cycles_sum[offset] = cyclesSub[key];
            }
        }
    }

    return cycles_sum;
}

static void
print_cycles(FILE* fplog, double c2t, const char* name, int nnodes, int nthreads, int ncalls, double c_sum, double tot)
{
    char   nnodes_str[STRLEN];
    char   nthreads_str[STRLEN];
    char   ncalls_str[STRLEN];
    double wallt;
    double percentage = (tot > 0.) ? (100. * c_sum / tot) : 0.;

    if (c_sum > 0)
    {
        if (ncalls > 0)
        {
            snprintf(ncalls_str, sizeof(ncalls_str), "%10d", ncalls);
            if (nnodes < 0)
            {
                snprintf(nnodes_str, sizeof(nnodes_str), "N/A");
            }
            else
            {
                snprintf(nnodes_str, sizeof(nnodes_str), "%4d", nnodes);
            }
            if (nthreads < 0)
            {
                snprintf(nthreads_str, sizeof(nthreads_str), "N/A");
            }
            else
            {
                snprintf(nthreads_str, sizeof(nthreads_str), "%4d", nthreads);
            }
        }
        else
        {
            nnodes_str[0]   = 0;
            nthreads_str[0] = 0;
            ncalls_str[0]   = 0;
        }
        /* Convert the cycle count to wallclock time for this task */
        wallt = c_sum * c2t;

        fprintf(fplog,
                " %-22.22s %4s %4s %10s  %10.3f %14.3f %5.1f\n",
                name,
                nnodes_str,
                nthreads_str,
                ncalls_str,
                wallt,
                c_sum * 1e-9,
                percentage);
    }
}

static void print_gputimes(FILE* fplog, const char* name, int n, double t, double tot_t)
{
    char num[11];
    char avg_perf[11];

    if (n > 0)
    {
        snprintf(num, sizeof(num), "%10d", n);
        snprintf(avg_perf, sizeof(avg_perf), "%10.3f", t / n);
    }
    else
    {
        sprintf(num, "          ");
        sprintf(avg_perf, "          ");
    }
    if (t != tot_t && tot_t > 0)
    {
        fprintf(fplog, " %-29s %10s%12.3f   %s   %5.1f\n", name, num, t / 1000, avg_perf, 100 * t / tot_t);
    }
    else
    {
        fprintf(fplog, " %-29s %10s%12.3f   %s   %5.1f\n", name, "", t / 1000, avg_perf, 100.0);
    }
}

static void print_header(FILE* fplog, int nrank_pp, int nth_pp, int nrank_pme, int nth_pme)
{
    int nrank_tot = nrank_pp + nrank_pme;
    if (0 == nrank_pme)
    {
        fprintf(fplog, "On %d MPI rank%s", nrank_tot, nrank_tot == 1 ? "" : "s");
        if (nth_pp > 1)
        {
            fprintf(fplog, ", each using %d OpenMP threads", nth_pp);
        }
        /* Don't report doing PP+PME, because we can't tell here if
         * this is RF, etc. */
    }
    else
    {
        fprintf(fplog, "On %d MPI rank%s doing PP", nrank_pp, nrank_pp == 1 ? "" : "s");
        if (nth_pp > 1)
        {
            fprintf(fplog, ",%s using %d OpenMP threads", nrank_pp > 1 ? " each" : "", nth_pp);
        }
        fprintf(fplog, ", and\non %d MPI rank%s doing PME", nrank_pme, nrank_pme == 1 ? "" : "s");
        if (nth_pme > 1)
        {
            fprintf(fplog, ",%s using %d OpenMP threads", nrank_pme > 1 ? " each" : "", nth_pme);
        }
    }

    fprintf(fplog, "\n\n");
    fprintf(fplog, " Activity:              Num   Num      Call    Wall time         Giga-Cycles\n");
    fprintf(fplog,
            "                        Ranks Threads  Count      (s)         total sum    %%\n");
}


void wallcycle_print(FILE*                            fplog,
                     const gmx::MDLogger&             mdlog,
                     int                              nnodes,
                     int                              npme,
                     int                              nth_pp,
                     int                              nth_pme,
                     double                           realtime,
                     gmx_wallcycle*                   wc,
                     const WallcycleCounts&           cyc_sum,
                     const gmx_wallclock_gpu_nbnxn_t* gpu_nbnxn_t,
                     const gmx_wallclock_gpu_pme_t*   gpu_pme_t)
{
    double      tot, tot_for_pp, tot_for_rest, tot_cpu_overlap, gpu_cpu_ratio;
    double      c2t, c2t_pp, c2t_pme = 0;
    int         npp, nth_tot;
    char        buf[STRLEN];
    const char* hline =
            "--------------------------------------------------------------------------------";

    if (wc == nullptr)
    {
        return;
    }

    GMX_ASSERT(nth_pp > 0, "Number of particle-particle threads must be >0");
    GMX_ASSERT(nth_pme > 0, "Number of PME threads must be >0");
    GMX_ASSERT(nnodes > 0, "Number of nodes must be >0");
    GMX_ASSERT(npme >= 0, "Number of PME nodes cannot be negative");
    npp = nnodes - npme;
    /* npme is the number of PME-only ranks used, and we always do PP work */
    GMX_ASSERT(npp > 0, "Number of particle-particle nodes must be >0");

    nth_tot = npp * nth_pp + npme * nth_pme;

    /* When using PME-only nodes, the next line is valid for both
       PP-only and PME-only nodes because they started ewcRUN at the
       same time. */
    tot        = cyc_sum[static_cast<int>(WallCycleCounter::Run)];
    tot_for_pp = 0;

    if (tot <= 0.0)
    {
        /* TODO This is heavy handed, but until someone reworks the
           code so that it is provably robust with respect to
           non-positive values for all possible timer and cycle
           counters, there is less value gained from printing whatever
           timing data might still be sensible for some non-CI
           run, than is lost from diagnosing CI FP exceptions on
           runs about whose execution time we don't care. */
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendTextFormatted(
                        "WARNING: A total of %f CPU cycles was recorded, so mdrun cannot print a "
                        "time accounting",
                        tot);
        return;
    }

    if (wc->haveInvalidCount)
    {
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendText(
                        "NOTE: Detected invalid cycle counts, probably because threads moved "
                        "between CPU cores that do not have synchronized cycle counters. Will not "
                        "print the cycle accounting.");
        return;
    }


    /* Conversion factor from cycles to seconds */
    c2t    = realtime / tot;
    c2t_pp = c2t * nth_tot / static_cast<double>(npp * nth_pp);
    if (npme > 0)
    {
        c2t_pme = c2t * nth_tot / static_cast<double>(npme * nth_pme);
    }
    else
    {
        c2t_pme = 0;
    }

    fprintf(fplog, "\n      R E A L   C Y C L E   A N D   T I M E   A C C O U N T I N G\n\n");

    print_header(fplog, npp, nth_pp, npme, nth_pme);

    fprintf(fplog, "%s\n", hline);
    gmx::EnumerationWrapper<WallCycleCounter> iter;
    for (auto key = gmx::EnumerationIterator<WallCycleCounter>(WallCycleCounter::Domdec);
         key != iter.end();
         ++key)
    {

        if (is_pme_subcounter(*key))
        {
            /* Do not count these at all */
        }
        else if (npme > 0 && is_pme_counter(*key))
        {
            /* Print timing information for PME-only nodes, but add an
             * asterisk so the reader of the table can know that the
             * walltimes are not meant to add up. The asterisk still
             * fits in the required maximum of 19 characters. */
            std::string message = gmx::formatString("%s *", enumValuetoString(*key));
            print_cycles(fplog,
                         c2t_pme,
                         message.c_str(),
                         npme,
                         nth_pme,
                         wc->wcc[*key].n,
                         cyc_sum[static_cast<int>(*key)],
                         tot);
        }
        else
        {
            /* Print timing information when it is for a PP or PP+PME
               node */
            print_cycles(fplog,
                         c2t_pp,
                         enumValuetoString(*key),
                         npp,
                         nth_pp,
                         wc->wcc[*key].n,
                         cyc_sum[static_cast<int>(*key)],
                         tot);
            tot_for_pp += cyc_sum[static_cast<int>(*key)];
        }
    }
    if (!wc->wcc_all.empty())
    {
        for (auto i : keysOf(wc->wcc))
        {
            const int countI = static_cast<int>(i);
            for (auto j : keysOf(wc->wcc))
            {
                const int countJ = static_cast<int>(j);
                snprintf(buf, 20, "%-9.9s %-9.9s", enumValuetoString(i), enumValuetoString(j));
                print_cycles(fplog,
                             c2t_pp,
                             buf,
                             npp,
                             nth_pp,
                             wc->wcc_all[countI * sc_numWallCycleCounters + countJ].n,
                             wc->wcc_all[countI * sc_numWallCycleCounters + countJ].c,
                             tot);
            }
        }
    }
    tot_for_rest = tot * npp * nth_pp / static_cast<double>(nth_tot);
    print_cycles(fplog, c2t_pp, "Rest", npp, nth_pp, -1, tot_for_rest - tot_for_pp, tot);
    fprintf(fplog, "%s\n", hline);
    print_cycles(fplog, c2t, "Total", npp, nth_pp, -1, tot, tot);
    fprintf(fplog, "%s\n", hline);

    if (npme > 0)
    {
        fprintf(fplog,
                "(*) Note that with separate PME ranks, the walltime column actually sums to\n"
                "    twice the total reported, but the cycle count total and %% are correct.\n"
                "%s\n",
                hline);
    }

    if (wc->wcc[WallCycleCounter::PmeMesh].n > 0 || wc->wcc[WallCycleCounter::PmeGpuMesh].n > 0)
    {
        // A workaround to not print breakdown when no subcounters were recorded.
        // TODO: figure out and record PME GPU counters (what to do with the waiting ones?)
        std::vector<WallCycleCounter> validPmeSubcounterIndices;
        for (auto key = gmx::EnumerationIterator<WallCycleCounter>(WallCycleCounter::Domdec);
             key != iter.end();
             key++)
        {
            if (is_pme_subcounter(*key) && wc->wcc[*key].n > 0)
            {
                validPmeSubcounterIndices.push_back(*key);
            }
        }

        if (!validPmeSubcounterIndices.empty())
        {
            fprintf(fplog, " Breakdown of PME mesh activities\n");
            fprintf(fplog, "%s\n", hline);
            for (auto i : validPmeSubcounterIndices)
            {
                print_cycles(fplog,
                             npme > 0 ? c2t_pme : c2t_pp,
                             enumValuetoString(i),
                             npme > 0 ? npme : npp,
                             nth_pme,
                             wc->wcc[i].n,
                             cyc_sum[static_cast<int>(i)],
                             tot);
            }
            fprintf(fplog, "%s\n", hline);
        }
    }

    // NOLINTNEXTLINE(readability-misleading-indentation)
    if constexpr (sc_useCycleSubcounters)
    {
        fprintf(fplog, " Breakdown of PP / PME activities\n");
        fprintf(fplog, "%s\n", hline);
        for (auto key : keysOf(wc->wcsc))
        {
            print_cycles(fplog,
                         c2t_pp,
                         enumValuetoString(key),
                         npp,
                         nth_pp,
                         wc->wcsc[key].n,
                         cyc_sum[sc_numWallCycleCounters + static_cast<int>(key)],
                         tot);
        }
        fprintf(fplog, "%s\n", hline);
    }

    /* print GPU timing summary */
    double tot_gpu = 0.0;
    if (gpu_pme_t)
    {
        for (auto key : keysOf(gpu_pme_t->timing))
        {
            tot_gpu += gpu_pme_t->timing[key].t;
        }
    }
    if (gpu_nbnxn_t)
    {
        const char* k_log_str[2][2] = { { "Nonbonded F kernel", "Nonbonded F+ene k." },
                                        { "Nonbonded F+prune k.", "Nonbonded F+ene+prune k." } };
        tot_gpu += gpu_nbnxn_t->pl_h2d_t + gpu_nbnxn_t->nb_h2d_t + gpu_nbnxn_t->nb_d2h_t;

        /* add up the kernel timings */
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                tot_gpu += gpu_nbnxn_t->ktime[i][j].t;
            }
        }
        tot_gpu += gpu_nbnxn_t->pruneTime.t;

        tot_cpu_overlap = wc->wcc[WallCycleCounter::Force].c;
        if (wc->wcc[WallCycleCounter::PmeMesh].n > 0)
        {
            tot_cpu_overlap += wc->wcc[WallCycleCounter::PmeMesh].c;
        }
        tot_cpu_overlap *= realtime * 1000 / tot; /* convert s to ms */

        fprintf(fplog, "\n GPU timings\n%s\n", hline);
        fprintf(fplog,
                " Computing:                         Count  Wall t (s)      ms/step       %c\n",
                '%');
        fprintf(fplog, "%s\n", hline);
        print_gputimes(fplog, "Pair list H2D", gpu_nbnxn_t->pl_h2d_c, gpu_nbnxn_t->pl_h2d_t, tot_gpu);
        print_gputimes(fplog, "X / q H2D", gpu_nbnxn_t->nb_c, gpu_nbnxn_t->nb_h2d_t, tot_gpu);

        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                if (gpu_nbnxn_t->ktime[i][j].c)
                {
                    print_gputimes(fplog,
                                   k_log_str[i][j],
                                   gpu_nbnxn_t->ktime[i][j].c,
                                   gpu_nbnxn_t->ktime[i][j].t,
                                   tot_gpu);
                }
            }
        }
        if (gpu_pme_t)
        {
            for (auto key : keysOf(gpu_pme_t->timing))
            {
                if (gpu_pme_t->timing[key].c)
                {
                    print_gputimes(fplog,
                                   enumValuetoString(key),
                                   gpu_pme_t->timing[key].c,
                                   gpu_pme_t->timing[key].t,
                                   tot_gpu);
                }
            }
        }
        if (gpu_nbnxn_t->pruneTime.c)
        {
            print_gputimes(fplog, "Pruning kernel", gpu_nbnxn_t->pruneTime.c, gpu_nbnxn_t->pruneTime.t, tot_gpu);
        }
        print_gputimes(fplog, "F D2H", gpu_nbnxn_t->nb_c, gpu_nbnxn_t->nb_d2h_t, tot_gpu);
        fprintf(fplog, "%s\n", hline);
        print_gputimes(fplog, "Total ", gpu_nbnxn_t->nb_c, tot_gpu, tot_gpu);
        fprintf(fplog, "%s\n", hline);
        if (gpu_nbnxn_t->dynamicPruneTime.c)
        {
            /* We print the dynamic pruning kernel timings after a separator
             * and avoid adding it to tot_gpu as this is not in the force
             * overlap. We print the fraction as relative to the rest.
             */
            print_gputimes(fplog,
                           "*Dynamic pruning",
                           gpu_nbnxn_t->dynamicPruneTime.c,
                           gpu_nbnxn_t->dynamicPruneTime.t,
                           tot_gpu);
            fprintf(fplog, "%s\n", hline);
        }
        gpu_cpu_ratio = tot_gpu / tot_cpu_overlap;
        if (gpu_nbnxn_t->nb_c > 0 && wc->wcc[WallCycleCounter::Force].n > 0)
        {
            fprintf(fplog,
                    "\nAverage per-step force GPU/CPU evaluation time ratio: %.3f ms/%.3f ms = "
                    "%.3f\n",
                    tot_gpu / gpu_nbnxn_t->nb_c,
                    tot_cpu_overlap / wc->wcc[WallCycleCounter::Force].n,
                    gpu_cpu_ratio);
        }

        /* only print notes related to CPU-GPU load balance with PME */
        if (wc->wcc[WallCycleCounter::PmeMesh].n > 0)
        {
            fprintf(fplog, "For optimal resource utilization this ratio should be close to 1\n");

            /* print note if the imbalance is high with PME case in which
             * CPU-GPU load balancing is possible */
            if (gpu_cpu_ratio < 0.8 || gpu_cpu_ratio > 1.25)
            {
                /* Only the sim main calls this function, so always print to stderr */
                if (gpu_cpu_ratio < 0.8)
                {
                    if (npp > 1)
                    {
                        /* The user could have used -notunepme,
                         * but we currently can't check that here.
                         */
                        GMX_LOG(mdlog.warning)
                                .asParagraph()
                                .appendText(
                                        "NOTE: The CPU has >25% more load than the GPU. This "
                                        "imbalance wastes\n"
                                        "      GPU resources. Maybe the domain decomposition "
                                        "limits the PME tuning.\n"
                                        "      In that case, try setting the DD grid manually "
                                        "(-dd) or lowering -dds.");
                    }
                    else
                    {
                        /* We should not end up here, unless the box is
                         * too small for increasing the cut-off for PME tuning.
                         */
                        GMX_LOG(mdlog.warning)
                                .asParagraph()
                                .appendText(
                                        "NOTE: The CPU has >25% more load than the GPU. This "
                                        "imbalance wastes\n"
                                        "      GPU resources.");
                    }
                }
                if (gpu_cpu_ratio > 1.25)
                {
                    GMX_LOG(mdlog.warning)
                            .asParagraph()
                            .appendText(
                                    "NOTE: The GPU has >25% more load than the CPU. This imbalance "
                                    "wastes\n"
                                    "      CPU resources.");
                }
            }
        }
    }

    if (wc->wc_barrier)
    {
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendText(
                        "MPI_Barrier was called before each cycle start/stop\n"
                        "call, so timings are not those of real runs.");
    }

    if (wc->wcc[WallCycleCounter::NbXFBufOps].n > 0
        && (cyc_sum[static_cast<int>(WallCycleCounter::Domdec)] > tot * 0.1
            || cyc_sum[static_cast<int>(WallCycleCounter::NS)] > tot * 0.1))
    {
        /* Only the sim main calls this function, so always print to stderr */
        if (wc->wcc[WallCycleCounter::Domdec].n == 0)
        {
            GMX_LOG(mdlog.warning)
                    .asParagraph()
                    .appendTextFormatted(
                            "NOTE: %d %% of the run time was spent in pair search,\n"
                            "      you might want to increase nstlist (this has no effect on "
                            "accuracy)\n",
                            gmx::roundToInt(100 * cyc_sum[static_cast<int>(WallCycleCounter::NS)] / tot));
        }
        else
        {
            GMX_LOG(mdlog.warning)
                    .asParagraph()
                    .appendTextFormatted(
                            "NOTE: %d %% of the run time was spent in domain decomposition,\n"
                            "      %d %% of the run time was spent in pair search,\n"
                            "      you might want to increase nstlist (this has no effect on "
                            "accuracy)\n",
                            gmx::roundToInt(100 * cyc_sum[static_cast<int>(WallCycleCounter::Domdec)] / tot),
                            gmx::roundToInt(100 * cyc_sum[static_cast<int>(WallCycleCounter::NS)] / tot));
        }
    }

    if (cyc_sum[static_cast<int>(WallCycleCounter::MoveE)] > tot * 0.05)
    {
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendTextFormatted(
                        "NOTE: %d %% of the run time was spent communicating energies,\n"
                        "      you might want to increase some nst* mdp options\n",
                        gmx::roundToInt(100 * cyc_sum[static_cast<int>(WallCycleCounter::MoveE)] / tot));
    }
}

int64_t wcycle_get_reset_counters(gmx_wallcycle* wc)
{
    if (wc == nullptr)
    {
        return -1;
    }
    return wc->reset_counters;
}

void wcycle_set_reset_counters(gmx_wallcycle* wc, int64_t reset_counters)
{
    if (wc == nullptr)
    {
        return;
    }
    wc->reset_counters = reset_counters;
}
