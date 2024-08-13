/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
/*! \internal \file
 *
 * \brief This file contains function definitions necessary for
 * managing automatic load balance of PME calculations (Coulomb and
 * LJ).
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_ewald
 */
#include "gmxpre.h"

#include "pme_load_balancing.h"

#include <cassert>
#include <cmath>

#include <algorithm>

#include "gromacs/domdec/dlb.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_network.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/partition.h"
#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/fft/calcgrid.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/dispersioncorrection.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/nbnxm/gpu_data_mgmt.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/strconvert.h"

#include "pme_internal.h"
#include "pme_pp.h"

/*! \brief Parameters and settings for one PP-PME setup */
struct pme_setup_t
{
    real rcut_coulomb;         /**< Coulomb cut-off                              */
    real rlistOuter;           /**< cut-off for the outer pair-list              */
    real rlistInner;           /**< cut-off for the inner pair-list              */
    real spacing;              /**< (largest) PME grid spacing                   */
    ivec grid;                 /**< the PME grid dimensions                      */
    real grid_efficiency;      /**< ineffiency factor for non-uniform grids <= 1 */
    real ewaldcoeff_q;         /**< Electrostatic Ewald coefficient            */
    real ewaldcoeff_lj;        /**< LJ Ewald coefficient, only for the call to send_switchgrid */
    struct gmx_pme_t* pmedata; /**< the data structure used in the PME code      */
    int               count;   /**< number of times this setup has been timed    */
    double            cycles;  /**< the fastest time for this setup in cycles    */
};

/*! \brief After 50 nstlist periods of not observing imbalance: never tune PME */
const int PMETunePeriod = 50;
/*! \brief Trigger PME load balancing at more than 5% PME overload */
const real loadBalanceTriggerFactor = 1.05;
/*! \brief Scale the grid by a most at factor 1.7.
 *
 * This still leaves room for about 4-4.5x decrease in grid spacing while limiting the cases where
 * large imbalance leads to extreme cutoff scaling for marginal benefits.
 *
 * This should help to avoid:
 *   - large increase in power consumption for little performance gain
 *   - increasing communication volume
 *   - limiting DLB
 */
const real c_maxSpacingScaling = 1.7;
/*! \brief In the initial scan, step by grids that are at least a factor 0.8 coarser */
const real gridpointsScaleFactor = 0.8;
/*! \brief In the initial scan, try to skip grids with uneven x/y/z spacing,
 * checking if the "efficiency" is more than 5% worse than the previous grid.
 */
const real relativeEfficiencyFactor = 1.05;
/*! \brief Rerun until a run is 12% slower setups than the fastest run so far */
const real maxRelativeSlowdownAccepted = 1.12;
/*! \brief If setups get more than 2% faster, do another round to avoid
 * choosing a slower setup due to acceleration or fluctuations.
 */
const real maxFluctuationAccepted = 1.02;

//! \brief Number of nstlist long tuning intervals to skip before starting
//         load-balancing at the beginning of the run.
const int c_numFirstTuningIntervalSkip = 5;
//! \brief Number of nstlist long tuning intervals to skip before starting
//         load-balancing at the beginning of the run with separate PME ranks. */
const int c_numFirstTuningIntervalSkipWithSepPme = 3;
//! \brief Number of nstlist long tuning intervals to skip after switching to a new setting
//         during balancing.
const int c_numPostSwitchTuningIntervalSkip = 1;
//! \brief Number of seconds to delay the tuning at startup to allow processors clocks to ramp up.
const double c_startupTimeDelay = 5.0;

/*! \brief Enumeration whose values describe the effect limiting the load balancing */
enum class PmeLoadBalancingLimit : int
{
    No,
    Box,
    DD,
    PmeGrid,
    MaxScaling,
    Count
};

/*! \brief Descriptive strings for PmeLoadBalancingLimit \c enumValue */
static const char* enumValueToString(PmeLoadBalancingLimit enumValue)
{
    constexpr gmx::EnumerationArray<PmeLoadBalancingLimit, const char*> pmeLoadBalancingLimitNames = {
        "no",
        "box size",
        "domain decomposition",
        "PME grid restriction",
        "maximum allowed grid scaling"
    };
    return pmeLoadBalancingLimitNames[enumValue];
}

struct pme_load_balancing_t
{
    gmx_bool bSepPMERanks;  /**< do we have separate PME ranks? */
    gmx_bool bActive;       /**< is PME tuning active? */
    int64_t  step_rel_stop; /**< stop the tuning after this value of step_rel */
    gmx_bool bTriggerOnDLB; /**< trigger balancing only on DD DLB */
    gmx_bool bBalance;      /**< are we in the balancing phase, i.e. trying different setups? */
    int      nstage;        /**< the current maximum number of stages */
    bool     startupTimeDelayElapsed; /**< Has the c_startupTimeDelay elapsed indicating that the balancing can start. */

    real                     cut_spacing;        /**< the minimum cutoff / PME grid spacing ratio */
    real                     rcut_vdw;           /**< Vdw cutoff (does not change) */
    real                     rcut_coulomb_start; /**< Initial electrostatics cutoff */
    real                     rbufOuter_coulomb;  /**< the outer pairlist buffer size */
    real                     rbufOuter_vdw;      /**< the outer pairlist buffer size */
    real                     rbufInner_coulomb;  /**< the inner pairlist buffer size */
    real                     rbufInner_vdw;      /**< the inner pairlist buffer size */
    matrix                   box_start;          /**< the initial simulation box */
    std::vector<pme_setup_t> setup;              /**< the PME+cutoff setups */
    int                      cur;                /**< the index (in setup) of the current setup */
    int                      fastest;            /**< index of the fastest setup up till now */
    int                      lower_limit;        /**< don't go below this setup index */
    int                      start;    /**< start of setup index range to consider in stage>0 */
    int                      end;      /**< end   of setup index range to consider in stage>0 */
    PmeLoadBalancingLimit    elimited; /**< was the balancing limited, uses enum above */
    CutoffScheme             cutoff_scheme; /**< Verlet or group cut-offs */

    int stage; /**< the current stage */

    int    cycles_n;  /**< step cycle counter cumulative count */
    double cycles_c;  /**< step cycle counter cumulative cycles */
    double startTime; /**< time stamp when the balancing was started on the main rank (relative to the UNIX epoch start).*/
};

/* TODO The code in this file should call this getter, rather than
 * read bActive anywhere */
bool pme_loadbal_is_active(const pme_load_balancing_t* pme_lb)
{
    return pme_lb != nullptr && pme_lb->bActive;
}

// TODO Return a unique_ptr to pme_load_balancing_t
void pme_loadbal_init(pme_load_balancing_t**         pme_lb_p,
                      t_commrec*                     cr,
                      const gmx::MDLogger&           mdlog,
                      const t_inputrec&              ir,
                      const matrix                   box,
                      const interaction_const_t&     ic,
                      const gmx::nonbonded_verlet_t& nbv,
                      gmx_pme_t*                     pmedata,
                      gmx_bool                       bUseGPU)
{

    pme_load_balancing_t* pme_lb;

    // Note that we don't (yet) support PME load balancing with LJ-PME only.
    GMX_RELEASE_ASSERT(usingPme(ir.coulombtype),
                       "pme_loadbal_init called without PME electrostatics");
    // To avoid complexity, we require a single cut-off with PME for q+LJ.
    // This is checked by grompp, but it doesn't hurt to check again.
    GMX_RELEASE_ASSERT(!(usingPme(ir.coulombtype) && usingLJPme(ir.vdwtype) && ir.rcoulomb != ir.rvdw),
                       "With Coulomb and LJ PME, rcoulomb should be equal to rvdw");

    pme_lb = new pme_load_balancing_t;

    pme_lb->bSepPMERanks = !thisRankHasDuty(cr, DUTY_PME);

    /* Initially we turn on balancing directly on based on PP/PME imbalance */
    pme_lb->bTriggerOnDLB = FALSE;

    /* Any number of stages >= 2 is supported */
    pme_lb->nstage = 2;

    pme_lb->cutoff_scheme = ir.cutoff_scheme;

    pme_lb->rbufOuter_coulomb = nbv.pairlistOuterRadius() - ic.rcoulomb;
    pme_lb->rbufOuter_vdw     = nbv.pairlistOuterRadius() - ic.rvdw;
    pme_lb->rbufInner_coulomb = nbv.pairlistInnerRadius() - ic.rcoulomb;
    pme_lb->rbufInner_vdw     = nbv.pairlistInnerRadius() - ic.rvdw;

    /* Scale box with Ewald wall factor; note that we pmedata->boxScaler
     * can't always usedd as it's not available with separate PME ranks.
     */
    EwaldBoxZScaler boxScaler(inputrecPbcXY2Walls(&ir), ir.wall_ewald_zfac);
    boxScaler.scaleBox(box, pme_lb->box_start);

    pme_lb->setup.resize(1);

    pme_lb->rcut_vdw           = ic.rvdw;
    pme_lb->rcut_coulomb_start = ir.rcoulomb;

    pme_lb->cur                    = 0;
    pme_lb->setup[0].rcut_coulomb  = ic.rcoulomb;
    pme_lb->setup[0].rlistOuter    = nbv.pairlistOuterRadius();
    pme_lb->setup[0].rlistInner    = nbv.pairlistInnerRadius();
    pme_lb->setup[0].grid[XX]      = ir.nkx;
    pme_lb->setup[0].grid[YY]      = ir.nky;
    pme_lb->setup[0].grid[ZZ]      = ir.nkz;
    pme_lb->setup[0].ewaldcoeff_q  = ic.ewaldcoeff_q;
    pme_lb->setup[0].ewaldcoeff_lj = ic.ewaldcoeff_lj;

    if (!pme_lb->bSepPMERanks)
    {
        GMX_RELEASE_ASSERT(pmedata, "On ranks doing both PP and PME we need a valid pmedata object");
        pme_lb->setup[0].pmedata = pmedata;
    }

    pme_lb->setup[0].spacing = getGridSpacingFromBox(pme_lb->box_start, pme_lb->setup[0].grid);

    if (ir.fourier_spacing > 0)
    {
        pme_lb->cut_spacing = ir.rcoulomb / ir.fourier_spacing;
    }
    else
    {
        pme_lb->cut_spacing = ir.rcoulomb / pme_lb->setup[0].spacing;
    }

    pme_lb->stage = 0;

    pme_lb->fastest     = 0;
    pme_lb->lower_limit = 0;
    pme_lb->start       = 0;
    pme_lb->end         = 0;
    pme_lb->elimited    = PmeLoadBalancingLimit::No;

    pme_lb->cycles_n = 0;
    pme_lb->cycles_c = 0;
    // only main ranks do timing
    if (!PAR(cr) || (haveDDAtomOrdering(*cr) && DDMAIN(cr->dd)))
    {
        pme_lb->startTime = gmx_gettime();
    }

    if (!wallcycle_have_counter())
    {
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendText(
                        "NOTE: Cycle counters unsupported or not enabled in kernel. Cannot use "
                        "PME-PP balancing.");
    }

    /* Tune with GPUs and/or separate PME ranks.
     * When running only on a CPU without PME ranks, PME tuning will only help
     * with small numbers of atoms in the cut-off sphere.
     */
    pme_lb->bActive = (wallcycle_have_counter() && (bUseGPU || pme_lb->bSepPMERanks));

    /* With GPUs and no separate PME ranks we can't measure the PP/PME
     * imbalance, so we start balancing right away.
     * Otherwise we only start balancing after we observe imbalance.
     */
    pme_lb->bBalance = (pme_lb->bActive && (bUseGPU && !pme_lb->bSepPMERanks));

    pme_lb->step_rel_stop = PMETunePeriod * ir.nstlist;

    /* Delay DD load balancing when GPUs are used */
    if (pme_lb->bActive && haveDDAtomOrdering(*cr) && cr->dd->nnodes > 1 && bUseGPU)
    {
        /* Lock DLB=auto to off (does nothing when DLB=yes/no.
         * With GPUs and separate PME nodes, we want to first
         * do PME tuning without DLB, since DLB might limit
         * the cut-off, which never improves performance.
         * We allow for DLB + PME tuning after a first round of tuning.
         */
        dd_dlb_lock(cr->dd);
        if (dd_dlb_is_locked(cr->dd))
        {
            GMX_LOG(mdlog.warning)
                    .asParagraph()
                    .appendText("NOTE: DLB will not turn on during the first phase of PME tuning");
        }
    }

    *pme_lb_p = pme_lb;
}

/*! \brief Try to increase the cutoff during load balancing */
static gmx_bool pme_loadbal_increase_cutoff(pme_load_balancing_t* pme_lb, int pme_order, const gmx_domdec_t* dd)
{
    real fac, sp;
    real tmpr_coulomb, tmpr_vdw;
    int  d;
    bool grid_ok;

    /* Try to add a new setup with next larger cut-off to the list */
    pme_setup_t set;

    set.pmedata = nullptr;

    NumPmeDomains numPmeDomains = getNumPmeDomains(dd);

    fac = 1;
    do
    {
        /* Avoid infinite while loop, which can occur at the minimum grid size.
         * Note that in practice load balancing will stop before this point.
         * The factor 2.1 allows for the extreme case in which only grids
         * of powers of 2 are allowed (the current code supports more grids).
         */
        if (fac > 2.1)
        {
            return FALSE;
        }

        fac *= 1.01;
        clear_ivec(set.grid);
        sp = calcFftGrid(nullptr,
                         pme_lb->box_start,
                         fac * pme_lb->setup[pme_lb->cur].spacing,
                         minimalPmeGridSize(pme_order),
                         &set.grid[XX],
                         &set.grid[YY],
                         &set.grid[ZZ]);

        /* As here we can't easily check if one of the PME ranks
         * uses threading, we do a conservative grid check.
         * This means we can't use pme_order or less grid lines
         * per PME rank along x, which is not a strong restriction.
         */
        grid_ok = gmx_pme_check_restrictions(
                pme_order, set.grid[XX], set.grid[YY], set.grid[ZZ], numPmeDomains.x, numPmeDomains.y, 0, false, true, false);
    } while (sp <= 1.001 * pme_lb->setup[pme_lb->cur].spacing || !grid_ok);

    set.rcut_coulomb = pme_lb->cut_spacing * sp;
    if (set.rcut_coulomb < pme_lb->rcut_coulomb_start)
    {
        /* This is unlikely, but can happen when e.g. continuing from
         * a checkpoint after equilibration where the box shrank a lot.
         * We want to avoid rcoulomb getting smaller than rvdw
         * and there might be more issues with decreasing rcoulomb.
         */
        set.rcut_coulomb = pme_lb->rcut_coulomb_start;
    }

    if (pme_lb->cutoff_scheme == CutoffScheme::Verlet)
    {
        /* Never decrease the Coulomb and VdW list buffers */
        set.rlistOuter = std::max(set.rcut_coulomb + pme_lb->rbufOuter_coulomb,
                                  pme_lb->rcut_vdw + pme_lb->rbufOuter_vdw);
        set.rlistInner = std::max(set.rcut_coulomb + pme_lb->rbufInner_coulomb,
                                  pme_lb->rcut_vdw + pme_lb->rbufInner_vdw);
    }
    else
    {
        /* TODO Remove these lines and pme_lb->cutoff_scheme */
        tmpr_coulomb = set.rcut_coulomb + pme_lb->rbufOuter_coulomb;
        tmpr_vdw     = pme_lb->rcut_vdw + pme_lb->rbufOuter_vdw;
        /* Two (known) bugs with cutoff-scheme=group here:
         * - This modification of rlist results in incorrect DD comunication.
         * - We should set fr->bTwinRange = (fr->rlistlong > fr->rlist).
         */
        set.rlistOuter = std::min(tmpr_coulomb, tmpr_vdw);
        set.rlistInner = set.rlistOuter;
    }

    set.spacing = sp;
    /* The grid efficiency is the size wrt a grid with uniform x/y/z spacing */
    set.grid_efficiency = 1;
    for (d = 0; d < DIM; d++)
    {
        set.grid_efficiency *= (set.grid[d] * sp) / norm(pme_lb->box_start[d]);
    }
    /* The Ewald coefficient is inversly proportional to the cut-off */
    set.ewaldcoeff_q = pme_lb->setup[0].ewaldcoeff_q * pme_lb->setup[0].rcut_coulomb / set.rcut_coulomb;
    /* We set ewaldcoeff_lj in set, even when LJ-PME is not used */
    set.ewaldcoeff_lj = pme_lb->setup[0].ewaldcoeff_lj * pme_lb->setup[0].rcut_coulomb / set.rcut_coulomb;

    set.count  = 0;
    set.cycles = 0;

    if (debug)
    {
        fprintf(debug,
                "PME loadbal: grid %d %d %d, coulomb cutoff %f\n",
                set.grid[XX],
                set.grid[YY],
                set.grid[ZZ],
                set.rcut_coulomb);
    }
    pme_lb->setup.push_back(set);
    return TRUE;
}

/*! \brief Print the PME grid */
static void print_grid(FILE* fp_err, FILE* fp_log, const char* pre, const char* desc, const pme_setup_t* set, double cycles)
{
    auto buf = gmx::formatString("%-11s%10s pme grid %d %d %d, coulomb cutoff %.3f",
                                 pre,
                                 desc,
                                 set->grid[XX],
                                 set->grid[YY],
                                 set->grid[ZZ],
                                 set->rcut_coulomb);
    if (cycles >= 0)
    {
        buf += gmx::formatString(": %.1f M-cycles", cycles * 1e-6);
    }
    if (fp_err != nullptr)
    {
        fprintf(fp_err, "\r%s\n", buf.c_str());
        fflush(fp_err);
    }
    if (fp_log != nullptr)
    {
        fprintf(fp_log, "%s\n", buf.c_str());
    }
}

/*! \brief Return the index of the last setup used in PME load balancing */
static int pme_loadbal_end(pme_load_balancing_t* pme_lb)
{
    /* In the initial stage only n is set; end is not set yet */
    if (pme_lb->end > 0)
    {
        return pme_lb->end;
    }
    else
    {
        return pme_lb->setup.size();
    }
}

/*! \brief Print descriptive string about what limits PME load balancing */
static void print_loadbal_limited(FILE* fp_err, FILE* fp_log, int64_t step, pme_load_balancing_t* pme_lb)
{
    auto buf = gmx::formatString(
            "step %4s: the %s limits the PME load balancing to a coulomb cut-off of %.3f",
            gmx::int64ToString(step).c_str(),
            enumValueToString(pme_lb->elimited),
            pme_lb->setup[pme_loadbal_end(pme_lb) - 1].rcut_coulomb);
    if (fp_err != nullptr)
    {
        fprintf(fp_err, "\r%s\n", buf.c_str());
        fflush(fp_err);
    }
    if (fp_log != nullptr)
    {
        fprintf(fp_log, "%s\n", buf.c_str());
    }
}

/*! \brief Switch load balancing to stage 1
 *
 * In this stage, only reasonably fast setups are run again. */
static void switch_to_stage1(pme_load_balancing_t* pme_lb)
{
    /* Increase start until we find a setup that is not slower than
     * maxRelativeSlowdownAccepted times the fastest setup.
     */
    pme_lb->start = pme_lb->lower_limit;
    while (pme_lb->start + 1 < gmx::ssize(pme_lb->setup)
           && (pme_lb->setup[pme_lb->start].count == 0
               || pme_lb->setup[pme_lb->start].cycles
                          > pme_lb->setup[pme_lb->fastest].cycles * maxRelativeSlowdownAccepted))
    {
        pme_lb->start++;
    }
    /* While increasing start, we might have skipped setups that we did not
     * time during stage 0. We want to extend the range for stage 1 to include
     * any skipped setups that lie between setups that were measured to be
     * acceptably fast and too slow.
     */
    while (pme_lb->start > pme_lb->lower_limit && pme_lb->setup[pme_lb->start - 1].count == 0)
    {
        pme_lb->start--;
    }

    /* Decrease end only with setups that we timed and that are slow. */
    pme_lb->end = pme_lb->setup.size();
    if (pme_lb->setup[pme_lb->end - 1].count > 0
        && pme_lb->setup[pme_lb->end - 1].cycles
                   > pme_lb->setup[pme_lb->fastest].cycles * maxRelativeSlowdownAccepted)
    {
        pme_lb->end--;
    }

    pme_lb->stage = 1;

    /* Next we want to choose setup pme_lb->end-1, but as we will decrease
     * pme_lb->cur by one right after returning, we set cur to end.
     */
    pme_lb->cur = pme_lb->end;
}

/*! \brief Process the timings and try to adjust the PME grid and Coulomb cut-off
 *
 * The adjustment is done to generate a different non-bonded PP and PME load.
 * With separate PME ranks (PP and PME on different processes) or with
 * a GPU (PP on GPU, PME on CPU), PP and PME run on different resources
 * and changing the load will affect the load balance and performance.
 * The total time for a set of integration steps is monitored and a range
 * of grid/cut-off setups is scanned. After calling pme_load_balance many
 * times and acquiring enough statistics, the best performing setup is chosen.
 * Here we try to take into account fluctuations and changes due to external
 * factors as well as DD load balancing.
 */
static void pme_load_balance(pme_load_balancing_t*          pme_lb,
                             t_commrec*                     cr,
                             FILE*                          fp_err,
                             FILE*                          fp_log,
                             const gmx::MDLogger&           mdlog,
                             const t_inputrec&              ir,
                             const matrix                   box,
                             gmx::ArrayRef<const gmx::RVec> x,
                             double                         cycles,
                             interaction_const_t*           ic,
                             gmx::nonbonded_verlet_t*       nbv,
                             struct gmx_pme_t**             pmedata,
                             int64_t                        step)
{
    gmx_bool     OK;
    pme_setup_t* set;
    double       cycles_fast;
    char         buf[STRLEN], sbuf[22];

    if (PAR(cr))
    {
        gmx_sumd(1, &cycles, cr);
        cycles /= cr->nnodes;
    }

    set = &pme_lb->setup[pme_lb->cur];
    set->count++;

    /* Skip the first c_numPostSwitchTuningIntervalSkip cycles because the first step
     * after a switch is much slower due to allocation and/or caching effects.
     */
    if (set->count % (c_numPostSwitchTuningIntervalSkip + 1) != 0)
    {
        return;
    }

    sprintf(buf, "step %4s: ", gmx_step_str(step, sbuf));
    print_grid(fp_err, fp_log, buf, "timed with", set, cycles);

    GMX_RELEASE_ASSERT(set->count > c_numPostSwitchTuningIntervalSkip, "We should skip cycles");
    if (set->count == (c_numPostSwitchTuningIntervalSkip + 1))
    {
        set->cycles = cycles;
    }
    else
    {
        if (cycles * maxFluctuationAccepted < set->cycles && pme_lb->stage == pme_lb->nstage - 1)
        {
            /* The performance went up a lot (due to e.g. DD load balancing).
             * Add a stage, keep the minima, but rescan all setups.
             */
            pme_lb->nstage++;

            if (debug)
            {
                fprintf(debug,
                        "The performance for grid %d %d %d went from %.3f to %.1f M-cycles, this "
                        "is more than %f\n"
                        "Increased the number stages to %d"
                        " and ignoring the previous performance\n",
                        set->grid[XX],
                        set->grid[YY],
                        set->grid[ZZ],
                        set->cycles * 1e-6,
                        cycles * 1e-6,
                        maxFluctuationAccepted,
                        pme_lb->nstage);
            }
        }
        set->cycles = std::min(set->cycles, cycles);
    }

    if (set->cycles < pme_lb->setup[pme_lb->fastest].cycles)
    {
        pme_lb->fastest = pme_lb->cur;

        if (haveDDAtomOrdering(*cr))
        {
            /* We found a new fastest setting, ensure that with subsequent
             * shorter cut-off's the dynamic load balancing does not make
             * the use of the current cut-off impossible. This solution is
             * a trade-off, as the PME load balancing and DD domain size
             * load balancing can interact in complex ways.
             * With the Verlet kernels, DD load imbalance will usually be
             * mainly due to bonded interaction imbalance, which will often
             * quickly push the domain boundaries beyond the limit for the
             * optimal, PME load balanced, cut-off. But it could be that
             * better overal performance can be obtained with a slightly
             * shorter cut-off and better DD load balancing.
             */
            set_dd_dlb_max_cutoff(cr, pme_lb->setup[pme_lb->fastest].rlistOuter);
        }
    }
    cycles_fast = pme_lb->setup[pme_lb->fastest].cycles;

    /* Check in stage 0 if we should stop scanning grids.
     * Stop when the time is more than maxRelativeSlowDownAccepted longer than the fastest.
     */
    if (pme_lb->stage == 0 && pme_lb->cur > 0
        && cycles > pme_lb->setup[pme_lb->fastest].cycles * maxRelativeSlowdownAccepted)
    {
        pme_lb->setup.resize(pme_lb->cur + 1);
        /* Done with scanning, go to stage 1 */
        switch_to_stage1(pme_lb);
    }

    if (pme_lb->stage == 0)
    {
        int gridsize_start;

        gridsize_start = set->grid[XX] * set->grid[YY] * set->grid[ZZ];

        do
        {
            if (pme_lb->cur + 1 < gmx::ssize(pme_lb->setup))
            {
                /* We had already generated the next setup */
                OK = TRUE;
            }
            else
            {
                /* Find the next setup */
                OK = pme_loadbal_increase_cutoff(pme_lb, ir.pme_order, cr->dd);

                if (!OK)
                {
                    pme_lb->elimited = PmeLoadBalancingLimit::PmeGrid;
                }
            }

            if (OK
                && pme_lb->setup[pme_lb->cur + 1].spacing > c_maxSpacingScaling * pme_lb->setup[0].spacing)
            {
                OK               = FALSE;
                pme_lb->elimited = PmeLoadBalancingLimit::MaxScaling;
            }

            if (OK && ir.pbcType != PbcType::No)
            {
                OK = (gmx::square(pme_lb->setup[pme_lb->cur + 1].rlistOuter)
                      <= max_cutoff2(ir.pbcType, box));
                if (!OK)
                {
                    pme_lb->elimited = PmeLoadBalancingLimit::Box;
                }
            }

            if (OK)
            {
                pme_lb->cur++;

                if (haveDDAtomOrdering(*cr))
                {
                    const bool checkGpuDdLimitation = true;
                    OK                              = change_dd_cutoff(
                            cr, box, x, pme_lb->setup[pme_lb->cur].rlistOuter, checkGpuDdLimitation);
                    if (!OK)
                    {
                        /* Failed: do not use this setup */
                        pme_lb->cur--;
                        pme_lb->elimited = PmeLoadBalancingLimit::DD;
                    }
                }
            }
            if (!OK)
            {
                /* We hit the upper limit for the cut-off,
                 * the setup should not go further than cur.
                 */
                pme_lb->setup.resize(pme_lb->cur + 1);
                print_loadbal_limited(fp_err, fp_log, step, pme_lb);
                /* Switch to the next stage */
                switch_to_stage1(pme_lb);
            }
        } while (OK
                 && !(pme_lb->setup[pme_lb->cur].grid[XX] * pme_lb->setup[pme_lb->cur].grid[YY]
                                      * pme_lb->setup[pme_lb->cur].grid[ZZ]
                              < gridsize_start * gridpointsScaleFactor
                      && pme_lb->setup[pme_lb->cur].grid_efficiency
                                 < pme_lb->setup[pme_lb->cur - 1].grid_efficiency * relativeEfficiencyFactor));
    }

    if (pme_lb->stage > 0 && pme_lb->end == 1)
    {
        pme_lb->cur   = pme_lb->lower_limit;
        pme_lb->stage = pme_lb->nstage;
    }
    else if (pme_lb->stage > 0 && pme_lb->end > 1)
    {
        /* If stage = nstage-1:
         *   scan over all setups, rerunning only those setups
         *   which are not much slower than the fastest
         * else:
         *   use the next setup
         * Note that we loop backward to minimize the risk of the cut-off
         * getting limited by DD DLB, since the DLB cut-off limit is set
         * to the fastest PME setup.
         */
        do
        {
            if (pme_lb->cur > pme_lb->start)
            {
                pme_lb->cur--;
            }
            else
            {
                pme_lb->stage++;

                pme_lb->cur = pme_lb->end - 1;
            }
        } while (pme_lb->stage == pme_lb->nstage - 1 && pme_lb->setup[pme_lb->cur].count > 0
                 && pme_lb->setup[pme_lb->cur].cycles > cycles_fast * maxRelativeSlowdownAccepted);

        if (pme_lb->stage == pme_lb->nstage)
        {
            /* We are done optimizing, use the fastest setup we found */
            pme_lb->cur = pme_lb->fastest;
        }
    }

    if (haveDDAtomOrdering(*cr) && pme_lb->stage > 0)
    {
        const bool checkGpuDdLimitation = true;
        OK = change_dd_cutoff(cr, box, x, pme_lb->setup[pme_lb->cur].rlistOuter, checkGpuDdLimitation);
        if (!OK)
        {
            /* For some reason the chosen cut-off is incompatible with DD.
             * We should continue scanning a more limited range of cut-off's.
             */
            if (pme_lb->cur > 1 && pme_lb->stage == pme_lb->nstage)
            {
                /* stage=nstage says we're finished, but we should continue
                 * balancing, so we set back stage which was just incremented.
                 */
                pme_lb->stage--;
            }
            if (pme_lb->cur <= pme_lb->fastest)
            {
                /* This should not happen, as we set limits on the DLB bounds.
                 * But we implement a complete failsafe solution anyhow.
                 */
                GMX_LOG(mdlog.warning)
                        .asParagraph()
                        .appendTextFormatted(
                                "The fastest PP/PME load balancing setting (cutoff %.3d nm) is no "
                                "longer available due to DD DLB or box size limitations",
                                pme_lb->fastest);
                pme_lb->fastest = pme_lb->lower_limit;
                pme_lb->start   = pme_lb->lower_limit;
            }
            /* Limit the range to below the current cut-off, scan from start */
            pme_lb->end      = pme_lb->cur;
            pme_lb->cur      = pme_lb->start;
            pme_lb->elimited = PmeLoadBalancingLimit::DD;
            print_loadbal_limited(fp_err, fp_log, step, pme_lb);
        }
    }

    /* Change the Coulomb cut-off and the PME grid */

    set = &pme_lb->setup[pme_lb->cur];

    ic->rcoulomb = set->rcut_coulomb;
    nbv->changePairlistRadii(set->rlistOuter, set->rlistInner);
    ic->ewaldcoeff_q = set->ewaldcoeff_q;
    /* TODO: centralize the code that sets the potentials shifts */
    if (ic->coulomb_modifier == InteractionModifiers::PotShift)
    {
        GMX_RELEASE_ASSERT(ic->rcoulomb != 0, "Cutoff radius cannot be zero");
        ic->sh_ewald = std::erfc(ic->ewaldcoeff_q * ic->rcoulomb) / ic->rcoulomb;
    }
    if (usingLJPme(ic->vdwtype))
    {
        /* We have PME for both Coulomb and VdW, set rvdw equal to rcoulomb */
        ic->rvdw          = set->rcut_coulomb;
        ic->ewaldcoeff_lj = set->ewaldcoeff_lj;
        if (ic->vdw_modifier == InteractionModifiers::PotShift)
        {
            real crc2;

            ic->dispersion_shift.cpot = -1.0 / gmx::power6(static_cast<double>(ic->rvdw));
            ic->repulsion_shift.cpot  = -1.0 / gmx::power12(static_cast<double>(ic->rvdw));
            crc2                      = gmx::square(ic->ewaldcoeff_lj * ic->rvdw);
            ic->sh_lj_ewald =
                    (std::exp(-crc2) * (1 + crc2 + 0.5 * crc2 * crc2) - 1) / gmx::power6(ic->rvdw);
        }
    }

    /* We always re-initialize the tables whether they are used or not */
    init_interaction_const_tables(nullptr, ic, set->rlistOuter, ir.tabext);

    gmx::gpu_pme_loadbal_update_param(nbv, *ic);

    if (!pme_lb->bSepPMERanks)
    {
        /* FIXME:
         * CPU PME keeps a list of allocated pmedata's, that's why pme_lb->setup[pme_lb->cur].pmedata is not always nullptr.
         * GPU PME, however, currently needs the gmx_pme_reinit always called on load balancing
         * (pme_gpu_reinit might be not sufficiently decoupled from gmx_pme_init).
         * This can lead to a lot of reallocations for PME GPU.
         * Would be nicer if the allocated grid list was hidden within a single pmedata structure.
         */
        if ((pme_lb->setup[pme_lb->cur].pmedata == nullptr)
            || pme_gpu_task_enabled(pme_lb->setup[pme_lb->cur].pmedata))
        {
            gmx_pme_t* newPmeData;
            // Generate a new PME data structure, copying part of the old pointers.
            gmx_pme_reinit(
                    &newPmeData, cr, pme_lb->setup[0].pmedata, &ir, set->grid, set->ewaldcoeff_q, set->ewaldcoeff_lj);
            // Destroy the old structure. Must be done after gmx_pme_reinit in case pme_lb->cur is 0.
            if (set->pmedata != nullptr)
            {
                gmx_pme_destroy(set->pmedata, false);
            }
            set->pmedata = newPmeData;
        }
        *pmedata = set->pmedata;
    }
    else
    {
        /* Tell our PME-only rank to switch grid */
        gmx_pme_send_switchgrid(cr, set->grid, set->ewaldcoeff_q, set->ewaldcoeff_lj);
    }

    if (debug)
    {
        print_grid(nullptr, debug, "", "switched to", set, -1);
    }

    if (pme_lb->stage == pme_lb->nstage)
    {
        print_grid(fp_err, fp_log, "", "optimal", set, -1);
    }
}

/*! \brief Prepare for another round of PME load balancing
 *
 * \param[in,out] pme_lb       Pointer to PME load balancing struct
 * \param[in]     bDlbUnlocked TRUE is DLB was locked and is now unlocked
 *
 * If the conditions (e.g. DLB off/on, CPU/GPU throttling etc.) changed,
 * the PP/PME balance might change and re-balancing can improve performance.
 * This function adds 2 stages and adjusts the considered setup range.
 */
static void continue_pme_loadbal(pme_load_balancing_t* pme_lb, gmx_bool bDlbUnlocked)
{
    /* Add 2 tuning stages, keep the detected end of the setup range */
    pme_lb->nstage += 2;
    if (bDlbUnlocked && pme_lb->bSepPMERanks)
    {
        /* With separate PME ranks, DLB should always lower the PP load and
         * can only increase the PME load (more communication and imbalance),
         * so we only need to scan longer cut-off's.
         */
        pme_lb->lower_limit = pme_lb->cur;
    }
    pme_lb->start = pme_lb->lower_limit;
}

void pme_loadbal_do(pme_load_balancing_t*          pme_lb,
                    t_commrec*                     cr,
                    FILE*                          fp_err,
                    FILE*                          fp_log,
                    const gmx::MDLogger&           mdlog,
                    const t_inputrec&              ir,
                    t_forcerec*                    fr,
                    const matrix                   box,
                    gmx::ArrayRef<const gmx::RVec> x,
                    gmx_wallcycle*                 wcycle,
                    int64_t                        step,
                    int64_t                        step_rel,
                    gmx_bool*                      bPrinting,
                    bool                           useGpuPmePpCommunication)
{
    int    n_prev;
    double cycles_prev;

    assert(pme_lb != nullptr);

    if (!pme_lb->bActive)
    {
        return;
    }

    n_prev      = pme_lb->cycles_n;
    cycles_prev = pme_lb->cycles_c;
    wallcycle_get(wcycle, WallCycleCounter::Step, &pme_lb->cycles_n, &pme_lb->cycles_c);

    /* Before the first step we haven't done any steps yet.
     * Also handle cases where ir.init_step % ir.nstlist != 0.
     * We also want to skip a number of steps and seconds while
     * the CPU and GPU, when used, performance stabilizes.
     */
    if (!PAR(cr) || (haveDDAtomOrdering(*cr) && DDMAIN(cr->dd)))
    {
        pme_lb->startupTimeDelayElapsed = (gmx_gettime() - pme_lb->startTime < c_startupTimeDelay);
    }
    if (haveDDAtomOrdering(*cr))
    {
        dd_bcast(cr->dd, sizeof(bool), &pme_lb->startupTimeDelayElapsed);
    }

    if (pme_lb->cycles_n == 0 || step_rel < c_numFirstTuningIntervalSkip * ir.nstlist
        || pme_lb->startupTimeDelayElapsed)
    {
        *bPrinting = FALSE;
        return;
    }
    /* Sanity check, we expect nstlist cycle counts */
    if (pme_lb->cycles_n - n_prev != ir.nstlist)
    {
        /* We could return here, but it's safer to issue an error and quit */
        gmx_incons("pme_loadbal_do called at an interval != nstlist");
    }

    /* PME grid + cut-off optimization with GPUs or PME ranks */
    if (!pme_lb->bBalance && pme_lb->bSepPMERanks)
    {
        if (pme_lb->bTriggerOnDLB)
        {
            pme_lb->bBalance = dd_dlb_is_on(cr->dd);
        }
        /* We should ignore the first timing to avoid timing allocation
         * overhead. And since the PME load balancing is called just
         * before DD repartitioning, the ratio returned by dd_pme_f_ratio
         * is not over the last nstlist steps, but the nstlist steps before
         * that. So the first useful ratio is available at step_rel=3*nstlist.
         */
        else if (step_rel >= c_numFirstTuningIntervalSkipWithSepPme * ir.nstlist)
        {
            GMX_ASSERT(haveDDAtomOrdering(*cr), "Domain decomposition should be active here");
            if (DDMAIN(cr->dd))
            {
                /* If PME rank load is too high, start tuning. If
                   PME-PP direct GPU communication is active,
                   unconditionally start tuning since ratio will be
                   unreliable due to CPU-GPU asynchronicity in codepath */
                pme_lb->bBalance = useGpuPmePpCommunication
                                           ? true
                                           : (dd_pme_f_ratio(cr->dd) >= loadBalanceTriggerFactor);
            }
            dd_bcast(cr->dd, sizeof(gmx_bool), &pme_lb->bBalance);
        }

        pme_lb->bActive = (pme_lb->bBalance || step_rel <= pme_lb->step_rel_stop);
    }

    /* The location in the code of this balancing termination is strange.
     * You would expect to have it after the call to pme_load_balance()
     * below, since there pme_lb->stage is updated.
     * But when terminating directly after deciding on and selecting the
     * optimal setup, DLB will turn on right away if it was locked before.
     * This might be due to PME reinitialization. So we check stage here
     * to allow for another nstlist steps with DLB locked to stabilize
     * the performance.
     */
    if (pme_lb->bBalance && pme_lb->stage == pme_lb->nstage)
    {
        pme_lb->bBalance = FALSE;

        if (haveDDAtomOrdering(*cr) && dd_dlb_is_locked(cr->dd))
        {
            /* Unlock the DLB=auto, DLB is allowed to activate */
            dd_dlb_unlock(cr->dd);
            GMX_LOG(mdlog.warning)
                    .asParagraph()
                    .appendText("NOTE: DLB can now turn on, when beneficial");

            /* We don't deactivate the tuning yet, since we will balance again
             * after DLB gets turned on, if it does within PMETune_period.
             */
            continue_pme_loadbal(pme_lb, TRUE);
            pme_lb->bTriggerOnDLB = TRUE;
            pme_lb->step_rel_stop = step_rel + PMETunePeriod * ir.nstlist;
        }
        else
        {
            /* We're completely done with PME tuning */
            pme_lb->bActive = FALSE;
        }

        if (haveDDAtomOrdering(*cr))
        {
            /* Set the cut-off limit to the final selected cut-off,
             * so we don't have artificial DLB limits.
             * This also ensures that we won't disable the currently
             * optimal setting during a second round of PME balancing.
             */
            set_dd_dlb_max_cutoff(cr, fr->nbv->pairlistOuterRadius());
        }
    }

    if (pme_lb->bBalance)
    {
        /* We might not have collected nstlist steps in cycles yet,
         * since init_step might not be a multiple of nstlist,
         * but the first data collected is skipped anyhow.
         */
        pme_load_balance(pme_lb,
                         cr,
                         fp_err,
                         fp_log,
                         mdlog,
                         ir,
                         box,
                         x,
                         pme_lb->cycles_c - cycles_prev,
                         fr->ic.get(),
                         fr->nbv.get(),
                         &fr->pmedata,
                         step);

        /* Update deprecated rlist in forcerec to stay in sync with fr->nbv */
        fr->rlist = fr->nbv->pairlistOuterRadius();

        if (ir.eDispCorr != DispersionCorrectionType::No)
        {
            fr->dispersionCorrection->setParameters(*fr->ic);
        }
    }

    if (!pme_lb->bBalance && (!pme_lb->bSepPMERanks || step_rel > pme_lb->step_rel_stop))
    {
        /* We have just deactivated the balancing and we're not measuring PP/PME
         * imbalance during the first steps of the run: deactivate the tuning.
         */
        pme_lb->bActive = FALSE;
    }

    if (!(pme_lb->bActive) && haveDDAtomOrdering(*cr) && dd_dlb_is_locked(cr->dd))
    {
        /* Make sure DLB is allowed when we deactivate PME tuning */
        dd_dlb_unlock(cr->dd);
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendText("NOTE: DLB can now turn on, when beneficial");
    }

    *bPrinting = pme_lb->bBalance;
}

/*! \brief Return product of the number of PME grid points in each dimension */
static int pme_grid_points(const pme_setup_t* setup)
{
    return setup->grid[XX] * setup->grid[YY] * setup->grid[ZZ];
}

/*! \brief Print one load-balancing setting */
static void print_pme_loadbal_setting(FILE* fplog, const char* name, const pme_setup_t* setup)
{
    fprintf(fplog,
            "   %-7s %6.3f nm %6.3f nm     %3d %3d %3d   %5.3f nm  %5.3f nm\n",
            name,
            setup->rcut_coulomb,
            setup->rlistInner,
            setup->grid[XX],
            setup->grid[YY],
            setup->grid[ZZ],
            setup->spacing,
            1 / setup->ewaldcoeff_q);
}

/*! \brief Print all load-balancing settings */
static void print_pme_loadbal_settings(pme_load_balancing_t* pme_lb,
                                       FILE*                 fplog,
                                       const gmx::MDLogger&  mdlog,
                                       gmx_bool              bNonBondedOnGPU)
{
    double pp_ratio, grid_ratio;
    real   pp_ratio_temporary;

    pp_ratio_temporary = pme_lb->setup[pme_lb->cur].rlistInner / pme_lb->setup[0].rlistInner;
    pp_ratio           = gmx::power3(pp_ratio_temporary);
    grid_ratio         = pme_grid_points(&pme_lb->setup[pme_lb->cur])
                 / static_cast<double>(pme_grid_points(&pme_lb->setup[0]));

    fprintf(fplog, "\n");
    fprintf(fplog, "       P P   -   P M E   L O A D   B A L A N C I N G\n");
    fprintf(fplog, "\n");
    /* Here we only warn when the optimal setting is the last one */
    if (pme_lb->elimited != PmeLoadBalancingLimit::No && pme_lb->cur == pme_loadbal_end(pme_lb) - 1)
    {
        fprintf(fplog,
                " NOTE: The PP/PME load balancing was limited by the %s,\n",
                enumValueToString(pme_lb->elimited));
        fprintf(fplog, "       you might not have reached a good load balance.\n");
        if (pme_lb->elimited == PmeLoadBalancingLimit::DD)
        {
            fprintf(fplog, "       Try different mdrun -dd settings or lower the -dds value.\n");
        }
        fprintf(fplog, "\n");
    }
    fprintf(fplog, " PP/PME load balancing changed the cut-off and PME settings:\n");
    fprintf(fplog, "           particle-particle                    PME\n");
    fprintf(fplog, "            rcoulomb  rlist            grid      spacing   1/beta\n");
    print_pme_loadbal_setting(fplog, "initial", &pme_lb->setup[0]);
    print_pme_loadbal_setting(fplog, "final", &pme_lb->setup[pme_lb->cur]);
    fprintf(fplog, " cost-ratio           %4.2f             %4.2f\n", pp_ratio, grid_ratio);
    fprintf(fplog, " (note that these numbers concern only part of the total PP and PME load)\n");

    if (pp_ratio > 1.5 && !bNonBondedOnGPU)
    {
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendText(
                        "NOTE: PME load balancing increased the non-bonded workload by more than "
                        "50%.\n"
                        "      For better performance, use (more) PME ranks (mdrun -npme),\n"
                        "      or if you are beyond the scaling limit, use fewer total ranks (or "
                        "nodes).");
    }
    else
    {
        fprintf(fplog, "\n");
    }
}

void pme_loadbal_done(pme_load_balancing_t* pme_lb, FILE* fplog, const gmx::MDLogger& mdlog, gmx_bool bNonBondedOnGPU)
{
    if (fplog != nullptr && (pme_lb->cur > 0 || pme_lb->elimited != PmeLoadBalancingLimit::No))
    {
        print_pme_loadbal_settings(pme_lb, fplog, mdlog, bNonBondedOnGPU);
    }
    for (int i = 0; i < gmx::ssize(pme_lb->setup); i++)
    {
        // current element is stored in forcerec and free'd in Mdrunner::mdruner, together with shared data
        if (i != pme_lb->cur)
        {
            gmx_pme_destroy(pme_lb->setup[i].pmedata, false);
        }
    }

    delete pme_lb;
}
