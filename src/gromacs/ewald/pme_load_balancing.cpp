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
 * \brief This file implements the PmeLoadBalancing class
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
#include "gromacs/mdlib/dispersioncorrection.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/simulation_workload.h"
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
#include "gromacs/utility/mpicomm.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/vec.h"

#include "pme_internal.h"
#include "pme_pp.h"

namespace gmx
{

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
    constexpr EnumerationArray<PmeLoadBalancingLimit, const char*> pmeLoadBalancingLimitNames = {
        "no",
        "box size",
        "domain decomposition",
        "PME grid restriction",
        "maximum allowed grid scaling"
    };
    return pmeLoadBalancingLimitNames[enumValue];
}

//! Struct for holding information on initial cut-off distances and the box
struct Cutoffs
{
    real   cutoffDivSpacing;   /**< the minimum cutoff / PME grid spacing ratio */
    real   rcut_vdw;           /**< Vdw cutoff (does not change) */
    real   rcut_coulomb_start; /**< Initial electrostatics cutoff */
    real   rbufOuter_coulomb;  /**< the outer pairlist buffer size */
    real   rbufOuter_vdw;      /**< the outer pairlist buffer size */
    real   rbufInner_coulomb;  /**< the inner pairlist buffer size */
    real   rbufInner_vdw;      /**< the inner pairlist buffer size */
    matrix startBox;           /**< the initial simulation box */
};

/*! \brief Impl class for PmeLoadBalancing
 *
 * The states the algorithm can reside in:
 *
 * - isActive_ == false: inactive, either from the beginning or after tuning has finished
 *
 * - isActive_ == true:
 *   - isInBalancingPhase == false: only timing and checking for imbalance,
 *                                  sufficient imbalance will trigger the balancing phase
 *   - isInBalancingPhase == true:
 *     - stage == 0: generate and step coarsely through a range of grid to get an initial grid range
 *                   and rough timings
 *     - stage > 0:  time grids within a sub-range of interest, this sub-range can be extended
 *                   when the timings change a lot from the initial ones
 *
 *   Will, permanently, move to isActive_ = false when no imbalance is observed over many steps
 *   or when the balancing has finished.
 */
class PmeLoadBalancing::Impl
{
public:
    Impl(gmx_domdec_t*              dd,
         const MDLogger&            mdlog,
         const t_inputrec&          ir,
         const matrix               box,
         const interaction_const_t& ic,
         const nonbonded_verlet_t&  nbv,
         gmx_pme_t*                 pmedata,
         const SimulationWorkload&  simulationWork);

    ~Impl();

    bool isActive() const { return isActive_; }

    bool isPrintingLoad() const { return isInBalancingPhase_; }

    //! Returns one index past the last setup that will be covered during balancing
    int setupEnd() const
    {
        /* In the initial stage only n is set; end is not set yet */
        if (endSetup_ > 0)
        {
            return endSetup_;
        }
        else
        {
            return setups_.size();
        }
    }

    //! Attempts to increase the cutoff, returns true when a setup has been added to the list
    bool increaseCutoff();

    /*! \brief Switch load balancing to stage 1
     *
     * In this stage, only sufficiently fast setups are run again.
     */
    void switchToStage1();

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
    void balance(FILE*                fp_err,
                 const matrix         box,
                 ArrayRef<const RVec> x,
                 double               cycles,
                 interaction_const_t* ic,
                 nonbonded_verlet_t*  nbv,
                 gmx_pme_t**          pmedata,
                 int64_t              step);

    /*! \brief Prepare for another round of PME load balancing
     *
     * \param[in]     dlbWasUnlocked  Pass true when DLB was locked and is now unlocked
     *
     * If the conditions (e.g. DLB off/on, CPU/GPU throttling etc.) changed,
     * the PP/PME balance might change and re-balancing can improve performance.
     * This function adds 2 stages and adjusts the considered setup range.
     */
    void addTwoStages(bool dlbWasUnlocked);

    void addCycles(FILE*                fp_err,
                   t_forcerec*          fr,
                   const matrix         box,
                   ArrayRef<const RVec> x,
                   const gmx_wallcycle* wcycle,
                   int64_t              step,
                   int64_t              step_rel);

    void printSettings() const;

private:
    const bool haveSepPMERanks_;          /**< do we have separate PME ranks? */
    const bool useGpuForNonbondeds_;      /**< do we use a GPU for, at least, the non-bondes? */
    const bool useGpuPmePpCommunication_; /**< do we perform PME-PP communication on GPUs? */

    bool    isActive_;        /**< is PME tuning active? */
    int64_t stepRelStop_;     /**< stop the tuning after this value of step_rel */
    bool    triggerOnDLB_;    /**< trigger balancing only on DD DLB */
    bool isInBalancingPhase_; /**< are we in the balancing phase, i.e. trying different setups? */
    int  numStages_;          /**< the current maximum number of stages */
    bool startupTimeDelayElapsed_; /**< Has the c_startupTimeDelay elapsed indicating that the balancing can start. */

    const Cutoffs cutoffs_; /**< Information on and settings for cutoffs */

    const t_inputrec& ir_; /**< Reference to the input record */

    std::vector<pme_setup_t> setups_;       /**< the PME+cutoff setups */
    int                      currentSetup_; /**< the index (in setup) of the current setup */
    int                      fastestSetup_; /**< index of the fastest setup up till now */
    int                      lowerLimit_; /**< don't go below this setup index (seems not needed) */
    int                      startSetup_; /**< start of setup index range to consider in stage>0 */
    int                      endSetup_;   /**< end   of setup index range to consider in stage>0 */
    PmeLoadBalancingLimit    limited_;    /**< was the balancing limited, uses enum above */

    int stage_; /**< the current stage */

    int    cyclesNumber_;  /**< step cycle counter cumulative count */
    double cyclesCounter_; /**< step cycle counter cumulative cycles */
    double startTime_; /**< time stamp when the balancing was started on the main rank (relative to the UNIX epoch start).*/

    gmx_domdec_t& dd_; /**< Reference to the domain decomposition object */

    const MDLogger& mdlog_; /**< Reference to the mdlogger */
};

//! Returns whether it might be useful to do PME tuning
static bool pmeTuningIsUseful(const SimulationWorkload& simulationWork)
{
    /* Tune with GPUs and/or separate PME ranks.
     * When running only on a CPU without PME ranks, PME tuning will only help
     * with small numbers of atoms in the cut-off sphere.
     * Disable PME tuning with GPU PME decomposition.
     */
    return wallcycle_have_counter()
           && (simulationWork.useCpuNonbonded || simulationWork.haveSeparatePmeRank)
           && !simulationWork.useGpuPmeDecomposition;
}

//! Computes and returns initially set cutoffs and box
static Cutoffs getCutoffs(const t_inputrec&          ir,
                          const matrix               box,
                          const interaction_const_t& ic,
                          const nonbonded_verlet_t&  nbv)
{
    Cutoffs cutoffs;

    cutoffs.rbufOuter_coulomb = nbv.pairlistOuterRadius() - ic.coulomb.cutoff;
    cutoffs.rbufOuter_vdw     = nbv.pairlistOuterRadius() - ic.vdw.cutoff;
    cutoffs.rbufInner_coulomb = nbv.pairlistInnerRadius() - ic.coulomb.cutoff;
    cutoffs.rbufInner_vdw     = nbv.pairlistInnerRadius() - ic.vdw.cutoff;

    /* Scale box with Ewald wall factor; note that we pmedata->boxScaler
     * can't always usedd as it's not available with separate PME ranks.
     */
    EwaldBoxZScaler boxScaler(inputrecPbcXY2Walls(&ir), ir.wall_ewald_zfac);
    boxScaler.scaleBox(box, cutoffs.startBox);

    cutoffs.rcut_vdw           = ic.vdw.cutoff;
    cutoffs.rcut_coulomb_start = ir.rcoulomb;

    real fourierSpacing = ir.fourier_spacing;
    if (fourierSpacing == 0)
    {
        fourierSpacing = getGridSpacingFromBox(cutoffs.startBox, IVec(ir.nkx, ir.nky, ir.nkz));
    }
    cutoffs.cutoffDivSpacing = ir.rcoulomb / fourierSpacing;

    return cutoffs;
}

PmeLoadBalancing::Impl::Impl(gmx_domdec_t*              dd,
                             const MDLogger&            mdlog,
                             const t_inputrec&          ir,
                             const matrix               box,
                             const interaction_const_t& ic,
                             const nonbonded_verlet_t&  nbv,
                             gmx_pme_t*                 pmedata,
                             const SimulationWorkload&  simulationWork) :
    haveSepPMERanks_(simulationWork.haveSeparatePmeRank),
    useGpuForNonbondeds_(simulationWork.useGpuNonbonded),
    useGpuPmePpCommunication_(simulationWork.useGpuPmePpCommunication),
    isActive_(pmeTuningIsUseful(simulationWork)),
    cutoffs_(getCutoffs(ir, box, ic, nbv)),
    ir_(ir),
    dd_(*dd),
    mdlog_(mdlog)
{
    // Note that we don't (yet) support PME load balancing with LJ-PME only.
    GMX_RELEASE_ASSERT(usingPme(ir.coulombtype),
                       "pme_loadbal_init called without PME electrostatics");
    // To avoid complexity, we require a single cut-off with PME for q+LJ.
    // This is checked by grompp, but it doesn't hurt to check again.
    GMX_RELEASE_ASSERT(!(usingPme(ir.coulombtype) && usingLJPme(ir.vdwtype) && ir.rcoulomb != ir.rvdw),
                       "With Coulomb and LJ PME, rcoulomb should be equal to rvdw");

    /* Initially we turn on balancing directly on based on PP/PME imbalance */
    triggerOnDLB_ = false;

    /* Any number of stages >= 2 is supported */
    numStages_ = 2;

    setups_.resize(1);

    currentSetup_            = 0;
    setups_[0].rcut_coulomb  = ic.coulomb.cutoff;
    setups_[0].rlistOuter    = nbv.pairlistOuterRadius();
    setups_[0].rlistInner    = nbv.pairlistInnerRadius();
    setups_[0].grid[XX]      = ir.nkx;
    setups_[0].grid[YY]      = ir.nky;
    setups_[0].grid[ZZ]      = ir.nkz;
    setups_[0].ewaldcoeff_q  = ic.coulomb.ewaldCoeff;
    setups_[0].ewaldcoeff_lj = ic.vdw.ewaldCoeff;

    if (!haveSepPMERanks_)
    {
        GMX_RELEASE_ASSERT(pmedata, "On ranks doing both PP and PME we need a valid pmedata object");
        setups_[0].pmedata = pmedata;
    }

    setups_[0].spacing = getGridSpacingFromBox(cutoffs_.startBox, setups_[0].grid);

    stage_ = 0;

    fastestSetup_ = 0;
    lowerLimit_   = 0;
    startSetup_   = 0;
    endSetup_     = 0;
    limited_      = PmeLoadBalancingLimit::No;

    cyclesNumber_  = 0;
    cyclesCounter_ = 0;
    // only main ranks do timing
    if (pmedata == nullptr || !pmedata->simulationIsParallel
        || (pmedata->haveDDAtomOrdering && DDMAIN(dd)))
    {
        startTime_ = gmx_gettime();
    }

    if (!wallcycle_have_counter())
    {
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendText(
                        "NOTE: Cycle counters unsupported or not enabled in kernel. Cannot use "
                        "PME-PP balancing.");
    }

    /* With GPUs and no separate PME ranks we can't measure the PP/PME
     * imbalance, so we start balancing right away.
     * Otherwise we only start balancing after we observe imbalance.
     */
    isInBalancingPhase_ = (isActive_ && (useGpuForNonbondeds_ && !haveSepPMERanks_));

    stepRelStop_ = PMETunePeriod * ir_.nstlist;

    /* Delay DD load balancing when GPUs are used */
    if (isActive_ && dd != nullptr && dd->nnodes > 1 && useGpuForNonbondeds_)
    {
        /* Lock DLB=auto to off (does nothing when DLB=yes/no.
         * With GPUs and separate PME nodes, we want to first
         * do PME tuning without DLB, since DLB might limit
         * the cut-off, which never improves performance.
         * We allow for DLB + PME tuning after a first round of tuning.
         */
        dd_dlb_lock(dd);
        if (dd_dlb_is_locked(dd))
        {
            GMX_LOG(mdlog.warning)
                    .asParagraph()
                    .appendText("NOTE: DLB will not turn on during the first phase of PME tuning");
        }
    }
}

PmeLoadBalancing::Impl::~Impl()
{
    for (int i = 0; i < gmx::ssize(setups_); i++)
    {
        // current element is stored in forcerec and free'd in Mdrunner::mdrunner, together with shared data
        if (i != currentSetup_)
        {
            gmx_pme_destroy(setups_[i].pmedata, false);
        }
    }
}

/*! \brief Return product of the number of PME grid points in each dimension */
static int numPmeGridPoints(const pme_setup_t& setup)
{
    return setup.grid[XX] * setup.grid[YY] * setup.grid[ZZ];
}

bool PmeLoadBalancing::Impl::increaseCutoff()
{
    /* Try to add a new setup with next larger cut-off to the list */
    pme_setup_t set;

    set.pmedata = nullptr;

    const NumPmeDomains numPmeDomains = getNumPmeDomains(&dd_);

    const pme_setup_t& currentSetup = setups_[currentSetup_];

    real fac = 1;
    real sp;
    bool grid_ok;
    do
    {
        /* Avoid infinite while loop, which can occur at the minimum grid size.
         * Note that in practice load balancing will stop before this point.
         * The factor 2.1 allows for the extreme case in which only grids
         * of powers of 2 are allowed (the current code supports more grids).
         */
        if (fac > 2.1)
        {
            return false;
        }

        fac *= 1.01;
        clear_ivec(set.grid);
        sp = calcFftGrid(nullptr,
                         cutoffs_.startBox,
                         fac * currentSetup.spacing,
                         minimalPmeGridSize(ir_.pme_order),
                         &set.grid[XX],
                         &set.grid[YY],
                         &set.grid[ZZ]);

        /* As here we can't easily check if one of the PME ranks
         * uses threading, we do a conservative grid check.
         * This means we can't use PME order or less grid lines
         * per PME rank along x, which is not a strong restriction.
         */
        grid_ok = gmx_pme_check_restrictions(ir_.pme_order,
                                             set.grid[XX],
                                             set.grid[YY],
                                             set.grid[ZZ],
                                             numPmeDomains.x,
                                             numPmeDomains.y,
                                             0,
                                             false,
                                             true,
                                             false);
    } while (sp <= 1.001 * currentSetup.spacing || !grid_ok);

    set.rcut_coulomb = cutoffs_.cutoffDivSpacing * sp;
    if (set.rcut_coulomb < cutoffs_.rcut_coulomb_start)
    {
        /* This is unlikely, but can happen when e.g. continuing from
         * a checkpoint after equilibration where the box shrank a lot.
         * We want to avoid rcoulomb getting smaller than rvdw
         * and there might be more issues with decreasing rcoulomb.
         */
        set.rcut_coulomb = cutoffs_.rcut_coulomb_start;
    }

    /* Never decrease the Coulomb and VdW list buffers */
    set.rlistOuter = std::max(set.rcut_coulomb + cutoffs_.rbufOuter_coulomb,
                              cutoffs_.rcut_vdw + cutoffs_.rbufOuter_vdw);
    set.rlistInner = std::max(set.rcut_coulomb + cutoffs_.rbufInner_coulomb,
                              cutoffs_.rcut_vdw + cutoffs_.rbufInner_vdw);

    set.spacing = sp;
    /* The grid efficiency is the size wrt a grid with uniform x/y/z spacing */
    set.grid_efficiency = 1;
    for (int d = 0; d < DIM; d++)
    {
        set.grid_efficiency *= (set.grid[d] * sp) / norm(cutoffs_.startBox[d]);
    }
    /* The Ewald coefficient is inversly proportional to the cut-off */
    set.ewaldcoeff_q = setups_[0].ewaldcoeff_q * setups_[0].rcut_coulomb / set.rcut_coulomb;
    /* We set ewaldcoeff_lj in set, even when LJ-PME is not used */
    set.ewaldcoeff_lj = setups_[0].ewaldcoeff_lj * setups_[0].rcut_coulomb / set.rcut_coulomb;

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

    setups_.push_back(set);

    return true;
}

/*! \brief Print the PME grid */
static void printGrid(FILE*              fp_err,
                      const MDLogger&    mdlog,
                      const char*        pre,
                      const char*        desc,
                      const pme_setup_t& set,
                      double             cycles)
{
    auto buf = formatString("%-11s%10s pme grid %d %d %d, coulomb cutoff %.3f",
                            pre,
                            desc,
                            set.grid[XX],
                            set.grid[YY],
                            set.grid[ZZ],
                            set.rcut_coulomb);
    if (cycles >= 0)
    {
        buf += formatString(": %.1f M-cycles", cycles * 1e-6);
    }
    if (fp_err != nullptr)
    {
        fprintf(fp_err, "\r%s\n", buf.c_str());
        std::fflush(fp_err);
    }
    GMX_LOG(mdlog.info).appendText(buf);
}

/*! \brief Print descriptive string about what limits PME load balancing */
static void printLoadBalLimited(FILE*                       fp_err,
                                const MDLogger&             mdlog,
                                int64_t                     step,
                                const PmeLoadBalancingLimit limited,
                                const pme_setup_t&          setup)
{
    auto buf = formatString(
            "step %4s: the %s limits the PME load balancing to a coulomb cut-off of %.3f",
            int64ToString(step).c_str(),
            enumValueToString(limited),
            setup.rcut_coulomb);
    if (fp_err != nullptr)
    {
        fprintf(fp_err, "\r%s\n", buf.c_str());
        std::fflush(fp_err);
    }
    GMX_LOG(mdlog.info).appendText(buf);
}

void PmeLoadBalancing::Impl::switchToStage1()
{
    /* Increase start until we find a setup that is not slower than
     * maxRelativeSlowdownAccepted times the fastest setup.
     */
    startSetup_ = lowerLimit_;
    while (startSetup_ + 1 < gmx::ssize(setups_)
           && (setups_[startSetup_].count == 0
               || setups_[startSetup_].cycles > setups_[fastestSetup_].cycles * maxRelativeSlowdownAccepted))
    {
        startSetup_++;
    }
    /* While increasing start, we might have skipped setups that we did not
     * time during stage 0. We want to extend the range for stage 1 to include
     * any skipped setups that lie between setups that were measured to be
     * acceptably fast and too slow.
     */
    while (startSetup_ > lowerLimit_ && setups_[startSetup_].count == 0)
    {
        startSetup_--;
    }

    /* Decrease end only with setups that we timed and that are slow. */
    endSetup_ = setups_.size();
    if (setups_[endSetup_ - 1].count > 0
        && setups_[endSetup_ - 1].cycles > setups_[fastestSetup_].cycles * maxRelativeSlowdownAccepted)
    {
        endSetup_--;
    }

    stage_ = 1;

    /* Next we want to choose setup endSetup_-1, but as we will decrease
     * currentSetup_ by one right after returning, we set cur to endSetup_.
     */
    currentSetup_ = endSetup_;
}

//! Updates all the mdrun machinery for \p setup, setup->pmedata might be updated
static void applySetup(pme_setup_t*         setup,
                       gmx_pme_t*           pmedataOfSetup0,
                       const t_inputrec&    ir,
                       interaction_const_t* ic,
                       nonbonded_verlet_t*  nbv,
                       gmx_domdec_t*        dd)
{
    ic->coulomb.cutoff = setup->rcut_coulomb;
    nbv->changePairlistRadii(setup->rlistOuter, setup->rlistInner);
    ic->coulomb.ewaldCoeff = setup->ewaldcoeff_q;
    /* TODO: centralize the code that sets the potentials shifts */
    if (ic->coulomb.modifier == InteractionModifiers::PotShift)
    {
        GMX_RELEASE_ASSERT(ic->coulomb.cutoff != 0, "Cutoff radius cannot be zero");
        ic->coulomb.ewaldShift =
                std::erfc(ic->coulomb.ewaldCoeff * ic->coulomb.cutoff) / ic->coulomb.cutoff;
    }
    if (usingLJPme(ic->vdw.type))
    {
        /* We have PME for both Coulomb and VdW, set rvdw equal to rcoulomb */
        ic->vdw.cutoff     = setup->rcut_coulomb;
        ic->vdw.ewaldCoeff = setup->ewaldcoeff_lj;
        if (ic->vdw.modifier == InteractionModifiers::PotShift)
        {
            ic->vdw.dispersionShift.cpot = -1.0 / power6(static_cast<double>(ic->vdw.cutoff));
            ic->vdw.repulsionShift.cpot  = -1.0 / power12(static_cast<double>(ic->vdw.cutoff));
            real crc2                    = square(ic->vdw.ewaldCoeff * ic->vdw.cutoff);
            ic->vdw.ewaldShift =
                    (std::exp(-crc2) * (1 + crc2 + 0.5 * crc2 * crc2) - 1) / power6(ic->vdw.cutoff);
        }
    }

    /* We always re-initialize the tables whether they are used or not */
    init_interaction_const_tables(nullptr, ic, setup->rlistOuter, ir.tabext);

    gpu_pme_loadbal_update_param(nbv, *ic);

    if (dd->hasPmeDuty)
    {
        /* FIXME:
         * CPU PME keeps a list of allocated pmedata's, that's why setups_[currentSetup_].pmedata is not always nullptr.
         * GPU PME, however, currently needs the gmx_pme_reinit always called on load balancing
         * (pme_gpu_reinit might be not sufficiently decoupled from gmx_pme_init).
         * This can lead to a lot of reallocations for PME GPU.
         * Would be nicer if the allocated grid list was hidden within a single pmedata structure.
         */
        if (setup->pmedata == nullptr || pme_gpu_task_enabled(setup->pmedata))
        {
            gmx_pme_t* newPmeData;
            // Generate a new PME data structure, copying part of the old pointers.
            gmx_pme_reinit(
                    &newPmeData, dd, pmedataOfSetup0, &ir, setup->grid, setup->ewaldcoeff_q, setup->ewaldcoeff_lj);
            // Destroy the old structure. Must be done after gmx_pme_reinit in case currenSetup_==0.
            if (setup->pmedata != nullptr)
            {
                gmx_pme_destroy(setup->pmedata, false);
            }
            setup->pmedata = newPmeData;
        }
    }
    else
    {
        /* Tell our PME-only rank to switch grid */
        gmx_pme_send_switchgrid(*dd, setup->grid, setup->ewaldcoeff_q, setup->ewaldcoeff_lj);
    }
}

/*! \brief Checks the cycles and adds them to \p *set
 *
 * \returns whether we should increase the number of stages by 1.
 */
static bool processCycles(FILE*           fp_err,
                          const MDLogger& mdlog,
                          const double    cycles,
                          const int64_t   step,
                          const bool      isInLastStage,
                          pme_setup_t*    set)
{
    bool increaseNumStages = false;

    const auto buf = gmx::formatString("step %4" PRId64 ": ", step);
    printGrid(fp_err, mdlog, buf.c_str(), "timed with", *set, cycles);

    GMX_RELEASE_ASSERT(set->count > c_numPostSwitchTuningIntervalSkip, "We should skip cycles");
    if (set->count == (c_numPostSwitchTuningIntervalSkip + 1))
    {
        set->cycles = cycles;
    }
    else
    {
        if (cycles * maxFluctuationAccepted < set->cycles && isInLastStage)
        {
            /* The performance went up a lot (due to e.g. DD load balancing).
             * Add a stage, keep the minima, but rescan all setups.
             */
            increaseNumStages = true;

            if (debug)
            {
                fprintf(debug,
                        "The performance for grid %d %d %d went from %.3f to %.1f M-cycles, this "
                        "is more than %f\n"
                        "Will increase the number stages to by 1"
                        " and ignoring the previous performance\n",
                        set->grid[XX],
                        set->grid[YY],
                        set->grid[ZZ],
                        set->cycles * 1e-6,
                        cycles * 1e-6,
                        maxFluctuationAccepted);
            }
        }
        set->cycles = std::min(set->cycles, cycles);
    }

    return increaseNumStages;
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
void PmeLoadBalancing::Impl::balance(FILE*                fp_err,
                                     const matrix         box,
                                     ArrayRef<const RVec> x,
                                     double               cycles,
                                     interaction_const_t* ic,
                                     nonbonded_verlet_t*  nbv,
                                     gmx_pme_t**          pmedata,
                                     int64_t              step)
{
    if (dd_.nnodes > 1)
    {
        // Average the cycles over all PP ranks
        dd_.mpiComm().sumReduce(1, &cycles);
        cycles /= dd_.mpiComm().size();
    }

    setups_[currentSetup_].count++;

    /* Skip the first c_numPostSwitchTuningIntervalSkip cycles because the first step
     * after a switch is much slower due to allocation and/or caching effects.
     */
    if (setups_[currentSetup_].count % (c_numPostSwitchTuningIntervalSkip + 1) != 0)
    {
        return;
    }

    const bool increaseNumStages = processCycles(
            fp_err, mdlog_, cycles, step, stage_ == numStages_ - 1, &setups_[currentSetup_]);

    if (increaseNumStages)
    {
        numStages_++;
    }

    if (setups_[currentSetup_].cycles < setups_[fastestSetup_].cycles)
    {
        fastestSetup_ = currentSetup_;

        if (haveDDAtomOrdering(&dd_))
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
            set_dd_dlb_max_cutoff(&dd_, setups_[fastestSetup_].rlistOuter);
        }
    }
    const double cyclesFastest = setups_[fastestSetup_].cycles;

    /* Check in stage 0 if we should stop scanning grids.
     * Stop when the time is more than maxRelativeSlowDownAccepted longer than the fastest.
     */
    if (stage_ == 0 && currentSetup_ > 0 && cycles > setups_[fastestSetup_].cycles * maxRelativeSlowdownAccepted)
    {
        setups_.resize(currentSetup_ + 1);
        /* Done with scanning, go to stage 1 */
        switchToStage1();
    }

    if (stage_ == 0)
    {
        const int gridsizeStart = numPmeGridPoints(setups_[currentSetup_]);

        bool continueToNextSetup;
        do
        {
            bool haveNextSetup;

            if (currentSetup_ + 1 < gmx::ssize(setups_))
            {
                /* We had already generated the next setup */
                haveNextSetup = true;
            }
            else
            {
                /* Find the next setup, when possible */
                haveNextSetup = increaseCutoff();

                if (!haveNextSetup)
                {
                    limited_ = PmeLoadBalancingLimit::PmeGrid;
                }
            }

            if (haveNextSetup
                && setups_[currentSetup_ + 1].spacing > c_maxSpacingScaling * setups_[0].spacing)
            {
                haveNextSetup = false;
                limited_      = PmeLoadBalancingLimit::MaxScaling;
            }

            if (haveNextSetup && ir_.pbcType != PbcType::No)
            {
                haveNextSetup =
                        (square(setups_[currentSetup_ + 1].rlistOuter) <= max_cutoff2(ir_.pbcType, box));
                if (!haveNextSetup)
                {
                    limited_ = PmeLoadBalancingLimit::Box;
                }
            }

            if (haveNextSetup)
            {
                currentSetup_++;

                if (haveDDAtomOrdering(&dd_))
                {
                    const bool checkGpuDdLimitation = true;
                    haveNextSetup                   = change_dd_cutoff(
                            &dd_, box, x, setups_[currentSetup_].rlistOuter, checkGpuDdLimitation);
                    if (!haveNextSetup)
                    {
                        /* Failed: do not use this setup */
                        currentSetup_--;
                        limited_ = PmeLoadBalancingLimit::DD;
                    }
                }
            }
            if (!haveNextSetup)
            {
                /* We hit the upper limit for the cut-off,
                 * the setup should not go further than currentSetup_.
                 */
                setups_.resize(currentSetup_ + 1);
                printLoadBalLimited(fp_err, mdlog_, step, limited_, setups_[setupEnd() - 1]);
                /* Switch to the next stage */
                switchToStage1();
            }

            continueToNextSetup = haveNextSetup;
            if (continueToNextSetup)
            {
                const auto cs = setups_[currentSetup_];

                // Skip this setup when it is not sufficiently coarser or more efficient
                continueToNextSetup =
                        !(cs.grid[XX] * cs.grid[YY] * cs.grid[ZZ] < gridsizeStart * gridpointsScaleFactor
                          && cs.grid_efficiency < setups_[currentSetup_ - 1].grid_efficiency
                                                          * relativeEfficiencyFactor);
            }
        } while (continueToNextSetup);
    }

    if (stage_ > 0 && endSetup_ == 1)
    {
        currentSetup_ = lowerLimit_;
        stage_        = numStages_;
    }
    else if (stage_ > 0 && endSetup_ > 1)
    {
        /* If stage_ = numStages_-1:
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
            if (currentSetup_ > startSetup_)
            {
                currentSetup_--;
            }
            else
            {
                stage_++;

                currentSetup_ = endSetup_ - 1;
            }
        } while (stage_ == numStages_ - 1 && setups_[currentSetup_].count > 0
                 && setups_[currentSetup_].cycles > cyclesFastest * maxRelativeSlowdownAccepted);

        if (stage_ == numStages_)
        {
            /* We are done optimizing, use the fastest setup we found */
            currentSetup_ = fastestSetup_;
        }
    }

    if (haveDDAtomOrdering(&dd_) && stage_ > 0)
    {
        const bool checkGpuDdLimitation = true;

        const bool cutoffIsAllowed =
                change_dd_cutoff(&dd_, box, x, setups_[currentSetup_].rlistOuter, checkGpuDdLimitation);
        if (!cutoffIsAllowed)
        {
            /* For some reason the chosen cut-off is incompatible with DD.
             * We should continue scanning a more limited range of cut-off's.
             */
            if (currentSetup_ > 1 && stage_ == numStages_)
            {
                /* stage_=numStages_ says we're finished, but we should continue
                 * balancing, so we set back stage_ which was just incremented.
                 */
                stage_--;
            }
            if (currentSetup_ <= fastestSetup_)
            {
                /* This should not happen, as we set limits on the DLB bounds.
                 * But we implement a complete failsafe solution anyhow.
                 */
                GMX_LOG(mdlog_.warning)
                        .asParagraph()
                        .appendTextFormatted(
                                "The fastest PP/PME load balancing setting (cutoff %.3f nm) is no "
                                "longer available due to DD DLB or box size limitations",
                                setups_[fastestSetup_].rcut_coulomb);
                fastestSetup_ = lowerLimit_;
                startSetup_   = lowerLimit_;
            }
            /* Limit the range to below the current cut-off, scan from start */
            endSetup_     = currentSetup_;
            currentSetup_ = startSetup_;
            limited_      = PmeLoadBalancingLimit::DD;
            printLoadBalLimited(fp_err, mdlog_, step, limited_, setups_[setupEnd() - 1]);
        }
    }

    pme_setup_t& setup = setups_[currentSetup_];

    /* Change the Coulomb cut-off and the PME grid */
    applySetup(&setup, setups_[0].pmedata, ir_, ic, nbv, &dd_);

    if (!haveSepPMERanks_)
    {
        *pmedata = setup.pmedata;
    }

    if (stage_ == numStages_)
    {
        printGrid(fp_err, mdlog_, "", "optimal", setup, -1);
    }
}

void PmeLoadBalancing::Impl::addTwoStages(const bool dlbWasUnlocked)
{
    /* Add 2 tuning stages, keep the detected end of the setup range */
    numStages_ += 2;
    if (dlbWasUnlocked && haveSepPMERanks_)
    {
        /* With separate PME ranks, DLB should always lower the PP load and
         * can only increase the PME load (more communication and imbalance),
         * so we only need to scan longer cut-off's.
         */
        lowerLimit_ = currentSetup_;
    }
    startSetup_ = lowerLimit_;
}

void PmeLoadBalancing::Impl::addCycles(FILE*                fp_err,
                                       t_forcerec*          fr,
                                       const matrix         box,
                                       ArrayRef<const RVec> x,
                                       const gmx_wallcycle* wcycle,
                                       int64_t              step,
                                       int64_t              step_rel)
{
    if (!isActive_)
    {
        return;
    }

    const int    cyclesNumberPrev  = cyclesNumber_;
    const double cyclesCounterPrev = cyclesCounter_;
    wallcycle_get(wcycle, WallCycleCounter::Step, &cyclesNumber_, &cyclesCounter_);

    /* Before the first step we haven't done any steps yet.
     * Also handle cases where ir.init_step % nstlist != 0.
     * We also want to skip a number of steps and seconds while
     * the CPU and GPU, when used, performance stabilizes.
     */
    if (haveDDAtomOrdering(&dd_) && DDMAIN(&dd_))
    {
        startupTimeDelayElapsed_ = (gmx_gettime() - startTime_ < c_startupTimeDelay);
    }
    if (haveDDAtomOrdering(&dd_))
    {
        dd_bcast(&dd_, sizeof(bool), &startupTimeDelayElapsed_);
    }

    if (cyclesNumber_ == 0 || step_rel < c_numFirstTuningIntervalSkip * ir_.nstlist || startupTimeDelayElapsed_)
    {
        return;
    }
    /* Sanity check, we expect nstlist cycle counts */
    if (cyclesNumber_ - cyclesNumberPrev != ir_.nstlist)
    {
        /* We could return here, but it's safer to issue an error and quit */
        gmx_incons("pme_loadbal_do called at an interval != nstlist");
    }

    /* PME grid + cut-off optimization with GPUs or PME ranks */
    if (!isInBalancingPhase_ && haveSepPMERanks_)
    {
        if (triggerOnDLB_)
        {
            isInBalancingPhase_ = dd_dlb_is_on(&dd_);
        }
        /* We should ignore the first timing to avoid timing allocation
         * overhead. And since the PME load balancing is called just
         * before DD repartitioning, the ratio returned by dd_pme_f_ratio
         * is not over the last nstlist steps, but the nstlist steps before
         * that. So the first useful ratio is available at step_rel=3*nstlist.
         */
        else if (step_rel >= c_numFirstTuningIntervalSkipWithSepPme * ir_.nstlist)
        {
            GMX_ASSERT(haveDDAtomOrdering(&dd_), "Domain decomposition should be active here");
            if (DDMAIN(&dd_))
            {
                /* If PME rank load is too high, start tuning. If
                   PME-PP direct GPU communication is active,
                   unconditionally start tuning since ratio will be
                   unreliable due to CPU-GPU asynchronicity in codepath */
                isInBalancingPhase_ = useGpuPmePpCommunication_
                                              ? true
                                              : (dd_pme_f_ratio(&dd_) >= loadBalanceTriggerFactor);
            }
            dd_bcast(&dd_, sizeof(bool), &isInBalancingPhase_);
        }

        isActive_ = (isInBalancingPhase_ || step_rel <= stepRelStop_);
    }

    /* The location in the code of this balancing termination is strange.
     * You would expect to have it after the call to pme_load_balance()
     * below, since there stage_ is updated.
     * But when terminating directly after deciding on and selecting the
     * optimal setup, DLB will turn on right away if it was locked before.
     * This might be due to PME reinitialization. So we check stage here
     * to allow for another nstlist steps with DLB locked to stabilize
     * the performance.
     */
    if (isInBalancingPhase_ && stage_ == numStages_)
    {
        isInBalancingPhase_ = false;

        if (haveDDAtomOrdering(&dd_) && dd_dlb_is_locked(&dd_))
        {
            /* Unlock the DLB=auto, DLB is allowed to activate */
            dd_dlb_unlock(&dd_);
            GMX_LOG(mdlog_.warning)
                    .asParagraph()
                    .appendText("NOTE: DLB can now turn on, when beneficial");

            /* We don't deactivate the tuning yet, since we will balance again
             * after DLB gets turned on, if it does within PMETune_period.
             */
            addTwoStages(true);
            triggerOnDLB_ = true;
            stepRelStop_  = step_rel + PMETunePeriod * ir_.nstlist;
        }
        else
        {
            /* We're completely done with PME tuning */
            isActive_ = false;
        }

        if (haveDDAtomOrdering(&dd_))
        {
            /* Set the cut-off limit to the final selected cut-off,
             * so we don't have artificial DLB limits.
             * This also ensures that we won't disable the currently
             * optimal setting during a second round of PME balancing.
             */
            set_dd_dlb_max_cutoff(&dd_, fr->nbv->pairlistOuterRadius());
        }
    }

    if (isInBalancingPhase_)
    {
        /* We might not have collected nstlist steps in cycles yet,
         * since init_step might not be a multiple of nstlist,
         * but the first data collected is skipped anyhow.
         */
        balance(fp_err, box, x, cyclesCounter_ - cyclesCounterPrev, fr->ic.get(), fr->nbv.get(), &fr->pmedata, step);

        /* Update deprecated rlist in forcerec to stay in sync with fr->nbv */
        fr->rlist = fr->nbv->pairlistOuterRadius();

        if (fr->dispersionCorrection)
        {
            fr->dispersionCorrection->setParameters(*fr->ic);
        }
    }

    if (!isInBalancingPhase_ && (!haveSepPMERanks_ || step_rel > stepRelStop_))
    {
        /* We have just deactivated the balancing and we're not measuring PP/PME
         * imbalance during the first steps of the run: deactivate the tuning.
         */
        isActive_ = false;
    }

    if (!isActive_ && haveDDAtomOrdering(&dd_) && dd_dlb_is_locked(&dd_))
    {
        /* Make sure DLB is allowed when we deactivate PME tuning */
        dd_dlb_unlock(&dd_);
        GMX_LOG(mdlog_.warning)
                .asParagraph()
                .appendText("NOTE: DLB can now turn on, when beneficial");
    }
}

/*! \brief Print one load-balancing setting */
static void printLoadBalSetup(const MDLogger& mdlog, const char* name, const pme_setup_t& setup)
{
    GMX_LOG(mdlog.info)
            .appendTextFormatted("   %-7s %6.3f nm %6.3f nm     %3d %3d %3d   %5.3f nm  %5.3f nm",
                                 name,
                                 setup.rcut_coulomb,
                                 setup.rlistInner,
                                 setup.grid[XX],
                                 setup.grid[YY],
                                 setup.grid[ZZ],
                                 setup.spacing,
                                 1 / setup.ewaldcoeff_q);
}

/*! \brief Print all load-balancing settings */
static void printLoadBalSettings(const PmeLoadBalancingLimit limited,
                                 const bool                  currentSetupIsLastSetup,
                                 const pme_setup_t&          currentSetup,
                                 const pme_setup_t&          originalSetup,
                                 const bool                  useGpuForNonbondeds,
                                 const MDLogger&             mdlog)
{
    double pp_ratio, grid_ratio;
    real   pp_ratio_temporary;

    pp_ratio_temporary = currentSetup.rlistInner / originalSetup.rlistInner;
    pp_ratio           = power3(pp_ratio_temporary);
    grid_ratio = numPmeGridPoints(currentSetup) / static_cast<double>(numPmeGridPoints(originalSetup));

    GMX_LOG(mdlog.info)
            .appendText(
                    "\n"
                    "       P P   -   P M E   L O A D   B A L A N C I N G\n"
                    "\n");
    /* Here we only warn when the optimal setting is the last one */
    if (limited != PmeLoadBalancingLimit::No && currentSetupIsLastSetup)
    {
        GMX_LOG(mdlog.info)
                .appendTextFormatted(
                        " NOTE: The PP/PME load balancing was limited by the %s,\n"
                        "       you might not have reached a good load balance.",
                        enumValueToString(limited));
        if (limited == PmeLoadBalancingLimit::DD)
        {
            GMX_LOG(mdlog.info)
                    .appendText("       Try different mdrun -dd settings or lower the -dds value.");
        }
        // Add empty line
        GMX_LOG(mdlog.info).appendText("");
    }
    GMX_LOG(mdlog.info)
            .appendText(
                    " PP/PME load balancing changed the cut-off and PME settings:\n"
                    "           particle-particle                    PME\n"

                    "            rcoulomb  rlist            grid      spacing   1/beta");
    printLoadBalSetup(mdlog, "initial", originalSetup);
    printLoadBalSetup(mdlog, "final", currentSetup);
    GMX_LOG(mdlog.info).appendTextFormatted(" cost-ratio           %4.2f             %4.2f", pp_ratio, grid_ratio);
    GMX_LOG(mdlog.info)
            .appendText(
                    " (note that these numbers concern only part of the total PP and PME load)");

    if (pp_ratio > 1.5 && !useGpuForNonbondeds)
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
        // Add empty line
        GMX_LOG(mdlog.info).appendText("");
    }
}

void PmeLoadBalancing::Impl::printSettings() const
{
    printLoadBalSettings(limited_,
                         currentSetup_ == setupEnd() - 1,
                         setups_[currentSetup_],
                         setups_[0],
                         useGpuForNonbondeds_,
                         mdlog_);
}

PmeLoadBalancing::PmeLoadBalancing(gmx_domdec_t*              dd,
                                   const MDLogger&            mdlog,
                                   const t_inputrec&          ir,
                                   const matrix               box,
                                   const interaction_const_t& ic,
                                   const nonbonded_verlet_t&  nbv,
                                   gmx_pme_t*                 pmedata,
                                   const SimulationWorkload&  simulationWork) :
    impl_(std::make_unique<Impl>(dd, mdlog, ir, box, ic, nbv, pmedata, simulationWork))
{
}

PmeLoadBalancing::~PmeLoadBalancing() = default;

bool PmeLoadBalancing::isActive() const
{
    return impl_->isActive();
}

bool PmeLoadBalancing::isPrintingLoad() const
{
    return impl_->isPrintingLoad();
}

void PmeLoadBalancing::addCycles(FILE*                fp_err,
                                 t_forcerec*          fr,
                                 const matrix         box,
                                 ArrayRef<const RVec> x,
                                 const gmx_wallcycle* wcycle,
                                 int64_t              step,
                                 int64_t              step_rel)
{
    impl_->addCycles(fp_err, fr, box, x, wcycle, step, step_rel);
}

void PmeLoadBalancing::printSettings() const
{
    impl_->printSettings();
}

} // namespace gmx
