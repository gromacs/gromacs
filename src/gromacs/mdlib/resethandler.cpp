/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * \brief
 * Defines the reset handler class.
 *
 * \author Pascal Merz <pascal.merz@colorado.edu>
 * \ingroup module_mdlib
 */

#include "gmxpre.h"

#include "resethandler.h"

#include <cinttypes>

#include <filesystem>

#include "gromacs/domdec/domdec.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme_load_balancing.h"
#include "gromacs/ewald/pme_pp.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/mdrunutility/printtime.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/nbnxm/gpu_data_mgmt.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

/*! \brief Convert signed char (as used by SimulationSignal) to ResetSignal enum
 *
 * Expected values are
 *   \p sig == 0 -- no signal
 *   \p sig >= 1 -- signal received
 */
static inline ResetSignal convertToResetSignal(signed char sig)
{
    GMX_ASSERT(sig >= 0, "Unexpected reset signal < 0 received");
    return sig >= 1 ? ResetSignal::doResetCounters : ResetSignal::noSignal;
}

ResetHandler::ResetHandler(compat::not_null<SimulationSignal*> signal,
                           bool                                simulationsShareState,
                           int64_t                             nsteps,
                           bool                                isMain,
                           bool                                resetHalfway,
                           real                                maximumHoursToRun,
                           const MDLogger&                     mdlog,
                           gmx_wallcycle*                      wcycle,
                           gmx_walltime_accounting_t           walltime_accounting) :
    signal_(*signal), rankCanSetSignal_(false), simulationNeedsReset_(false), maximumHoursToRun_(maximumHoursToRun)
{
    if (simulationsShareState)
    {
        signal_.isLocal = false;
    }
    if (resetHalfway)
    {
        GMX_LOG(mdlog.info)
                .asParagraph()
                .appendText(
                        "The -resethway functionality is deprecated, and may be removed in a "
                        "future version.");
        if (nsteps > 0)
        {
            /* Signal to reset the counters half the simulation steps. */
            wcycle_set_reset_counters(wcycle, nsteps / 2);
        }
        simulationNeedsReset_ = true;

        if (isMain && (maximumHoursToRun > 0))
        {
            rankCanSetSignal_ = true;
        }
    }
    else if (wcycle_get_reset_counters(wcycle) > 0)
    {
        simulationNeedsReset_ = true;
    }
    else
    {
        // if no reset is happening, this will always be valid
        walltime_accounting_set_valid_finish(walltime_accounting);
    }
}

bool ResetHandler::setSignalImpl(gmx_walltime_accounting_t walltime_accounting)
{
    const double secondsSinceStart = walltime_accounting_get_time_since_start(walltime_accounting);
    if (secondsSinceStart > maximumHoursToRun_ * 60.0 * 60.0 * 0.495)
    {
        /* Set flag that will communicate the signal to all ranks in the simulation */
        signal_.sig = static_cast<signed char>(ResetSignal::doResetCounters);
        /* Let handler know that we did signal a reset */
        return true;
    }
    /* Let handler know that we did not signal a reset */
    return false;
}

bool ResetHandler::resetCountersImpl(int64_t                     step,
                                     int64_t                     step_rel,
                                     const MDLogger&             mdlog,
                                     FILE*                       fplog,
                                     const t_commrec*            cr,
                                     nonbonded_verlet_t*         nbv,
                                     t_nrnb*                     nrnb,
                                     const gmx_pme_t*            pme,
                                     const pme_load_balancing_t* pme_loadbal,
                                     gmx_wallcycle*              wcycle,
                                     gmx_walltime_accounting_t   walltime_accounting)
{
    /* Reset either if signal has been passed, or if reset step has been reached */
    if (convertToResetSignal(signal_.set) == ResetSignal::doResetCounters
        || step_rel == wcycle_get_reset_counters(wcycle))
    {
        if (pme_loadbal_is_active(pme_loadbal))
        {
            /* Do not permit counter reset while PME load
             * balancing is active. The only purpose for resetting
             * counters is to measure reliable performance data,
             * and that can't be done before balancing
             * completes.
             *
             * TODO consider fixing this by delaying the reset
             * until after load balancing completes,
             * e.g. https://gerrit.gromacs.org/#/c/4964/2 */
            gmx_fatal(FARGS,
                      "PME tuning was still active when attempting to "
                      "reset mdrun counters at step %" PRId64
                      ". Try "
                      "resetting counters later in the run, e.g. with gmx "
                      "mdrun -resetstep.",
                      step);
        }

        char sbuf[STEPSTRSIZE];

        /* Reset all the counters related to performance over the run */
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendTextFormatted("step %s: resetting all time and cycle counters",
                                     gmx_step_str(step, sbuf));

        if (nbv && nbv->useGpu())
        {
            gpu_reset_timings(nbv);
        }

        if (pme_gpu_task_enabled(pme))
        {
            pme_gpu_reset_timings(pme);
        }

        if ((nbv && nbv->useGpu()) || pme_gpu_task_enabled(pme))
        {
            resetGpuProfiler();
        }

        wallcycle_stop(wcycle, WallCycleCounter::Run);
        wallcycle_reset_all(wcycle);
        if (haveDDAtomOrdering(*cr))
        {
            reset_dd_statistics_counters(cr->dd);
        }
        clear_nrnb(nrnb);
        wallcycle_start(wcycle, WallCycleCounter::Run);
        walltime_accounting_reset_time(walltime_accounting, step);
        print_date_and_time(fplog, cr->nodeid, "Restarted time", gmx_gettime());

        wcycle_set_reset_counters(wcycle, -1);
        if (!thisRankHasDuty(cr, DUTY_PME))
        {
            /* Tell our PME node to reset its counters */
            gmx_pme_send_resetcounters(cr, step);
        }
        /* Reset can only happen once, so clear the triggering flag. */
        signal_.set = static_cast<signed char>(ResetSignal::noSignal);
        /* We have done a reset, so the finish will be valid. */
        walltime_accounting_set_valid_finish(walltime_accounting);
        /* Let handler know that we handled a reset */
        return true;
    }

    /* Let handler know that we did not handle a reset */
    return false;
}

} // namespace gmx
