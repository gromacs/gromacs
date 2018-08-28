/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Defines the reset handler class.
 *
 * \author Pascal Merz <pascal.merz@colorado.edu>
 * \ingroup module_mdlib
 */

#include "gmxpre.h"

#include "resethandler.h"

#include "gromacs/domdec/domdec.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/mdlib/nbnxn_gpu_data_mgmt.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"

using namespace gmx;

ResetHandler::ResetHandler(
        gmx::SimulationSignal      *sig,
        bool                        needSync,
        const t_inputrec           *ir,
        const t_commrec            *cr,
        const MdrunOptions         &mdrunOptions,
        const gmx::MDLogger        &mdlog,
        FILE                       *fplog,
        const int64_t              &step,
        const int64_t              &step_rel,
        nonbonded_verlet_t         *nbv,
        t_nrnb                     *nrnb,
        const gmx_pme_t            *pme,
        const pme_load_balancing_t *pme_loadbal,
        gmx_walltime_accounting_t   walltime_accounting,
        gmx_wallcycle_t             wcycle) :
    maximumHoursToRun(mdrunOptions.maximumHoursToRun),
    signal(sig),
    step(step),
    step_rel(step_rel),
    mdlog(mdlog),
    fplog(fplog),
    cr(cr),
    nbv(nbv),
    nrnb(nrnb),
    pme(pme),
    pme_loadbal(pme_loadbal),
    wcycle(wcycle),
    walltime_accounting(walltime_accounting)

{
    if (needSync)
    {
        signal->isLocal = false;
    }
    if (mdrunOptions.timingOptions.resetHalfway)
    {
        GMX_LOG(mdlog.info).asParagraph().
            appendText(
                "The -resethway functionality is deprecated, and may be removed in a future version.");
        if (ir->nsteps > 0)
        {
            /* Signal to reset the counters half the simulation steps. */
            wcycle_set_reset_counters(wcycle, ir->nsteps / 2);
        }
        doHandle = true;

        if (MASTER(cr) && (mdrunOptions.maximumHoursToRun > 0))
        {
            doSet = true;
        }
    }
    else if (wcycle_get_reset_counters(wcycle) > 0)
    {
        doHandle = true;
    }
    else
    {
        // if no reset is happening, this will always be valid
        walltime_accounting_set_valid_finish(walltime_accounting);
    }
}

bool ResetHandler::setSignalImpl()
{
    const double secondsSinceStart = walltime_accounting_get_time_since_start(walltime_accounting);
    if (secondsSinceStart > maximumHoursToRun * 60.0 * 60.0 * 0.495)
    {
        /* Set flag that will communicate the signal to all ranks in the simulation */
        signal->sig = 1;
        /* Let helper know that we did signal a reset */
        return true;
    }
    /* Let helper know that we did not signal a reset */
    return false;
}

bool ResetHandler::handleSignalImpl()
{
    /* Reset either if signal has been passed,  */
    if (signal->set != 0 || step_rel == wcycle_get_reset_counters(wcycle))
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
            gmx_fatal(FARGS, "PME tuning was still active when attempting to "
                      "reset mdrun counters at step %" PRId64 ". Try "
                      "resetting counters later in the run, e.g. with gmx "
                      "mdrun -resetstep.", step);
        }

        char sbuf[STEPSTRSIZE];

        /* Reset all the counters related to performance over the run */
        GMX_LOG(mdlog.warning).asParagraph().appendTextFormatted(
                "step %s: resetting all time and cycle counters",
                gmx_step_str(step, sbuf));

        if (use_GPU(nbv))
        {
            nbnxn_gpu_reset_timings(nbv);
        }

        if (pme_gpu_task_enabled(pme))
        {
            pme_gpu_reset_timings(pme);
        }

        if (use_GPU(nbv) || pme_gpu_task_enabled(pme))
        {
            resetGpuProfiler();
        }

        wallcycle_stop(wcycle, ewcRUN);
        wallcycle_reset_all(wcycle);
        if (DOMAINDECOMP(cr))
        {
            reset_dd_statistics_counters(cr->dd);
        }
        init_nrnb(nrnb);
        wallcycle_start(wcycle, ewcRUN);
        walltime_accounting_reset_time(walltime_accounting, step);
        print_date_and_time(fplog, cr->nodeid, "Restarted time", gmx_gettime());

        wcycle_set_reset_counters(wcycle, -1);
        if (!thisRankHasDuty(cr, DUTY_PME))
        {
            /* Tell our PME node to reset its counters */
            gmx_pme_send_resetcounters(cr, step);
        }
        /* Reset can only happen once, so clear the triggering flag. */
        signal->set = 0;
        /* We have done a reset, so the finish will be valid. */
        walltime_accounting_set_valid_finish(walltime_accounting);
        /* Let helper know that we handled a reset */
        return true;
    }

    /* Let helper know that we did not handle a reset */
    return false;
}
