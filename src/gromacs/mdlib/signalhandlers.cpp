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

#include <gromacs/utility/fatalerror.h>
#include <gromacs/ewald/pme-load-balancing.h>
#include <gromacs/ewald/pme.h>
#include <gromacs/gmxlib/nrnb.h>
#include <gromacs/domdec/domdec.h>
#include <gromacs/gpu_utils/gpu_utils.h>
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/mdtypes/inputrec.h"
#include "signalhandlers.h"
#include "sim_util.h"
#include "nbnxn_gpu_data_mgmt.h"

gmx::CheckpointSignalHandler::CheckpointSignalHandler(gmx::SimulationSignal   * sig,
                                                      const t_inputrec* const   ir,
                                                      const t_commrec* const    cr,
                                                      const MdrunOptions       &mdrunOptions,
                                                      const gmx_bool           &bGStat,
                                                      const gmx_bool           &bNS,
                                                      const gmx_bool           &bLastStep,
                                                      const gmx_int64_t        &step,
                                                      gmx_bool                 &bCPT,
                                                      gmx_walltime_accounting_t walltime_accounting) :
    SignalHandler(sig),
    nstlist(ir->nstlist),
    init_step(ir->init_step),
    cpt_period(mdrunOptions.checkpointOptions.period),
    isMaster(MASTER(cr)),
    isParallel(PAR(cr)),
    writeConfout(bool(mdrunOptions.writeConfout)),
    isRerun(bool(mdrunOptions.rerun)),
    bGStat(bGStat),
    bNS(bNS),
    bLastStep(bLastStep),
    step(step),
    bCPT(bCPT),
    walltime_accounting(walltime_accounting) {}

void gmx::CheckpointSignalHandler::set_signal()
{
    /* In parallel we only have to check for checkpointing in steps
     * where we do global communication,
     *  otherwise the other nodes don't know.
     */
    const double secondsSinceStart = walltime_accounting_get_time_since_start(walltime_accounting);
    if (isMaster && ((bGStat || !isParallel) &&
                     cpt_period >= 0 &&
                     (cpt_period == 0 ||
                      secondsSinceStart >= nchkpt * cpt_period * 60.0)) &&
        signal->set == 0)
    {
        signal->sig = 1;
        nchkpt++;
    }
}

void gmx::CheckpointSignalHandler::handle_signal()
{
    /* We write a checkpoint at this MD step when:
     * either at an NS step when we signalled through gs,
     * or at the last step (but not when we do not want confout),
     * but never at the first step or with rerun.
     */
    bCPT = (((signal->set && (bNS || nstlist == 0)) ||
             (bLastStep && writeConfout)) &&
            (step > init_step && !isRerun));
    if (bCPT)
    {
        signal->set = 0;
    }
}

gmx::StopSignalHandler::StopSignalHandler(gmx::SimulationSignal   * sig,
                                          const t_inputrec* const   ir,
                                          const t_commrec* const    cr,
                                          const MdrunOptions       &mdrunOptions,
                                          int                       nstSignalComm,
                                          FILE                    * fplog,
                                          const gmx_bool           &bNS,
                                          const gmx_int64_t        &step,
                                          gmx_bool                 &bLastStep,
                                          gmx_walltime_accounting_t walltime_accounting) :
    SignalHandler(sig),
    maximumHoursToRun(mdrunOptions.maximumHoursToRun),
    isMaster(MASTER(cr)),
    reproducible(bool(mdrunOptions.reproducible)),
    nstSignalComm(nstSignalComm),
    nstlist(ir->nstlist),
    bNS(bNS),
    step(step),
    bLastStep(bLastStep),
    fplog(fplog),
    walltime_accounting(walltime_accounting)
{}

void gmx::StopSignalHandler::set_signal()
{
    const double secondsSinceStart = walltime_accounting_get_time_since_start(walltime_accounting);

    /* Check whether everything is still allright */
    if (((int) gmx_get_stop_condition() > handled_stop_condition)
#if GMX_THREAD_MPI
        && isMaster
#endif
        )
    {
        int nsteps_stop = -1;

        /* this just makes signals[].sig compatible with the hack
           of sending signals around by MPI_Reduce together with
           other floats */
        if ((gmx_get_stop_condition() == gmx_stop_cond_next_ns) ||
            (reproducible && gmx_get_stop_condition() == gmx_stop_cond_next))
        {
            /* We need at least two global communication steps to pass
             * around the signal. We stop at a pair-list creation step
             * to allow for exact continuation, when possible.
             */
            signal->sig = 1;
            nsteps_stop = std::max(nstlist, 2 * nstSignalComm);
        }
        else if (gmx_get_stop_condition() == gmx_stop_cond_next)
        {
            /* Stop directly after the next global communication step.
             * This breaks exact continuation.
             */
            signal->sig = -1;
            nsteps_stop = nstSignalComm + 1;
        }
        if (fplog)
        {
            fprintf(fplog,
                    "\n\nReceived the %s signal, stopping within %d steps\n\n",
                    gmx_get_signal_name(), nsteps_stop);
            fflush(fplog);
        }
        fprintf(stderr,
                "\n\nReceived the %s signal, stopping within %d steps\n\n",
                gmx_get_signal_name(), nsteps_stop);
        fflush(stderr);
        handled_stop_condition = (int) gmx_get_stop_condition();
    }
    else if (isMaster && (bNS || nstlist <= 0) &&
             (maximumHoursToRun > 0 &&
              secondsSinceStart > maximumHoursToRun * 60.0 * 60.0 * 0.99) &&
             signal->sig == 0 && signal->set == 0)
    {
        /* Signal to terminate the run */
        signal->sig = 1;
        char sbuf[STEPSTRSIZE];
        if (fplog)
        {
            fprintf(fplog, "\nStep %s: Run time exceeded %.3f hours, will terminate the run\n",
                    gmx_step_str(step, sbuf), maximumHoursToRun * 0.99);
        }
        fprintf(stderr, "\nStep %s: Run time exceeded %.3f hours, will terminate the run\n",
                gmx_step_str(step, sbuf), maximumHoursToRun * 0.99);
    }
}

void gmx::StopSignalHandler::handle_signal()
{
    /* < 0 means stop at next step, > 0 means stop at next NS step */
    if (signal->set < 0 ||
        (signal->set > 0 && ( bNS || nstlist == 0)))
    {
        bLastStep = TRUE;
    }
}

gmx::ResetSignalHandler::ResetSignalHandler(gmx::SimulationSignal           * sig,
                                            const t_inputrec* const           ir,
                                            const t_commrec* const            cr,
                                            const MdrunOptions               &mdrunOptions,
                                            const gmx_int64_t                &step,
                                            const gmx_int64_t                &step_rel,
                                            gmx_walltime_accounting_t         walltime_accounting,
                                            gmx_wallcycle_t                   wcycle,
                                            const gmx::MDLogger              &mdlog,
                                            FILE* const                       fplog,
                                            nonbonded_verlet_t* const         nbv,
                                            const gmx_pme_t* const            pme,
                                            t_nrnb* const                     nrnb,
                                            const pme_load_balancing_t* const pme_loadbal) :
    SignalHandler(sig),
    maximumHoursToRun(mdrunOptions.maximumHoursToRun),
    isMaster(MASTER(cr)),
    cr(cr),
    step(step),
    step_rel(step_rel),
    walltime_accounting(walltime_accounting),
    wcycle(wcycle),
    mdlog(mdlog),
    fplog(fplog),
    nbv(nbv),
    pme(pme),
    nrnb(nrnb),
    pme_loadbal(pme_loadbal)
{
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
        /* Signal to reset the counters halfway the simulation time.
         * Handled by master only
         * */
        bResetCountersHalfMaxH = isMaster &&
            (mdrunOptions.maximumHoursToRun > 0);
    }
}

void gmx::ResetSignalHandler::set_signal()
{
    if (bResetCountersHalfMaxH)
    {
        const double secondsSinceStart = walltime_accounting_get_time_since_start(walltime_accounting);
        if (secondsSinceStart > maximumHoursToRun * 60.0 * 60.0 * 0.495)
        {
            /* Set flag that will communicate the signal to all ranks in the simulation */
            signal->sig = 1;
            /* If mdrun -maxh -resethway was active, it can only trigger once */
            bResetCountersHalfMaxH = false;
        }
    }
    if (step_rel == wcycle_get_reset_counters(wcycle))
    {
        // Counter reset based on step does not need communication
        signal->set = 1;
    }
}

void gmx::ResetSignalHandler::handle_signal()
{
    /* If it is time to reset counters, set a flag that remains
       true until counters actually get reset */
    if (signal->set != 0)
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
                      "reset mdrun counters at step %" GMX_PRId64 ". Try "
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
    }
}

bool gmx::ResetSignalHandler::valid_finish() const
{
    return (step_rel >= wcycle_get_reset_counters(wcycle) &&
            signal->sig == 0 && signal->set == 0 &&
            !bResetCountersHalfMaxH);
}
