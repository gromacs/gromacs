//
// Created by Pascal Merz on 4/21/18.
//

#include "gromacs/mdlib/trajectory_writing.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/commrec.h"
#include "writetrajectoryelement.h"

void gmx::WriteTrajectoryElement::run()
{
    /* Now we have the energies and forces corresponding to the
     * coordinates at time t. We must output all of this before
     * the update.
     */
    do_md_trajectory_writing(integrator.fplog, integrator.cr, integrator.nfile, integrator.fnm, dataManager.step, dataManager.step_rel, dataManager.t,
                             integrator.inputrec, dataManager.state, integrator.state_global, integrator.observablesHistory,
                             integrator.top_global, integrator.fr,
                             dataManager.outf, dataManager.mdebin, dataManager.ekind, dataManager.f,
                             &dataManager.nchkpt,
                             dataManager.bCPT, dataManager.bRerunMD, dataManager.bLastStep,
                             integrator.mdrunOptions.writeConfout,
                             dataManager.bSumEkinhOld);

    dataManager.elapsed_time = walltime_accounting_get_current_elapsed_time(integrator.walltime_accounting);

    // TODO: does this belong here?
    /* Check whether everything is still allright */
    if (((int)gmx_get_stop_condition() > dataManager.handled_stop_condition)
#if GMX_THREAD_MPI
        && MASTER(integrator.cr)
#endif
        )
    {
        int nsteps_stop = -1;

        /* this just makes signals[].sig compatible with the hack
           of sending signals around by MPI_Reduce together with
           other floats */
        if ((gmx_get_stop_condition() == gmx_stop_cond_next_ns) ||
            (integrator.mdrunOptions.reproducible &&
             gmx_get_stop_condition() == gmx_stop_cond_next))
        {
            /* We need at least two global communication steps to pass
             * around the signal. We stop at a pair-list creation step
             * to allow for exact continuation, when possible.
             */
            dataManager.signals[eglsSTOPCOND].sig = 1;
            nsteps_stop                           = std::max(integrator.inputrec->nstlist, 2*dataManager.nstSignalComm);
        }
        else if (gmx_get_stop_condition() == gmx_stop_cond_next)
        {
            /* Stop directly after the next global communication step.
             * This breaks exact continuation.
             */
            dataManager.signals[eglsSTOPCOND].sig = -1;
            nsteps_stop                           = dataManager.nstSignalComm + 1;
        }
        if (integrator.fplog)
        {
            fprintf(integrator.fplog,
                    "\n\nReceived the %s signal, stopping within %d steps\n\n",
                    gmx_get_signal_name(), nsteps_stop);
            fflush(integrator.fplog);
        }
        fprintf(stderr,
                "\n\nReceived the %s signal, stopping within %d steps\n\n",
                gmx_get_signal_name(), nsteps_stop);
        fflush(stderr);
        dataManager.handled_stop_condition = (int)gmx_get_stop_condition();
    }
    else if (MASTER(integrator.cr) && (dataManager.bNS || integrator.inputrec->nstlist <= 0) &&
             (dataManager.max_hours > 0 && dataManager.elapsed_time > dataManager.max_hours*60.0*60.0*0.99) &&
             dataManager.signals[eglsSTOPCOND].sig == 0 && dataManager.signals[eglsSTOPCOND].set == 0)
    {
        /* Signal to terminate the run */
        dataManager.signals[eglsSTOPCOND].sig = 1;
        if (integrator.fplog)
        {
            fprintf(integrator.fplog, "\nStep %s: Run time exceeded %.3f hours, will terminate the run\n", gmx_step_str(dataManager.step, dataManager.sbuf), dataManager.max_hours*0.99);
        }
        fprintf(stderr, "\nStep %s: Run time exceeded %.3f hours, will terminate the run\n", gmx_step_str(dataManager.step, dataManager.sbuf), dataManager.max_hours*0.99);
    }

    if (dataManager.bResetCountersHalfMaxH && MASTER(integrator.cr) &&
        dataManager.elapsed_time > dataManager.max_hours*60.0*60.0*0.495)
    {
        /* Set flag that will communicate the signal to all ranks in the simulation */
        dataManager.signals[eglsRESETCOUNTERS].sig = 1;
    }

    /* In parallel we only have to check for checkpointing in steps
     * where we do global communication,
     *  otherwise the other nodes don't know.
     */
    const real cpt_period = integrator.mdrunOptions.checkpointOptions.period;
    if (MASTER(integrator.cr) && ((dataManager.bGStat || !PAR(integrator.cr)) &&
                                  cpt_period >= 0 &&
                                  (cpt_period == 0 ||
                                   dataManager.elapsed_time >= dataManager.nchkpt*cpt_period*60.0)) &&
        dataManager.signals[eglsCHKPT].set == 0)
    {
        dataManager.signals[eglsCHKPT].sig = 1;
    }
}
