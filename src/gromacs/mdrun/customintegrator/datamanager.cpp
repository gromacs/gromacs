//
// integrator.created by Pascal Merz on 4/21/18.
//
#include "gmxpre.h"

#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include <cmath>

#include <algorithm>
#include <memory>

#include "thread_mpi/threads.h"

#include "gromacs/awh/awh.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/compat/make_unique.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_network.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/essentialdynamics/edsam.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme-load-balancing.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/imd/imd.h"
#include "gromacs/listed-forces/manage-threading.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/compute_io.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/deform.h"
#include "gromacs/mdlib/ebin.h"
#include "gromacs/mdlib/expanded.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/md_support.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/mdebin.h"
#include "gromacs/mdlib/mdoutf.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/mdsetup.h"
#include "gromacs/mdlib/membed.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_gpu_data_mgmt.h"
#include "gromacs/mdlib/ns.h"
#include "gromacs/mdlib/repl_ex.h"
#include "gromacs/mdlib/shellfc.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdlib/simulationsignal.h"
#include "gromacs/mdlib/tgroup.h"
#include "gromacs/mdlib/trajectory_writing.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdlib/vcm.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdtypes/awh-history.h"
#include "gromacs/mdtypes/awh-params.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/df_history.h"
#include "gromacs/mdtypes/energyhistory.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/observableshistory.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/swap/swapcoords.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#include "gromacs/mdrun/integrator.h"

#include "datamanager.h"

gmx::DataManager::DataManager(Integrator &integrator_ref) :
    integrator(integrator_ref),
    continuationOptions(integrator_ref.mdrunOptions.continuationOptions)
{

    t_inputrec       *ir   = integrator.inputrec;

    //
    // MISSING FEATURES
    //
    if (opt2bSet("-ei", integrator.nfile, integrator.fnm) ||
        integrator.observablesHistory->edsamHistory != nullptr)
    {
        gmx_fatal(FARGS, "Essential dynamics is not implemented for custom integrators.");
    }

    if (ir->eSwapCoords != eswapNO)
    {
        gmx_fatal(FARGS, "Ion swapping is not implemented for custom integrators.");
    }

    /* Check for polarizable models and flexible constraints */
    shellfc = init_shell_flexcon(integrator.fplog,
                                 integrator.top_global,
                                 n_flexible_constraints(integrator.constr),
                                 ir->nstcalcenergy,
                                 DOMAINDECOMP(integrator.cr));
    if (shellfc)
    {
        gmx_fatal(FARGS, "Shell particles are not implemented for custom integrators.");
    }

    if (ir->bIMD)
    {
        gmx_fatal(FARGS, "Interactive MD is not implemented for custom integrators.");
    }

    if (ir->bExpanded)
    {
        gmx_fatal(FARGS, "Expanded ensemble is not implemented for custom integrators.");
    }

    /* Initialize AWH and restore state from history in checkpoint if needed. */
    if (ir->bDoAwh)
    {
        gmx_fatal(FARGS, "AWH is not implemented for custom integrators.");
    }

    const bool useReplicaExchange = (integrator.replExParams.exchangeInterval > 0);
    if (useReplicaExchange)
    {
        gmx_fatal(FARGS, "Replica exchange is not implemented for custom integrators.");
    }
    if (integrator.vsite)
    {
        gmx_fatal(FARGS, "Virtual sites are not implemented for custom integrators.");
    }

    if (ir->efep != efepNO)
    {
        gmx_fatal(FARGS, "Free energy calculation is not implemented for custom integrators.");
    }

    if (ir->bSimTemp)
    {
        gmx_fatal(FARGS, "Simulated tempering is not implemented for custom integrators.");
    }
    //
    // END MISSING FEATURES
    //

    if (!integrator.mdrunOptions.writeConfout)
    {
        // This is on by default, and the main known use case for
        // turning it off is for convenience in benchmarking, which is
        // something that should not show up in the general user
        // interface.
        GMX_LOG(integrator.mdlog.info).asParagraph().
            appendText("The -noconfout functionality is deprecated, and may be removed in a future version.");
    }
    if (integrator.mdrunOptions.timingOptions.resetHalfway)
    {
        GMX_LOG(integrator.mdlog.info).asParagraph().
            appendText("The -resethway functionality is deprecated, and may be removed in a future version.");
        if (ir->nsteps > 0)
        {
            /* Signal to reset the counters half the simulation steps. */
            wcycle_set_reset_counters(integrator.wcycle, ir->nsteps/2);
        }
        /* Signal to reset the counters halfway the simulation time. */
        bResetCountersHalfMaxH = (integrator.mdrunOptions.maximumHoursToRun > 0);
    }

    nstglobalcomm = integrator.mdrunOptions.globalCommunicationInterval;

    nstglobalcomm   = check_nstglobalcomm(integrator.mdlog, nstglobalcomm, ir);
    bGStatEveryStep = (nstglobalcomm == 1);

    groups = &integrator.top_global->groups;

    /* Initial values */
    init_md(integrator.fplog, integrator.cr, integrator.outputProvider, ir, integrator.oenv, integrator.mdrunOptions,
            &t, &t0, integrator.state_global, lam0,
            integrator.nrnb, integrator.top_global, &upd,
            integrator.nfile, integrator.fnm, &outf, &mdebin,
            force_vir, shake_vir, mu_tot, &bSimAnn, &vcm, integrator.wcycle);

    clear_mat(total_vir);
    clear_mat(pres);
    /* Energy terms and groups */
    snew(enerd, 1);
    init_enerdata(integrator.top_global->groups.grps[egcENER].nr, ir->fepvals->n_lambda,
                  enerd);

    /* Kinetic energy data */
    snew(ekind, 1);
    init_ekindata(integrator.fplog, integrator.top_global, &(ir->opts), ekind);
    /* Copy the cos acceleration to the groups struct */
    ekind->cosacc.cos_accel = ir->cos_accel;

    gstat = global_stat_init(ir);

    if (inputrecDeform(ir))
    {
        tMPI_Thread_mutex_lock(&deform_init_box_mutex);
        set_deform_reference_box(upd,
                                 deform_init_init_step_tpx,
                                 deform_init_box_tpx);
        tMPI_Thread_mutex_unlock(&deform_init_box_mutex);
    }

    {
        double io = compute_io(ir, integrator.top_global->natoms, groups, mdebin->ebin->nener, 1);
        if ((io > 2000) && MASTER(integrator.cr))
        {
            fprintf(stderr,
                    "\nWARNING: This run will generate roughly %.0f Mb of data\n\n",
                    io);
        }
    }

    // Local state only becomes valid now.

    if (DOMAINDECOMP(integrator.cr))
    {
        top = dd_init_local_top(integrator.top_global);

        stateInstance = gmx::compat::make_unique<t_state>();
        state         = stateInstance.get();
        dd_init_local_state(integrator.cr->dd, integrator.state_global, state);
    }
    else
    {
        state_change_natoms(integrator.state_global, integrator.state_global->natoms);
        /* We need to allocate one element extra, since we might use
         * (unaligned) 4-wide SIMD loads to access rvec entries.
         */
        f.resize(gmx::paddedRVecVectorSize(integrator.state_global->natoms));
        /* Copy the pointer to the global state */
        state = integrator.state_global;

        snew(top, 1);
        mdAlgorithmsSetupAtomData(integrator.cr, ir, integrator.top_global, top, integrator.fr,
                                  &graph, integrator.mdAtoms, integrator.vsite, shellfc);

        update_realloc(upd, state->natoms);
    }

    if (DOMAINDECOMP(integrator.cr))
    {
        /* Distribute the charge groups over the nodes from the master node */
        dd_partition_system(integrator.fplog, ir->init_step, integrator.cr, TRUE, 1,
                            integrator.state_global, integrator.top_global, ir,
                            state, &f, integrator.mdAtoms, top, integrator.fr,
                            integrator.vsite, integrator.constr,
                            integrator.nrnb, nullptr, FALSE);
        shouldCheckNumberOfBondedInteractions = true;
        update_realloc(upd, state->natoms);
    }

    //auto
    mdatoms = integrator.mdAtoms->mdatoms();

    // NOTE: The global state is no longer used at this point.
    // But integrator.state_global is still used as temporary storage space for writing
    // the global state to file and potentially for replica exchange.
    // (Global topology should persist.)

    update_mdatoms(mdatoms, state->lambda[efptMASS]);

    startingFromCheckpoint = continuationOptions.startedFromCheckpoint;


    if (MASTER(integrator.cr))
    {
        if (startingFromCheckpoint)
        {
            /* Update mdebin with energy history if appending to output files */
            if (continuationOptions.appendFiles)
            {
                restore_energyhistory_from_state(mdebin, integrator.observablesHistory->energyHistory.get());
            }
            else if (integrator.observablesHistory->energyHistory.get() != nullptr)
            {
                /* We might have read an energy history from checkpoint.
                 * As we are not appending, we want to restart the statistics.
                 * Free the allocated memory and reset the counts.
                 */
                integrator.observablesHistory->energyHistory = {};
            }
        }
        if (integrator.observablesHistory->energyHistory.get() == nullptr)
        {
            integrator.observablesHistory->energyHistory = std::unique_ptr<energyhistory_t>(new energyhistory_t {});
        }
        /* Set the initial energy history in state by updating once */
        update_energyhistory(integrator.observablesHistory->energyHistory.get(), mdebin);
    }

    /* Initialize constraints */
    if (integrator.constr && !DOMAINDECOMP(integrator.cr))
    {
        set_constraints(integrator.constr, top, ir, mdatoms, integrator.cr);
    }

    /* PME tuning is only supported in the Verlet scheme, with PME for
     * Coulomb. It is not supported with only LJ PME, or for
     * reruns. */
    bPMETune = (integrator.mdrunOptions.tunePme && EEL_PME(integrator.fr->ic->eeltype) &&
                !integrator.mdrunOptions.reproducible && ir->cutoff_scheme != ecutsGROUP);
    // Setup exported to element

    if (!ir->bContinuation)
    {
        if (state->flags & (1 << estV))
        {
            /* Set the velocities of integrator.vsites, shells and frozen atoms to zero */
            for (i = 0; i < mdatoms->homenr; i++)
            {
                if (mdatoms->ptype[i] == eptVSite ||
                    mdatoms->ptype[i] == eptShell)
                {
                    clear_rvec(state->v[i]);
                }
                else if (mdatoms->cFREEZE)
                {
                    for (m = 0; m < DIM; m++)
                    {
                        if (ir->opts.nFreeze[mdatoms->cFREEZE[i]][m])
                        {
                            state->v[i][m] = 0;
                        }
                    }
                }
            }
        }

        if (integrator.constr)
        {
            /* Constrain the initial coordinates and velocities */
            do_constrain_first(integrator.fplog, integrator.constr, ir, mdatoms, state,
                               integrator.cr, integrator.ms, integrator.nrnb, integrator.fr, top);
        }
    }

    /* Be REALLY careful about what flags you set here. You CANNOT assume
     * this is the first step, since we might be restarting from a checkpoint,
     * and in that case we should not do any modifications to the state.
     */
    bStopCM = (ir->comm_mode != ecmNO && !ir->bContinuation);

    if (continuationOptions.haveReadEkin)
    {
        restore_ekinstate_from_state(integrator.cr, ekind, &integrator.state_global->ekinstate);
    }

    cglo_flags = (CGLO_INITIALIZATION | CGLO_TEMPERATURE | CGLO_GSTAT
                  | (continuationOptions.haveReadEkin ? CGLO_READEKIN : 0));

    bSumEkinhOld = FALSE;
    /* To minimize communication, compute_globals computes the COM velocity
     * and the kinetic energy for the velocities without COM motion removed.
     * Thus to get the kinetic energy without the COM contribution, we need
     * to call compute_globals twice.
     */
    for (int cgloIteration = 0; cgloIteration < (bStopCM ? 2 : 1); cgloIteration++)
    {
        int cglo_flags_iteration = cglo_flags;
        if (bStopCM && cgloIteration == 0)
        {
            cglo_flags_iteration |= CGLO_STOPCM;
            cglo_flags_iteration &= ~CGLO_TEMPERATURE;
        }
        compute_globals(integrator.fplog, gstat, integrator.cr, ir, integrator.fr, ekind, state, mdatoms, integrator.nrnb, vcm,
                        nullptr, enerd, force_vir, shake_vir, total_vir, pres, mu_tot,
                        integrator.constr, &nullSignaller, state->box,
                        &totalNumberOfBondedInteractions, &bSumEkinhOld, cglo_flags_iteration
                        | (shouldCheckNumberOfBondedInteractions ? CGLO_CHECK_NUMBER_OF_BONDED_INTERACTIONS : 0));
    }
    checkNumberOfBondedInteractions(integrator.fplog, integrator.cr, totalNumberOfBondedInteractions,
                                    integrator.top_global, top, state,
                                    &shouldCheckNumberOfBondedInteractions);

    /* Calculate the initial half step temperature, and save the ekinh_old */
    if (!continuationOptions.startedFromCheckpoint)
    {
        for (i = 0; (i < ir->opts.ngtc); i++)
        {
            copy_mat(ekind->tcstat[i].ekinh, ekind->tcstat[i].ekinh_old);
        }
    }

    if (MASTER(integrator.cr))
    {
        if (!ir->bContinuation)
        {
            if (integrator.constr && ir->eConstrAlg == econtLINCS)
            {
                fprintf(integrator.fplog,
                        "RMS relative constraint deviation after constraining: %.2e\n",
                        constr_rmsd(integrator.constr));
            }
            if (EI_STATE_VELOCITY(ir->eI))
            {
                real temp = enerd->term[F_TEMP];
                {
                    /* Result of Ekin averaged over velocities of -half
                     * and +half step, while we only have -half step here.
                     */
                    temp *= 2;
                }
                fprintf(integrator.fplog, "Initial temperature: %g K\n", temp);
            }
        }

        // TODO: Clean up indenting
        {
            char tbuf[20];
            fprintf(stderr, "starting mdrun '%s'\n",
                    *(integrator.top_global->name));
            if (ir->nsteps >= 0)
            {
                sprintf(tbuf, "%8.1f", (ir->init_step+ir->nsteps)*ir->delta_t);
            }
            else
            {
                sprintf(tbuf, "%s", "infinite");
            }
            if (ir->init_step > 0)
            {
                fprintf(stderr, "%s steps, %s ps (continuing from step %s, %8.1f ps).\n",
                        gmx_step_str(ir->init_step+ir->nsteps, sbuf), tbuf,
                        gmx_step_str(ir->init_step, sbuf2),
                        ir->init_step*ir->delta_t);
            }
            else
            {
                fprintf(stderr, "%s steps, %s ps.\n",
                        gmx_step_str(ir->nsteps, sbuf), tbuf);
            }
        }
        fprintf(integrator.fplog, "\n");
    }

    walltime_accounting_start(integrator.walltime_accounting);
    wallcycle_start(integrator.wcycle, ewcRUN);
    print_start(integrator.fplog, integrator.cr, integrator.walltime_accounting, "mdrun");

    max_hours        = integrator.mdrunOptions.maximumHoursToRun;
    bFirstStep       = TRUE;
    /* Skip the first Nose-Hoover integration when we get the state from tpx */
    bInitStep        = !startingFromCheckpoint;
    bSumEkinhOld     = FALSE;

    nstSignalComm         = nstglobalcomm;

    step     = ir->init_step;
    step_rel = 0;

    /* and stop now if we should */
    bLastStep = (bLastStep || (ir->nsteps >= 0 && step_rel > ir->nsteps));
}

void gmx::DataManager::preStep()
{
    t_inputrec       *ir   = integrator.inputrec;
    /* Determine if this is a neighbor search step */
    bNStList = (ir->nstlist > 0  && step % ir->nstlist == 0);


    wallcycle_start(integrator.wcycle, ewcSTEP);

    bLastStep = (step_rel == ir->nsteps);
    t         = t0 + step*ir->delta_t;

    /* Stop Center of Mass motion */
    bStopCM = (ir->comm_mode != ecmNO && do_per_step(step, ir->nstcomm));

    /* Determine whether or not to do Neighbour Searching */
    bNS = (bFirstStep || bNStList);

    /* < 0 means stop at next step, > 0 means stop at next NS step */
    if ( (signals[eglsSTOPCOND].set < 0) ||
         ( (signals[eglsSTOPCOND].set > 0 ) && ( bNS || ir->nstlist == 0)))
    {
        bLastStep = TRUE;
    }

    /* do_log triggers energy and virial calculation. Because this leads
     * to different code paths, forces can be different. Thus for exact
     * continuation we should avoid extra log output.
     * Note that the || bLastStep can result in non-exact continuation
     * beyond the last step. But we don't consider that to be an issue.
     */
    do_log     = do_per_step(step, ir->nstlog) || (bFirstStep && !startingFromCheckpoint) || bLastStep;
    do_verbose = integrator.mdrunOptions.verbose &&
        (step % integrator.mdrunOptions.verboseStepPrintInterval == 0 || bFirstStep || bLastStep);

    if (MASTER(integrator.cr) && do_log)
    {
        print_ebin_header(integrator.fplog, step, t); /* can we improve the information printed here? */
    }

    clear_mat(force_vir);

    /* We write a checkpoint at this MD step when:
     * either at an NS step when we signalled through gs,
     * or at the last step (but not when we do not want confout),
     * but never at the first step or with rerun.
     */
    bCPT = (((signals[eglsCHKPT].set && (bNS || ir->nstlist == 0)) ||
             (bLastStep && integrator.mdrunOptions.writeConfout)) &&
            step > ir->init_step);
    if (bCPT)
    {
        signals[eglsCHKPT].set = 0;
    }

    /* Determine the energy and pressure:
     * at nstcalcenergy steps and at energy output steps (set below).
     */
    {
        bCalcEnerStep = do_per_step(step, ir->nstcalcenergy);
        bCalcVir      = bCalcEnerStep ||
            (ir->epc != epcNO && do_per_step(step, ir->nstpcouple));
    }
    bCalcEner = bCalcEnerStep;

    do_ene = (do_per_step(step, ir->nstenergy) || bLastStep);

    if (do_ene || do_log)
    {
        bCalcVir  = TRUE;
        bCalcEner = TRUE;
    }

    /* Do we need global communication ? */
    bGStat = (bCalcVir || bCalcEner || bStopCM ||
              do_per_step(step, nstglobalcomm));
}

void gmx::DataManager::postStep()
{
    t_inputrec       *ir   = integrator.inputrec;

    /* Print the remaining wall clock time for the run */
    if (isMasterSimMasterRank(integrator.ms, integrator.cr) &&
        (do_verbose || gmx_got_usr_signal()) &&
        !bPMETunePrinting)
    {
        print_time(stderr, integrator.walltime_accounting, step, ir, integrator.cr);
    }

    bFirstStep             = FALSE;
    bInitStep              = FALSE;
    startingFromCheckpoint = false;

    /* #######  SET VARIABLES FOR NEXT ITERATION IF THEY STILL NEED IT ###### */
    /* With all integrators, except VV, we need to retain the pressure
     * at the current step for coupling at the next step.
     */
    if ((state->flags & (1<<estPRES_PREV)) &&
        (bGStatEveryStep ||
         (ir->nstpcouple > 0 && step % ir->nstpcouple == 0)))
    {
        /* Store the pressure in t_state for pressure coupling
         * at the next MD step.
         */
        copy_mat(pres, state->pres_prev);
    }

    /* #######  END SET VARIABLES FOR NEXT ITERATION ###### */

    if ( (integrator.membed != nullptr) && (!bLastStep) )
    {
        rescale_membed(step_rel, integrator.membed, as_rvec_array(integrator.state_global->x.data()));
    }

    cycles = wallcycle_stop(integrator.wcycle, ewcSTEP);
    if (DOMAINDECOMP(integrator.cr) && integrator.wcycle)
    {
        dd_cycles_add(integrator.cr->dd, cycles, ddCyclStep);
    }

    step++;
    step_rel++;

    /* TODO make a counter-reset module */
    /* If it is time to reset counters, set a flag that remains
       true until counters actually get reset */
    if (step_rel == wcycle_get_reset_counters(integrator.wcycle) ||
        signals[eglsRESETCOUNTERS].set != 0)
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
        reset_all_counters(integrator.fplog, integrator.mdlog, integrator.cr, step, &step_rel, ir, integrator.wcycle, integrator.nrnb, integrator.walltime_accounting,
                           use_GPU(integrator.fr->nbv) ? integrator.fr->nbv : nullptr, integrator.fr->pmedata);
        wcycle_set_reset_counters(integrator.wcycle, -1);
        if (!thisRankHasDuty(integrator.cr, DUTY_PME))
        {
            /* Tell our PME node to reset its counters */
            gmx_pme_send_resetcounters(integrator.cr, step);
        }
        /* Correct max_hours for the elapsed time */
        max_hours                -= elapsed_time/(60.0*60.0);
        /* If mdrun -maxh -resethway was active, it can only trigger once */
        bResetCountersHalfMaxH    = FALSE; /* TODO move this to where signals[eglsRESETCOUNTERS].sig is set */
        /* Reset can only happen once, so clear the triggering flag. */
        signals[eglsRESETCOUNTERS].set = 0;
    }
}


gmx::DataManager::~DataManager()
{
    t_inputrec       *ir   = integrator.inputrec;

    /* Closing TNG files can include compressing data. Therefore it is good to do that
     * before stopping the time measurements. */
    mdoutf_tng_close(outf);

    /* Stop measuring walltime */
    walltime_accounting_end(integrator.walltime_accounting);

    if (!thisRankHasDuty(integrator.cr, DUTY_PME))
    {
        /* Tell the PME only node to finish */
        gmx_pme_send_finish(integrator.cr);
    }

    if (MASTER(integrator.cr))
    {
        if (ir->nstcalcenergy > 0)
        {
            print_ebin(mdoutf_get_fp_ene(outf), FALSE, FALSE, FALSE,
                       integrator.fplog, step, t,
                       eprAVER, mdebin, integrator.fcd,
                       groups, &(ir->opts), ir->awh);
        }
    }
    done_mdebin(mdebin);
    done_mdoutf(outf);

    done_shellfc(integrator.fplog, shellfc, step_rel);

    // Clean up swapcoords
    if (ir->eSwapCoords != eswapNO)
    {
        finish_swapcoords(ir->swap);
    }

    walltime_accounting_set_nsteps_done(integrator.walltime_accounting, step_rel);
    if (step_rel >= wcycle_get_reset_counters(integrator.wcycle) &&
        signals[eglsRESETCOUNTERS].set == 0 &&
        !bResetCountersHalfMaxH)
    {
        walltime_accounting_set_valid_finish(integrator.walltime_accounting);
    }

    destroy_enerdata(enerd);
    sfree(enerd);
    mdAlgorithmsTearDownAtomData(integrator.fr->bonded_threading, top);
    integrator.fr->bonded_threading = nullptr;
    sfree(top);
}

//! Resets all the counters.
void gmx::DataManager::reset_all_counters(FILE *fplog, const gmx::MDLogger &mdlog, t_commrec *cr,
                                          gmx_int64_t step,
                                          gmx_int64_t *step_rel, t_inputrec *ir,
                                          gmx_wallcycle_t wcycle, t_nrnb *nrnb,
                                          gmx_walltime_accounting_t walltime_accounting,
                                          struct nonbonded_verlet_t *nbv,
                                          struct gmx_pme_t *pme)
{
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
    ir->init_step += *step_rel;
    ir->nsteps    -= *step_rel;
    *step_rel      = 0;
    wallcycle_start(wcycle, ewcRUN);
    walltime_accounting_start(walltime_accounting);
    print_date_and_time(fplog, cr->nodeid, "Restarted time", gmx_gettime());
}
