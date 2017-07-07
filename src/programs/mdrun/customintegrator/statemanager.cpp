#include "statemanager.h"

#include "gromacs/domdec/domdec.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme-load-balancing.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/imd/imd.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/compute_io.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/mdebin.h"
#include "gromacs/mdlib/md_support.h"
#include "gromacs/mdlib/mdoutf.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/mdsetup.h"
#include "gromacs/mdlib/nbnxn_gpu_data_mgmt.h"
#include "gromacs/mdlib/shellfc.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdlib/tgroup.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/energyhistory.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/swap/swapcoords.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#include "programs/mdrun/deform.h"
#include "programs/mdrun/membed.h"
#include "programs/mdrun/repl_ex.h"
#include "helperFunctionsTemporary.h"

void reset_all_counters(FILE *fplog, const gmx::MDLogger &mdlog, t_commrec *cr,
                               gmx_int64_t step,
                               gmx_int64_t *step_rel, t_inputrec *ir,
                               gmx_wallcycle_t wcycle, t_nrnb *nrnb,
                               gmx_walltime_accounting_t walltime_accounting,
                               struct nonbonded_verlet_t *nbv)
{
    char sbuf[STEPSTRSIZE];

    /* Reset all the counters related to performance over the run */
    GMX_LOG(mdlog.warning).asParagraph().appendTextFormatted(
            "step %s: resetting all time and cycle counters",
            gmx_step_str(step, sbuf));

    if (use_GPU(nbv))
    {
        nbnxn_gpu_reset_timings(nbv);
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

void StateManager::loopSetup()
{
	    /* Check for special mdrun options */
    if (Flags & MD_RESETCOUNTERSHALFWAY)
    {
        if (ir->nsteps > 0)
        {
            /* Signal to reset the counters half the simulation steps. */
            wcycle_set_reset_counters(wcycle, ir->nsteps/2);
        }
        /* Signal to reset the counters halfway the simulation time. */
        bResetCountersHalfMaxH = (max_hours > 0);
    }

    /* md-vv uses averaged full step velocities for T-control
       md uses averaged half step kinetic energies to determine temperature unless defined otherwise by GMX_EKIN_AVE_VEL; */
    bTrotter = (EI_VV(ir->eI) && (inputrecNptTrotter(ir) || inputrecNphTrotter(ir) || inputrecNvtTrotter(ir)));

    nstglobalcomm   = check_nstglobalcomm(*mdlog, nstglobalcomm, ir);
    bGStatEveryStep = (nstglobalcomm == 1);

    groups = &top_global->groups;

    if (ir->eSwapCoords != eswapNO)
    {
        /* Initialize ion swapping code */
        init_swapcoords(fplog, bVerbose, ir, opt2fn_master("-swap", nfile, fnm, cr),
                        top_global, as_rvec_array(state_global->x.data()), state_global->box, &state_global->swapstate, cr, oenv, Flags);
    }

    /* Initial values */
    init_md(fplog, cr, ir, oenv, &t, &t0, &state_global->lambda,
            &(state_global->fep_state), lam0,
            nrnb, top_global, &upd,
            nfile, fnm, &outf, &mdebin,
            force_vir, shake_vir, mu_tot, &bSimAnn, &vcm, Flags, wcycle);

    clear_mat(total_vir);
    clear_mat(pres);
    /* Energy terms and groups */
    snew(enerd, 1);
    init_enerdata(top_global->groups.grps[egcENER].nr, ir->fepvals->n_lambda,
                  enerd);
    if (!DOMAINDECOMP(cr))
    {
        f.resize(top_global->natoms + 1);
    }

    /* Kinetic energy data */
    snew(ekind, 1);
    init_ekindata(fplog, top_global, &(ir->opts), ekind);
    /* Copy the cos acceleration to the groups struct */
    ekind->cosacc.cos_accel = ir->cos_accel;

    gstat = global_stat_init(ir);

    /* Check for polarizable models and flexible constraints */
    shellfc = init_shell_flexcon(fplog,
                                 top_global, n_flexible_constraints(constr),
                                 ir->nstcalcenergy, DOMAINDECOMP(cr));

    if (shellfc && ir->nstcalcenergy != 1)
    {
        gmx_fatal(FARGS, "You have nstcalcenergy set to a value (%d) that is different from 1.\nThis is not supported in combinations with shell particles.\nPlease make a new tpr file.", ir->nstcalcenergy);
    }
    if (shellfc && DOMAINDECOMP(cr))
    {
        gmx_fatal(FARGS, "Shell particles are not implemented with domain decomposition, use a single rank");
    }

    if (inputrecDeform(ir))
    {
        tMPI_Thread_mutex_lock(&deform_init_box_mutex);
        set_deform_reference_box(upd,
                                 deform_init_init_step_tpx,
                                 deform_init_box_tpx);
        tMPI_Thread_mutex_unlock(&deform_init_box_mutex);
    }

    {
        double io = compute_io(ir, top_global->natoms, groups, mdebin->ebin->nener, 1);
        if ((io > 2000) && MASTER(cr))
        {
            fprintf(stderr,
                    "\nWARNING: This run will generate roughly %.0f Mb of data\n\n",
                    io);
        }
    }

    if (DOMAINDECOMP(cr))
    {
        top = dd_init_local_top(top_global);

        stateInstance = std::unique_ptr<t_state>(new t_state {});
        state         = stateInstance.get();
        dd_init_local_state(cr->dd, state_global, state);
    }
    else
    {
        /* Copy the pointer to the global state */
        state = state_global;

        snew(top, 1);
        mdAlgorithmsSetupAtomData(cr, ir, top_global, top, fr,
                                  &graph, mdatoms, vsite, shellfc);

        update_realloc(upd, state->natoms);
    }

    /* Set up interactive MD (IMD) */
    init_IMD(ir, cr, top_global, fplog, ir->nstcalcenergy, as_rvec_array(state_global->x.data()),
             nfile, fnm, oenv, imdport, Flags);

    if (DOMAINDECOMP(cr))
    {
        /* Distribute the charge groups over the nodes from the master node */
        dd_partition_system(fplog, ir->init_step, cr, TRUE, 1,
                            state_global, top_global, ir,
                            state, &f, mdatoms, top, fr,
                            vsite, constr,
                            nrnb, nullptr, FALSE);
        shouldCheckNumberOfBondedInteractions = true;
        update_realloc(upd, state->natoms);
    }

    update_mdatoms(mdatoms, state->lambda[efptMASS]);

    startingFromCheckpoint = Flags & MD_STARTFROMCPT;

    if (ir->bExpanded)
    {
        init_expanded_ensemble(startingFromCheckpoint, ir, state->dfhist);
    }

    if (MASTER(cr))
    {
        if (startingFromCheckpoint)
        {
            /* Update mdebin with energy history if appending to output files */
            if (Flags & MD_APPENDFILES)
            {
                restore_energyhistory_from_state(mdebin, energyHistory);
            }
            else
            {
                /* We might have read an energy history from checkpoint,
                 * free the allocated memory and reset the counts.
                 */
                *energyHistory = {};
            }
        }
        /* Set the initial energy history in state by updating once */
        update_energyhistory(energyHistory, mdebin);
    }

    /* Initialize constraints */
    if (constr && !DOMAINDECOMP(cr))
    {
        set_constraints(constr, top, ir, mdatoms, cr);
    }

    if (repl_ex_nst > 0 && MASTER(cr))
    {
        repl_ex = init_replica_exchange(fplog, cr->ms, state_global, ir,
                                        repl_ex_nst, repl_ex_nex, repl_ex_seed);
    }

    /* PME tuning is only supported with PME for Coulomb. Is is not supported
     * with only LJ PME.
     */
    bPMETune = ((Flags & MD_TUNEPME) && EEL_PME(fr->eeltype) &&
                !(Flags & MD_REPRODUCIBLE));
    if (bPMETune)
    {
        pme_loadbal_init(&pme_loadbal, cr, *mdlog, ir, state->box,
                         fr->ic, fr->pmedata, use_GPU(fr->nbv),
                         &bPMETunePrinting);
    }

    if (!ir->bContinuation)
    {
        if (state->flags & (1 << estV))
        {
            /* Set the velocities of vsites, shells and frozen atoms to zero */
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

        if (vsite)
        {
            /* Construct the virtual sites for the initial configuration */
            construct_vsites(vsite, as_rvec_array(state->x.data()), ir->delta_t, nullptr,
                             top->idef.iparams, top->idef.il,
                             fr->ePBC, fr->bMolPBC, cr, state->box);
        }
    }

    if (ir->efep != efepNO)
    {
        /* Set free energy calculation frequency as the greatest common
         * denominator of nstdhdl and repl_ex_nst. */
        nstfep = ir->fepvals->nstdhdl;
        if (ir->bExpanded)
        {
            nstfep = gmx_greatest_common_divisor(ir->expandedvals->nstexpanded, nstfep);
        }
        if (repl_ex_nst > 0)
        {
            nstfep = gmx_greatest_common_divisor(repl_ex_nst, nstfep);
        }
    }

    /* Be REALLY careful about what flags you set here. You CANNOT assume
     * this is the first step, since we might be restarting from a checkpoint,
     * and in that case we should not do any modifications to the state.
     */
    bStopCM = (ir->comm_mode != ecmNO && !ir->bContinuation);

    if (Flags & MD_READ_EKIN)
    {
        restore_ekinstate_from_state(cr, ekind, &state_global->ekinstate);
    }

    cglo_flags = (CGLO_TEMPERATURE | CGLO_GSTAT
                  | (bStopCM ? CGLO_STOPCM : 0)
                  | (EI_VV(ir->eI) ? CGLO_PRESSURE : 0)
                  | (EI_VV(ir->eI) ? CGLO_CONSTRAINT : 0)
                  | ((Flags & MD_READ_EKIN) ? CGLO_READEKIN : 0));

    bSumEkinhOld = FALSE;
    compute_globals(fplog, gstat, cr, ir, fr, ekind, state, mdatoms, nrnb, vcm,
                    nullptr, enerd, force_vir, shake_vir, total_vir, pres, mu_tot,
                    constr, &nullSignaller, state->box,
                    &totalNumberOfBondedInteractions, &bSumEkinhOld, cglo_flags
                    | (shouldCheckNumberOfBondedInteractions ? CGLO_CHECK_NUMBER_OF_BONDED_INTERACTIONS : 0));
    checkNumberOfBondedInteractions(fplog, cr, totalNumberOfBondedInteractions,
                                    top_global, top, state,
                                    &shouldCheckNumberOfBondedInteractions);
    
    /* Calculate the initial half step temperature, and save the ekinh_old */
    if (!(Flags & MD_STARTFROMCPT))
    {
        for (i = 0; (i < ir->opts.ngtc); i++)
        {
            copy_mat(ekind->tcstat[i].ekinh, ekind->tcstat[i].ekinh_old);
        }
    }
    if (ir->eI != eiVV)
    {
        enerd->term[F_TEMP] *= 2; /* result of averages being done over previous and current step,
                                     and there is no previous step */
    }

    /* need to make an initiation call to get the Trotter variables set, as well as other constants for non-trotter
       temperature control */
    trotter_seq = init_npt_vars(ir, state, &MassQ, bTrotter);

    if (MASTER(cr))
    {
        if (constr && !ir->bContinuation && ir->eConstrAlg == econtLINCS)
        {
            fprintf(fplog,
                    "RMS relative constraint deviation after constraining: %.2e\n",
                    constr_rmsd(constr));
        }
        if (EI_STATE_VELOCITY(ir->eI))
        {
            fprintf(fplog, "Initial temperature: %g K\n", enerd->term[F_TEMP]);
        }

        {
            char tbuf[20];
            fprintf(stderr, "starting mdrun '%s'\n",
                    *(top_global->name));
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
        fprintf(fplog, "\n");
    }

    walltime_accounting_start(walltime_accounting);
    wallcycle_start(wcycle, ewcRUN);
    print_start(fplog, cr, walltime_accounting, "mdrun");

    /* safest point to do file checkpointing is here.  More general point would be immediately before integrator call */
#ifdef GMX_FAHCORE
    chkpt_ret = fcCheckPointParallel( cr->nodeid,
                                      NULL, 0);
    if (chkpt_ret == 0)
    {
        gmx_fatal( 3, __FILE__, __LINE__, "Checkpoint error on step %d\n", 0 );
    }
#endif

    /* Initial values for booleans, that might later change during integration loop */
    bFirstStep = TRUE;
    /* Skip the first Nose-Hoover integration when we get the state from tpx */
    bInitStep        = !startingFromCheckpoint || EI_VV(ir->eI);
    bSumEkinhOld     = FALSE;
    bExchanged       = FALSE;
    bNeedRepartition = FALSE;
    // TODO This implementation of ensemble orientation restraints is nasty because
    // a user can't just do multi-sim with single-sim orientation restraints.
    bUsingEnsembleRestraints = (fcd->disres.nsystems > 1) || (cr->ms && fcd->orires.nr);

    {
        // Replica exchange and ensemble restraints need all
        // simulations to remain synchronized, so they need
        // checkpoints and stop conditions to act on the same step, so
        // the propagation of such signals must take place between
        // simulations, not just within simulations.
        bool checkpointIsLocal    = (repl_ex_nst <= 0) && !bUsingEnsembleRestraints;
        bool stopConditionIsLocal = (repl_ex_nst <= 0) && !bUsingEnsembleRestraints;
        bool resetCountersIsLocal = true;
        signals[eglsCHKPT]         = gmx::SimulationSignal(checkpointIsLocal);
        signals[eglsSTOPCOND]      = gmx::SimulationSignal(stopConditionIsLocal);
        signals[eglsRESETCOUNTERS] = gmx::SimulationSignal(resetCountersIsLocal);
    }

    step     = ir->init_step;
    step_rel = 0;

    // TODO extract this to new multi-simulation module
    if (MASTER(cr) && MULTISIM(cr) && (repl_ex_nst <= 0 ))
    {
        if (!multisim_int_all_are_equal(cr->ms, ir->nsteps))
        {
            GMX_LOG(mdlog->warning).appendText(
                    "Note: The number of steps is not consistent across multi simulations,\n"
                    "but we are proceeding anyway!");
        }
        if (!multisim_int_all_are_equal(cr->ms, ir->init_step))
        {
            GMX_LOG(mdlog->warning).appendText(
                    "Note: The initial step is not consistent across multi simulations,\n"
                    "but we are proceeding anyway!");
        }
    }

    /* and stop now if we should */
    bLastStep = (bLastStep || (ir->nsteps >= 0 && step_rel > ir->nsteps));
    copy_mat(state->box, lastbox);
}

void StateManager::stepSetup()
{
	    /* Determine if this is a neighbor search step */
    bNStList = (ir->nstlist > 0  && step % ir->nstlist == 0);

    if (bPMETune && bNStList)
    {
        /* PME grid + cut-off optimization with GPUs or PME nodes */
        pme_loadbal_do(pme_loadbal, cr,
                       (bVerbose && MASTER(cr)) ? stderr : nullptr,
                       fplog, *mdlog,
                       ir, fr, state,
                       wcycle,
                       step, step_rel,
                       &bPMETunePrinting);
    }

    wallcycle_start(wcycle, ewcSTEP);

    {
        bLastStep = (step_rel == ir->nsteps);
        t         = t0 + step*ir->delta_t;
    }

    // TODO Refactor this, so that nstfep does not need a default value of zero
    if (ir->efep != efepNO || ir->bSimTemp)
    {
        /* find and set the current lambdas. */

        set_current_lambdas(step, ir->fepvals, FALSE, nullptr, state_global, state, lam0);
        bDoDHDL      = do_per_step(step, ir->fepvals->nstdhdl);
        bDoFEP       = ((ir->efep != efepNO) && do_per_step(step, nstfep));
        bDoExpanded  = (do_per_step(step, ir->expandedvals->nstexpanded)
                        && (ir->bExpanded) && (step > 0) && (!startingFromCheckpoint));
    }

    bDoReplEx = ((repl_ex_nst > 0) && (step > 0) && !bLastStep &&
                 do_per_step(step, repl_ex_nst));

    if (bSimAnn)
    {
        update_annealing_target_temp(ir, t, upd);
    }

    /* Stop Center of Mass motion */
    bStopCM = (ir->comm_mode != ecmNO && do_per_step(step, ir->nstcomm));

    {
        /* Determine whether or not to do Neighbour Searching */
        bNS = (bFirstStep || bNStList || bExchanged || bNeedRepartition);
    }

    /* < 0 means stop at next step, > 0 means stop at next NS step */
    if ( (signals[eglsSTOPCOND].set < 0) ||
         ( (signals[eglsSTOPCOND].set > 0 ) && ( bNS || ir->nstlist == 0)))
    {
        bLastStep = TRUE;
    }

    /* Determine whether or not to update the Born radii if doing GB */
    bBornRadii = bFirstStep;
    if (ir->implicit_solvent && (step % ir->nstgbradii == 0))
    {
        bBornRadii = TRUE;
    }

    /* do_log triggers energy and virial calculation. Because this leads
     * to different code paths, forces can be different. Thus for exact
     * continuation we should avoid extra log output.
     * Note that the || bLastStep can result in non-exact continuation
     * beyond the last step. But we don't consider that to be an issue.
     */
    do_log     = do_per_step(step, ir->nstlog) || (bFirstStep && !startingFromCheckpoint) || bLastStep;
    do_verbose = bVerbose &&
        (step % stepout == 0 || bFirstStep || bLastStep);

    if (bNS && !(bFirstStep && ir->bContinuation))
    {
        {
            bMasterState = FALSE;
            /* Correct the new box if it is too skewed */
            if (inputrecDynamicBox(ir))
            {
                if (correct_box(fplog, step, state->box, graph))
                {
                    bMasterState = TRUE;
                }
            }
            if (DOMAINDECOMP(cr) && bMasterState)
            {
                dd_collect_state(cr->dd, state, state_global);
            }
        }

        if (DOMAINDECOMP(cr))
        {
            /* Repartition the domain decomposition */
            dd_partition_system(fplog, step, cr,
                                bMasterState, nstglobalcomm,
                                state_global, top_global, ir,
                                state, &f, mdatoms, top, fr,
                                vsite, constr,
                                nrnb, wcycle,
                                do_verbose && !bPMETunePrinting);
            shouldCheckNumberOfBondedInteractions = true;
            update_realloc(upd, state->natoms);
        }
    }

    if (MASTER(cr) && do_log)
    {
        print_ebin_header(fplog, step, t); /* can we improve the information printed here? */
    }

    if (ir->efep != efepNO)
    {
        update_mdatoms(mdatoms, state->lambda[efptMASS]);
    }

    if (bExchanged)
    {

        /* We need the kinetic energy at minus the half step for determining
         * the full step kinetic energy and possibly for T-coupling.*/
        /* This may not be quite working correctly yet . . . . */
        compute_globals(fplog, gstat, cr, ir, fr, ekind, state, mdatoms, nrnb, vcm,
                        wcycle, enerd, nullptr, nullptr, nullptr, nullptr, mu_tot,
                        constr, &nullSignaller, state->box,
                        &totalNumberOfBondedInteractions, &bSumEkinhOld,
                        CGLO_GSTAT | CGLO_TEMPERATURE | CGLO_CHECK_NUMBER_OF_BONDED_INTERACTIONS);
        checkNumberOfBondedInteractions(fplog, cr, totalNumberOfBondedInteractions,
                                        top_global, top, state,
                                        &shouldCheckNumberOfBondedInteractions);
    }

    /* We write a checkpoint at this MD step when:
     * either at an NS step when we signalled through gs,
     * or at the last step (but not when we do not want confout),
     * but never at the first step.
     */
    bCPT = (((signals[eglsCHKPT].set && (bNS || ir->nstlist == 0)) ||
             (bLastStep && (Flags & MD_CONFOUT))) &&
            step > ir->init_step);
    if (bCPT)
    {
        signals[eglsCHKPT].set = 0;
    }

    /* Determine the energy and pressure:
     * at nstcalcenergy steps and at energy output steps (set below).
     */
    if (EI_VV(ir->eI) && (!bInitStep))
    {
        /* for vv, the first half of the integration actually corresponds
           to the previous step.  bCalcEner is only required to be evaluated on the 'next' step,
           but the virial needs to be calculated on both the current step and the 'next' step. Future
           reorganization may be able to get rid of one of the bCalcVir=TRUE steps. */

        /* TODO: This is probably not what we want, we will write to energy file one step after nstcalcenergy steps. */
        bCalcEnerStep = do_per_step(step - 1, ir->nstcalcenergy);
        bCalcVir      = bCalcEnerStep ||
            (ir->epc != epcNO && (do_per_step(step, ir->nstpcouple) || do_per_step(step-1, ir->nstpcouple)));
    }
    else
    {
        bCalcEnerStep = do_per_step(step, ir->nstcalcenergy);
        bCalcVir      = bCalcEnerStep ||
            (ir->epc != epcNO && do_per_step(step, ir->nstpcouple));
    }
    bCalcEner = bCalcEnerStep;

    do_ene = (do_per_step(step, ir->nstenergy) || bLastStep);

    if (do_ene || do_log || bDoReplEx)
    {
        bCalcVir  = TRUE;
        bCalcEner = TRUE;
    }

    /* Do we need global communication ? */
    bGStat = (bCalcVir || bCalcEner || bStopCM ||
              do_per_step(step, nstglobalcomm) ||
              (EI_VV(ir->eI) && inputrecNvtTrotter(ir) && do_per_step(step-1, nstglobalcomm)));

    force_flags = (GMX_FORCE_STATECHANGED |
                   (inputrecDynamicBox(ir) ? GMX_FORCE_DYNAMICBOX : 0) |
                   GMX_FORCE_ALLFORCES |
                   (bCalcVir ? GMX_FORCE_VIRIAL : 0) |
                   (bCalcEner ? GMX_FORCE_ENERGY : 0) |
                   (bDoFEP ? GMX_FORCE_DHDL : 0)
                   );
}

void StateManager::stepTeardown()
{
    if (bDoExpanded)
    {
        /* Have to do this part _after_ outputting the logfile and the edr file */
        /* Gets written into the state at the beginning of next loop*/
        state->fep_state = lamnew;
    }
    /* Print the remaining wall clock time for the run */
    if (MULTIMASTER(cr) &&
        (do_verbose || gmx_got_usr_signal()) &&
        !bPMETunePrinting)
    {
        if (shellfc)
        {
            fprintf(stderr, "\n");
        }
        print_time(stderr, walltime_accounting, step, ir, cr);
    }

    /* Ion/water position swapping.
     * Not done in last step since trajectory writing happens before this call
     * in the MD loop and exchanges would be lost anyway. */
    bNeedRepartition = FALSE;
    if ((ir->eSwapCoords != eswapNO) && (step > 0) && !bLastStep &&
        do_per_step(step, ir->swap->nstswap))
    {
        bNeedRepartition = do_swapcoords(cr, step, t, ir, wcycle,
                                         as_rvec_array(state->x.data()),
                                         state->box,
                                         MASTER(cr) && bVerbose, FALSE);

        if (bNeedRepartition && DOMAINDECOMP(cr))
        {
            dd_collect_state(cr->dd, state, state_global);
        }
    }

    /* Replica exchange */
    bExchanged = FALSE;
    if (bDoReplEx)
    {
        bExchanged = replica_exchange(fplog, cr, repl_ex,
                                      state_global, enerd,
                                      state, step, t);
    }

    if ( (bExchanged || bNeedRepartition) && DOMAINDECOMP(cr) )
    {
        dd_partition_system(fplog, step, cr, TRUE, 1,
                            state_global, top_global, ir,
                            state, &f, mdatoms, top, fr,
                            vsite, constr,
                            nrnb, wcycle, FALSE);
        shouldCheckNumberOfBondedInteractions = true;
        update_realloc(upd, state->natoms);
    }

    bFirstStep             = FALSE;
    bInitStep              = FALSE;
    startingFromCheckpoint = FALSE;

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

    if ( (membed != nullptr) && (!bLastStep) )
    {
        rescale_membed(step_rel, membed, as_rvec_array(state_global->x.data()));
    }

    cycles = wallcycle_stop(wcycle, ewcSTEP);
    if (DOMAINDECOMP(cr) && wcycle)
    {
        dd_cycles_add(cr->dd, cycles, ddCyclStep);
    }

    {
        /* increase the MD step number */
        step++;
        step_rel++;
    }

    /* TODO make a counter-reset module */
    /* If it is time to reset counters, set a flag that remains
       true until counters actually get reset */
    if (step_rel == wcycle_get_reset_counters(wcycle) ||
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
        reset_all_counters(fplog, *mdlog, cr, step, &step_rel, ir, wcycle, nrnb, walltime_accounting,
                           use_GPU(fr->nbv) ? fr->nbv : nullptr);
        wcycle_set_reset_counters(wcycle, -1);
        if (!(cr->duty & DUTY_PME))
        {
            /* Tell our PME node to reset its counters */
            gmx_pme_send_resetcounters(cr, step);
        }
        /* Correct max_hours for the elapsed time */
        max_hours                -= elapsed_time/(60.0*60.0);
        /* If mdrun -maxh -resethway was active, it can only trigger once */
        bResetCountersHalfMaxH    = FALSE; /* TODO move this to where signals[eglsRESETCOUNTERS].sig is set */
        /* Reset can only happen once, so clear the triggering flag. */
        signals[eglsRESETCOUNTERS].set = 0;
    }

    /* If bIMD is TRUE, the master updates the IMD energy record and sends positions to VMD client */
    IMD_prep_energies_send_positions(ir->bIMD && MASTER(cr), bIMDstep, ir->imd, enerd, step, bCalcEner, wcycle);
}

void StateManager::loopTeardown()
{
	/* Closing TNG files can include compressing data. Therefore it is good to do that
	 * before stopping the time measurements. */
	mdoutf_tng_close(outf);

	/* Stop measuring walltime */
	walltime_accounting_end(walltime_accounting);

	if (!(cr->duty & DUTY_PME))
	{
	    /* Tell the PME only node to finish */
	    gmx_pme_send_finish(cr);
	}

	if (MASTER(cr))
	{
	    if (ir->nstcalcenergy > 0)
	    {
	        print_ebin(mdoutf_get_fp_ene(outf), FALSE, FALSE, FALSE, fplog, step, t,
	                   eprAVER, mdebin, fcd, groups, &(ir->opts));
	    }
	}

	done_mdoutf(outf, ir);

	if (bPMETune)
	{
	    pme_loadbal_done(pme_loadbal, fplog, *mdlog, use_GPU(fr->nbv));
	}

	done_shellfc(fplog, shellfc, step_rel);

	if (repl_ex_nst > 0 && MASTER(cr))
	{
	    print_replica_exchange_statistics(fplog, repl_ex);
	}

	// Clean up swapcoords
	if (ir->eSwapCoords != eswapNO)
	{
	    finish_swapcoords(ir->swap);
	}

	/* IMD cleanup, if bIMD is TRUE. */
	IMD_finalize(ir->bIMD, ir->imd);

	walltime_accounting_set_nsteps_done(walltime_accounting, step_rel);
}
