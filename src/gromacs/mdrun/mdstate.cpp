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
/*! \libinternal \file
 *
 * \brief This file defines the state of an MD run
 *
 * \author Pascal Merz <pascal.merz@colorado.edu>
 * \inlibraryapi
 */

#include "mdstate.h"

using namespace gmx;

SimulationSetup::SimulationSetup(std::shared_ptr<Integrator> &integrator) :
        bGStat(false),
        bCalcVir(false),
        bCalcEnerStep(false),
        bCalcEner(false),
        bNS(false),
        bNStList(false),
        bSimAnn(false),
        bLastStep(false),
        bDoDHDL(false),
        bDoFEP(false),
        bDoExpanded(false),
        do_ene(false),
        do_log(false),
        do_verbose(false),
        bRerunWarnNoV(false),
        bMasterState(false),
        bDoReplEx(false),
        bTemp(false),
        bPres(false),
        bPMETunePrinting(false),
        bIMDstep(false),
        resetCountersIsLocal(true)
{
    const gmx_multisim_t            *ms = integrator->ms;
    const MDLogger                  &mdlog = integrator->mdlog;
    const MdrunOptions              &mdrunOptions = integrator->mdrunOptions;
    t_inputrec                      *ir = integrator->inputrec;
    const t_fcdata                  *fcd = integrator->fcd;
    const t_forcerec                *fr = integrator->fr;
    const ReplicaExchangeParameters &replExParams = integrator->replExParams;

    /* md-vv uses averaged full step velocities for T-control
       md-vv-avek uses averaged half step velocities for T-control (but full step ekin for P control)
       md uses averaged half step kinetic energies to determine temperature unless defined otherwise by GMX_EKIN_AVE_VEL; */
    bTrotter = (EI_VV(ir->eI) && (inputrecNptTrotter(ir) || inputrecNphTrotter(ir) || inputrecNvtTrotter(ir)));

    bRerunMD      = mdrunOptions.rerun;
    nstglobalcomm = mdrunOptions.globalCommunicationInterval;

    if (bRerunMD)
    {
    /* Since we don't know if the frames read are related in any way,
     * rebuild the neighborlist at every step.
     */
    ir->nstlist       = 1;
    ir->nstcalcenergy = 1;
    nstglobalcomm     = 1;
    }

    nstglobalcomm   = check_nstglobalcomm(mdlog, nstglobalcomm, ir);
    bGStatEveryStep = (nstglobalcomm == 1);

    startingFromCheckpoint = mdrunOptions.continuationOptions.startedFromCheckpoint;

    useReplicaExchange = (replExParams.exchangeInterval > 0);

    /* PME tuning is only supported in the Verlet scheme, with PME for
     * Coulomb. It is not supported with only LJ PME, or for
     * reruns. */
    bPMETune = (mdrunOptions.tunePme && EEL_PME(fr->ic->eeltype) && !bRerunMD &&
                !mdrunOptions.reproducible && ir->cutoff_scheme != ecutsGROUP);

    if (ir->efep != efepNO)
    {
    /* Set free energy calculation frequency as the greatest common
     * denominator of nstdhdl and repl_ex_nst. */
    nstfep = ir->fepvals->nstdhdl;
    if (ir->bExpanded)
    {
    nstfep = gmx_greatest_common_divisor(ir->expandedvals->nstexpanded, nstfep);
    }
    if (useReplicaExchange)
    {
    nstfep = gmx_greatest_common_divisor(replExParams.exchangeInterval, nstfep);
    }
    }

    /* Be REALLY careful about what flags you set here. You CANNOT assume
     * this is the first step, since we might be restarting from a checkpoint,
     * and in that case we should not do any modifications to the state.
     */
    bStopCM = (ir->comm_mode != ecmNO && !ir->bContinuation);


    if (bRerunMD)
    {
    if (getenv("GMX_FORCE_UPDATE"))
    {
    bForceUpdate = true;
    }
    }

    bFirstStep       = true;
    /* Skip the first Nose-Hoover integration when we get the state from tpx */
    bInitStep        = !startingFromCheckpoint || EI_VV(ir->eI);
    bSumEkinhOld     = true;
    bExchanged       = true;
    bNeedRepartition = true;

    simulationsShareState = true;
    nstSignalComm         = nstglobalcomm;
    {
        // TODO This implementation of ensemble orientation restraints is nasty because
        // a user can't just do multi-sim with single-sim orientation restraints.
        bool usingEnsembleRestraints = (fcd->disres.nsystems > 1) || ((ms != nullptr) && (fcd->orires.nr != 0));
        bool awhUsesMultiSim         = (ir->bDoAwh && ir->awhParams->shareBiasMultisim && (ms != nullptr));

        // Replica exchange, ensemble restraints and AWH need all
        // simulations to remain synchronized, so they need
        // checkpoints and stop conditions to act on the same step, so
        // the propagation of such signals must take place between
        // simulations, not just within simulations.
        // TODO: Make algorithm initializers set these flags.
        simulationsShareState     = useReplicaExchange || usingEnsembleRestraints || awhUsesMultiSim;

        if (simulationsShareState)
        {
        // Inter-simulation signal communication does not need to happen
        // often, so we use a minimum of 200 steps to reduce overhead.
        const int c_minimumInterSimulationSignallingInterval = 200;
        nstSignalComm = ((c_minimumInterSimulationSignallingInterval + nstglobalcomm - 1)/nstglobalcomm)*nstglobalcomm;
        }
    }

    step     = ir->init_step;
    step_rel = 0;

    // TODO extract this to new multi-simulation module
    if (MASTER(cr) && isMultiSim(ms) && !useReplicaExchange)
    {
        if (!multisim_int_all_are_equal(ms, ir->nsteps))
        {
            GMX_LOG(mdlog.warning).appendText(
                    "Note: The number of steps is not consistent across multi simulations,\n"
                    "but we are proceeding anyway!");
        }
        if (!multisim_int_all_are_equal(ms, ir->init_step))
        {
            GMX_LOG(mdlog.warning).appendText(
                    "Note: The initial step is not consistent across multi simulations,\n"
                    "but we are proceeding anyway!");
        }
    }
}

void MDState::init(
        std::shared_ptr<Integrator>      &integrator,
        std::shared_ptr<SimulationSetup> &setup)
{
    // Data from integrator (using alias here to reduce number of changed lines(
    FILE *fplog = integrator->fplog;
    t_commrec *cr = integrator->cr;
    const gmx_multisim_t *ms = integrator->ms;
    const MDLogger &mdlog = integrator->mdlog;
    const int nfile = integrator->nfile;
    const t_filenm *fnm = integrator->fnm;
    const gmx_output_env_t *oenv = integrator->oenv;
    const MdrunOptions &mdrunOptions = integrator->mdrunOptions;
    gmx_vsite_t *vsite = integrator->vsite;
    Constraints *constr = integrator->constr;
    // gmx_enfrot                      *enforcedRotation;
    BoxDeformation *deform = integrator->deform;
    IMDOutputProvider *outputProvider = integrator->outputProvider;
    t_inputrec *ir = integrator->inputrec;
    const gmx_mtop_t *top_global = integrator->top_global;
    const t_fcdata *fcd = integrator->fcd;
    t_state *state_global = integrator->state_global;
    ObservablesHistory *observablesHistory = integrator->observablesHistory;
    MDAtoms *mdAtoms = integrator->mdAtoms;
    t_nrnb *nrnb = integrator->nrnb;
    gmx_wallcycle *wcycle = integrator->wcycle;
    t_forcerec *fr = integrator->fr;
    const ReplicaExchangeParameters &replExParams = integrator->replExParams;
    // gmx_membed_t                    *membed;
    gmx_walltime_accounting *walltime_accounting = integrator->walltime_accounting;

    // Data from setup
    bool &bSimAnn = setup->bSimAnn;
    bool &startingFromCheckpoint = setup->startingFromCheckpoint;
    bool &bPMETune = setup->bPMETune;
    bool &bPMETunePrinting = setup->bPMETunePrinting;
    bool &bRerunMD = setup->bRerunMD;
    bool &bStopCM = setup->bStopCM;
    bool &bSumEkinhOld = setup->bSumEkinhOld;
    bool &bTrotter = setup->bTrotter;
    bool &bLastStep = setup->bLastStep;

    if (opt2bSet("-ei", nfile, fnm) || observablesHistory->edsamHistory != nullptr) {
        /* Initialize essential dynamics sampling */
        ed = init_edsam(opt2fn_null("-ei", nfile, fnm), opt2fn("-eo", nfile, fnm),
                        top_global,
                        ir, cr, constr,
                        state_global, observablesHistory,
                        oenv, mdrunOptions.continuationOptions.appendFiles);
    }

    /* Initial values */
    init_md(fplog, cr, outputProvider, ir, oenv, mdrunOptions,
            t.write(), t0.write(), state_global, lam0,
            nrnb, top_global, upd.write(), deform,
            nfile, fnm, outf.write(), mdebin.write(),
            force_vir, shake_vir, total_vir, pres,
            mu_tot, &bSimAnn, vcm.write(), wcycle);

    /* Energy terms and groups */
    init_enerdata(top_global->groups.grps[egcENER].nr, ir->fepvals->n_lambda,
                  enerd.write().get());

    /* Kinetic energy data */
    init_ekindata(fplog, top_global, &(ir->opts), ekind.write().get());
    /* Copy the cos acceleration to the groups struct */
    ekind.write()->cosacc.cos_accel = ir->cos_accel;

    gstat = global_stat_init(ir);

    /* Check for polarizable models and flexible constraints */
    shellfc = init_shell_flexcon(fplog,
                                 top_global, constr ? constr->numFlexibleConstraints() : 0,
                                 ir->nstcalcenergy, DOMAINDECOMP(cr));

    /* Set up interactive MD (IMD) */
    init_IMD(ir, cr, ms, top_global, fplog, ir->nstcalcenergy,
             MASTER(cr) ? as_rvec_array(state_global->x.data()) : nullptr,
             nfile, fnm, oenv, mdrunOptions);

    if (DOMAINDECOMP(cr)) {
        top = dd_init_local_top(top_global);

        stateInstance = compat::make_unique<t_state>();
        state = stateInstance.get();
        dd_init_local_state(cr->dd, state_global, state.write().get());

        /* Distribute the charge groups over the nodes from the master node */
        dd_partition_system(fplog, ir->init_step, cr, TRUE, 1,
                            state_global, top_global, ir,
                            state.write().get(), &f, mdAtoms, top.write().get(), fr,
                            vsite, constr,
                            nrnb, nullptr, FALSE);
        shouldCheckNumberOfBondedInteractions = true;
        update_realloc(upd.write().get(), state->natoms);
    } else {
        state_change_natoms(state_global, state_global->natoms);
        /* We need to allocate one element extra, since we might use
         * (unaligned) 4-wide SIMD loads to access rvec entries.
         */
        f.resize(gmx::paddedRVecVectorSize(state_global->natoms));
        /* Copy the pointer to the global state */
        state = state_global;

        // TODO - get rid of this hack (even uglier than the rest...)
        t_graph *tempGraph = graph.write().get();
        mdAlgorithmsSetupAtomData(cr, ir, top_global, top.write().get(), fr,
                                  &tempGraph, mdAtoms, constr, vsite, shellfc.write().get());
        graph = tempGraph;

        update_realloc(upd.write().get(), state->natoms);
    }

    mdatoms = mdAtoms->mdatoms();

    // NOTE: The global state is no longer used at this point.
    // But state_global is still used as temporary storage space for writing
    // the global state to file and potentially for replica exchange.
    // (Global topology should persist.)

    update_mdatoms(mdatoms.write().get(), state->lambda[efptMASS]);

    const ContinuationOptions &continuationOptions = mdrunOptions.continuationOptions;

    if (ir->bExpanded) {
        init_expanded_ensemble(startingFromCheckpoint, ir, state->dfhist);
    }

    if (MASTER(cr)) {
        if (startingFromCheckpoint) {
            /* Update mdebin with energy history if appending to output files */
            if (continuationOptions.appendFiles) {
                restore_energyhistory_from_state(mdebin.write().get(), observablesHistory->energyHistory.get());
            } else if (observablesHistory->energyHistory != nullptr) {
                /* We might have read an energy history from checkpoint.
                 * As we are not appending, we want to restart the statistics.
                 * Free the allocated memory and reset the counts.
                 */
                observablesHistory->energyHistory = {};
            }
        }
        if (observablesHistory->energyHistory == nullptr) {
            observablesHistory->energyHistory = compat::make_unique<energyhistory_t>();
        }
        /* Set the initial energy history in state by updating once */
        update_energyhistory(observablesHistory->energyHistory.get(), mdebin.get());
    }

    // TODO: Remove this by converting AWH into a ForceProvider
    awh = prepareAwhModule(fplog, *ir, state_global, cr, ms, startingFromCheckpoint,
                           shellfc != nullptr,
                           opt2fn("-awh", nfile, fnm), ir->pull_work);

    const bool useReplicaExchange = (replExParams.exchangeInterval > 0);
    if (useReplicaExchange && MASTER(cr)) {
        repl_ex = init_replica_exchange(fplog, ms, top_global->natoms, ir,
                                        replExParams);
    }

    if (bPMETune) {
        pme_loadbal_init(&pme_loadbal, cr, mdlog, ir, state->box,
                         fr->ic, fr->nbv->listParams.get(), fr->pmedata, use_GPU(fr->nbv),
                         &bPMETunePrinting);
    }

    if (!ir->bContinuation && !bRerunMD) {
        if (state->flags & (1 << estV)) {
            /* Set the velocities of vsites, shells and frozen atoms to zero */
            for (int i = 0; i < mdatoms->homenr; i++) {
                if (mdatoms->ptype[i] == eptVSite ||
                    mdatoms->ptype[i] == eptShell) {
                    clear_rvec(state.write()->v[i]);
                } else if (mdatoms->cFREEZE) {
                    for (int m = 0; m < DIM; m++) {
                        if (ir->opts.nFreeze[mdatoms->cFREEZE[i]][m]) {
                            state.write()->v[i][m] = 0;
                        }
                    }
                }
            }
        }

        if (constr) {
            /* Constrain the initial coordinates and velocities */
            do_constrain_first(fplog, constr, ir, mdatoms.write().get(), state.write().get());
        }
        if (vsite) {
            /* Construct the virtual sites for the initial configuration */
            construct_vsites(vsite, as_rvec_array(state.write()->x.data()), ir->delta_t, nullptr,
                             top->idef.iparams, top->idef.il,
                             fr->ePBC, fr->bMolPBC, cr, state->box);
        }
    }

    if (continuationOptions.haveReadEkin) {
        restore_ekinstate_from_state(cr, ekind.write().get(), &state_global->ekinstate);
    }

    int cglo_flags = (CGLO_INITIALIZATION | CGLO_TEMPERATURE | CGLO_GSTAT
                      | (EI_VV(ir->eI) ? CGLO_PRESSURE : 0)
                      | (EI_VV(ir->eI) ? CGLO_CONSTRAINT : 0)
                      | (continuationOptions.haveReadEkin ? CGLO_READEKIN : 0));

    /* To minimize communication, compute_globals computes the COM velocity
     * and the kinetic energy for the velocities without COM motion removed.
     * Thus to get the kinetic energy without the COM contribution, we need
     * to call compute_globals twice.
     */
    SimulationSignaller nullSignaller(nullptr, nullptr, nullptr, false, false);
    for (int cgloIteration = 0; cgloIteration < (bStopCM ? 2 : 1); cgloIteration++) {
        int cglo_flags_iteration = cglo_flags;
        if (bStopCM && cgloIteration == 0) {
            cglo_flags_iteration |= CGLO_STOPCM;
            cglo_flags_iteration &= ~CGLO_TEMPERATURE;
        }
        compute_globals(fplog, gstat, cr, ir, fr, ekind, state, mdatoms, nrnb, vcm,
                        nullptr, enerd, force_vir, shake_vir, total_vir, pres, mu_tot,
                        constr, &nullSignaller, state->box,
                        &totalNumberOfBondedInteractions, &bSumEkinhOld, cglo_flags_iteration
                                                                         | (shouldCheckNumberOfBondedInteractions
                                                                            ? CGLO_CHECK_NUMBER_OF_BONDED_INTERACTIONS
                                                                            : 0));
    }
    checkNumberOfBondedInteractions(fplog, cr, totalNumberOfBondedInteractions,
                                    top_global, top, state,
                                    &shouldCheckNumberOfBondedInteractions);
    if (ir->eI == eiVVAK) {
        /* a second call to get the half step temperature initialized as well */
        /* we do the same call as above, but turn the pressure off -- internally to
           compute_globals, this is recognized as a velocity verlet half-step
           kinetic energy calculation.  This minimized excess variables, but
           perhaps loses some logic?*/

        compute_globals(fplog, gstat, cr, ir, fr, ekind, state, mdatoms, nrnb, vcm,
                        nullptr, enerd, force_vir, shake_vir, total_vir, pres, mu_tot,
                        constr, &nullSignaller, state->box,
                        nullptr, &bSumEkinhOld,
                        cglo_flags & ~CGLO_PRESSURE);
    }

    /* Calculate the initial half step temperature, and save the ekinh_old */
    if (!continuationOptions.startedFromCheckpoint) {
        for (int i = 0; (i < ir->opts.ngtc); i++) {
            copy_mat(ekind->tcstat[i].ekinh, ekind->tcstat[i].ekinh_old);
        }
    }

    /* need to make an initiation call to get the Trotter variables set, as well as other constants for non-trotter
       temperature control */
    trotter_seq = init_npt_vars(ir, state.write().get(), &MassQ, bTrotter);

    if (MASTER(cr)) {
        if (!ir->bContinuation) {
            if (constr && ir->eConstrAlg == econtLINCS) {
                fprintf(fplog,
                        "RMS relative constraint deviation after constraining: %.2e\n",
                        constr->rmsd());
            }
            if (EI_STATE_VELOCITY(ir->eI)) {
                real temp = enerd->term[F_TEMP];
                if (ir->eI != eiVV) {
                    /* Result of Ekin averaged over velocities of -half
                     * and +half step, while we only have -half step here.
                     */
                    temp *= 2;
                }
                fprintf(fplog, "Initial temperature: %g K\n", temp);
            }
        }

        if (bRerunMD) {
            fprintf(stderr, "starting md rerun '%s', reading coordinates from"
                            " input trajectory '%s'\n\n",
                    *(top_global->name), opt2fn("-rerun", nfile, fnm));
            if (mdrunOptions.verbose) {
                fprintf(stderr, "Calculated time to finish depends on nsteps from "
                                "run input file,\nwhich may not correspond to the time "
                                "needed to process input trajectory.\n\n");
            }
        } else {
            char tbuf[20];
            char sbuf[STEPSTRSIZE], sbuf2[STEPSTRSIZE];
            fprintf(stderr, "starting mdrun '%s'\n",
                    *(top_global->name));
            if (ir->nsteps >= 0) {
                sprintf(tbuf, "%8.1f", (ir->init_step + ir->nsteps) * ir->delta_t);
            } else {
                sprintf(tbuf, "%s", "infinite");
            }
            if (ir->init_step > 0) {
                fprintf(stderr, "%s steps, %s ps (continuing from step %s, %8.1f ps).\n",
                        gmx_step_str(ir->init_step + ir->nsteps, sbuf), tbuf,
                        gmx_step_str(ir->init_step, sbuf2),
                        ir->init_step * ir->delta_t);
            } else {
                fprintf(stderr, "%s steps, %s ps.\n",
                        gmx_step_str(ir->nsteps, sbuf), tbuf);
            }
        }
        fprintf(fplog, "\n");
    }

    walltime_accounting_start_time(walltime_accounting);
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

    /***********************************************************
     *
     *             Loop over MD steps
     *
     ************************************************************/

    /* if rerunMD then read coordinates and velocities from input trajectory */
    if (bRerunMD) {
        rerun_fr.natoms = 0;
        if (MASTER(cr)) {
            bLastStep = !read_first_frame(oenv, &status,
                                          opt2fn("-rerun", nfile, fnm),
                                          &rerun_fr, TRX_NEED_X | TRX_READ_V);
            if (rerun_fr.natoms != top_global->natoms) {
                gmx_fatal(FARGS,
                          "Number of atoms in trajectory (%d) does not match the "
                          "run input file (%d)\n",
                          rerun_fr.natoms, top_global->natoms);
            }
            if (ir->ePBC != epbcNONE) {
                if (!rerun_fr.bBox) {
                    gmx_fatal(FARGS,
                              "Rerun trajectory frame step %ld time %f does not contain a box, while pbc is used",
                              rerun_fr.step, rerun_fr.time);
                }
                if (max_cutoff2(ir->ePBC, rerun_fr.box) < gmx::square(fr->rlist)) {
                    gmx_fatal(FARGS, "Rerun trajectory frame step %ld time %f has too small box dimensions",
                              rerun_fr.step, rerun_fr.time);
                }
            }
        }

        if (PAR(cr)) {
            rerun_parallel_comm(cr, &rerun_fr, &bLastStep);
        }

        if (ir->ePBC != epbcNONE) {
            /* Set the shift vectors.
             * Necessary here when have a static box different from the tpr box.
             */
            calc_shifts(rerun_fr.box, fr->shift_vec);
        }
    }

    DdOpenBalanceRegionBeforeForceComputation ddOpenBalanceRegion = (DOMAINDECOMP(cr)
                                                                     ? DdOpenBalanceRegionBeforeForceComputation::yes
                                                                     : DdOpenBalanceRegionBeforeForceComputation::no);
    DdCloseBalanceRegionAfterForceComputation ddCloseBalanceRegion = (DOMAINDECOMP(cr)
                                                                      ? DdCloseBalanceRegionAfterForceComputation::yes
                                                                      : DdCloseBalanceRegionAfterForceComputation::no);

    // TODO extract this to new multi-simulation module
    if (MASTER(cr) && isMultiSim(ms) && !useReplicaExchange) {
        if (!multisim_int_all_are_equal(ms, ir->nsteps)) {
            GMX_LOG(mdlog.warning).appendText(
                        "Note: The number of steps is not consistent across multi simulations,\n"
                        "but we are proceeding anyway!");
        }
        if (!multisim_int_all_are_equal(ms, ir->init_step)) {
            GMX_LOG(mdlog.warning).appendText(
                        "Note: The initial step is not consistent across multi simulations,\n"
                        "but we are proceeding anyway!");
        }
    }

}
