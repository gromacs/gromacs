/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011,2012,2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
 *
 * \brief Implements the integrator for normal molecular dynamics simulations
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_mdrun
 */
#include "gmxpre.h"

#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <memory>

#include "gromacs/awh/awh.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/compat/make_unique.h"
#include "gromacs/domdec/collect.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_network.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/partition.h"
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
#include "gromacs/mdlib/checkpointhandler.h"
#include "gromacs/mdlib/compute_io.h"
#include "gromacs/mdlib/constr.h"
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
#include "gromacs/mdlib/resethandler.h"
#include "gromacs/mdlib/shellfc.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdlib/simulationsignal.h"
#include "gromacs/mdlib/stophandler.h"
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
#include "gromacs/mdtypes/pullhistory.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/output.h"
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

#include "integrator.h"
#include "replicaexchange.h"

#if GMX_FAHCORE
#include "corewrap.h"
#endif

using gmx::SimulationSignaller;

void gmx::Integrator::do_md()
{
    // TODO Historically, the EM and MD "integrators" used different
    // names for the t_inputrec *parameter, but these must have the
    // same name, now that it's a member of a struct. We use this ir
    // alias to avoid a large ripple of nearly useless changes.
    // t_inputrec is being replaced by IMdpOptionsProvider, so this
    // will go away eventually.
    t_inputrec             *ir   = inputrec;
    gmx_mdoutf             *outf = nullptr;
    int64_t                 step, step_rel;
    double                  t, t0, lam0[efptNR];
    gmx_bool                bGStatEveryStep, bGStat, bCalcVir, bCalcEnerStep, bCalcEner;
    gmx_bool                bNS, bNStList, bSimAnn, bStopCM,
                            bFirstStep, bInitStep, bLastStep = FALSE;
    gmx_bool                bDoDHDL = FALSE, bDoFEP = FALSE, bDoExpanded = FALSE;
    gmx_bool                do_ene, do_log, do_verbose;
    gmx_bool                bMasterState;
    int                     force_flags, cglo_flags;
    tensor                  force_vir, shake_vir, total_vir, tmp_vir, pres;
    int                     i, m;
    rvec                    mu_tot;
    t_vcm                  *vcm;
    matrix                  parrinellorahmanMu, M;
    gmx_repl_ex_t           repl_ex = nullptr;
    gmx_localtop_t         *top;
    t_mdebin               *mdebin   = nullptr;
    gmx_enerdata_t         *enerd;
    PaddedVector<gmx::RVec> f {};
    gmx_global_stat_t       gstat;
    gmx_update_t           *upd   = nullptr;
    t_graph                *graph = nullptr;
    gmx_groups_t           *groups;
    gmx_ekindata_t         *ekind;
    gmx_shellfc_t          *shellfc;
    gmx_bool                bSumEkinhOld, bDoReplEx, bExchanged, bNeedRepartition;
    gmx_bool                bTemp, bPres, bTrotter;
    real                    dvdl_constr;
    rvec                   *cbuf        = nullptr;
    int                     cbuf_nalloc = 0;
    matrix                  lastbox;
    int                     lamnew  = 0;
    /* for FEP */
    int                     nstfep = 0;
    double                  cycles;
    real                    saved_conserved_quantity = 0;
    real                    last_ekin                = 0;
    t_extmass               MassQ;
    int                   **trotter_seq;
    char                    sbuf[STEPSTRSIZE], sbuf2[STEPSTRSIZE];

    /* PME load balancing data for GPU kernels */
    pme_load_balancing_t *pme_loadbal      = nullptr;
    gmx_bool              bPMETune         = FALSE;
    gmx_bool              bPMETunePrinting = FALSE;

    /* Interactive MD */
    gmx_bool          bIMDstep = FALSE;

    /* Domain decomposition could incorrectly miss a bonded
       interaction, but checking for that requires a global
       communication stage, which does not otherwise happen in DD
       code. So we do that alongside the first global energy reduction
       after a new DD is made. These variables handle whether the
       check happens, and the result it returns. */
    bool              shouldCheckNumberOfBondedInteractions = false;
    int               totalNumberOfBondedInteractions       = -1;

    SimulationSignals signals;
    // Most global communnication stages don't propagate mdrun
    // signals, and will use this object to achieve that.
    SimulationSignaller nullSignaller(nullptr, nullptr, nullptr, false, false);

    if (!mdrunOptions.writeConfout)
    {
        // This is on by default, and the main known use case for
        // turning it off is for convenience in benchmarking, which is
        // something that should not show up in the general user
        // interface.
        GMX_LOG(mdlog.info).asParagraph().
            appendText("The -noconfout functionality is deprecated, and may be removed in a future version.");
    }

    /* md-vv uses averaged full step velocities for T-control
       md-vv-avek uses averaged half step velocities for T-control (but full step ekin for P control)
       md uses averaged half step kinetic energies to determine temperature unless defined otherwise by GMX_EKIN_AVE_VEL; */
    bTrotter = (EI_VV(ir->eI) && (inputrecNptTrotter(ir) || inputrecNphTrotter(ir) || inputrecNvtTrotter(ir)));

    const bool bRerunMD      = false;
    int        nstglobalcomm = mdrunOptions.globalCommunicationInterval;

    nstglobalcomm   = check_nstglobalcomm(mdlog, nstglobalcomm, ir, cr);
    bGStatEveryStep = (nstglobalcomm == 1);

    groups = &top_global->groups;

    std::unique_ptr<EssentialDynamics> ed = nullptr;
    if (opt2bSet("-ei", nfile, fnm) || observablesHistory->edsamHistory != nullptr)
    {
        /* Initialize essential dynamics sampling */
        ed = init_edsam(mdlog,
                        opt2fn_null("-ei", nfile, fnm), opt2fn("-eo", nfile, fnm),
                        top_global,
                        ir, cr, constr,
                        state_global, observablesHistory,
                        oenv, mdrunOptions.continuationOptions.appendFiles);
    }

    /* Initial values */
    init_md(fplog, cr, outputProvider, ir, oenv, mdrunOptions,
            &t, &t0, state_global, lam0,
            nrnb, top_global, &upd, deform,
            nfile, fnm, &outf, &mdebin,
            force_vir, shake_vir, total_vir, pres, mu_tot, &bSimAnn, &vcm, wcycle);

    /* Energy terms and groups */
    snew(enerd, 1);
    init_enerdata(top_global->groups.grps[egcENER].nr, ir->fepvals->n_lambda,
                  enerd);

    /* Kinetic energy data */
    snew(ekind, 1);
    init_ekindata(fplog, top_global, &(ir->opts), ekind);
    /* Copy the cos acceleration to the groups struct */
    ekind->cosacc.cos_accel = ir->cos_accel;

    gstat = global_stat_init(ir);

    /* Check for polarizable models and flexible constraints */
    shellfc = init_shell_flexcon(fplog,
                                 top_global, constr ? constr->numFlexibleConstraints() : 0,
                                 ir->nstcalcenergy, DOMAINDECOMP(cr));

    {
        double io = compute_io(ir, top_global->natoms, groups, mdebin->ebin->nener, 1);
        if ((io > 2000) && MASTER(cr))
        {
            fprintf(stderr,
                    "\nWARNING: This run will generate roughly %.0f Mb of data\n\n",
                    io);
        }
    }

    /* Set up interactive MD (IMD) */
    init_IMD(ir, cr, ms, top_global, fplog, ir->nstcalcenergy,
             MASTER(cr) ? state_global->x.rvec_array() : nullptr,
             nfile, fnm, oenv, mdrunOptions);

    // Local state only becomes valid now.
    std::unique_ptr<t_state> stateInstance;
    t_state *                state;

    if (DOMAINDECOMP(cr))
    {
        top = dd_init_local_top(top_global);

        stateInstance = compat::make_unique<t_state>();
        state         = stateInstance.get();
        dd_init_local_state(cr->dd, state_global, state);

        /* Distribute the charge groups over the nodes from the master node */
        dd_partition_system(fplog, mdlog, ir->init_step, cr, TRUE, 1,
                            state_global, top_global, ir,
                            state, &f, mdAtoms, top, fr,
                            vsite, constr,
                            nrnb, nullptr, FALSE);
        shouldCheckNumberOfBondedInteractions = true;
        update_realloc(upd, state->natoms);
    }
    else
    {
        state_change_natoms(state_global, state_global->natoms);
        f.resizeWithPadding(state_global->natoms);
        /* Copy the pointer to the global state */
        state = state_global;

        snew(top, 1);
        mdAlgorithmsSetupAtomData(cr, ir, top_global, top, fr,
                                  &graph, mdAtoms, constr, vsite, shellfc);

        update_realloc(upd, state->natoms);
    }

    auto mdatoms = mdAtoms->mdatoms();

    // NOTE: The global state is no longer used at this point.
    // But state_global is still used as temporary storage space for writing
    // the global state to file and potentially for replica exchange.
    // (Global topology should persist.)

    update_mdatoms(mdatoms, state->lambda[efptMASS]);

    const ContinuationOptions &continuationOptions    = mdrunOptions.continuationOptions;
    bool                       startingFromCheckpoint = continuationOptions.startedFromCheckpoint;

    if (ir->bExpanded)
    {
        /* Check nstexpanded here, because the grompp check was broken */
        if (ir->expandedvals->nstexpanded % ir->nstcalcenergy != 0)
        {
            gmx_fatal(FARGS, "With expanded ensemble, nstexpanded should be a multiple of nstcalcenergy");
        }
        init_expanded_ensemble(startingFromCheckpoint, ir, state->dfhist);
    }

    if (MASTER(cr))
    {
        if (startingFromCheckpoint)
        {
            /* Update mdebin with energy history if appending to output files */
            if (continuationOptions.appendFiles)
            {
                /* If no history is available (because a checkpoint is from before
                 * it was written) make a new one later, otherwise restore it.
                 */
                if (observablesHistory->energyHistory)
                {
                    restore_energyhistory_from_state(mdebin, observablesHistory->energyHistory.get());
                }
            }
            else if (observablesHistory->energyHistory)
            {
                /* We might have read an energy history from checkpoint.
                 * As we are not appending, we want to restart the statistics.
                 * Free the allocated memory and reset the counts.
                 */
                observablesHistory->energyHistory = {};
                /* We might have read a pull history from checkpoint.
                 * We will still want to keep the statistics, so that the files
                 * can be joined and still be meaningful.
                 * This means that observablesHistory->pullHistory
                 * should not be reset.
                 */
            }
        }
        if (!observablesHistory->energyHistory)
        {
            observablesHistory->energyHistory = compat::make_unique<energyhistory_t>();
        }
        if (!observablesHistory->pullHistory)
        {
            observablesHistory->pullHistory = compat::make_unique<PullHistory>();
        }
        /* Set the initial energy history in state by updating once */
        update_energyhistory(observablesHistory->energyHistory.get(), mdebin);
    }

    preparePrevStepPullCom(ir, mdatoms, state, state_global, cr, startingFromCheckpoint);

    // TODO: Remove this by converting AWH into a ForceProvider
    auto awh = prepareAwhModule(fplog, *ir, state_global, cr, ms, startingFromCheckpoint,
                                shellfc != nullptr,
                                opt2fn("-awh", nfile, fnm), ir->pull_work);

    const bool useReplicaExchange = (replExParams.exchangeInterval > 0);
    if (useReplicaExchange && MASTER(cr))
    {
        repl_ex = init_replica_exchange(fplog, ms, top_global->natoms, ir,
                                        replExParams);
    }
    /* PME tuning is only supported in the Verlet scheme, with PME for
     * Coulomb. It is not supported with only LJ PME. */
    bPMETune = (mdrunOptions.tunePme && EEL_PME(fr->ic->eeltype) &&
                !mdrunOptions.reproducible && ir->cutoff_scheme != ecutsGROUP);
    if (bPMETune)
    {
        pme_loadbal_init(&pme_loadbal, cr, mdlog, *ir, state->box,
                         *fr->ic, *fr->nbv->listParams, fr->pmedata, use_GPU(fr->nbv),
                         &bPMETunePrinting);
    }

    if (!ir->bContinuation)
    {
        if (state->flags & (1 << estV))
        {
            auto v = makeArrayRef(state->v);
            /* Set the velocities of vsites, shells and frozen atoms to zero */
            for (i = 0; i < mdatoms->homenr; i++)
            {
                if (mdatoms->ptype[i] == eptVSite ||
                    mdatoms->ptype[i] == eptShell)
                {
                    clear_rvec(v[i]);
                }
                else if (mdatoms->cFREEZE)
                {
                    for (m = 0; m < DIM; m++)
                    {
                        if (ir->opts.nFreeze[mdatoms->cFREEZE[i]][m])
                        {
                            v[i][m] = 0;
                        }
                    }
                }
            }
        }

        if (constr)
        {
            /* Constrain the initial coordinates and velocities */
            do_constrain_first(fplog, constr, ir, mdatoms, state);
        }
        if (vsite)
        {
            /* Construct the virtual sites for the initial configuration */
            construct_vsites(vsite, state->x.rvec_array(), ir->delta_t, nullptr,
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

    if (continuationOptions.haveReadEkin)
    {
        restore_ekinstate_from_state(cr, ekind, &state_global->ekinstate);
    }

    cglo_flags = (CGLO_INITIALIZATION | CGLO_TEMPERATURE | CGLO_GSTAT
                  | (EI_VV(ir->eI) ? CGLO_PRESSURE : 0)
                  | (EI_VV(ir->eI) ? CGLO_CONSTRAINT : 0)
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
        compute_globals(fplog, gstat, cr, ir, fr, ekind, state, mdatoms, nrnb, vcm,
                        nullptr, enerd, force_vir, shake_vir, total_vir, pres, mu_tot,
                        constr, &nullSignaller, state->box,
                        &totalNumberOfBondedInteractions, &bSumEkinhOld, cglo_flags_iteration
                        | (shouldCheckNumberOfBondedInteractions ? CGLO_CHECK_NUMBER_OF_BONDED_INTERACTIONS : 0));
    }
    checkNumberOfBondedInteractions(mdlog, cr, totalNumberOfBondedInteractions,
                                    top_global, top, state,
                                    &shouldCheckNumberOfBondedInteractions);
    if (ir->eI == eiVVAK)
    {
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
    if (!continuationOptions.startedFromCheckpoint)
    {
        for (i = 0; (i < ir->opts.ngtc); i++)
        {
            copy_mat(ekind->tcstat[i].ekinh, ekind->tcstat[i].ekinh_old);
        }
    }

    /* need to make an initiation call to get the Trotter variables set, as well as other constants for non-trotter
       temperature control */
    trotter_seq = init_npt_vars(ir, state, &MassQ, bTrotter);

    if (MASTER(cr))
    {
        if (!ir->bContinuation)
        {
            if (constr && ir->eConstrAlg == econtLINCS)
            {
                fprintf(fplog,
                        "RMS relative constraint deviation after constraining: %.2e\n",
                        constr->rmsd());
            }
            if (EI_STATE_VELOCITY(ir->eI))
            {
                real temp = enerd->term[F_TEMP];
                if (ir->eI != eiVV)
                {
                    /* Result of Ekin averaged over velocities of -half
                     * and +half step, while we only have -half step here.
                     */
                    temp *= 2;
                }
                fprintf(fplog, "Initial temperature: %g K\n", temp);
            }
        }

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
        fprintf(fplog, "\n");
    }

    walltime_accounting_start_time(walltime_accounting);
    wallcycle_start(wcycle, ewcRUN);
    print_start(fplog, cr, walltime_accounting, "mdrun");

#if GMX_FAHCORE
    /* safest point to do file checkpointing is here.  More general point would be immediately before integrator call */
    int chkpt_ret = fcCheckPointParallel( cr->nodeid,
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

    bFirstStep       = TRUE;
    /* Skip the first Nose-Hoover integration when we get the state from tpx */
    bInitStep        = !startingFromCheckpoint || EI_VV(ir->eI);
    bSumEkinhOld     = FALSE;
    bExchanged       = FALSE;
    bNeedRepartition = FALSE;

    bool simulationsShareState = false;
    int  nstSignalComm         = nstglobalcomm;
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
        simulationsShareState = useReplicaExchange || usingEnsembleRestraints || awhUsesMultiSim;

        if (simulationsShareState)
        {
            // Inter-simulation signal communication does not need to happen
            // often, so we use a minimum of 200 steps to reduce overhead.
            const int c_minimumInterSimulationSignallingInterval = 200;
            nstSignalComm = ((c_minimumInterSimulationSignallingInterval + nstglobalcomm - 1)/nstglobalcomm)*nstglobalcomm;
        }
    }

    auto stopHandler = stopHandlerBuilder->getStopHandlerMD(
                compat::not_null<SimulationSignal*>(&signals[eglsSTOPCOND]), simulationsShareState,
                MASTER(cr), ir->nstlist, mdrunOptions.reproducible, nstSignalComm,
                mdrunOptions.maximumHoursToRun, ir->nstlist == 0, fplog, step, bNS, walltime_accounting);

    auto checkpointHandler = compat::make_unique<CheckpointHandler>(
                compat::make_not_null<SimulationSignal*>(&signals[eglsCHKPT]),
                simulationsShareState, ir->nstlist == 0, MASTER(cr),
                mdrunOptions.writeConfout, mdrunOptions.checkpointOptions.period);

    const bool resetCountersIsLocal = true;
    auto       resetHandler         = compat::make_unique<ResetHandler>(
                compat::make_not_null<SimulationSignal*>(&signals[eglsRESETCOUNTERS]), !resetCountersIsLocal,
                ir->nsteps, MASTER(cr), mdrunOptions.timingOptions.resetHalfway,
                mdrunOptions.maximumHoursToRun, mdlog, wcycle, walltime_accounting);

    DdOpenBalanceRegionBeforeForceComputation ddOpenBalanceRegion   = (DOMAINDECOMP(cr) ? DdOpenBalanceRegionBeforeForceComputation::yes : DdOpenBalanceRegionBeforeForceComputation::no);
    DdCloseBalanceRegionAfterForceComputation ddCloseBalanceRegion  = (DOMAINDECOMP(cr) ? DdCloseBalanceRegionAfterForceComputation::yes : DdCloseBalanceRegionAfterForceComputation::no);

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

    /* and stop now if we should */
    bLastStep = (bLastStep || (ir->nsteps >= 0 && step_rel > ir->nsteps));
    while (!bLastStep)
    {

        /* Determine if this is a neighbor search step */
        bNStList = (ir->nstlist > 0  && step % ir->nstlist == 0);

        if (bPMETune && bNStList)
        {
            /* PME grid + cut-off optimization with GPUs or PME nodes */
            pme_loadbal_do(pme_loadbal, cr,
                           (mdrunOptions.verbose && MASTER(cr)) ? stderr : nullptr,
                           fplog, mdlog,
                           *ir, fr, *state,
                           wcycle,
                           step, step_rel,
                           &bPMETunePrinting);
        }

        wallcycle_start(wcycle, ewcSTEP);

        bLastStep = (step_rel == ir->nsteps);
        t         = t0 + step*ir->delta_t;

        // TODO Refactor this, so that nstfep does not need a default value of zero
        if (ir->efep != efepNO || ir->bSimTemp)
        {
            /* find and set the current lambdas */
            setCurrentLambdasLocal(step, ir->fepvals, lam0, state);

            bDoDHDL      = do_per_step(step, ir->fepvals->nstdhdl);
            bDoFEP       = ((ir->efep != efepNO) && do_per_step(step, nstfep));
            bDoExpanded  = (do_per_step(step, ir->expandedvals->nstexpanded)
                            && (ir->bExpanded) && (step > 0) && (!startingFromCheckpoint));
        }

        bDoReplEx = (useReplicaExchange && (step > 0) && !bLastStep &&
                     do_per_step(step, replExParams.exchangeInterval));

        if (bSimAnn)
        {
            update_annealing_target_temp(ir, t, upd);
        }

        /* Stop Center of Mass motion */
        bStopCM = (ir->comm_mode != ecmNO && do_per_step(step, ir->nstcomm));

        /* Determine whether or not to do Neighbour Searching */
        bNS = (bFirstStep || bNStList || bExchanged || bNeedRepartition);

        /* Note that the stopHandler will cause termination at nstglobalcomm
         * steps. Since this concides with nstcalcenergy, nsttcouple and/or
         * nstpcouple steps, we have computed the half-step kinetic energy
         * of the previous step and can always output energies at the last step.
         */
        bLastStep = bLastStep || stopHandler->stoppingAfterCurrentStep(bNS);

        /* do_log triggers energy and virial calculation. Because this leads
         * to different code paths, forces can be different. Thus for exact
         * continuation we should avoid extra log output.
         * Note that the || bLastStep can result in non-exact continuation
         * beyond the last step. But we don't consider that to be an issue.
         */
        do_log     = do_per_step(step, ir->nstlog) || (bFirstStep && !startingFromCheckpoint) || bLastStep;
        do_verbose = mdrunOptions.verbose &&
            (step % mdrunOptions.verboseStepPrintInterval == 0 || bFirstStep || bLastStep);

        if (bNS && !(bFirstStep && ir->bContinuation))
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

            if (DOMAINDECOMP(cr))
            {
                /* Repartition the domain decomposition */
                dd_partition_system(fplog, mdlog, step, cr,
                                    bMasterState, nstglobalcomm,
                                    state_global, top_global, ir,
                                    state, &f, mdAtoms, top, fr,
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
            checkNumberOfBondedInteractions(mdlog, cr, totalNumberOfBondedInteractions,
                                            top_global, top, state,
                                            &shouldCheckNumberOfBondedInteractions);
        }
        clear_mat(force_vir);

        checkpointHandler->decideIfCheckpointingThisStep(bNS, bFirstStep, bLastStep);

        /* Determine the energy and pressure:
         * at nstcalcenergy steps and at energy output steps (set below).
         */
        if (EI_VV(ir->eI) && (!bInitStep))
        {
            bCalcEnerStep = do_per_step(step, ir->nstcalcenergy);
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
                       ((inputrecDynamicBox(ir)) ? GMX_FORCE_DYNAMICBOX : 0) |
                       GMX_FORCE_ALLFORCES |
                       (bCalcVir ? GMX_FORCE_VIRIAL : 0) |
                       (bCalcEner ? GMX_FORCE_ENERGY : 0) |
                       (bDoFEP ? GMX_FORCE_DHDL : 0)
                       );

        if (shellfc)
        {
            /* Now is the time to relax the shells */
            relax_shell_flexcon(fplog, cr, ms, mdrunOptions.verbose,
                                enforcedRotation, step,
                                ir, bNS, force_flags, top,
                                constr, enerd, fcd,
                                state, f.arrayRefWithPadding(), force_vir, mdatoms,
                                nrnb, wcycle, graph, groups,
                                shellfc, fr, ppForceWorkload, t, mu_tot,
                                vsite,
                                ddOpenBalanceRegion, ddCloseBalanceRegion);
        }
        else
        {
            /* The AWH history need to be saved _before_ doing force calculations where the AWH bias is updated
               (or the AWH update will be performed twice for one step when continuing). It would be best to
               call this update function from do_md_trajectory_writing but that would occur after do_force.
               One would have to divide the update_awh function into one function applying the AWH force
               and one doing the AWH bias update. The update AWH bias function could then be called after
               do_md_trajectory_writing (then containing update_awh_history).
               The checkpointing will in the future probably moved to the start of the md loop which will
               rid of this issue. */
            if (awh && checkpointHandler->isCheckpointingStep() && MASTER(cr))
            {
                awh->updateHistory(state_global->awhHistory.get());
            }

            /* The coordinates (x) are shifted (to get whole molecules)
             * in do_force.
             * This is parallellized as well, and does communication too.
             * Check comments in sim_util.c
             */
            do_force(fplog, cr, ms, ir, awh.get(), enforcedRotation,
                     step, nrnb, wcycle, top, groups,
                     state->box, state->x.arrayRefWithPadding(), &state->hist,
                     f.arrayRefWithPadding(), force_vir, mdatoms, enerd, fcd,
                     state->lambda, graph,
                     fr, ppForceWorkload, vsite, mu_tot, t, ed ? ed->getLegacyED() : nullptr,
                     (bNS ? GMX_FORCE_NS : 0) | force_flags,
                     ddOpenBalanceRegion, ddCloseBalanceRegion);
        }

        if (EI_VV(ir->eI) && !startingFromCheckpoint)
        /*  ############### START FIRST UPDATE HALF-STEP FOR VV METHODS############### */
        {
            rvec *vbuf = nullptr;

            wallcycle_start(wcycle, ewcUPDATE);
            if (ir->eI == eiVV && bInitStep)
            {
                /* if using velocity verlet with full time step Ekin,
                 * take the first half step only to compute the
                 * virial for the first step. From there,
                 * revert back to the initial coordinates
                 * so that the input is actually the initial step.
                 */
                snew(vbuf, state->natoms);
                copy_rvecn(state->v.rvec_array(), vbuf, 0, state->natoms); /* should make this better for parallelizing? */
            }
            else
            {
                /* this is for NHC in the Ekin(t+dt/2) version of vv */
                trotter_update(ir, step, ekind, enerd, state, total_vir, mdatoms, &MassQ, trotter_seq, ettTSEQ1);
            }

            update_coords(step, ir, mdatoms, state, f.arrayRefWithPadding(), fcd,
                          ekind, M, upd, etrtVELOCITY1,
                          cr, constr);

            wallcycle_stop(wcycle, ewcUPDATE);
            constrain_velocities(step, nullptr,
                                 state,
                                 shake_vir,
                                 constr,
                                 bCalcVir, do_log, do_ene);
            wallcycle_start(wcycle, ewcUPDATE);
            /* if VV, compute the pressure and constraints */
            /* For VV2, we strictly only need this if using pressure
             * control, but we really would like to have accurate pressures
             * printed out.
             * Think about ways around this in the future?
             * For now, keep this choice in comments.
             */
            /*bPres = (ir->eI==eiVV || inputrecNptTrotter(ir)); */
            /*bTemp = ((ir->eI==eiVV &&(!bInitStep)) || (ir->eI==eiVVAK && inputrecNptTrotter(ir)));*/
            bPres = TRUE;
            bTemp = ((ir->eI == eiVV && (!bInitStep)) || (ir->eI == eiVVAK));
            if (bCalcEner && ir->eI == eiVVAK)
            {
                bSumEkinhOld = TRUE;
            }
            /* for vv, the first half of the integration actually corresponds to the previous step.
               So we need information from the last step in the first half of the integration */
            if (bGStat || do_per_step(step-1, nstglobalcomm))
            {
                wallcycle_stop(wcycle, ewcUPDATE);
                compute_globals(fplog, gstat, cr, ir, fr, ekind, state, mdatoms, nrnb, vcm,
                                wcycle, enerd, force_vir, shake_vir, total_vir, pres, mu_tot,
                                constr, &nullSignaller, state->box,
                                &totalNumberOfBondedInteractions, &bSumEkinhOld,
                                (bGStat ? CGLO_GSTAT : 0)
                                | (bCalcEner ? CGLO_ENERGY : 0)
                                | (bTemp ? CGLO_TEMPERATURE : 0)
                                | (bPres ? CGLO_PRESSURE : 0)
                                | (bPres ? CGLO_CONSTRAINT : 0)
                                | (bStopCM ? CGLO_STOPCM : 0)
                                | (shouldCheckNumberOfBondedInteractions ? CGLO_CHECK_NUMBER_OF_BONDED_INTERACTIONS : 0)
                                | CGLO_SCALEEKIN
                                );
                /* explanation of above:
                   a) We compute Ekin at the full time step
                   if 1) we are using the AveVel Ekin, and it's not the
                   initial step, or 2) if we are using AveEkin, but need the full
                   time step kinetic energy for the pressure (always true now, since we want accurate statistics).
                   b) If we are using EkinAveEkin for the kinetic energy for the temperature control, we still feed in
                   EkinAveVel because it's needed for the pressure */
                checkNumberOfBondedInteractions(mdlog, cr, totalNumberOfBondedInteractions,
                                                top_global, top, state,
                                                &shouldCheckNumberOfBondedInteractions);
                wallcycle_start(wcycle, ewcUPDATE);
            }
            /* temperature scaling and pressure scaling to produce the extended variables at t+dt */
            if (!bInitStep)
            {
                if (bTrotter)
                {
                    m_add(force_vir, shake_vir, total_vir);     /* we need the un-dispersion corrected total vir here */
                    trotter_update(ir, step, ekind, enerd, state, total_vir, mdatoms, &MassQ, trotter_seq, ettTSEQ2);

                    /* TODO This is only needed when we're about to write
                     * a checkpoint, because we use it after the restart
                     * (in a kludge?). But what should we be doing if
                     * startingFromCheckpoint or bInitStep are true? */
                    if (inputrecNptTrotter(ir) || inputrecNphTrotter(ir))
                    {
                        copy_mat(shake_vir, state->svir_prev);
                        copy_mat(force_vir, state->fvir_prev);
                    }
                    if (inputrecNvtTrotter(ir) && ir->eI == eiVV)
                    {
                        /* update temperature and kinetic energy now that step is over - this is the v(t+dt) point */
                        enerd->term[F_TEMP] = sum_ekin(&(ir->opts), ekind, nullptr, (ir->eI == eiVV), FALSE);
                        enerd->term[F_EKIN] = trace(ekind->ekin);
                    }
                }
                else if (bExchanged)
                {
                    wallcycle_stop(wcycle, ewcUPDATE);
                    /* We need the kinetic energy at minus the half step for determining
                     * the full step kinetic energy and possibly for T-coupling.*/
                    /* This may not be quite working correctly yet . . . . */
                    compute_globals(fplog, gstat, cr, ir, fr, ekind, state, mdatoms, nrnb, vcm,
                                    wcycle, enerd, nullptr, nullptr, nullptr, nullptr, mu_tot,
                                    constr, &nullSignaller, state->box,
                                    nullptr, &bSumEkinhOld,
                                    CGLO_GSTAT | CGLO_TEMPERATURE);
                    wallcycle_start(wcycle, ewcUPDATE);
                }
            }
            /* if it's the initial step, we performed this first step just to get the constraint virial */
            if (ir->eI == eiVV && bInitStep)
            {
                copy_rvecn(vbuf, state->v.rvec_array(), 0, state->natoms);
                sfree(vbuf);
            }
            wallcycle_stop(wcycle, ewcUPDATE);
        }

        /* compute the conserved quantity */
        if (EI_VV(ir->eI))
        {
            saved_conserved_quantity = NPT_energy(ir, state, &MassQ);
            if (ir->eI == eiVV)
            {
                last_ekin = enerd->term[F_EKIN];
            }
            if ((ir->eDispCorr != edispcEnerPres) && (ir->eDispCorr != edispcAllEnerPres))
            {
                saved_conserved_quantity -= enerd->term[F_DISPCORR];
            }
            /* sum up the foreign energy and dhdl terms for vv.  currently done every step so that dhdl is correct in the .edr */
            if (ir->efep != efepNO)
            {
                sum_dhdl(enerd, state->lambda, ir->fepvals);
            }
        }

        /* ########  END FIRST UPDATE STEP  ############## */
        /* ########  If doing VV, we now have v(dt) ###### */
        if (bDoExpanded)
        {
            /* perform extended ensemble sampling in lambda - we don't
               actually move to the new state before outputting
               statistics, but if performing simulated tempering, we
               do update the velocities and the tau_t. */

            lamnew = ExpandedEnsembleDynamics(fplog, ir, enerd, state, &MassQ, state->fep_state, state->dfhist, step, state->v.rvec_array(), mdatoms);
            /* history is maintained in state->dfhist, but state_global is what is sent to trajectory and log output */
            if (MASTER(cr))
            {
                copy_df_history(state_global->dfhist, state->dfhist);
            }
        }

        /* Now we have the energies and forces corresponding to the
         * coordinates at time t. We must output all of this before
         * the update.
         */
        do_md_trajectory_writing(fplog, cr, nfile, fnm, step, step_rel, t,
                                 ir, state, state_global, observablesHistory,
                                 top_global, fr,
                                 outf, mdebin, ekind, f,
                                 checkpointHandler->isCheckpointingStep(),
                                 bRerunMD, bLastStep,
                                 mdrunOptions.writeConfout,
                                 bSumEkinhOld);
        /* Check if IMD step and do IMD communication, if bIMD is TRUE. */
        bIMDstep = do_IMD(ir->bIMD, step, cr, bNS, state->box, state->x.rvec_array(), ir, t, wcycle);

        /* kludge -- virial is lost with restart for MTTK NPT control. Must reload (saved earlier). */
        if (startingFromCheckpoint && (inputrecNptTrotter(ir) || inputrecNphTrotter(ir)))
        {
            copy_mat(state->svir_prev, shake_vir);
            copy_mat(state->fvir_prev, force_vir);
        }

        stopHandler->setSignal();
        resetHandler->setSignal(walltime_accounting);

        if (bGStat || !PAR(cr))
        {
            /* In parallel we only have to check for checkpointing in steps
             * where we do global communication,
             *  otherwise the other nodes don't know.
             */
            checkpointHandler->setSignal(walltime_accounting);
        }

        /* #########   START SECOND UPDATE STEP ################# */

        /* at the start of step, randomize or scale the velocities ((if vv. Restriction of Andersen controlled
           in preprocessing */

        if (ETC_ANDERSEN(ir->etc)) /* keep this outside of update_tcouple because of the extra info required to pass */
        {
            gmx_bool bIfRandomize;
            bIfRandomize = update_randomize_velocities(ir, step, cr, mdatoms, state->v, upd, constr);
            /* if we have constraints, we have to remove the kinetic energy parallel to the bonds */
            if (constr && bIfRandomize)
            {
                constrain_velocities(step, nullptr,
                                     state,
                                     tmp_vir,
                                     constr,
                                     bCalcVir, do_log, do_ene);
            }
        }
        /* Box is changed in update() when we do pressure coupling,
         * but we should still use the old box for energy corrections and when
         * writing it to the energy file, so it matches the trajectory files for
         * the same timestep above. Make a copy in a separate array.
         */
        copy_mat(state->box, lastbox);

        dvdl_constr = 0;

        wallcycle_start(wcycle, ewcUPDATE);
        /* UPDATE PRESSURE VARIABLES IN TROTTER FORMULATION WITH CONSTRAINTS */
        if (bTrotter)
        {
            trotter_update(ir, step, ekind, enerd, state, total_vir, mdatoms, &MassQ, trotter_seq, ettTSEQ3);
            /* We can only do Berendsen coupling after we have summed
             * the kinetic energy or virial. Since the happens
             * in global_state after update, we should only do it at
             * step % nstlist = 1 with bGStatEveryStep=FALSE.
             */
        }
        else
        {
            update_tcouple(step, ir, state, ekind, &MassQ, mdatoms);
            update_pcouple_before_coordinates(fplog, step, ir, state,
                                              parrinellorahmanMu, M,
                                              bInitStep);
        }

        if (EI_VV(ir->eI))
        {
            /* velocity half-step update */
            update_coords(step, ir, mdatoms, state, f.arrayRefWithPadding(), fcd,
                          ekind, M, upd, etrtVELOCITY2,
                          cr, constr);
        }

        /* Above, initialize just copies ekinh into ekin,
         * it doesn't copy position (for VV),
         * and entire integrator for MD.
         */

        if (ir->eI == eiVVAK)
        {
            /* We probably only need md->homenr, not state->natoms */
            if (state->natoms > cbuf_nalloc)
            {
                cbuf_nalloc = state->natoms;
                srenew(cbuf, cbuf_nalloc);
            }
            copy_rvecn(as_rvec_array(state->x.data()), cbuf, 0, state->natoms);
        }

        update_coords(step, ir, mdatoms, state, f.arrayRefWithPadding(), fcd,
                      ekind, M, upd, etrtPOSITION, cr, constr);
        wallcycle_stop(wcycle, ewcUPDATE);

        constrain_coordinates(step, &dvdl_constr, state,
                              shake_vir,
                              upd, constr,
                              bCalcVir, do_log, do_ene);
        update_sd_second_half(step, &dvdl_constr, ir, mdatoms, state,
                              cr, nrnb, wcycle, upd, constr, do_log, do_ene);
        finish_update(ir, mdatoms,
                      state, graph,
                      nrnb, wcycle, upd, constr);

        if (ir->bPull && ir->pull->bSetPbcRefToPrevStepCOM)
        {
            updatePrevStepPullCom(ir->pull_work, state);
        }

        if (ir->eI == eiVVAK)
        {
            /* erase F_EKIN and F_TEMP here? */
            /* just compute the kinetic energy at the half step to perform a trotter step */
            compute_globals(fplog, gstat, cr, ir, fr, ekind, state, mdatoms, nrnb, vcm,
                            wcycle, enerd, force_vir, shake_vir, total_vir, pres, mu_tot,
                            constr, &nullSignaller, lastbox,
                            nullptr, &bSumEkinhOld,
                            (bGStat ? CGLO_GSTAT : 0) | CGLO_TEMPERATURE
                            );
            wallcycle_start(wcycle, ewcUPDATE);
            trotter_update(ir, step, ekind, enerd, state, total_vir, mdatoms, &MassQ, trotter_seq, ettTSEQ4);
            /* now we know the scaling, we can compute the positions again again */
            copy_rvecn(cbuf, as_rvec_array(state->x.data()), 0, state->natoms);

            update_coords(step, ir, mdatoms, state, f.arrayRefWithPadding(), fcd,
                          ekind, M, upd, etrtPOSITION, cr, constr);
            wallcycle_stop(wcycle, ewcUPDATE);

            /* do we need an extra constraint here? just need to copy out of as_rvec_array(state->v.data()) to upd->xp? */
            /* are the small terms in the shake_vir here due
             * to numerical errors, or are they important
             * physically? I'm thinking they are just errors, but not completely sure.
             * For now, will call without actually constraining, constr=NULL*/
            finish_update(ir, mdatoms,
                          state, graph,
                          nrnb, wcycle, upd, nullptr);
        }
        if (EI_VV(ir->eI))
        {
            /* this factor or 2 correction is necessary
               because half of the constraint force is removed
               in the vv step, so we have to double it.  See
               the Redmine issue #1255.  It is not yet clear
               if the factor of 2 is exact, or just a very
               good approximation, and this will be
               investigated.  The next step is to see if this
               can be done adding a dhdl contribution from the
               rattle step, but this is somewhat more
               complicated with the current code. Will be
               investigated, hopefully for 4.6.3. However,
               this current solution is much better than
               having it completely wrong.
             */
            enerd->term[F_DVDL_CONSTR] += 2*dvdl_constr;
        }
        else
        {
            enerd->term[F_DVDL_CONSTR] += dvdl_constr;
        }

        if (vsite != nullptr)
        {
            wallcycle_start(wcycle, ewcVSITECONSTR);
            if (graph != nullptr)
            {
                shift_self(graph, state->box, state->x.rvec_array());
            }
            construct_vsites(vsite, state->x.rvec_array(), ir->delta_t, state->v.rvec_array(),
                             top->idef.iparams, top->idef.il,
                             fr->ePBC, fr->bMolPBC, cr, state->box);

            if (graph != nullptr)
            {
                unshift_self(graph, state->box, state->x.rvec_array());
            }
            wallcycle_stop(wcycle, ewcVSITECONSTR);
        }

        /* ############## IF NOT VV, Calculate globals HERE  ############ */
        /* With Leap-Frog we can skip compute_globals at
         * non-communication steps, but we need to calculate
         * the kinetic energy one step before communication.
         */
        {
            // Organize to do inter-simulation signalling on steps if
            // and when algorithms require it.
            const bool doInterSimSignal = (simulationsShareState && do_per_step(step, nstSignalComm));

            // With leap-frog we also need to compute the half-step
            // kinetic energy at the step before we need to write
            // the full-step kinetic energy
            const bool needEkinAtNextStep =
                (!EI_VV(ir->eI) && (do_per_step(step + 1, nstglobalcomm) ||
                                    step_rel + 1 == ir->nsteps));

            if (bGStat || needEkinAtNextStep || doInterSimSignal)
            {
                // Since we're already communicating at this step, we
                // can propagate intra-simulation signals. Note that
                // check_nstglobalcomm has the responsibility for
                // choosing the value of nstglobalcomm that is one way
                // bGStat becomes true, so we can't get into a
                // situation where e.g. checkpointing can't be
                // signalled.
                bool                doIntraSimSignal = true;
                SimulationSignaller signaller(&signals, cr, ms, doInterSimSignal, doIntraSimSignal);

                compute_globals(fplog, gstat, cr, ir, fr, ekind, state, mdatoms, nrnb, vcm,
                                wcycle, enerd, force_vir, shake_vir, total_vir, pres, mu_tot,
                                constr, &signaller,
                                lastbox,
                                &totalNumberOfBondedInteractions, &bSumEkinhOld,
                                (bGStat ? CGLO_GSTAT : 0)
                                | (!EI_VV(ir->eI) && bCalcEner ? CGLO_ENERGY : 0)
                                | (!EI_VV(ir->eI) && bStopCM ? CGLO_STOPCM : 0)
                                | (!EI_VV(ir->eI) ? CGLO_TEMPERATURE : 0)
                                | (!EI_VV(ir->eI) ? CGLO_PRESSURE : 0)
                                | CGLO_CONSTRAINT
                                | (shouldCheckNumberOfBondedInteractions ? CGLO_CHECK_NUMBER_OF_BONDED_INTERACTIONS : 0)
                                );
                checkNumberOfBondedInteractions(mdlog, cr, totalNumberOfBondedInteractions,
                                                top_global, top, state,
                                                &shouldCheckNumberOfBondedInteractions);
            }
        }

        /* #############  END CALC EKIN AND PRESSURE ################# */

        /* Note: this is OK, but there are some numerical precision issues with using the convergence of
           the virial that should probably be addressed eventually. state->veta has better properies,
           but what we actually need entering the new cycle is the new shake_vir value. Ideally, we could
           generate the new shake_vir, but test the veta value for convergence.  This will take some thought. */

        if (ir->efep != efepNO && !EI_VV(ir->eI))
        {
            /* Sum up the foreign energy and dhdl terms for md and sd.
               Currently done every step so that dhdl is correct in the .edr */
            sum_dhdl(enerd, state->lambda, ir->fepvals);
        }

        update_pcouple_after_coordinates(fplog, step, ir, mdatoms,
                                         pres, force_vir, shake_vir,
                                         parrinellorahmanMu,
                                         state, nrnb, upd);

        /* ################# END UPDATE STEP 2 ################# */
        /* #### We now have r(t+dt) and v(t+dt/2)  ############# */

        /* The coordinates (x) were unshifted in update */
        if (!bGStat)
        {
            /* We will not sum ekinh_old,
             * so signal that we still have to do it.
             */
            bSumEkinhOld = TRUE;
        }

        if (bCalcEner)
        {
            /* #########  BEGIN PREPARING EDR OUTPUT  ###########  */

            /* use the directly determined last velocity, not actually the averaged half steps */
            if (bTrotter && ir->eI == eiVV)
            {
                enerd->term[F_EKIN] = last_ekin;
            }
            enerd->term[F_ETOT] = enerd->term[F_EPOT] + enerd->term[F_EKIN];

            if (integratorHasConservedEnergyQuantity(ir))
            {
                if (EI_VV(ir->eI))
                {
                    enerd->term[F_ECONSERVED] = enerd->term[F_ETOT] + saved_conserved_quantity;
                }
                else
                {
                    enerd->term[F_ECONSERVED] = enerd->term[F_ETOT] + NPT_energy(ir, state, &MassQ);
                }
            }
            /* #########  END PREPARING EDR OUTPUT  ###########  */
        }

        /* Output stuff */
        if (MASTER(cr))
        {
            if (fplog && do_log && bDoExpanded)
            {
                /* only needed if doing expanded ensemble */
                PrintFreeEnergyInfoToFile(fplog, ir->fepvals, ir->expandedvals, ir->bSimTemp ? ir->simtempvals : nullptr,
                                          state_global->dfhist, state->fep_state, ir->nstlog, step);
            }
            if (bCalcEner)
            {
                upd_mdebin(mdebin, bDoDHDL, bCalcEnerStep,
                           t, mdatoms->tmass, enerd, state,
                           ir->fepvals, ir->expandedvals, lastbox,
                           shake_vir, force_vir, total_vir, pres,
                           ekind, mu_tot, constr);
            }
            else
            {
                upd_mdebin_step(mdebin);
            }

            gmx_bool do_dr  = do_per_step(step, ir->nstdisreout);
            gmx_bool do_or  = do_per_step(step, ir->nstorireout);

            print_ebin(mdoutf_get_fp_ene(outf), do_ene, do_dr, do_or, do_log ? fplog : nullptr,
                       step, t,
                       eprNORMAL, mdebin, fcd, groups, &(ir->opts), awh.get());

            if (ir->bPull)
            {
                pull_print_output(ir->pull_work, step, t);
            }

            if (do_per_step(step, ir->nstlog))
            {
                if (fflush(fplog) != 0)
                {
                    gmx_fatal(FARGS, "Cannot flush logfile - maybe you are out of disk space?");
                }
            }
        }
        if (bDoExpanded)
        {
            /* Have to do this part _after_ outputting the logfile and the edr file */
            /* Gets written into the state at the beginning of next loop*/
            state->fep_state = lamnew;
        }
        /* Print the remaining wall clock time for the run */
        if (isMasterSimMasterRank(ms, cr) &&
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
                                             MASTER(cr) && mdrunOptions.verbose,
                                             bRerunMD);

            if (bNeedRepartition && DOMAINDECOMP(cr))
            {
                dd_collect_state(cr->dd, state, state_global);
            }
        }

        /* Replica exchange */
        bExchanged = FALSE;
        if (bDoReplEx)
        {
            bExchanged = replica_exchange(fplog, cr, ms, repl_ex,
                                          state_global, enerd,
                                          state, step, t);
        }

        if ( (bExchanged || bNeedRepartition) && DOMAINDECOMP(cr) )
        {
            dd_partition_system(fplog, mdlog, step, cr, TRUE, 1,
                                state_global, top_global, ir,
                                state, &f, mdAtoms, top, fr,
                                vsite, constr,
                                nrnb, wcycle, FALSE);
            shouldCheckNumberOfBondedInteractions = true;
            update_realloc(upd, state->natoms);
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

        if ( (membed != nullptr) && (!bLastStep) )
        {
            rescale_membed(step_rel, membed, as_rvec_array(state_global->x.data()));
        }

        cycles = wallcycle_stop(wcycle, ewcSTEP);
        if (DOMAINDECOMP(cr) && wcycle)
        {
            dd_cycles_add(cr->dd, cycles, ddCyclStep);
        }

        /* increase the MD step number */
        step++;
        step_rel++;

        resetHandler->resetCounters(
                step, step_rel, mdlog, fplog, cr, (use_GPU(fr->nbv) ? fr->nbv : nullptr),
                nrnb, fr->pmedata, pme_loadbal, wcycle, walltime_accounting);

        /* If bIMD is TRUE, the master updates the IMD energy record and sends positions to VMD client */
        IMD_prep_energies_send_positions(ir->bIMD && MASTER(cr), bIMDstep, ir->imd, enerd, step, bCalcEner, wcycle);

    }
    /* End of main MD loop */

    /* Closing TNG files can include compressing data. Therefore it is good to do that
     * before stopping the time measurements. */
    mdoutf_tng_close(outf);

    /* Stop measuring walltime */
    walltime_accounting_end_time(walltime_accounting);

    if (!thisRankHasDuty(cr, DUTY_PME))
    {
        /* Tell the PME only node to finish */
        gmx_pme_send_finish(cr);
    }

    if (MASTER(cr))
    {
        if (ir->nstcalcenergy > 0)
        {
            print_ebin(mdoutf_get_fp_ene(outf), FALSE, FALSE, FALSE, fplog, step, t,
                       eprAVER, mdebin, fcd, groups, &(ir->opts), awh.get());
        }
    }
    done_mdebin(mdebin);
    done_mdoutf(outf);

    if (bPMETune)
    {
        pme_loadbal_done(pme_loadbal, fplog, mdlog, use_GPU(fr->nbv));
    }

    done_shellfc(fplog, shellfc, step_rel);

    if (useReplicaExchange && MASTER(cr))
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

    destroy_enerdata(enerd);
    sfree(enerd);
    sfree(top);
}
