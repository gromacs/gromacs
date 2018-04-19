/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011,2012,2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include <cmath>

#include <algorithm>
#include <memory>
#include <gromacs/mdrun/customintegrator/datamanager.h>
#include <gromacs/mdrun/customintegrator/loopelement.h>
#include <gromacs/mdrun/customintegrator/pmeloadbalanceelement.h>
#include <gromacs/mdrun/customintegrator/ddelement.h>
#include <gromacs/mdrun/customintegrator/forceelement.h>
#include <gromacs/mdrun/customintegrator/updatefusedelement.h>
#include <gromacs/mdrun/customintegrator/constrainelement.h>
#include <gromacs/mdrun/customintegrator/computeglobalselement.h>
#include <gromacs/mdrun/customintegrator/writetrajectoryelement.h>
#include <gromacs/mdrun/customintegrator/writeenergyelement.h>
#include <iostream>

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

#include "integrator.h"

#ifdef GMX_FAHCORE
#include "corewrap.h"
#endif

using gmx::SimulationSignaller;

//! Resets all the counters.
static void reset_all_counters(FILE *fplog, const gmx::MDLogger &mdlog, t_commrec *cr,
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

void gmx::Integrator::do_custom()
{
    // TODO Historically, the EM and MD "integrators" used different
    // names for the t_inputrec *parameter, but these must have the
    // same name, now that it's a member of a struct. We use this ir
    // alias to avoid a large ripple of nearly useless changes.
    // t_inputrec is being replaced by IMdpOptionsProvider, so this
    // will go away eventually.
    t_inputrec       *ir   = inputrec;
    gmx_mdoutf       *outf = nullptr;
    gmx_int64_t       step, step_rel;
    double            elapsed_time;
    double            t, t0, lam0[efptNR];
    gmx_bool          bGStatEveryStep, bGStat, bCalcVir, bCalcEnerStep, bCalcEner;
    gmx_bool          bNS, bNStList, bSimAnn, bStopCM,
                      bFirstStep, bInitStep, bLastStep = FALSE;
    gmx_bool          bDoDHDL = FALSE, bDoFEP = FALSE;
    gmx_bool          do_ene, do_log, do_verbose,
                      bCPT;
    gmx_bool          bMasterState;
    int               force_flags, cglo_flags;
    tensor            force_vir, shake_vir, total_vir, pres;
    int               i, m;
    rvec              mu_tot;
    t_vcm            *vcm;
    matrix            parrinellorahmanMu, M;
    int               nchkpt  = 1;
    gmx_localtop_t   *top;
    t_mdebin         *mdebin   = nullptr;
    gmx_enerdata_t   *enerd;
    PaddedRVecVector  f {};
    gmx_global_stat_t gstat;
    gmx_update_t     *upd   = nullptr;
    t_graph          *graph = nullptr;
    gmx_groups_t     *groups;
    gmx_ekindata_t   *ekind;
    gmx_shellfc_t    *shellfc;
    gmx_bool          bSumEkinhOld;
    gmx_bool          bResetCountersHalfMaxH = FALSE;
    real              dvdl_constr;
    matrix            lastbox;
    /* for FEP */
    double            cycles;
    t_extmass         MassQ;
    char              sbuf[STEPSTRSIZE], sbuf2[STEPSTRSIZE];
    int               handled_stop_condition = gmx_stop_cond_none; /* compare to get_stop_condition*/


    /* PME load balancing data for GPU kernels */
    pme_load_balancing_t *pme_loadbal      = nullptr;
    gmx_bool              bPMETune         = FALSE;
    gmx_bool              bPMETunePrinting = FALSE;

#ifdef GMX_FAHCORE
    /* Temporary addition for FAHCORE checkpointing */
    int chkpt_ret;
#endif
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

    //
    // MISSING FEATURES
    //
    gmx_edsam *ed = nullptr;
    if (opt2bSet("-ei", nfile, fnm) || observablesHistory->edsamHistory != nullptr)
    {
        gmx_fatal(FARGS, "Essential dynamics is not implemented for custom integrators.");
    }

    if (ir->eSwapCoords != eswapNO)
    {
        gmx_fatal(FARGS, "Ion swapping is not implemented for custom integrators.");
    }

    /* Check for polarizable models and flexible constraints */
    shellfc = init_shell_flexcon(fplog,
                                 top_global, n_flexible_constraints(constr),
                                 ir->nstcalcenergy, DOMAINDECOMP(cr));
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

    const bool useReplicaExchange = (replExParams.exchangeInterval > 0);
    if (useReplicaExchange)
    {
        gmx_fatal(FARGS, "Replica exchange is not implemented for custom integrators.");
    }
    if (vsite)
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


    if (!mdrunOptions.writeConfout)
    {
        // This is on by default, and the main known use case for
        // turning it off is for convenience in benchmarking, which is
        // something that should not show up in the general user
        // interface.
        GMX_LOG(mdlog.info).asParagraph().
            appendText("The -noconfout functionality is deprecated, and may be removed in a future version.");
    }
    if (mdrunOptions.timingOptions.resetHalfway)
    {
        GMX_LOG(mdlog.info).asParagraph().
            appendText("The -resethway functionality is deprecated, and may be removed in a future version.");
        if (ir->nsteps > 0)
        {
            /* Signal to reset the counters half the simulation steps. */
            wcycle_set_reset_counters(wcycle, ir->nsteps/2);
        }
        /* Signal to reset the counters halfway the simulation time. */
        bResetCountersHalfMaxH = (mdrunOptions.maximumHoursToRun > 0);
    }

    const gmx_bool bRerunMD      = FALSE;
    int            nstglobalcomm = mdrunOptions.globalCommunicationInterval;

    nstglobalcomm   = check_nstglobalcomm(mdlog, nstglobalcomm, ir);
    bGStatEveryStep = (nstglobalcomm == 1);

    groups = &top_global->groups;

    /* Initial values */
    init_md(fplog, cr, outputProvider, ir, oenv, mdrunOptions,
            &t, &t0, state_global, lam0,
            nrnb, top_global, &upd,
            nfile, fnm, &outf, &mdebin,
            force_vir, shake_vir, mu_tot, &bSimAnn, &vcm, wcycle);

    clear_mat(total_vir);
    clear_mat(pres);
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

    // Local state only becomes valid now.
    std::unique_ptr<t_state> stateInstance;
    t_state *                state;

    if (DOMAINDECOMP(cr))
    {
        top = dd_init_local_top(top_global);

        stateInstance = gmx::compat::make_unique<t_state>();
        state         = stateInstance.get();
        dd_init_local_state(cr->dd, state_global, state);
    }
    else
    {
        state_change_natoms(state_global, state_global->natoms);
        /* We need to allocate one element extra, since we might use
         * (unaligned) 4-wide SIMD loads to access rvec entries.
         */
        f.resize(gmx::paddedRVecVectorSize(state_global->natoms));
        /* Copy the pointer to the global state */
        state = state_global;

        snew(top, 1);
        mdAlgorithmsSetupAtomData(cr, ir, top_global, top, fr,
                                  &graph, mdAtoms, vsite, shellfc);

        update_realloc(upd, state->natoms);
    }

    if (DOMAINDECOMP(cr))
    {
        /* Distribute the charge groups over the nodes from the master node */
        dd_partition_system(fplog, ir->init_step, cr, TRUE, 1,
                            state_global, top_global, ir,
                            state, &f, mdAtoms, top, fr,
                            vsite, constr,
                            nrnb, nullptr, FALSE);
        shouldCheckNumberOfBondedInteractions = true;
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


    if (MASTER(cr))
    {
        if (startingFromCheckpoint)
        {
            /* Update mdebin with energy history if appending to output files */
            if (continuationOptions.appendFiles)
            {
                restore_energyhistory_from_state(mdebin, observablesHistory->energyHistory.get());
            }
            else if (observablesHistory->energyHistory.get() != nullptr)
            {
                /* We might have read an energy history from checkpoint.
                 * As we are not appending, we want to restart the statistics.
                 * Free the allocated memory and reset the counts.
                 */
                observablesHistory->energyHistory = {};
            }
        }
        if (observablesHistory->energyHistory.get() == nullptr)
        {
            observablesHistory->energyHistory = std::unique_ptr<energyhistory_t>(new energyhistory_t {});
        }
        /* Set the initial energy history in state by updating once */
        update_energyhistory(observablesHistory->energyHistory.get(), mdebin);
    }

    /* Initialize constraints */
    if (constr && !DOMAINDECOMP(cr))
    {
        set_constraints(constr, top, ir, mdatoms, cr);
    }

    /* PME tuning is only supported in the Verlet scheme, with PME for
     * Coulomb. It is not supported with only LJ PME, or for
     * reruns. */
    bPMETune = (mdrunOptions.tunePme && EEL_PME(fr->ic->eeltype) &&
                !mdrunOptions.reproducible && ir->cutoff_scheme != ecutsGROUP);
    if (bPMETune)
    {
        pme_loadbal_init(&pme_loadbal, cr, mdlog, ir, state->box,
                         fr->ic, fr->nbv->listParams.get(), fr->pmedata, use_GPU(fr->nbv),
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

        if (constr)
        {
            /* Constrain the initial coordinates and velocities */
            do_constrain_first(fplog, constr, ir, mdatoms, state,
                               cr, ms, nrnb, fr, top);
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
    checkNumberOfBondedInteractions(fplog, cr, totalNumberOfBondedInteractions,
                                    top_global, top, state,
                                    &shouldCheckNumberOfBondedInteractions);

    /* Calculate the initial half step temperature, and save the ekinh_old */
    if (!continuationOptions.startedFromCheckpoint)
    {
        for (i = 0; (i < ir->opts.ngtc); i++)
        {
            copy_mat(ekind->tcstat[i].ekinh, ekind->tcstat[i].ekinh_old);
        }
    }

    if (MASTER(cr))
    {
        if (!ir->bContinuation)
        {
            if (constr && ir->eConstrAlg == econtLINCS)
            {
                fprintf(fplog,
                        "RMS relative constraint deviation after constraining: %.2e\n",
                        constr_rmsd(constr));
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
                fprintf(fplog, "Initial temperature: %g K\n", temp);
            }
        }

        // TODO: Clean up indenting
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

    /***********************************************************
     *
     *             Loop over MD steps
     *
     ************************************************************/

    /* Loop over MD steps or if rerunMD to end of input trajectory,
     * or, if max_hours>0, until max_hours is reached.
     */
    real max_hours   = mdrunOptions.maximumHoursToRun;
    bFirstStep       = TRUE;
    /* Skip the first Nose-Hoover integration when we get the state from tpx */
    bInitStep        = !startingFromCheckpoint;
    bSumEkinhOld     = FALSE;

    int  nstSignalComm         = nstglobalcomm;

    DdOpenBalanceRegionBeforeForceComputation ddOpenBalanceRegion   = (DOMAINDECOMP(cr) ? DdOpenBalanceRegionBeforeForceComputation::yes : DdOpenBalanceRegionBeforeForceComputation::no);
    DdCloseBalanceRegionAfterForceComputation ddCloseBalanceRegion  = (DOMAINDECOMP(cr) ? DdCloseBalanceRegionAfterForceComputation::yes : DdCloseBalanceRegionAfterForceComputation::no);

    step     = ir->init_step;
    step_rel = 0;

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
                           ir, fr, state,
                           wcycle,
                           step, step_rel,
                           &bPMETunePrinting);
        }

        wallcycle_start(wcycle, ewcSTEP);

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
        do_verbose = mdrunOptions.verbose &&
            (step % mdrunOptions.verboseStepPrintInterval == 0 || bFirstStep || bLastStep);

        if (bNS && !(bFirstStep && ir->bContinuation))
        {
            // TODO: Clean up indenting
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

        clear_mat(force_vir);

        /* We write a checkpoint at this MD step when:
         * either at an NS step when we signalled through gs,
         * or at the last step (but not when we do not want confout),
         * but never at the first step or with rerun.
         */
        bCPT = (((signals[eglsCHKPT].set && (bNS || ir->nstlist == 0)) ||
                 (bLastStep && mdrunOptions.writeConfout)) &&
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

        force_flags = (GMX_FORCE_STATECHANGED |
                       ((inputrecDynamicBox(ir)) ? GMX_FORCE_DYNAMICBOX : 0) |
                       GMX_FORCE_ALLFORCES |
                       (bCalcVir ? GMX_FORCE_VIRIAL : 0) |
                       (bCalcEner ? GMX_FORCE_ENERGY : 0) |
                       (bDoFEP ? GMX_FORCE_DHDL : 0)
                       );

        {
            /* The coordinates (x) are shifted (to get whole molecules)
             * in do_force.
             * This is parallellized as well, and does communication too.
             * Check comments in sim_util.c
             */
            do_force(fplog, cr, ms, ir, step, nrnb, wcycle, top, groups,
                     state->box, state->x, &state->hist,
                     f, force_vir, mdatoms, enerd, fcd,
                     state->lambda, graph,
                     fr, vsite, mu_tot, t, ed,
                     (bNS ? GMX_FORCE_NS : 0) | force_flags,
                     ddOpenBalanceRegion, ddCloseBalanceRegion);
        }

        /* ########  END FIRST UPDATE STEP  ############## */
        /* ########  If doing VV, we now have v(dt) ###### */

        /* Now we have the energies and forces corresponding to the
         * coordinates at time t. We must output all of this before
         * the update.
         */
        do_md_trajectory_writing(fplog, cr, nfile, fnm, step, step_rel, t,
                                 ir, state, state_global, observablesHistory,
                                 top_global, fr,
                                 outf, mdebin, ekind, f,
                                 &nchkpt,
                                 bCPT, bRerunMD, bLastStep,
                                 mdrunOptions.writeConfout,
                                 bSumEkinhOld);

        elapsed_time = walltime_accounting_get_current_elapsed_time(walltime_accounting);

        /* Check whether everything is still allright */
        if (((int)gmx_get_stop_condition() > handled_stop_condition)
#if GMX_THREAD_MPI
            && MASTER(cr)
#endif
            )
        {
            int nsteps_stop = -1;

            /* this just makes signals[].sig compatible with the hack
               of sending signals around by MPI_Reduce together with
               other floats */
            if ((gmx_get_stop_condition() == gmx_stop_cond_next_ns) ||
                (mdrunOptions.reproducible &&
                 gmx_get_stop_condition() == gmx_stop_cond_next))
            {
                /* We need at least two global communication steps to pass
                 * around the signal. We stop at a pair-list creation step
                 * to allow for exact continuation, when possible.
                 */
                signals[eglsSTOPCOND].sig = 1;
                nsteps_stop               = std::max(ir->nstlist, 2*nstSignalComm);
            }
            else if (gmx_get_stop_condition() == gmx_stop_cond_next)
            {
                /* Stop directly after the next global communication step.
                 * This breaks exact continuation.
                 */
                signals[eglsSTOPCOND].sig = -1;
                nsteps_stop               = nstSignalComm + 1;
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
            handled_stop_condition = (int)gmx_get_stop_condition();
        }
        else if (MASTER(cr) && (bNS || ir->nstlist <= 0) &&
                 (max_hours > 0 && elapsed_time > max_hours*60.0*60.0*0.99) &&
                 signals[eglsSTOPCOND].sig == 0 && signals[eglsSTOPCOND].set == 0)
        {
            /* Signal to terminate the run */
            signals[eglsSTOPCOND].sig = 1;
            if (fplog)
            {
                fprintf(fplog, "\nStep %s: Run time exceeded %.3f hours, will terminate the run\n", gmx_step_str(step, sbuf), max_hours*0.99);
            }
            fprintf(stderr, "\nStep %s: Run time exceeded %.3f hours, will terminate the run\n", gmx_step_str(step, sbuf), max_hours*0.99);
        }

        if (bResetCountersHalfMaxH && MASTER(cr) &&
            elapsed_time > max_hours*60.0*60.0*0.495)
        {
            /* Set flag that will communicate the signal to all ranks in the simulation */
            signals[eglsRESETCOUNTERS].sig = 1;
        }

        /* In parallel we only have to check for checkpointing in steps
         * where we do global communication,
         *  otherwise the other nodes don't know.
         */
        const real cpt_period = mdrunOptions.checkpointOptions.period;
        if (MASTER(cr) && ((bGStat || !PAR(cr)) &&
                           cpt_period >= 0 &&
                           (cpt_period == 0 ||
                            elapsed_time >= nchkpt*cpt_period*60.0)) &&
            signals[eglsCHKPT].set == 0)
        {
            signals[eglsCHKPT].sig = 1;
        }

        /* #########   START SECOND UPDATE STEP ################# */

        /* Box is changed in update() when we do pressure coupling,
         * but we should still use the old box for energy corrections and when
         * writing it to the energy file, so it matches the trajectory files for
         * the same timestep above. Make a copy in a separate array.
         */
        copy_mat(state->box, lastbox);

        dvdl_constr = 0;

        // TODO: Clean up indenting
        {
            wallcycle_start(wcycle, ewcUPDATE);
            {
                update_tcouple(step, ir, state, ekind, &MassQ, mdatoms);
                update_pcouple_before_coordinates(fplog, step, ir, state,
                                                  parrinellorahmanMu, M,
                                                  bInitStep);
            }

            /* Above, initialize just copies ekinh into ekin,
             * it doesn't copy position (for VV),
             * and entire integrator for MD.
             */

            update_coords(step, ir, mdatoms, state, f, fcd,
                          ekind, M, upd, etrtPOSITION, cr, constr);
            wallcycle_stop(wcycle, ewcUPDATE);

            constrain_coordinates(step, &dvdl_constr, ir, mdatoms, state,
                                  fr->bMolPBC,
                                  &top->idef, shake_vir,
                                  cr, ms, nrnb, wcycle, upd, constr,
                                  bCalcVir, do_log, do_ene);
            update_sd_second_half(step, &dvdl_constr, ir, mdatoms, state,
                                  fr->bMolPBC, f, &top->idef, cr, ms,
                                  nrnb, wcycle, upd, constr, do_log, do_ene);
            finish_update(ir, mdatoms,
                          state, graph,
                          nrnb, wcycle, upd, constr);
            {
                enerd->term[F_DVDL_CONSTR] += dvdl_constr;
            }
        }

        /* ############## IF NOT VV, Calculate globals HERE  ############ */
        /* With Leap-Frog we can skip compute_globals at
         * non-communication steps, but we need to calculate
         * the kinetic energy one step before communication.
         */
        {
            // Organize to do inter-simulation signalling on steps if
            // and when algorithms require it.
            bool doInterSimSignal = false;

            if (bGStat || (do_per_step(step+1, nstglobalcomm)) || doInterSimSignal)
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
                                | CGLO_ENERGY
                                | (bStopCM ? CGLO_STOPCM : 0)
                                | CGLO_TEMPERATURE
                                | CGLO_PRESSURE
                                | CGLO_CONSTRAINT
                                | (shouldCheckNumberOfBondedInteractions ? CGLO_CHECK_NUMBER_OF_BONDED_INTERACTIONS : 0)
                                );
                checkNumberOfBondedInteractions(fplog, cr, totalNumberOfBondedInteractions,
                                                top_global, top, state,
                                                &shouldCheckNumberOfBondedInteractions);
            }
        }

        /* #############  END CALC EKIN AND PRESSURE ################# */

        /* Note: this is OK, but there are some numerical precision issues with using the convergence of
           the virial that should probably be addressed eventually. state->veta has better properies,
           but what we actually need entering the new cycle is the new shake_vir value. Ideally, we could
           generate the new shake_vir, but test the veta value for convergence.  This will take some thought. */

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
            enerd->term[F_ETOT] = enerd->term[F_EPOT] + enerd->term[F_EKIN];

            if (integratorHasConservedEnergyQuantity(ir))
            {
                {
                    enerd->term[F_ECONSERVED] = enerd->term[F_ETOT] + NPT_energy(ir, state, &MassQ);
                }
            }
            /* #########  END PREPARING EDR OUTPUT  ###########  */
        }

        /* Output stuff */
        if (MASTER(cr))
        {
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
                       eprNORMAL, mdebin, fcd, groups, &(ir->opts), ir->awh);

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
        /* Print the remaining wall clock time for the run */
        if (isMasterSimMasterRank(ms, cr) &&
            (do_verbose || gmx_got_usr_signal()) &&
            !bPMETunePrinting)
        {
            print_time(stderr, walltime_accounting, step, ir, cr);
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

        step++;
        step_rel++;

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
            reset_all_counters(fplog, mdlog, cr, step, &step_rel, ir, wcycle, nrnb, walltime_accounting,
                               use_GPU(fr->nbv) ? fr->nbv : nullptr, fr->pmedata);
            wcycle_set_reset_counters(wcycle, -1);
            if (!thisRankHasDuty(cr, DUTY_PME))
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

    }
    /* End of main MD loop */

    /* Closing TNG files can include compressing data. Therefore it is good to do that
     * before stopping the time measurements. */
    mdoutf_tng_close(outf);

    /* Stop measuring walltime */
    walltime_accounting_end(walltime_accounting);

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
                       eprAVER, mdebin, fcd, groups, &(ir->opts), ir->awh);
        }
    }
    done_mdebin(mdebin);
    done_mdoutf(outf);

    if (bPMETune)
    {
        pme_loadbal_done(pme_loadbal, fplog, mdlog, use_GPU(fr->nbv));
    }

    done_shellfc(fplog, shellfc, step_rel);

    // Clean up swapcoords
    if (ir->eSwapCoords != eswapNO)
    {
        finish_swapcoords(ir->swap);
    }

    walltime_accounting_set_nsteps_done(walltime_accounting, step_rel);
    if (step_rel >= wcycle_get_reset_counters(wcycle) &&
        signals[eglsRESETCOUNTERS].set == 0 &&
        !bResetCountersHalfMaxH)
    {
        walltime_accounting_set_valid_finish(walltime_accounting);
    }

    destroy_enerdata(enerd);
    sfree(enerd);
    mdAlgorithmsTearDownAtomData(fr->bonded_threading, top);
    fr->bonded_threading = nullptr;
    sfree(top);
}


void gmx::Integrator::do_custom_2()
{
    /* This creates a NVE integrator, having an outer loop handling
     * parallelization (load balancing, DD) and pair listing, and an
     * inner loop doing the actual propagation and trajectory writing.
     *
     * Note that trajectory writing should not be in the inner loop,
     * but in an outer loop depending on the chosen writing frequency.
     */

    DataManager data(*this);
    LoopElement outerLoop(*this, data);

    int         nreps_inner = 1;
    if (inputrec->nstlist > 0)
    {
        nreps_inner = inputrec->nstlist;

        if (data.bPMETune)
        {
            outerLoop.add(new PMELoadBalanceElement(*this, data));
        }

        outerLoop.add(new DDElement(*this, data));
    }

    auto innerLoop = new LoopElement(*this, data, nreps_inner, true);
    innerLoop->add(new ForceElement(*this, data));
    innerLoop->add(new WriteTrajectoryElement(*this, data));
    innerLoop->add(new UpdateFusedElement(*this, data));
    if (constr)
    {
        innerLoop->add(new ConstrainElement(*this, data));
        //innerLoop->add(constr->makeFusedConstrainer(state, &state->readX, &state->readV, &state->writeX))
    }
    innerLoop->add(new ComputeGlobalsElement(*this, data));
    innerLoop->add(new WriteEnergyElement(*this, data));

    outerLoop.add(innerLoop);

    while (!data.bLastStep)
    {
        outerLoop.run();
    }

}
