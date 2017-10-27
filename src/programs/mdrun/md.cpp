/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011,2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include "md.h"

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
#include "gromacs/fda/FDA.h"
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
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_gpu_data_mgmt.h"
#include "gromacs/mdlib/ns.h"
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

#include "deform.h"
#include "membed.h"
#include "repl_ex.h"

#ifdef GMX_FAHCORE
#include "corewrap.h"
#endif

using gmx::SimulationSignaller;

/*! \brief Check whether bonded interactions are missing, if appropriate
 *
 * \param[in]    fplog                                  Log file pointer
 * \param[in]    cr                                     Communication object
 * \param[in]    totalNumberOfBondedInteractions        Result of the global reduction over the number of bonds treated in each domain
 * \param[in]    top_global                             Global topology for the error message
 * \param[in]    top_local                              Local topology for the error message
 * \param[in]    state                                  Global state for the error message
 * \param[inout] shouldCheckNumberOfBondedInteractions  Whether we should do the check.
 *
 * \return Nothing, except that shouldCheckNumberOfBondedInteractions
 * is always set to false after exit.
 */
static void checkNumberOfBondedInteractions(FILE *fplog, t_commrec *cr, int totalNumberOfBondedInteractions,
                                            gmx_mtop_t *top_global, gmx_localtop_t *top_local, t_state *state,
                                            bool *shouldCheckNumberOfBondedInteractions)
{
    if (*shouldCheckNumberOfBondedInteractions)
    {
        if (totalNumberOfBondedInteractions != cr->dd->nbonded_global)
        {
            dd_print_missing_interactions(fplog, cr, totalNumberOfBondedInteractions, top_global, top_local, state); // Does not return
        }
        *shouldCheckNumberOfBondedInteractions = false;
    }
}

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

/*! \brief Copy the state from \p rerunFrame to \p globalState and, if requested, construct vsites
 *
 * \param[in]     rerunFrame      The trajectory frame to compute energy/forces for
 * \param[in,out] globalState     The global state container
 * \param[in]     constructVsites When true, vsite coordinates are constructed
 * \param[in]     vsite           Vsite setup, can be nullptr when \p constructVsites = false
 * \param[in]     idef            Topology parameters, used for constructing vsites
 * \param[in]     timeStep        Time step, used for constructing vsites
 * \param[in]     forceRec        Force record, used for constructing vsites
 * \param[in,out] graph           The molecular graph, used for constructing vsites when != nullptr
 * \param[in,out] warnWhenNoV     When true, issue a warning when no velocities are present in \p rerunFrame; is set to false when a warning was issued
 */
static void prepareRerunState(const t_trxframe  &rerunFrame,
                              t_state           *globalState,
                              bool               constructVsites,
                              const gmx_vsite_t *vsite,
                              const t_idef      &idef,
                              double             timeStep,
                              const t_forcerec  &forceRec,
                              t_graph           *graph,
                              gmx_bool          *warnWhenNoV)
{
    for (int i = 0; i < globalState->natoms; i++)
    {
        copy_rvec(rerunFrame.x[i], globalState->x[i]);
    }
    if (rerunFrame.bV)
    {
        for (int i = 0; i < globalState->natoms; i++)
        {
            copy_rvec(rerunFrame.v[i], globalState->v[i]);
        }
    }
    else
    {
        for (int i = 0; i < globalState->natoms; i++)
        {
            clear_rvec(globalState->v[i]);
        }
        if (*warnWhenNoV)
        {
            fprintf(stderr, "\nWARNING: Some frames do not contain velocities.\n"
                    "         Ekin, temperature and pressure are incorrect,\n"
                    "         the virial will be incorrect when constraints are present.\n"
                    "\n");
            *warnWhenNoV = FALSE;
        }
    }
    copy_mat(rerunFrame.box, globalState->box);

    if (constructVsites)
    {
        GMX_ASSERT(vsite, "Need valid vsite for constructing vsites");

        if (graph)
        {
            /* Following is necessary because the graph may get out of sync
             * with the coordinates if we only have every N'th coordinate set
             */
            mk_mshift(nullptr, graph, forceRec.ePBC, globalState->box, as_rvec_array(globalState->x.data()));
            shift_self(graph, globalState->box, as_rvec_array(globalState->x.data()));
        }
        construct_vsites(vsite, as_rvec_array(globalState->x.data()), timeStep, as_rvec_array(globalState->v.data()),
                         idef.iparams, idef.il,
                         forceRec.ePBC, forceRec.bMolPBC, nullptr, globalState->box);
        if (graph)
        {
            unshift_self(graph, globalState->box, as_rvec_array(globalState->x.data()));
        }
    }
}

/*! \libinternal
    \copydoc integrator_t (FILE *fplog, t_commrec *cr, const gmx::MDLogger &mdlog,
                           int nfile, const t_filenm fnm[],
                           const gmx_output_env_t *oenv,
                           const MdrunOptions &mdrunOptions,
                           gmx_vsite_t *vsite, gmx_constr_t constr,
                           gmx::IMDOutputProvider *outputProvider,
                           t_inputrec *inputrec,
                           gmx_mtop_t *top_global, t_fcdata *fcd,
                           t_state *state_global,
                           t_mdatoms *mdatoms,
                           t_nrnb *nrnb, gmx_wallcycle_t wcycle,
                           t_forcerec *fr,
                           const ReplicaExchangeParameters &replExParams,
                           gmx_membed_t *membed,
                           gmx_walltime_accounting_t walltime_accounting)
 */
double gmx::do_md(FILE *fplog, t_commrec *cr, const gmx::MDLogger &mdlog,
                  int nfile, const t_filenm fnm[],
                  const gmx_output_env_t *oenv,
                  const MdrunOptions &mdrunOptions,
                  gmx_vsite_t *vsite, gmx_constr_t constr,
                  gmx::IMDOutputProvider *outputProvider,
                  t_inputrec *ir,
                  gmx_mtop_t *top_global,
                  t_fcdata *fcd,
                  t_state *state_global,
                  ObservablesHistory *observablesHistory,
                  gmx::MDAtoms *mdAtoms,
                  t_nrnb *nrnb, gmx_wallcycle_t wcycle,
                  t_forcerec *fr,
                  const ReplicaExchangeParameters &replExParams,
                  gmx_membed_t *membed,
                  gmx_walltime_accounting_t walltime_accounting)
{
    gmx_mdoutf_t    outf = nullptr;
    gmx_int64_t     step, step_rel;
    double          elapsed_time;
    double          t, t0, lam0[efptNR];
    gmx_bool        bGStatEveryStep, bGStat, bCalcVir, bCalcEnerStep, bCalcEner;
    gmx_bool        bNS, bNStList, bSimAnn, bStopCM,
                    bFirstStep, bInitStep, bLastStep = FALSE,
                    bBornRadii, bUsingEnsembleRestraints;
    gmx_bool          bDoDHDL = FALSE, bDoFEP = FALSE, bDoExpanded = FALSE;
    gmx_bool          do_ene, do_log, do_verbose, bRerunWarnNoV = TRUE,
                      bForceUpdate = FALSE, bCPT;
    gmx_bool          bMasterState;
    int               force_flags, cglo_flags;
    tensor            force_vir, shake_vir, total_vir, tmp_vir, pres;
    int               i, m;
    t_trxstatus      *status;
    rvec              mu_tot;
    t_vcm            *vcm;
    matrix            parrinellorahmanMu, M;
    t_trxframe        rerun_fr;
    gmx_repl_ex_t     repl_ex = nullptr;
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
    gmx_bool          bSumEkinhOld, bDoReplEx, bExchanged, bNeedRepartition;
    gmx_bool          bResetCountersHalfMaxH = FALSE;
    gmx_bool          bTemp, bPres, bTrotter;
    real              dvdl_constr;
    rvec             *cbuf        = nullptr;
    int               cbuf_nalloc = 0;
    matrix            lastbox;
    int               lamnew  = 0;
    /* for FEP */
    int               nstfep = 0;
    double            cycles;
    real              saved_conserved_quantity = 0;
    real              last_ekin                = 0;
    t_extmass         MassQ;
    int             **trotter_seq;
    char              sbuf[STEPSTRSIZE], sbuf2[STEPSTRSIZE];
    int               handled_stop_condition = gmx_stop_cond_none; /* compare to get_stop_condition*/


    /* PME load balancing data for GPU kernels */
    pme_load_balancing_t *pme_loadbal      = nullptr;
    gmx_bool              bPMETune         = FALSE;
    gmx_bool              bPMETunePrinting = FALSE;

    /* Interactive MD */
    gmx_bool          bIMDstep = FALSE;

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
    SimulationSignaller nullSignaller(nullptr, nullptr, false, false);

    if (mdrunOptions.timingOptions.resetHalfway)
    {
        if (ir->nsteps > 0)
        {
            /* Signal to reset the counters half the simulation steps. */
            wcycle_set_reset_counters(wcycle, ir->nsteps/2);
        }
        /* Signal to reset the counters halfway the simulation time. */
        bResetCountersHalfMaxH = (mdrunOptions.maximumHoursToRun > 0);
    }

    /* md-vv uses averaged full step velocities for T-control
       md-vv-avek uses averaged half step velocities for T-control (but full step ekin for P control)
       md uses averaged half step kinetic energies to determine temperature unless defined otherwise by GMX_EKIN_AVE_VEL; */
    bTrotter = (EI_VV(ir->eI) && (inputrecNptTrotter(ir) || inputrecNphTrotter(ir) || inputrecNvtTrotter(ir)));

    const gmx_bool bRerunMD      = mdrunOptions.rerun;
    int            nstglobalcomm = mdrunOptions.globalCommunicationInterval;

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

    if (bRerunMD)
    {
        ir->nstxout_compressed = 0;
    }
    groups = &top_global->groups;

    gmx_edsam *ed = nullptr;
    if (opt2bSet("-ei", nfile, fnm) || observablesHistory->edsamHistory != nullptr)
    {
        /* Initialize essential dynamics sampling */
        ed = init_edsam(opt2fn_null("-ei", nfile, fnm), opt2fn("-eo", nfile, fnm),
                        top_global,
                        ir, cr, constr,
                        state_global, observablesHistory,
                        oenv, mdrunOptions.continuationOptions.appendFiles);
    }

    if (ir->eSwapCoords != eswapNO)
    {
        /* Initialize ion swapping code */
        init_swapcoords(fplog, ir, opt2fn_master("-swap", nfile, fnm, cr),
                        top_global,
                        state_global, observablesHistory,
                        cr, oenv, mdrunOptions);
    }

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
    if (shellfc && ir->bDoAwh)
    {
        gmx_fatal(FARGS, "AWH biasing does not support shell particles.");
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
#ifdef BUILD_WITH_FDA
        fda::FDASettings fda_settings = fr->fda->get_settings();
        top = gmx_mtop_generate_local_top(top_global, ir->efep != efepNO, &fda_settings);
#endif

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

    /* Set up interactive MD (IMD) */
    init_IMD(ir, cr, top_global, fplog, ir->nstcalcenergy, MASTER(cr) ? as_rvec_array(state_global->x.data()) : nullptr,
             nfile, fnm, oenv, mdrunOptions);

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

    if (ir->bExpanded)
    {
        init_expanded_ensemble(startingFromCheckpoint, ir, state->dfhist);
    }

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

    /* Initialize AWH and restore state from history in checkpoint if needed. */
    if (ir->bDoAwh)
    {
        ir->awh = new gmx::Awh(fplog, *ir, cr, *ir->awhParams, opt2fn("-awh", nfile, fnm), ir->pull_work);

        if (startingFromCheckpoint)
        {
            /* Restore the AWH history read from checkpoint */
            ir->awh->restoreStateFromHistory(MASTER(cr) ? state_global->awhHistory.get() : nullptr);
        }
        else if (MASTER(cr))
        {
            /* Initialize the AWH history here */
            state_global->awhHistory = ir->awh->initHistoryFromState();
        }
    }

    const bool useReplicaExchange = (replExParams.exchangeInterval > 0);
    if (useReplicaExchange && MASTER(cr))
    {
        repl_ex = init_replica_exchange(fplog, cr->ms, top_global->natoms, ir,
                                        replExParams);
    }

    /* PME tuning is only supported in the Verlet scheme, with PME for
     * Coulomb. It is not supported with only LJ PME, or for
     * reruns. */
    bPMETune = (mdrunOptions.tunePme && EEL_PME(fr->ic->eeltype) && !bRerunMD &&
                !mdrunOptions.reproducible && ir->cutoff_scheme != ecutsGROUP);
    if (bPMETune)
    {
        pme_loadbal_init(&pme_loadbal, cr, mdlog, ir, state->box,
                         fr->ic, fr->nbv->listParams.get(), fr->pmedata, use_GPU(fr->nbv),
                         &bPMETunePrinting);
    }

    if (!ir->bContinuation && !bRerunMD)
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
                               cr, nrnb, fr, top);
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

    cglo_flags = (CGLO_TEMPERATURE | CGLO_GSTAT
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
    checkNumberOfBondedInteractions(fplog, cr, totalNumberOfBondedInteractions,
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
                        constr_rmsd(constr));
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

        if (bRerunMD)
        {
            fprintf(stderr, "starting md rerun '%s', reading coordinates from"
                    " input trajectory '%s'\n\n",
                    *(top_global->name), opt2fn("-rerun", nfile, fnm));
            if (mdrunOptions.verbose)
            {
                fprintf(stderr, "Calculated time to finish depends on nsteps from "
                        "run input file,\nwhich may not correspond to the time "
                        "needed to process input trajectory.\n\n");
            }
        }
        else
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

    /* if rerunMD then read coordinates and velocities from input trajectory */
    if (bRerunMD)
    {
        if (getenv("GMX_FORCE_UPDATE"))
        {
            bForceUpdate = TRUE;
        }

        rerun_fr.natoms = 0;
        if (MASTER(cr))
        {
            bLastStep = !read_first_frame(oenv, &status,
                                          opt2fn("-rerun", nfile, fnm),
                                          &rerun_fr, TRX_NEED_X | TRX_READ_V);
            if (rerun_fr.natoms != top_global->natoms)
            {
                gmx_fatal(FARGS,
                          "Number of atoms in trajectory (%d) does not match the "
                          "run input file (%d)\n",
                          rerun_fr.natoms, top_global->natoms);
            }
            if (ir->ePBC != epbcNONE)
            {
                if (!rerun_fr.bBox)
                {
                    gmx_fatal(FARGS, "Rerun trajectory frame step %d time %f does not contain a box, while pbc is used", rerun_fr.step, rerun_fr.time);
                }
                if (max_cutoff2(ir->ePBC, rerun_fr.box) < gmx::square(fr->rlist))
                {
                    gmx_fatal(FARGS, "Rerun trajectory frame step %d time %f has too small box dimensions", rerun_fr.step, rerun_fr.time);
                }
            }
        }

        if (PAR(cr))
        {
            rerun_parallel_comm(cr, &rerun_fr, &bLastStep);
        }

        if (ir->ePBC != epbcNONE)
        {
            /* Set the shift vectors.
             * Necessary here when have a static box different from the tpr box.
             */
            calc_shifts(rerun_fr.box, fr->shift_vec);
        }
    }

    /* Loop over MD steps or if rerunMD to end of input trajectory,
     * or, if max_hours>0, until max_hours is reached.
     */
    real max_hours   = mdrunOptions.maximumHoursToRun;
    bFirstStep       = TRUE;
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
        bool checkpointIsLocal    = !useReplicaExchange && !bUsingEnsembleRestraints;
        bool stopConditionIsLocal = !useReplicaExchange && !bUsingEnsembleRestraints;
        bool resetCountersIsLocal = true;
        signals[eglsCHKPT]         = SimulationSignal(checkpointIsLocal);
        signals[eglsSTOPCOND]      = SimulationSignal(stopConditionIsLocal);
        signals[eglsRESETCOUNTERS] = SimulationSignal(resetCountersIsLocal);
    }

    DdOpenBalanceRegionBeforeForceComputation ddOpenBalanceRegion   = (DOMAINDECOMP(cr) ? DdOpenBalanceRegionBeforeForceComputation::yes : DdOpenBalanceRegionBeforeForceComputation::no);
    DdCloseBalanceRegionAfterForceComputation ddCloseBalanceRegion  = (DOMAINDECOMP(cr) ? DdCloseBalanceRegionAfterForceComputation::yes : DdCloseBalanceRegionAfterForceComputation::no);

    step     = ir->init_step;
    step_rel = 0;

    // TODO extract this to new multi-simulation module
    if (MASTER(cr) && MULTISIM(cr) && !useReplicaExchange)
    {
        if (!multisim_int_all_are_equal(cr->ms, ir->nsteps))
        {
            GMX_LOG(mdlog.warning).appendText(
                    "Note: The number of steps is not consistent across multi simulations,\n"
                    "but we are proceeding anyway!");
        }
        if (!multisim_int_all_are_equal(cr->ms, ir->init_step))
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
                           ir, fr, state,
                           wcycle,
                           step, step_rel,
                           &bPMETunePrinting);
        }

        wallcycle_start(wcycle, ewcSTEP);

        if (bRerunMD)
        {
            if (rerun_fr.bStep)
            {
                step     = rerun_fr.step;
                step_rel = step - ir->init_step;
            }
            if (rerun_fr.bTime)
            {
                t = rerun_fr.time;
            }
            else
            {
                t = step;
            }
        }
        else
        {
            bLastStep = (step_rel == ir->nsteps);
            t         = t0 + step*ir->delta_t;
        }

        // TODO Refactor this, so that nstfep does not need a default value of zero
        if (ir->efep != efepNO || ir->bSimTemp)
        {
            /* find and set the current lambdas.  If rerunning, we either read in a state, or a lambda value,
               requiring different logic. */
            if (bRerunMD)
            {
                if (MASTER(cr))
                {
                    setCurrentLambdasRerun(step, ir->fepvals, &rerun_fr, lam0, state_global);
                }
            }
            else
            {
                setCurrentLambdasLocal(step, ir->fepvals, lam0, state);
            }
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

        if (bRerunMD && MASTER(cr))
        {
            const bool constructVsites = (vsite && mdrunOptions.rerunConstructVsites);
            if (constructVsites && DOMAINDECOMP(cr))
            {
                gmx_fatal(FARGS, "Vsite recalculation with -rerun is not implemented with domain decomposition, use a single rank");
            }
            prepareRerunState(rerun_fr, state_global, constructVsites, vsite, top->idef, ir->delta_t, *fr, graph, &bRerunWarnNoV);
        }

        /* Stop Center of Mass motion */
        bStopCM = (ir->comm_mode != ecmNO && do_per_step(step, ir->nstcomm));

        if (bRerunMD)
        {
            /* for rerun MD always do Neighbour Searching */
            bNS      = (bFirstStep || ir->nstlist != 0);
        }
        else
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
        do_log     = do_per_step(step, ir->nstlog) || (bFirstStep && !startingFromCheckpoint) || bLastStep || bRerunMD;
        do_verbose = mdrunOptions.verbose &&
            (step % mdrunOptions.verboseStepPrintInterval == 0 || bFirstStep || bLastStep || bRerunMD);

        if (bNS && !(bFirstStep && ir->bContinuation && !bRerunMD))
        {
            if (bRerunMD)
            {
                bMasterState = TRUE;
            }
            else
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

        if (ir->efep != efepNO)
        {
            update_mdatoms(mdatoms, state->lambda[efptMASS]);
        }

        if ((bRerunMD && rerun_fr.bV) || bExchanged)
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
        clear_mat(force_vir);

        /* We write a checkpoint at this MD step when:
         * either at an NS step when we signalled through gs,
         * or at the last step (but not when we do not want confout),
         * but never at the first step or with rerun.
         */
        bCPT = (((signals[eglsCHKPT].set && (bNS || ir->nstlist == 0)) ||
                 (bLastStep && mdrunOptions.writeConfout)) &&
                step > ir->init_step && !bRerunMD);
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

        do_ene = (do_per_step(step, ir->nstenergy) || bLastStep || bRerunMD);

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
                       ((inputrecDynamicBox(ir) || bRerunMD) ? GMX_FORCE_DYNAMICBOX : 0) |
                       GMX_FORCE_ALLFORCES |
                       (bCalcVir ? GMX_FORCE_VIRIAL : 0) |
                       (bCalcEner ? GMX_FORCE_ENERGY : 0) |
                       (bDoFEP ? GMX_FORCE_DHDL : 0)
                       );

        if (shellfc)
        {
            /* Now is the time to relax the shells */
            relax_shell_flexcon(fplog, cr, mdrunOptions.verbose, step,
                                ir, bNS, force_flags, top,
                                constr, enerd, fcd,
                                state, &f, force_vir, mdatoms,
                                nrnb, wcycle, graph, groups,
                                shellfc, fr, bBornRadii, t, mu_tot,
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
            if (ir->bDoAwh && bCPT && MASTER(cr))
            {
                ir->awh->updateHistory(state_global->awhHistory.get());
            }

            /* The coordinates (x) are shifted (to get whole molecules)
             * in do_force.
             * This is parallellized as well, and does communication too.
             * Check comments in sim_util.c
             */
            do_force(fplog, cr, ir, step, nrnb, wcycle, top, groups,
                     state->box, state->x, &state->hist,
                     f, force_vir, mdatoms, enerd, fcd,
                     state->lambda, graph,
                     fr, vsite, mu_tot, t, ed, bBornRadii,
                     (bNS ? GMX_FORCE_NS : 0) | force_flags,
                     ddOpenBalanceRegion, ddCloseBalanceRegion);
        }

        if (EI_VV(ir->eI) && !startingFromCheckpoint && !bRerunMD)
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
                copy_rvecn(as_rvec_array(state->v.data()), vbuf, 0, state->natoms); /* should make this better for parallelizing? */
            }
            else
            {
                /* this is for NHC in the Ekin(t+dt/2) version of vv */
                trotter_update(ir, step, ekind, enerd, state, total_vir, mdatoms, &MassQ, trotter_seq, ettTSEQ1);
            }

            update_coords(fplog, step, ir, mdatoms, state, f, fcd,
                          ekind, M, upd, etrtVELOCITY1,
                          cr, constr);

            if (!bRerunMD || rerun_fr.bV || bForceUpdate)         /* Why is rerun_fr.bV here?  Unclear. */
            {
                wallcycle_stop(wcycle, ewcUPDATE);
                update_constraints(fplog, step, nullptr, ir, mdatoms,
                                   state, fr->bMolPBC, graph, f,
                                   &top->idef, shake_vir,
                                   cr, nrnb, wcycle, upd, constr,
                                   TRUE, bCalcVir);
                wallcycle_start(wcycle, ewcUPDATE);
            }
            else if (graph)
            {
                /* Need to unshift here if a do_force has been
                   called in the previous step */
                unshift_self(graph, state->box, as_rvec_array(state->x.data()));
            }
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
                                | CGLO_ENERGY
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
                checkNumberOfBondedInteractions(fplog, cr, totalNumberOfBondedInteractions,
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
                copy_rvecn(vbuf, as_rvec_array(state->v.data()), 0, state->natoms);
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
            if (ir->efep != efepNO && !bRerunMD)
            {
                sum_dhdl(enerd, state->lambda, ir->fepvals);
            }
        }

        /* ########  END FIRST UPDATE STEP  ############## */
        /* ########  If doing VV, we now have v(dt) ###### */

#ifdef BUILD_WITH_FDA
        fr->fda->save_and_write_scalar_time_averages(state->x, top_global);
#endif

        if (bDoExpanded)
        {
            /* perform extended ensemble sampling in lambda - we don't
               actually move to the new state before outputting
               statistics, but if performing simulated tempering, we
               do update the velocities and the tau_t. */

            lamnew = ExpandedEnsembleDynamics(fplog, ir, enerd, state, &MassQ, state->fep_state, state->dfhist, step, as_rvec_array(state->v.data()), mdatoms);
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
                                 &nchkpt,
                                 bCPT, bRerunMD, bLastStep,
                                 mdrunOptions.writeConfout,
                                 bSumEkinhOld);
        /* Check if IMD step and do IMD communication, if bIMD is TRUE. */
        bIMDstep = do_IMD(ir->bIMD, step, cr, bNS, state->box, as_rvec_array(state->x.data()), ir, t, wcycle);

        /* kludge -- virial is lost with restart for MTTK NPT control. Must reload (saved earlier). */
        if (startingFromCheckpoint && (inputrecNptTrotter(ir) || inputrecNphTrotter(ir)))
        {
            copy_mat(state->svir_prev, shake_vir);
            copy_mat(state->fvir_prev, force_vir);
        }

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
                nsteps_stop               = std::max(ir->nstlist, 2*nstglobalcomm);
            }
            else if (gmx_get_stop_condition() == gmx_stop_cond_next)
            {
                /* Stop directly after the next global communication step.
                 * This breaks exact continuation.
                 */
                signals[eglsSTOPCOND].sig = -1;
                nsteps_stop               = nstglobalcomm + 1;
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

        /* at the start of step, randomize or scale the velocities ((if vv. Restriction of Andersen controlled
           in preprocessing */

        if (ETC_ANDERSEN(ir->etc)) /* keep this outside of update_tcouple because of the extra info required to pass */
        {
            gmx_bool bIfRandomize;
            bIfRandomize = update_randomize_velocities(ir, step, cr, mdatoms, state, upd, constr);
            /* if we have constraints, we have to remove the kinetic energy parallel to the bonds */
            if (constr && bIfRandomize)
            {
                update_constraints(fplog, step, nullptr, ir, mdatoms,
                                   state, fr->bMolPBC, graph, f,
                                   &top->idef, tmp_vir,
                                   cr, nrnb, wcycle, upd, constr,
                                   TRUE, bCalcVir);
            }
        }
        /* Box is changed in update() when we do pressure coupling,
         * but we should still use the old box for energy corrections and when
         * writing it to the energy file, so it matches the trajectory files for
         * the same timestep above. Make a copy in a separate array.
         */
        copy_mat(state->box, lastbox);

        dvdl_constr = 0;

        if (!bRerunMD || rerun_fr.bV || bForceUpdate)
        {
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
                update_coords(fplog, step, ir, mdatoms, state, f, fcd,
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

            update_coords(fplog, step, ir, mdatoms, state, f, fcd,
                          ekind, M, upd, etrtPOSITION, cr, constr);
            wallcycle_stop(wcycle, ewcUPDATE);

            update_constraints(fplog, step, &dvdl_constr, ir, mdatoms, state,
                               fr->bMolPBC, graph, f,
                               &top->idef, shake_vir,
                               cr, nrnb, wcycle, upd, constr,
                               FALSE, bCalcVir);

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

                update_coords(fplog, step, ir, mdatoms, state, f, fcd,
                              ekind, M, upd, etrtPOSITION, cr, constr);
                wallcycle_stop(wcycle, ewcUPDATE);

                /* do we need an extra constraint here? just need to copy out of as_rvec_array(state->v.data()) to upd->xp? */
                /* are the small terms in the shake_vir here due
                 * to numerical errors, or are they important
                 * physically? I'm thinking they are just errors, but not completely sure.
                 * For now, will call without actually constraining, constr=NULL*/
                update_constraints(fplog, step, nullptr, ir, mdatoms,
                                   state, fr->bMolPBC, graph, f,
                                   &top->idef, tmp_vir,
                                   cr, nrnb, wcycle, upd, nullptr,
                                   FALSE, bCalcVir);
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
        }
        else if (graph)
        {
            /* Need to unshift here */
            unshift_self(graph, state->box, as_rvec_array(state->x.data()));
        }

        if (vsite != nullptr)
        {
            wallcycle_start(wcycle, ewcVSITECONSTR);
            if (graph != nullptr)
            {
                shift_self(graph, state->box, as_rvec_array(state->x.data()));
            }
            construct_vsites(vsite, as_rvec_array(state->x.data()), ir->delta_t, as_rvec_array(state->v.data()),
                             top->idef.iparams, top->idef.il,
                             fr->ePBC, fr->bMolPBC, cr, state->box);

            if (graph != nullptr)
            {
                unshift_self(graph, state->box, as_rvec_array(state->x.data()));
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
            bool doInterSimSignal = (!bFirstStep && bDoReplEx) || bUsingEnsembleRestraints;

            if (bGStat || (!EI_VV(ir->eI) && do_per_step(step+1, nstglobalcomm)) || doInterSimSignal)
            {
                // Since we're already communicating at this step, we
                // can propagate intra-simulation signals. Note that
                // check_nstglobalcomm has the responsibility for
                // choosing the value of nstglobalcomm that is one way
                // bGStat becomes true, so we can't get into a
                // situation where e.g. checkpointing can't be
                // signalled.
                bool                doIntraSimSignal = true;
                SimulationSignaller signaller(&signals, cr, doInterSimSignal, doIntraSimSignal);

                compute_globals(fplog, gstat, cr, ir, fr, ekind, state, mdatoms, nrnb, vcm,
                                wcycle, enerd, force_vir, shake_vir, total_vir, pres, mu_tot,
                                constr, &signaller,
                                lastbox,
                                &totalNumberOfBondedInteractions, &bSumEkinhOld,
                                (bGStat ? CGLO_GSTAT : 0)
                                | (!EI_VV(ir->eI) || bRerunMD ? CGLO_ENERGY : 0)
                                | (!EI_VV(ir->eI) && bStopCM ? CGLO_STOPCM : 0)
                                | (!EI_VV(ir->eI) ? CGLO_TEMPERATURE : 0)
                                | (!EI_VV(ir->eI) || bRerunMD ? CGLO_PRESSURE : 0)
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

        if (ir->efep != efepNO && (!EI_VV(ir->eI) || bRerunMD))
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
                                             bRerunMD ? rerun_fr.x   : as_rvec_array(state->x.data()),
                                             bRerunMD ? rerun_fr.box : state->box,
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
            bExchanged = replica_exchange(fplog, cr, repl_ex,
                                          state_global, enerd,
                                          state, step, t);
        }

        if ( (bExchanged || bNeedRepartition) && DOMAINDECOMP(cr) )
        {
            dd_partition_system(fplog, step, cr, TRUE, 1,
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

        if (bRerunMD)
        {
            if (MASTER(cr))
            {
                /* read next frame from input trajectory */
                bLastStep = !read_next_frame(oenv, status, &rerun_fr);
            }

            if (PAR(cr))
            {
                rerun_parallel_comm(cr, &rerun_fr, &bLastStep);
            }
        }

        cycles = wallcycle_stop(wcycle, ewcSTEP);
        if (DOMAINDECOMP(cr) && wcycle)
        {
            dd_cycles_add(cr->dd, cycles, ddCyclStep);
        }

        if (!bRerunMD || !rerun_fr.bStep)
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

        /* If bIMD is TRUE, the master updates the IMD energy record and sends positions to VMD client */
        IMD_prep_energies_send_positions(ir->bIMD && MASTER(cr), bIMDstep, ir->imd, enerd, step, bCalcEner, wcycle);

    }
    /* End of main MD loop */

    /* Closing TNG files can include compressing data. Therefore it is good to do that
     * before stopping the time measurements. */
    mdoutf_tng_close(outf);

    /* Stop measuring walltime */
    walltime_accounting_end(walltime_accounting);

    if (bRerunMD && MASTER(cr))
    {
        close_trx(status);
    }

#ifdef BUILD_WITH_FDA
    fr->fda->write_scalar_time_averages();
#endif

    if (!thisRankHasDuty(cr, DUTY_PME))
    {
        /* Tell the PME only node to finish */
        gmx_pme_send_finish(cr);
    }

    if (MASTER(cr))
    {
        if (ir->nstcalcenergy > 0 && !bRerunMD)
        {
            print_ebin(mdoutf_get_fp_ene(outf), FALSE, FALSE, FALSE, fplog, step, t,
                       eprAVER, mdebin, fcd, groups, &(ir->opts), ir->awh);
        }
    }

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

    if (ir->bDoAwh)
    {
        delete ir->awh;
    }

    // Clean up swapcoords
    if (ir->eSwapCoords != eswapNO)
    {
        finish_swapcoords(ir->swap);
    }

    /* Do essential dynamics cleanup if needed. Close .edo file */
    done_ed(&ed);

    /* IMD cleanup, if bIMD is TRUE. */
    IMD_finalize(ir->bIMD, ir->imd);

    walltime_accounting_set_nsteps_done(walltime_accounting, step_rel);
    if (step_rel >= wcycle_get_reset_counters(wcycle) &&
        signals[eglsRESETCOUNTERS].set == 0 &&
        !bResetCountersHalfMaxH)
    {
        walltime_accounting_set_valid_finish(walltime_accounting);
    }

    return 0;
}
