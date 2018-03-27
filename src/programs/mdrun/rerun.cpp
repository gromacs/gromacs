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
#include "gmxpre.h"

#include "rerun.h"

#include "config.h"

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
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/imd/imd.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/compute_io.h"
#include "gromacs/mdlib/ebin.h"
#include "gromacs/mdlib/expanded.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/md_support.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/mdebin.h"
#include "gromacs/mdlib/mdoutf.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/mdsetup.h"
#include "gromacs/mdlib/ns.h"
#include "gromacs/mdlib/shellfc.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdlib/simulationsignal.h"
#include "gromacs/mdlib/trajectory_writing.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdtypes/awh-history.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/energyhistory.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/observableshistory.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/swap/swapcoords.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"

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
                              t_graph           *graph)
{
    for (int i = 0; i < globalState->natoms; i++)
    {
        copy_rvec(rerunFrame.x[i], globalState->x[i]);
    }
    for (int i = 0; i < globalState->natoms; i++)
    {
        clear_rvec(globalState->v[i]);
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
    \copydoc integrator_t (FILE *fplog, t_commrec *cr,
                           const gmx_multisim_t *ms,
                           const gmx::MDLogger &mdlog,
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
double gmx::do_rerun(FILE *fplog, t_commrec *cr,
                     const gmx_multisim_t *ms,
                     const gmx::MDLogger &mdlog,
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
    gmx_mdoutf_t      outf = nullptr;
    gmx_int64_t       step, step_rel;
    double            elapsed_time;
    double            t, t0, lam0[efptNR];
    gmx_bool          bSimAnn, bLastStep = FALSE;
    gmx_bool          bDoDHDL = FALSE, bDoFEP = FALSE;
    gmx_bool          do_verbose;
    int               force_flags, cglo_flags;
    tensor            force_vir, shake_vir, total_vir, pres;
    t_trxstatus      *status;
    rvec              mu_tot;
    t_vcm            *vcm;
    t_trxframe        rerun_fr;
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
    /* for FEP */
    int               nstfep = 0;
    double            cycles;
    char              sbuf[STEPSTRSIZE];
    int               handled_stop_condition = gmx_stop_cond_none; /* compare to get_stop_condition */

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

    if (constr) {
        GMX_LOG(mdlog.info).asParagraph().
                appendText("NOTE: Rerun does not recalculate constraints.");
    }

    /* Settings for rerun */
    ir->nstlist       = 1;
    ir->nstcalcenergy = 1;
    int nstglobalcomm = 1;

    ir->nstxout_compressed = 0;
    groups = &top_global->groups;
    // no kinetic energy (but bool passed to compute_globals and trajectory writing)
    gmx_bool bSumEkinhOld = FALSE;
    // no essential dynamics (but object needs to be passed to force calc)
    gmx_edsam *ed = nullptr;

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

    if (ir->bExpanded)
    {
        gmx_fatal(FARGS, "Expanded ensemble not supported by rerun.");
    }

    if (MASTER(cr))
    {
        if (observablesHistory->energyHistory.get() == nullptr)
        {
            observablesHistory->energyHistory = std::unique_ptr<energyhistory_t>(new energyhistory_t {});
        }
        /* Set the initial energy history in state by updating once */
        update_energyhistory(observablesHistory->energyHistory.get(), mdebin);
    }

    if (ir->efep != efepNO)
    {
        /* Set free energy calculation frequency as the greatest common
         * denominator of nstdhdl and repl_ex_nst. */
        nstfep = ir->fepvals->nstdhdl;
    }

    /* removing this 'CGLO_TEMPERATURE' leads to error using DD:
     *
    A list of missing interactions:
              exclusions of      0 missing      1
    -------------------------------------------------------
    Program:     gmx mdrun, version 2019-dev-20180305-27c754e-dirty
    Source file: src/gromacs/domdec/domdec_topology.cpp (line 437)
    MPI rank:    0 (out of 20)

    Fatal error:
    1 of the 0 bonded interactions could not be calculated because some atoms
    involved moved further apart than the multi-body cut-off distance (0.606042
    nm) or the two-body cut-off distance (1.196 nm), see option -rdd, for pairs
    and tabulated bonds also see option -ddcheck
    */

    cglo_flags = (CGLO_INITIALIZATION | CGLO_TEMPERATURE | CGLO_GSTAT);

    /* To minimize communication, compute_globals computes the COM velocity
     * and the kinetic energy for the velocities without COM motion removed.
     * Thus to get the kinetic energy without the COM contribution, we need
     * to call compute_globals twice.
     */
    compute_globals(fplog, gstat, cr, ir, fr, ekind, state, mdatoms, nrnb, vcm,
                    nullptr, enerd, force_vir, shake_vir, total_vir, pres, mu_tot,
                    constr, &nullSignaller, state->box,
                    &totalNumberOfBondedInteractions, &bSumEkinhOld,
                    cglo_flags | (shouldCheckNumberOfBondedInteractions ? CGLO_CHECK_NUMBER_OF_BONDED_INTERACTIONS : 0));
    checkNumberOfBondedInteractions(fplog, cr, totalNumberOfBondedInteractions,
                                    top_global, top, state,
                                    &shouldCheckNumberOfBondedInteractions);

    if (MASTER(cr))
    {
        fprintf(stderr, "starting md rerun '%s', reading coordinates from"
                        " input trajectory '%s'\n\n",
                *(top_global->name), opt2fn("-rerun", nfile, fnm));
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
    rerun_fr.natoms = 0;
    if (MASTER(cr))
    {
        bLastStep = !read_first_frame(oenv, &status,
                                      opt2fn("-rerun", nfile, fnm),
                                      &rerun_fr, TRX_NEED_X | TRX_READ_V);
        if (rerun_fr.bV)
        {
            GMX_LOG(mdlog.warning).asParagraph().appendText(
                    "NOTE: Rerun trajectory contains velocities. "
                            "Rerun does only evaluate potential energy and forces. "
                            "The velocities will be ignored.");
            // re-read first frame to have following frames ignoring velocities
            read_first_frame(oenv, &status,
                             opt2fn("-rerun", nfile, fnm),
                             &rerun_fr, TRX_NEED_X);
        }
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

    /* Loop over MD steps or if rerunMD to end of input trajectory,
     * or, if max_hours>0, until max_hours is reached.
     */
    real max_hours   = mdrunOptions.maximumHoursToRun;

    DdOpenBalanceRegionBeforeForceComputation ddOpenBalanceRegion   = (DOMAINDECOMP(cr) ? DdOpenBalanceRegionBeforeForceComputation::yes : DdOpenBalanceRegionBeforeForceComputation::no);
    DdCloseBalanceRegionAfterForceComputation ddCloseBalanceRegion  = (DOMAINDECOMP(cr) ? DdCloseBalanceRegionAfterForceComputation::yes : DdCloseBalanceRegionAfterForceComputation::no);

    step     = ir->init_step;
    step_rel = 0;

    /* and stop now if we should */
    bLastStep = (bLastStep || (ir->nsteps >= 0 && step_rel > ir->nsteps));
    while (!bLastStep)
    {

        wallcycle_start(wcycle, ewcSTEP);

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

        // TODO Refactor this, so that nstfep does not need a default value of zero
        if (ir->efep != efepNO || ir->bSimTemp)
        {
            if (MASTER(cr))
            {
                setCurrentLambdasRerun(step, ir->fepvals, &rerun_fr, lam0, state_global);
            }

            bDoDHDL      = do_per_step(step, ir->fepvals->nstdhdl);
            bDoFEP       = ((ir->efep != efepNO) && do_per_step(step, nstfep));
        }

        if (bSimAnn)
        {
            update_annealing_target_temp(ir, t, upd);
        }

        if (MASTER(cr)) {
            const bool constructVsites = (vsite && mdrunOptions.rerunConstructVsites);
            if (constructVsites && DOMAINDECOMP(cr)) {
                gmx_fatal(FARGS,
                          "Vsite recalculation with -rerun is not implemented with domain decomposition, use a single rank");
            }
            prepareRerunState(rerun_fr, state_global, constructVsites, vsite, top->idef, ir->delta_t, *fr, graph);
        }

        /* < 0 means stop at next step, > 0 means stop at next NS step */
        if ((signals[eglsSTOPCOND].set < 0) || (signals[eglsSTOPCOND].set > 0 ))
        {
            bLastStep = TRUE;
        }

        do_verbose = mdrunOptions.verbose;

        {
            gmx_bool bMasterState = TRUE;

            if (DOMAINDECOMP(cr))
            {
                /* Repartition the domain decomposition */
                dd_partition_system(fplog, step, cr,
                                    bMasterState, nstglobalcomm,
                                    state_global, top_global, ir,
                                    state, &f, mdAtoms, top, fr,
                                    vsite, constr,
                                    nrnb, wcycle,
                                    do_verbose);
                shouldCheckNumberOfBondedInteractions = true;
                update_realloc(upd, state->natoms);
            }
        }

        if (MASTER(cr))
        {
            print_ebin_header(fplog, step, t); /* can we improve the information printed here? */
        }

        if (ir->efep != efepNO)
        {
            update_mdatoms(mdatoms, state->lambda[efptMASS]);
        }

        clear_mat(force_vir);

        force_flags = (GMX_FORCE_STATECHANGED |
                       GMX_FORCE_DYNAMICBOX |
                       GMX_FORCE_ALLFORCES |
                       GMX_FORCE_VIRIAL |
                       GMX_FORCE_ENERGY |
                       (bDoFEP ? GMX_FORCE_DHDL : 0)
        );

        if (shellfc)
        {
            /* Now is the time to relax the shells */
            gmx_bool bNS = TRUE;
            relax_shell_flexcon(fplog, cr, ms, mdrunOptions.verbose, step,
                                ir, bNS, force_flags, top,
                                constr, enerd, fcd,
                                state, &f, force_vir, mdatoms,
                                nrnb, wcycle, graph, groups,
                                shellfc, fr, t, mu_tot,
                                vsite,
                                ddOpenBalanceRegion, ddCloseBalanceRegion);
        }
        else
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
                     GMX_FORCE_NS | force_flags,
                     ddOpenBalanceRegion, ddCloseBalanceRegion);
        }

        /* ########  END FIRST UPDATE STEP  ############## */
        /* ########  If doing VV, we now have v(dt) ###### */

        /* Now we have the energies and forces corresponding to the
         * coordinates at time t. We must output all of this before
         * the update.
         */
        {
            // dummy variables for rerun
            int nchkpt = 0;
            const gmx_bool bCPT = FALSE;
            const gmx_bool bRerunMD = TRUE;
            do_md_trajectory_writing(fplog, cr, nfile, fnm, step, step_rel, t,
                                     ir, state, state_global, observablesHistory,
                                     top_global, fr,
                                     outf, mdebin, ekind, f,
                                     &nchkpt,
                                     bCPT, bRerunMD, bLastStep,
                                     mdrunOptions.writeConfout,
                                     bSumEkinhOld);
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
        else if (MASTER(cr) &&
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

        /* #########   START SECOND UPDATE STEP ################# */

        if (graph)
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
            bool doInterSimSignal = FALSE;

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
                                state->box,
                                &totalNumberOfBondedInteractions, &bSumEkinhOld,
                                CGLO_GSTAT | CGLO_ENERGY
                                | CGLO_PRESSURE  // needed to compute virial
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

        if (ir->efep != efepNO)
        {
            /* Sum up the foreign energy and dhdl terms for md and sd.
               Currently done every step so that dhdl is correct in the .edr */
            sum_dhdl(enerd, state->lambda, ir->fepvals);
        }

        /* ################# END UPDATE STEP 2 ################# */
        /* #### We now have r(t+dt) and v(t+dt/2)  ############# */

        /* #########  BEGIN PREPARING EDR OUTPUT  ###########  */

        enerd->term[F_ETOT] = enerd->term[F_EPOT] + enerd->term[F_EKIN];

        /* #########  END PREPARING EDR OUTPUT  ###########  */

        /* Output stuff */
        if (MASTER(cr))
        {
            gmx_bool bCalcEnerStep = TRUE;
            upd_mdebin(mdebin, bDoDHDL, bCalcEnerStep,
                       t, mdatoms->tmass, enerd, state,
                       ir->fepvals, ir->expandedvals, state->box,
                       shake_vir, force_vir, total_vir, pres,
                       ekind, mu_tot, constr);

            gmx_bool do_ene = TRUE;
            gmx_bool do_log = TRUE;
            gmx_bool do_dr  = do_per_step(step, ir->nstdisreout);
            gmx_bool do_or  = do_per_step(step, ir->nstorireout);

            print_ebin(mdoutf_get_fp_ene(outf), do_ene, do_dr, do_or, do_log ? fplog : nullptr,
                       step, t,
                       eprNORMAL, mdebin, fcd, groups, &(ir->opts), ir->awh);

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
            (do_verbose || gmx_got_usr_signal()))
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
        gmx_bool bNeedRepartition = FALSE;
        if ((ir->eSwapCoords != eswapNO) && (step > 0) && !bLastStep &&
            do_per_step(step, ir->swap->nstswap))
        {
            bNeedRepartition = do_swapcoords(cr, step, t, ir, wcycle,
                                             rerun_fr.x,
                                             rerun_fr.box,
                                             MASTER(cr) && mdrunOptions.verbose,
                                             TRUE);

            if (bNeedRepartition && DOMAINDECOMP(cr))
            {
                dd_collect_state(cr->dd, state, state_global);
            }
        }

        if (bNeedRepartition && DOMAINDECOMP(cr))
        {
            dd_partition_system(fplog, step, cr, TRUE, 1,
                                state_global, top_global, ir,
                                state, &f, mdAtoms, top, fr,
                                vsite, constr,
                                nrnb, wcycle, FALSE);
            shouldCheckNumberOfBondedInteractions = true;
            update_realloc(upd, state->natoms);
        }

        if (MASTER(cr))
        {
            /* read next frame from input trajectory */
            bLastStep = !read_next_frame(oenv, status, &rerun_fr);
        }

        if (PAR(cr))
        {
            rerun_parallel_comm(cr, &rerun_fr, &bLastStep);
        }

        cycles = wallcycle_stop(wcycle, ewcSTEP);
        if (DOMAINDECOMP(cr) && wcycle)
        {
            dd_cycles_add(cr->dd, cycles, ddCyclStep);
        }

        if (!rerun_fr.bStep) {
            /* increase the MD step number */
            step++;
            step_rel++;
        }
    }
    /* End of main MD loop */

    /* Closing TNG files can include compressing data. Therefore it is good to do that
     * before stopping the time measurements. */
    mdoutf_tng_close(outf);

    /* Stop measuring walltime */
    walltime_accounting_end(walltime_accounting);

    if (MASTER(cr))
    {
        close_trx(status);
    }

    if (!thisRankHasDuty(cr, DUTY_PME))
    {
        /* Tell the PME only node to finish */
        gmx_pme_send_finish(cr);
    }

    done_ebin(mdebin->ebin);
    done_mdoutf(outf);

    done_shellfc(fplog, shellfc, step_rel);

    if (ir->bDoAwh)
    {
        delete ir->awh;
    }

    // Clean up swapcoords
    if (ir->eSwapCoords != eswapNO)
    {
        finish_swapcoords(ir->swap);
    }

    walltime_accounting_set_nsteps_done(walltime_accounting, step_rel);
    if (step_rel >= wcycle_get_reset_counters(wcycle) &&
        signals[eglsRESETCOUNTERS].set == 0)
    {
        walltime_accounting_set_valid_finish(walltime_accounting);
    }

    return 0;
}
