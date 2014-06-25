/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "smalloc.h"
#include "sysstuff.h"
#include "vec.h"
#include "statutil.h"
#include "vcm.h"
#include "mdebin.h"
#include "nrnb.h"
#include "calcmu.h"
#include "index.h"
#include "vsite.h"
#include "update.h"
#include "ns.h"
#include "trnio.h"
#include "xtcio.h"
#include "mdrun.h"
#include "md_support.h"
#include "md_logging.h"
#include "confio.h"
#include "network.h"
#include "pull.h"
#include "xvgr.h"
#include "physics.h"
#include "names.h"
#include "xmdrun.h"
#include "ionize.h"
#include "disre.h"
#include "orires.h"
#include "pme.h"
#include "mdatoms.h"
#include "repl_ex.h"
#include "qmmm.h"
#include "mpelogging.h"
#include "domdec.h"
#include "domdec_network.h"
#include "partdec.h"
#include "topsort.h"
#include "coulomb.h"
#include "constr.h"
#include "shellfc.h"
#include "compute_io.h"
#include "mvdata.h"
#include "checkpoint.h"
#include "mtop_util.h"
#include "sighandler.h"
#include "txtdump.h"
#include "string2.h"
#include "pme_loadbal.h"
#include "bondf.h"
#include "membed.h"
#include "types/nlistheuristics.h"
#include "types/iteratedconstraints.h"
#include "nbnxn_cuda_data_mgmt.h"

#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif
#ifdef GMX_THREAD_MPI
#include "tmpi.h"
#endif

#ifdef GMX_FAHCORE
#include "corewrap.h"
#endif

static void reset_all_counters(FILE *fplog, t_commrec *cr,
                               gmx_large_int_t step,
                               gmx_large_int_t *step_rel, t_inputrec *ir,
                               gmx_wallcycle_t wcycle, t_nrnb *nrnb,
                               gmx_runtime_t *runtime,
                               nbnxn_cuda_ptr_t cu_nbv)
{
    char sbuf[STEPSTRSIZE];

    /* Reset all the counters related to performance over the run */
    md_print_warn(cr, fplog, "step %s: resetting all time and cycle counters\n",
                  gmx_step_str(step, sbuf));

    if (cu_nbv)
    {
        nbnxn_cuda_reset_timings(cu_nbv);
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
    runtime_start(runtime);
    print_date_and_time(fplog, cr->nodeid, "Restarted time", runtime);
}

double do_md(FILE *fplog, t_commrec *cr, int nfile, const t_filenm fnm[],
             const output_env_t oenv, gmx_bool bVerbose, gmx_bool bCompact,
             int nstglobalcomm,
             gmx_vsite_t *vsite, gmx_constr_t constr,
             int stepout, t_inputrec *ir,
             gmx_mtop_t *top_global,
             t_fcdata *fcd,
             t_state *state_global,
             t_mdatoms *mdatoms,
             t_nrnb *nrnb, gmx_wallcycle_t wcycle,
             gmx_edsam_t ed, t_forcerec *fr,
             int repl_ex_nst, int repl_ex_nex, int repl_ex_seed, gmx_membed_t membed,
             real cpt_period, real max_hours,
             const char *deviceOptions,
             unsigned long Flags,
             gmx_runtime_t *runtime)
{
    gmx_mdoutf_t   *outf;
    gmx_large_int_t step, step_rel;
    double          run_time;
    double          t, t0, lam0[efptNR];
    gmx_bool        bGStatEveryStep, bGStat, bCalcVir, bCalcEner;
    gmx_bool        bNS, bNStList, bSimAnn, bStopCM, bRerunMD, bNotLastFrame = FALSE,
                    bFirstStep, bStateFromCP, bStateFromTPX, bInitStep, bLastStep,
                    bBornRadii, bStartingFromCpt;
    gmx_bool          bDoDHDL = FALSE, bDoFEP = FALSE, bDoExpanded = FALSE;
    gmx_bool          do_ene, do_log, do_verbose, bRerunWarnNoV = TRUE,
                      bForceUpdate = FALSE, bCPT;
    int               mdof_flags;
    gmx_bool          bMasterState;
    int               force_flags, cglo_flags;
    tensor            force_vir, shake_vir, total_vir, tmp_vir, pres;
    int               i, m;
    t_trxstatus      *status;
    rvec              mu_tot;
    t_vcm            *vcm;
    t_state          *bufstate = NULL;
    matrix           *scale_tot, pcoupl_mu, M, ebox;
    gmx_nlheur_t      nlh;
    t_trxframe        rerun_fr;
    gmx_repl_ex_t     repl_ex = NULL;
    int               nchkpt  = 1;
    gmx_localtop_t   *top;
    t_mdebin         *mdebin = NULL;
    t_state          *state    = NULL;
    rvec             *f_global = NULL;
    int               n_xtc    = -1;
    rvec             *x_xtc    = NULL;
    gmx_enerdata_t   *enerd;
    rvec             *f = NULL;
    gmx_global_stat_t gstat;
    gmx_update_t      upd   = NULL;
    t_graph          *graph = NULL;
    globsig_t         gs;
    gmx_rng_t         mcrng = NULL;
    gmx_bool          bFFscan;
    gmx_groups_t     *groups;
    gmx_ekindata_t   *ekind, *ekind_save;
    gmx_shellfc_t     shellfc;
    int               count, nconverged = 0;
    real              timestep = 0;
    double            tcount   = 0;
    gmx_bool          bIonize  = FALSE;
    gmx_bool          bTCR     = FALSE, bConverged = TRUE, bOK, bSumEkinhOld, bExchanged;
    gmx_bool          bAppend;
    gmx_bool          bResetCountersHalfMaxH = FALSE;
    gmx_bool          bVV, bIterativeCase, bFirstIterate, bTemp, bPres, bTrotter;
    gmx_bool          bUpdateDoLR;
    real              mu_aver = 0, dvdl_constr;
    int               a0, a1, gnx = 0, ii;
    atom_id          *grpindex = NULL;
    char             *grpname;
    t_coupl_rec      *tcr     = NULL;
    rvec             *xcopy   = NULL, *vcopy = NULL, *cbuf = NULL;
    matrix            boxcopy = {{0}}, lastbox;
    tensor            tmpvir;
    real              fom, oldfom, veta_save, pcurr, scalevir, tracevir;
    real              vetanew = 0;
    int               lamnew  = 0;
    /* for FEP */
    int               nstfep;
    real              rate;
    double            cycles;
    real              saved_conserved_quantity = 0;
    real              last_ekin                = 0;
    int               iter_i;
    t_extmass         MassQ;
    int             **trotter_seq;
    char              sbuf[STEPSTRSIZE], sbuf2[STEPSTRSIZE];
    int               handled_stop_condition = gmx_stop_cond_none; /* compare to get_stop_condition*/
    gmx_iterate_t     iterate;
    gmx_large_int_t   multisim_nsteps = -1;                        /* number of steps to do  before first multisim
                                                                      simulation stops. If equal to zero, don't
                                                                      communicate any more between multisims.*/
    /* PME load balancing data for GPU kernels */
    pme_load_balancing_t pme_loadbal = NULL;
    double               cycles_pmes;
    gmx_bool             bPMETuneTry = FALSE, bPMETuneRunning = FALSE;

#ifdef GMX_FAHCORE
    /* Temporary addition for FAHCORE checkpointing */
    int chkpt_ret;
#endif

    /* Check for special mdrun options */
    bRerunMD = (Flags & MD_RERUN);
    bIonize  = (Flags & MD_IONIZE);
    bFFscan  = (Flags & MD_FFSCAN);
    bAppend  = (Flags & MD_APPENDFILES);
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
       md-vv-avek uses averaged half step velocities for T-control (but full step ekin for P control)
       md uses averaged half step kinetic energies to determine temperature unless defined otherwise by GMX_EKIN_AVE_VEL; */
    bVV = EI_VV(ir->eI);
    if (bVV) /* to store the initial velocities while computing virial */
    {
        snew(cbuf, top_global->natoms);
    }
    /* all the iteratative cases - only if there are constraints */
    bIterativeCase = ((IR_NPH_TROTTER(ir) || IR_NPT_TROTTER(ir)) && (constr) && (!bRerunMD));
    gmx_iterate_init(&iterate, FALSE); /* The default value of iterate->bIterationActive is set to
                                          false in this step.  The correct value, true or false,
                                          is set at each step, as it depends on the frequency of temperature
                                          and pressure control.*/
    bTrotter = (bVV && (IR_NPT_TROTTER(ir) || IR_NPH_TROTTER(ir) || IR_NVT_TROTTER(ir)));

    if (bRerunMD)
    {
        /* Since we don't know if the frames read are related in any way,
         * rebuild the neighborlist at every step.
         */
        ir->nstlist       = 1;
        ir->nstcalcenergy = 1;
        nstglobalcomm     = 1;
    }

    check_ir_old_tpx_versions(cr, fplog, ir, top_global);

    nstglobalcomm   = check_nstglobalcomm(fplog, cr, nstglobalcomm, ir);
    bGStatEveryStep = (nstglobalcomm == 1);

    if (!bGStatEveryStep && ir->nstlist == -1 && fplog != NULL)
    {
        fprintf(fplog,
                "To reduce the energy communication with nstlist = -1\n"
                "the neighbor list validity should not be checked at every step,\n"
                "this means that exact integration is not guaranteed.\n"
                "The neighbor list validity is checked after:\n"
                "  <n.list life time> - 2*std.dev.(n.list life time)  steps.\n"
                "In most cases this will result in exact integration.\n"
                "This reduces the energy communication by a factor of 2 to 3.\n"
                "If you want less energy communication, set nstlist > 3.\n\n");
    }

    if (bRerunMD || bFFscan)
    {
        ir->nstxtcout = 0;
    }
    groups = &top_global->groups;

    /* Initial values */
    init_md(fplog, cr, ir, oenv, &t, &t0, state_global->lambda,
            &(state_global->fep_state), lam0,
            nrnb, top_global, &upd,
            nfile, fnm, &outf, &mdebin,
            force_vir, shake_vir, mu_tot, &bSimAnn, &vcm, state_global, Flags);

    clear_mat(total_vir);
    clear_mat(pres);
    /* Energy terms and groups */
    snew(enerd, 1);
    init_enerdata(top_global->groups.grps[egcENER].nr, ir->fepvals->n_lambda,
                  enerd);
    if (DOMAINDECOMP(cr))
    {
        f = NULL;
    }
    else
    {
        snew(f, top_global->natoms);
    }

    /* Kinetic energy data */
    snew(ekind, 1);
    init_ekindata(fplog, top_global, &(ir->opts), ekind);
    /* needed for iteration of constraints */
    snew(ekind_save, 1);
    init_ekindata(fplog, top_global, &(ir->opts), ekind_save);
    /* Copy the cos acceleration to the groups struct */
    ekind->cosacc.cos_accel = ir->cos_accel;

    gstat = global_stat_init(ir);
    debug_gmx();

    /* Check for polarizable models and flexible constraints */
    shellfc = init_shell_flexcon(fplog,
                                 top_global, n_flexible_constraints(constr),
                                 (ir->bContinuation ||
                                  (DOMAINDECOMP(cr) && !MASTER(cr))) ?
                                 NULL : state_global->x);
    if (shellfc && ir->nstcalcenergy != 1)
    {
        gmx_fatal(FARGS, "You have nstcalcenergy set to a value (%d) that is different from 1.\nThis is not supported in combinations with shell particles.\nPlease make a new tpr file.", ir->nstcalcenergy);
    }
    if (shellfc && DOMAINDECOMP(cr))
    {
        gmx_fatal(FARGS, "In order to run parallel simulations with shells you need to use the -pd flag to mdrun.");
    }
    if (shellfc && ir->eI == eiNM)
    {
        /* Currently shells don't work with Normal Modes */
        gmx_fatal(FARGS, "Normal Mode analysis is not supported with shells.\nIf you'd like to help with adding support, we have an open discussion at http://redmine.gromacs.org/issues/879\n");
    }

    if (vsite && ir->eI == eiNM)
    {
        /* Currently virtual sites don't work with Normal Modes */
        gmx_fatal(FARGS, "Normal Mode analysis is not supported with virtual sites.\nIf you'd like to help with adding support, we have an open discussion at http://redmine.gromacs.org/issues/879\n");
    }

    if (DEFORM(*ir))
    {
#ifdef GMX_THREAD_MPI
        tMPI_Thread_mutex_lock(&deform_init_box_mutex);
#endif
        set_deform_reference_box(upd,
                                 deform_init_init_step_tpx,
                                 deform_init_box_tpx);
#ifdef GMX_THREAD_MPI
        tMPI_Thread_mutex_unlock(&deform_init_box_mutex);
#endif
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

        snew(state, 1);
        dd_init_local_state(cr->dd, state_global, state);

        if (DDMASTER(cr->dd) && ir->nstfout)
        {
            snew(f_global, state_global->natoms);
        }
    }
    else
    {
        if (PAR(cr))
        {
            /* Initialize the particle decomposition and split the topology */
            top = split_system(fplog, top_global, ir, cr);

            pd_cg_range(cr, &fr->cg0, &fr->hcg);
            pd_at_range(cr, &a0, &a1);
        }
        else
        {
            top = gmx_mtop_generate_local_top(top_global, ir);

            a0 = 0;
            a1 = top_global->natoms;
        }

        forcerec_set_excl_load(fr, top, cr);

        state    = partdec_init_local_state(cr, state_global);
        f_global = f;

        atoms2md(top_global, ir, 0, NULL, a0, a1-a0, mdatoms);

        if (vsite)
        {
            set_vsite_top(vsite, top, mdatoms, cr);
        }

        if (ir->ePBC != epbcNONE && !fr->bMolPBC)
        {
            graph = mk_graph(fplog, &(top->idef), 0, top_global->natoms, FALSE, FALSE);
        }

        if (shellfc)
        {
            make_local_shells(cr, mdatoms, shellfc);
        }

        setup_bonded_threading(fr, &top->idef);

        if (ir->pull && PAR(cr))
        {
            dd_make_local_pull_groups(NULL, ir->pull, mdatoms);
        }
    }

    if (DOMAINDECOMP(cr))
    {
        /* Distribute the charge groups over the nodes from the master node */
        dd_partition_system(fplog, ir->init_step, cr, TRUE, 1,
                            state_global, top_global, ir,
                            state, &f, mdatoms, top, fr,
                            vsite, shellfc, constr,
                            nrnb, wcycle, FALSE);

    }

    update_mdatoms(mdatoms, state->lambda[efptMASS]);

    if (opt2bSet("-cpi", nfile, fnm))
    {
        bStateFromCP = gmx_fexist_master(opt2fn_master("-cpi", nfile, fnm, cr), cr);
    }
    else
    {
        bStateFromCP = FALSE;
    }

    if (ir->bExpanded)
    {
        init_expanded_ensemble(bStateFromCP,ir,&mcrng,&state->dfhist);
    }

    if (MASTER(cr))
    {
        if (bStateFromCP)
        {
            /* Update mdebin with energy history if appending to output files */
            if (Flags & MD_APPENDFILES)
            {
                restore_energyhistory_from_state(mdebin, &state_global->enerhist);
            }
            else
            {
                /* We might have read an energy history from checkpoint,
                 * free the allocated memory and reset the counts.
                 */
                done_energyhistory(&state_global->enerhist);
                init_energyhistory(&state_global->enerhist);
            }
        }
        /* Set the initial energy history in state by updating once */
        update_energyhistory(&state_global->enerhist, mdebin);
    }

    if ((state->flags & (1<<estLD_RNG)) && (Flags & MD_READ_RNG))
    {
        /* Set the random state if we read a checkpoint file */
        set_stochd_state(upd, state);
    }

    if (state->flags & (1<<estMC_RNG))
    {
        set_mc_state(mcrng, state);
    }

    /* Initialize constraints */
    if (constr)
    {
        if (!DOMAINDECOMP(cr))
        {
            set_constraints(constr, top, ir, mdatoms, cr);
        }
    }

    /* Check whether we have to GCT stuff */
    bTCR = ftp2bSet(efGCT, nfile, fnm);
    if (bTCR)
    {
        if (MASTER(cr))
        {
            fprintf(stderr, "Will do General Coupling Theory!\n");
        }
        gnx = top_global->mols.nr;
        snew(grpindex, gnx);
        for (i = 0; (i < gnx); i++)
        {
            grpindex[i] = i;
        }
    }

    if (repl_ex_nst > 0)
    {
        /* We need to be sure replica exchange can only occur
         * when the energies are current */
        check_nst_param(fplog, cr, "nstcalcenergy", ir->nstcalcenergy,
                        "repl_ex_nst", &repl_ex_nst);
        /* This check needs to happen before inter-simulation
         * signals are initialized, too */
    }
    if (repl_ex_nst > 0 && MASTER(cr))
    {
        repl_ex = init_replica_exchange(fplog, cr->ms, state_global, ir,
                                        repl_ex_nst, repl_ex_nex, repl_ex_seed);
    }

    /* PME tuning is only supported with GPUs or PME nodes and not with rerun.
     * With perturbed charges with soft-core we should not change the cut-off.
     */
    if ((Flags & MD_TUNEPME) &&
        EEL_PME(fr->eeltype) &&
        ( (fr->cutoff_scheme == ecutsVERLET && fr->nbv->bUseGPU) || !(cr->duty & DUTY_PME)) &&
        !(ir->efep != efepNO && mdatoms->nChargePerturbed > 0 && ir->fepvals->bScCoul) &&
        !bRerunMD)
    {
        pme_loadbal_init(&pme_loadbal, ir, state->box, fr->ic, fr->pmedata);
        cycles_pmes = 0;
        if (cr->duty & DUTY_PME)
        {
            /* Start tuning right away, as we can't measure the load */
            bPMETuneRunning = TRUE;
        }
        else
        {
            /* Separate PME nodes, we can measure the PP/PME load balance */
            bPMETuneTry = TRUE;
        }
    }

    if (!ir->bContinuation && !bRerunMD)
    {
        if (mdatoms->cFREEZE && (state->flags & (1<<estV)))
        {
            /* Set the velocities of frozen particles to zero */
            for (i = mdatoms->start; i < mdatoms->start+mdatoms->homenr; i++)
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

        if (constr)
        {
            /* Constrain the initial coordinates and velocities */
            do_constrain_first(fplog, constr, ir, mdatoms, state, f,
                               graph, cr, nrnb, fr, top, shake_vir);
        }
        if (vsite)
        {
            /* Construct the virtual sites for the initial configuration */
            construct_vsites(fplog, vsite, state->x, nrnb, ir->delta_t, NULL,
                             top->idef.iparams, top->idef.il,
                             fr->ePBC, fr->bMolPBC, graph, cr, state->box);
        }
    }

    debug_gmx();

    /* set free energy calculation frequency as the minimum
       greatest common denominator of nstdhdl, nstexpanded, and repl_ex_nst*/
    nstfep = ir->fepvals->nstdhdl;
    if (ir->bExpanded)
    {
        nstfep = gmx_greatest_common_divisor(ir->fepvals->nstdhdl,nstfep);
    }
    if (repl_ex_nst > 0)
    {
        nstfep = gmx_greatest_common_divisor(repl_ex_nst,nstfep);
    }

    /* I'm assuming we need global communication the first time! MRS */
    cglo_flags = (CGLO_TEMPERATURE | CGLO_GSTAT
                  | ((ir->comm_mode != ecmNO) ? CGLO_STOPCM : 0)
                  | (bVV ? CGLO_PRESSURE : 0)
                  | (bVV ? CGLO_CONSTRAINT : 0)
                  | (bRerunMD ? CGLO_RERUNMD : 0)
                  | ((Flags & MD_READ_EKIN) ? CGLO_READEKIN : 0));

    bSumEkinhOld = FALSE;
    compute_globals(fplog, gstat, cr, ir, fr, ekind, state, state_global, mdatoms, nrnb, vcm,
                    NULL, enerd, force_vir, shake_vir, total_vir, pres, mu_tot,
                    constr, NULL, FALSE, state->box,
                    top_global, &pcurr, top_global->natoms, &bSumEkinhOld, cglo_flags);
    if (ir->eI == eiVVAK)
    {
        /* a second call to get the half step temperature initialized as well */
        /* we do the same call as above, but turn the pressure off -- internally to
           compute_globals, this is recognized as a velocity verlet half-step
           kinetic energy calculation.  This minimized excess variables, but
           perhaps loses some logic?*/

        compute_globals(fplog, gstat, cr, ir, fr, ekind, state, state_global, mdatoms, nrnb, vcm,
                        NULL, enerd, force_vir, shake_vir, total_vir, pres, mu_tot,
                        constr, NULL, FALSE, state->box,
                        top_global, &pcurr, top_global->natoms, &bSumEkinhOld,
                        cglo_flags &~(CGLO_STOPCM | CGLO_PRESSURE));
    }

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

    /* if using an iterative algorithm, we need to create a working directory for the state. */
    if (bIterativeCase)
    {
        bufstate = init_bufstate(state);
    }
    if (bFFscan)
    {
        snew(xcopy, state->natoms);
        snew(vcopy, state->natoms);
        copy_rvecn(state->x, xcopy, 0, state->natoms);
        copy_rvecn(state->v, vcopy, 0, state->natoms);
        copy_mat(state->box, boxcopy);
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
                    constr_rmsd(constr, FALSE));
        }
        if (EI_STATE_VELOCITY(ir->eI))
        {
            fprintf(fplog, "Initial temperature: %g K\n", enerd->term[F_TEMP]);
        }
        if (bRerunMD)
        {
            fprintf(stderr, "starting md rerun '%s', reading coordinates from"
                    " input trajectory '%s'\n\n",
                    *(top_global->name), opt2fn("-rerun", nfile, fnm));
            if (bVerbose)
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

    print_start(fplog, cr, runtime, "mdrun");
    runtime_start(runtime);
    wallcycle_start(wcycle, ewcRUN);

    /* safest point to do file checkpointing is here.  More general point would be immediately before integrator call */
#ifdef GMX_FAHCORE
    chkpt_ret = fcCheckPointParallel( cr->nodeid,
                                      NULL, 0);
    if (chkpt_ret == 0)
    {
        gmx_fatal( 3, __FILE__, __LINE__, "Checkpoint error on step %d\n", 0 );
    }
#endif

    debug_gmx();
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
            bNotLastFrame = read_first_frame(oenv, &status,
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
                if (max_cutoff2(ir->ePBC, rerun_fr.box) < sqr(fr->rlistlong))
                {
                    gmx_fatal(FARGS, "Rerun trajectory frame step %d time %f has too small box dimensions", rerun_fr.step, rerun_fr.time);
                }
            }
        }

        if (PAR(cr))
        {
            rerun_parallel_comm(cr, &rerun_fr, &bNotLastFrame);
        }

        if (ir->ePBC != epbcNONE)
        {
            /* Set the shift vectors.
             * Necessary here when have a static box different from the tpr box.
             */
            calc_shifts(rerun_fr.box, fr->shift_vec);
        }
    }

    /* loop over MD steps or if rerunMD to end of input trajectory */
    bFirstStep = TRUE;
    /* Skip the first Nose-Hoover integration when we get the state from tpx */
    bStateFromTPX    = !bStateFromCP;
    bInitStep        = bFirstStep && (bStateFromTPX || bVV);
    bStartingFromCpt = (Flags & MD_STARTFROMCPT) && bInitStep;
    bLastStep        = FALSE;
    bSumEkinhOld     = FALSE;
    bExchanged       = FALSE;

    init_global_signals(&gs, cr, ir, repl_ex_nst);

    step     = ir->init_step;
    step_rel = 0;

    if (ir->nstlist == -1)
    {
        init_nlistheuristics(&nlh, bGStatEveryStep, step);
    }

    if (MULTISIM(cr) && (repl_ex_nst <= 0 ))
    {
        /* check how many steps are left in other sims */
        multisim_nsteps = get_multisim_nsteps(cr, ir->nsteps);
    }


    /* and stop now if we should */
    bLastStep = (bRerunMD || (ir->nsteps >= 0 && step_rel > ir->nsteps) ||
                 ((multisim_nsteps >= 0) && (step_rel >= multisim_nsteps )));
    while (!bLastStep || (bRerunMD && bNotLastFrame))
    {

        wallcycle_start(wcycle, ewcSTEP);

        GMX_MPE_LOG(ev_timestep1);

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

        if (ir->efep != efepNO || ir->bSimTemp)
        {
            /* find and set the current lambdas.  If rerunning, we either read in a state, or a lambda value,
               requiring different logic. */

            set_current_lambdas(step, ir->fepvals, bRerunMD, &rerun_fr, state_global, state, lam0);
            bDoDHDL      = do_per_step(step, ir->fepvals->nstdhdl);
            bDoFEP       = (do_per_step(step, nstfep) && (ir->efep != efepNO));
            bDoExpanded  = (do_per_step(step, ir->expandedvals->nstexpanded)
                            && (ir->bExpanded) && (step > 0) && (!bStartingFromCpt));
        }

        if (bSimAnn)
        {
            update_annealing_target_temp(&(ir->opts), t);
        }

        if (bRerunMD)
        {
            if (!(DOMAINDECOMP(cr) && !MASTER(cr)))
            {
                for (i = 0; i < state_global->natoms; i++)
                {
                    copy_rvec(rerun_fr.x[i], state_global->x[i]);
                }
                if (rerun_fr.bV)
                {
                    for (i = 0; i < state_global->natoms; i++)
                    {
                        copy_rvec(rerun_fr.v[i], state_global->v[i]);
                    }
                }
                else
                {
                    for (i = 0; i < state_global->natoms; i++)
                    {
                        clear_rvec(state_global->v[i]);
                    }
                    if (bRerunWarnNoV)
                    {
                        fprintf(stderr, "\nWARNING: Some frames do not contain velocities.\n"
                                "         Ekin, temperature and pressure are incorrect,\n"
                                "         the virial will be incorrect when constraints are present.\n"
                                "\n");
                        bRerunWarnNoV = FALSE;
                    }
                }
            }
            copy_mat(rerun_fr.box, state_global->box);
            copy_mat(state_global->box, state->box);

            if (vsite && (Flags & MD_RERUN_VSITE))
            {
                if (DOMAINDECOMP(cr))
                {
                    gmx_fatal(FARGS, "Vsite recalculation with -rerun is not implemented for domain decomposition, use particle decomposition");
                }
                if (graph)
                {
                    /* Following is necessary because the graph may get out of sync
                     * with the coordinates if we only have every N'th coordinate set
                     */
                    mk_mshift(fplog, graph, fr->ePBC, state->box, state->x);
                    shift_self(graph, state->box, state->x);
                }
                construct_vsites(fplog, vsite, state->x, nrnb, ir->delta_t, state->v,
                                 top->idef.iparams, top->idef.il,
                                 fr->ePBC, fr->bMolPBC, graph, cr, state->box);
                if (graph)
                {
                    unshift_self(graph, state->box, state->x);
                }
            }
        }

        /* Stop Center of Mass motion */
        bStopCM = (ir->comm_mode != ecmNO && do_per_step(step, ir->nstcomm));

        /* Copy back starting coordinates in case we're doing a forcefield scan */
        if (bFFscan)
        {
            for (ii = 0; (ii < state->natoms); ii++)
            {
                copy_rvec(xcopy[ii], state->x[ii]);
                copy_rvec(vcopy[ii], state->v[ii]);
            }
            copy_mat(boxcopy, state->box);
        }

        if (bRerunMD)
        {
            /* for rerun MD always do Neighbour Searching */
            bNS      = (bFirstStep || ir->nstlist != 0);
            bNStList = bNS;
        }
        else
        {
            /* Determine whether or not to do Neighbour Searching and LR */
            bNStList = (ir->nstlist > 0  && step % ir->nstlist == 0);

            bNS = (bFirstStep || bExchanged || bNStList || bDoFEP ||
                   (ir->nstlist == -1 && nlh.nabnsb > 0));

            if (bNS && ir->nstlist == -1)
            {
                set_nlistheuristics(&nlh, bFirstStep || bExchanged || bDoFEP, step);
            }
        }

        /* check whether we should stop because another simulation has
           stopped. */
        if (MULTISIM(cr))
        {
            if ( (multisim_nsteps >= 0) &&  (step_rel >= multisim_nsteps)  &&
                 (multisim_nsteps != ir->nsteps) )
            {
                if (bNS)
                {
                    if (MASTER(cr))
                    {
                        fprintf(stderr,
                                "Stopping simulation %d because another one has finished\n",
                                cr->ms->sim);
                    }
                    bLastStep         = TRUE;
                    gs.sig[eglsCHKPT] = 1;
                }
            }
        }

        /* < 0 means stop at next step, > 0 means stop at next NS step */
        if ( (gs.set[eglsSTOPCOND] < 0) ||
             ( (gs.set[eglsSTOPCOND] > 0) && (bNStList || ir->nstlist == 0) ) )
        {
            bLastStep = TRUE;
        }

        /* Determine whether or not to update the Born radii if doing GB */
        bBornRadii = bFirstStep;
        if (ir->implicit_solvent && (step % ir->nstgbradii == 0))
        {
            bBornRadii = TRUE;
        }

        do_log     = do_per_step(step, ir->nstlog) || bFirstStep || bLastStep;
        do_verbose = bVerbose &&
            (step % stepout == 0 || bFirstStep || bLastStep);

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
                if (DYNAMIC_BOX(*ir))
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
                wallcycle_start(wcycle, ewcDOMDEC);
                dd_partition_system(fplog, step, cr,
                                    bMasterState, nstglobalcomm,
                                    state_global, top_global, ir,
                                    state, &f, mdatoms, top, fr,
                                    vsite, shellfc, constr,
                                    nrnb, wcycle,
                                    do_verbose && !bPMETuneRunning);
                wallcycle_stop(wcycle, ewcDOMDEC);
                /* If using an iterative integrator, reallocate space to match the decomposition */
            }
        }

        if (MASTER(cr) && do_log && !bFFscan)
        {
            print_ebin_header(fplog, step, t, state->lambda[efptFEP]); /* can we improve the information printed here? */
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
            compute_globals(fplog, gstat, cr, ir, fr, ekind, state, state_global, mdatoms, nrnb, vcm,
                            wcycle, enerd, NULL, NULL, NULL, NULL, mu_tot,
                            constr, NULL, FALSE, state->box,
                            top_global, &pcurr, top_global->natoms, &bSumEkinhOld,
                            CGLO_RERUNMD | CGLO_GSTAT | CGLO_TEMPERATURE);
        }
        clear_mat(force_vir);

        /* Ionize the atoms if necessary */
        if (bIonize)
        {
            ionize(fplog, oenv, mdatoms, top_global, t, ir, state->x, state->v,
                   mdatoms->start, mdatoms->start+mdatoms->homenr, state->box, cr);
        }

        /* Update force field in ffscan program */
        if (bFFscan)
        {
            if (update_forcefield(fplog,
                                  nfile, fnm, fr,
                                  mdatoms->nr, state->x, state->box))
            {
                gmx_finalize_par();

                exit(0);
            }
        }

        GMX_MPE_LOG(ev_timestep2);

        /* We write a checkpoint at this MD step when:
         * either at an NS step when we signalled through gs,
         * or at the last step (but not when we do not want confout),
         * but never at the first step or with rerun.
         */
        bCPT = (((gs.set[eglsCHKPT] && (bNS || ir->nstlist == 0)) ||
                 (bLastStep && (Flags & MD_CONFOUT))) &&
                step > ir->init_step && !bRerunMD);
        if (bCPT)
        {
            gs.set[eglsCHKPT] = 0;
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

            bCalcEner = do_per_step(step-1, ir->nstcalcenergy);
            bCalcVir  = bCalcEner ||
                (ir->epc != epcNO && (do_per_step(step, ir->nstpcouple) || do_per_step(step-1, ir->nstpcouple)));
        }
        else
        {
            bCalcEner = do_per_step(step, ir->nstcalcenergy);
            bCalcVir  = bCalcEner ||
                (ir->epc != epcNO && do_per_step(step, ir->nstpcouple));
        }

        /* Do we need global communication ? */
        bGStat = (bCalcVir || bCalcEner || bStopCM ||
                  do_per_step(step, nstglobalcomm) || (bVV && IR_NVT_TROTTER(ir) && do_per_step(step-1, nstglobalcomm)) ||
                  (ir->nstlist == -1 && !bRerunMD && step >= nlh.step_nscheck));

        do_ene = (do_per_step(step, ir->nstenergy) || bLastStep);

        if (do_ene || do_log)
        {
            bCalcVir  = TRUE;
            bCalcEner = TRUE;
            bGStat    = TRUE;
        }

        /* these CGLO_ options remain the same throughout the iteration */
        cglo_flags = ((bRerunMD ? CGLO_RERUNMD : 0) |
                      (bGStat ? CGLO_GSTAT : 0)
                      );

        force_flags = (GMX_FORCE_STATECHANGED |
                       ((DYNAMIC_BOX(*ir) || bRerunMD) ? GMX_FORCE_DYNAMICBOX : 0) |
                       GMX_FORCE_ALLFORCES |
                       GMX_FORCE_SEPLRF |
                       (bCalcVir ? GMX_FORCE_VIRIAL : 0) |
                       (bCalcEner ? GMX_FORCE_ENERGY : 0) |
                       (bDoFEP ? GMX_FORCE_DHDL : 0)
                       );

        if (fr->bTwinRange)
        {
            if (do_per_step(step, ir->nstcalclr))
            {
                force_flags |= GMX_FORCE_DO_LR;
            }
        }

        if (shellfc)
        {
            /* Now is the time to relax the shells */
            count = relax_shell_flexcon(fplog, cr, bVerbose, bFFscan ? step+1 : step,
                                        ir, bNS, force_flags,
                                        bStopCM, top, top_global,
                                        constr, enerd, fcd,
                                        state, f, force_vir, mdatoms,
                                        nrnb, wcycle, graph, groups,
                                        shellfc, fr, bBornRadii, t, mu_tot,
                                        state->natoms, &bConverged, vsite,
                                        outf->fp_field);
            tcount += count;

            if (bConverged)
            {
                nconverged++;
            }
        }
        else
        {
            /* The coordinates (x) are shifted (to get whole molecules)
             * in do_force.
             * This is parallellized as well, and does communication too.
             * Check comments in sim_util.c
             */
            do_force(fplog, cr, ir, step, nrnb, wcycle, top, top_global, groups,
                     state->box, state->x, &state->hist,
                     f, force_vir, mdatoms, enerd, fcd,
                     state->lambda, graph,
                     fr, vsite, mu_tot, t, outf->fp_field, ed, bBornRadii,
                     (bNS ? GMX_FORCE_NS : 0) | force_flags);
        }

        GMX_BARRIER(cr->mpi_comm_mygroup);

        if (bTCR)
        {
            mu_aver = calc_mu_aver(cr, state->x, mdatoms->chargeA,
                                   mu_tot, &top_global->mols, mdatoms, gnx, grpindex);
        }

        if (bTCR && bFirstStep)
        {
            tcr = init_coupling(fplog, nfile, fnm, cr, fr, mdatoms, &(top->idef));
            fprintf(fplog, "Done init_coupling\n");
            fflush(fplog);
        }

        if (bVV && !bStartingFromCpt && !bRerunMD)
        /*  ############### START FIRST UPDATE HALF-STEP FOR VV METHODS############### */
        {
            if (ir->eI == eiVV && bInitStep)
            {
                /* if using velocity verlet with full time step Ekin,
                 * take the first half step only to compute the
                 * virial for the first step. From there,
                 * revert back to the initial coordinates
                 * so that the input is actually the initial step.
                 */
                copy_rvecn(state->v, cbuf, 0, state->natoms); /* should make this better for parallelizing? */
            }
            else
            {
                /* this is for NHC in the Ekin(t+dt/2) version of vv */
                trotter_update(ir, step, ekind, enerd, state, total_vir, mdatoms, &MassQ, trotter_seq, ettTSEQ1);
            }

            /* If we are using twin-range interactions where the long-range component
             * is only evaluated every nstcalclr>1 steps, we should do a special update
             * step to combine the long-range forces on these steps.
             * For nstcalclr=1 this is not done, since the forces would have been added
             * directly to the short-range forces already.
             *
             * TODO Remove various aspects of VV+twin-range in master
             * branch, because VV integrators did not ever support
             * twin-range multiple time stepping with constraints.
             */
            bUpdateDoLR = (fr->bTwinRange && do_per_step(step, ir->nstcalclr));

            update_coords(fplog, step, ir, mdatoms, state, fr->bMolPBC,
                          f, bUpdateDoLR, fr->f_twin, bCalcVir ? &fr->vir_twin_constr : NULL, fcd,
                          ekind, M, wcycle, upd, bInitStep, etrtVELOCITY1,
                          cr, nrnb, constr, &top->idef);

            if (bIterativeCase && do_per_step(step-1, ir->nstpcouple) && !bInitStep)
            {
                gmx_iterate_init(&iterate, TRUE);
            }
            /* for iterations, we save these vectors, as we will be self-consistently iterating
               the calculations */

            /*#### UPDATE EXTENDED VARIABLES IN TROTTER FORMULATION */

            /* save the state */
            if (iterate.bIterationActive)
            {
                copy_coupling_state(state, bufstate, ekind, ekind_save, &(ir->opts));
            }

            bFirstIterate = TRUE;
            while (bFirstIterate || iterate.bIterationActive)
            {
                if (iterate.bIterationActive)
                {
                    copy_coupling_state(bufstate, state, ekind_save, ekind, &(ir->opts));
                    if (bFirstIterate && bTrotter)
                    {
                        /* The first time through, we need a decent first estimate
                           of veta(t+dt) to compute the constraints.  Do
                           this by computing the box volume part of the
                           trotter integration at this time. Nothing else
                           should be changed by this routine here.  If
                           !(first time), we start with the previous value
                           of veta.  */

                        veta_save = state->veta;
                        trotter_update(ir, step, ekind, enerd, state, total_vir, mdatoms, &MassQ, trotter_seq, ettTSEQ0);
                        vetanew     = state->veta;
                        state->veta = veta_save;
                    }
                }

                bOK = TRUE;
                if (!bRerunMD || rerun_fr.bV || bForceUpdate)     /* Why is rerun_fr.bV here?  Unclear. */
                {
                    update_constraints(fplog, step, NULL, ir, ekind, mdatoms,
                                       state, fr->bMolPBC, graph, f,
                                       &top->idef, shake_vir, NULL,
                                       cr, nrnb, wcycle, upd, constr,
                                       bInitStep, TRUE, bCalcVir, vetanew);

                    if (bCalcVir && bUpdateDoLR && ir->nstcalclr > 1)
                    {
                        /* Correct the virial for multiple time stepping */
                        m_sub(shake_vir, fr->vir_twin_constr, shake_vir);
                    }

                    if (!bOK && !bFFscan)
                    {
                        gmx_fatal(FARGS, "Constraint error: Shake, Lincs or Settle could not solve the constrains");
                    }

                }
                else if (graph)
                {
                    /* Need to unshift here if a do_force has been
                       called in the previous step */
                    unshift_self(graph, state->box, state->x);
                }

                /* if VV, compute the pressure and constraints */
                /* For VV2, we strictly only need this if using pressure
                 * control, but we really would like to have accurate pressures
                 * printed out.
                 * Think about ways around this in the future?
                 * For now, keep this choice in comments.
                 */
                /*bPres = (ir->eI==eiVV || IR_NPT_TROTTER(ir)); */
                /*bTemp = ((ir->eI==eiVV &&(!bInitStep)) || (ir->eI==eiVVAK && IR_NPT_TROTTER(ir)));*/
                bPres = TRUE;
                bTemp = ((ir->eI == eiVV && (!bInitStep)) || (ir->eI == eiVVAK));
                if (bCalcEner && ir->eI == eiVVAK)  /*MRS:  7/9/2010 -- this still doesn't fix it?*/
                {
                    bSumEkinhOld = TRUE;
                }
                /* for vv, the first half of the integration actually corresponds to the previous step.
                   So we need information from the last step in the first half of the integration */
                if (bGStat || do_per_step(step-1, nstglobalcomm))
                {
                    compute_globals(fplog, gstat, cr, ir, fr, ekind, state, state_global, mdatoms, nrnb, vcm,
                                    wcycle, enerd, force_vir, shake_vir, total_vir, pres, mu_tot,
                                    constr, NULL, FALSE, state->box,
                                    top_global, &pcurr, top_global->natoms, &bSumEkinhOld,
                                    cglo_flags
                                    | CGLO_ENERGY
                                    | (bTemp ? CGLO_TEMPERATURE : 0)
                                    | (bPres ? CGLO_PRESSURE : 0)
                                    | (bPres ? CGLO_CONSTRAINT : 0)
                                    | ((iterate.bIterationActive) ? CGLO_ITERATE : 0)
                                    | (bFirstIterate ? CGLO_FIRSTITERATE : 0)
                                    | CGLO_SCALEEKIN
                                    );
                    /* explanation of above:
                       a) We compute Ekin at the full time step
                       if 1) we are using the AveVel Ekin, and it's not the
                       initial step, or 2) if we are using AveEkin, but need the full
                       time step kinetic energy for the pressure (always true now, since we want accurate statistics).
                       b) If we are using EkinAveEkin for the kinetic energy for the temperature control, we still feed in
                       EkinAveVel because it's needed for the pressure */
                }
                /* temperature scaling and pressure scaling to produce the extended variables at t+dt */
                if (!bInitStep)
                {
                    if (bTrotter)
                    {
                        m_add(force_vir, shake_vir, total_vir); /* we need the un-dispersion corrected total vir here */
                        trotter_update(ir, step, ekind, enerd, state, total_vir, mdatoms, &MassQ, trotter_seq, ettTSEQ2);
                    }
                    else
                    {
                        if (bExchanged)
                        {

                            /* We need the kinetic energy at minus the half step for determining
                             * the full step kinetic energy and possibly for T-coupling.*/
                            /* This may not be quite working correctly yet . . . . */
                            compute_globals(fplog, gstat, cr, ir, fr, ekind, state, state_global, mdatoms, nrnb, vcm,
                                            wcycle, enerd, NULL, NULL, NULL, NULL, mu_tot,
                                            constr, NULL, FALSE, state->box,
                                            top_global, &pcurr, top_global->natoms, &bSumEkinhOld,
                                            CGLO_RERUNMD | CGLO_GSTAT | CGLO_TEMPERATURE);
                        }
                    }
                }

                if (iterate.bIterationActive &&
                    done_iterating(cr, fplog, step, &iterate, bFirstIterate,
                                   state->veta, &vetanew))
                {
                    break;
                }
                bFirstIterate = FALSE;
            }

            if (bTrotter && !bInitStep)
            {
                copy_mat(shake_vir, state->svir_prev);
                copy_mat(force_vir, state->fvir_prev);
                if (IR_NVT_TROTTER(ir) && ir->eI == eiVV)
                {
                    /* update temperature and kinetic energy now that step is over - this is the v(t+dt) point */
                    enerd->term[F_TEMP] = sum_ekin(&(ir->opts), ekind, NULL, (ir->eI == eiVV), FALSE, FALSE);
                    enerd->term[F_EKIN] = trace(ekind->ekin);
                }
            }
            /* if it's the initial step, we performed this first step just to get the constraint virial */
            if (bInitStep && ir->eI == eiVV)
            {
                copy_rvecn(cbuf, state->v, 0, state->natoms);
            }

            GMX_MPE_LOG(ev_timestep1);
        }

        /* MRS -- now done iterating -- compute the conserved quantity */
        if (bVV)
        {
            saved_conserved_quantity = compute_conserved_from_auxiliary(ir, state, &MassQ);
            if (ir->eI == eiVV)
            {
                last_ekin = enerd->term[F_EKIN];
            }
            if ((ir->eDispCorr != edispcEnerPres) && (ir->eDispCorr != edispcAllEnerPres))
            {
                saved_conserved_quantity -= enerd->term[F_DISPCORR];
            }
            /* sum up the foreign energy and dhdl terms for vv.  currently done every step so that dhdl is correct in the .edr */
            if (!bRerunMD)
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

            lamnew = ExpandedEnsembleDynamics(fplog, ir, enerd, state, &MassQ, state->fep_state, &state->dfhist, step, mcrng, state->v, mdatoms);
            /* history is maintained in state->dfhist, but state_global is what is sent to trajectory and log output */
            copy_df_history(&state_global->dfhist,&state->dfhist);
        }
        /* ################## START TRAJECTORY OUTPUT ################# */

        /* Now we have the energies and forces corresponding to the
         * coordinates at time t. We must output all of this before
         * the update.
         * for RerunMD t is read from input trajectory
         */
        GMX_MPE_LOG(ev_output_start);

        mdof_flags = 0;
        if (do_per_step(step, ir->nstxout))
        {
            mdof_flags |= MDOF_X;
        }
        if (do_per_step(step, ir->nstvout))
        {
            mdof_flags |= MDOF_V;
        }
        if (do_per_step(step, ir->nstfout))
        {
            mdof_flags |= MDOF_F;
        }
        if (do_per_step(step, ir->nstxtcout))
        {
            mdof_flags |= MDOF_XTC;
        }
        if (bCPT)
        {
            mdof_flags |= MDOF_CPT;
        }
        ;

#if defined(GMX_FAHCORE) || defined(GMX_WRITELASTSTEP)
        if (bLastStep)
        {
            /* Enforce writing positions and velocities at end of run */
            mdof_flags |= (MDOF_X | MDOF_V);
        }
#endif
#ifdef GMX_FAHCORE
        if (MASTER(cr))
        {
            fcReportProgress( ir->nsteps, step );
        }

#if defined(__native_client__)
        fcCheckin(MASTER(cr));
#endif

        /* sync bCPT and fc record-keeping */
        if (bCPT && MASTER(cr))
        {
            fcRequestCheckPoint();
        }
#endif

        if (mdof_flags != 0)
        {
            wallcycle_start(wcycle, ewcTRAJ);
            if (bCPT)
            {
                if (state->flags & (1<<estLD_RNG))
                {
                    get_stochd_state(upd, state);
                }
                if (state->flags  & (1<<estMC_RNG))
                {
                    get_mc_state(mcrng, state);
                }
                if (MASTER(cr))
                {
                    if (bSumEkinhOld)
                    {
                        state_global->ekinstate.bUpToDate = FALSE;
                    }
                    else
                    {
                        update_ekinstate(&state_global->ekinstate, ekind);
                        state_global->ekinstate.bUpToDate = TRUE;
                    }
                    update_energyhistory(&state_global->enerhist, mdebin);
                }
            }
            write_traj(fplog, cr, outf, mdof_flags, top_global,
                       step, t, state, state_global, f, f_global, &n_xtc, &x_xtc);
            if (bCPT)
            {
                nchkpt++;
                bCPT = FALSE;
            }
            debug_gmx();
            if (bLastStep && step_rel == ir->nsteps &&
                (Flags & MD_CONFOUT) && MASTER(cr) &&
                !bRerunMD && !bFFscan)
            {
                /* x and v have been collected in write_traj,
                 * because a checkpoint file will always be written
                 * at the last step.
                 */
                fprintf(stderr, "\nWriting final coordinates.\n");
                if (fr->bMolPBC)
                {
                    /* Make molecules whole only for confout writing */
                    do_pbc_mtop(fplog, ir->ePBC, state->box, top_global, state_global->x);
                }
                write_sto_conf_mtop(ftp2fn(efSTO, nfile, fnm),
                                    *top_global->name, top_global,
                                    state_global->x, state_global->v,
                                    ir->ePBC, state->box);
                debug_gmx();
            }
            wallcycle_stop(wcycle, ewcTRAJ);
        }
        GMX_MPE_LOG(ev_output_finish);

        /* kludge -- virial is lost with restart for NPT control. Must restart */
        if (bStartingFromCpt && bVV)
        {
            copy_mat(state->svir_prev, shake_vir);
            copy_mat(state->fvir_prev, force_vir);
        }
        /*  ################## END TRAJECTORY OUTPUT ################ */

        /* Determine the wallclock run time up till now */
        run_time = gmx_gettime() - (double)runtime->real;

        /* Check whether everything is still allright */
        if (((int)gmx_get_stop_condition() > handled_stop_condition)
#ifdef GMX_THREAD_MPI
            && MASTER(cr)
#endif
            )
        {
            /* this is just make gs.sig compatible with the hack
               of sending signals around by MPI_Reduce with together with
               other floats */
            if (gmx_get_stop_condition() == gmx_stop_cond_next_ns)
            {
                gs.sig[eglsSTOPCOND] = 1;
            }
            if (gmx_get_stop_condition() == gmx_stop_cond_next)
            {
                gs.sig[eglsSTOPCOND] = -1;
            }
            /* < 0 means stop at next step, > 0 means stop at next NS step */
            if (fplog)
            {
                fprintf(fplog,
                        "\n\nReceived the %s signal, stopping at the next %sstep\n\n",
                        gmx_get_signal_name(),
                        gs.sig[eglsSTOPCOND] == 1 ? "NS " : "");
                fflush(fplog);
            }
            fprintf(stderr,
                    "\n\nReceived the %s signal, stopping at the next %sstep\n\n",
                    gmx_get_signal_name(),
                    gs.sig[eglsSTOPCOND] == 1 ? "NS " : "");
            fflush(stderr);
            handled_stop_condition = (int)gmx_get_stop_condition();
        }
        else if (MASTER(cr) && (bNS || ir->nstlist <= 0) &&
                 (max_hours > 0 && run_time > max_hours*60.0*60.0*0.99) &&
                 gs.sig[eglsSTOPCOND] == 0 && gs.set[eglsSTOPCOND] == 0)
        {
            /* Signal to terminate the run */
            gs.sig[eglsSTOPCOND] = 1;
            if (fplog)
            {
                fprintf(fplog, "\nStep %s: Run time exceeded %.3f hours, will terminate the run\n", gmx_step_str(step, sbuf), max_hours*0.99);
            }
            fprintf(stderr, "\nStep %s: Run time exceeded %.3f hours, will terminate the run\n", gmx_step_str(step, sbuf), max_hours*0.99);
        }

        if (bResetCountersHalfMaxH && MASTER(cr) &&
            run_time > max_hours*60.0*60.0*0.495)
        {
            gs.sig[eglsRESETCOUNTERS] = 1;
        }

        if (ir->nstlist == -1 && !bRerunMD)
        {
            /* When bGStatEveryStep=FALSE, global_stat is only called
             * when we check the atom displacements, not at NS steps.
             * This means that also the bonded interaction count check is not
             * performed immediately after NS. Therefore a few MD steps could
             * be performed with missing interactions.
             * But wrong energies are never written to file,
             * since energies are only written after global_stat
             * has been called.
             */
            if (step >= nlh.step_nscheck)
            {
                nlh.nabnsb = natoms_beyond_ns_buffer(ir, fr, &top->cgs,
                                                     nlh.scale_tot, state->x);
            }
            else
            {
                /* This is not necessarily true,
                 * but step_nscheck is determined quite conservatively.
                 */
                nlh.nabnsb = 0;
            }
        }

        /* In parallel we only have to check for checkpointing in steps
         * where we do global communication,
         *  otherwise the other nodes don't know.
         */
        if (MASTER(cr) && ((bGStat || !PAR(cr)) &&
                           cpt_period >= 0 &&
                           (cpt_period == 0 ||
                            run_time >= nchkpt*cpt_period*60.0)) &&
            gs.set[eglsCHKPT] == 0)
        {
            gs.sig[eglsCHKPT] = 1;
        }

        /* at the start of step, randomize or scale the velocities (trotter done elsewhere) */
        if (EI_VV(ir->eI))
        {
            if (!bInitStep)
            {
                update_tcouple(fplog, step, ir, state, ekind, wcycle, upd, &MassQ, mdatoms);
            }
            if (ETC_ANDERSEN(ir->etc)) /* keep this outside of update_tcouple because of the extra info required to pass */
            {
                gmx_bool bIfRandomize;
                bIfRandomize = update_randomize_velocities(ir, step, mdatoms, state, upd, &top->idef, constr, DOMAINDECOMP(cr));
                /* if we have constraints, we have to remove the kinetic energy parallel to the bonds */
                if (constr && bIfRandomize)
                {
                    update_constraints(fplog, step, NULL, ir, ekind, mdatoms,
                                       state, fr->bMolPBC, graph, f,
                                       &top->idef, tmp_vir, NULL,
                                       cr, nrnb, wcycle, upd, constr,
                                       bInitStep, TRUE, bCalcVir, vetanew);
                }
            }
        }

        if (bIterativeCase && do_per_step(step, ir->nstpcouple))
        {
            gmx_iterate_init(&iterate, TRUE);
            /* for iterations, we save these vectors, as we will be redoing the calculations */
            copy_coupling_state(state, bufstate, ekind, ekind_save, &(ir->opts));
        }

        bFirstIterate = TRUE;
        while (bFirstIterate || iterate.bIterationActive)
        {
            /* We now restore these vectors to redo the calculation with improved extended variables */
            if (iterate.bIterationActive)
            {
                copy_coupling_state(bufstate, state, ekind_save, ekind, &(ir->opts));
            }

            /* We make the decision to break or not -after- the calculation of Ekin and Pressure,
               so scroll down for that logic */

            /* #########   START SECOND UPDATE STEP ################# */
            GMX_MPE_LOG(ev_update_start);
            /* Box is changed in update() when we do pressure coupling,
             * but we should still use the old box for energy corrections and when
             * writing it to the energy file, so it matches the trajectory files for
             * the same timestep above. Make a copy in a separate array.
             */
            copy_mat(state->box, lastbox);

            bOK = TRUE;
            dvdl_constr = 0;

            if (!(bRerunMD && !rerun_fr.bV && !bForceUpdate))
            {
                wallcycle_start(wcycle, ewcUPDATE);
                /* UPDATE PRESSURE VARIABLES IN TROTTER FORMULATION WITH CONSTRAINTS */
                if (bTrotter)
                {
                    if (iterate.bIterationActive)
                    {
                        if (bFirstIterate)
                        {
                            scalevir = 1;
                        }
                        else
                        {
                            /* we use a new value of scalevir to converge the iterations faster */
                            scalevir = tracevir/trace(shake_vir);
                        }
                        msmul(shake_vir, scalevir, shake_vir);
                        m_add(force_vir, shake_vir, total_vir);
                        clear_mat(shake_vir);
                    }
                    trotter_update(ir, step, ekind, enerd, state, total_vir, mdatoms, &MassQ, trotter_seq, ettTSEQ3);
                    /* We can only do Berendsen coupling after we have summed
                     * the kinetic energy or virial. Since the happens
                     * in global_state after update, we should only do it at
                     * step % nstlist = 1 with bGStatEveryStep=FALSE.
                     */
                }
                else
                {
                    update_tcouple(fplog, step, ir, state, ekind, wcycle, upd, &MassQ, mdatoms);
                    update_pcouple(fplog, step, ir, state, pcoupl_mu, M, wcycle,
                                   upd, bInitStep);
                }

                if (bVV)
                {
                    bUpdateDoLR = (fr->bTwinRange && do_per_step(step, ir->nstcalclr));

                    /* velocity half-step update */
                    update_coords(fplog, step, ir, mdatoms, state, fr->bMolPBC, f,
                                  bUpdateDoLR, fr->f_twin, bCalcVir ? &fr->vir_twin_constr : NULL, fcd,
                                  ekind, M, wcycle, upd, FALSE, etrtVELOCITY2,
                                  cr, nrnb, constr, &top->idef);
                }

                /* Above, initialize just copies ekinh into ekin,
                 * it doesn't copy position (for VV),
                 * and entire integrator for MD.
                 */

                if (ir->eI == eiVVAK)
                {
                    copy_rvecn(state->x, cbuf, 0, state->natoms);
                }
                bUpdateDoLR = (fr->bTwinRange && do_per_step(step, ir->nstcalclr));

                update_coords(fplog, step, ir, mdatoms, state, fr->bMolPBC, f,
                              bUpdateDoLR, fr->f_twin, bCalcVir ? &fr->vir_twin_constr : NULL, fcd,
                              ekind, M, wcycle, upd, bInitStep, etrtPOSITION, cr, nrnb, constr, &top->idef);
                wallcycle_stop(wcycle, ewcUPDATE);

                update_constraints(fplog, step, &dvdl_constr, ir, ekind, mdatoms, state,
                                   fr->bMolPBC, graph, f,
                                   &top->idef, shake_vir, force_vir,
                                   cr, nrnb, wcycle, upd, constr,
                                   bInitStep, FALSE, bCalcVir, state->veta);

                if (bCalcVir && bUpdateDoLR && ir->nstcalclr > 1)
                {
                    /* Correct the virial for multiple time stepping */
                    m_sub(shake_vir, fr->vir_twin_constr, shake_vir);
                }

                if (ir->eI == eiVVAK)
                {
                    /* erase F_EKIN and F_TEMP here? */
                    /* just compute the kinetic energy at the half step to perform a trotter step */
                    compute_globals(fplog, gstat, cr, ir, fr, ekind, state, state_global, mdatoms, nrnb, vcm,
                                    wcycle, enerd, force_vir, shake_vir, total_vir, pres, mu_tot,
                                    constr, NULL, FALSE, lastbox,
                                    top_global, &pcurr, top_global->natoms, &bSumEkinhOld,
                                    cglo_flags | CGLO_TEMPERATURE
                                    );
                    wallcycle_start(wcycle, ewcUPDATE);
                    trotter_update(ir, step, ekind, enerd, state, total_vir, mdatoms, &MassQ, trotter_seq, ettTSEQ4);
                    /* now we know the scaling, we can compute the positions again again */
                    copy_rvecn(cbuf, state->x, 0, state->natoms);

                    bUpdateDoLR = (fr->bTwinRange && do_per_step(step, ir->nstcalclr));

                    update_coords(fplog, step, ir, mdatoms, state, fr->bMolPBC, f,
                                  bUpdateDoLR, fr->f_twin, bCalcVir ? &fr->vir_twin_constr : NULL, fcd,
                                  ekind, M, wcycle, upd, bInitStep, etrtPOSITION, cr, nrnb, constr, &top->idef);
                    wallcycle_stop(wcycle, ewcUPDATE);

                    /* do we need an extra constraint here? just need to copy out of state->v to upd->xp? */
                    /* are the small terms in the shake_vir here due
                     * to numerical errors, or are they important
                     * physically? I'm thinking they are just errors, but not completely sure.
                     * For now, will call without actually constraining, constr=NULL*/
                    update_constraints(fplog, step, NULL, ir, ekind, mdatoms,
                                       state, fr->bMolPBC, graph, f,
                                       &top->idef, tmp_vir, force_vir,
                                       cr, nrnb, wcycle, upd, NULL,
                                       bInitStep, FALSE, bCalcVir,
                                       state->veta);
                }
                if (!bOK && !bFFscan)
                {
                    gmx_fatal(FARGS, "Constraint error: Shake, Lincs or Settle could not solve the constrains");
                }

                if (fr->bSepDVDL && fplog && do_log)
                {
                    fprintf(fplog, sepdvdlformat, "Constraint dV/dl", 0.0, dvdl_constr);
                }
                if (bVV)
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
                unshift_self(graph, state->box, state->x);
            }

            GMX_BARRIER(cr->mpi_comm_mygroup);
            GMX_MPE_LOG(ev_update_finish);

            if (vsite != NULL)
            {
                wallcycle_start(wcycle, ewcVSITECONSTR);
                if (graph != NULL)
                {
                    shift_self(graph, state->box, state->x);
                }
                construct_vsites(fplog, vsite, state->x, nrnb, ir->delta_t, state->v,
                                 top->idef.iparams, top->idef.il,
                                 fr->ePBC, fr->bMolPBC, graph, cr, state->box);

                if (graph != NULL)
                {
                    unshift_self(graph, state->box, state->x);
                }
                wallcycle_stop(wcycle, ewcVSITECONSTR);
            }

            /* ############## IF NOT VV, Calculate globals HERE, also iterate constraints  ############ */
            /* With Leap-Frog we can skip compute_globals at
             * non-communication steps, but we need to calculate
             * the kinetic energy one step before communication.
             */
            if (bGStat || (!EI_VV(ir->eI) && do_per_step(step+1, nstglobalcomm)))
            {
                if (ir->nstlist == -1 && bFirstIterate)
                {
                    gs.sig[eglsNABNSB] = nlh.nabnsb;
                }
                compute_globals(fplog, gstat, cr, ir, fr, ekind, state, state_global, mdatoms, nrnb, vcm,
                                wcycle, enerd, force_vir, shake_vir, total_vir, pres, mu_tot,
                                constr,
                                bFirstIterate ? &gs : NULL,
                                (step_rel % gs.nstms == 0) &&
                                (multisim_nsteps < 0 || (step_rel < multisim_nsteps)),
                                lastbox,
                                top_global, &pcurr, top_global->natoms, &bSumEkinhOld,
                                cglo_flags
                                | (!EI_VV(ir->eI) || bRerunMD ? CGLO_ENERGY : 0)
                                | (!EI_VV(ir->eI) && bStopCM ? CGLO_STOPCM : 0)
                                | (!EI_VV(ir->eI) ? CGLO_TEMPERATURE : 0)
                                | (!EI_VV(ir->eI) || bRerunMD ? CGLO_PRESSURE : 0)
                                | (iterate.bIterationActive ? CGLO_ITERATE : 0)
                                | (bFirstIterate ? CGLO_FIRSTITERATE : 0)
                                | CGLO_CONSTRAINT
                                );
                if (ir->nstlist == -1 && bFirstIterate)
                {
                    nlh.nabnsb         = gs.set[eglsNABNSB];
                    gs.set[eglsNABNSB] = 0;
                }
            }
            /* bIterate is set to keep it from eliminating the old ekin kinetic energy terms */
            /* #############  END CALC EKIN AND PRESSURE ################# */

            /* Note: this is OK, but there are some numerical precision issues with using the convergence of
               the virial that should probably be addressed eventually. state->veta has better properies,
               but what we actually need entering the new cycle is the new shake_vir value. Ideally, we could
               generate the new shake_vir, but test the veta value for convergence.  This will take some thought. */

            if (iterate.bIterationActive &&
                done_iterating(cr, fplog, step, &iterate, bFirstIterate,
                               trace(shake_vir), &tracevir))
            {
                break;
            }
            bFirstIterate = FALSE;
        }

        if (!bVV || bRerunMD)
        {
            /* sum up the foreign energy and dhdl terms for md and sd. currently done every step so that dhdl is correct in the .edr */
            sum_dhdl(enerd, state->lambda, ir->fepvals);
        }
        update_box(fplog, step, ir, mdatoms, state, graph, f,
                   ir->nstlist == -1 ? &nlh.scale_tot : NULL, pcoupl_mu, nrnb, wcycle, upd, bInitStep, FALSE);

        /* ################# END UPDATE STEP 2 ################# */
        /* #### We now have r(t+dt) and v(t+dt/2)  ############# */

        /* The coordinates (x) were unshifted in update */
        if (bFFscan && (shellfc == NULL || bConverged))
        {
            if (print_forcefield(fplog, enerd->term, mdatoms->homenr,
                                 f, NULL, xcopy,
                                 &(top_global->mols), mdatoms->massT, pres))
            {
                gmx_finalize_par();

                fprintf(stderr, "\n");
                exit(0);
            }
        }
        if (!bGStat)
        {
            /* We will not sum ekinh_old,
             * so signal that we still have to do it.
             */
            bSumEkinhOld = TRUE;
        }

        if (bTCR)
        {
            /* Only do GCT when the relaxation of shells (minimization) has converged,
             * otherwise we might be coupling to bogus energies.
             * In parallel we must always do this, because the other sims might
             * update the FF.
             */

            /* Since this is called with the new coordinates state->x, I assume
             * we want the new box state->box too. / EL 20040121
             */
            do_coupling(fplog, oenv, nfile, fnm, tcr, t, step, enerd->term, fr,
                        ir, MASTER(cr),
                        mdatoms, &(top->idef), mu_aver,
                        top_global->mols.nr, cr,
                        state->box, total_vir, pres,
                        mu_tot, state->x, f, bConverged);
            debug_gmx();
        }

        /* #########  BEGIN PREPARING EDR OUTPUT  ###########  */

        /* use the directly determined last velocity, not actually the averaged half steps */
        if (bTrotter && ir->eI == eiVV)
        {
            enerd->term[F_EKIN] = last_ekin;
        }
        enerd->term[F_ETOT] = enerd->term[F_EPOT] + enerd->term[F_EKIN];

        if (bVV)
        {
            enerd->term[F_ECONSERVED] = enerd->term[F_ETOT] + saved_conserved_quantity;
        }
        else
        {
            enerd->term[F_ECONSERVED] = enerd->term[F_ETOT] + compute_conserved_from_auxiliary(ir, state, &MassQ);
        }
        /* Check for excessively large energies */
        if (bIonize)
        {
#ifdef GMX_DOUBLE
            real etot_max = 1e200;
#else
            real etot_max = 1e30;
#endif
            if (fabs(enerd->term[F_ETOT]) > etot_max)
            {
                fprintf(stderr, "Energy too large (%g), giving up\n",
                        enerd->term[F_ETOT]);
            }
        }
        /* #########  END PREPARING EDR OUTPUT  ###########  */

        /* Time for performance */
        if (((step % stepout) == 0) || bLastStep)
        {
            runtime_upd_proc(runtime);
        }

        /* Output stuff */
        if (MASTER(cr))
        {
            gmx_bool do_dr, do_or;

            if (fplog && do_log && bDoExpanded)
            {
                /* only needed if doing expanded ensemble */
                PrintFreeEnergyInfoToFile(fplog, ir->fepvals, ir->expandedvals, ir->bSimTemp ? ir->simtempvals : NULL,
                                          &state_global->dfhist, state->fep_state, ir->nstlog, step);
            }
            if (!(bStartingFromCpt && (EI_VV(ir->eI))))
            {
                if (bCalcEner)
                {
                    upd_mdebin(mdebin, bDoDHDL, TRUE,
                               t, mdatoms->tmass, enerd, state,
                               ir->fepvals, ir->expandedvals, lastbox,
                               shake_vir, force_vir, total_vir, pres,
                               ekind, mu_tot, constr);
                }
                else
                {
                    upd_mdebin_step(mdebin);
                }

                do_dr  = do_per_step(step, ir->nstdisreout);
                do_or  = do_per_step(step, ir->nstorireout);

                print_ebin(outf->fp_ene, do_ene, do_dr, do_or, do_log ? fplog : NULL,
                           step, t,
                           eprNORMAL, bCompact, mdebin, fcd, groups, &(ir->opts));
            }
            if (ir->ePull != epullNO)
            {
                pull_print_output(ir->pull, step, t);
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

        /* Remaining runtime */
        if (MULTIMASTER(cr) && (do_verbose || gmx_got_usr_signal()) && !bPMETuneRunning)
        {
            if (shellfc)
            {
                fprintf(stderr, "\n");
            }
            print_time(stderr, runtime, step, ir, cr);
        }

        /* Replica exchange */
        bExchanged = FALSE;
        if ((repl_ex_nst > 0) && (step > 0) && !bLastStep &&
            do_per_step(step, repl_ex_nst))
        {
            bExchanged = replica_exchange(fplog, cr, repl_ex,
                                          state_global, enerd,
                                          state, step, t);

            if (bExchanged && DOMAINDECOMP(cr))
            {
                dd_partition_system(fplog, step, cr, TRUE, 1,
                                    state_global, top_global, ir,
                                    state, &f, mdatoms, top, fr,
                                    vsite, shellfc, constr,
                                    nrnb, wcycle, FALSE);
            }
        }

        bFirstStep       = FALSE;
        bInitStep        = FALSE;
        bStartingFromCpt = FALSE;

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

        if ( (membed != NULL) && (!bLastStep) )
        {
            rescale_membed(step_rel, membed, state_global->x);
        }

        if (bRerunMD)
        {
            if (MASTER(cr))
            {
                /* read next frame from input trajectory */
                bNotLastFrame = read_next_frame(oenv, status, &rerun_fr);
            }

            if (PAR(cr))
            {
                rerun_parallel_comm(cr, &rerun_fr, &bNotLastFrame);
            }
        }

        if (!bRerunMD || !rerun_fr.bStep)
        {
            /* increase the MD step number */
            step++;
            step_rel++;
        }

        cycles = wallcycle_stop(wcycle, ewcSTEP);
        if (DOMAINDECOMP(cr) && wcycle)
        {
            dd_cycles_add(cr->dd, cycles, ddCyclStep);
        }

        if (bPMETuneRunning || bPMETuneTry)
        {
            /* PME grid + cut-off optimization with GPUs or PME nodes */

            /* Count the total cycles over the last steps */
            cycles_pmes += cycles;

            /* We can only switch cut-off at NS steps */
            if (step % ir->nstlist == 0)
            {
                /* PME grid + cut-off optimization with GPUs or PME nodes */
                if (bPMETuneTry)
                {
                    if (DDMASTER(cr->dd))
                    {
                        /* PME node load is too high, start tuning */
                        bPMETuneRunning = (dd_pme_f_ratio(cr->dd) >= 1.05);
                    }
                    dd_bcast(cr->dd, sizeof(gmx_bool), &bPMETuneRunning);

                    if (bPMETuneRunning || step_rel > ir->nstlist*50)
                    {
                        bPMETuneTry     = FALSE;
                    }
                }
                if (bPMETuneRunning)
                {
                    /* init_step might not be a multiple of nstlist,
                     * but the first cycle is always skipped anyhow.
                     */
                    bPMETuneRunning =
                        pme_load_balance(pme_loadbal, cr,
                                         (bVerbose && MASTER(cr)) ? stderr : NULL,
                                         fplog,
                                         ir, state, cycles_pmes,
                                         fr->ic, fr->nbv, &fr->pmedata,
                                         step);

                    /* Update constants in forcerec/inputrec to keep them in sync with fr->ic */
                    fr->ewaldcoeff = fr->ic->ewaldcoeff;
                    fr->rlist      = fr->ic->rlist;
                    fr->rlistlong  = fr->ic->rlistlong;
                    fr->rcoulomb   = fr->ic->rcoulomb;
                    fr->rvdw       = fr->ic->rvdw;
                }
                cycles_pmes = 0;
            }
        }

        if (step_rel == wcycle_get_reset_counters(wcycle) ||
            gs.set[eglsRESETCOUNTERS] != 0)
        {
            /* Reset all the counters related to performance over the run */
            reset_all_counters(fplog, cr, step, &step_rel, ir, wcycle, nrnb, runtime,
                               fr->nbv != NULL && fr->nbv->bUseGPU ? fr->nbv->cu_nbv : NULL);
            wcycle_set_reset_counters(wcycle, -1);
            if (!(cr->duty & DUTY_PME))
            {
                /* Tell our PME node to reset its counters */
                gmx_pme_send_resetcounters(cr, step);
            }
            /* Correct max_hours for the elapsed time */
            max_hours                -= run_time/(60.0*60.0);
            bResetCountersHalfMaxH    = FALSE;
            gs.set[eglsRESETCOUNTERS] = 0;
        }

    }
    /* End of main MD loop */
    debug_gmx();

    /* Stop the time */
    runtime_end(runtime);

    if (bRerunMD && MASTER(cr))
    {
        close_trj(status);
    }

    if (!(cr->duty & DUTY_PME))
    {
        /* Tell the PME only node to finish */
        gmx_pme_send_finish(cr);
    }

    if (MASTER(cr))
    {
        if (ir->nstcalcenergy > 0 && !bRerunMD)
        {
            print_ebin(outf->fp_ene, FALSE, FALSE, FALSE, fplog, step, t,
                       eprAVER, FALSE, mdebin, fcd, groups, &(ir->opts));
        }
    }

    done_mdoutf(outf);

    debug_gmx();

    if (ir->nstlist == -1 && nlh.nns > 0 && fplog)
    {
        fprintf(fplog, "Average neighborlist lifetime: %.1f steps, std.dev.: %.1f steps\n", nlh.s1/nlh.nns, sqrt(nlh.s2/nlh.nns - sqr(nlh.s1/nlh.nns)));
        fprintf(fplog, "Average number of atoms that crossed the half buffer length: %.1f\n\n", nlh.ab/nlh.nns);
    }

    if (pme_loadbal != NULL)
    {
        pme_loadbal_done(pme_loadbal, cr, fplog,
                         fr->nbv != NULL && fr->nbv->bUseGPU);
    }

    if (shellfc && fplog)
    {
        fprintf(fplog, "Fraction of iterations that converged:           %.2f %%\n",
                (nconverged*100.0)/step_rel);
        fprintf(fplog, "Average number of force evaluations per MD step: %.2f\n\n",
                tcount/step_rel);
    }

    if (repl_ex_nst > 0 && MASTER(cr))
    {
        print_replica_exchange_statistics(fplog, repl_ex);
    }

    runtime->nsteps_done = step_rel;

    return 0;
}
