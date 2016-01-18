/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011,2012,2013,2014,2015,2016, by the GROMACS development team, led by
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "thread_mpi/threads.h"

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_network.h"
#include "gromacs/ewald/pme-load-balancing.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/fileio/filenm.h"
#include "gromacs/fileio/mdoutf.h"
#include "gromacs/fileio/trajectory_writing.h"
#include "gromacs/fileio/trx.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/imd/imd.h"
#include "gromacs/legacyheaders/constr.h"
#include "gromacs/legacyheaders/ebin.h"
#include "gromacs/legacyheaders/force.h"
#include "gromacs/legacyheaders/md_logging.h"
#include "gromacs/legacyheaders/md_support.h"
#include "gromacs/legacyheaders/mdatoms.h"
#include "gromacs/legacyheaders/mdebin.h"
#include "gromacs/legacyheaders/mdrun.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/nrnb.h"
#include "gromacs/legacyheaders/ns.h"
#include "gromacs/legacyheaders/shellfc.h"
#include "gromacs/legacyheaders/sighandler.h"
#include "gromacs/legacyheaders/sim_util.h"
#include "gromacs/legacyheaders/tgroup.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/update.h"
#include "gromacs/legacyheaders/vcm.h"
#include "gromacs/legacyheaders/vsite.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/legacyheaders/types/constr.h"
#include "gromacs/legacyheaders/types/enums.h"
#include "gromacs/legacyheaders/types/fcdata.h"
#include "gromacs/legacyheaders/types/force_flags.h"
#include "gromacs/legacyheaders/types/forcerec.h"
#include "gromacs/legacyheaders/types/group.h"
#include "gromacs/legacyheaders/types/inputrec.h"
#include "gromacs/legacyheaders/types/interaction_const.h"
#include "gromacs/legacyheaders/types/mdatom.h"
#include "gromacs/legacyheaders/types/membedt.h"
#include "gromacs/legacyheaders/types/nrnb.h"
#include "gromacs/legacyheaders/types/oenv.h"
#include "gromacs/legacyheaders/types/shellfc.h"
#include "gromacs/legacyheaders/types/state.h"
#include "gromacs/listed-forces/manage-threading.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/compute_io.h"
#include "gromacs/mdlib/mdrun_signalling.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_gpu_data_mgmt.h"
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
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#include "deform.h"
#include "membed.h"
#include "repl_ex.h"

#ifdef GMX_FAHCORE
#include "corewrap.h"
#endif

static void reset_all_counters(FILE *fplog, t_commrec *cr,
                               gmx_int64_t step,
                               gmx_int64_t *step_rel, t_inputrec *ir,
                               gmx_wallcycle_t wcycle, t_nrnb *nrnb,
                               gmx_walltime_accounting_t walltime_accounting,
                               struct nonbonded_verlet_t *nbv)
{
    char sbuf[STEPSTRSIZE];

    /* Reset all the counters related to performance over the run */
    md_print_warn(cr, fplog, "step %s: resetting all time and cycle counters\n",
                  gmx_step_str(step, sbuf));

    if (use_GPU(nbv))
    {
        nbnxn_gpu_reset_timings(nbv);
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
             int imdport,
             unsigned long Flags,
             gmx_walltime_accounting_t walltime_accounting)
{
    gmx_mdoutf_t    outf = NULL;
    gmx_int64_t     step, step_rel;
    double          elapsed_time;
    double          t, t0, lam0[efptNR];
    gmx_bool        bGStatEveryStep, bGStat, bCalcVir, bCalcEnerStep, bCalcEner;
    gmx_bool        bNS, bNStList, bSimAnn, bStopCM, bRerunMD, bNotLastFrame = FALSE,
                    bFirstStep, bStateFromCP, bStateFromTPX, bInitStep, bLastStep,
                    bBornRadii, bStartingFromCpt;
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
    matrix            pcoupl_mu, M;
    t_trxframe        rerun_fr;
    gmx_repl_ex_t     repl_ex = NULL;
    int               nchkpt  = 1;
    gmx_localtop_t   *top;
    t_mdebin         *mdebin   = NULL;
    t_state          *state    = NULL;
    rvec             *f_global = NULL;
    gmx_enerdata_t   *enerd;
    rvec             *f = NULL;
    gmx_global_stat_t gstat;
    gmx_update_t      upd   = NULL;
    t_graph          *graph = NULL;
    gmx_signalling_t  gs;
    gmx_groups_t     *groups;
    gmx_ekindata_t   *ekind;
    gmx_shellfc_t     shellfc;
    int               count, nconverged = 0;
    double            tcount                 = 0;
    gmx_bool          bConverged             = TRUE, bSumEkinhOld, bDoReplEx, bExchanged, bNeedRepartition;
    gmx_bool          bResetCountersHalfMaxH = FALSE;
    gmx_bool          bVV, bTemp, bPres, bTrotter;
    gmx_bool          bUpdateDoLR;
    real              dvdl_constr;
    rvec             *cbuf        = NULL;
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
    gmx_int64_t       multisim_nsteps        = -1;                 /* number of steps to do  before first multisim
                                                                          simulation stops. If equal to zero, don't
                                                                          communicate any more between multisims.*/
    /* PME load balancing data for GPU kernels */
    pme_load_balancing_t *pme_loadbal      = NULL;
    gmx_bool              bPMETune         = FALSE;
    gmx_bool              bPMETunePrinting = FALSE;

    /* Interactive MD */
    gmx_bool          bIMDstep = FALSE;

#ifdef GMX_FAHCORE
    /* Temporary addition for FAHCORE checkpointing */
    int chkpt_ret;
#endif

    /* Check for special mdrun options */
    bRerunMD = (Flags & MD_RERUN);
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
    bVV      = EI_VV(ir->eI);
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

    if (bRerunMD)
    {
        ir->nstxout_compressed = 0;
    }
    groups = &top_global->groups;

    /* Initial values */
    init_md(fplog, cr, ir, oenv, &t, &t0, state_global->lambda,
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
        gmx_fatal(FARGS, "Shell particles are not implemented with domain decomposition, use a single rank");
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

        snew(state, 1);
        dd_init_local_state(cr->dd, state_global, state);

        if (DDMASTER(cr->dd) && ir->nstfout)
        {
            snew(f_global, state_global->natoms);
        }
    }
    else
    {
        top = gmx_mtop_generate_local_top(top_global, ir);

        forcerec_set_excl_load(fr, top);

        state    = serial_init_local_state(state_global);
        f_global = f;

        atoms2md(top_global, ir, 0, NULL, top_global->natoms, mdatoms);

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
    }

    /* Set up interactive MD (IMD) */
    init_IMD(ir, cr, top_global, fplog, ir->nstcalcenergy, state_global->x,
             nfile, fnm, oenv, imdport, Flags);

    if (DOMAINDECOMP(cr))
    {
        /* Distribute the charge groups over the nodes from the master node */
        dd_partition_system(fplog, ir->init_step, cr, TRUE, 1,
                            state_global, top_global, ir,
                            state, &f, mdatoms, top, fr,
                            vsite, shellfc, constr,
                            nrnb, NULL, FALSE);

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
        init_expanded_ensemble(bStateFromCP, ir, &state->dfhist);
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
     * with only LJ PME, or for reruns.
     */
    bPMETune = ((Flags & MD_TUNEPME) && EEL_PME(fr->eeltype) && !bRerunMD &&
                !(Flags & MD_REPRODUCIBLE));
    if (bPMETune)
    {
        pme_loadbal_init(&pme_loadbal, cr, fplog, ir, state->box,
                         fr->ic, fr->pmedata, use_GPU(fr->nbv),
                         &bPMETunePrinting);
    }

    if (!ir->bContinuation && !bRerunMD)
    {
        if (mdatoms->cFREEZE && (state->flags & (1<<estV)))
        {
            /* Set the velocities of frozen particles to zero */
            for (i = 0; i < mdatoms->homenr; i++)
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
            do_constrain_first(fplog, constr, ir, mdatoms, state,
                               cr, nrnb, fr, top);
        }
        if (vsite)
        {
            /* Construct the virtual sites for the initial configuration */
            construct_vsites(vsite, state->x, ir->delta_t, NULL,
                             top->idef.iparams, top->idef.il,
                             fr->ePBC, fr->bMolPBC, cr, state->box);
        }
    }

    debug_gmx();

    if (IR_TWINRANGE(*ir) && repl_ex_nst % ir->nstcalclr != 0)
    {
        /* We should exchange at nstcalclr steps to get correct integration */
        gmx_fatal(FARGS, "The replica exchange period (%d) is not divisible by nstcalclr (%d)", repl_ex_nst, ir->nstcalclr);
    }

    if (ir->efep != efepNO)
    {
        /* Set free energy calculation frequency as the greatest common
         * denominator of nstdhdl and repl_ex_nst.
         * Check for nstcalclr with twin-range, since we need the long-range
         * contribution to the free-energy at the correct (nstcalclr) steps.
         */
        nstfep = ir->fepvals->nstdhdl;
        if (ir->bExpanded)
        {
            if (IR_TWINRANGE(*ir) &&
                ir->expandedvals->nstexpanded % ir->nstcalclr != 0)
            {
                gmx_fatal(FARGS, "nstexpanded should be divisible by nstcalclr");
            }
            nstfep = gmx_greatest_common_divisor(ir->expandedvals->nstexpanded, nstfep);
        }
        if (repl_ex_nst > 0)
        {
            nstfep = gmx_greatest_common_divisor(repl_ex_nst, nstfep);
        }
        /* We checked divisibility of repl_ex_nst and nstcalclr above */
        if (IR_TWINRANGE(*ir) && nstfep % ir->nstcalclr != 0)
        {
            gmx_incons("nstfep not divisible by nstcalclr");
        }
    }

    /* Be REALLY careful about what flags you set here. You CANNOT assume
     * this is the first step, since we might be restarting from a checkpoint,
     * and in that case we should not do any modifications to the state.
     */
    bStopCM = (ir->comm_mode != ecmNO && !ir->bContinuation);

    cglo_flags = (CGLO_TEMPERATURE | CGLO_GSTAT
                  | (bStopCM ? CGLO_STOPCM : 0)
                  | (bVV ? CGLO_PRESSURE : 0)
                  | (bVV ? CGLO_CONSTRAINT : 0)
                  | (bRerunMD ? CGLO_RERUNMD : 0)
                  | ((Flags & MD_READ_EKIN) ? CGLO_READEKIN : 0));

    bSumEkinhOld = FALSE;
    compute_globals(fplog, gstat, cr, ir, fr, ekind, state, state_global, mdatoms, nrnb, vcm,
                    NULL, enerd, force_vir, shake_vir, total_vir, pres, mu_tot,
                    constr, NULL, FALSE, state->box,
                    top_global, &bSumEkinhOld, cglo_flags);
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
                        top_global, &bSumEkinhOld,
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
    bSumEkinhOld     = FALSE;
    bExchanged       = FALSE;
    bNeedRepartition = FALSE;

    init_global_signals(&gs, cr, ir, repl_ex_nst);

    step     = ir->init_step;
    step_rel = 0;

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

        /* Determine if this is a neighbor search step */
        bNStList = (ir->nstlist > 0  && step % ir->nstlist == 0);

        if (bPMETune && bNStList)
        {
            /* PME grid + cut-off optimization with GPUs or PME nodes */
            pme_loadbal_do(pme_loadbal, cr,
                           (bVerbose && MASTER(cr)) ? stderr : NULL,
                           fplog,
                           ir, fr, state, wcycle,
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

        if (ir->efep != efepNO || ir->bSimTemp)
        {
            /* find and set the current lambdas.  If rerunning, we either read in a state, or a lambda value,
               requiring different logic. */

            set_current_lambdas(step, ir->fepvals, bRerunMD, &rerun_fr, state_global, state, lam0);
            bDoDHDL      = do_per_step(step, ir->fepvals->nstdhdl);
            bDoFEP       = ((ir->efep != efepNO) && do_per_step(step, nstfep));
            bDoExpanded  = (do_per_step(step, ir->expandedvals->nstexpanded)
                            && (ir->bExpanded) && (step > 0) && (!bStartingFromCpt));
        }

        bDoReplEx = ((repl_ex_nst > 0) && (step > 0) && !bLastStep &&
                     do_per_step(step, repl_ex_nst));

        if (bSimAnn)
        {
            update_annealing_target_temp(&(ir->opts), t);
        }

        if (bRerunMD)
        {
            if (!DOMAINDECOMP(cr) || MASTER(cr))
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
                    gmx_fatal(FARGS, "Vsite recalculation with -rerun is not implemented with domain decomposition, use a single rank");
                }
                if (graph)
                {
                    /* Following is necessary because the graph may get out of sync
                     * with the coordinates if we only have every N'th coordinate set
                     */
                    mk_mshift(fplog, graph, fr->ePBC, state->box, state->x);
                    shift_self(graph, state->box, state->x);
                }
                construct_vsites(vsite, state->x, ir->delta_t, state->v,
                                 top->idef.iparams, top->idef.il,
                                 fr->ePBC, fr->bMolPBC, cr, state->box);
                if (graph)
                {
                    unshift_self(graph, state->box, state->x);
                }
            }
        }

        /* Stop Center of Mass motion */
        bStopCM = (ir->comm_mode != ecmNO && do_per_step(step, ir->nstcomm));

        if (bRerunMD)
        {
            /* for rerun MD always do Neighbour Searching */
            bNS      = (bFirstStep || ir->nstlist != 0);
            bNStList = bNS;
        }
        else
        {
            /* Determine whether or not to do Neighbour Searching and LR */
            bNS = (bFirstStep || bNStList || bExchanged || bNeedRepartition);
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

        /* do_log triggers energy and virial calculation. Because this leads
         * to different code paths, forces can be different. Thus for exact
         * continuation we should avoid extra log output.
         * Note that the || bLastStep can result in non-exact continuation
         * beyond the last step. But we don't consider that to be an issue.
         */
        do_log     = do_per_step(step, ir->nstlog) || (bFirstStep && !bStateFromCP) || bLastStep;
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
                dd_partition_system(fplog, step, cr,
                                    bMasterState, nstglobalcomm,
                                    state_global, top_global, ir,
                                    state, &f, mdatoms, top, fr,
                                    vsite, shellfc, constr,
                                    nrnb, wcycle,
                                    do_verbose && !bPMETunePrinting);
            }
        }

        if (MASTER(cr) && do_log)
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
                            top_global, &bSumEkinhOld,
                            CGLO_RERUNMD | CGLO_GSTAT | CGLO_TEMPERATURE);
        }
        clear_mat(force_vir);

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
                  (bVV && IR_NVT_TROTTER(ir) && do_per_step(step-1, nstglobalcomm)));

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
            count = relax_shell_flexcon(fplog, cr, bVerbose, step,
                                        ir, bNS, force_flags,
                                        top,
                                        constr, enerd, fcd,
                                        state, f, force_vir, mdatoms,
                                        nrnb, wcycle, graph, groups,
                                        shellfc, fr, bBornRadii, t, mu_tot,
                                        &bConverged, vsite,
                                        mdoutf_get_fp_field(outf));
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
            do_force(fplog, cr, ir, step, nrnb, wcycle, top, groups,
                     state->box, state->x, &state->hist,
                     f, force_vir, mdatoms, enerd, fcd,
                     state->lambda, graph,
                     fr, vsite, mu_tot, t, mdoutf_get_fp_field(outf), ed, bBornRadii,
                     (bNS ? GMX_FORCE_NS : 0) | force_flags);
        }

        if (bVV && !bStartingFromCpt && !bRerunMD)
        /*  ############### START FIRST UPDATE HALF-STEP FOR VV METHODS############### */
        {
            rvec *vbuf = NULL;

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
                copy_rvecn(state->v, vbuf, 0, state->natoms); /* should make this better for parallelizing? */
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
                          ekind, M, upd, bInitStep, etrtVELOCITY1,
                          cr, nrnb, constr, &top->idef);

            if (!bRerunMD || rerun_fr.bV || bForceUpdate)         /* Why is rerun_fr.bV here?  Unclear. */
            {
                wallcycle_stop(wcycle, ewcUPDATE);
                update_constraints(fplog, step, NULL, ir, mdatoms,
                                   state, fr->bMolPBC, graph, f,
                                   &top->idef, shake_vir,
                                   cr, nrnb, wcycle, upd, constr,
                                   TRUE, bCalcVir);
                wallcycle_start(wcycle, ewcUPDATE);
                if (bCalcVir && bUpdateDoLR && ir->nstcalclr > 1)
                {
                    /* Correct the virial for multiple time stepping */
                    m_sub(shake_vir, fr->vir_twin_constr, shake_vir);
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
            if (bCalcEner && ir->eI == eiVVAK)
            {
                bSumEkinhOld = TRUE;
            }
            /* for vv, the first half of the integration actually corresponds to the previous step.
               So we need information from the last step in the first half of the integration */
            if (bGStat || do_per_step(step-1, nstglobalcomm))
            {
                wallcycle_stop(wcycle, ewcUPDATE);
                compute_globals(fplog, gstat, cr, ir, fr, ekind, state, state_global, mdatoms, nrnb, vcm,
                                wcycle, enerd, force_vir, shake_vir, total_vir, pres, mu_tot,
                                constr, NULL, FALSE, state->box,
                                top_global, &bSumEkinhOld,
                                cglo_flags
                                | CGLO_ENERGY
                                | (bTemp ? CGLO_TEMPERATURE : 0)
                                | (bPres ? CGLO_PRESSURE : 0)
                                | (bPres ? CGLO_CONSTRAINT : 0)
                                | (bStopCM ? CGLO_STOPCM : 0)
                                | CGLO_SCALEEKIN
                                );
                /* explanation of above:
                   a) We compute Ekin at the full time step
                   if 1) we are using the AveVel Ekin, and it's not the
                   initial step, or 2) if we are using AveEkin, but need the full
                   time step kinetic energy for the pressure (always true now, since we want accurate statistics).
                   b) If we are using EkinAveEkin for the kinetic energy for the temperature control, we still feed in
                   EkinAveVel because it's needed for the pressure */
                wallcycle_start(wcycle, ewcUPDATE);
            }
            /* temperature scaling and pressure scaling to produce the extended variables at t+dt */
            if (!bInitStep)
            {
                if (bTrotter)
                {
                    m_add(force_vir, shake_vir, total_vir);     /* we need the un-dispersion corrected total vir here */
                    trotter_update(ir, step, ekind, enerd, state, total_vir, mdatoms, &MassQ, trotter_seq, ettTSEQ2);
                }
                else
                {
                    if (bExchanged)
                    {
                        wallcycle_stop(wcycle, ewcUPDATE);
                        /* We need the kinetic energy at minus the half step for determining
                         * the full step kinetic energy and possibly for T-coupling.*/
                        /* This may not be quite working correctly yet . . . . */
                        compute_globals(fplog, gstat, cr, ir, fr, ekind, state, state_global, mdatoms, nrnb, vcm,
                                        wcycle, enerd, NULL, NULL, NULL, NULL, mu_tot,
                                        constr, NULL, FALSE, state->box,
                                        top_global, &bSumEkinhOld,
                                        CGLO_RERUNMD | CGLO_GSTAT | CGLO_TEMPERATURE);
                        wallcycle_start(wcycle, ewcUPDATE);
                    }
                }
            }
            if (bTrotter && !bInitStep)
            {
                copy_mat(shake_vir, state->svir_prev);
                copy_mat(force_vir, state->fvir_prev);
                if (IR_NVT_TROTTER(ir) && ir->eI == eiVV)
                {
                    /* update temperature and kinetic energy now that step is over - this is the v(t+dt) point */
                    enerd->term[F_TEMP] = sum_ekin(&(ir->opts), ekind, NULL, (ir->eI == eiVV), FALSE);
                    enerd->term[F_EKIN] = trace(ekind->ekin);
                }
            }
            /* if it's the initial step, we performed this first step just to get the constraint virial */
            if (ir->eI == eiVV && bInitStep)
            {
                copy_rvecn(vbuf, state->v, 0, state->natoms);
                sfree(vbuf);
            }
            wallcycle_stop(wcycle, ewcUPDATE);
        }

        /* compute the conserved quantity */
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
            if (ir->efep != efepNO && !bRerunMD)
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

            lamnew = ExpandedEnsembleDynamics(fplog, ir, enerd, state, &MassQ, state->fep_state, &state->dfhist, step, state->v, mdatoms);
            /* history is maintained in state->dfhist, but state_global is what is sent to trajectory and log output */
            copy_df_history(&state_global->dfhist, &state->dfhist);
        }

        /* Now we have the energies and forces corresponding to the
         * coordinates at time t. We must output all of this before
         * the update.
         */
        do_md_trajectory_writing(fplog, cr, nfile, fnm, step, step_rel, t,
                                 ir, state, state_global, top_global, fr,
                                 outf, mdebin, ekind, f, f_global,
                                 &nchkpt,
                                 bCPT, bRerunMD, bLastStep, (Flags & MD_CONFOUT),
                                 bSumEkinhOld);
        /* Check if IMD step and do IMD communication, if bIMD is TRUE. */
        bIMDstep = do_IMD(ir->bIMD, step, cr, bNS, state->box, state->x, ir, t, wcycle);

        /* kludge -- virial is lost with restart for MTTK NPT control. Must reload (saved earlier). */
        if (bStartingFromCpt && bTrotter)
        {
            copy_mat(state->svir_prev, shake_vir);
            copy_mat(state->fvir_prev, force_vir);
        }

        elapsed_time = walltime_accounting_get_current_elapsed_time(walltime_accounting);

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
                 (max_hours > 0 && elapsed_time > max_hours*60.0*60.0*0.99) &&
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
            elapsed_time > max_hours*60.0*60.0*0.495)
        {
            /* Set flag that will communicate the signal to all ranks in the simulation */
            gs.sig[eglsRESETCOUNTERS] = 1;
        }

        /* In parallel we only have to check for checkpointing in steps
         * where we do global communication,
         *  otherwise the other nodes don't know.
         */
        if (MASTER(cr) && ((bGStat || !PAR(cr)) &&
                           cpt_period >= 0 &&
                           (cpt_period == 0 ||
                            elapsed_time >= nchkpt*cpt_period*60.0)) &&
            gs.set[eglsCHKPT] == 0)
        {
            gs.sig[eglsCHKPT] = 1;
        }

        /* at the start of step, randomize or scale the velocities ((if vv. Restriction of Andersen controlled
           in preprocessing */

        if (ETC_ANDERSEN(ir->etc)) /* keep this outside of update_tcouple because of the extra info required to pass */
        {
            gmx_bool bIfRandomize;
            bIfRandomize = update_randomize_velocities(ir, step, cr, mdatoms, state, upd, constr);
            /* if we have constraints, we have to remove the kinetic energy parallel to the bonds */
            if (constr && bIfRandomize)
            {
                update_constraints(fplog, step, NULL, ir, mdatoms,
                                   state, fr->bMolPBC, graph, f,
                                   &top->idef, tmp_vir,
                                   cr, nrnb, wcycle, upd, constr,
                                   TRUE, bCalcVir);
            }
        }
        /* #########   START SECOND UPDATE STEP ################# */
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
                update_pcouple(fplog, step, ir, state, pcoupl_mu, M, bInitStep);
            }

            if (bVV)
            {
                bUpdateDoLR = (fr->bTwinRange && do_per_step(step, ir->nstcalclr));

                /* velocity half-step update */
                update_coords(fplog, step, ir, mdatoms, state, fr->bMolPBC, f,
                              bUpdateDoLR, fr->f_twin, bCalcVir ? &fr->vir_twin_constr : NULL, fcd,
                              ekind, M, upd, FALSE, etrtVELOCITY2,
                              cr, nrnb, constr, &top->idef);
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
                copy_rvecn(state->x, cbuf, 0, state->natoms);
            }
            bUpdateDoLR = (fr->bTwinRange && do_per_step(step, ir->nstcalclr));

            update_coords(fplog, step, ir, mdatoms, state, fr->bMolPBC, f,
                          bUpdateDoLR, fr->f_twin, bCalcVir ? &fr->vir_twin_constr : NULL, fcd,
                          ekind, M, upd, bInitStep, etrtPOSITION, cr, nrnb, constr, &top->idef);
            wallcycle_stop(wcycle, ewcUPDATE);

            update_constraints(fplog, step, &dvdl_constr, ir, mdatoms, state,
                               fr->bMolPBC, graph, f,
                               &top->idef, shake_vir,
                               cr, nrnb, wcycle, upd, constr,
                               FALSE, bCalcVir);

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
                                top_global, &bSumEkinhOld,
                                cglo_flags | CGLO_TEMPERATURE
                                );
                wallcycle_start(wcycle, ewcUPDATE);
                trotter_update(ir, step, ekind, enerd, state, total_vir, mdatoms, &MassQ, trotter_seq, ettTSEQ4);
                /* now we know the scaling, we can compute the positions again again */
                copy_rvecn(cbuf, state->x, 0, state->natoms);

                bUpdateDoLR = (fr->bTwinRange && do_per_step(step, ir->nstcalclr));

                update_coords(fplog, step, ir, mdatoms, state, fr->bMolPBC, f,
                              bUpdateDoLR, fr->f_twin, bCalcVir ? &fr->vir_twin_constr : NULL, fcd,
                              ekind, M, upd, bInitStep, etrtPOSITION, cr, nrnb, constr, &top->idef);
                wallcycle_stop(wcycle, ewcUPDATE);

                /* do we need an extra constraint here? just need to copy out of state->v to upd->xp? */
                /* are the small terms in the shake_vir here due
                 * to numerical errors, or are they important
                 * physically? I'm thinking they are just errors, but not completely sure.
                 * For now, will call without actually constraining, constr=NULL*/
                update_constraints(fplog, step, NULL, ir, mdatoms,
                                   state, fr->bMolPBC, graph, f,
                                   &top->idef, tmp_vir,
                                   cr, nrnb, wcycle, upd, NULL,
                                   FALSE, bCalcVir);
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

        if (vsite != NULL)
        {
            wallcycle_start(wcycle, ewcVSITECONSTR);
            if (graph != NULL)
            {
                shift_self(graph, state->box, state->x);
            }
            construct_vsites(vsite, state->x, ir->delta_t, state->v,
                             top->idef.iparams, top->idef.il,
                             fr->ePBC, fr->bMolPBC, cr, state->box);

            if (graph != NULL)
            {
                unshift_self(graph, state->box, state->x);
            }
            wallcycle_stop(wcycle, ewcVSITECONSTR);
        }

        /* ############## IF NOT VV, Calculate globals HERE  ############ */
        /* With Leap-Frog we can skip compute_globals at
         * non-communication steps, but we need to calculate
         * the kinetic energy one step before communication.
         */
        if (bGStat || (!EI_VV(ir->eI) && do_per_step(step+1, nstglobalcomm)))
        {
            compute_globals(fplog, gstat, cr, ir, fr, ekind, state, state_global, mdatoms, nrnb, vcm,
                            wcycle, enerd, force_vir, shake_vir, total_vir, pres, mu_tot,
                            constr, &gs,
                            (step_rel % gs.nstms == 0) &&
                            (multisim_nsteps < 0 || (step_rel < multisim_nsteps)),
                            lastbox,
                            top_global, &bSumEkinhOld,
                            cglo_flags
                            | (!EI_VV(ir->eI) || bRerunMD ? CGLO_ENERGY : 0)
                            | (!EI_VV(ir->eI) && bStopCM ? CGLO_STOPCM : 0)
                            | (!EI_VV(ir->eI) ? CGLO_TEMPERATURE : 0)
                            | (!EI_VV(ir->eI) || bRerunMD ? CGLO_PRESSURE : 0)
                            | CGLO_CONSTRAINT
                            );
        }

        /* #############  END CALC EKIN AND PRESSURE ################# */

        /* Note: this is OK, but there are some numerical precision issues with using the convergence of
           the virial that should probably be addressed eventually. state->veta has better properies,
           but what we actually need entering the new cycle is the new shake_vir value. Ideally, we could
           generate the new shake_vir, but test the veta value for convergence.  This will take some thought. */

        if (ir->efep != efepNO && (!bVV || bRerunMD))
        {
            /* Sum up the foreign energy and dhdl terms for md and sd.
               Currently done every step so that dhdl is correct in the .edr */
            sum_dhdl(enerd, state->lambda, ir->fepvals);
        }
        update_box(fplog, step, ir, mdatoms, state, f,
                   pcoupl_mu, nrnb, upd);

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
        /* #########  END PREPARING EDR OUTPUT  ###########  */

        /* Output stuff */
        if (MASTER(cr))
        {
            if (fplog && do_log && bDoExpanded)
            {
                /* only needed if doing expanded ensemble */
                PrintFreeEnergyInfoToFile(fplog, ir->fepvals, ir->expandedvals, ir->bSimTemp ? ir->simtempvals : NULL,
                                          &state_global->dfhist, state->fep_state, ir->nstlog, step);
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

            print_ebin(mdoutf_get_fp_ene(outf), do_ene, do_dr, do_or, do_log ? fplog : NULL,
                       step, t,
                       eprNORMAL, bCompact, mdebin, fcd, groups, &(ir->opts));

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
                                             bRerunMD ? rerun_fr.x   : state->x,
                                             bRerunMD ? rerun_fr.box : state->box,
                                             top_global, MASTER(cr) && bVerbose, bRerunMD);

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
                                vsite, shellfc, constr,
                                nrnb, wcycle, FALSE);
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
            gs.set[eglsRESETCOUNTERS] != 0)
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
            reset_all_counters(fplog, cr, step, &step_rel, ir, wcycle, nrnb, walltime_accounting,
                               use_GPU(fr->nbv) ? fr->nbv : NULL);
            wcycle_set_reset_counters(wcycle, -1);
            if (!(cr->duty & DUTY_PME))
            {
                /* Tell our PME node to reset its counters */
                gmx_pme_send_resetcounters(cr, step);
            }
            /* Correct max_hours for the elapsed time */
            max_hours                -= elapsed_time/(60.0*60.0);
            /* If mdrun -maxh -resethway was active, it can only trigger once */
            bResetCountersHalfMaxH    = FALSE; /* TODO move this to where gs.sig[eglsRESETCOUNTERS] is set */
            /* Reset can only happen once, so clear the triggering flag. */
            gs.set[eglsRESETCOUNTERS] = 0;
        }

        /* If bIMD is TRUE, the master updates the IMD energy record and sends positions to VMD client */
        IMD_prep_energies_send_positions(ir->bIMD && MASTER(cr), bIMDstep, ir->imd, enerd, step, bCalcEner, wcycle);

    }
    /* End of main MD loop */
    debug_gmx();

    /* Closing TNG files can include compressing data. Therefore it is good to do that
     * before stopping the time measurements. */
    mdoutf_tng_close(outf);

    /* Stop measuring walltime */
    walltime_accounting_end(walltime_accounting);

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
            print_ebin(mdoutf_get_fp_ene(outf), FALSE, FALSE, FALSE, fplog, step, t,
                       eprAVER, FALSE, mdebin, fcd, groups, &(ir->opts));
        }
    }

    done_mdoutf(outf);
    debug_gmx();

    if (bPMETune)
    {
        pme_loadbal_done(pme_loadbal, cr, fplog, use_GPU(fr->nbv));
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

    /* IMD cleanup, if bIMD is TRUE. */
    IMD_finalize(ir->bIMD, ir->imd);

    walltime_accounting_set_nsteps_done(walltime_accounting, step_rel);

    return 0;
}
