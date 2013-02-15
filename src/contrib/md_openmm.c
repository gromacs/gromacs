/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2010, The GROMACS development team,
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

#include <signal.h>
#include <stdlib.h>

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
#include "qmmm.h"
#include "mpelogging.h"
#include "domdec.h"
#include "partdec.h"
#include "topsort.h"
#include "coulomb.h"
#include "constr.h"
#include "compute_io.h"
#include "mvdata.h"
#include "checkpoint.h"
#include "mtop_util.h"
#include "sighandler.h"
#include "genborn.h"
#include "string2.h"
#include "copyrite.h"
#include "membed.h"

#ifdef GMX_THREAD_MPI
#include "tmpi.h"
#endif

/* include even when OpenMM not used to force compilation of do_md_openmm */
#include "openmm_wrapper.h"

double do_md_openmm(FILE *fplog,t_commrec *cr,int nfile,const t_filenm fnm[],
                    const output_env_t oenv, gmx_bool bVerbose,gmx_bool bCompact,
                    int nstglobalcomm,
                    gmx_vsite_t *vsite,gmx_constr_t constr,
                    int stepout,t_inputrec *ir,
                    gmx_mtop_t *top_global,
                    t_fcdata *fcd,
                    t_state *state_global,
                    t_mdatoms *mdatoms,
                    t_nrnb *nrnb,gmx_wallcycle_t wcycle,
                    gmx_edsam_t ed,t_forcerec *fr,
                    int repl_ex_nst, int repl_ex_nex, int repl_ex_seed,
                    gmx_membed_t membed,
                    real cpt_period,real max_hours,
                    const char *deviceOptions,
                    unsigned long Flags,
                    gmx_runtime_t *runtime)
{
    gmx_mdoutf_t *outf;
    gmx_large_int_t step,step_rel;
    double     run_time;
    double     t,t0,lam0;
    gmx_bool       bSimAnn,
    bFirstStep,bStateFromTPX,bLastStep,bStartingFromCpt;
    gmx_bool       bInitStep=TRUE;
    gmx_bool       do_ene,do_log, do_verbose,
    bX,bV,bF,bCPT;
    tensor     force_vir,shake_vir,total_vir,pres;
    int        i,m;
    int        mdof_flags;
    rvec       mu_tot;
    t_vcm      *vcm;
    int        nchkpt=1;
    gmx_localtop_t *top;
    t_mdebin *mdebin;
    t_state    *state=NULL;
    rvec       *f_global=NULL;
    int        n_xtc=-1;
    rvec       *x_xtc=NULL;
    gmx_enerdata_t *enerd;
    rvec       *f=NULL;
    gmx_global_stat_t gstat;
    gmx_update_t upd=NULL;
    t_graph    *graph=NULL;
    globsig_t   gs;

    gmx_groups_t *groups;
    gmx_ekindata_t *ekind, *ekind_save;
    gmx_bool        bAppend;
    int         a0,a1;
    matrix      lastbox;
    real        reset_counters=0,reset_counters_now=0;
    char        sbuf[STEPSTRSIZE],sbuf2[STEPSTRSIZE];
    int         handled_stop_condition=gmx_stop_cond_none; 

    const char *ommOptions = NULL;
    void   *openmmData;

#ifdef GMX_DOUBLE
    /* Checks in cmake should prevent the compilation in double precision
     * with OpenMM, but just to be sure we check here.
     */
    gmx_fatal(FARGS,"Compilation was performed in double precision, but OpenMM only supports single precision. If you want to use to OpenMM, compile in single precision.");
#endif

    bAppend  = (Flags & MD_APPENDFILES);
    check_ir_old_tpx_versions(cr,fplog,ir,top_global);

    groups = &top_global->groups;

    /* Initial values */
    init_md(fplog,cr,ir,oenv,&t,&t0,state_global->lambda,
            &(state_global->fep_state),&lam0,
            nrnb,top_global,&upd,
            nfile,fnm,&outf,&mdebin,
            force_vir,shake_vir,mu_tot,&bSimAnn,&vcm,state_global,Flags);

    clear_mat(total_vir);
    clear_mat(pres);
    /* Energy terms and groups */
    snew(enerd,1);
    init_enerdata(top_global->groups.grps[egcENER].nr,ir->fepvals->n_lambda,
                  enerd);
    snew(f,top_global->natoms);

    /* Kinetic energy data */
    snew(ekind,1);
    init_ekindata(fplog,top_global,&(ir->opts),ekind);
    /* needed for iteration of constraints */
    snew(ekind_save,1);
    init_ekindata(fplog,top_global,&(ir->opts),ekind_save);
    /* Copy the cos acceleration to the groups struct */
    ekind->cosacc.cos_accel = ir->cos_accel;

    gstat = global_stat_init(ir);
    debug_gmx();

    {
        double io = compute_io(ir,top_global->natoms,groups,mdebin->ebin->nener,1);
        if ((io > 2000) && MASTER(cr))
            fprintf(stderr,
                    "\nWARNING: This run will generate roughly %.0f Mb of data\n\n",
                    io);
    }

    top = gmx_mtop_generate_local_top(top_global,ir);

    a0 = 0;
    a1 = top_global->natoms;

    state = partdec_init_local_state(cr,state_global);
    f_global = f;

    atoms2md(top_global,ir,0,NULL,a0,a1-a0,mdatoms);

    if (vsite)
    {
        set_vsite_top(vsite,top,mdatoms,cr);
    }

    if (ir->ePBC != epbcNONE && !ir->bPeriodicMols)
    {
        graph = mk_graph(fplog,&(top->idef),0,top_global->natoms,FALSE,FALSE);
    }

    update_mdatoms(mdatoms,state->lambda[efptMASS]);

    if (deviceOptions[0]=='\0')
    {
        /* empty options, which should default to OpenMM in this build */
        ommOptions=deviceOptions;
    }
    else
    {
        if (gmx_strncasecmp(deviceOptions,"OpenMM",6)!=0)
        {
            gmx_fatal(FARGS, "This Gromacs version currently only works with OpenMM. Use -device \"OpenMM:<options>\"");
        }
        else
        {
            ommOptions=strchr(deviceOptions,':');
            if (NULL!=ommOptions)
            {
                /* Increase the pointer to skip the colon */
                ommOptions++;
            }
        }
    }

    openmmData = openmm_init(fplog, ommOptions, ir, top_global, top, mdatoms, fr, state);
    please_cite(fplog,"Friedrichs2009");

    if (MASTER(cr))
    {
        /* Update mdebin with energy history if appending to output files */
        if ( Flags & MD_APPENDFILES )
        {
            restore_energyhistory_from_state(mdebin,&state_global->enerhist);
        }
        /* Set the initial energy history in state to zero by updating once */
        update_energyhistory(&state_global->enerhist,mdebin);
    }

    if (constr)
    {
        set_constraints(constr,top,ir,mdatoms,cr);
    }

    if (!ir->bContinuation)
    {
        if (mdatoms->cFREEZE && (state->flags & (1<<estV)))
        {
            /* Set the velocities of frozen particles to zero */
            for (i=mdatoms->start; i<mdatoms->start+mdatoms->homenr; i++)
            {
                for (m=0; m<DIM; m++)
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
            do_constrain_first(fplog,constr,ir,mdatoms,state,f,
                               graph,cr,nrnb,fr,top,shake_vir);
        }
        if (vsite)
        {
            /* Construct the virtual sites for the initial configuration */
            construct_vsites(fplog,vsite,state->x,nrnb,ir->delta_t,NULL,
                             top->idef.iparams,top->idef.il,
                             fr->ePBC,fr->bMolPBC,graph,cr,state->box);
        }
    }

    debug_gmx();

    if (MASTER(cr))
    {
        char tbuf[20];
        fprintf(stderr,"starting mdrun '%s'\n",
                *(top_global->name));
        if (ir->nsteps >= 0)
        {
            sprintf(tbuf,"%8.1f",(ir->init_step+ir->nsteps)*ir->delta_t);
        }
        else
        {
            sprintf(tbuf,"%s","infinite");
        }
        if (ir->init_step > 0)
        {
            fprintf(stderr,"%s steps, %s ps (continuing from step %s, %8.1f ps).\n",
                    gmx_step_str(ir->init_step+ir->nsteps,sbuf),tbuf,
                    gmx_step_str(ir->init_step,sbuf2),
                    ir->init_step*ir->delta_t);
        }
        else
        {
            fprintf(stderr,"%s steps, %s ps.\n",
                    gmx_step_str(ir->nsteps,sbuf),tbuf);
        }
    }

    fprintf(fplog,"\n");

    /* Set and write start time */
    runtime_start(runtime);
    print_date_and_time(fplog,cr->nodeid,"Started mdrun",runtime);
    wallcycle_start(wcycle,ewcRUN);
    if (fplog)
        fprintf(fplog,"\n");

    /* safest point to do file checkpointing is here.  More general point would be immediately before integrator call */

    debug_gmx();
    /***********************************************************
     *
     *             Loop over MD steps
     *
     ************************************************************/

    /* loop over MD steps or if rerunMD to end of input trajectory */
    bFirstStep = TRUE;
    /* Skip the first Nose-Hoover integration when we get the state from tpx */
    bStateFromTPX = !opt2bSet("-cpi",nfile,fnm);
    bInitStep = bFirstStep && bStateFromTPX;
    bStartingFromCpt = (Flags & MD_STARTFROMCPT) && bInitStep;
    bLastStep = FALSE;

    init_global_signals(&gs,cr,ir,repl_ex_nst);

    step = ir->init_step;
    step_rel = 0;

    while (!bLastStep)
    {
        wallcycle_start(wcycle,ewcSTEP);

        GMX_MPE_LOG(ev_timestep1);

        bLastStep = (step_rel == ir->nsteps);
        t = t0 + step*ir->delta_t;

        if (gs.set[eglsSTOPCOND] != 0)
        {
            bLastStep = TRUE;
        }

        do_log = do_per_step(step,ir->nstlog) || bFirstStep || bLastStep;
        do_verbose = bVerbose &&
                     (step % stepout == 0 || bFirstStep || bLastStep);

        if (MASTER(cr) && do_log)
        {
            print_ebin_header(fplog,step,t,state->lambda[efptFEP]);
        }

        clear_mat(force_vir);
        GMX_MPE_LOG(ev_timestep2);

        /* We write a checkpoint at this MD step when:
         * either when we signalled through gs (in OpenMM NS works different),
         * or at the last step (but not when we do not want confout),
         * but never at the first step.
         */
        bCPT = ((gs.set[eglsCHKPT] ||
                 (bLastStep && (Flags & MD_CONFOUT))) &&
                step > ir->init_step );
        if (bCPT)
        {
            gs.set[eglsCHKPT] = 0;
        }

        /* Now we have the energies and forces corresponding to the
         * coordinates at time t. We must output all of this before
         * the update.
         * for RerunMD t is read from input trajectory
         */
        GMX_MPE_LOG(ev_output_start);

        mdof_flags = 0;
        if (do_per_step(step,ir->nstxout))
        {
            mdof_flags |= MDOF_X;
        }
        if (do_per_step(step,ir->nstvout))
        {
            mdof_flags |= MDOF_V;
        }
        if (do_per_step(step,ir->nstfout))
        {
            mdof_flags |= MDOF_F;
        }
        if (do_per_step(step,ir->nstxtcout))
        {
            mdof_flags |= MDOF_XTC;
        }
        if (bCPT)
        {
            mdof_flags |= MDOF_CPT;
        };
        do_ene = (do_per_step(step,ir->nstenergy) || bLastStep);

        if (mdof_flags != 0 || do_ene || do_log)
        {
            wallcycle_start(wcycle,ewcTRAJ);
            bF = (mdof_flags & MDOF_F);
            bX = (mdof_flags & (MDOF_X | MDOF_XTC | MDOF_CPT));
            bV = (mdof_flags & (MDOF_V | MDOF_CPT));

            openmm_copy_state(openmmData, state, &t, f, enerd, bX, bV, bF, do_ene);

            upd_mdebin(mdebin,FALSE,TRUE,
                       t,mdatoms->tmass,enerd,state,ir->fepvals,ir->expandedvals,lastbox,
                       shake_vir,force_vir,total_vir,pres,
                       ekind,mu_tot,constr);
            print_ebin(outf->fp_ene,do_ene,FALSE,FALSE,do_log?fplog:NULL,
                       step,t,
                       eprNORMAL,bCompact,mdebin,fcd,groups,&(ir->opts));
            write_traj(fplog,cr,outf,mdof_flags,top_global,
                       step,t,state,state_global,f,f_global,&n_xtc,&x_xtc);
            if (bCPT)
            {
                nchkpt++;
                bCPT = FALSE;
            }
            debug_gmx();
            if (bLastStep && step_rel == ir->nsteps &&
                    (Flags & MD_CONFOUT) && MASTER(cr))
            {
                /* x and v have been collected in write_traj,
                 * because a checkpoint file will always be written
                 * at the last step.
                 */
                fprintf(stderr,"\nWriting final coordinates.\n");
                if (ir->ePBC != epbcNONE && !ir->bPeriodicMols)
                {
                    /* Make molecules whole only for confout writing */
                    do_pbc_mtop(fplog,ir->ePBC,state->box,top_global,state_global->x);
                }
                write_sto_conf_mtop(ftp2fn(efSTO,nfile,fnm),
                                    *top_global->name,top_global,
                                    state_global->x,state_global->v,
                                    ir->ePBC,state->box);
                debug_gmx();
            }
            wallcycle_stop(wcycle,ewcTRAJ);
        }
        GMX_MPE_LOG(ev_output_finish);


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
            /* NOTE: this only works for serial code. For code that allows
               MPI nodes to propagate their condition, see kernel/md.c*/
            if ( gmx_get_stop_condition() == gmx_stop_cond_next_ns )
                gs.set[eglsSTOPCOND]=1;
            if ( gmx_get_stop_condition() == gmx_stop_cond_next )
                gs.set[eglsSTOPCOND]=1;
            /* < 0 means stop at next step, > 0 means stop at next NS step */
            if (fplog)
            {
                fprintf(fplog,
                        "\n\nReceived the %s signal, stopping at the next %sstep\n\n",
                        gmx_get_signal_name(),
                        gs.sig[eglsSTOPCOND]==1 ? "NS " : "");
                fflush(fplog);
            }
            fprintf(stderr,
                    "\n\nReceived the %s signal, stopping at the next %sstep\n\n",
                    gmx_get_signal_name(),
                    gs.sig[eglsSTOPCOND]==1 ? "NS " : "");
            fflush(stderr);
            handled_stop_condition=(int)gmx_get_stop_condition();
        }
        else if (MASTER(cr) &&
                 (max_hours > 0 && run_time > max_hours*60.0*60.0*0.99) &&
                 gs.set[eglsSTOPCOND] == 0)
        {
            /* Signal to terminate the run */
            gs.set[eglsSTOPCOND] = 1;
            if (fplog)
            {
                fprintf(fplog,"\nStep %s: Run time exceeded %.3f hours, will terminate the run\n",gmx_step_str(step,sbuf),max_hours*0.99);
            }
            fprintf(stderr, "\nStep %s: Run time exceeded %.3f hours, will terminate the run\n",gmx_step_str(step,sbuf),max_hours*0.99);
        }

        /* checkpoints */
        if (MASTER(cr) && (cpt_period >= 0 &&
                           (cpt_period == 0 ||
                            run_time >= nchkpt*cpt_period*60.0)) &&
                gs.set[eglsCHKPT] == 0)
        {
            gs.set[eglsCHKPT] = 1;
        }

        /* Time for performance */
        if (((step % stepout) == 0) || bLastStep)
        {
            runtime_upd_proc(runtime);
        }

        if (do_per_step(step,ir->nstlog))
        {
            if (fflush(fplog) != 0)
            {
                gmx_fatal(FARGS,"Cannot flush logfile - maybe you are out of disk space?");
            }
        }

        /* Remaining runtime */
        if (MULTIMASTER(cr) && (do_verbose || gmx_got_usr_signal() ))
        {
            print_time(stderr,runtime,step,ir,cr);
        }

        bFirstStep = FALSE;
        bInitStep = FALSE;
        bStartingFromCpt = FALSE;
        step++;
        step_rel++;

        openmm_take_one_step(openmmData);
    }
    /* End of main MD loop */
    debug_gmx();

    /* Stop the time */
    runtime_end(runtime);

    if (MASTER(cr))
    {
        if (ir->nstcalcenergy > 0) 
        {
            print_ebin(outf->fp_ene,FALSE,FALSE,FALSE,fplog,step,t,
                       eprAVER,FALSE,mdebin,fcd,groups,&(ir->opts));
        }
    }

    openmm_cleanup(fplog, openmmData);

    done_mdoutf(outf);

    debug_gmx();

    runtime->nsteps_done = step_rel;

    return 0;
}
