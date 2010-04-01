#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <signal.h>
#include <stdlib.h>

#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
/* _isnan() */
#include <float.h>
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
#include "dihre.h"
#include "pppm.h"
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

#ifdef GMX_THREADS
#include "tmpi.h"
#endif

/* include even when OpenMM not used to force compilation of do_md_openmm */
#include "openmm_wrapper.h"


enum { eglsNABNSB, eglsCHKPT, eglsTERM, eglsRESETCOUNTERS, eglsNR };

typedef struct {
    int nstms;       /* The frequency for intersimulation communication */
    int sig[eglsNR]; /* The signal set by one process in do_md */
    int set[eglsNR]; /* The communicated signal, equal for all processes */
} globsig_t;


static int multisim_min(const gmx_multisim_t *ms,int nmin,int n)
{
    int  *buf;
    bool bPos,bEqual;
    int  s,d;

    snew(buf,ms->nsim);
    buf[ms->sim] = n;
    gmx_sumi_sim(ms->nsim,buf,ms);
    bPos   = TRUE;
    bEqual = TRUE;
    for(s=0; s<ms->nsim; s++)
    {
        bPos   = bPos   && (buf[s] > 0);
        bEqual = bEqual && (buf[s] == buf[0]);
    }
    if (bPos)
    {
        if (bEqual)
        {
            nmin = min(nmin,buf[0]);
        }
        else
        {
            /* Find the least common multiple */
            for(d=2; d<nmin; d++)
            {
                s = 0;
                while (s < ms->nsim && d % buf[s] == 0)
                {
                    s++;
                }
                if (s == ms->nsim)
                {
                    /* We found the LCM and it is less than nmin */
                    nmin = d;
                    break;
                }
            }
        }
    }
    sfree(buf);

    return nmin;
}

static int multisim_nstsimsync(const t_commrec *cr,
                               const t_inputrec *ir,int repl_ex_nst)
{
    int nmin;

    if (MASTER(cr))
    {
        nmin = INT_MAX;
        nmin = multisim_min(cr->ms,nmin,ir->nstlist);
        nmin = multisim_min(cr->ms,nmin,ir->nstcalcenergy);
        nmin = multisim_min(cr->ms,nmin,repl_ex_nst);
        if (nmin == INT_MAX)
        {
            gmx_fatal(FARGS,"Can not find an appropriate interval for inter-simulation communication, since nstlist, nstcalcenergy and -replex are all <= 0");
        }
        /* Avoid inter-simulation communication at every (second) step */
        if (nmin <= 2)
        {
            nmin = 10;
        }
    }

    gmx_bcast(sizeof(int),&nmin,cr);

    return nmin;
}

static void init_global_signals(globsig_t *gs,const t_commrec *cr,
                                const t_inputrec *ir,int repl_ex_nst)
{
    int i;

    gs->nstms = 1;

    for(i=0; i<eglsNR; i++)
    {
        gs->sig[i] = 0;
        gs->set[i] = 0;
    }
}


double do_md_openmm(FILE *fplog,t_commrec *cr,int nfile,const t_filenm fnm[],
             const output_env_t oenv, bool bVerbose,bool bCompact,
             int nstglobalcomm,
             gmx_vsite_t *vsite,gmx_constr_t constr,
             int stepout,t_inputrec *ir,
             gmx_mtop_t *top_global,
             t_fcdata *fcd,
             t_state *state_global,
             t_mdatoms *mdatoms,
             t_nrnb *nrnb,gmx_wallcycle_t wcycle,
             gmx_edsam_t ed,t_forcerec *fr,
             int repl_ex_nst,int repl_ex_seed,
             real cpt_period,real max_hours,
             const char *deviceOptions,
             unsigned long Flags,
             gmx_runtime_t *runtime)
{
    gmx_mdoutf_t *outf;
    gmx_large_int_t step,step_rel;
    double     run_time;
    double     t,t0,lam0;
    bool       bSimAnn,
               bFirstStep,bStateFromTPX,bLastStep,bStartingFromCpt;
    bool       bInitStep=TRUE;
    bool       bNEMD,do_ene,do_log, do_verbose,
               bX,bV,bF,bXTC,bCPT;
    tensor     force_vir,shake_vir,total_vir,pres;
    int        i,m;
    int        mdof_flags;
    rvec       mu_tot;
    t_vcm      *vcm;
    int        nchkpt=1;
    gmx_localtop_t *top;	
    t_mdebin *mdebin=NULL;
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
    bool        bAppend;
    int         a0,a1;
    matrix      lastbox;
    real        reset_counters=0,reset_counters_now=0;
    char        sbuf[STEPSTRSIZE],sbuf2[STEPSTRSIZE];
    int         handledSignal=-1; /* compare to last_signal_recvd */
    bool        bHandledSignal=FALSE;
 
    const char *ommOptions = NULL;
    void   *openmmData;

    bAppend  = (Flags & MD_APPENDFILES);
    check_ir_old_tpx_versions(cr,fplog,ir,top_global);

    groups = &top_global->groups;

    /* Initial values */
    init_md(fplog,cr,ir,oenv,&t,&t0,&state_global->lambda,&lam0,
            nrnb,top_global,&upd,
            nfile,fnm,&outf,&mdebin,
            force_vir,shake_vir,mu_tot,&bNEMD,&bSimAnn,&vcm,state_global,Flags);

    clear_mat(total_vir);
    clear_mat(pres);
    /* Energy terms and groups */
    snew(enerd,1);
    init_enerdata(top_global->groups.grps[egcENER].nr,ir->n_flambda,enerd);
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

        if (vsite) {
            set_vsite_top(vsite,top,mdatoms,cr);
        }

        if (ir->ePBC != epbcNONE && !ir->bPeriodicMols) {
            graph = mk_graph(fplog,&(top->idef),0,top_global->natoms,FALSE,FALSE);
        }

    update_mdatoms(mdatoms,state->lambda);

                                                                                                                                                                                              
      if(deviceOptions[0]=='\0')                                                                                                                                                              
      {                                                                                                                                                                                       
          /* empty options, which should default to OpenMM in this build */                                                                                                                   
          ommOptions=deviceOptions;                                                                                                                                                           
      }                                                                                                                                                                                       
      else                                                                                                                                                                                    
      {                                                                                                                                                                                       
          if(gmx_strncasecmp(deviceOptions,"OpenMM",6)!=0)                                                                                                                                    
          {                                                                                                                                                                                   
              gmx_fatal(FARGS, "This Gromacs version currently only works with OpenMM. Use -device \"OpenMM:<options>\"");                                                                    
          }                                                                                                                                                                                   
          else                                                                                                                                                                                
          {                                                                                                                                                                                   
              ommOptions=strchr(deviceOptions,':');                                                                                                                                           
              if(NULL!=ommOptions)                                                                                                                                                            
              {                                                                                                                                                                               
                  /* Increase the pointer to skip the colon */                                                                                                                                
                  ommOptions++;                                     
            }                                                                                                                                                                               
          }                                                                                                                                                                                   
      }                                                                                                                                                                                       

      openmmData = openmm_init(fplog, ommOptions, cr, ir, top_global, top, mdatoms, fr, state);                                                                                         


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

    if (constr) {
            set_constraints(constr,top,ir,mdatoms,cr);
    }

    if (!ir->bContinuation)
    {
        if (mdatoms->cFREEZE && (state->flags & (1<<estV)))
        {
            /* Set the velocities of frozen particles to zero */
            for(i=mdatoms->start; i<mdatoms->start+mdatoms->homenr; i++)
            {
                for(m=0; m<DIM; m++)
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
        fprintf(fplog,"Initial temperature: %g K\n",enerd->term[F_TEMP]);
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

    bLastStep = ((ir->nsteps >= 0 && step_rel > ir->nsteps));

    while (!bLastStep) {

        wallcycle_start(wcycle,ewcSTEP);

        GMX_MPE_LOG(ev_timestep1);

        bLastStep = (step_rel == ir->nsteps);
        t = t0 + step*ir->delta_t;

        if (gs.set[eglsTERM] != 0 )
        {
            bLastStep = TRUE;
        }

        do_log = do_per_step(step,ir->nstlog) || bFirstStep || bLastStep;
        do_verbose = bVerbose &&
                  (step % stepout == 0 || bFirstStep || bLastStep);

        if (MASTER(cr) && do_log)
        {
            print_ebin_header(fplog,step,t,state->lambda);
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
        if (do_per_step(step,ir->nstxout)) { mdof_flags |= MDOF_X; }
        if (do_per_step(step,ir->nstvout)) { mdof_flags |= MDOF_V; }
        if (do_per_step(step,ir->nstfout)) { mdof_flags |= MDOF_F; }
        if (do_per_step(step,ir->nstxtcout)) { mdof_flags |= MDOF_XTC; }
        if (bCPT) { mdof_flags |= MDOF_CPT; };
        bX   = MDOF_X || (bLastStep && ir->nstxout);
        bV   = MDOF_V || (bLastStep && ir->nstvout);
        bF   = MDOF_F || (bLastStep && ir->nstfout);
        bXTC = MDOF_XTC || (bLastStep && ir->nstxtcout);
        do_ene = do_per_step(step,ir->nstenergy) || (bLastStep && ir->nstenergy);

      if( bX || bXTC || bV ){                                                                                                                                                                 
        wallcycle_start(wcycle,ewcTRAJ);                                                                                                                                                              
        openmm_copy_state(openmmData, state, &t, f, enerd, bX||bXTC, bV, 0, 0);                                                                                                               
        wallcycle_stop(wcycle,ewcTRAJ);                                                                                                                                                       
      }                                                                                                                                                                                       

      openmm_take_one_step(openmmData);                                                                                                                                                       

      if (bX || bV || bF || bXTC || do_ene) {
        wallcycle_start(wcycle,ewcTRAJ);
        if( bF || do_ene ){                                                                                                                                                                   
           openmm_copy_state(openmmData, state, &t, f, enerd, 0, 0, bF, do_ene);                                                                                                              
      }                                                                                                                                                                                     
     upd_mdebin(mdebin, NULL,TRUE,
		 t,mdatoms->tmass,enerd,state,lastbox,
		 shake_vir,force_vir,total_vir,pres,
		 ekind,mu_tot,constr);
         print_ebin(outf->fp_ene,do_ene,FALSE,FALSE,fplog,step,t,
                       eprAVER,FALSE,mdebin,fcd,groups,&(ir->opts));

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
        if ((bGotStopNextStepSignal || bGotStopNextNSStepSignal) && 
            (handledSignal!=last_signal_number_recvd) &&
            MASTERTHREAD(cr))
        {
            if (bGotStopNextStepSignal)
            {
                gs.set[eglsTERM] = 1;
            }
            else
            {
                gs.set[eglsTERM] = -1;
            }
            if (fplog)
            {
                fprintf(fplog,
                        "\n\nReceived the %s signal, stopping at the next %sstep\n\n",
                        signal_name[last_signal_number_recvd], 
                        gs.set[eglsTERM]==-1 ? "NS " : "");
                fflush(fplog);
            }
            fprintf(stderr,
                    "\n\nReceived the %s signal, stopping at the next %sstep\n\n",
                    signal_name[last_signal_number_recvd], 
                    gs.set[eglsTERM]==-1 ? "NS " : "");
            fflush(stderr);
            handledSignal=last_signal_number_recvd;
        }
        else if (MASTER(cr) &&
                 (max_hours > 0 && run_time > max_hours*60.0*60.0*0.99) &&
                 gs.set[eglsTERM] == 0)
        {
            /* Signal to terminate the run */
            gs.set[eglsTERM] = 1;
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
                if(fflush(fplog) != 0)
                {
                    gmx_fatal(FARGS,"Cannot flush logfile - maybe you are out of quota?");
                }
            }


        /* Remaining runtime */
        if (MULTIMASTER(cr) && do_verbose) 
        {
            print_time(stderr,runtime,step,ir,cr);
        }


        bFirstStep = FALSE;
        bInitStep = FALSE;
        bStartingFromCpt = FALSE;
        step++;
        step_rel++;



    }
    /* End of main MD loop */
    debug_gmx();
    
    openmm_cleanup(fplog, openmmData);

    /* Stop the time */
    runtime_end(runtime);
    

    done_mdoutf(outf);

    debug_gmx();

    runtime->nsteps_done = step_rel;
    
    return 0;
}
