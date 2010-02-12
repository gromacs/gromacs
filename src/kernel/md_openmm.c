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
#include "repl_ex.h"
#include "qmmm.h"
#include "mpelogging.h"
#include "domdec.h"
#include "partdec.h"
#include "topsort.h"
#include "coulomb.h"
#include "constr.h"
#include "shellfc.h"
#include "compute_io.h"
#include "mvdata.h"
#include "checkpoint.h"
#include "mtop_util.h"
#include "genborn.h"
#include "string2.h"

#ifdef GMX_THREADS
#include "tmpi.h"
#endif

/* include even when OpenMM not used to force compilation of do_md_openmm */
#include "openmm_wrapper.h"

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
    int        fp_trn=0,fp_xtc=0;
    ener_file_t fp_ene=NULL;
    gmx_large_int_t step,step_rel;
    const char *fn_cpt;
    FILE       *fp_dhdl=NULL,*fp_field=NULL;
    double     run_time;
    double     t,t0,lam0;
    bool       bGStatEveryStep,bGStat,bNstEner,bCalcPres,bCalcEner;
    bool       bNS,bNStList,bSimAnn,bStopCM,bRerunMD,bNotLastFrame=FALSE,
               bFirstStep,bStateFromTPX,bInitStep,bLastStep,
               bBornRadii,bStartingFromCpt;
    bool       bDoDHDL=FALSE;
    bool       bNEMD,do_ene,do_log,do_verbose,bRerunWarnNoV=TRUE,
               bForceUpdate=FALSE,bX,bV,bF,bXTC,bCPT=FALSE;
    bool       bMasterState;
    int        force_flags,cglo_flags;
    tensor     force_vir,shake_vir,total_vir,tmp_vir,pres;
    int        i,m,status;
    rvec       mu_tot;
    t_vcm      *vcm;
    t_state    *bufstate=NULL;   
    matrix     *scale_tot,pcoupl_mu,M,ebox;
//    gmx_nlheur_t nlh;
    t_trxframe rerun_fr;
    gmx_repl_ex_t repl_ex=NULL;
    int        nchkpt=1;
    /* Booleans (disguised as a reals) to checkpoint and terminate mdrun */  
    real       chkpt=0,terminate=0,terminate_now=0;

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

    bool        bFFscan;
    gmx_groups_t *groups;
    gmx_ekindata_t *ekind, *ekind_save;
    gmx_shellfc_t shellfc;
    int         count,nconverged=0;
    real        timestep=0;
    double      tcount=0;
    bool        bIonize=FALSE;
    bool        bTCR=FALSE,bConverged=TRUE,bOK,bSumEkinhOld,bExchanged;
    bool        bAppend;
    bool        bResetCountersHalfMaxH=FALSE;
    bool        bVV,bIterations,bIterate,bFirstIterate,bTemp,bPres,bTrotter;
    real        temp0,mu_aver=0,dvdl;
    int         a0,a1,gnx=0,ii;
    atom_id     *grpindex=NULL;
    char        *grpname;
    t_coupl_rec *tcr=NULL;
    rvec        *xcopy=NULL,*vcopy=NULL,*cbuf=NULL;
    matrix      boxcopy={{0}},lastbox;
	tensor      tmpvir;
	real        fom,oldfom,veta_save,pcurr,scalevir,tracevir;
	real        vetanew = 0;
    double      cycles;
    real        reset_counters=0,reset_counters_now=0;
	real        last_conserved = 0;
    real        last_ekin = 0;
	int         iter_i;
	t_extmass   MassQ;
    int         **trotter_seq; 
    char        sbuf[22],sbuf2[22];
    bool        bHandledSignal=FALSE;
//    gmx_iterate_t iterate;
 
    const char *ommOptions = NULL;
    void   *openmmData;

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
            wcycle_set_reset_counters(wcycle,ir->nsteps/2);
        }
        /* Signal to reset the counters halfway the simulation time. */
        bResetCountersHalfMaxH = (max_hours > 0);
    }

    /* md-vv uses averaged full step velocities for T-control 
       md-vv2 uses averaged half step velocities for T-control (but full step ekin for P control)
       md uses averaged half step kinetic energies to determine temperature unless defined otherwise by GMX_EKIN_AVE_VEL; */
    bVV = EI_VV(ir->eI);
    if (bVV) /* to store the initial velocities while computing virial */
    {
        snew(cbuf,top_global->natoms);
    }
    /* all the iteratative cases - only if there are constraints */ 
    bIterations = ((IR_NPT_TROTTER(ir)) && (constr) && (!bRerunMD));
    bTrotter = (bVV && (IR_NPT_TROTTER(ir) || (IR_NVT_TROTTER(ir))));        
    
    if (bRerunMD)
    {
        /* Since we don't know if the frames read are related in any way,
         * rebuild the neighborlist at every step.
         */
        ir->nstlist       = 1;
        ir->nstcalcenergy = 1;
        nstglobalcomm     = 1;
    }

    check_ir_old_tpx_versions(cr,fplog,ir,top_global);

//****    nstglobalcomm = check_nstglobalcomm(fplog,cr,nstglobalcomm,ir);
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
    init_md(fplog,cr,ir,oenv,&t,&t0,&state_global->lambda,&lam0,
            nrnb,top_global,&upd,
            nfile,fnm,&fp_trn,&fp_xtc,&fp_ene,&fn_cpt,
            &fp_dhdl,&fp_field,&mdebin,
            force_vir,shake_vir,mu_tot,&bNEMD,&bSimAnn,&vcm,state_global,Flags);

    clear_mat(total_vir);
    clear_mat(pres);
    /* Energy terms and groups */
    snew(enerd,1);
    init_enerdata(top_global->groups.grps[egcENER].nr,ir->n_flambda,enerd);
    if (DOMAINDECOMP(cr))
    {
        f = NULL;
    }
    else
    {
        snew(f,top_global->natoms);
    }

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

    /* Check for polarizable models and flexible constraints */
    shellfc = init_shell_flexcon(fplog,
                                 top_global,n_flexible_constraints(constr),
                                 (ir->bContinuation || 
                                  (DOMAINDECOMP(cr) && !MASTER(cr))) ?
                                 NULL : state_global->x);

    if (DEFORM(*ir))
    {
#ifdef GMX_THREADS
        tMPI_Thread_mutex_lock(&deform_init_box_mutex);
#endif
        set_deform_reference_box(upd,
                                 deform_init_init_step_tpx,
                                 deform_init_box_tpx);
#ifdef GMX_THREADS
        tMPI_Thread_mutex_unlock(&deform_init_box_mutex);
#endif
    }

    {
        double io = compute_io(ir,top_global->natoms,groups,mdebin->ebin->nener,1);
        if ((io > 2000) && MASTER(cr))
            fprintf(stderr,
                    "\nWARNING: This run will generate roughly %.0f Mb of data\n\n",
                    io);
    }

    if (DOMAINDECOMP(cr)) {
        top = dd_init_local_top(top_global);

        snew(state,1);
        dd_init_local_state(cr->dd,state_global,state);

        if (DDMASTER(cr->dd) && ir->nstfout) {
            snew(f_global,state_global->natoms);
        }
    } else {
        if (PAR(cr)) {
            /* Initialize the particle decomposition and split the topology */
            top = split_system(fplog,top_global,ir,cr);

            pd_cg_range(cr,&fr->cg0,&fr->hcg);
            pd_at_range(cr,&a0,&a1);
        } else {
            top = gmx_mtop_generate_local_top(top_global,ir);

            a0 = 0;
            a1 = top_global->natoms;
        }

        state = partdec_init_local_state(cr,state_global);
        f_global = f;

        atoms2md(top_global,ir,0,NULL,a0,a1-a0,mdatoms);

        if (vsite) {
            set_vsite_top(vsite,top,mdatoms,cr);
        }

        if (ir->ePBC != epbcNONE && !ir->bPeriodicMols) {
            graph = mk_graph(fplog,&(top->idef),0,top_global->natoms,FALSE,FALSE);
        }

        if (shellfc) {
            make_local_shells(cr,mdatoms,shellfc);
        }

        if (ir->pull && PAR(cr)) {
            dd_make_local_pull_groups(NULL,ir->pull,mdatoms);
        }
    }

    if (DOMAINDECOMP(cr))
    {
        /* Distribute the charge groups over the nodes from the master node */
        dd_partition_system(fplog,ir->init_step,cr,TRUE,1,
                            state_global,top_global,ir,
                            state,&f,mdatoms,top,fr,
                            vsite,shellfc,constr,
                            nrnb,wcycle,FALSE);
    }

    /* If not DD, copy gb data */
    if(ir->implicit_solvent && !DOMAINDECOMP(cr))
    {
        make_local_gb(cr,fr->born,ir->gb_algorithm);
    }

    update_mdatoms(mdatoms,state->lambda);
//#ifdef GMX_OPENMM                                                                                                                                                                           
                                                                                                                                                                                              
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
//#endif             

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

    if ((state->flags & (1<<estLD_RNG)) && (Flags & MD_READ_RNG)) {
        /* Set the random state if we read a checkpoint file */
        set_stochd_state(upd,state);
    }

    /* Initialize constraints */
    if (constr) {
        if (!DOMAINDECOMP(cr))
            set_constraints(constr,top,ir,mdatoms,cr);
    }

    /* Check whether we have to GCT stuff */
    bTCR = ftp2bSet(efGCT,nfile,fnm);
    if (bTCR) {
        if (MASTER(cr)) {
            fprintf(stderr,"Will do General Coupling Theory!\n");
        }
        gnx = top_global->mols.nr;
        snew(grpindex,gnx);
        for(i=0; (i<gnx); i++) {
            grpindex[i] = i;
        }
    }

    if (repl_ex_nst > 0 && MASTER(cr))
        repl_ex = init_replica_exchange(fplog,cr->ms,state_global,ir,
                                        repl_ex_nst,repl_ex_seed);

    if (!ir->bContinuation && !bRerunMD)
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
  
    /* I'm assuming we need global communication the first time! MRS */
//*************
//  cglo_flags = (CGLO_TEMPERATURE | CGLO_GSTAT
//                  | (bVV ? CGLO_PRESSURE:0)
//                  | (bVV ? CGLO_CONSTRAINT:0)
//                  | (bRerunMD ? CGLO_RERUNMD:0)
//                  | ((Flags & MD_READ_EKIN) ? CGLO_READEKIN:0));
//    
//    bSumEkinhOld = FALSE;
//    compute_globals(fplog,gstat,cr,ir,fr,ekind,state,state_global,mdatoms,nrnb,vcm,
//                    wcycle,enerd,force_vir,shake_vir,total_vir,pres,mu_tot,
//                    constr,NULL,NULL,NULL,NULL,NULL,
//                    NULL,state->box,
//                    top_global,&pcurr,top_global->natoms,&bSumEkinhOld,cglo_flags);
//    if (ir->eI == eiVVAK) {
//        /* a second call to get the half step temperature initialized as well */ 
//        /* we do the same call as above, but turn the pressure off -- internally, this 
//           is recognized as a velocity verlet half-step kinetic energy calculation.
//           This minimized excess variables, but perhaps loses some logic?*/
//        
//        compute_globals(fplog,gstat,cr,ir,fr,ekind,state,state_global,mdatoms,nrnb,vcm,
//                        wcycle,enerd,force_vir,shake_vir,total_vir,pres,mu_tot,
//                        constr,NULL,NULL,NULL,NULL,NULL,NULL,state->box,
//                        top_global,&pcurr,top_global->natoms,&bSumEkinhOld,
//                        cglo_flags &~ CGLO_PRESSURE);
//    }
    
    /* Calculate the initial half step temperature, and save the ekinh_old */
    if (!(Flags & MD_STARTFROMCPT)) 
    {
        for(i=0; (i<ir->opts.ngtc); i++) 
        {
            copy_mat(ekind->tcstat[i].ekinh,ekind->tcstat[i].ekinh_old);
        } 
    }
    temp0 = enerd->term[F_TEMP];
    
    /* if using an iterative algorithm, we need to create a working directory for the state. */
    if (bIterations) 
    {
            bufstate = init_bufstate(state);
    }
    if (bFFscan) 
    {
        snew(xcopy,state->natoms);
        snew(vcopy,state->natoms);
        copy_rvecn(state->x,xcopy,0,state->natoms);
        copy_rvecn(state->v,vcopy,0,state->natoms);
        copy_mat(state->box,boxcopy);
    } 
    
    /* need to make an initiation call to get the Trotter variables set, as well as other constants for non-trotter
       temperature control */
    trotter_seq = init_npt_vars(ir,state,&MassQ,bTrotter);
    
    if (MASTER(cr))
    {
        if (constr && !ir->bContinuation && ir->eConstrAlg == econtLINCS)
        {
            fprintf(fplog,
                    "RMS relative constraint deviation after constraining: %.2e\n",
                    constr_rmsd(constr,FALSE));
        }
        if (!EI_VV(ir->eI)) 
        {
            enerd->term[F_TEMP] *= 2; /* result of averages being done over previous and current step,
                                         and there is no previous step */
        }
        fprintf(fplog,"Initial temperature: %g K\n",enerd->term[F_TEMP]);
        if (bRerunMD)
        {
            fprintf(stderr,"starting md rerun '%s', reading coordinates from"
                    " input trajectory '%s'\n\n",
                    *(top_global->name),opt2fn("-rerun",nfile,fnm));
            if (bVerbose)
            {
                fprintf(stderr,"Calculated time to finish depends on nsteps from "
                        "run input file,\nwhich may not correspond to the time "
                        "needed to process input trajectory.\n\n");
            }
        }
        else
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
    }

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

    /* if rerunMD then read coordinates and velocities from input trajectory */
    if (bRerunMD)
    {
        if (getenv("GMX_FORCE_UPDATE"))
        {
            bForceUpdate = TRUE;
        }

        bNotLastFrame = read_first_frame(oenv,&status,
                                         opt2fn("-rerun",nfile,fnm),
                                         &rerun_fr,TRX_NEED_X | TRX_READ_V);
        if (rerun_fr.natoms != top_global->natoms)
        {
            gmx_fatal(FARGS,
                      "Number of atoms in trajectory (%d) does not match the "
                      "run input file (%d)\n",
                      rerun_fr.natoms,top_global->natoms);
        }
        if (ir->ePBC != epbcNONE)
        {
            if (!rerun_fr.bBox)
            {
                gmx_fatal(FARGS,"Rerun trajectory frame step %d time %f does not contain a box, while pbc is used",rerun_fr.step,rerun_fr.time);
            }
            if (max_cutoff2(ir->ePBC,rerun_fr.box) < sqr(fr->rlistlong))
            {
                gmx_fatal(FARGS,"Rerun trajectory frame step %d time %f has too small box dimensions",rerun_fr.step,rerun_fr.time);
            }

            /* Set the shift vectors.
             * Necessary here when have a static box different from the tpr box.
             */
            calc_shifts(rerun_fr.box,fr->shift_vec);
        }
    }

    /* loop over MD steps or if rerunMD to end of input trajectory */
    bFirstStep = TRUE;
    /* Skip the first Nose-Hoover integration when we get the state from tpx */
    bStateFromTPX = !opt2bSet("-cpi",nfile,fnm);
    bInitStep = bFirstStep && (bStateFromTPX || bVV);
    bStartingFromCpt = (Flags & MD_STARTFROMCPT) && bInitStep;
    bLastStep = FALSE;
    bSumEkinhOld = FALSE;
    bExchanged = FALSE;

    step = ir->init_step;
    step_rel = 0;

//    if (ir->nstlist == -1)
//    {
//        init_nlistheuristics(&nlh,bGStatEveryStep,step);
//    }

    bLastStep = (bRerunMD || (ir->nsteps >= 0 && step_rel > ir->nsteps));
    while (!bLastStep || (bRerunMD && bNotLastFrame)) {

        wallcycle_start(wcycle,ewcSTEP);

        GMX_MPE_LOG(ev_timestep1);

        /* Now we have the energies and forces corresponding to the 
         * coordinates at time t. We must output all of this before
         * the update.
         * for RerunMD t is read from input trajectory
         */
        GMX_MPE_LOG(ev_output_start);

        bX   = do_per_step(step,ir->nstxout) || (bLastStep && ir->nstxout);
        bV   = do_per_step(step,ir->nstvout) || (bLastStep && ir->nstvout);
        bF   = do_per_step(step,ir->nstfout) || (bLastStep && ir->nstfout);
        bXTC = do_per_step(step,ir->nstxtcout) || (bLastStep && ir->nstxtcout);
        do_ene = do_per_step(step,ir->nstenergy) || (bLastStep && ir->nstenergy);

//#ifdef GMX_OPENMM

      if( bX || bXTC || bV ){                                                                                                                                                                 
        wallcycle_start(wcycle,ewcTRAJ);                                                                                                                                                              
        openmm_copy_state(openmmData, state, &t, f, enerd, bX||bXTC, bV, 0, 0);                                                                                                               
        wallcycle_stop(wcycle,ewcTRAJ);                                                                                                                                                       
      }                                                                                                                                                                                       

      openmm_take_one_step(openmmData);                                                                                                                                                       
      bLastStep = (step_rel == ir->nsteps);
      if (bX || bV || bF || bXTC || do_ene) {
        wallcycle_start(wcycle,ewcTRAJ);
        if( bF || do_ene ){                                                                                                                                                                   
           openmm_copy_state(openmmData, state, &t, f, enerd, 0, 0, bF, do_ene);                                                                                                              
      }                                                                                                                                                                                     
      upd_mdebin(mdebin,fp_dhdl,bGStatEveryStep && !bRerunMD,
		 t,mdatoms->tmass,enerd,state,lastbox,
		 shake_vir,force_vir,total_vir,pres,
		 ekind,mu_tot,constr);
      print_ebin(fp_ene,do_ene,FALSE,FALSE,do_log?fplog:NULL,step,t,
		 eprNORMAL,bCompact,mdebin,fcd,groups,&(ir->opts)); // XXX check this do_log stuff!
//#else

            write_traj(fplog,cr,fp_trn,bX,bV,bF,fp_xtc,bXTC,ir->xtcprec,
                       fn_cpt,bCPT,top_global,ir->eI,ir->simulation_part,
                       step,t,state,state_global,f,f_global,&n_xtc,&x_xtc);
            debug_gmx();
            if (bLastStep && step_rel == ir->nsteps &&
                (Flags & MD_CONFOUT) && MASTER(cr) &&
                !bRerunMD && !bFFscan)
            {
                /* x and v have been collected in write_traj */
                fprintf(stderr,"\nWriting final coordinates.\n");
                if (ir->ePBC != epbcNONE && !ir->bPeriodicMols &&
                    DOMAINDECOMP(cr))
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
        
        bFirstStep = FALSE;
        bInitStep = FALSE;
        step++;
        step_rel++;
    }
    /* End of main MD loop */
    debug_gmx();
    
//#ifdef GMX_OPENMM
    openmm_cleanup(fplog, openmmData);
//#endif

    /* Stop the time */
    runtime_end(runtime);
    
    if (bRerunMD)
    {
        close_trj(status);
    }
    
    if (!(cr->duty & DUTY_PME))
    {
        /* Tell the PME only node to finish */
        gmx_pme_finish(cr);
    }
    
    if (MASTER(cr))
    {
        if (bGStatEveryStep && !bRerunMD) 
        {
            print_ebin(fp_ene,FALSE,FALSE,FALSE,fplog,step,t,
                       eprAVER,FALSE,mdebin,fcd,groups,&(ir->opts));
        }
        close_enx(fp_ene);
        if (ir->nstxtcout)
        {
            close_xtc(fp_xtc);
        }
        close_trn(fp_trn);
        if (fp_dhdl)
        {
            gmx_fio_fclose(fp_dhdl);
        }
        if (fp_field)
        {
            gmx_fio_fclose(fp_field);
        }
    }
    debug_gmx();

//    if (ir->nstlist == -1 && nlh.nns > 0 && fplog)
//    {
//        fprintf(fplog,"Average neighborlist lifetime: %.1f steps, std.dev.: %.1f steps\n",nlh.s1/nlh.nns,sqrt(nlh.s2/nlh.nns - sqr(nlh.s1/nlh.nns)));
//        fprintf(fplog,"Average number of atoms that crossed the half buffer length: %.1f\n\n",nlh.ab/nlh.nns);
//    }
    
    if (shellfc && fplog)
    {
        fprintf(fplog,"Fraction of iterations that converged:           %.2f %%\n",
                (nconverged*100.0)/step_rel);
        fprintf(fplog,"Average number of force evaluations per MD step: %.2f\n\n",
                tcount/step_rel);
    }
    
    if (repl_ex_nst > 0 && MASTER(cr))
    {
        print_replica_exchange_statistics(fplog,repl_ex);
    }
    
    runtime->nsteps_done = step_rel;
    
    return 0;
}

