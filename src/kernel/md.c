/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
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
#include "confio.h"
#include "network.h"
#include "pull.h"
#include "xvgr.h"
#include "physics.h"
#include "names.h"
#include "xmdrun.h"
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

#ifdef GMX_MPI
#include <mpi.h>
#endif

/* The following two variables and the signal_handler function
 * is used from pme.c as well 
 */
extern bool bGotTermSignal, bGotUsr1Signal;

static RETSIGTYPE signal_handler(int n)
{
  switch (n) {
  case SIGTERM:
    bGotTermSignal = TRUE;
    break;
  case SIGUSR1:
    bGotUsr1Signal = TRUE;
    break;
  }
}

typedef struct { 
  gmx_integrator_t *func;
} gmx_intp_t;

/* The array should match the eI array in include/types/enums.h */
const gmx_intp_t integrator[eiNR] = { {do_md}, {do_steep}, {do_cg}, {do_md}, {do_md}, {do_nm}, {do_lbfgs}, {do_tpi}, {do_tpi}, {do_md} };

void mdrunner(FILE *fplog,t_commrec *cr,int nfile,t_filenm fnm[],
	      bool bVerbose,bool bCompact,
	      ivec ddxyz,int dd_node_order,real rdd,real rconstr,
	      real dlb_scale,char *ddcsx,char *ddcsy,char *ddcsz,
	      int nstepout,gmx_edsam_t ed,int repl_ex_nst,int repl_ex_seed,
	      real pforce,real cpt_period,real max_hours,
	      unsigned long Flags)
{
  double     nodetime=0,realtime,t;
  t_inputrec *inputrec;
  t_state    *state=NULL;
  matrix     box;
  rvec       *buf=NULL,*f=NULL;
  real       tmpr1,tmpr2;
  real       *ener=NULL;
  t_nrnb     *nrnb;
  gmx_mtop_t *mtop=NULL;
  t_groups   *grps=NULL;
  t_mdatoms  *mdatoms=NULL;
  t_forcerec *fr=NULL;
  t_fcdata   *fcd=NULL;
  real       ewaldcoeff=0;
  gmx_pme_t  *pmedata=NULL;
  time_t     start_t=0;
  gmx_vsite_t *vsite=NULL;
  gmx_constr_t constr;
  int        i,m,nChargePerturbed=-1,status,nalloc;
  char       *gro;
  gmx_wallcycle_t wcycle;
  bool       bReadRNG;

  snew(inputrec,1);
  snew(mtop,1);
  
  if (bVerbose && SIMMASTER(cr)) 
    fprintf(stderr,"Getting Loaded...\n");
  
  if (PAR(cr)) {
    /* The master thread on the master node reads from disk, 
     * then distributes everything to the other processors.
     */
    snew(state,1);
    init_parallel(fplog,ftp2fn(efTPX,nfile,fnm),cr,
		  inputrec,mtop,state,
		  SIMMASTER(cr) ? LIST_SCALARS | LIST_INPUTREC : 0);
  } else {
    /* Read a file for a single processor */
    snew(state,1);
    init_single(fplog,inputrec,ftp2fn(efTPX,nfile,fnm),mtop,state);
  }
  if (!EEL_PME(inputrec->coulombtype) || (Flags & MD_PARTDEC)) {
    cr->npmenodes = 0;
  }
  
  /* NMR restraints must be initialized before read_checkpoint,
   * since with time averaging the history is added to t_state.
   * For proper consistency check we therefore need to extend
   * t_state here.
   * So the PME-only nodes (if present) will also initialize
   * the distance restraints.
   */
  snew(fcd,1);

  /* This needs to be called before read_checkpoint to extend the state */
  init_disres(fplog,mtop,inputrec,cr,Flags & MD_PARTDEC,fcd,state);

  if (gmx_mtop_ftype_count(mtop,F_ORIRES) > 0) {
    if (PAR(cr) && !(Flags & MD_PARTDEC)) {
      gmx_fatal(FARGS,"Orientation restraints do not work (yet) with domain decomposition, use particle decomposition (mdrun option -pd)");
    }
    /* Orientation restraints */
    if (MASTER(cr)) {
      init_orires(fplog,mtop,state->x,inputrec,cr->ms,&(fcd->orires),state);
    }
  }
  
  if (DEFORM(*inputrec)) {
    /* Store the deform reference box before reading the checkpoint */
    set_deform_reference_box(inputrec->init_step,state->box);
  }

  if (opt2bSet("-cpi",nfile,fnm)) {
    if (SIMMASTER(cr)) {
      /* Read the state from the checkpoint file */
      read_checkpoint(opt2fn("-cpi",nfile,fnm),fplog,
		      cr,Flags & MD_PARTDEC,ddxyz,
		      inputrec->eI,&i,&t,state,&bReadRNG);
    }
    if (PAR(cr)) {
      gmx_bcast(sizeof(i),&i,cr);
      gmx_bcast(sizeof(bReadRNG),&bReadRNG,cr);
      gmx_bcast(DIM*sizeof(ddxyz[0]),ddxyz,cr);
    }
    inputrec->bContinuation = TRUE;
    inputrec->init_step += i;
    inputrec->nsteps    -= i;
    if (bReadRNG) {
      Flags |= MD_READ_RNG;
    }
  }

  if (SIMMASTER(cr)) {
    copy_mat(state->box,box);
  }
  if (PAR(cr)) {
    gmx_bcast(sizeof(box),box,cr);
  }
    
  if (bVerbose && SIMMASTER(cr))
    fprintf(stderr,"Loaded with Money\n\n");

  if (PAR(cr) && !((Flags & MD_PARTDEC) || EI_TPI(inputrec->eI))) {
    cr->dd = init_domain_decomposition(fplog,cr,ddxyz,rdd,rconstr,
				       Flags & MD_DLB,dlb_scale,
				       ddcsx,ddcsy,ddcsz,
				       mtop,box,inputrec);
    
    make_dd_communicators(fplog,cr,dd_node_order);

    /* Set overallocation to avoid frequent reallocation of arrays */
    set_over_alloc_dd(TRUE);
  } else {
    cr->duty = (DUTY_PP | DUTY_PME);

    if (inputrec->ePBC == epbcSCREW)
      gmx_fatal(FARGS,"pbc=%s is only implemented with domain decomposition",
		epbc_names[inputrec->ePBC]);
  }

  if (PAR(cr)) {
    /* After possible communicator splitting in make_dd_communicators.
     * we can set up the intra/inter node communication.
     */
    gmx_setup_nodecomm(fplog,cr);
  }

  wcycle = wallcycle_init(fplog,cr);

  snew(nrnb,1);
  if (cr->duty & DUTY_PP) {
    /* Initiate everything (snew sets to zero!) */
    snew(ener,F_NRE);
    snew(grps,1);

    /* For domain decomposition we allocate dynamically
     * in dd_partition_system.
     */
    if (DOMAINDECOMP(cr)) {
      bcast_state_setup(cr,state);
    } else {
      if (PAR(cr)) {
	if (!MASTER(cr)) {
	  snew(state,1);
	}
	bcast_state(cr,state,TRUE);
      }

      snew(buf,mtop->natoms);
      snew(f,mtop->natoms);
    }
    
    /* Group stuff (energies etc) */
    init_t_groups(fplog,mtop,&(inputrec->opts),grps);
    /* Copy the cos acceleration to the groups struct */
    grps->cosacc.cos_accel = inputrec->cos_accel;
    
    /* Dihedral Restraints */
    if (gmx_mtop_ftype_count(mtop,F_DIHRES) > 0) {
      init_dihres(fplog,mtop,inputrec,fcd);
    }
    
    /* Initiate forcerecord */
    fr = mk_forcerec();
    init_forcerec(fplog,fr,fcd,inputrec,mtop,cr,box,FALSE,
		  opt2fn("-table",nfile,fnm),opt2fn("-tablep",nfile,fnm),
		  opt2fn("-tableb",nfile,fnm),FALSE,pforce);
    fr->bSepDVDL = ((Flags & MD_SEPPOT) == MD_SEPPOT);
    
    /* Initialize QM-MM */
    if(fr->bQMMM)
      init_QMMMrec(cr,box,mtop,inputrec,fr);
    
    /* Initialize the mdatoms structure.
     * mdatoms is not filled with atom data,
     * as this can not be done now with domain decomposition.
     */
    mdatoms = init_mdatoms(fplog,mtop,inputrec->efep!=efepNO);
    
    /* Initialize the virtual site communication */
    vsite = init_vsite(mtop,cr);

    calc_shifts(box,fr->shift_vec);

    /* With periodic molecules the charge groups should be whole at start up
     * and the virtual sites should not be far from their proper positions.
     */
    if (!inputrec->bContinuation && MASTER(cr) &&
	!(inputrec->ePBC != epbcNONE && inputrec->bPeriodicMols)) {
      /* Make molecules whole at start of run */
      if (fr->ePBC != epbcNONE)  {
	do_pbc_first_mtop(fplog,inputrec->ePBC,box,fr,mtop,state->x);
      }
      if (vsite) {
	/* Correct initial vsite positions are required
	 * for the initial distribution in the domain decomposition
	 * and for the initial shell prediction.
	 */
	construct_vsites_mtop(fplog,vsite,mtop,state->x);
      }
    }

    /* Initiate PPPM if necessary */
    if (fr->eeltype == eelPPPM) {
      if (mdatoms->nChargePerturbed)
	gmx_fatal(FARGS,"Free energy with %s is not implemented",
		  eel_names[fr->eeltype]);
      status = gmx_pppm_init(fplog,cr,FALSE,TRUE,box,
			     getenv("GMXGHAT"),inputrec, (Flags & MD_REPRODUCIBLE));
      if (status != 0)
	gmx_fatal(FARGS,"Error %d initializing PPPM",status);
    }

    if (EEL_PME(fr->eeltype)) {
      ewaldcoeff = fr->ewaldcoeff;
      pmedata = &fr->pmedata;
    } else {
      pmedata = NULL;
    }
  } else {
    /* This is a PME only node */
    
    /* We don't need the state */
    done_state(state);

    ewaldcoeff = calc_ewaldcoeff(inputrec->rcoulomb, inputrec->ewald_rtol);
    snew(pmedata,1);
  }

  /* Initiate PME if necessary,
   * either on all nodes or on dedicated PME nodes only. */
  if (EEL_PME(inputrec->coulombtype)) {
    if (mdatoms) {
      nChargePerturbed = mdatoms->nChargePerturbed;
    }
    if (cr->npmenodes > 0) {
      /* The PME only nodes need to know nChargePerturbed */
      gmx_bcast_sim(sizeof(nChargePerturbed),&nChargePerturbed,cr);
    }
    if (cr->duty & DUTY_PME) {
      status = gmx_pme_init(pmedata,cr,inputrec,
			    mtop ? mtop->natoms : 0,nChargePerturbed,
			    (Flags & MD_REPRODUCIBLE));
      if (status != 0)
	gmx_fatal(FARGS,"Error %d initializing PME",status);
    }
  }
  
  if (integrator[inputrec->eI].func == do_md) {
    /* Turn on signal handling on all nodes */
    /*
     * (A user signal from the PME nodes (if any)
     * is communicated to the PP nodes.
     */
    if (getenv("GMX_NO_TERM") == NULL) {
      if (debug)
	fprintf(debug,"Installing signal handler for SIGTERM\n");
      signal(SIGTERM,signal_handler);
    }
    if (getenv("GMX_NO_USR1") == NULL) {
      if (debug)
	fprintf(debug,"Installing signal handler for SIGUSR1\n");
      signal(SIGUSR1,signal_handler);
    }
  }

  if (cr->duty & DUTY_PP) {
    if (inputrec->ePull != epullNO) {
      /* Initialize pull code */
      init_pull(fplog,inputrec,nfile,fnm,mtop,cr,
		EI_DYNAMICS(inputrec->eI) && MASTER(cr));
    }

    constr = init_constraints(fplog,mtop,inputrec,ed,state,cr);

    if (DOMAINDECOMP(cr)) {
      dd_make_reverse_top(fplog,cr->dd,mtop,vsite,constr,inputrec,
			  Flags & MD_DDBONDCHECK);

      set_dd_parameters(fplog,cr->dd,dlb_scale,inputrec,fr,box);
     
      setup_dd_grid(fplog,cr->dd);
    }

    /* Now do whatever the user wants us to do (how flexible...) */
    start_t = integrator[inputrec->eI].func(fplog,cr,nfile,fnm,
					    bVerbose,bCompact,
					    vsite,constr,
					    nstepout,inputrec,grps,mtop,
					    ener,fcd,state,f,buf,
					    mdatoms,nrnb,wcycle,ed,fr,
					    repl_ex_nst,repl_ex_seed,
					    cpt_period,max_hours,
					    Flags);

    if (inputrec->ePull != epullNO)
      finish_pull(fplog,inputrec->pull);
  } else {
    /* do PME only */
    gmx_pmeonly(*pmedata,cr,nrnb,wcycle,ewaldcoeff,FALSE);
  }
 
  /* Some timing stats */  
  if (MASTER(cr)) {
    realtime=difftime(time(NULL),start_t);
    if ((nodetime=node_time()) == 0)
      nodetime=realtime;
  }
  else 
    realtime=0;

  wallcycle_stop(wcycle,ewcRUN);
    
  /* Finish up, write some stuff
   * if rerunMD, don't write last frame again 
   */
  finish_run(fplog,cr,ftp2fn(efSTO,nfile,fnm),
	     inputrec,nrnb,wcycle,nodetime,realtime,inputrec->nsteps,
	     EI_DYNAMICS(inputrec->eI) && !MULTISIM(cr));
  
  /* Does what it says */  
  print_date_and_time(fplog,cr->nodeid,"Finished mdrun");
}

time_t do_md(FILE *fplog,t_commrec *cr,int nfile,t_filenm fnm[],
	     bool bVerbose,bool bCompact,
	     gmx_vsite_t *vsite,gmx_constr_t constr,
	     int stepout,t_inputrec *ir,t_groups *grps,
	     gmx_mtop_t *top_global,
	     real ener[],t_fcdata *fcd,
	     t_state *state_global,rvec f[],
	     rvec buf[],t_mdatoms *mdatoms,
	     t_nrnb *nrnb,gmx_wallcycle_t wcycle,
	     gmx_edsam_t ed,t_forcerec *fr,
	     int repl_ex_nst,int repl_ex_seed,
	     real cpt_period,real max_hours,
	     unsigned long Flags)
{
  t_mdebin   *mdebin;
  int        fp_ene=0,fp_trn=0,fp_xtc=0,step,step_rel;
  char       *fn_cpt;
  FILE       *fp_dgdl=NULL,*fp_field=NULL;
  time_t     start_t;
  real       t,t0,lam0;
  bool       bGStat,bNS,bSimAnn,bStopCM,bRerunMD,bNotLastFrame=FALSE,
             bFirstStep,bStateFromTPX,bLastStep;
  bool       bNEMD,do_log,do_verbose,bRerunWarnNoV=TRUE,
	     bForceUpdate=FALSE,bX,bV,bF,bXTC,bCPT;
  bool       bMasterState;
  tensor     force_vir,shake_vir,total_vir,pres,ekin;
  int        i,m,status;
  rvec       mu_tot;
  t_vcm      *vcm;
  int        nabnsb,ns_step=0,nns=0;
  double     ns_s1=0,ns_s2=0,ns_ab=0;
  matrix     *scale_tot;
  t_trxframe rerun_fr;
  gmx_repl_ex_t repl_ex=NULL;
  int        nchkpt=1;
  /* Booleans (disguised as a reals) to checkpoint and terminate mdrun */  
  real       chkpt=0,terminate=0,terminate_now=0;

  t_topology *top;
  t_state    *state=NULL;
  rvec       *f_global=NULL;
  gmx_stochd_t sd=NULL;
  t_graph    *graph=NULL;

  bool        bFFscan;
  gmx_groups_t *groups;
  gmx_shellfc_t shellfc;
  int         count,nconverged=0;
  real        timestep=0;
  double      tcount=0;
  bool        bHaveConstr=FALSE,bIonize=FALSE,bGlas=FALSE;
  bool        bTCR=FALSE,bConverged=TRUE,bOK,bSumEkinhOld,bExchanged;
  real        temp0,mu_aver=0,dvdl;
  int         a0,a1,gnx,ii;
  atom_id     *grpindex;
  char        *grpname;
  t_coupl_rec *tcr=NULL;
  rvec        *xcopy=NULL,*vcopy=NULL;
  matrix      boxcopy,lastbox;
  double      cycles;

  /* Check for special mdrun options */
  bRerunMD = (Flags & MD_RERUN);
  bIonize  = (Flags & MD_IONIZE);
  bGlas    = (Flags & MD_GLAS);
  bFFscan  = (Flags & MD_FFSCAN);
  bGStat   = !(Flags & MD_NOGSTAT);

  if (!bGStat) {
    if (EI_DYNAMICS(ir->eI)) {
      if (fplog) {
	fprintf(fplog,"\nWill not sum the energies at every step,\n"
		"therefore the energy file does not contain exact averages and fluctuations.\n\n");
      }
      if (ir->etc != etcNO || ir->epc != epcNO) {
	if (fplog) {
	  fprintf(fplog,"WARNING:\nThe temperature and/or pressure for scaling will only be updated every nstlist (%d) steps\n\n",ir->nstlist);
	}
      }
      if (ir->comm_mode != ecmNO && ir->nstcomm != ir->nstlist) {
	if (fplog) {
	  fprintf(fplog,"WARNING:\nBecause of the no energy summing option setting nstcomm (was %d) to nstlist (%d)\n\n",ir->nstcomm,ir->nstlist);
	}
      }
    } else {
      char *warn="\nWARNING:\nNo energy summing can only be used with dynamics, ignoring this option\n";
      fprintf(stderr,"%s\n",warn);
      if (fplog)
	fprintf(fplog,"%s\n",warn);
      bGStat = TRUE;
    }
  }

  if (bRerunMD || bFFscan)
    ir->nstxtcout = 0;

  groups = &top_global->groups;

  /* Initial values */
  init_md(fplog,cr,ir,&t,&t0,&state_global->lambda,&lam0,
	  nrnb,top_global,&sd,
	  nfile,fnm,&fp_trn,&fp_xtc,&fp_ene,&fn_cpt,
	  &fp_dgdl,&fp_field,&mdebin,
	  force_vir,shake_vir,mu_tot,&bNEMD,&bSimAnn,&vcm);

  debug_gmx();

  /* Check for polarizable models and flexible constraints */
  shellfc = init_shell_flexcon(fplog,
			       top_global,n_flexible_constraints(constr),
			       (ir->bContinuation || 
				(DOMAINDECOMP(cr) && !MASTER(cr))) ?
			       NULL : state_global->x);

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
      graph = mk_graph(fplog,&(top->idef),0,top->atoms.nr,FALSE,FALSE);
    }

    if (shellfc) {
      make_local_shells(cr,mdatoms,shellfc);
    }
  }

  if (DOMAINDECOMP(cr)) {
    /* Distribute the charge groups over the nodes from the master node */
    dd_partition_system(fplog,ir->init_step,cr,TRUE,
			state_global,top_global,ir,
			state,&f,&buf,mdatoms,top,fr,vsite,shellfc,constr,
			nrnb,wcycle,FALSE);
  }

  update_mdatoms(mdatoms,state->lambda);

  if (sd && (Flags & MD_READ_RNG)) {
    /* Set the random state if we read a checkpoint file */
    set_stochd_state(sd,state);
  }

  /* Initialize constraints */
  if (constr) {
    if (!DOMAINDECOMP(cr))
      set_constraints(constr,top,ir,mdatoms,NULL);
    bHaveConstr = TRUE;
  }

  gnx = top->mols.nr;
  snew(grpindex,gnx);
  for(i=0; (i<gnx); i++)
    grpindex[i] = i;

  /* Check whether we have to GCT stuff */
  bTCR = ftp2bSet(efGCT,nfile,fnm);
  if (MASTER(cr) && bTCR)
    fprintf(stderr,"Will do General Coupling Theory!\n");

  if (repl_ex_nst > 0 && MASTER(cr))
    repl_ex = init_replica_exchange(fplog,cr->ms,state_global,ir,
				    repl_ex_nst,repl_ex_seed);

  if (bHaveConstr && !ir->bContinuation && !bRerunMD) {
    do_shakefirst(fplog,constr,ir,mdatoms,state,buf,f,
		  graph,cr,nrnb,grps,fr,top);

    if (vsite) {
      construct_vsites(fplog,vsite,state->x,nrnb,ir->delta_t,NULL,
		       top->idef.iparams,top->idef.il,
		       fr->ePBC,fr->bMolPBC,graph,cr,state->box);
    }
  }

  debug_gmx();

  /* Compute initial EKin for all.. */
  if (grps->cosacc.cos_accel == 0)
    calc_ke_part(state->v,&(ir->opts),mdatoms,grps,nrnb,state->lambda);
  else
    calc_ke_part_visc(state->box,state->x,state->v,&(ir->opts),
		      mdatoms,grps,nrnb,state->lambda);
  debug_gmx();

  if (PAR(cr))
  {
    GMX_MPE_LOG(ev_global_stat_start);
       
    global_stat(fplog,cr,ener,force_vir,shake_vir,mu_tot,
		ir,grps,FALSE,constr,vcm,NULL,NULL,&terminate);

    GMX_MPE_LOG(ev_global_stat_finish);
  }
  debug_gmx();
  
  /* Calculate the initial half step temperature */
  temp0 = sum_ekin(TRUE,&(ir->opts),grps,ekin,NULL);

  debug_gmx();
   
  /* Initiate data for the special cases */
  if (bFFscan) {
    snew(xcopy,state->natoms);
    snew(vcopy,state->natoms);
    for(ii=0; (ii<state->natoms); ii++) {
      copy_rvec(state->x[ii],xcopy[ii]);
      copy_rvec(state->v[ii],vcopy[ii]);
    }
    copy_mat(state->box,boxcopy);
  } 

  if (MASTER(cr)) {
    if (constr && !ir->bContinuation && ir->eConstrAlg == econtLINCS)
      fprintf(fplog,
	      "RMS relative constraint deviation after constraining: %.2e\n",
	      constr_rmsd(constr,FALSE));
    fprintf(fplog,"Initial temperature: %g K\n",temp0);
    if (bRerunMD) {
      fprintf(stderr,"starting md rerun '%s', reading coordinates from"
	      " input trajectory '%s'\n\n",
	      *(top->name),opt2fn("-rerun",nfile,fnm));
      if (bVerbose)
	fprintf(stderr,"Calculated time to finish depends on nsteps from "
		"run input file,\nwhich may not correspond to the time "
		"needed to process input trajectory.\n\n");
    } else {
      fprintf(stderr,"starting mdrun '%s'\n%d steps, %8.1f ps.\n",
	      *(top->name),ir->nsteps,ir->nsteps*ir->delta_t);
    }
    fprintf(fplog,"\n");
  }

  if (ir->nstlist == -1)
    snew(scale_tot,1);
  else
    scale_tot = NULL;

  /* Write start time */
  start_t=print_date_and_time(fplog,cr->nodeid,"Started mdrun");
  wallcycle_start(wcycle,ewcRUN);
  if (fplog)
    fprintf(fplog,"\n");

  /* Set the node time counter to 0 after initialisation */
  start_time();
  debug_gmx();
  /***********************************************************
   *
   *             Loop over MD steps 
   *
   ************************************************************/
  
  /* if rerunMD then read coordinates and velocities from input trajectory */
  if (bRerunMD) {
    if (getenv("GMX_FORCE_UPDATE"))
      bForceUpdate = TRUE;

    bNotLastFrame = read_first_frame(&status,opt2fn("-rerun",nfile,fnm),
				     &rerun_fr,TRX_NEED_X | TRX_READ_V);
    if (rerun_fr.natoms != mdatoms->nr)
      gmx_fatal(FARGS,"Number of atoms in trajectory (%d) does not match the "
		  "run input file (%d)\n",rerun_fr.natoms,mdatoms->nr);
    if (ir->ePBC != epbcNONE) {
      if (!rerun_fr.bBox)
	gmx_fatal(FARGS,"Rerun trajectory frame step %d time %f does not contain a box, while pbc is used",rerun_fr.step,rerun_fr.time);
      if (max_cutoff2(ir->ePBC,rerun_fr.box) < sqr(fr->rlistlong))
	gmx_fatal(FARGS,"Rerun trajectory frame step %d time %f has too small box dimensions",rerun_fr.step,rerun_fr.time);
    }
  }

  /* loop over MD steps or if rerunMD to end of input trajectory */
  bFirstStep = TRUE;
  /* Skip the first Nose-Hoover integration when we get the state from tpx */
  bStateFromTPX = !opt2bSet("-cpi",nfile,fnm);
  bLastStep = FALSE;
  bSumEkinhOld = FALSE,
  bExchanged = FALSE;
  nabnsb = 1;
  step = ir->init_step;
  step_rel = 0;

  bLastStep = (bRerunMD || step_rel > ir->nsteps);
  while (!bLastStep || (bRerunMD && bNotLastFrame)) {

    wallcycle_start(wcycle,ewcSTEP);

    GMX_MPE_LOG(ev_timestep1);

    if (bRerunMD) {
      if (rerun_fr.bStep) {
	step = rerun_fr.step;
	step_rel = step - ir->init_step;
      }
      if (rerun_fr.bTime)
	t = rerun_fr.time;
      else
	t = step;
    } else {
      bLastStep = (step_rel == ir->nsteps);

      t = t0 + step*ir->delta_t;
    }
    
    if (ir->efep != efepNO) {
      if (bRerunMD && rerun_fr.bLambda && (ir->delta_lambda!=0))
	state->lambda = rerun_fr.lambda;
      else
	state->lambda = lam0 + step*ir->delta_lambda;
    }
    
    if (bSimAnn) 
      update_annealing_target_temp(&(ir->opts),t);
    
    if (bRerunMD) {
      for(i=0; i<mdatoms->nr; i++)
	copy_rvec(rerun_fr.x[i],state->x[i]);
      if (rerun_fr.bV)
	for(i=0; i<mdatoms->nr; i++)
	  copy_rvec(rerun_fr.v[i],state->v[i]);
      else {
	for(i=0; i<mdatoms->nr; i++)
	    clear_rvec(state->v[i]);
	if (bRerunWarnNoV) {
	  fprintf(stderr,"\nWARNING: Some frames do not contain velocities.\n"
		  "         Ekin, temperature and pressure are incorrect,\n"
		  "         the virial will be incorrect when constraints are present.\n"
		  "\n");
	  bRerunWarnNoV = FALSE;
	}
      }
      copy_mat(rerun_fr.box,state->box);

      if (vsite) {
	if (graph) {
	  /* Following is necessary because the graph may get out of sync
	   * with the coordinates if we only have every N'th coordinate set
	   */
	  mk_mshift(fplog,graph,fr->ePBC,state->box,state->x);
	  shift_self(graph,state->box,state->x);
	}
	construct_vsites(fplog,vsite,state->x,nrnb,ir->delta_t,state->v,
			 top->idef.iparams,top->idef.il,
			 fr->ePBC,fr->bMolPBC,graph,cr,state->box);
	if (graph)
	  unshift_self(graph,state->box,state->x);
      }
    }

    /* Stop Center of Mass motion */
    bStopCM = do_per_step(step,abs(ir->nstcomm));

    /* Copy back starting coordinates in case we're doing a forcefield scan */
    if (bFFscan) {
      for(ii=0; (ii<state->natoms); ii++) {
	copy_rvec(xcopy[ii],state->x[ii]);
	copy_rvec(vcopy[ii],state->v[ii]);
      }
      copy_mat(boxcopy,state->box);
    }
    
    bNS = bFirstStep;
    if (bRerunMD) {
      /* for rerun MD always do Neighbour Searching */
      if (ir->nstlist != 0)
	bNS = TRUE;
    } else {
      /* Determine whether or not to do Neighbour Searching */
      if (bExchanged || (ir->nstlist>0 && (step % ir->nstlist==0))) {
	bNS = TRUE;
      } else if (ir->nstlist == -1) {
	bNS = (nabnsb > 0);
	if (bNS) {
	  if (debug)
	    fprintf(debug,"%d atoms beyond ns buffer, updating neighbor list after %d steps\n",nabnsb,step-ns_step);
	  nns++;
	  ns_s1 += step - ns_step;
	  ns_s2 += sqr(step - ns_step);
	  ns_ab += nabnsb;
	  ns_step = step;
	  /* Initialize the cumulative coordinate scaling matrix */
	  clear_mat(*scale_tot);
	  for(ii=0; ii<DIM; ii++)
	    (*scale_tot)[ii][ii] = 1.0;
	}
      }
    } 

    if (terminate_now > 0 || (terminate_now < 0 && bNS))
      bLastStep = TRUE;

    do_log = do_per_step(step,ir->nstlog) || bFirstStep || bLastStep;
    do_verbose = bVerbose && (step % stepout == 0 || bFirstStep || bLastStep);

    if (bNS && !(bFirstStep && ir->bContinuation)) {
      bMasterState = FALSE;
      /* Correct the new box if it is too skewed */
      if (DYNAMIC_BOX(*ir) && !bRerunMD) {
	if (correct_box(fplog,step,state->box,graph))
	  bMasterState = TRUE;
      }
      if (DOMAINDECOMP(cr) && bMasterState)
	dd_collect_state(cr->dd,state,state_global);
      
      if (DOMAINDECOMP(cr)) {
	/* Repartition the domain decomposition */
	wallcycle_start(wcycle,ewcDOMDEC);
	dd_partition_system(fplog,step,cr,bMasterState,
			    state_global,top_global,ir,
			    state,&f,&buf,mdatoms,top,fr,vsite,shellfc,constr,
			    nrnb,wcycle,do_verbose);
	wallcycle_stop(wcycle,ewcDOMDEC);
      }
    }

    if (MASTER(cr) && do_log && !bFFscan)
      print_ebin_header(fplog,step,t,state->lambda);

    if (ir->efep != efepNO)
      update_mdatoms(mdatoms,state->lambda); 
    
    clear_mat(force_vir);
    
    /* Ionize the atoms if necessary */
    if (bIonize)
      ionize(fplog,mdatoms,top->atoms.atomname,t,ir,state->x,state->v,
	     mdatoms->start,mdatoms->start+mdatoms->homenr,state->box,cr);
      
    /* Update force field in ffscan program */
    if (bFFscan) {
      if (update_forcefield(fplog,
			    nfile,fnm,fr,mdatoms->nr,state->x,state->box)) {
	if (gmx_parallel_env)
	  gmx_finalize(cr);
	exit(0);
      }
    }

    GMX_MPE_LOG(ev_timestep2);

    if (shellfc) {
      /* Now is the time to relax the shells */
      count=relax_shell_flexcon(fplog,cr,bVerbose,bFFscan ? step+1 : step,
				ir,bNS,bStopCM,top,constr,ener,fcd,
				state,f,buf,force_vir,mdatoms,
				nrnb,wcycle,graph,groups,grps,
				shellfc,fr,t,mu_tot,
				state->natoms,&bConverged,vsite,
				fp_field);
      tcount+=count;
      
      if (bConverged)
	nconverged++;
    }
    else {
      /* The coordinates (x) are shifted (to get whole molecules) in do_force
       * This is parallellized as well, and does communication too. 
       * Check comments in sim_util.c
       */
      do_force(fplog,cr,ir,step,nrnb,wcycle,top,groups,grps,
	       state->box,state->x,&state->hist,
	       f,buf,force_vir,mdatoms,ener,fcd,
	       state->lambda,graph,
	       fr,vsite,mu_tot,t,fp_field,ed,
	       GMX_FORCE_STATECHANGED | (bNS ? GMX_FORCE_NS : 0) |
	       GMX_FORCE_ALLFORCES);
    }

    GMX_BARRIER(cr->mpi_comm_mygroup);

    if (bTCR)
      mu_aver = calc_mu_aver(cr,state->x,mdatoms->chargeA,
			     mu_tot,top,mdatoms,gnx,grpindex);
    if (bGlas)
      do_glas(fplog,mdatoms->start,mdatoms->homenr,state->x,f,
	      fr,mdatoms,top->idef.atnr,ir,ener);
    
    if (bTCR && bFirstStep) {
      tcr=init_coupling(fplog,nfile,fnm,cr,fr,mdatoms,&(top->idef));
      fprintf(fplog,"Done init_coupling\n"); 
      fflush(fplog);
    }

    /* Now we have the energies and forces corresponding to the 
     * coordinates at time t. We must output all of this before
     * the update.
     * for RerunMD t is read from input trajectory
     */

    GMX_MPE_LOG(ev_output_start);

    bX   = do_per_step(step,ir->nstxout);
    bV   = do_per_step(step,ir->nstvout);
    bF   = do_per_step(step,ir->nstfout);
    bXTC = do_per_step(step,ir->nstxtcout);
    if ((bNS || bLastStep) && (step > ir->init_step) && !bRerunMD) {
      bCPT = (chkpt > 0 || bLastStep);
      chkpt = 0;
    } else {
      bCPT = FALSE;
    }

    if (bX || bV || bF || bXTC || bCPT) {
      wallcycle_start(wcycle,ewcTRAJ);
      if (sd && bCPT) {
	get_stochd_state(sd,state);
      }
      write_traj(fplog,cr,fp_trn,bX,bV,bF,fp_xtc,bXTC,ir->xtcprec,fn_cpt,bCPT,
		 top_global,ir->eI,step,t,state,state_global,f,f_global);
      debug_gmx();

      if (bLastStep && step_rel == ir->nsteps &&
	  (Flags & MD_CONFOUT) && MASTER(cr) &&
	  !bRerunMD && !bFFscan) {
	/* x and v have been collected in write_traj */
	fprintf(stderr,"\nWriting final coordinates.\n");
	write_sto_conf_mtop(ftp2fn(efSTO,nfile,fnm),
			    *top_global->name,top_global,
			    state_global->x,state_global->v,
			    ir->ePBC,state->box);
	debug_gmx();
      }
      wallcycle_stop(wcycle,ewcTRAJ);
    }
    GMX_MPE_LOG(ev_output_finish);

    clear_mat(shake_vir);

    if (bFFscan)
      clear_rvecs(state->natoms,buf);

    /* Box is changed in update() when we do pressure coupling,
     * but we should still use the old box for energy corrections and when
     * writing it to the energy file, so it matches the trajectory files for
     * the same timestep above. Make a copy in a separate array.
     */
    copy_mat(state->box,lastbox);
 
    
    GMX_MPE_LOG(ev_update_start);
    /* This is also parallellized, but check code in update.c */
    /* bOK = update(nsb->natoms,START(nsb),HOMENR(nsb),step,state->lambda,&ener[F_DVDL], */
    bOK = TRUE;
    if (!bRerunMD || rerun_fr.bV || bForceUpdate) {
      wallcycle_start(wcycle,ewcUPDATE);
      dvdl = 0;
      /* We can only do Berendsen coupling after we have summed the kinetic
       * energy or virial. Since the happens in global_state after update,
       * we should only do it at step % nstlist = 1 with bGStat=FALSE.
       */
      update(fplog,step,&dvdl,ir,mdatoms,state,graph,f,buf,fcd,
	     top,grps,shake_vir,scale_tot,
	     cr,nrnb,wcycle,sd,constr,bHaveConstr,
	     bNEMD,bFirstStep && bStateFromTPX);
      if (fr->bSepDVDL && fplog && do_log) {
	fprintf(fplog,sepdvdlformat,"Constraint",0.0,dvdl);
      }
      ener[F_DGDL_CON] += dvdl;
      wallcycle_stop(wcycle,ewcUPDATE);
    } else if (graph) {
      /* Need to unshift here */
      unshift_self(graph,state->box,state->x);
    }

    GMX_BARRIER(cr->mpi_comm_mygroup);
    GMX_MPE_LOG(ev_update_finish);

    if (!bOK && !bFFscan)
      gmx_fatal(FARGS,"Constraint error: Shake, Lincs or Settle could not solve the constrains");

    if (vsite) {
      wallcycle_start(wcycle,ewcVSITECONSTR);
      if (graph)
	shift_self(graph,state->box,state->x);
      
      construct_vsites(fplog,vsite,state->x,nrnb,ir->delta_t,state->v,
		       top->idef.iparams,top->idef.il,
		       fr->ePBC,fr->bMolPBC,graph,cr,state->box);
      
      if (graph)
	unshift_self(graph,state->box,state->x);
      wallcycle_stop(wcycle,ewcVSITECONSTR);
    }

    /* Non-equilibrium MD: this is parallellized, but only does communication
     * when there really is NEMD.
     */
    if (PAR(cr) && bNEMD) 
      accumulate_u(cr,&(ir->opts),grps);
      
    debug_gmx();
    if (grps->cosacc.cos_accel == 0)
      calc_ke_part(state->v,&(ir->opts),mdatoms,grps,nrnb,state->lambda);
    else
      calc_ke_part_visc(state->box,state->x,state->v,&(ir->opts),
			mdatoms,grps,nrnb,state->lambda);

    /* since we use the new coordinates in calc_ke_part_visc, we should use
     * the new box too. Still, won't this be offset by one timestep in the
     * energy file? / EL 20040121
     */ 

    debug_gmx();
    /* Calculate center of mass velocity if necessary, also parallellized */
    if (bStopCM && !bFFscan && !bRerunMD)
      calc_vcm_grp(fplog,mdatoms->start,mdatoms->homenr,mdatoms,
		   state->x,state->v,vcm);
    
    /* Check whether everything is still allright */    
    if (bGotTermSignal || bGotUsr1Signal) {
      if (bGotTermSignal || ir->nstlist == 0)
	terminate = 1;
      else
	terminate = -1;
      if (!PAR(cr))
	terminate_now = terminate;
      if (fplog) {
	fprintf(fplog,
		"\n\nReceived the %s signal, stopping at the next %sstep\n\n",
		bGotTermSignal ? "TERM" : "USR1",terminate==-1 ? "NS " : "");
	fflush(fplog);
      }
      fprintf(stderr,
	      "\n\nReceived the %s signal, stopping at the next %sstep\n\n",
	      bGotTermSignal ? "TERM" : "USR1",terminate==-1 ? "NS " : "");
      fflush(stderr);
      bGotTermSignal = FALSE;
      bGotUsr1Signal = FALSE;
    } else if (MASTER(cr) && (bNS || ir->nstlist <= 0) &&
	       (max_hours > 0 && (time(NULL) - start_t) > max_hours*24*60.0*0.99) &&
	       terminate == 0) {
      /* Signal to terminate the run */
      terminate = (ir->nstlist == 0 ? 1 : -1);
      if (!PAR(cr))
	terminate_now = terminate;
     if (fplog)
	fprintf(fplog,"\nStep %d: Run time exceeded %.3f hours, will terminate the run\n",step,max_hours*0.99);
      fprintf(stderr, "\nStep %d: Run time exceeded %.3f hours, will terminate the run\n",step,max_hours*0.99);
    }
    
    if (ir->nstlist == -1 && !bRerunMD) {
      nabnsb = natoms_beyond_ns_buffer(ir,fr,&top->cgs,*scale_tot,state->x);
    }

    if (MASTER(cr) && (cpt_period >= 0 &&
		       (cpt_period == 0 || 
			time(NULL) - start_t >= nchkpt*cpt_period*60.0))) {
      if (chkpt == 0) {
	nchkpt++;
      }
      chkpt = 1;
    }

    if (bGStat || bNS || do_per_step(step,ir->nstenergy) || bCPT || do_log) {
      if (PAR(cr)) {
	wallcycle_start(wcycle,ewcMoveE);
	/* Globally (over all NODEs) sum energy, virial etc. 
	 * This includes communication 
	 */
	global_stat(fplog,cr,ener,force_vir,shake_vir,mu_tot,
		    ir,grps,bSumEkinhOld,constr,vcm,
		    ir->nstlist==-1 ? &nabnsb : NULL,&chkpt,&terminate);
	if (terminate != 0) {
	  terminate_now = terminate;
	  terminate = 0;
	}
	
	wallcycle_stop(wcycle,ewcMoveE);
	bSumEkinhOld = FALSE;
      } else {
	bSumEkinhOld = TRUE;
      }

      /* This is just for testing. Nothing is actually done to Ekin
       * since that would require extra communication.
       */
      if (!bNEMD && debug && (vcm->nr > 0)) {
	correct_ekin(debug,mdatoms->start,mdatoms->start+mdatoms->homenr,
		     state->v,vcm->group_p[0],
		     mdatoms->massT,mdatoms->tmass,ekin);
      }
    
      /* Do center of mass motion removal */
      if (bStopCM && !bFFscan && !bRerunMD) {
	check_cm_grp(fplog,vcm,1);
	do_stopcm_grp(fplog,mdatoms->start,mdatoms->homenr,mdatoms->cVCM,
		      state->x,state->v,vcm);
	inc_nrnb(nrnb,eNR_STOPCM,mdatoms->homenr);
	/*
	  calc_vcm_grp(fplog,START(nsb),HOMENR(nsb),mdatoms->massT,x,v,vcm);
	  check_cm_grp(fplog,vcm);
	  do_stopcm_grp(fplog,START(nsb),HOMENR(nsb),x,v,vcm);
	  check_cm_grp(fplog,vcm);
	*/
      }
      
      /* Add force and shake contribution to the virial */
      m_add(force_vir,shake_vir,total_vir);
      
      /* Calculate the amplitude of the cosine velocity profile */
      grps->cosacc.vcos = grps->cosacc.mvcos/mdatoms->tmass;
      
      /* Sum the kinetic energies of the groups & calc temp */
      ener[F_TEMP] = sum_ekin(bRerunMD,&(ir->opts),grps,ekin,
			      &(ener[F_DKDL]));
      ener[F_EKIN] = trace(ekin);
      
      /* Calculate pressure and apply LR correction if PPPM is used.
       * Use the box from last timestep since we already called update().
       */
      ener[F_PRES] = calc_pres(fr->ePBC,ir->nwall,lastbox,ekin,total_vir,pres,
			       (fr->eeltype==eelPPPM)?ener[F_COUL_RECIP]:0.0);
      
      /* Calculate long range corrections to pressure and energy */
      if (bTCR || bFFscan) {
	set_avcsixtwelve(fplog,fr,top_global);
      }
      
      /* Calculate long range corrections to pressure and energy */
      calc_dispcorr(fplog,ir,fr,step,top_global->natoms,
		    lastbox,state->lambda,
		    pres,total_vir,ener);
      
      ener[F_ETOT] = ener[F_EPOT] + ener[F_EKIN];
      
      switch (ir->etc) {
      case etcNO:
      case etcBERENDSEN:
	break;
      case etcNOSEHOOVER:
	ener[F_ECONSERVED] =
	  ener[F_ETOT] + nosehoover_energy(&(ir->opts),grps,
					   state->nosehoover_xi,
					   state->therm_integral);
	break;
      case etcVRESCALE:
	ener[F_ECONSERVED] =
	  ener[F_ETOT] + vrescale_energy(&(ir->opts),state->therm_integral);
	break;
      }
      
      if ((state->flags & (1<<estPRES_PREV)) &&
	  (bGStat || ir->nstlist<=1 || step%ir->nstlist==0)) {
	/* Store the pressure in t_state for pressure coupling
	 * at the next MD step.
	 */
	copy_mat(pres,state->pres_prev);
      }

      /* Check for excessively large energies */
      if (bIonize) {
#ifdef GMX_DOUBLE
	real etot_max = 1e200;
#else
	real etot_max = 1e30;
#endif
	if (fabs(ener[F_ETOT]) > etot_max) {
	  fprintf(stderr,"Energy too large (%g), giving up\n",ener[F_ETOT]);
	  break;
	}
      }
    }
    
    /* The coordinates (x) were unshifted in update */
    if (bFFscan && (shellfc==NULL || bConverged)) {
      if (print_forcefield(fplog,ener,mdatoms->homenr,f,buf,xcopy,
			   &(top->mols),mdatoms->massT,pres)) {
	if (gmx_parallel_env)
	  gmx_finalize(cr);
	fprintf(stderr,"\n");
	exit(0);
      }
    }
    
    if (bTCR) {
      /* Only do GCT when the relaxation of shells (minimization) has converged,
       * otherwise we might be coupling to bogus energies. 
       * In parallel we must always do this, because the other sims might
       * update the FF.
       */
      
      /* Since this is called with the new coordinates state->x, I assume
       * we want the new box state->box too. / EL 20040121
       */
      do_coupling(fplog,nfile,fnm,tcr,t,step,ener,fr,
		  ir,MASTER(cr),
		  mdatoms,&(top->idef),mu_aver,
		  top->mols.nr,cr,
		  state->box,total_vir,pres,
		  mu_tot,state->x,f,bConverged);
      debug_gmx();
    }

    /* Time for performance */
    if (((step % stepout) == 0) || bLastStep)
      update_time();

    /* Output stuff */
    if (MASTER(cr)) {
      bool do_ene,do_dr,do_or;
      
      upd_mdebin(mdebin,fp_dgdl,bGStat,
		 mdatoms->tmass,step_rel,t,ener,state,lastbox,
		 shake_vir,force_vir,total_vir,pres,grps,mu_tot,constr);
      do_ene = do_per_step(step,ir->nstenergy) || bCPT || bFirstStep || bLastStep;
      do_dr  = do_per_step(step,ir->nstdisreout);
      do_or  = do_per_step(step,ir->nstorireout);
      print_ebin(fp_ene,do_ene,do_dr,do_or,do_log?fplog:NULL,
		 step,step_rel,t,
		 eprNORMAL,bCompact,mdebin,fcd,groups,&(ir->opts));
      if (ir->ePull != epullNO)
	pull_print_output(ir->pull,step,t);

      if (bVerbose)
	fflush(fplog);
    }
    
    /* Remaining runtime */
    if (MULTIMASTER(cr) && do_verbose) {
      if (shellfc)
	fprintf(stderr,"\n");
      print_time(stderr,start_t,step,ir);
    }

    /* Replica exchange */
    bExchanged = FALSE;
    if ((repl_ex_nst > 0) && (step > 0) && !bLastStep &&
	do_per_step(step,repl_ex_nst))
      bExchanged = replica_exchange(fplog,cr,repl_ex,state_global,ener,
				    state,step,t);
    if (bExchanged && PAR(cr)) {
      if (DOMAINDECOMP(cr)) {
	dd_partition_system(fplog,step,cr,TRUE,
			    state_global,top_global,ir,
			    state,&f,&buf,mdatoms,top,fr,vsite,shellfc,constr,
			    nrnb,wcycle,FALSE);
      } else {
	bcast_state(cr,state,FALSE);
      }
    }
    
    bFirstStep = FALSE;
    
    if (bRerunMD) 
      /* read next frame from input trajectory */
      bNotLastFrame = read_next_frame(status,&rerun_fr);

    if (!bRerunMD || !rerun_fr.bStep) {
      /* increase the MD step number */
      step++;
      step_rel++;
    }

    cycles = wallcycle_stop(wcycle,ewcSTEP);
    if (DOMAINDECOMP(cr) && wcycle)
      dd_cycles_add(cr->dd,cycles,ddCyclStep);
  }
  /* End of main MD loop */
  debug_gmx();

  if (bRerunMD)
    close_trj(status);

  if (!(cr->duty & DUTY_PME)) {
    /* Tell the PME only node to finish */
    gmx_pme_finish(cr);
  }
	  
  if (MASTER(cr)) {
    if (bGStat) {
      print_ebin(fp_ene,FALSE,FALSE,FALSE,fplog,step,step_rel,t,
		 eprAVER,FALSE,mdebin,fcd,groups,&(ir->opts));
      print_ebin(fp_ene,FALSE,FALSE,FALSE,fplog,step,step_rel,t,
		 eprRMS,FALSE,mdebin,fcd,groups,&(ir->opts));
    }
    close_enx(fp_ene);
    if (ir->nstxtcout)
      close_xtc(fp_xtc);
    close_trn(fp_trn);
    if (fp_dgdl)
      fclose(fp_dgdl);
    if (fp_field)
      fclose(fp_field);
  }
  debug_gmx();
  
  if (ir->nstlist == -1 && nns>0 && fplog) {
    fprintf(fplog,"Average neighborlist lifetime: %.1f steps, std.dev.: %.1f steps\n",ns_s1/nns,sqrt(ns_s2/nns - sqr(ns_s1/nns)));
    fprintf(fplog,"Average number of atoms that crossed the half buffer length: %.1f\n\n",ns_ab/nns);
  }

  if (shellfc && fplog) {
    fprintf(fplog,"Fraction of iterations that converged:           %.2f %%\n",
	    (nconverged*100.0)/step_rel);
    fprintf(fplog,"Average number of force evaluations per MD step: %.2f\n\n",
	    tcount/step_rel);
  }

  if (repl_ex_nst > 0 && MASTER(cr)) {
    print_replica_exchange_statistics(fplog,repl_ex);
  }
    
  return start_t;
}
