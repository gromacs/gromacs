/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
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

#ifdef GMX_MPI
#include <mpi.h>
#endif

#ifdef GMX_FAHCORE
#include "corewrap.h"
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
#ifdef HAVE_SIGUSR1
  case SIGUSR1:
    bGotUsr1Signal = TRUE;
    break;
#endif
  }
}

typedef struct { 
  gmx_integrator_t *func;
} gmx_intp_t;

/* The array should match the eI array in include/types/enums.h */
const gmx_intp_t integrator[eiNR] = { {do_md}, {do_steep}, {do_cg}, {do_md}, {do_md}, {do_nm}, {do_lbfgs}, {do_tpi}, {do_tpi}, {do_md} };

int  mdrunner(FILE *fplog,t_commrec *cr,int nfile,t_filenm fnm[],
	      bool bVerbose,bool bCompact,
	      ivec ddxyz,int dd_node_order,real rdd,real rconstr,
	      char *dddlb_opt,real dlb_scale,
	      char *ddcsx,char *ddcsy,char *ddcsz,
	      int nstepout,gmx_edsam_t ed,int repl_ex_nst,int repl_ex_seed,
	      real pforce,real cpt_period,real max_hours,
	      unsigned long Flags)
{
  double     nodetime=0,realtime;
  t_inputrec *inputrec;
  t_state    *state=NULL;
  matrix     box;
  rvec       *buf=NULL,*f=NULL;
  real       tmpr1,tmpr2;
  t_nrnb     *nrnb;
  gmx_mtop_t *mtop=NULL;
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
  bool       bReadRNG,bReadEkin;
  int        list;
  gmx_step_t nsteps_done;
  int        rc;
  gmx_genborn_t *born;

  snew(inputrec,1);
  snew(mtop,1);
  	
  if (bVerbose && SIMMASTER(cr)) 
    fprintf(stderr,"Getting Loaded...\n");
 
  if (Flags & MD_APPENDFILES) 
  {
	  fplog = NULL;
  }
	
  if (PAR(cr)) {
    /* The master thread on the master node reads from disk, 
     * then distributes everything to the other processors.
     */
	
    list = (SIMMASTER(cr) && !(Flags & MD_APPENDFILES)) ?  (LIST_SCALARS | LIST_INPUTREC) : 0;
		  
    snew(state,1);
    init_parallel(fplog,ftp2fn(efTPX,nfile,fnm),cr,
		  inputrec,mtop,state,list);
	  
  } else {
    /* Read a file for a single processor */
    snew(state,1);
    init_single(fplog,inputrec,ftp2fn(efTPX,nfile,fnm),mtop,state);
  }
  if (!EEL_PME(inputrec->coulombtype) || (Flags & MD_PARTDEC)) {
    cr->npmenodes = 0;
  }
 
  #ifdef GMX_FAHCORE
    fcRegisterSteps(inputrec->nsteps,inputrec->init_step);
  #endif
 
  /* NMR restraints must be initialized before load_checkpoint,
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
    if (SIMMASTER(cr)) {
      copy_mat(state->box,box);
    }
    if (PAR(cr)) {
      gmx_bcast(sizeof(box),box,cr);
    }
    set_deform_reference_box(inputrec->init_step,box);
  }

  if (opt2bSet("-cpi",nfile,fnm)) 
  {
      /* Check if checkpoint file exists before doing continuation.
       * This way we can use identical input options for the first and subsequent runs...
       */
      if( gmx_fexist(opt2fn("-cpi",nfile,fnm)) )
      {
	  load_checkpoint(opt2fn("-cpi",nfile,fnm),fplog,
			  cr,Flags & MD_PARTDEC,ddxyz,
			  inputrec,state,&bReadRNG,&bReadEkin,
			  (Flags & MD_APPENDFILES));
	  
	  if (bReadRNG)
	  {
	      Flags |= MD_READ_RNG;
	  }
	  if (bReadEkin)
	  {
	      Flags |= MD_READ_EKIN;
	  }
      }
  }
  
  if (MASTER(cr) && (Flags & MD_APPENDFILES))
  {
      fplog = gmx_log_open(ftp2fn(efLOG,nfile,fnm),cr, !(Flags & MD_SEPPOT) ,Flags);
  }
  
  if (SIMMASTER(cr)) 
  {
      copy_mat(state->box,box);
  }
	
  if (PAR(cr)) 
  {
      gmx_bcast(sizeof(box),box,cr);
  }
	
  if (bVerbose && SIMMASTER(cr))
    fprintf(stderr,"Loaded with Money\n\n");

  if (PAR(cr) && !((Flags & MD_PARTDEC) || EI_TPI(inputrec->eI))) {
    cr->dd = init_domain_decomposition(fplog,cr,Flags,ddxyz,rdd,rconstr,
				       dddlb_opt,dlb_scale,
				       ddcsx,ddcsy,ddcsz,
				       mtop,inputrec,
				       box,state->x);
    
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
    
    /* Initialise GB-stuff */
    if(fr->bGB)
      init_gb(&born,cr,fr,inputrec,mtop,state->x,inputrec->rgbradii,inputrec->gb_algorithm);

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
	do_pbc_first_mtop(fplog,inputrec->ePBC,box,mtop,state->x);
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
#ifdef HAVE_SIGUSR1
    if (getenv("GMX_NO_USR1") == NULL) {
        if (debug)
            fprintf(debug,"Installing signal handler for SIGUSR1\n");
        signal(SIGUSR1,signal_handler);
    }
#endif
  }

  if (cr->duty & DUTY_PP) {
    if (inputrec->ePull != epullNO) {
      /* Initialize pull code */
      init_pull(fplog,inputrec,nfile,fnm,mtop,cr,
		EI_DYNAMICS(inputrec->eI) && MASTER(cr),Flags);
    }

    constr = init_constraints(fplog,mtop,inputrec,ed,state,cr);

    if (DOMAINDECOMP(cr)) {
      dd_init_bondeds(fplog,cr->dd,mtop,vsite,constr,inputrec,
		      Flags & MD_DDBONDCHECK,fr->cginfo_global);

      set_dd_parameters(fplog,cr->dd,dlb_scale,inputrec,fr,box);
     
      setup_dd_grid(fplog,cr->dd);
    }

    /* Now do whatever the user wants us to do (how flexible...) */
    start_t = integrator[inputrec->eI].func(fplog,cr,nfile,fnm,
					    bVerbose,bCompact,
					    vsite,constr,
					    nstepout,inputrec,mtop,
					    fcd,state,f,buf,
					    mdatoms,nrnb,wcycle,ed,fr,
					    born,
					    repl_ex_nst,repl_ex_seed,
					    cpt_period,max_hours,
					    Flags,
					    &nsteps_done);

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
	     inputrec,nrnb,wcycle,nodetime,realtime,nsteps_done,
	     EI_DYNAMICS(inputrec->eI) && !MULTISIM(cr));
  
  /* Does what it says */  
  print_date_and_time(fplog,cr->nodeid,"Finished mdrun");

	/* Close logfile already here if we were appending to it */
	if (MASTER(cr) && (Flags & MD_APPENDFILES))
	{
		gmx_log_close(fplog);
	}	

  if(bGotTermSignal)
  {
    rc = 1;
  }
  else if(bGotUsr1Signal)
  {
    rc = 2;
  }
  else
  {
    rc = 0;
  }
	
  return rc;

}

time_t do_md(FILE *fplog,t_commrec *cr,int nfile,t_filenm fnm[],
	     bool bVerbose,bool bCompact,
	     gmx_vsite_t *vsite,gmx_constr_t constr,
	     int stepout,t_inputrec *ir,
	     gmx_mtop_t *top_global,
	     t_fcdata *fcd,
	     t_state *state_global,rvec f[],
	     rvec buf[],t_mdatoms *mdatoms,
	     t_nrnb *nrnb,gmx_wallcycle_t wcycle,
	     gmx_edsam_t ed,t_forcerec *fr, gmx_genborn_t *born,
	     int repl_ex_nst,int repl_ex_seed,
	     real cpt_period,real max_hours,
	     unsigned long Flags,
	     gmx_step_t *nsteps_done)
{
  int        fp_ene=0,fp_trn=0,fp_xtc=0;
  gmx_step_t step,step_rel;
  char       *fn_cpt;
  FILE       *fp_dhdl=NULL,*fp_field=NULL;
  time_t     start_t;
  double     run_time;
  double     t,t0,lam0;
  bool       bGStatEveryStep,bGStat,bCalcPres;
  bool       bNS,bSimAnn,bStopCM,bRerunMD,bNotLastFrame=FALSE,
             bFirstStep,bStateFromTPX,bLastStep,bBornRadii;
  bool       bDoDHDL=FALSE;
  bool       bNEMD,do_ene,do_log,do_verbose,bRerunWarnNoV=TRUE,
	         bForceUpdate=FALSE,bX,bV,bF,bXTC,bCPT;
  bool       bMasterState;
  tensor     force_vir,shake_vir,total_vir,pres,ekin;
  int        i,m,status;
  rvec       mu_tot;
  t_vcm      *vcm;
  gmx_step_t step_ns=0,step_nscheck=0,nns=0,ns_lt;
  int        nabnsb=0;
  double     ns_s1=0,ns_s2=0,ns_ab=0,ns_lt_runav=0,ns_lt_runav2=0;
  matrix     *scale_tot;
  t_trxframe rerun_fr;
  gmx_repl_ex_t repl_ex=NULL;
  int        nchkpt=1;
  /* Booleans (disguised as a reals) to checkpoint and terminate mdrun */  
  real       chkpt=0,terminate=0,terminate_now=0;

  gmx_localtop_t *top;	
  t_mdebin *mdebin=NULL;
  t_state    *state=NULL;
  rvec       *f_global=NULL;
  gmx_enerdata_t *enerd;
  gmx_stochd_t sd=NULL;
  t_graph    *graph=NULL;

  bool        bFFscan;
  gmx_groups_t *groups;
  gmx_ekindata_t *ekind;
  gmx_shellfc_t shellfc;
  int         count,nconverged=0;
  real        timestep=0;
  double      tcount=0;
  bool        bIonize=FALSE;
  bool        bTCR=FALSE,bConverged=TRUE,bOK,bSumEkinhOld,bExchanged;
  bool        bAppend;
  real        temp0,mu_aver=0,dvdl;
  int         a0,a1,gnx=0,ii;
  atom_id     *grpindex=NULL;
  char        *grpname;
  t_coupl_rec *tcr=NULL;
  rvec        *xcopy=NULL,*vcopy=NULL;
  matrix      boxcopy,lastbox;
  double      cycles;
  char        sbuf[22],sbuf2[22];
#ifdef GMX_FAHCORE
  /* Temporary addition for FAHCORE checkpointing */
  int chkpt_ret;
#endif

  /* Check for special mdrun options */
  bRerunMD = (Flags & MD_RERUN);
  bIonize  = (Flags & MD_IONIZE);
  bFFscan  = (Flags & MD_FFSCAN);
  bGStatEveryStep = !(Flags & MD_NOGSTAT);
  bAppend  = (Flags & MD_APPENDFILES);

  /* If we do reruns, the step numbers in the output energy frames cannot be
   * used for averages (since energies are only calculated for trajectory frames).
   * By turning of bGStatEveryStep we force g_energy to use the actual energy frame 
   * contents for the averages instead.
   */
  if(bRerunMD)
  {
      bGStatEveryStep = FALSE;
  }
	
  if (!bGStatEveryStep && !EI_DYNAMICS(ir->eI)) {
    char *warn="\nWARNING:\nNo energy summing can only be used with dynamics, ignoring this option\n";
    fprintf(stderr,"%s\n",warn);
    if (fplog)
      fprintf(fplog,"%s\n",warn);
    bGStatEveryStep = TRUE;
  }
  if (!bGStatEveryStep) {
    if (fplog) {
      fprintf(fplog,"\nWill not sum the energies at every step,\n"
	      "therefore the energy file does not contain exact averages and fluctuations.\n\n");
      if (ir->etc != etcNO || ir->epc != epcNO) {
	fprintf(fplog,"WARNING:\nThe temperature and/or pressure for scaling will only be updated every nstlist (%d) steps\n\n",ir->nstlist);
      }
      if (ir->nstlist == -1) {
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
    }
    if (ir->comm_mode != ecmNO && ir->nstcomm == 1) {
      if (fplog) {
	fprintf(fplog,"WARNING:\nWe should not remove the COM motion every step with option -nosum,\n");
      }
      if (ir->nstlist > 0) {
	ir->nstcomm = ir->nstlist;
	if (fplog) {
	  fprintf(fplog,"setting nstcomm to nstlist (%d)\n\n",ir->nstcomm);
	}
      } else {
	ir->nstcomm = 100;
	if (fplog) {
	  fprintf(fplog,"setting nstcomm to %d\n\n",ir->nstcomm);
	}
      }
    }
  }

  if (bRerunMD || bFFscan)
    ir->nstxtcout = 0;

  groups = &top_global->groups;

  /* Initial values */
  init_md(fplog,cr,ir,&t,&t0,&state_global->lambda,&lam0,
	  nrnb,top_global,&sd,
	  nfile,fnm,&fp_trn,&fp_xtc,&fp_ene,&fn_cpt,
	  &fp_dhdl,&fp_field,&mdebin,
	  force_vir,shake_vir,mu_tot,&bNEMD,&bSimAnn,&vcm,Flags);

  /* Energy terms and groups */
  snew(enerd,1);
  init_enerdata(top_global->groups.grps[egcENER].nr,ir->n_flambda,enerd);
  /* Kinetic energy data */
  snew(ekind,1);
  init_ekindata(fplog,top_global,&(ir->opts),ekind);
  /* Copy the cos acceleration to the groups struct */
  ekind->cosacc.cos_accel = ir->cos_accel;

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
      graph = mk_graph(fplog,&(top->idef),0,top_global->natoms,FALSE,FALSE);
    }

    if (shellfc) {
      make_local_shells(cr,mdatoms,shellfc);
    }

    if (ir->pull && PAR(cr)) {
      dd_make_local_pull_groups(NULL,ir->pull,mdatoms);
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
	
	
  if (sd && (Flags & MD_READ_RNG)) {
    /* Set the random state if we read a checkpoint file */
    set_stochd_state(sd,state);
  }

  /* Initialize constraints */
  if (constr) {
    if (!DOMAINDECOMP(cr))
      set_constraints(constr,top,ir,mdatoms,NULL);
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
            for(i=mdatoms->start; i<mdatoms->homenr; i++)
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
            do_constrain_first(fplog,constr,ir,mdatoms,state,
                               graph,cr,nrnb,fr,&top->idef);
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

  if (Flags & MD_READ_EKIN)
  {
      restore_ekinstate_from_state(cr,ekind,&state_global->ekinstate);
  }
  else
  {
      /* Compute initial EKin for all.. */
      if (ekind->cosacc.cos_accel == 0) {
          calc_ke_part(state->v,&(ir->opts),mdatoms,ekind,nrnb,state->lambda);
      } else {
	  calc_ke_part_visc(state->box,state->x,state->v,&(ir->opts),
			    mdatoms,ekind,nrnb,state->lambda);
      }
      debug_gmx();

      if (PAR(cr))
      {
	  GMX_MPE_LOG(ev_global_stat_start);
	  
	  global_stat(fplog,cr,enerd,force_vir,shake_vir,mu_tot,
		      ir,ekind,FALSE,constr,vcm,NULL,NULL,&terminate,
		      top_global,state);
	  
	  GMX_MPE_LOG(ev_global_stat_finish);
      }
       debug_gmx();
  }
  
  /* Calculate the initial half step temperature */
  temp0 = sum_ekin(TRUE,&(ir->opts),ekind,ekin,NULL);

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

    if (MASTER(cr))
    {
        if (constr && !ir->bContinuation && ir->eConstrAlg == econtLINCS)
        {
            fprintf(fplog,
                    "RMS relative constraint deviation after constraining: %.2e\n",
                    constr_rmsd(constr,FALSE));
        }
        fprintf(fplog,"Initial temperature: %g K\n",temp0);
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
            fprintf(stderr,"starting mdrun '%s'\n",
                    *(top_global->name));
            if (ir->init_step > 0)
            {
                fprintf(stderr,"%s steps, %8.1f ps (continuing from step %s, %8.1f ps).\n",
                        gmx_step_str(ir->init_step+ir->nsteps,sbuf),
                        (ir->init_step+ir->nsteps)*ir->delta_t,
                        gmx_step_str(ir->init_step,sbuf2),
                        ir->init_step*ir->delta_t);
            }
            else
            {
                fprintf(stderr,"%s steps, %8.1f ps.\n",
                        gmx_step_str(ir->nsteps,sbuf),ir->nsteps*ir->delta_t);
            }
        }
        fprintf(fplog,"\n");
    }

  if (ir->nstlist == -1) {
    snew(scale_tot,1);
  } else {
    scale_tot = NULL;
  }

  /* Write start time */
  start_t=print_date_and_time(fplog,cr->nodeid,"Started mdrun");
  wallcycle_start(wcycle,ewcRUN);
  if (fplog)
    fprintf(fplog,"\n");

  /* safest point to do file checkpointing is here.  More general point would be immediately before integrator call */
  #ifdef GMX_FAHCORE
	chkpt_ret=fcCheckPointParallel( (MASTER(cr)==0),
             NULL,0);
  if ( chkpt_ret == 0 ) 
        gmx_fatal( 3,__FILE__,__LINE__, "Checkpoint error on step %d\n", 0 );
  #endif
			 
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
    if (rerun_fr.natoms != top_global->natoms)
      gmx_fatal(FARGS,"Number of atoms in trajectory (%d) does not match the "
		"run input file (%d)\n",rerun_fr.natoms,top_global->natoms);
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

        if (ir->efep != efepNO)
        {
            if (bRerunMD && rerun_fr.bLambda && (ir->delta_lambda!=0))
            {
                state->lambda = rerun_fr.lambda;
            }
            else
            {
                state->lambda = lam0 + step*ir->delta_lambda;
            }
            bDoDHDL = do_per_step(step,ir->nstdhdl);
        }
        
        if (bSimAnn) 
        {
            update_annealing_target_temp(&(ir->opts),t);
        }
        
        if (bRerunMD)
        {
            for(i=0; i<mdatoms->nr; i++)
            {
                copy_rvec(rerun_fr.x[i],state->x[i]);
            }
            if (rerun_fr.bV)
            {
                for(i=0; i<mdatoms->nr; i++)
                {
                    copy_rvec(rerun_fr.v[i],state->v[i]);
                }
            }
            else
            {
                for(i=0; i<mdatoms->nr; i++)
                {
                    clear_rvec(state->v[i]);
                }
                if (bRerunWarnNoV)
                {
                    fprintf(stderr,"\nWARNING: Some frames do not contain velocities.\n"
                            "         Ekin, temperature and pressure are incorrect,\n"
                            "         the virial will be incorrect when constraints are present.\n"
                            "\n");
                    bRerunWarnNoV = FALSE;
                }
            }
            copy_mat(rerun_fr.box,state->box);

            if (vsite)
            {
                if (graph)
                {
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
                {
                    unshift_self(graph,state->box,state->x);
                }
            }
        }
        
        /* Stop Center of Mass motion */
        bStopCM = (ir->comm_mode != ecmNO && do_per_step(step,ir->nstcomm));
        
        /* Copy back starting coordinates in case we're doing a forcefield scan */
        if (bFFscan)
        {
            for(ii=0; (ii<state->natoms); ii++)
            {
                copy_rvec(xcopy[ii],state->x[ii]);
                copy_rvec(vcopy[ii],state->v[ii]);
            }
            copy_mat(boxcopy,state->box);
        }
        
        bNS = bFirstStep;
        if (bRerunMD)
        {
            /* for rerun MD always do Neighbour Searching */
            if (ir->nstlist != 0)
            {
                bNS = TRUE;
            }
        }
        else
        {
            /* Determine whether or not to do Neighbour Searching */
            if (bExchanged || (ir->nstlist > 0 && (step % ir->nstlist == 0)))
            {
                bNS = TRUE;
            }
            else if (ir->nstlist == -1)
            {
                bNS = (bFirstStep || nabnsb > 0);
                if (bNS)
                {
                    if (bFirstStep)
                    {
                        ns_lt_runav = 0;
                        step_nscheck = step;
                    }
                    else
                    {
                        /* Determine the neighbor list life time */
                        ns_lt = step - step_ns;
                        if (debug)
                        {
                            fprintf(debug,"%d atoms beyond ns buffer, updating neighbor list after %s steps\n",nabnsb,gmx_step_str(ns_lt,sbuf));
                        }
                        nns++;
                        ns_s1 += ns_lt;
                        ns_s2 += ns_lt*ns_lt;
                        ns_ab += nabnsb;
                        if (ns_lt_runav == 0)
                        {
                            ns_lt_runav  = ns_lt;
                            /* Initialize the fluctuation average
                             * such that at startup we check after 0 steps.
                             */
                            ns_lt_runav2 = sqr(ns_lt/2.0);
                        }
                        /* Running average with 0.9 gives an exp. history of 9.5 */
                        ns_lt_runav2 = 0.9*ns_lt_runav2 + 0.1*sqr(ns_lt_runav - ns_lt);
                        ns_lt_runav  = 0.9*ns_lt_runav  + 0.1*ns_lt;
                        if (bGStatEveryStep)
                        {
                            /* Always check the nlist validity */
                            step_nscheck = step;
                        }
                        else
                        {
                            /* We check after:  <life time> - 2*sigma
                             * The factor 2 is quite conservative,
                             * but we assume that with nstlist=-1 the user
                             * prefers exact integration over performance.
                             */
                            step_nscheck = step
                                + (int)(ns_lt_runav - 2.0*sqrt(ns_lt_runav2)) - 1;
                        }
                        if (debug)
                        {
                            fprintf(debug,"nlist life time %s run av. %4.1f sig %3.1f check %s check with -nosum %d\n",
                                    gmx_step_str(ns_lt,sbuf),ns_lt_runav,sqrt(ns_lt_runav2),
                                    gmx_step_str(step_nscheck-step+1,sbuf2),
                                    (int)(ns_lt_runav - 2.0*sqrt(ns_lt_runav2)));
                        }
                    }
                    step_ns = step;
                    /* Initialize the cumulative coordinate scaling matrix */
                    clear_mat(*scale_tot);
                    for(ii=0; ii<DIM; ii++)
                    {
                        (*scale_tot)[ii][ii] = 1.0;
                    }
                }
            }
        } 
        
        if (terminate_now > 0 || (terminate_now < 0 && bNS))
        {
            bLastStep = TRUE;
        }

        /* Determine whether or not to update the Born radii if doing GB */
        bBornRadii=bFirstStep;
        if(ir->implicit_solvent && (step % ir->nstgbradii==0))
        {
            bBornRadii=TRUE;
        }
        do_log = do_per_step(step,ir->nstlog) || bFirstStep || bLastStep;
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
                    if (correct_box(fplog,step,state->box,graph))
                    {
                        bMasterState = TRUE;
                    }
                }
                if (DOMAINDECOMP(cr) && bMasterState)
                {
                    dd_collect_state(cr->dd,state,state_global);
                }
            }
            
            if (DOMAINDECOMP(cr))
            {
                /* Repartition the domain decomposition */
                wallcycle_start(wcycle,ewcDOMDEC);
                dd_partition_system(fplog,step,cr,bMasterState,
                                    state_global,top_global,ir,
                                    state,&f,&buf,mdatoms,top,fr,
                                    vsite,shellfc,constr,
                                    nrnb,wcycle,do_verbose);
                wallcycle_stop(wcycle,ewcDOMDEC);
            }
        }
        
        if (MASTER(cr) && do_log && !bFFscan)
        {
            print_ebin_header(fplog,step,t,state->lambda);
        }
        
        if (ir->efep != efepNO)
        {
            update_mdatoms(mdatoms,state->lambda); 
        }
        
        clear_mat(force_vir);
        
        /* Ionize the atoms if necessary */
        if (bIonize)
        {
            ionize(fplog,mdatoms,top_global,t,ir,state->x,state->v,
                   mdatoms->start,mdatoms->start+mdatoms->homenr,state->box,cr);
        }
        
        /* Update force field in ffscan program */
        if (bFFscan)
        {
            if (update_forcefield(fplog,
                                  nfile,fnm,fr,
                                  mdatoms->nr,state->x,state->box)) {
                if (gmx_parallel_env)
                {
                    gmx_finalize(cr);
                }
                exit(0);
            }
        }
        
        GMX_MPE_LOG(ev_timestep2);
        
        if ((bNS || bLastStep) && (step > ir->init_step) && !bRerunMD)
        {
            bCPT = (chkpt > 0 || bLastStep);
            if (bCPT)
            {
                chkpt = 0;
            }
        }
        else
        {
            bCPT = FALSE;
        }

        /* Determine the pressure:
         * always when we want exact averages in the energy file,
         * at ns steps when we have pressure coupling,
         * otherwise only at energy output steps (set below).
         */
        bCalcPres = (bGStatEveryStep || (ir->epc != epcNO && bNS));
        
        /* Do we need global communication ? */
        bGStat = (bGStatEveryStep || bStopCM || bNS ||
                  (ir->nstlist == -1 && !bRerunMD && step >= step_nscheck));
        
        do_ene = (do_per_step(step,ir->nstenergy) || bLastStep);
        
        if (do_ene || do_log)
        {
            bCalcPres = TRUE;
            bGStat    = TRUE;
        }
        
        if (shellfc)
        {
            /* Now is the time to relax the shells */
            count=relax_shell_flexcon(fplog,cr,bVerbose,bFFscan ? step+1 : step,
                                      ir,bNS,bStopCM,top,top_global,
                                      constr,enerd,fcd,
                                      state,f,buf,force_vir,mdatoms,
                                      nrnb,wcycle,graph,groups,
                                      shellfc,fr,born,bBornRadii,t,mu_tot,
                                      state->natoms,&bConverged,vsite,
                                      fp_field);
            tcount+=count;
            
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
            do_force(fplog,cr,ir,step,nrnb,wcycle,top,top_global,groups,
                     state->box,state->x,&state->hist,
                     f,buf,force_vir,mdatoms,enerd,fcd,
                     state->lambda,graph,
                     fr,vsite,mu_tot,t,fp_field,ed,born,bBornRadii,
                     GMX_FORCE_STATECHANGED | (bNS ? GMX_FORCE_NS : 0) |
                     GMX_FORCE_ALLFORCES | (bCalcPres ? GMX_FORCE_VIRIAL : 0) |
                     (bDoDHDL ? GMX_FORCE_DHDL : 0));
        }
        
        GMX_BARRIER(cr->mpi_comm_mygroup);
        
        if (bTCR)
        {
            mu_aver = calc_mu_aver(cr,state->x,mdatoms->chargeA,
                                   mu_tot,&top_global->mols,mdatoms,gnx,grpindex);
        }
        
        if (bTCR && bFirstStep)
        {
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
        
#ifdef GMX_FAHCORE
        bX = bX || bLastStep; /*enforce writing positions and velocities at end of run */
        bV = bV || bLastStep;
        bCPT = bCPT || fcCheckPointPending();  
        if (bCPT)
            fcRequestCheckPoint(); /* sync bCPT and fc record-keeping */
        
        if (MASTER(cr))
            fcReportProgress( ir->nsteps, step );
#endif
        
        if (bX || bV || bF || bXTC || bCPT)
        {
            wallcycle_start(wcycle,ewcTRAJ);
            if (bCPT)
            {
                if (sd)
                {
                    get_stochd_state(sd,state);
                }
                if (MASTER(cr))
                {
                    update_ekinstate(&state_global->ekinstate,ekind);
                    update_energyhistory(&state_global->enerhist,mdebin);
                }
            }
            write_traj(fplog,cr,fp_trn,bX,bV,bF,fp_xtc,bXTC,ir->xtcprec,fn_cpt,bCPT,
                       top_global,ir->eI,ir->simulation_part,step,t,state,state_global,f,f_global);
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
        
        clear_mat(shake_vir);
        
        if (bFFscan)
        {
            clear_rvecs(state->natoms,buf);
        }
        
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
        if (!bRerunMD || rerun_fr.bV || bForceUpdate)
        {
            wallcycle_start(wcycle,ewcUPDATE);
            dvdl = 0;
            /* We can only do Berendsen coupling after we have summed
             * the kinetic energy or virial.
             * Since the happens in global_state after update,
             * we should only do it at step % nstlist = 1
             * with bGStatEveryStep=FALSE.
             */
            update(fplog,step,&dvdl,ir,mdatoms,state,graph,f,buf,fcd,
                   &top->idef,ekind,scale_tot,
                   cr,nrnb,wcycle,sd,constr,bCalcPres,shake_vir,
                   bNEMD,bFirstStep && bStateFromTPX);
            if (fr->bSepDVDL && fplog && do_log)
            {
                fprintf(fplog,sepdvdlformat,"Constraint",0.0,dvdl);
            }
            enerd->term[F_DHDL_CON] += dvdl;
            wallcycle_stop(wcycle,ewcUPDATE);
        }
        else if (graph)
        {
            /* Need to unshift here */
            unshift_self(graph,state->box,state->x);
        }
        
        GMX_BARRIER(cr->mpi_comm_mygroup);
        GMX_MPE_LOG(ev_update_finish);
        
        if (!bOK && !bFFscan)
        {
            gmx_fatal(FARGS,"Constraint error: Shake, Lincs or Settle could not solve the constrains");
        }
        
        if (vsite != NULL)
        {
            wallcycle_start(wcycle,ewcVSITECONSTR);
            if (graph != NULL)
            {
                shift_self(graph,state->box,state->x);
            }
            construct_vsites(fplog,vsite,state->x,nrnb,ir->delta_t,state->v,
                             top->idef.iparams,top->idef.il,
                             fr->ePBC,fr->bMolPBC,graph,cr,state->box);
            
            if (graph != NULL)
            {
                unshift_self(graph,state->box,state->x);
            }
            wallcycle_stop(wcycle,ewcVSITECONSTR);
        }
        
        /* Non-equilibrium MD:  this is parallellized,
         * but only does communication when there really is NEMD.
         */
        if (PAR(cr) && bNEMD) 
        {
            accumulate_u(cr,&(ir->opts),ekind);
        }
        
        debug_gmx();
        if (ekind->cosacc.cos_accel == 0)
        {
            calc_ke_part(state->v,&(ir->opts),mdatoms,ekind,nrnb,state->lambda);
        }
        else
        {
            calc_ke_part_visc(state->box,state->x,state->v,&(ir->opts),
                              mdatoms,ekind,nrnb,state->lambda);
        }
        
        /* since we use the new coordinates in calc_ke_part_visc, we should use
         * the new box too. Still, won't this be offset by one timestep in the
         * energy file? / EL 20040121
         */ 
        
        debug_gmx();
        /* Calculate center of mass velocity if necessary, also parallellized */
        if (bStopCM && !bFFscan && !bRerunMD)
        {
            calc_vcm_grp(fplog,mdatoms->start,mdatoms->homenr,mdatoms,
                         state->x,state->v,vcm);
        }
        
        /* Determine the wallclock run time up till now */
        run_time = (double)time(NULL) - (double)start_t;
        
        /* Check whether everything is still allright */    
        if (bGotTermSignal || bGotUsr1Signal)
        {
            if (bGotTermSignal || ir->nstlist == 0)
            {
                terminate = 1;
            }
            else
            {
                terminate = -1;
            }
            if (!PAR(cr))
            {
                terminate_now = terminate;
            }
            if (fplog)
            {
                fprintf(fplog,
                        "\n\nReceived the %s signal, stopping at the next %sstep\n\n",
                        bGotTermSignal ? "TERM" : "USR1",
                        terminate==-1 ? "NS " : "");
                fflush(fplog);
            }
            fprintf(stderr,
                    "\n\nReceived the %s signal, stopping at the next %sstep\n\n",
                    bGotTermSignal ? "TERM" : "USR1",
                    terminate==-1 ? "NS " : "");
            fflush(stderr);
            bGotTermSignal = FALSE;
            bGotUsr1Signal = FALSE;
        }
        else if (MASTER(cr) && (bNS || ir->nstlist <= 0) &&
                 (max_hours > 0 && run_time > max_hours*60.0*60.0*0.99) &&
                 terminate == 0)
        {
            /* Signal to terminate the run */
            terminate = (ir->nstlist == 0 ? 1 : -1);
            if (!PAR(cr))
            {
                terminate_now = terminate;
            }
            if (fplog)
            {
                fprintf(fplog,"\nStep %s: Run time exceeded %.3f hours, will terminate the run\n",gmx_step_str(step,sbuf),max_hours*0.99);
            }
            fprintf(stderr, "\nStep %s: Run time exceeded %.3f hours, will terminate the run\n",gmx_step_str(step,sbuf),max_hours*0.99);
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
            if (step >= step_nscheck)
            {
                nabnsb = natoms_beyond_ns_buffer(ir,fr,&top->cgs,
                                                 *scale_tot,state->x);
            }
            else
            {
                /* This is not necessarily true,
                 * but step_nscheck is determined quite conservatively.
                 */
                nabnsb = 0;
            }
        }
        
        /* In parallel we only have to check for checkpointing in steps
         * where we do global communication,
         *  otherwise the other nodes don't know.
         */
        if (MASTER(cr) && ((bGStat || !PAR(cr)) &&
                           cpt_period >= 0 &&
                           (cpt_period == 0 || 
                            run_time >= nchkpt*cpt_period*60.0)))
        {
            if (chkpt == 0)
            {
                nchkpt++;
            }
            chkpt = 1;
        }
        
        if (!bGStat)
        {
            /* We will not sum ekinh_old,
             * so signal that we still have to do it.
             */
            bSumEkinhOld = TRUE;
        }
        else
        {
            if (PAR(cr))
            {
                wallcycle_start(wcycle,ewcMoveE);
                /* Globally (over all NODEs) sum energy, virial etc. 
                 * This includes communication 
                 */
                global_stat(fplog,cr,enerd,force_vir,shake_vir,mu_tot,
                            ir,ekind,bSumEkinhOld,constr,vcm,
                            ir->nstlist==-1 ? &nabnsb : NULL,&chkpt,&terminate,
                            top_global, state);
                if (terminate != 0)
                {
                    terminate_now = terminate;
                    terminate = 0;
                }
                
                wallcycle_stop(wcycle,ewcMoveE);
                bSumEkinhOld = FALSE;
            }
            
            /* This is just for testing. Nothing is actually done to Ekin
             * since that would require extra communication.
             */
            if (!bNEMD && debug && (vcm->nr > 0))
            {
                correct_ekin(debug,
                             mdatoms->start,mdatoms->start+mdatoms->homenr,
                             state->v,vcm->group_p[0],
                             mdatoms->massT,mdatoms->tmass,ekin);
            }
            
            /* Do center of mass motion removal */
            if (bStopCM && !bFFscan && !bRerunMD)
            {
                check_cm_grp(fplog,vcm,ir,1);
                do_stopcm_grp(fplog,mdatoms->start,mdatoms->homenr,mdatoms->cVCM,
                              state->x,state->v,vcm);
                inc_nrnb(nrnb,eNR_STOPCM,mdatoms->homenr);
                /*
                  calc_vcm_grp(fplog,START(nsb),HOMENR(nsb),mdatoms->massT,x,v,vcm);
                  check_cm_grp(fplog,vcm,ir);
                  do_stopcm_grp(fplog,START(nsb),HOMENR(nsb),x,v,vcm);
                  check_cm_grp(fplog,vcm,ir);
                */
            }
            
            /* Add force and shake contribution to the virial */
            m_add(force_vir,shake_vir,total_vir);
            
            /* Calculate the amplitude of the cosine velocity profile */
            ekind->cosacc.vcos = ekind->cosacc.mvcos/mdatoms->tmass;
            
            /* Sum the kinetic energies of the groups & calc temp */
            enerd->term[F_TEMP] = sum_ekin(bRerunMD,&(ir->opts),ekind,ekin,
                                           &(enerd->term[F_DKDL]));
            enerd->term[F_EKIN] = trace(ekin);
            
            /* Calculate pressure and apply LR correction if PPPM is used.
             * Use the box from last timestep since we already called update().
             */
            enerd->term[F_PRES] =
                calc_pres(fr->ePBC,ir->nwall,lastbox,ekin,total_vir,pres,
                          (fr->eeltype==eelPPPM)?enerd->term[F_COUL_RECIP]:0.0);
            
            /* Calculate long range corrections to pressure and energy */
            if (bTCR || bFFscan)
            {
                set_avcsixtwelve(fplog,fr,top_global);
            }
            
            /* Calculate long range corrections to pressure and energy */
            calc_dispcorr(fplog,ir,fr,step,top_global->natoms,
                          lastbox,state->lambda,pres,total_vir,enerd);
            
            sum_dhdl(enerd,state->lambda,ir);
            
            enerd->term[F_ETOT] = enerd->term[F_EPOT] + enerd->term[F_EKIN];
            
            switch (ir->etc) {
            case etcNO:
            case etcBERENDSEN:
                break;
            case etcNOSEHOOVER:
                enerd->term[F_ECONSERVED] =
                    enerd->term[F_ETOT] +
                    nosehoover_energy(&(ir->opts),ekind,
                                      state->nosehoover_xi,
                                      state->therm_integral);
                break;
            case etcVRESCALE:
                enerd->term[F_ECONSERVED] =
                    enerd->term[F_ETOT] +
                    vrescale_energy(&(ir->opts),
                                    state->therm_integral);
                break;
            }
            
            /* Complicated conditional when bGStatEveryStep=FALSE.
             * We can not just use bCalcPres, since then the simulation results
             * would depend on nstenergy and nstlog or step_nscheck.
             */
            if ((state->flags & (1<<estPRES_PREV)) &&
                (bGStatEveryStep || bNS || (ir->nstlist == 0 && bCalcPres)))
            {
                /* Store the pressure in t_state for pressure coupling
                 * at the next MD step.
                 */
                copy_mat(pres,state->pres_prev);
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
                    fprintf(stderr,"Energy too large (%g), giving up\n",
                            enerd->term[F_ETOT]);
                    break;
                }
            }
        }
        
        /* The coordinates (x) were unshifted in update */
        if (bFFscan && (shellfc==NULL || bConverged))
        {
            if (print_forcefield(fplog,enerd->term,mdatoms->homenr,f,buf,xcopy,
                                 &(top_global->mols),mdatoms->massT,pres))
            {
                if (gmx_parallel_env)
                    gmx_finalize(cr);
                fprintf(stderr,"\n");
                exit(0);
            }
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
            do_coupling(fplog,nfile,fnm,tcr,t,step,enerd->term,fr,
                        ir,MASTER(cr),
                        mdatoms,&(top->idef),mu_aver,
                        top_global->mols.nr,cr,
                        state->box,total_vir,pres,
                        mu_tot,state->x,f,bConverged);
            debug_gmx();
        }
        
        /* Time for performance */
        if (((step % stepout) == 0) || bLastStep)
        {
            update_time();
        }
        
        /* Output stuff */
        if (MASTER(cr))
        {
            bool do_dr,do_or;
            
            upd_mdebin(mdebin,bDoDHDL ? fp_dhdl : NULL,bGStatEveryStep,
                       t,mdatoms->tmass,enerd,state,lastbox,
                       shake_vir,force_vir,total_vir,pres,
                       ekind,mu_tot,constr);
            
            do_dr  = do_per_step(step,ir->nstdisreout);
            do_or  = do_per_step(step,ir->nstorireout);
            
            print_ebin(fp_ene,do_ene,do_dr,do_or,do_log?fplog:NULL,step,t,
                       eprNORMAL,bCompact,mdebin,fcd,groups,&(ir->opts));
      
            if (ir->ePull != epullNO)
            {
                pull_print_output(ir->pull,step,t);
            }
            
            if (do_per_step(step,ir->nstlog))
            {
                if(fflush(fplog) != 0)
                {
                    gmx_fatal(FARGS,"Cannot flush logfile - maybe you are out of quota?");
                }
            }
        }
        
    
        /* Remaining runtime */
        if (MULTIMASTER(cr) && do_verbose)
        {
            if (shellfc)
            {
                fprintf(stderr,"\n");
            }
            print_time(stderr,start_t,step,ir);
        }
        
        /* Replica exchange */
        bExchanged = FALSE;
        if ((repl_ex_nst > 0) && (step > 0) && !bLastStep &&
            do_per_step(step,repl_ex_nst))
        {
            bExchanged = replica_exchange(fplog,cr,repl_ex,
                                          state_global,enerd->term,
                                          state,step,t);
        }
        if (bExchanged && PAR(cr))
        {
            if (DOMAINDECOMP(cr))
            {
                dd_partition_system(fplog,step,cr,TRUE,
                                    state_global,top_global,ir,
                                    state,&f,&buf,mdatoms,top,fr,
                                    vsite,shellfc,constr,
                                    nrnb,wcycle,FALSE);
            }
            else
            {
                bcast_state(cr,state,FALSE);
            }
        }
        
        bFirstStep = FALSE;
        
        if (bRerunMD) 
        {
            /* read next frame from input trajectory */
            bNotLastFrame = read_next_frame(status,&rerun_fr);
        }
        
        if (!bRerunMD || !rerun_fr.bStep)
        {
            /* increase the MD step number */
            step++;
            step_rel++;
        }
        
        cycles = wallcycle_stop(wcycle,ewcSTEP);
        if (DOMAINDECOMP(cr) && wcycle)
        {
            dd_cycles_add(cr->dd,cycles,ddCyclStep);
        }
    }
    /* End of main MD loop */
    debug_gmx();
    
    if (bRerunMD)
    {
        close_trj(status);
    }
    
    if (!(cr->duty & DUTY_PME))
    {
        /* Tell the PME only node to finish */
        gmx_pme_finish(cr);
    }
    
    if (MASTER(cr)) {
        if (bGStatEveryStep) {
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
        
    if (ir->nstlist == -1 && nns>0 && fplog)
    {
        fprintf(fplog,"Average neighborlist lifetime: %.1f steps, std.dev.: %.1f steps\n",ns_s1/nns,sqrt(ns_s2/nns - sqr(ns_s1/nns)));
        fprintf(fplog,"Average number of atoms that crossed the half buffer length: %.1f\n\n",ns_ab/nns);
    }
    
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
    
    *nsteps_done = step_rel;
    
    return start_t;
}
