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
#include "dummies.h"
#include "update.h"
#include "ns.h"
#include "trnio.h"
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

volatile bool bGotTermSignal = FALSE, bGotUsr1Signal = FALSE;

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



void mdrunner(t_commrec *cr,t_commrec *mcr,int nfile,t_filenm fnm[],
	      bool bVerbose,bool bCompact,
	      int nDlb,int nstepout,t_edsamyn *edyn,
	      unsigned long Flags)
{
  double     nodetime=0,realtime;
  t_parm     *parm;
  t_state    *state;
  rvec       *buf,*f,*vold,*vt,box_size;
  real       tmpr1,tmpr2;
  real       *ener;
  t_nrnb     *nrnb;
  t_nsborder *nsb;
  t_topology *top;
  t_groups   *grps;
  t_graph    *graph;
  t_mdatoms  *mdatoms;
  t_forcerec *fr;
  t_fcdata   *fcd;
  time_t     start_t=0;
  bool       bDummies,bParDummies;
  t_comm_dummies dummycomm;
  int        i,m;
  char       *gro;
  
  /* Initiate everything (snew sets to zero!) */
  snew(ener,F_NRE);
  snew(fcd,1);
  snew(nsb,1);
  snew(top,1);
  snew(grps,1);
  snew(parm,1);
  snew(state,1);
  snew(nrnb,cr->nnodes);
  
  if (bVerbose && MASTER(cr)) 
    fprintf(stderr,"Getting Loaded...\n");

  if (PAR(cr)) {
    /* The master thread on the master node reads from disk, then passes everything
     * around the ring, and finally frees the stuff
     */
    if (MASTER(cr)) 
      distribute_parts(cr->left,cr->right,cr->nodeid,cr->nnodes,parm,
		       ftp2fn(efTPX,nfile,fnm),nDlb);
    
    /* Every node (including the master) reads the data from the ring */
    init_parts(stdlog,cr,
	       parm,top,state,&mdatoms,nsb,
	       MASTER(cr) ? LIST_SCALARS | LIST_PARM : 0,
	       &bParDummies,&dummycomm);
  } else {
    /* Read it up... */
    init_single(stdlog,parm,ftp2fn(efTPX,nfile,fnm),top,state,&mdatoms,nsb);
    bParDummies=FALSE;
  }
  if (parm->ir.eI == eiSD) {
    /* Is not read from TPR yet, so we allocate space here */
    snew(state->sd_X,nsb->natoms);
  }
  snew(buf,nsb->natoms);
  snew(f,nsb->natoms);
  snew(vt,nsb->natoms);
  snew(vold,nsb->natoms);

  if (bVerbose && MASTER(cr))
    fprintf(stderr,"Loaded with Money\n\n");
  
  /* Index numbers for parallellism... */
  nsb->nodeid      = cr->nodeid;
  top->idef.nodeid = cr->nodeid;

  /* Group stuff (energies etc) */
  init_groups(stdlog,mdatoms,&(parm->ir.opts),grps);
  /* Copy the cos acceleration to the groups struct */
  grps->cosacc.cos_accel = parm->ir.cos_accel;
  
  /* Periodicity stuff */  
  if (parm->ir.ePBC == epbcXYZ) {
    graph=mk_graph(&(top->idef),top->atoms.nr,FALSE,FALSE);
    if (debug)
      p_graph(debug,"Initial graph",graph);
  }
  else
    graph = NULL;
    
  /* Distance Restraints */
  init_disres(stdlog,top->idef.il[F_DISRES].nr,top->idef.il[F_DISRES].iatoms,
	      top->idef.iparams,&(parm->ir),mcr,fcd);

  /* Orientation restraints */
  init_orires(stdlog,top->idef.il[F_ORIRES].nr,top->idef.il[F_ORIRES].iatoms,
	      top->idef.iparams,state->x,mdatoms,&(parm->ir),mcr,
	      &(fcd->orires));

  /* Dihedral Restraints */
  init_dihres(stdlog,top->idef.il[F_DIHRES].nr,top->idef.il[F_DIHRES].iatoms,
	      top->idef.iparams,&(parm->ir),fcd);

  /* check if there are dummies */
  bDummies=FALSE;
  for(i=0; (i<F_NRE) && !bDummies; i++)
    bDummies = ((interaction_function[i].flags & IF_DUMMY) && 
		(top->idef.il[i].nr > 0));

  /* Initiate forcerecord */
  fr = mk_forcerec();
  init_forcerec(stdlog,fr,&(parm->ir),top,cr,mdatoms,nsb,state->box,FALSE,
		opt2fn("-table",nfile,fnm),FALSE);
  fr->bSepDVDL = ((Flags & MD_SEPDVDL) == MD_SEPDVDL);
  /* Initiate box */
  for(m=0; (m<DIM); m++)
    box_size[m]=state->box[m][m];
    
  /* Initiate PPPM if necessary */
  if (fr->eeltype == eelPPPM)
    init_pppm(stdlog,cr,nsb,FALSE,TRUE,box_size,getenv("GMXGHAT"),&parm->ir);
  if (fr->eeltype == eelPME)
    (void) init_pme(stdlog,cr,parm->ir.nkx,parm->ir.nky,parm->ir.nkz,
		    parm->ir.pme_order,HOMENR(nsb),
		    parm->ir.bOptFFT,parm->ir.ewald_geometry);
  
  /* Now do whatever the user wants us to do (how flexible...) */
  switch (parm->ir.eI) {
  case eiMD:
  case eiSD:
  case eiBD:
    start_t=do_md(stdlog,cr,mcr,nfile,fnm,
		  bVerbose,bCompact,bDummies,
		  bParDummies ? &dummycomm : NULL,
		  nstepout,parm,grps,top,ener,fcd,state,vold,vt,f,buf,
		  mdatoms,nsb,nrnb,graph,edyn,fr,box_size,Flags);
    break;
  case eiCG:
    start_t=do_cg(stdlog,nfile,fnm,parm,top,grps,nsb,
		  state,f,buf,mdatoms,parm->ekin,ener,fcd,
		  nrnb,bVerbose,bDummies,
		  bParDummies ? &dummycomm : NULL,
		  cr,mcr,graph,fr,box_size);
    break;
  case eiLBFGS:
    start_t=do_lbfgs(stdlog,nfile,fnm,parm,top,grps,nsb,
		     state,f,buf,mdatoms,parm->ekin,ener,fcd,
		     nrnb,bVerbose,bDummies,
		     bParDummies ? &dummycomm : NULL,
		     cr,mcr,graph,fr,box_size);
    break;
  case eiSteep:
    start_t=do_steep(stdlog,nfile,fnm,parm,top,grps,nsb,
		     state,f,buf,mdatoms,parm->ekin,ener,fcd,
		     nrnb,bVerbose,bDummies,
		     bParDummies ? &dummycomm : NULL,
		     cr,mcr,graph,fr,box_size);
    break;
  case eiNM:
    start_t=do_nm(stdlog,cr,nfile,fnm,
		  bVerbose,bCompact,nstepout,parm,grps,
		  top,ener,fcd,state,vold,vt,f,buf,
		  mdatoms,nsb,nrnb,graph,edyn,fr,box_size);
    break;
  default:
    fatal_error(0,"Invalid integrator (%d)...\n",parm->ir.eI);
  }
  
    /* Some timing stats */  
  if (MASTER(cr)) {
    realtime=difftime(time(NULL),start_t);
    if ((nodetime=node_time()) == 0)
      nodetime=realtime;
  }
  else 
    realtime=0;
    
  /* Convert back the atoms */
  md2atoms(mdatoms,&(top->atoms),TRUE);
  
  /* Finish up, write some stuff
   * if rerunMD, don't write last frame again 
   */
  finish_run(stdlog,cr,ftp2fn(efSTO,nfile,fnm),
	     nsb,top,parm,nrnb,nodetime,realtime,parm->ir.nsteps,
	     parm->ir.eI==eiMD || parm->ir.eI==eiSD || parm->ir.eI==eiBD);
  
  /* Does what it says */  
  print_date_and_time(stdlog,cr->nodeid,"Finished mdrun");
}

time_t do_md(FILE *log,t_commrec *cr,t_commrec *mcr,int nfile,t_filenm fnm[],
	     bool bVerbose,bool bCompact,
	     bool bDummies, t_comm_dummies *dummycomm,
	     int stepout,t_parm *parm,t_groups *grps,t_topology *top,
	     real ener[],t_fcdata *fcd,
	     t_state *state,rvec vold[],rvec vt[],rvec f[],
	     rvec buf[],t_mdatoms *mdatoms,t_nsborder *nsb,t_nrnb nrnb[],
	     t_graph *graph,t_edsamyn *edyn,t_forcerec *fr,rvec box_size,
	     unsigned long Flags)
{
  t_mdebin   *mdebin;
  int        fp_ene=0,fp_trn=0,step,step_rel;
  FILE       *fp_dgdl=NULL,*fp_field=NULL;
  time_t     start_t;
  real       t,t0,lam0;
  bool       bNS,bSimAnn,bStopCM,bTYZ,bRerunMD,bNotLastFrame=FALSE,
             bFirstStep,bLastStep,bNEMD,do_log,bRerunWarnNoV=TRUE,
	     bFullPBC;
  tensor     force_vir,shake_vir;
  t_nrnb     mynrnb;
  char       *traj,*xtc_traj; /* normal and compressed trajectory filename */
  int        i,m,status;
  rvec       mu_tot;
  rvec       *xx,*vv,*ff;
  t_vcm      *vcm;
  t_trxframe rerun_fr;
  t_pull     pulldata; /* for pull code */
  /* A boolean (disguised as a real) to terminate mdrun */  
  real       terminate=0;
 
  /* XMDRUN stuff: shell, general coupling etc. */
  bool        bFFscan;
  int         nshell,nflexcon,nshell_flexcon_tot,count,nconverged=0;
  t_shell     *shells=NULL;
  real        timestep=0;
  double      tcount=0;
  bool        bShell_FlexCon,bIonize=FALSE,bGlas=FALSE;
  bool        bTCR=FALSE,bConverged=TRUE,bOK;
  real        mu_aver=0,fmax;
  int         gnx,ii;
  atom_id     *grpindex;
  char        *grpname;
  t_coupl_rec *tcr=NULL;
  rvec        *xcopy=NULL,*vcopy=NULL;
  matrix      boxcopy,lastbox;
  /* End of XMDRUN stuff */

  /* Turn on signal handling */
  signal(SIGTERM,signal_handler);
  signal(SIGUSR1,signal_handler);

  /* Check for special mdrun options */
  bRerunMD = (Flags & MD_RERUN)  == MD_RERUN;
  bIonize  = (Flags & MD_IONIZE) == MD_IONIZE;
  bGlas    = (Flags & MD_GLAS)   == MD_GLAS;
  bFFscan  = (Flags & MD_FFSCAN) == MD_FFSCAN;
  
  /* Initial values */
  init_md(cr,&parm->ir,state->box,&t,&t0,&state->lambda,&lam0,
	  &mynrnb,&bTYZ,top,
	  nfile,fnm,&traj,&xtc_traj,&fp_ene,&fp_dgdl,&fp_field,&mdebin,grps,
	  force_vir,shake_vir,mdatoms,mu_tot,&bNEMD,&bSimAnn,&vcm,nsb);
  debug_gmx();

  /* Check for full periodicity calculations */
  bFullPBC = (parm->ir.ePBC == epbcFULL);  
  
  /* Check for polarizable models */
  shells     = init_shells(log,START(nsb),HOMENR(nsb),&top->idef,
			   mdatoms,&nshell);

  /* Check for flexible constraints */
  nflexcon = count_flexible_constraints(log,fr,&top->idef);

  nshell_flexcon_tot = nshell + nflexcon;
  if (PAR(cr))
    gmx_sumi(1,&nshell_flexcon_tot,cr);
  bShell_FlexCon = nshell_flexcon_tot > 0;
  
  gnx = top->blocks[ebMOLS].nr;
  snew(grpindex,gnx);
  for(i=0; (i<gnx); i++)
    grpindex[i] = i;

  /* Check whether we have to GCT stuff */
  bTCR = ftp2bSet(efGCT,nfile,fnm);
  if (MASTER(cr) && bTCR)
    fprintf(stderr,"Will do General Coupling Theory!\n");

  /* Remove periodicity */  
  if (fr->ePBC != epbcNONE)
    do_pbc_first(log,state->box,box_size,fr,graph,state->x);
  debug_gmx();

  /* Initialize pull code */
  init_pull(log,nfile,fnm,&pulldata,state->x,mdatoms,parm->ir.opts.nFreeze,
	    state->box,START(nsb),HOMENR(nsb),cr);
    
  if (!parm->ir.bUncStart && !bRerunMD) 
    do_shakefirst(log,bTYZ,ener,parm,nsb,mdatoms,state,vold,buf,f,
		  graph,cr,&mynrnb,grps,fr,top,edyn,&pulldata);
  debug_gmx();
    
  /* Compute initial EKin for all.. */
  if (grps->cosacc.cos_accel == 0)
    calc_ke_part(TRUE,parm->ir.eI==eiSD,
		 START(nsb),HOMENR(nsb),vold,state->v,vt,&(parm->ir.opts),
		 mdatoms,grps,&mynrnb,state->lambda,&ener[F_DVDLKIN]);
  else
    calc_ke_part_visc(TRUE,START(nsb),HOMENR(nsb),
		      state->box,state->x,vold,state->v,vt,&(parm->ir.opts),
		      mdatoms,grps,&mynrnb,state->lambda,&ener[F_DVDLKIN]);
  debug_gmx();
	
  if (PAR(cr)) 
    global_stat(log,cr,ener,force_vir,shake_vir,
		&(parm->ir.opts),grps,&mynrnb,nrnb,vcm,&terminate);
  debug_gmx();
  
  /* Calculate Temperature coupling parameters lambda */
  ener[F_TEMP] = sum_ekin(&(parm->ir.opts),grps,parm->ekin,bTYZ);
  /*
  if(parm->ir.etc==etcBERENDSEN)
    berendsen_tcoupl(&(parm->ir.opts),grps,
		     parm->ir.delta_t);
  else if(parm->ir.etc==etcNOSEHOOVER)
    nosehoover_tcoupl(&(parm->ir.opts),grps,
		      parm->ir.delta_t);
  */
  debug_gmx();
  
  /* Initiate data for the special cases */
  if (bFFscan) {
    snew(xcopy,nsb->natoms);
    snew(vcopy,nsb->natoms);
    for(ii=0; (ii<nsb->natoms); ii++) {
      copy_rvec(state->x[ii],xcopy[ii]);
      copy_rvec(state->v[ii],vcopy[ii]);
    }
    copy_mat(state->box,boxcopy);
  } 
  /* Write start time and temperature */
  start_t=print_date_and_time(log,cr->nodeid,"Started mdrun");
  
  if (MASTER(cr)) {
    fprintf(log,"Initial temperature: %g K\n",ener[F_TEMP]);
    if (bRerunMD) {
      fprintf(stderr,"starting md rerun '%s', reading coordinates from"
	      " input trajectory '%s'\n\n",
	      *(top->name),opt2fn("-rerun",nfile,fnm));
      if (bVerbose)
	fprintf(stderr,"Calculated time to finish depends on nsteps from "
		"run input file,\nwhich may not correspond to the time "
		"needed to process input trajectory.\n\n");
    } else
      fprintf(stderr,"starting mdrun '%s'\n%d steps, %8.1f ps.\n\n",
	      *(top->name),parm->ir.nsteps,parm->ir.nsteps*parm->ir.delta_t);
  }

  /* Initialize values for invmass, etc. */
  update_mdatoms(mdatoms,state->lambda,TRUE);

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
    bNotLastFrame = read_first_frame(&status,opt2fn("-rerun",nfile,fnm),
				     &rerun_fr,TRX_NEED_X | TRX_READ_V);
    if (rerun_fr.natoms != mdatoms->nr)
      fatal_error(0,"Number of atoms in trajectory (%d) does not match the "
		  "run input file (%d)\n",rerun_fr.natoms,mdatoms->nr);
  } 
  
  /* loop over MD steps or if rerunMD to end of input trajectory */
  bFirstStep = TRUE;
  bLastStep = FALSE;
  step = parm->ir.init_step;
  step_rel = 0;
  while ((!bRerunMD && (step_rel <= parm->ir.nsteps)) ||  
	 (bRerunMD && bNotLastFrame)) {
    
    if (bRerunMD) {
      if (rerun_fr.bStep) {
	step = rerun_fr.step;
	step_rel = step - parm->ir.init_step;
      }
      if (rerun_fr.bTime)
	t = rerun_fr.time;
      else
	t = step;
    } else {
      bLastStep = (step_rel == parm->ir.nsteps);

      t = t0 + step*parm->ir.delta_t;
    }
    
    do_log = do_per_step(step,parm->ir.nstlog) || bLastStep;

    if (parm->ir.efep != efepNO) {
      if (bRerunMD && rerun_fr.bLambda && (parm->ir.delta_lambda!=0))
	state->lambda = rerun_fr.lambda;
      else
	state->lambda = lam0 + step*parm->ir.delta_lambda;
    }
    if (MASTER(cr) && do_log && !bFFscan)
      print_ebin_header(log,step,t,state->lambda);

    if (bSimAnn) 
      update_annealing_target_temp(&(parm->ir.opts),t);
    
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
      
      /* for rerun MD always do Neighbour Searching */
      bNS = ((parm->ir.nstlist!=0) || bFirstStep);
    } else {
      /* Determine whether or not to do Neighbour Searching */
      bNS = ((parm->ir.nstlist && (step % parm->ir.nstlist==0)) || bFirstStep);
    }
    
    /* Stop Center of Mass motion */
    bStopCM = do_per_step(step,abs(parm->ir.nstcomm));

    /* Copy back starting coordinates in case we're doing a forcefield scan */
    if (bFFscan) {
      for(ii=0; (ii<nsb->natoms); ii++) {
	copy_rvec(xcopy[ii],state->x[ii]);
	copy_rvec(vcopy[ii],state->v[ii]);
      }
      copy_mat(boxcopy,state->box);
    }
    
    if (bDummies) {
      if (graph) {
	/* Following is necessary because the graph may get out of sync
	 * with the coordinates if we only have every N'th coordinate set
	 */
	if (bRerunMD)
	  mk_mshift(log,graph,state->box,state->x);
	shift_self(graph,state->box,state->x);
      }
      construct_dummies(log,state->x,&mynrnb,parm->ir.delta_t,state->v,
			&top->idef,graph,cr,state->box,dummycomm);
      
      if (graph)
	unshift_self(graph,state->box,state->x);
    }
     
    debug_gmx();
    
    /* Set values for invmass etc. This routine not parallellized, but hardly
     * ever used, only when doing free energy calculations.
     */
    if(parm->ir.efep != efepNO)
      update_mdatoms(mdatoms,state->lambda,FALSE); 
    
    clear_mat(force_vir);
    
    /* Ionize the atoms if necessary */
    if (bIonize)
      ionize(log,mdatoms,top->atoms.atomname,t,&parm->ir,state->x,state->v,
	     START(nsb),START(nsb)+HOMENR(nsb),state->box,cr);
      
    /* Update force field in ffscan program */
    if (bFFscan) 
      update_forcefield(nfile,fnm,fr,mdatoms->nr,state->x,state->box);
    
    if (bShell_FlexCon) {
      /* Now is the time to relax the shells */
      count=relax_shells(log,cr,mcr,bVerbose,bFFscan ? step+1 : step,
			 parm,bNS,bStopCM,top,ener,fcd,
			 state,vold,vt,f,buf,mdatoms,nsb,&mynrnb,graph,
			 grps,force_vir,
			 nshell,shells,nflexcon,fr,traj,t,mu_tot,
			 nsb->natoms,&bConverged,bDummies,dummycomm,
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
      do_force(log,cr,mcr,parm,nsb,force_vir,step,&mynrnb,top,grps,
	       state->box,state->x,f,buf,mdatoms,ener,fcd,bVerbose && !PAR(cr),
	       state->lambda,graph,bNS,FALSE,fr,mu_tot,FALSE,t,fp_field);
    }
   
    if (bTCR)
      mu_aver=calc_mu_aver(cr,nsb,state->x,mdatoms->chargeA,mu_tot,top,mdatoms,
			   gnx,grpindex);
    
    if (bGlas)
      do_glas(log,START(nsb),HOMENR(nsb),state->x,f,
	      fr,mdatoms,top->idef.atnr,&parm->ir,ener);
    
    if (bTCR && bFirstStep) {
      tcr=init_coupling(log,nfile,fnm,cr,fr,mdatoms,&(top->idef));
      fprintf(log,"Done init_coupling\n"); 
      fflush(log);
    }
    
    /* Now we have the energies and forces corresponding to the 
     * coordinates at time t. We must output all of this before
     * the update.
     * for RerunMD t is read from input trajectory
     */
    if (bDummies) 
      spread_dummy_f(log,state->x,f,&mynrnb,&top->idef,dummycomm,cr);
      
    /* Calculation of the virial must be done after dummies!    */
    /* Question: Is it correct to do the PME forces after this? */
    calc_virial(log,START(nsb),HOMENR(nsb),state->x,f,
		force_vir,fr->vir_el_recip,graph,state->box,&mynrnb,fr);
		  
    /* Spread the LR force on dummy particle to the other particles... 
     * This is parallellized. MPI communication is performed
     * if the constructing atoms aren't local.
     */
    if (bDummies && fr->bEwald) 
      spread_dummy_f(log,state->x,fr->f_el_recip,&mynrnb,&top->idef,
		     dummycomm,cr);
    
    sum_lrforces(f,fr,START(nsb),HOMENR(nsb));

    xx = (do_per_step(step,parm->ir.nstxout) || bLastStep) ? state->x : NULL;
    vv = (do_per_step(step,parm->ir.nstvout) || bLastStep) ? state->v : NULL;
    ff = (do_per_step(step,parm->ir.nstfout)) ? f : NULL;

    fp_trn = write_traj(log,cr,traj,nsb,step,t,state->lambda,
			nrnb,nsb->natoms,xx,vv,ff,state->box);
    debug_gmx();
    
    /* don't write xtc and last structure for rerunMD */
    if (!bRerunMD && !bFFscan) {
      if (do_per_step(step,parm->ir.nstxtcout))
	write_xtc_traj(log,cr,xtc_traj,nsb,mdatoms,
		       step,t,state->x,state->box,parm->ir.xtcprec);
      if (bLastStep && MASTER(cr)) {
	fprintf(stderr,"Writing final coordinates.\n");
	write_sto_conf(ftp2fn(efSTO,nfile,fnm),
		       *top->name, &(top->atoms),state->x,state->v,state->box);
      }
      debug_gmx();
    }
    clear_mat(shake_vir);
    
    /* Afm and Umbrella type pulling happens before the update, 
     * other types in update 
     */
    if (pulldata.bPull && 
	(pulldata.runtype == eAfm || pulldata.runtype == eUmbrella))
      pull(&pulldata,state->x,f,state->box,top,parm->ir.delta_t,step,t,
	   mdatoms,START(nsb),HOMENR(nsb),cr); 
    
    if (bFFscan)
      clear_rvecs(nsb->natoms,buf);

    /* Box is changed in update() when we do pressure coupling,
     * but we should still use the old box for energy corrections and when
     * writing it to the energy file, so it matches the trajectory files for
     * the same timestep above. Make a copy in a separate array.
     */
    copy_mat(state->box,lastbox);
    
    /* This is also parallellized, but check code in update.c */
    /* bOK = update(nsb->natoms,START(nsb),HOMENR(nsb),step,state->lambda,&ener[F_DVDL], */
    bOK = TRUE;
    if (!rerun_fr.bV)
      update(nsb->natoms,START(nsb),HOMENR(nsb),step,&ener[F_DVDL],
	     parm,mdatoms,state,graph,f,buf,vold,
	     top,grps,shake_vir,cr,&mynrnb,bTYZ,edyn,&pulldata,bNEMD,
	     TRUE,bFirstStep,NULL);
    else {
      /* Need to unshift here */
      if ((parm->ir.ePBC == epbcXYZ) && (graph->nnodes > 0))
	unshift_self(graph,state->box,state->x);
    }
    if (!bOK && !bFFscan)
      fatal_error(0,"Constraint error: Shake, Lincs or Settle could not solve the constrains");

    /* Correct the new box if it is too skewed */
    if ((parm->ir.epc != epcNO) && !bRerunMD)
      correct_box(state->box,fr,graph);
    /* (un)shifting should NOT be done after this,
     * since the box vectors might have changed
     */

    /* Non-equilibrium MD: this is parallellized, but only does communication
     * when there really is NEMD.
     */
    if (PAR(cr) && bNEMD) 
      accumulate_u(cr,&(parm->ir.opts),grps);
      
    debug_gmx();
    if (grps->cosacc.cos_accel == 0)
      calc_ke_part(FALSE,parm->ir.eI==eiSD,
		   START(nsb),HOMENR(nsb),vold,state->v,vt,&(parm->ir.opts),
		   mdatoms,grps,&mynrnb,state->lambda,&ener[F_DVDLKIN]);
    else
      calc_ke_part_visc(FALSE,START(nsb),HOMENR(nsb),
			state->box,state->x,vold,state->v,vt,&(parm->ir.opts),
			mdatoms,grps,&mynrnb,state->lambda,&ener[F_DVDLKIN]);

    /* since we use the new coordinates in calc_ke_part_visc, we should use
     * the new box too. Still, won't this be offset by one timestep in the
     * energy file? / EL 20040121
     */ 
    
    debug_gmx();
    /* Calculate center of mass velocity if necessary, also parallellized */
    if (bStopCM && !bFFscan && !bRerunMD)
      calc_vcm_grp(log,START(nsb),HOMENR(nsb),mdatoms->massT,
		   state->x,state->v,vcm);

    /* Check whether everything is still allright */    
    if (bGotTermSignal || bGotUsr1Signal) {
      if (bGotTermSignal)
	terminate = 1;
      else
	terminate = -1;
      fprintf(log,"\n\nReceived the %s signal\n\n",
	      bGotTermSignal ? "TERM" : "USR1");
      fflush(log);
      if (MASTER(cr)) {
	fprintf(stderr,"\n\nReceived the %s signal\n\n",
	      bGotTermSignal ? "TERM" : "USR1");
	fflush(stderr);
      }
      bGotTermSignal = FALSE;
      bGotUsr1Signal = FALSE;
    }

    if (PAR(cr)) {
      /* Globally (over all NODEs) sum energy, virial etc. 
       * This includes communication 
       */
      global_stat(log,cr,ener,force_vir,shake_vir,
		  &(parm->ir.opts),grps,&mynrnb,nrnb,vcm,&terminate);
      /* Correct for double counting energies, should be moved to 
       * global_stat 
       */
      if (fr->bTwinRange && !bNS) 
	for(i=0; (i<grps->estat.nn); i++) {
	  grps->estat.ee[egCOULLR][i] /= cr->nnodes;
	  grps->estat.ee[egLJLR][i]   /= cr->nnodes;
	}
    }
    else
      cp_nrnb(&(nrnb[0]),&mynrnb);
    
    /* This is just for testing. Nothing is actually done to Ekin
     * since that would require extra communication.
     */
    if (!bNEMD && debug && (vcm->nr > 0))
      correct_ekin(debug,START(nsb),START(nsb)+HOMENR(nsb),state->v,
		   vcm->group_p[0],
		   mdatoms->massT,mdatoms->tmass,parm->ekin);
    
    if ((terminate != 0) && (step - parm->ir.init_step < parm->ir.nsteps)) {
      if (terminate<0 && parm->ir.nstxout)
	/* this is the USR1 signal and we are writing x to trr, 
	   stop at next x frame in trr */
	parm->ir.nsteps =
	  (step/parm->ir.nstxout + 1) * parm->ir.nstxout - parm->ir.init_step;
      else
	parm->ir.nsteps = step+1;
      fprintf(log,"\nSetting nsteps to %d\n\n",parm->ir.nsteps);
      fflush(log);
      if (MASTER(cr)) {
	fprintf(stderr,"\nSetting nsteps to %d\n\n",parm->ir.nsteps);
	fflush(stderr);
      }
      /* erase the terminate signal */
      terminate = 0;
    }
      
    /* Do center of mass motion removal */
    if (bStopCM && !bFFscan && !bRerunMD) {
      check_cm_grp(log,vcm);
      do_stopcm_grp(log,START(nsb),HOMENR(nsb),state->x,state->v,vcm);
      inc_nrnb(&mynrnb,eNR_STOPCM,HOMENR(nsb));
      /*
      calc_vcm_grp(log,START(nsb),HOMENR(nsb),mdatoms->massT,x,v,vcm);
      check_cm_grp(log,vcm);
      do_stopcm_grp(log,START(nsb),HOMENR(nsb),x,v,vcm);
      check_cm_grp(log,vcm);
      */
    }
        
    /* Add force and shake contribution to the virial */
    m_add(force_vir,shake_vir,parm->vir);
  
    /* Sum the potential energy terms from group contributions */
    sum_epot(&(parm->ir.opts),grps,ener);

    /* Calculate the amplitude of the cosine velocity profile */
    grps->cosacc.vcos = grps->cosacc.mvcos/mdatoms->tmass;

    /* Sum the kinetic energies of the groups & calc temp */
    ener[F_TEMP]=sum_ekin(&(parm->ir.opts),grps,parm->ekin,bTYZ);
    ener[F_EKIN]=trace(parm->ekin);
    ener[F_ETOT]=ener[F_EPOT]+ener[F_EKIN];
    
    /* Check for excessively large energies */
    if (bIonize && fabs(ener[F_ETOT]) > 1e18) {
      fprintf(stderr,"Energy too large (%g), giving up\n",ener[F_ETOT]);
      break;
    }
      
    /* Calculate Temperature coupling parameters lambda and adjust
     * target temp when doing simulated annealing
     */
    /*
    if(parm->ir.etc==etcBERENDSEN)
      berendsen_tcoupl(&(parm->ir.opts),grps,parm->ir.delta_t);
    else if(parm->ir.etc==etcNOSEHOOVER)
      nosehoover_tcoupl(&(parm->ir.opts),grps,parm->ir.delta_t);
    */

    /* Calculate pressure and apply LR correction if PPPM is used.
     * Use the box from last timestep since we already called update().
     */
    calc_pres(fr->ePBC,lastbox,parm->ekin,parm->vir,parm->pres,
	      (fr->eeltype==eelPPPM) ? ener[F_COUL_RECIP] : 0.0);
    
    /* Calculate long range corrections to pressure and energy */
    if (bTCR || bFFscan)
      set_avcsixtwelve(log,fr,mdatoms,&top->atoms.excl);
    
    /* Calculate long range corrections to pressure and energy */
    calc_dispcorr(log,parm->ir.eDispCorr,
		  fr,mdatoms->nr,lastbox,parm->pres,parm->vir,ener);

    /* The coordinates (x) were unshifted in update */
    if (bFFscan && (!bShell_FlexCon || bConverged))
      print_forcefield(log,ener,HOMENR(nsb),f,buf,xcopy,
		       &(top->blocks[ebMOLS]),mdatoms->massT,
		       parm->pres); 
    
    if (bTCR) {
      /* Only do GCT when the relaxation of shells (minimization) has converged,
       * otherwise we might be coupling to bogus energies. 
       * In parallel we must always do this, because the other sims might
       * update the FF.
       */
      
      /* Since this is called with the new coordinates state->x, I assume
       * we want the new box state->box too. / EL 20040121
       */
      do_coupling(log,nfile,fnm,tcr,t,step,ener,fr,
		  &(parm->ir),MASTER(cr),
		  mdatoms,&(top->idef),mu_aver,
		  top->blocks[ebMOLS].nr,(mcr!=NULL) ? mcr : cr,
		  state->box,parm->vir,parm->pres,
		  mu_tot,state->x,f,bConverged);
      debug_gmx();
    }

    /* Time for performance */
    if (((step % stepout) == 0) || bLastStep)
      update_time();

    /* Output stuff */
    if (MASTER(cr)) {
      bool do_ene,do_dr,do_or,do_dihr;
      
      upd_mdebin(mdebin,fp_dgdl,mdatoms->tmass,step_rel,t,ener,state,lastbox,
		 shake_vir,force_vir,parm->vir,parm->pres,grps,mu_tot);
      do_ene = do_per_step(step,parm->ir.nstenergy) || bLastStep;
      do_dr  = do_per_step(step,parm->ir.nstdisreout) || bLastStep;
      do_or  = do_per_step(step,parm->ir.nstorireout) || bLastStep;
      do_dihr= do_per_step(step,parm->ir.nstdihreout) || bLastStep;
      print_ebin(fp_ene,do_ene,do_dr,do_or,do_dihr,do_log?log:NULL,step,t,
		 eprNORMAL,bCompact,mdebin,fcd,&(top->atoms),&(parm->ir.opts));
      if (bVerbose)
	fflush(log);
    }
    
    /* Remaining runtime */
    if (MASTER(cr) && bVerbose && ( ((step % stepout)==0) || bLastStep)) {
      if (bShell_FlexCon)
	fprintf(stderr,"\n");
      print_time(stderr,start_t,step,&parm->ir);
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
  }
  /* End of main MD loop */
  debug_gmx();

  /* Dump the NODE time to the log file on each node */
  if (PAR(cr)) {
    double *ct,ctmax,ctsum;
    
    snew(ct,cr->nnodes);
    ct[cr->nodeid] = node_time();
    gmx_sumd(cr->nnodes,ct,cr);
    ctmax = ct[0];
    ctsum = ct[0];
    for(i=1; (i<cr->nodeid); i++) {
      ctmax = max(ctmax,ct[i]);
      ctsum += ct[i];
    }
    ctsum /= cr->nnodes;
    fprintf(log,"\nTotal NODE time on node %d: %g\n",cr->nodeid,ct[cr->nodeid]);
    fprintf(log,"Average NODE time: %g\n",ctsum);
    fprintf(log,"Load imbalance reduced performance to %3d%% of max\n",
	    (int) (100.0*ctmax/ctsum));
    sfree(ct);
  }
  if (MASTER(cr)) {
    print_ebin(fp_ene,FALSE,FALSE,FALSE,FALSE,log,step,t,
	       eprAVER,FALSE,mdebin,fcd,&(top->atoms),&(parm->ir.opts));
    print_ebin(fp_ene,FALSE,FALSE,FALSE,FALSE,log,step,t,
	       eprRMS,FALSE,mdebin,fcd,&(top->atoms),&(parm->ir.opts));
    close_enx(fp_ene);
    if (!bRerunMD && parm->ir.nstxtcout)
      close_xtc_traj();
    close_trn(fp_trn);
    if (fp_dgdl)
      fclose(fp_dgdl);
    if (fp_field)
      fclose(fp_field);
  }
  debug_gmx();

  if (bShell_FlexCon) {
    fprintf(log,"Fraction of iterations that converged:           %.2f %%\n",
	    (nconverged*100.0)/(parm->ir.nsteps+1));
    fprintf(log,"Average number of force evaluations per MD step: %.2f\n",
	    tcount/(parm->ir.nsteps+1));
  }

  return start_t;
}
