/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Gyas ROwers Mature At Cryogenic Speed
 */
static char *SRCID_md_c = "$Id$";
#include <signal.h>
#include <stdlib.h>
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "statutil.h"
#include "vcm.h"
#include "mdebin.h"
#include "nrnb.h"
#include "calcmu.h"
#include "dummies.h"
#include "update.h"
#include "trnio.h"
#include "mdrun.h"
#include "confio.h"
#include "network.h"
#include "sim_util.h"
#include "pull.h"
#include "physics.h"
#include "names.h"

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

time_t do_md(FILE *log,t_commrec *cr,int nfile,t_filenm fnm[],
	     bool bVerbose,bool bCompact,bool bDummies, t_comm_dummies *dummycomm,
	     int stepout,t_parm *parm,t_groups *grps,t_topology *top,real ener[],
	     rvec x[],rvec vold[],rvec v[],rvec vt[],rvec f[],
	     rvec buf[],t_mdatoms *mdatoms,t_nsborder *nsb,t_nrnb nrnb[],
	     t_graph *graph,t_edsamyn *edyn,t_forcerec *fr,rvec box_size,
	     unsigned long Flags)
{
  t_mdebin   *mdebin;
  int        fp_ene=0,fp_trn=0,step;
  FILE       *fp_dgdl=NULL;
  time_t     start_t;
  real       t,lambda,t0,lam0,SAfactor;
  bool       bNS,bStopCM,bStopRot,bTYZ,bRerunMD,bNotLastFrame=FALSE,
             bFirstStep,bLastStep,bNEMD,do_log,bRerunWarnNoV=TRUE,
	     bLateVir,bTweak;
  tensor     force_vir,pme_vir,shake_vir;
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
 
#ifdef XMDRUN
  /* Shell stuff */
  int         nshell,count,nconverged=0;
  t_shell     *shells;
  real        timestep;
  double      tcount=0;
  bool        bDynamicStep,bIonize,bMultiSim,bGlas;
  bool        bTCR,bConverged;
  real        mu_aver=0,fmax;
  int         gnx;
  atom_id     *grpindex;
  char        *grpname;
  t_coupl_rec *tcr=NULL;
#endif

  /* Turn on signal handling */
  signal(SIGTERM,signal_handler);
  signal(SIGUSR1,signal_handler);

  /* check if "-rerun" is used. */
  bRerunMD = ((Flags & MD_RERUN) == MD_RERUN);
  /* Initial values */
  init_md(cr,&parm->ir,parm->box,&t,&t0,&lambda,&lam0,&SAfactor,
	  &mynrnb,&bTYZ,top,
	  nfile,fnm,&traj,&xtc_traj,&fp_ene,&fp_dgdl,&mdebin,grps,
	  force_vir,pme_vir,shake_vir,mdatoms,mu_tot,&bNEMD,&vcm,nsb);
  debug_gmx();

  bLateVir = ((Flags & MD_LATEVIR) == MD_LATEVIR);
  bTweak   = ((Flags & MD_TWEAK) == MD_TWEAK);
  
  fprintf(log,"Late virial: %s Tweak PME virial: %s\n",
	  bool_names[bLateVir],bool_names[bTweak]);
  
#ifdef XMDRUN
  bDynamicStep = FALSE;
  bIonize      = (Flags & MD_IONIZE)   == MD_IONIZE;
  bMultiSim    = (Flags & MD_MULTISIM) == MD_MULTISIM;
  bGlas        = (Flags & MD_GLAS)     == MD_GLAS;
  timestep = parm->ir.delta_t;

  /* Check whether we have to do dipole stuff */
  if (ftp2bSet(efNDX,nfile,fnm))
    rd_index(ftp2fn(efNDX,nfile,fnm),1,&gnx,&grpindex,&grpname);
  else {
    gnx = top->blocks[ebMOLS].nr;
    snew(grpindex,gnx);
    for(i=0; (i<gnx); i++)
      grpindex[i] = i;
  }
  /* Check whether we have to GCT stuff */
  bTCR = ftp2bSet(efGCT,nfile,fnm);
  if (MASTER(cr) && bTCR)
    fprintf(stderr,"Will do General Coupling Theory!\n");

  fr->k_dirmin = parm->ir.userreal4;
#endif

  /* Remove periodicity */  
  if (fr->ePBC != epbcNONE)
    do_pbc_first(log,parm,box_size,fr,graph,x);
  debug_gmx();

  /* Initialize pull code */
  init_pull(log,nfile,fnm,&pulldata,x,mdatoms,parm->box);
  if (pulldata.bPull && cr->nnodes>1)
    fatal_error(0,"Can not pull in parallel");
    
  if (!parm->ir.bUncStart) 
    do_shakefirst(log,bTYZ,lambda,ener,parm,nsb,mdatoms,x,vold,buf,f,v,
		  graph,cr,&mynrnb,grps,fr,top,edyn,&pulldata);
  debug_gmx();
    
  /* Compute initial EKin for all.. */
  if (grps->cosacc.cos_accel == 0)
    calc_ke_part(TRUE,parm->ir.eI==eiSD,
		 START(nsb),HOMENR(nsb),vold,v,vt,&(parm->ir.opts),
		 mdatoms,grps,&mynrnb,lambda,&ener[F_DVDLKIN]);
  else
    calc_ke_part_visc(TRUE,START(nsb),HOMENR(nsb),
		      parm->box,x,vold,v,vt,&(parm->ir.opts),
		      mdatoms,grps,&mynrnb,lambda,&ener[F_DVDLKIN]);
  debug_gmx();
	
  if (PAR(cr)) 
    global_stat(log,cr,ener,force_vir,shake_vir,
		&(parm->ir.opts),grps,&mynrnb,nrnb,vcm,&terminate);
  debug_gmx();
  
  /* Calculate Temperature coupling parameters lambda */
  ener[F_TEMP] = sum_ekin(&(parm->ir.opts),grps,parm->ekin,bTYZ);
  if(parm->ir.etc==etcBERENDSEN)
    berendsen_tcoupl(&(parm->ir.opts),grps,
		     parm->ir.delta_t,SAfactor);
  else if(parm->ir.etc==etcNOSEHOOVER)
    nosehoover_tcoupl(&(parm->ir.opts),grps,
		      parm->ir.delta_t,SAfactor);
  debug_gmx();
  
#ifdef XMDRUN
  shells = init_shells(log,START(nsb),HOMENR(nsb),&top->idef,mdatoms,&nshell);
#endif
  
  /* Write start time and temperature */
  start_t=print_date_and_time(log,cr->nodeid,
#ifdef XMDRUN
			      "Started Xmdrun"
#else
			      "Started mdrun"
#endif		      
			      );
  
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
  step = 0;
  while ((!bRerunMD && (step<=parm->ir.nsteps)) ||  
	 (bRerunMD && bNotLastFrame)) {
    
    /*dynamic_load_balancing(bVerbose,cr,&(top->blocks[ebCGS]),
			   &(top->blocks[ebSBLOCKS]),NULL,x,parm->box,XX);
    */
    if (bRerunMD) {
      if (rerun_fr.bStep)
	step = rerun_fr.step;
      if (rerun_fr.bTime)
	t = rerun_fr.time;
      else
	t = step;
      for(i=0; i<mdatoms->nr; i++)
	copy_rvec(rerun_fr.x[i],x[i]);
      if (rerun_fr.bV)
	for(i=0; i<mdatoms->nr; i++)
	  copy_rvec(rerun_fr.v[i],v[i]);
      else {
	for(i=0; i<mdatoms->nr; i++)
	    clear_rvec(v[i]);
	if (bRerunWarnNoV) {
	  fprintf(stderr,"\nWARNING: Some frames do not contain velocities.\n"
		  "         Ekin, temperature and pressure are incorrect,\n"
		  "         the virial will be incorrect when constraints are present.\n"
		  "\n");
	  bRerunWarnNoV = FALSE;
	}
      }
      copy_mat(rerun_fr.box,parm->box);
      
      /* for rerun MD always do Neighbour Searching */
      bNS = ((parm->ir.nstlist!=0) || bFirstStep);
    } else
      /* Determine whether or not to do Neighbour Searching */
      bNS = ((parm->ir.nstlist && (step % parm->ir.nstlist==0)) || bFirstStep);
    
    bLastStep=(step==parm->ir.nsteps);
    
    do_log = do_per_step(step,parm->ir.nstlog) || bLastStep;
    
    /* Stop Center of Mass motion */
    get_cmparm(&parm->ir,step,&bStopCM,&bStopRot);

    if (bDummies) {
      /* Construct dummy particles. This is parallellized.
       * If the atoms constructing a dummy are not present
       * on the local processor we start by fetching them,
       * and afterwards distribute possible nonlocal 
       * constructed dummies. (This is checked for during
       * setup to avoid unnecessary communication).
       */
      if(dummycomm) 
	move_construct_x(dummycomm,x,cr);

      shift_self(graph,parm->box,x);
    
      construct_dummies(log,x,&mynrnb,parm->ir.delta_t,v,&top->idef);

      unshift_self(graph,parm->box,x);
      
      if(dummycomm)
	move_dummy_xv(dummycomm,x,v,cr);
    }
     
    debug_gmx();
    
    /* Set values for invmass etc. This routine not parallellized, but hardly
     * ever used, only when doing free energy calculations.
     */
    init_mdatoms(mdatoms,lambda,bFirstStep);
    
    clear_mat(force_vir);
    
#ifdef XMDRUN
    /* Ionize the atoms if necessary */    
    if (bIonize)
      ionize(log,mdatoms,top->atoms.atomname,t,&parm->ir,x,v,
	     START(nsb),START(nsb)+HOMENR(nsb),parm->box,cr);
    
    /* Now is the time to relax the shells */
    count=relax_shells(log,cr,bVerbose,step,parm,bNS,bStopCM,top,ener,
		       x,vold,v,vt,f,
		       buf,mdatoms,nsb,&mynrnb,graph,grps,force_vir,
		       pme_vir,nshell,shells,fr,traj,t,lambda,mu_tot,
		       nsb->natoms,parm->box,&bConverged);
    tcount+=count;
    
    if (bConverged)
      nconverged++;

    if (bIonize)
      ener[F_SR] += 
	electron_atom_interactions(log,mdatoms,&parm->ir,
				   START(nsb),START(nsb)+HOMENR(nsb),
				   x,v,f,parm->box);
    
#else

    /* The coordinates (x) are shifted (to get whole molecules) in do_force
     * This is parallellized as well, and does communication too. 
     * Check comments in sim_util.c
     */
    do_force(log,cr,parm,nsb,force_vir,pme_vir,step,&mynrnb,
	     top,grps,x,v,f,buf,mdatoms,ener,bVerbose && !PAR(cr),
	     lambda,graph,bNS,FALSE,fr,mu_tot,FALSE);
 
#endif
    /* HACK */
    if (bTweak)
      sum_lrforces(f,fr,START(nsb),HOMENR(nsb));
    /* HACK */
    if (!bLateVir)
      calc_virial(log,START(nsb),HOMENR(nsb),x,f,
		  force_vir,pme_vir,cr,graph,parm->box,&mynrnb,fr,bTweak);
   
#ifdef XMDRUN
    if (bTCR)
      mu_aver=calc_mu_aver(cr,nsb,x,mdatoms->chargeA,mu_tot,top,mdatoms,
			   gnx,grpindex);
    
    if (bGlas)
      do_glas(log,START(nsb),HOMENR(nsb),x,f,
	      fr,mdatoms,top->idef.atnr,&parm->ir,ener);
	         
    if (bTCR && (step == 0)) {
      tcr=init_coupling(log,nfile,fnm,cr,fr,mdatoms,&(top->idef));
      fprintf(log,"Done init_coupling\n"); 
      fflush(log);
    }
#endif


    /* Now we have the energies and forces corresponding to the 
     * coordinates at time t. We must output all of this before
     * the update.
     * for RerunMD t is read from input trajectory
     */
    if (!bRerunMD)
      t        = t0   + step*parm->ir.delta_t;
    
    if (parm->ir.efep != efepNO) {
      if (bRerunMD && rerun_fr.bLambda)
	lambda = rerun_fr.lambda;
      else
	lambda = lam0 + step*parm->ir.delta_lambda;
    }
    if (parm->ir.bSimAnn) {
      SAfactor = 1.0  - t/parm->ir.zero_temp_time;
      if (SAfactor < 0) 
	SAfactor = 0;
    }

    if (MASTER(cr) && do_log)
      print_ebin_header(log,step,t,lambda,SAfactor);
    
    if (bDummies) { 
      /* Spread the force on dummy particle to the other particles... 
       * This is parallellized. MPI communication is performed
       * if the constructing atoms aren't local.
       */
      if(dummycomm)
	move_dummy_f(dummycomm,f,cr);

      spread_dummy_f(log,x,f,buf,&mynrnb,&top->idef);

      if(dummycomm)
        move_construct_f(dummycomm,f,cr);
    }    
    
    if (bLateVir)
      calc_virial(log,START(nsb),HOMENR(nsb),x,f,
		  force_vir,pme_vir,cr,graph,parm->box,&mynrnb,fr,bTweak);
		  
    if (bDummies && fr->bEwald) { 
      /* Spread the LR force on dummy particle to the other particles... 
       * This is parallellized. MPI communication is performed
       * if the constructing atoms aren't local.
       */
      if(dummycomm)
	move_dummy_f(dummycomm,fr->f_pme,cr);

      spread_dummy_f(log,x,fr->f_pme,buf,&mynrnb,&top->idef);

      if(dummycomm)
        move_construct_f(dummycomm,fr->f_pme,cr);
    }    
    if (!bTweak)
      sum_lrforces(f,fr,START(nsb),HOMENR(nsb));
    
    xx = (do_per_step(step,parm->ir.nstxout) || bLastStep) ? x : NULL;
    vv = (do_per_step(step,parm->ir.nstvout) || bLastStep) ? v : NULL;
    ff = (do_per_step(step,parm->ir.nstfout)) ? f : NULL;

    fp_trn = write_traj(log,cr,traj,nsb,step,t,lambda,
			nrnb,nsb->natoms,xx,vv,ff,parm->box);
    debug_gmx();
    
    /* don't write xtc and last structure for rerunMD */
    if (!bRerunMD) {
      if (do_per_step(step,parm->ir.nstxtcout))
	write_xtc_traj(log,cr,xtc_traj,nsb,mdatoms,
		       step,t,x,parm->box,parm->ir.xtcprec);
      if (bLastStep && MASTER(cr)) {
	fprintf(stderr,"Writing final coordinates.\n");
	write_sto_conf(ftp2fn(efSTO,nfile,fnm),
		       *top->name, &(top->atoms),x,v,parm->box);
      }
      debug_gmx();
    }
    clear_mat(shake_vir);
    
    /* Afm and Umbrella type pulling happens before the update, 
	 other types in update */
    if (pulldata.bPull && 
	(pulldata.runtype == eAfm || pulldata.runtype == eUmbrella ||
	 pulldata.runtype == eTest))
      pull(&pulldata,x,f,parm->box,top,parm->ir.delta_t,step,
	   mdatoms->nr,mdatoms); 
    
#ifdef XMDRUN
      /* Check magnitude of the forces */
    fmax = f_max(cr->left,cr->right,cr->nnodes,START(nsb),
		 START(nsb)+HOMENR(nsb),f);
    debug_gmx();
    parm->ir.delta_t = timestep;
#endif
    
    /* This is also parallellized, but check code in update.c */
    update(nsb->natoms,START(nsb),HOMENR(nsb),step,lambda,&ener[F_DVDL],
	   parm,SAfactor,mdatoms,
           x,graph,f,buf,vold,vt,v,
	   top,grps,shake_vir,cr,&mynrnb,bTYZ,TRUE,edyn,&pulldata,bNEMD);
    /* The coordinates (x) were unshifted in update */

    /* Non-equilibrium MD: this is parallellized, but only does communication
     * when there really is NEMD.
     */
    if (PAR(cr) && bNEMD) 
      accumulate_u(cr,&(parm->ir.opts),grps);
      
    debug_gmx();
    if (grps->cosacc.cos_accel == 0)
      calc_ke_part(FALSE,parm->ir.eI==eiSD,
		   START(nsb),HOMENR(nsb),vold,v,vt,&(parm->ir.opts),
		   mdatoms,grps,&mynrnb,lambda,&ener[F_DVDLKIN]);
    else
      calc_ke_part_visc(FALSE,START(nsb),HOMENR(nsb),
			parm->box,x,vold,v,vt,&(parm->ir.opts),
			mdatoms,grps,&mynrnb,lambda,&ener[F_DVDLKIN]);

    debug_gmx();
    /* Calculate center of mass velocity if necessary, also parallellized */
    if (bStopCM)
      calc_vcm_grp(log,HOMENR(nsb),START(nsb),mdatoms->massT,v,vcm);

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
	  grps->estat.ee[egLR][i]   /= cr->nnodes;
	  grps->estat.ee[egLJLR][i] /= cr->nnodes;
	}
    }
    else
      cp_nrnb(&(nrnb[0]),&mynrnb);
    
    /* This is just for testing. Nothing is actually done to Ekin
     * since that would require extra communication.
     */
    if (!bNEMD && debug)
      correct_ekin(debug,START(nsb),START(nsb)+HOMENR(nsb),v,
		   vcm->group_mvcm[0],
		   mdatoms->massT,mdatoms->tmass,parm->ekin);
    
    if ((terminate != 0) && (step < parm->ir.nsteps)) {
      if (terminate<0 && parm->ir.nstxout)
	/* this is the USR1 signal and we are writing x to trr, 
	   stop at next x frame in trr */
	parm->ir.nsteps = (step / parm->ir.nstxout + 1) * parm->ir.nstxout;
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
    if (bStopCM) {
      check_cm_grp(log,vcm);
      do_stopcm_grp(log,HOMENR(nsb),START(nsb),v,vcm,mdatoms->invmass);
      inc_nrnb(&mynrnb,eNR_STOPCM,HOMENR(nsb));
    }
        
    /* Do fit to remove overall rotation */
    if (bStopRot) {
      /* this check is also in grompp.c, if it becomes obsolete here,
	 also remove it there */
      if (PAR(cr))
	  fatal_error(0,"Can not stop rotation about center of mass in a "
		      "parallel run\n");
      do_stoprot(log,top->atoms.nr,box_size,x,mdatoms->massT);
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
    
#ifdef XMDRUN
    /* Check for excessively large energies */
    if (fabs(ener[F_ETOT]) > 1e18) {
      fprintf(stderr,"Energy too large (%g), giving up\n",ener[F_ETOT]);
      break;
    }
#endif
      
    /* Calculate Temperature coupling parameters lambda and adjust
     * target temp when doing simulated annealing
     */
    if(parm->ir.etc==etcBERENDSEN)
      berendsen_tcoupl(&(parm->ir.opts),grps,parm->ir.delta_t,SAfactor);
    else if(parm->ir.etc==etcNOSEHOOVER)
      nosehoover_tcoupl(&(parm->ir.opts),grps,parm->ir.delta_t,SAfactor);

    /* Calculate pressure and apply LR correction if PPPM is used */
    calc_pres(fr->ePBC,parm->box,parm->ekin,parm->vir,parm->pres,
	      (fr->eeltype==eelPPPM) ? ener[F_LR] : 0.0);
    
#ifdef XMDRUN    
    /* Calculate long range corrections to pressure and energy */
    if (bTCR)
      set_avcsix(log,fr,mdatoms);
#endif

    /* Calculate long range corrections to pressure and energy */
    calc_dispcorr(log,parm->ir.eDispCorr,
		  fr,mdatoms->nr,parm->box,parm->pres,parm->vir,ener);

#ifdef XMDRUN
    /* Only do GCT when the relaxation of shells (minimization) has converged,
     * otherwise we might be coupling agains bogus energies. 
     * In parallel we must always do this, because the other sims might
     * update the FF.
     */
    if (bTCR)
      do_coupling(log,nfile,fnm,tcr,t,step,ener,fr,
		  &(parm->ir),MASTER(cr) || bMultiSim,
		  mdatoms,&(top->idef),mu_aver,
		  top->blocks[ebMOLS].nr,bMultiSim ? cr_msim : cr,
		  parm->box,parm->vir,parm->pres,
		  mu_tot,x,f,bConverged);
    debug_gmx();
#endif

    /* Time for performance */
    if (((step % 10) == 0) || bLastStep)
      update_time();

    /* Output stuff */
    if ( MASTER(cr) ) {
      bool do_ene,do_dr;
      
      upd_mdebin(mdebin,fp_dgdl,mdatoms->tmass,step,t,ener,parm->box,shake_vir,
		 force_vir,parm->vir,parm->pres,grps,mu_tot,
		 (parm->ir.etc==etcNOSEHOOVER));
      do_ene = do_per_step(step,parm->ir.nstenergy) || bLastStep;
      if (top->idef.il[F_DISRES].nr)
	do_dr = do_per_step(step,parm->ir.nstdisreout) || bLastStep;
      else
	do_dr = FALSE; 
      print_ebin(fp_ene,do_ene,do_dr,do_log?log:NULL,step,t,
		 eprNORMAL,bCompact,mdebin,&(top->atoms));
      if (bVerbose)
	fflush(log);
    }
    
    /* Remaining runtime */
    if (MASTER(cr) && bVerbose && ( ((step % stepout)==0) || bLastStep)) {
#ifdef XMDRUN
      if (nshell > 0)
	fprintf(stderr,"\n");
#endif
      print_time(stderr,start_t,step,&parm->ir);
    }
    
    bFirstStep = FALSE;

    if (bRerunMD) 
      /* read next frame from input trajectory */
      bNotLastFrame = read_next_frame(status,&rerun_fr);
    
    if (!bRerunMD || !rerun_fr.bStep)
      /* increase the MD step number */
      step++;
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
    print_ebin(fp_ene,FALSE,FALSE,log,step,t,
	       eprAVER,FALSE,mdebin,&(top->atoms));
    print_ebin(fp_ene,FALSE,FALSE,log,step,t,
	       eprRMS,FALSE,mdebin,&(top->atoms));
    close_enx(fp_ene);
    if (!bRerunMD && parm->ir.nstxtcout)
      close_xtc_traj();
    close_trn(fp_trn);
    if (fp_dgdl)
      fclose(fp_dgdl);
  }
  debug_gmx();

#ifdef XMDRUN
  fprintf(log,"Fraction of iterations that converged:           %.2f %%\n",
	  (nconverged*100.0)/(parm->ir.nsteps+1));
  fprintf(log,"Average number of force evaluations per MD step: %.2f\n",
	  tcount/(parm->ir.nsteps+1));
#endif

  return start_t;
}




