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
#include "rdgroup.h"
#include "dummies.h"
#include "update.h"
#include "ns.h"
#include "trnio.h"
#include "mdrun.h"
#include "confio.h"
#include "network.h"
#include "sim_util.h"
#include "pull.h"
#include "xvgr.h"
#include "physics.h"
#include "names.h"
#include "xmdrun.h"

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

static void init_md(t_commrec *cr,t_inputrec *ir,tensor box,real *t,real *t0,
		    real *lambda,real *lam0,real *SAfactor,
		    t_nrnb *mynrnb,bool *bTYZ,t_topology *top,
		    int nfile,t_filenm fnm[],char **traj,
		    char **xtc_traj,int *fp_ene,
		    FILE **fp_dgdl,t_mdebin **mdebin,t_groups *grps,
		    tensor force_vir,tensor pme_vir,
		    tensor shake_vir,t_mdatoms *mdatoms,rvec mu_tot,
		    bool *bNEMD,t_vcm **vcm,t_nsborder *nsb)
{
  bool bBHAM,b14,bLR,bLJLR;
  int  i;
  
  /* Initial values */
  *t = *t0       = ir->init_t;
  if (ir->efep != efepNO) {
    *lambda = *lam0 = ir->init_lambda;
  }
  else {
    *lambda = *lam0   = 0.0;
  } 
  if (ir->bSimAnn) {
    *SAfactor = 1.0 - *t0/ir->zero_temp_time;
    if (*SAfactor < 0) 
      *SAfactor = 0;
  } else
    *SAfactor     = 1.0;
    
  init_nrnb(mynrnb);
  
  /* Check Environment variables & other booleans */
  *bTYZ=getenv("TYZ") != NULL;
  set_pot_bools(ir,top,&bLR,&bLJLR,&bBHAM,&b14);
  
  if (nfile != -1) {
    /* Filenames */
    *traj     = ftp2fn(efTRN,nfile,fnm);
    *xtc_traj = ftp2fn(efXTC,nfile,fnm);
    
    if (MASTER(cr)) {
      *fp_ene = open_enx(ftp2fn(efENX,nfile,fnm),"w");
      if ((fp_dgdl != NULL) && ir->efep!=efepNO)
	*fp_dgdl =
	  xvgropen(opt2fn("-dgdl",nfile,fnm),
		   "dG/d\\8l\\4","Time (ps)",
		   "dG/d\\8l\\4 (kJ mol\\S-1\\N nm\\S-2\\N \\8l\\4\\S-1\\N)");
    } else
      *fp_ene = -1;

    *mdebin = init_mdebin(*fp_ene,grps,&(top->atoms),&(top->idef),
			  bLR,bLJLR,bBHAM,b14,ir->efep!=efepNO,ir->epc,
			  ir->eDispCorr,(TRICLINIC(ir->compress) || TRICLINIC(box)),
			  (ir->etc==etcNOSEHOOVER),cr);
  }
  
  /* Initiate variables */  
  clear_mat(force_vir);
  clear_mat(pme_vir);
  clear_mat(shake_vir);
  clear_rvec(mu_tot);
  
  /* Set initial values for invmass etc. */
  init_mdatoms(mdatoms,*lambda,TRUE);

  *vcm = init_vcm(stdlog,top,mdatoms,START(nsb),HOMENR(nsb),ir->nstcomm);
    
  debug_gmx();

  *bNEMD = (ir->opts.ngacc > 1) || (norm(ir->opts.acc[0]) > 0);

  if (ir->eI == eiSD)
    init_sd_consts(ir->opts.ngtc,ir->opts.tau_t,ir->delta_t);

}

time_t do_md(FILE *log,t_commrec *cr,t_commrec *mcr,int nfile,t_filenm fnm[],
	     bool bVerbose,bool bCompact,
	     bool bDummies, t_comm_dummies *dummycomm,
	     int stepout,t_parm *parm,t_groups *grps,t_topology *top,
	     real ener[],t_fcdata *fcd,
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
  bool       bNS,bStopCM,bTYZ,bRerunMD,bNotLastFrame=FALSE,
             bFirstStep,bLastStep,bNEMD,do_log,bRerunWarnNoV=TRUE;
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
 
  /* XMDRUN stuff: shell, general coupling etc. */
  bool        bFFscan;
  int         nshell,nshell_tot,count,nconverged=0;
  t_shell     *shells=NULL;
  real        timestep=0;
  double      tcount=0;
  bool        bShell,bIonize=FALSE,bGlas=FALSE;
  bool        bTCR=FALSE,bConverged=FALSE,bOK;
  real        mu_aver=0,fmax;
  int         gnx,ii;
  atom_id     *grpindex;
  char        *grpname;
  t_coupl_rec *tcr=NULL;
  rvec        *xcopy=NULL;
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
  init_md(cr,&parm->ir,parm->box,&t,&t0,&lambda,&lam0,&SAfactor,
	  &mynrnb,&bTYZ,top,
	  nfile,fnm,&traj,&xtc_traj,&fp_ene,&fp_dgdl,&mdebin,grps,
	  force_vir,pme_vir,shake_vir,mdatoms,mu_tot,&bNEMD,&vcm,nsb);
  debug_gmx();

  /* Check for polarizable models */
  shells     = init_shells(log,START(nsb),HOMENR(nsb),&top->idef,
			   mdatoms,&nshell);
  nshell_tot = nshell;
  if (PAR(cr))
    gmx_sumi(1,&nshell_tot,cr);
  bShell = nshell_tot > 0;
  
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
  
  /* Initiate data for the special cases */
  if (bFFscan) {
    snew(xcopy,nsb->natoms);
    for(ii=0; (ii<nsb->natoms); ii++)
      copy_rvec(x[ii],xcopy[ii]);
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
    bStopCM = do_per_step(step,abs(parm->ir.nstcomm));

    /* Copy back starting coordinates in case we're doing a forcefield scan */
    if (bFFscan) {
      for(ii=0; (ii<nsb->natoms); ii++)
	copy_rvec(xcopy[ii],x[ii]);
    }
    
    if (bDummies) {
      shift_self(graph,parm->box,x);
    
      construct_dummies(log,x,&mynrnb,parm->ir.delta_t,v,&top->idef,
			graph,cr,parm->box,dummycomm);
      
      unshift_self(graph,parm->box,x);
    }
     
    debug_gmx();
    
    /* Set values for invmass etc. This routine not parallellized, but hardly
     * ever used, only when doing free energy calculations.
     */
    init_mdatoms(mdatoms,lambda,bFirstStep);
    
    clear_mat(force_vir);
    
    /* Ionize the atoms if necessary */
    if (bIonize)
      ionize(log,mdatoms,top->atoms.atomname,t,&parm->ir,x,v,
	     START(nsb),START(nsb)+HOMENR(nsb),parm->box,cr);
      
    /* Update force field in ffscan program */
    if (bFFscan) 
      update_forcefield(nfile,fnm,fr);
    
    if (bShell) {
      /* Now is the time to relax the shells */
      count=relax_shells(log,cr,mcr,bVerbose,step,parm,bNS,bStopCM,top,
			 ener,fcd,x,vold,v,vt,f,buf,mdatoms,nsb,&mynrnb,graph,
			 grps,force_vir,pme_vir,bShell,
			 nshell,shells,fr,traj,t,lambda,mu_tot,
			 nsb->natoms,parm->box,&bConverged);
      tcount+=count;
      
      if (bConverged)
	nconverged++;
    }
    else {
      /* The coordinates (x) are shifted (to get whole molecules) in do_force
       * This is parallellized as well, and does communication too. 
       * Check comments in sim_util.c
       */
      do_force(log,cr,mcr,parm,nsb,force_vir,pme_vir,step,&mynrnb,
	       top,grps,x,v,f,buf,mdatoms,ener,fcd,bVerbose && !PAR(cr),
	       lambda,graph,bNS,FALSE,fr,mu_tot,FALSE);
    }
   
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

    if (MASTER(cr) && do_log && !bFFscan)
      print_ebin_header(log,step,t,lambda,SAfactor);
    
    if (bDummies) 
      spread_dummy_f(log,x,f,&mynrnb,&top->idef,dummycomm,cr);
      
    /* Calculation of the virial must be done after dummies!    */
    /* Question: Is it correct to do the PME forces after this? */
    calc_virial(log,START(nsb),HOMENR(nsb),x,f,
		force_vir,pme_vir,graph,parm->box,&mynrnb,fr,FALSE);
		  
    /* Spread the LR force on dummy particle to the other particles... 
     * This is parallellized. MPI communication is performed
     * if the constructing atoms aren't local.
     */
    if (bDummies && fr->bEwald) 
      spread_dummy_f(log,x,fr->f_pme,&mynrnb,&top->idef,dummycomm,cr);
    
    sum_lrforces(f,fr,START(nsb),HOMENR(nsb));

    xx = (do_per_step(step,parm->ir.nstxout) || bLastStep) ? x : NULL;
    vv = (do_per_step(step,parm->ir.nstvout) || bLastStep) ? v : NULL;
    ff = (do_per_step(step,parm->ir.nstfout)) ? f : NULL;

    fp_trn = write_traj(log,cr,traj,nsb,step,t,lambda,
			nrnb,nsb->natoms,xx,vv,ff,parm->box);
    debug_gmx();
    
    /* don't write xtc and last structure for rerunMD */
    if (!bRerunMD && !bFFscan) {
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
     * other types in update 
     */
    if (pulldata.bPull && 
	(pulldata.runtype == eAfm || pulldata.runtype == eUmbrella ||
	 pulldata.runtype == eTest))
      pull(&pulldata,x,f,parm->box,top,parm->ir.delta_t,step,
	   mdatoms->nr,mdatoms); 
    
    if (bFFscan)
      clear_rvecs(nsb->natoms,buf);
      
    /* This is also parallellized, but check code in update.c */
    /* bOK = update(nsb->natoms,START(nsb),HOMENR(nsb),step,lambda,&ener[F_DVDL], */
    bOK = TRUE;
    update(nsb->natoms,START(nsb),HOMENR(nsb),step,lambda,&ener[F_DVDL],
		 parm,SAfactor,mdatoms,x,graph,f,buf,vold,vt,v,
		 top,grps,shake_vir,cr,&mynrnb,bTYZ,TRUE,edyn,&pulldata,bNEMD);
    if (!bOK && !bFFscan)
      fatal_error(0,"Constraint error: Shake, Lincs or Settle could not solve the constrains");
    
    /* The coordinates (x) were unshifted in update */
    if (bFFscan && (!bShell || bConverged))
      print_forcefield(log,ener[F_EPOT],HOMENR(nsb),f,buf,xcopy,
		       &(top->blocks[ebMOLS]),mdatoms->massT); 
    
    if (parm->ir.epc!=epcNO)
      correct_box(parm->box,fr,graph);
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
		   START(nsb),HOMENR(nsb),vold,v,vt,&(parm->ir.opts),
		   mdatoms,grps,&mynrnb,lambda,&ener[F_DVDLKIN]);
    else
      calc_ke_part_visc(FALSE,START(nsb),HOMENR(nsb),
			parm->box,x,vold,v,vt,&(parm->ir.opts),
			mdatoms,grps,&mynrnb,lambda,&ener[F_DVDLKIN]);

    debug_gmx();
    /* Calculate center of mass velocity if necessary, also parallellized */
    if (bStopCM && !bFFscan)
      calc_vcm_grp(log,START(nsb),HOMENR(nsb),mdatoms->massT,x,v,vcm);

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
		   vcm->group_p[0],
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
    if (bStopCM && !bFFscan) {
      check_cm_grp(log,vcm);
      do_stopcm_grp(log,START(nsb),HOMENR(nsb),x,v,vcm);
      inc_nrnb(&mynrnb,eNR_STOPCM,HOMENR(nsb));
      calc_vcm_grp(log,START(nsb),HOMENR(nsb),mdatoms->massT,x,v,vcm);
      check_cm_grp(log,vcm);
      do_stopcm_grp(log,START(nsb),HOMENR(nsb),x,v,vcm);
      check_cm_grp(log,vcm);
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
    if(parm->ir.etc==etcBERENDSEN)
      berendsen_tcoupl(&(parm->ir.opts),grps,parm->ir.delta_t,SAfactor);
    else if(parm->ir.etc==etcNOSEHOOVER)
      nosehoover_tcoupl(&(parm->ir.opts),grps,parm->ir.delta_t,SAfactor);

    /* Calculate pressure and apply LR correction if PPPM is used */
    calc_pres(fr->ePBC,parm->box,parm->ekin,parm->vir,parm->pres,
	      (fr->eeltype==eelPPPM) ? ener[F_LR] : 0.0);
    
    /* Calculate long range corrections to pressure and energy */
    if (bTCR)
      set_avcsix(log,fr,mdatoms);
    
    /* Calculate long range corrections to pressure and energy */
    calc_dispcorr(log,parm->ir.eDispCorr,
		  fr,mdatoms->nr,parm->box,parm->pres,parm->vir,ener);

    if (bTCR) {
      /* Only do GCT when the relaxation of shells (minimization) has converged,
       * otherwise we might be coupling to bogus energies. 
       * In parallel we must always do this, because the other sims might
       * update the FF.
       */
      do_coupling(log,nfile,fnm,tcr,t,step,ener,fr,
		  &(parm->ir),MASTER(cr),
		  mdatoms,&(top->idef),mu_aver,
		  top->blocks[ebMOLS].nr,(mcr!=NULL) ? mcr : cr,
		  parm->box,parm->vir,parm->pres,
		  mu_tot,x,f,bConverged);
      debug_gmx();
    }

    /* Time for performance */
    if (((step % 10) == 0) || bLastStep)
      update_time();

    /* Output stuff */
    if (MASTER(cr) && !bFFscan) {
      bool do_ene,do_dr,do_or;
      
      upd_mdebin(mdebin,fp_dgdl,mdatoms->tmass,step,t,ener,parm->box,shake_vir,
		 force_vir,parm->vir,parm->pres,grps,mu_tot,
		 (parm->ir.etc==etcNOSEHOOVER));
      do_ene = do_per_step(step,parm->ir.nstenergy) || bLastStep;
      do_dr  = do_per_step(step,parm->ir.nstdisreout) || bLastStep;
      do_or  = do_per_step(step,parm->ir.nstorireout) || bLastStep;
      print_ebin(fp_ene,do_ene,do_dr,do_or,do_log?log:NULL,step,t,
		 eprNORMAL,bCompact,mdebin,fcd,&(top->atoms));
      if (bVerbose)
	fflush(log);
    }
    
    /* Remaining runtime */
    if (MASTER(cr) && bVerbose && ( ((step % stepout)==0) || bLastStep)) {
      if (bShell)
	fprintf(stderr,"\n");
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
    print_ebin(fp_ene,FALSE,FALSE,FALSE,log,step,t,
	       eprAVER,FALSE,mdebin,fcd,&(top->atoms));
    print_ebin(fp_ene,FALSE,FALSE,FALSE,log,step,t,
	       eprRMS,FALSE,mdebin,fcd,&(top->atoms));
    close_enx(fp_ene);
    if (!bRerunMD && parm->ir.nstxtcout)
      close_xtc_traj();
    close_trn(fp_trn);
    if (fp_dgdl)
      fclose(fp_dgdl);
  }
  debug_gmx();

  if (bShell) {
    fprintf(log,"Fraction of iterations that converged:           %.2f %%\n",
	    (nconverged*100.0)/(parm->ir.nsteps+1));
    fprintf(log,"Average number of force evaluations per MD step: %.2f\n",
	    tcount/(parm->ir.nsteps+1));
  }

  return start_t;
}
