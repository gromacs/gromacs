/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_md_c = "$Id$";

#include <signal.h>
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

volatile bool bGotTermSignal = FALSE, bGotUsr1Signal = FALSE;

static void signal_handler(int n)
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
	     bool bVerbose,bool bCompact,bool bDummies,int stepout,
	     t_parm *parm,t_groups *grps,t_topology *top,real ener[],
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
  bool       bNS,bStopCM,bStopRot,bTYZ,bRerunMD,bNotLastFrame=FALSE,bLastStep,
             bNEMD,do_log;
  tensor     force_vir,shake_vir;
  t_nrnb     mynrnb;
  char       *traj,*xtc_traj; /* normal and compressed trajectory filename */
  int        i,m;
  rvec       vcm,mu_tot;
  rvec       *xx,*vv,*ff;
  int        natoms=0,status;
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
  real        mu_aver,fmax;
  int         gnx;
  atom_id     *grpindex;
  char        *grpname;
  t_coupl_rec *tcr=NULL;
#endif
  
  /* Turn on signal handling */
  signal(SIGTERM,signal_handler);
  signal(SIGUSR1,signal_handler);

  /* check if "-rerun" is used. */
  bRerunMD = (Flags & MD_RERUN) == MD_RERUN;

  /* Initial values */
  init_md(cr,&parm->ir,&t,&t0,&lambda,&lam0,&SAfactor,&mynrnb,&bTYZ,top,
	  nfile,fnm,&traj,&xtc_traj,&fp_ene,&fp_dgdl,&mdebin,grps,vcm,
	  force_vir,shake_vir,mdatoms,mu_tot,&bNEMD);
  debug_gmx();
  
#ifdef XMDRUN
  bDynamicStep = FALSE;
  bIonize      = (Flags & MD_IONIZE)   == MD_IONIZE;
  bMultiSim    = (Flags & MD_MULTISIM) == MD_MULTISIM;
  bGlas        = (Flags & MD_GLAS)     == MD_GLAS;
  timestep = parm->ir.delta_t;

  /* Check whether we have to do dipole stuff */
  if (ftp2bSet(efNDX,nfile,fnm))
    rd_index(ftp2fn(efNDX,nfile,fnm),1,&gnx,&grpindex,&grpname);
  else 
    gnx = 0;
  
  /* Check whether we have to GCT stuff */
  bTCR = ftp2bSet(efGCT,nfile,fnm);
  if (MASTER(cr) && bTCR)
    fprintf(stderr,"Will do General Coupling Theory!\n");
  
#endif

  /* Remove periodicity */  
  if (fr->ePBC != epbcNONE)
    do_pbc_first(log,parm,box_size,fr,graph,x);
  debug_gmx();

  /* Initialize pull code */
  init_pull(log,nfile,fnm,&pulldata,x,mdatoms,parm->box); 
  
  if (!parm->ir.bUncStart) 
    do_shakefirst(log,bTYZ,lambda,ener,parm,nsb,mdatoms,x,vold,buf,f,v,
		  graph,cr,&mynrnb,grps,fr,top,edyn,&pulldata);
  debug_gmx();
    
  /* Compute initial EKin for all.. */
  if (grps->cosacc.cos_accel == 0)
    calc_ke_part(TRUE,START(nsb),HOMENR(nsb),vold,v,vt,&(parm->ir.opts),
		 mdatoms,grps,&mynrnb,lambda,&ener[F_DVDLKIN]);
  else
    calc_ke_part_visc(TRUE,START(nsb),HOMENR(nsb),
		      parm->box,x,vold,v,vt,&(parm->ir.opts),
		      mdatoms,grps,&mynrnb,lambda,&ener[F_DVDLKIN]);
  debug_gmx();
	
  if (PAR(cr)) 
    global_stat(log,cr,ener,force_vir,shake_vir,
		&(parm->ir.opts),grps,&mynrnb,nrnb,vcm,mu_tot,&terminate);
  clear_rvec(vcm);
  debug_gmx();
  
  /* Calculate Temperature coupling parameters lambda */
  ener[F_TEMP] = sum_ekin(&(parm->ir.opts),grps,parm->ekin,bTYZ);
  tcoupl(parm->ir.btc,&(parm->ir.opts),grps,parm->ir.delta_t,SAfactor);
  debug_gmx();
  
#ifdef XMDRUN
  shells = init_shells(log,START(nsb),HOMENR(nsb),&top->idef,mdatoms,&nshell);
#endif
  
  /* Write start time and temperature */
  start_t=print_date_and_time(log,cr->pid,
#ifdef XMDRUN
			      "Started Xmdrun"
#else
			      "Started mdrun"
#endif		      
			      );
  
  if (MASTER(cr)) {
    fprintf(log,"Initial temperature: %g K\n",ener[F_TEMP]);
    if (bRerunMD) 
      fprintf(stderr,"starting md rerun '%s', reading coordinates from"
	      " input trajectory '%s'\n\n",
	      *(top->name),opt2fn("-rerun",nfile,fnm));
    else
      fprintf(stderr,"starting mdrun '%s'\n%d steps, %8.1f ps.\n\n",
	      *(top->name),parm->ir.nsteps,parm->ir.nsteps*parm->ir.delta_t);
  }
  /* Set the cpu time counter to 0 after initialisation */
  start_time();
  debug_gmx();
  
  /***********************************************************
   *
   *             Loop over MD steps 
   *
   ************************************************************/
  
  /* if rerunMD then read coordinates from input trajectory,
   * set velocities to zero */
  if (bRerunMD) {
    int i;

    natoms=read_first_x(&status,opt2fn("-rerun",nfile,fnm),&t,&x,parm->box);
    if (natoms != mdatoms->nr)
      fatal_error(0,"Number of atoms in trajectory (%d) does not match the "
		  "run input file (%d)\n",natoms,mdatoms->nr);
    for(i=0;(i<natoms);i++) 
      clear_rvec(v[i]);

    bNotLastFrame = (natoms!=0);
  } 

  /* loop over MD steps or if rerunMD to end of input trajectory */
  for(step=0; ((!bRerunMD && (step<=parm->ir.nsteps)) || 
	       (bRerunMD && bNotLastFrame)); step++) {
	       
    bLastStep=(step==parm->ir.nsteps);

    do_log = do_per_step(step,parm->ir.nstlog) || bLastStep;
    
    if (bRerunMD)
      /* for rerun MD always do Neighbour Searching */
      bNS = ((parm->ir.nstlist!=0) || (step==0));
    else {
      /* Stop Center of Mass motion */
      get_cmparm(&parm->ir,step,&bStopCM,&bStopRot);
    
      /* Determine whether or not to do Neighbour Searching */
      bNS=((parm->ir.nstlist && ((step % parm->ir.nstlist)==0)) || (step==0));
    }
    
    if (bDummies) {
      /* Construct dummy particles. This is parallellized */
      shift_self(graph,parm->box,x);
      construct_dummies(log,x,&mynrnb,parm->ir.delta_t,v,&top->idef);
      unshift_self(graph,parm->box,x);
    }
    debug_gmx();
    
    /* Set values for invmass etc. This routine not parallellized, but hardly
     * ever used, only when doing free energy calculations.
     */
    init_mdatoms(mdatoms,lambda,(step==0));
    
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
		       nshell,shells,fr,traj,t,lambda,nsb->natoms,parm->box,
		       &bConverged);
    tcount+=count;
    
    if (bConverged)
      nconverged++;

#else

    /* The coordinates (x) are shifted (to get whole molecules) in do_force
     * This is parallellized as well, and does communication too. 
     * Check comments in sim_util.c
     */
    do_force(log,cr,parm,nsb,force_vir,step,&mynrnb,
	     top,grps,x,v,f,buf,mdatoms,ener,bVerbose && !PAR(cr),
	     lambda,graph,bNS,FALSE,fr);
    debug_gmx();
#ifdef DEBUG
    pr_rvecs(log,0,"force_vir",force_vir,DIM);
#endif     

#endif

    /* Calculate total (local) dipole moment if necessary. 
     * This is parallellized 
     */
    calc_mu(nsb,x,mdatoms->chargeT,mu_tot);
    
#ifdef XMDRUN
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
    
    if (parm->ir.efep != efepNO)
      lambda   = lam0 + step*parm->ir.delta_lambda;
    if (parm->ir.bSimAnn) {
      SAfactor = 1.0  - t/parm->ir.zero_temp_time;
      if (SAfactor < 0) 
	SAfactor = 0;
    }

    if (MASTER(cr) && do_log)
      print_ebin_header(log,step,t,lambda,SAfactor);
    
    if (bDummies)
      /* Spread the force on dummy particle to the other particles... 
       * This is parallellized
       */
      spread_dummy_f(log,x,f,&mynrnb,&top->idef);
    
    xx = (do_per_step(step,parm->ir.nstxout) || bLastStep) ? x : NULL;
    vv = (do_per_step(step,parm->ir.nstvout) || bLastStep) ? v : NULL;
    ff = (do_per_step(step,parm->ir.nstfout)) ? f : NULL;
    fp_trn = write_traj(log,cr,traj,nsb,step,t,lambda,
			nrnb,nsb->natoms,xx,vv,ff,parm->box);
    debug_gmx();
    
    /* for rerunMD, certain things don't have to be done */
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
      
      clear_mat(shake_vir);

      /* Afm and Umbrella type pulling happens before the update, 
	 other types in update */
      if (pulldata.bPull && 
	  (pulldata.runtype == eAfm || pulldata.runtype == eUmbrella ||
	   pulldata.runtype == eTest))
	pull(&pulldata,x,f,parm->box,top,parm->ir.delta_t,step,
	     natoms,mdatoms); 

#ifdef XMDRUN
      /* Check magnitude of the forces */
      fmax = f_max(cr->left,cr->right,cr->nprocs,START(nsb),
		   START(nsb)+HOMENR(nsb),f);
      debug_gmx();
      parm->ir.delta_t = timestep;
#endif

      /* This is also parallellized, but check code in update.c */
      update(nsb->natoms,START(nsb),HOMENR(nsb),step,lambda,&ener[F_DVDL],
	     &(parm->ir),mdatoms,x,graph,f,
	     buf,vold,v,vt,parm->pres,parm->box,
	     top,grps,shake_vir,cr,&mynrnb,bTYZ,TRUE,edyn,&pulldata,bNEMD);
      /* The coordinates (x) were unshifted in update */

#ifdef DEBUG
      pr_rvecs(log,0,"shake_vir",shake_vir,DIM);
#endif
      /* Non-equilibrium MD: this is parallellized, but only does communication
       * when there really is NEMD.
       */
      if (PAR(cr) && bNEMD) 
	accumulate_u(cr,&(parm->ir.opts),grps);
      
      debug_gmx();
      /* Calculate partial Kinetic Energy (for this processor) 
       * per group! Parallelized
       */
      if (grps->cosacc.cos_accel == 0)
	calc_ke_part(FALSE,START(nsb),HOMENR(nsb),vold,v,vt,&(parm->ir.opts),
		     mdatoms,grps,&mynrnb,lambda,&ener[F_DVDLKIN]);
      else
	calc_ke_part_visc(FALSE,START(nsb),HOMENR(nsb),
			  parm->box,x,vold,v,vt,&(parm->ir.opts),
			  mdatoms,grps,&mynrnb,lambda,&ener[F_DVDLKIN]);
      debug_gmx();
      /* Calculate center of mass velocity if necessary, also parallellized */
      if (bStopCM)
	calc_vcm(log,HOMENR(nsb),START(nsb),mdatoms->massT,v,vcm);
    }

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
      /* Globally (over all CPUs) sum energy, virial etc. 
       * This includes communication 
       */
      global_stat(log,cr,ener,force_vir,shake_vir,
		  &(parm->ir.opts),grps,&mynrnb,nrnb,vcm,mu_tot,&terminate);
      /* Correct for double counting energies, should be moved to 
       * global_stat 
       */
      if (fr->bTwinRange && !bNS) 
	for(i=0; (i<grps->estat.nn); i++) {
	  grps->estat.ee[egLR][i]   /= cr->nprocs;
	  grps->estat.ee[egLJLR][i] /= cr->nprocs;
	}
    }
    else
      cp_nrnb(&(nrnb[0]),&mynrnb);
    
    /* This is just for testing. Nothing is actually done to Ekin
     * since that would require extra communication.
     */
    if (!bNEMD && debug)
      correct_ekin(debug,START(nsb),START(nsb)+HOMENR(nsb),v,vcm,
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
    
    if (!bRerunMD) {
      /* Do center of mass motion removal */
      if (bStopCM) {
	check_cm(log,vcm,mdatoms->tmass);
	do_stopcm(log,HOMENR(nsb),START(nsb),v,vcm,
		  mdatoms->tmass,mdatoms->invmass);
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
    }
    
    /* Sum the potential energy terms from group contributions */
    sum_epot(&(parm->ir.opts),grps,ener);

    /* Calculate the amplitude of the cosine velocity profile */
    grps->cosacc.vcos = grps->cosacc.mvcos/mdatoms->tmass;

    if (!bRerunMD) {
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
      
      /* Calculate Temperature coupling parameters lambda */
      tcoupl(parm->ir.btc,&(parm->ir.opts),grps,parm->ir.delta_t,SAfactor);
    
      /* Calculate pressure and apply LR correction if PPPM is used */
      calc_pres(fr->ePBC,parm->box,parm->ekin,parm->vir,parm->pres,
		(fr->eeltype==eelPPPM) ? ener[F_LR] : 0.0);
    }
    
#ifdef XMDRUN    
    /* Calculate long range corrections to pressure and energy */
    if (bTCR)
      set_avcsix(log,fr,mdatoms);
#endif

    /* Calculate long range corrections to pressure and energy */
    calc_dispcorr(log,parm->ir.bDispCorr,
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
		 force_vir,parm->vir,parm->pres,grps,mu_tot);
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
    
    /* if rerunMD read next frame from input trajectory */
    if (bRerunMD) 
      bNotLastFrame = read_next_x(status,&t,natoms,x,parm->box);
  }
  /* End of main MD loop */
  debug_gmx();

  /* Dump the CPU time to the log file on each processor */
  if (PAR(cr)) {
    double *ct,ctmax,ctsum;
    
    snew(ct,cr->nprocs);
    ct[cr->pid] = cpu_time();
    gmx_sumd(cr->nprocs,ct,cr);
    ctmax = ct[0];
    ctsum = ct[0];
    for(i=1; (i<cr->pid); i++) {
      ctmax = max(ctmax,ct[i]);
      ctsum += ct[i];
    }
    ctsum /= cr->nprocs;
    fprintf(log,"\nTotal CPU time on processor %d: %g\n",cr->pid,ct[cr->pid]);
    fprintf(log,"Average CPU time: %g\n",ctsum);
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
