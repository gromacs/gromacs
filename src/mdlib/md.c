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

volatile bool bGotTermSignal = FALSE;

static void signal_handler(int n /* to keep the DEC compiler happy */)
{
  bGotTermSignal = TRUE;
}

time_t do_md(FILE *log,t_commrec *cr,int nfile,t_filenm fnm[],
	     bool bVerbose,bool bCompact,bool bDummies,int stepout,
	     t_parm *parm,t_groups *grps,t_topology *top,real ener[],
	     rvec x[],rvec vold[],rvec v[],rvec vt[],rvec f[],
	     rvec buf[],t_mdatoms *mdatoms,t_nsborder *nsb,t_nrnb nrnb[],
	     t_graph *graph,t_edsamyn *edyn,t_forcerec *fr,rvec box_size)
{
  t_mdebin   *mdebin;
  int        fp_ene=0,fp_trn=0,step;
  FILE       *fp_dgdl=NULL;
  time_t     start_t;
  real       t,lambda,t0,lam0,SAfactor;
  bool       bNS,bStopCM,bStopRot,bTYZ,bRerunMD,bNotLastFrame=FALSE,bLastStep;
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

  /* Turn on signal handling */
  signal(SIGTERM,signal_handler);

  /* check if "-rerun" is used. */
  bRerunMD = optRerunMDset(nfile,fnm);
  
  /* Initial values */
  init_md(cr,&parm->ir,&t,&t0,&lambda,&lam0,&SAfactor,&mynrnb,&bTYZ,top,
	  nfile,fnm,&traj,&xtc_traj,&fp_ene,&fp_dgdl,&mdebin,grps,vcm,
	  force_vir,shake_vir,mdatoms);
  debug_gmx();
  
  /* Remove periodicity */  
  if (fr->eBox != ebtNONE)
    do_pbc_first(log,parm,box_size,fr,graph,x);
  debug_gmx();

  /* Initialize pull code */
  init_pull(log,nfile,fnm,&pulldata,x,mdatoms,box_size); 
  
  if (!parm->ir.bUncStart) 
    do_shakefirst(log,bTYZ,lambda,ener,parm,nsb,mdatoms,x,vold,buf,f,v,
		  graph,cr,&mynrnb,grps,fr,top,edyn,&pulldata);
    
  /* Compute initial EKin for all.. */
  calc_ke_part(TRUE,START(nsb),HOMENR(nsb),vold,v,vt,&(parm->ir.opts),
	       mdatoms,grps,&mynrnb,lambda,&ener[F_DVDLKIN]);
	
  if (PAR(cr)) 
    global_stat(log,cr,ener,force_vir,shake_vir,
		&(parm->ir.opts),grps,&mynrnb,nrnb,vcm,mu_tot,&terminate);
  clear_rvec(vcm);
  debug_gmx();
  
  /* Calculate Temperature coupling parameters lambda */
  ener[F_TEMP] = sum_ekin(&(parm->ir.opts),grps,parm->ekin,bTYZ);
  tcoupl(parm->ir.btc,&(parm->ir.opts),grps,parm->ir.delta_t,SAfactor,0,
	 parm->ir.ntcmemory);
  debug_gmx();

  /* Write start time and temperature */
  start_t=print_date_and_time(log,cr->pid,"Started mdrun");
  if (MASTER(cr)) {
    fprintf(log,"Initial temperature: %g K\n",ener[F_TEMP]);
    if (bRerunMD) 
      fprintf(stderr,"starting md rerun '%s', reading coordinates from input trajectory '%s'\n\n",*(top->name),opt2fn("-rerun",nfile,fnm));
    else
      fprintf(stderr,"starting mdrun '%s'\n%d steps, %8.1f ps.\n\n",
	      *(top->name),parm->ir.nsteps,parm->ir.nsteps*parm->ir.delta_t);
  }
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
    for(i=0;(i<natoms);i++) 
      clear_rvec(v[i]);

    bNotLastFrame = (natoms!=0);
  } 

  /* loop over MD steps or if rerunMD to end of input trajectory */
  for(step=0; ((!bRerunMD && (step<=parm->ir.nsteps)) || 
	       (bRerunMD && bNotLastFrame)); step++) {
	       
    bLastStep=(step==parm->ir.nsteps);

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
      /* Construct dummy particles */
      shift_self(graph,fr->shift_vec,x);
      construct_dummies(log,x,&mynrnb,parm->ir.delta_t,v,&top->idef);
      unshift_self(graph,fr->shift_vec,x);
    }
    debug_gmx();
    
    /* Set values for invmass etc. */
    init_mdatoms(mdatoms,lambda,(step==0));

    /* Calculate total (local) dipole moment if necessary */    
    calc_mu(nsb,x,mdatoms->chargeT,mu_tot);
    
    /* Calc forces and virial
     * The coordinates (x) are shifted (to get whole molecules) in do_force
     */
    clear_mat(force_vir);
    do_force(log,cr,parm,nsb,force_vir,step,&mynrnb,
	     top,grps,x,v,f,buf,mdatoms,ener,bVerbose && !PAR(cr),
	     lambda,graph,bNS,FALSE,fr);
    debug_gmx();
#ifdef DEBUG
    pr_rvecs(log,0,"force_vir",force_vir,DIM);
#endif     
    /* Now we have the energies and forces corresponding to the 
     * coordinates at time t. We must output all of this before
     * the update.
     * for RerunMD t is read from input trajectory
     */
    if (!bRerunMD)
      t        = t0   + step*parm->ir.delta_t;
    
    if (parm->ir.bPert)
      lambda   = lam0 + step*parm->ir.delta_lambda;
    if (parm->ir.bSimAnn) {
      SAfactor = 1.0  - t/parm->ir.zero_temp_time;
      if (SAfactor < 0) 
	SAfactor = 0;
    }

    if (bDummies)
      /* Spread the force on dummy particle to the other particles... */
      spread_dummy_f(log,x,f,&mynrnb,&top->idef);
    
    xx = (do_per_step(step,parm->ir.nstxout) || bLastStep) ? x : NULL;
    vv = (do_per_step(step,parm->ir.nstvout) || bLastStep) ? v : NULL;
    ff = (do_per_step(step,parm->ir.nstfout) || bLastStep) ? f : NULL;
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

      update(nsb->natoms,START(nsb),HOMENR(nsb),step,lambda,&ener[F_DVDL],
	     &(parm->ir),mdatoms,x,graph,
	     fr->shift_vec,f,buf,vold,v,vt,parm->pres,parm->box,
	     top,grps,shake_vir,cr,&mynrnb,bTYZ,TRUE,edyn,&pulldata);
      /* The coordinates (x) were unshifted in update */
#ifdef DEBUG
      pr_rvecs(log,0,"shake_vir",shake_vir,DIM);
#endif
      if (PAR(cr)) 
	accumulate_u(cr,&(parm->ir.opts),grps);
      
      debug_gmx();
      /* Calculate partial Kinetic Energy (for this processor) 
       * per group!
       */
      calc_ke_part(FALSE,START(nsb),HOMENR(nsb),
		   vold,v,vt,&(parm->ir.opts),
		   mdatoms,grps,&mynrnb,
		   lambda,&ener[F_DVDL]);
      debug_gmx();
      if (bStopCM)
	calc_vcm(log,HOMENR(nsb),START(nsb),mdatoms->massT,v,vcm);
    }
    
    if (bGotTermSignal) {
      terminate = 1;
      fprintf(log,"\n\nReceived the TERM signal\n\n");
      if (MASTER(cr))
	fprintf(stderr,"\n\nReceived the TERM signal\n\n");
      bGotTermSignal = FALSE;
    }

    if (PAR(cr)) {
      /* Globally (over all CPUs) sum energy, virial etc. */
      global_stat(log,cr,ener,force_vir,shake_vir,
		  &(parm->ir.opts),grps,&mynrnb,nrnb,vcm,mu_tot,&terminate);
      if (fr->bTwinRange && !bNS) 
	for(i=0; (i<grps->estat.nn); i++) {
	  grps->estat.ee[egLR][i]   /= cr->nprocs;
	  grps->estat.ee[egLJLR][i] /= cr->nprocs;
	}
    }
    else
      cp_nrnb(&(nrnb[0]),&mynrnb);

    if ((terminate > 0) && (step < parm->ir.nsteps)) {
      parm->ir.nsteps = step+1;
      fprintf(log,"\nSetting nsteps to %d\n\n",parm->ir.nsteps);
      if (MASTER(cr))
	fprintf(stderr,"\nSetting nsteps to %d\n\n",parm->ir.nsteps);
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
      if (bStopRot)
	do_stoprot(log,top->atoms.nr,box_size,x,mdatoms->massT);
    
      /* Add force and shake contribution to the virial */
      m_add(force_vir,shake_vir,parm->vir);
    }
    
    /* Sum the potential energy terms from group contributions */
    sum_epot(&(parm->ir.opts),grps,ener);
    
    if (!bRerunMD) {
      /* Sum the kinetic energies of the groups & calc temp */
      ener[F_TEMP]=sum_ekin(&(parm->ir.opts),grps,parm->ekin,bTYZ);
      ener[F_EKIN]=trace(parm->ekin);
      ener[F_ETOT]=ener[F_EPOT]+ener[F_EKIN];
      
      /* Calculate Temperature coupling parameters lambda */
      tcoupl(parm->ir.btc,&(parm->ir.opts),grps,parm->ir.delta_t,SAfactor,
	     step,parm->ir.ntcmemory);
    
      /* Calculate pressure and apply LR correction if PPPM is used */
      calc_pres(fr->eBox,parm->box,parm->ekin,parm->vir,parm->pres,
		(fr->eeltype==eelPPPM) ? ener[F_LR] : 0.0);
    }
    
    /* Calculate long range corrections to pressure and energy */
    calc_dispcorr(log,parm->ir.bDispCorr,
		  fr,mdatoms->nr,parm->box,parm->pres,parm->vir,ener);
    
    if (MASTER(cr))
      upd_mdebin(mdebin,fp_dgdl,mdatoms->tmass,step,t,ener,parm->box,shake_vir,
		 force_vir,parm->vir,parm->pres,grps,mu_tot);
    debug_gmx();
    
    if ( MASTER(cr) ) {
      bool do_ene,do_log,do_dr;
      
      do_ene = do_per_step(step,parm->ir.nstenergy) || bLastStep;
      if (top->idef.il[F_DISRES].nr)
	do_dr = do_per_step(step,parm->ir.nstdisreout) || bLastStep;
      else
	do_dr = FALSE; 
      do_log = do_per_step(step,parm->ir.nstlog) || bLastStep;
      print_ebin(fp_ene,do_ene,do_dr,do_log?log:NULL,step,t,lambda,SAfactor,
		 eprNORMAL,bCompact,mdebin,grps,&(top->atoms));
      if (bVerbose)
	fflush(log);
    }
    
    /* Time for performance */
    if (((step % 10) == 0) || bLastStep)
      update_time();
    
    /* Remaining runtime */
    if (MASTER(cr) && bVerbose && ( ((step % stepout)==0) || bLastStep))
      print_time(stderr,start_t,step,&parm->ir);
    
    /* if rerunMD read next frame from input trajectory */
    if (bRerunMD) 
      bNotLastFrame = read_next_x(status,&t,natoms,x,parm->box);
  }
  /* End of main MD loop */
  debug_gmx();
  
  if (MASTER(cr)) {
    print_ebin(fp_ene,FALSE,FALSE,log,step,t,lambda,SAfactor,
	       eprAVER,FALSE,mdebin,grps,&(top->atoms));
    print_ebin(fp_ene,FALSE,FALSE,log,step,t,lambda,SAfactor,
	       eprRMS,FALSE,mdebin,grps,&(top->atoms));
    close_enx(fp_ene);
    if (!bRerunMD && parm->ir.nstxtcout)
      close_xtc_traj();
    close_trn(fp_trn);
    if (fp_dgdl)
      fclose(fp_dgdl);
  }
  debug_gmx();

  return start_t;
}
