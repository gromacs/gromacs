/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Giant Rising Ordinary Mutants for A Clerical Setup
 */
static char *SRCID_md_c = "$Id$";

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

time_t do_md(FILE *log,t_commrec *cr,int nfile,t_filenm fnm[],
	     bool bVerbose,bool bCompact,bool bDummies,int stepout,
	     t_parm *parm,t_groups *grps,t_topology *top,real ener[],
	     rvec x[],rvec vold[],rvec v[],rvec vt[],rvec f[],
	     rvec buf[],t_mdatoms *mdatoms,t_nsborder *nsb,t_nrnb nrnb[],
	     t_graph *graph,t_edsamyn *edyn,t_forcerec *fr,rvec box_size)
{
  t_mdebin   *mdebin;
  int        fp_ene,fp_trn,step;
  time_t     start_t;
  real       t,lambda,t0,lam0,SAfactor;
  bool       bNS,bStopCM,bStopRot,bTYZ,bRerunMD,bNotLastFrame,bLastStep;
  tensor     force_vir,shake_vir;
  t_nrnb     mynrnb;
  char       *traj,*xtc_traj; /* normal and compressed trajectory filename */
  int        i,m;
  rvec       vcm,mu_tot;
  rvec       *xx,*vv,*ff;
  int        natoms,status;
  
  /* check if "-rerun" is used. */
  bRerunMD = optRerunMDset(nfile,fnm);
  
  /* Initial values */
  init_md(cr,&parm->ir,&t,&t0,&lambda,&lam0,&SAfactor,&mynrnb,&bTYZ,top,
	  nfile,fnm,&traj,&xtc_traj,&fp_ene,&mdebin,grps,vcm,
	  force_vir,shake_vir,mdatoms);
  
  /* Remove periodicity */  
  if (parm->ir.eBox != ebtNONE)
    do_pbc_first(log,parm,box_size,fr,graph,x);
  
  if (!parm->ir.bUncStart) 
    do_shakefirst(log,bTYZ,lambda,ener,parm,nsb,mdatoms,x,vold,buf,f,v,
		  graph,cr,&mynrnb,grps,fr,top,edyn);
    
  /* Compute initial EKin for all.. */
  calc_ke_part(TRUE,START(nsb),HOMENR(nsb),vold,v,vt,&(parm->ir.opts),
	       mdatoms,grps,&mynrnb,lambda,&ener[F_DVDLKIN]);
	       
  if (PAR(cr)) 
    global_stat(log,cr,ener,force_vir,shake_vir,
		&(parm->ir.opts),grps,&mynrnb,nrnb,vcm,mu_tot);
  clear_rvec(vcm);
  
  /* Calculate Temperature coupling parameters lambda */
  ener[F_TEMP] = sum_ekin(&(parm->ir.opts),grps,parm->ekin,bTYZ);
  tcoupl(parm->ir.btc,&(parm->ir.opts),grps,parm->ir.delta_t,SAfactor,0,
	 parm->ir.ntcmemory);
  where();
  
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
    
    /* for rerun MD always do Neighbour Searching */
    if (bRerunMD) 
      bNS = TRUE;
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

    /* Set values for invmass etc. */
    init_mdatoms(mdatoms,lambda,(step==0));

    /* Calculate total dipole moment if necessary */    
    calc_mu(nsb,x,mdatoms->chargeT,mu_tot);
    
    /* Calc forces and virial */
    clear_mat(force_vir);
    do_force(log,cr,parm,nsb,force_vir,step,&mynrnb,
	     top,grps,x,v,f,buf,mdatoms,ener,bVerbose && !PAR(cr),
	     lambda,graph,bNS,FALSE,fr);
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
    
    if (do_per_step(step,parm->ir.nstxout) || bLastStep) xx=x; else xx=NULL;
    if (do_per_step(step,parm->ir.nstvout) || bLastStep) vv=v; else vv=NULL;
    if (do_per_step(step,parm->ir.nstfout)) ff=f; else ff=NULL;
    fp_trn = write_traj(log,cr,traj,nsb,step,t,lambda,
			nrnb,nsb->natoms,xx,vv,ff,parm->box);
    where();
    
    /* for rerunMD, certain things don't have to be done */
    if (!bRerunMD) {
      if (do_per_step(step,parm->ir.nstxtcout)) {
	write_xtc_traj(log,cr,xtc_traj,nsb,mdatoms,
		       step,t,x,parm->box,parm->ir.xtcprec);
	if (bLastStep) {
	  fprintf(stderr,"Writing final coordinates.\n");
	  write_sto_conf(ftp2fn(efSTO,nfile,fnm),
			 *top->name, &(top->atoms),x,v,parm->box);
	}
	where();
      }
      
      where();
      
      clear_mat(shake_vir);
      update(nsb->natoms,START(nsb),HOMENR(nsb),step,lambda,&ener[F_DVDL],
	     &(parm->ir),FALSE,mdatoms,x,graph,
	     fr->shift_vec,f,buf,vold,v,vt,parm->pres,parm->box,
	     top,grps,shake_vir,cr,&mynrnb,bTYZ,TRUE,edyn);
#ifdef DEBUG
      pr_rvecs(log,0,"shake_vir",shake_vir,DIM);
#endif
      if (PAR(cr)) 
	accumulate_u(cr,&(parm->ir.opts),grps);
      
      where();
      /* Calculate partial Kinetic Energy (for this processor) 
       * per group!
       */
      calc_ke_part(FALSE,START(nsb),HOMENR(nsb),
		   vold,v,vt,&(parm->ir.opts),
		   mdatoms,grps,&mynrnb,
		   lambda,&ener[F_DVDL]);
      where();
      if (bStopCM)
	calc_vcm(log,HOMENR(nsb),START(nsb),mdatoms->massT,v,vcm);
    }
    
    if (PAR(cr)) {
      /* Globally (over all CPUs) sum energy, virial etc. */
      global_stat(log,cr,ener,force_vir,shake_vir,
		  &(parm->ir.opts),grps,&mynrnb,nrnb,vcm,mu_tot);
      if (fr->bTwinRange && !bNS) 
	for(i=0; (i<grps->estat.nn); i++)
	  grps->estat.ee[egLR][i] /= cr->nprocs;
    }
    else
      cp_nrnb(&(nrnb[0]),&mynrnb);
    
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
    
      /* Calculate pressure ! */
      calc_pres(parm->box,parm->ekin,parm->vir,parm->pres,
		EEL_LR(fr->eeltype) ? ener[F_LR] : 0.0);
    }
    
    /* Calculate long range corrections to pressure and energy */
    calc_ljcorr(log,parm->ir.bLJcorr,
		fr,mdatoms->nr,parm->box,parm->pres,parm->vir,ener);
    
    if (MASTER(cr))
      upd_mdebin(mdebin,mdatoms->tmass,step,ener,parm->box,shake_vir,
		 force_vir,parm->vir,parm->pres,grps,mu_tot);
    where();
    
    if ( MASTER(cr) ) {
      bool do_ene,do_log;
      
      do_ene=do_per_step(step,parm->ir.nstenergy) || bLastStep;
      do_log=do_per_step(step,parm->ir.nstlog) || bLastStep;
      print_ebin(do_ene?fp_ene:-1,do_log?log:NULL,step,t,lambda,SAfactor,
		 eprNORMAL,bCompact,mdebin,grps,&(top->atoms));
    }
    if (bVerbose)
      fflush(log);
    
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
  
  if (MASTER(cr)) {
    print_ebin(-1,log,step,t,lambda,SAfactor,
	       eprAVER,FALSE,mdebin,grps,&(top->atoms));
    print_ebin(-1,log,step,t,lambda,SAfactor,
	       eprRMS,FALSE,mdebin,grps,&(top->atoms));
    close_enx(fp_ene);
    if (parm->ir.nstxtcout)
      close_xtc_traj();
    close_trn(fp_trn);
  }
  
  return start_t;
}
