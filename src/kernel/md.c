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

/*
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "sysstuff.h"
#include "string2.h"
#include "led.h"
#include "nrnb.h"
#include "network.h"
#include "confio.h"
#include "copyrite.h"
#include "smalloc.h"
#include "main.h"
#include "pbc.h"
#include "force.h"
#include "macros.h"
#include "names.h"
#include "fatal.h"
#include "txtdump.h"
#include "futil.h"
#include "typedefs.h"
#include "update.h"
#include "random.h"
#include "vec.h"
#include "filenm.h"
#include "statutil.h"
#include "tgroup.h"
#include "mdrun.h"
#include "vcm.h"
#include "ebin.h"
#include "mdebin.h"
#include "dummies.h"
#include "mdrun.h"
#include "physics.h"
#include "pppm.h"
#include "edsam.h"
#include "calcmu.h"
*/
#include "typedefs.h"
#include "vec.h"
#include "statutil.h"
#include "vcm.h"
#include "mdebin.h"
#include "nrnb.h"
#include "calcmu.h"
#include "dummies.h"
#include "pppm.h"
#include "update.h"
#include "mdrun.h"

time_t do_md(FILE *log,t_commrec *cr,int nfile,t_filenm fnm[],
	     bool bVerbose,bool bCompact,bool bDummies,int stepout,
	     t_parm *parm,t_groups *grps,
	     t_topology *top,real ener[],
	     rvec x[],rvec vold[],rvec v[],rvec vt[],rvec f[],
	     rvec buf[],t_mdatoms *mdatoms,
	     t_nsborder *nsb,t_nrnb nrnb[],
	     t_graph *graph,t_edsamyn *edyn,
	     t_forcerec *fr,rvec box_size)
{
  t_mdebin   *mdebin;
  int        fp_ene,step;
  time_t     start_t;
  real       t,lambda,t0,lam0,SAfactor;
  bool       bNS,bStopCM,bStopRot,bTYZ,bLR,bBHAM,b14,bRerunMD,bNotLastFrame;
  tensor     force_vir,shake_vir;
  t_nrnb     mynrnb;
  char       *traj,*enerfile;
  char       *xtc_traj; /* compressed trajectory filename */
  int        i,m;
  rvec       vcm,mu_tot;
  rvec       *xx,*vv,*ff;
  int        natoms,status;
  
  /* check if "-rerun" is used. */
  bRerunMD = optRerunMDset(nfile,fnm);
  
  /* Initial values */
  t0           = parm->ir.init_t;
  if (parm->ir.bPert) {
    lam0         = parm->ir.init_lambda;
    lambda       = lam0;
  }
  else {
    lam0   = 0.0;
    lambda = 0.0;
  } 
  if (parm->ir.bSimAnn) {
    SAfactor = 1.0  - t0/parm->ir.zero_temp_time;
    if (SAfactor < 0) 
      SAfactor = 0;
  } else
    SAfactor     = 1.0;
  
  /* Check Environment variables */
  bTYZ=getenv("TYZ") != NULL;
  
  init_nrnb(&mynrnb);
  
  calc_shifts(parm->box,box_size,fr->shift_vec,FALSE);
  
  fprintf(log,"Removing pbc first time\n");
  mk_mshift(log,graph,parm->box,x);
  shift_self(graph,fr->shift_vec,x);
  fprintf(log,"Done rmpbc\n");

  traj     = ftp2fn(efTRN,nfile,fnm);
  xtc_traj = ftp2fn(efXTC,nfile,fnm);
  where();
  
  bLR      = (parm->ir.rlong > parm->ir.rshort);
  bBHAM    = (top->idef.functype[0]==F_BHAM);
  b14      = (top->idef.il[F_LJ14].nr > 0);
  
  if (MASTER(cr)) 
    fp_ene=open_enx(ftp2fn(efENX,nfile,fnm),"w");
  else 
    fp_ene = -1;
  mdebin=init_mdebin(fp_ene,grps,&(top->atoms),bLR,bBHAM,b14);
  where();
  
  clear_rvec(vcm);
  
  /* Set initial values for invmass etc. */
  init_mdatoms(mdatoms,lambda,TRUE);
  where();
  
  if (!parm->ir.bUncStart) 
    do_shakefirst(log,bTYZ,lambda,ener,parm,nsb,mdatoms,x,vold,buf,f,v,
		  graph,cr,&mynrnb,grps,fr,top,edyn);
  where();
    
  /* Compute initial EKin for all.. */
  clear_mat(force_vir);
  clear_mat(shake_vir);
  calc_ke_part(TRUE,START(nsb),HOMENR(nsb),
	       vold,v,vt,&(parm->ir.opts),
	       mdatoms,grps,&mynrnb,
	       lambda,&ener[F_DVDLKIN]);
  if (PAR(cr)) 
    global_stat(log,cr,ener,force_vir,shake_vir,
		&(parm->ir.opts),grps,&mynrnb,nrnb,vcm,mu_tot);
  clear_rvec(vcm);
  
  /* Calculate Temperature coupling parameters lambda */
  ener[F_TEMP]=sum_ekin(&(parm->ir.opts),grps,parm->ekin,bTYZ);
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
    for(i=0;(i<natoms);i++) {
      v[i][XX]=0;
      v[i][YY]=0;
      v[i][ZZ]=0;
    }
    bNotLastFrame = (natoms!=0);
  } 

  
  /* loop over MD steps or if rerunMD to end of 
   * input trajectory */
  for(step=0; ((!bRerunMD && (step<parm->ir.nsteps)) || 
	       (bRerunMD && bNotLastFrame)); step++) {
    
    /* for rerun MD always do Neighbour Searching */
    if (bRerunMD) 
      bNS = TRUE;
    else {
      /* Stop Center of Mass motion */
      if (parm->ir.nstcomm == 0) {
	bStopCM=FALSE;
	bStopRot=FALSE;
      } else if (parm->ir.nstcomm > 0) {
	bStopCM=do_per_step(step,parm->ir.nstcomm);
	bStopRot=FALSE;
      } else {
	bStopCM=FALSE;
	bStopRot=do_per_step(step,-parm->ir.nstcomm);
      }
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
    init_mdatoms(mdatoms,lambda,FALSE);

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
    
    if (do_per_step(step,parm->ir.nstxout)) xx=x; else xx=NULL;
    if (do_per_step(step,parm->ir.nstvout)) vv=v; else vv=NULL;
    if (do_per_step(step,parm->ir.nstfout)) ff=f; else ff=NULL;
    
    write_traj(log,cr,traj,
	       nsb,step,t,lambda,nrnb,nsb->natoms,xx,vv,ff,parm->box);
    where();
    
    /* for rerunMD, certain things don't have to be done */
    if (!bRerunMD) {
      if (do_per_step(step,parm->ir.nstxtcout)) {
	write_xtc_traj(log,cr,xtc_traj,nsb,mdatoms,
		       step,t,x,parm->box,parm->ir.xtcprec);
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
	do_stopcm(log,HOMENR(nsb),START(nsb),v,vcm,mdatoms->tmass);
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
    
    upd_mdebin(mdebin,mdatoms->tmass,step,ener,parm->box,shake_vir,
	       force_vir,parm->vir,parm->pres,grps,mu_tot);
    
    where();
    
    if ( MASTER(cr) && (do_per_step(step,parm->ir.nstprint)) ) {
      print_ebin(fp_ene,log,step,t,lambda,SAfactor,
		 eprNORMAL,bCompact,mdebin,grps,&(top->atoms));
    }
    if (bVerbose)
      fflush(log);
    
    /* Time for performance */
    if ((step % 10) == 0)
      update_time();
      
    /* Remaining runtime */
    if (MASTER(cr) && bVerbose && ((step % stepout)==0))
      print_time(stderr,start_t,step,&parm->ir);
    
    /* if rerunMD read next frame from input trajectory */
    if (bRerunMD) 
      bNotLastFrame = read_next_x(status,&t,natoms,x,parm->box);
  }
  /* End of main MD loop */
  
  if (MASTER(cr)) {
    /* if RerunMD don't print energies of last frame again */
    if ( (parm->ir.nstprint > 1) && !bRerunMD )
      print_ebin(fp_ene,log,step-1,t,lambda,SAfactor,
		 eprNORMAL,bCompact,mdebin,grps,&(top->atoms));
    
    print_ebin(-1,log,step,t,lambda,SAfactor,
	       eprAVER,FALSE,mdebin,grps,&(top->atoms));
    print_ebin(-1,log,step,t,lambda,SAfactor,
	       eprRMS,FALSE,mdebin,grps,&(top->atoms));
    close_enx(fp_ene);
  }
  
  /* Construct dummy particles, for last output frame */
  construct_dummies(log,x,&mynrnb,parm->ir.delta_t,v,&top->idef);
    
  /*free_nslist(log);*/
  
  return start_t;
}
