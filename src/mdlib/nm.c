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
static char *SRCID_nm_c = "$Id$";

#include "typedefs.h"
#include "vec.h"
#include "smalloc.h"
#include "statutil.h"
#include "vcm.h"
#include "mdebin.h"
#include "nrnb.h"
#include "vcm.h"
#include "calcmu.h"
#include "dummies.h"
#include "pppm.h"
#include "update.h"
#include "mdrun.h"

time_t do_nm(FILE *log,t_commrec *cr,int nfile,t_filenm fnm[],
	     bool bVerbose,bool bCompact,int stepout,
	     t_parm *parm,t_groups *grps,
	     t_topology *top,real ener[],
	     rvec x[],rvec vold[],rvec v[],rvec vt[],rvec f[],
	     rvec buf[],t_mdatoms *mdatoms,
	     t_nsborder *nsb,t_nrnb nrnb[],
	     t_graph *graph,t_edsamyn *edyn,
	     t_forcerec *fr,rvec box_size)
{
  t_mdebin   *mdebin;
  int        fp_ene,step,nre;
  time_t     start_t;
  real       t,lambda,t0,lam0;
  bool       bNS,bStopCM,bTYZ,bLR,bLJLR,bBHAM,b14,bBox;
  tensor     force_vir,shake_vir;
  t_nrnb     mynrnb;
  int        i,m;
  rvec       mu_tot;
  t_vcm      *vcm;
  rvec       *xx,*vv,*ff;
  
  /* added with respect to mdrun */
  int        idum,jdum,kdum;
  real       der_range=1.0e-6,fmax;
  rvec       *dfdx;
  snew(dfdx,top->atoms.nr);

  vv=NULL;
  ff=NULL; 

  /* end nmrun differences */

  /* Initial values */
  t0           = parm->ir.init_t;
  lam0         = parm->ir.init_lambda;
  t            = t0;
  lambda       = lam0;

  /* Check Environment variables */
  bTYZ=getenv("TYZ") != NULL;
  
  init_nrnb(&mynrnb);

  bBox = (fr->ePBC != epbcNONE);
  if (bBox) {
    calc_shifts(parm->box,box_size,fr->shift_vec,FALSE);
    fprintf(log,"Removing pbc first time\n");
    mk_mshift(log,graph,parm->box,x);
    shift_self(graph,parm->box,x);
  }

  fp_ene=-1;
  set_pot_bools(&(parm->ir),top,&bLR,&bLJLR,&bBHAM,&b14);
  mdebin=init_mdebin(fp_ene,grps,&(top->atoms),&(top->idef),bLR,bLJLR,
		     bBHAM,b14,parm->ir.efep!=efepNO,parm->ir.epc,
		     parm->ir.bDispCorr,TRICLINIC(parm->ir.compress),(parm->ir.etc==etcNOSEHOOVER),cr);

  /* Compute initial EKin for all.. */
  calc_ke_part(TRUE,FALSE,0,top->atoms.nr,
               vold,v,vt,&(parm->ir.opts),
               mdatoms,grps,&mynrnb,
	       lambda,&ener[F_DVDL]);
  
  /* Calculate Temperature coupling parameters lambda */
  ener[F_TEMP]=sum_ekin(&(parm->ir.opts),grps,parm->ekin,bTYZ);
  if(parm->ir.etc==etcBERENDSEN)
    berendsen_tcoupl(&(parm->ir.opts),grps,parm->ir.delta_t,lam0);
  else if(parm->ir.etc==etcNOSEHOOVER)
    nosehoover_tcoupl(&(parm->ir.opts),grps,parm->ir.delta_t,lam0);

  where();
  
  /* Write start time and temperature */
  start_t=print_date_and_time(log,cr->nodeid,"Started nmrun");
  if (MASTER(cr)) {
    fprintf(stderr,"starting nmrun '%s'\n%d steps.\n\n",*(top->name),
            top->atoms.nr);
  }

  vcm = init_vcm(stdlog,top,mdatoms);  
  
  /* Call do_force once to make pairlist */
  clear_mat(force_vir);
    
  bNS=TRUE;
  do_force(log,cr,parm,nsb,force_vir,0,&mynrnb,
	   top,grps,x,v,f,buf,mdatoms,ener,bVerbose && !PAR(cr),
	   lambda,graph,bNS,FALSE,fr,mu_tot,FALSE);
  bNS=FALSE;
  if (bBox)
    /* Shift back the coordinates, since we're not calling update */
    unshift_self(graph,parm->box,x);

  
  /* if forces not small, warn user */
  fmax=f_max(cr->left,cr->right,nsb->nnodes,0,top->atoms.nr,f);
  fprintf(stderr,"Maximum force:%12.5e\n",fmax);
  if (fmax > 1.0e-3) {
    fprintf(stderr,"Maximum force probably not small enough to");
    fprintf(stderr," ensure that you are in a \nenergy well. ");
    fprintf(stderr,"Be aware that negative eigenvalues may occur");
    fprintf(stderr," when the\nresulting matrix is diagonalized.\n");
  }

  /***********************************************************
   *
   *      Loop over all pairs in matrix 
   *
   *      do_force called twice. Once with positive and 
   *      once with negative displacement 
   *
   ************************************************************/

  /* fudge nr of steps to nr of atoms */

  parm->ir.nsteps=top->atoms.nr;

  for (step=0; (step<parm->ir.nsteps); step++) {

    for (idum=0; (idum<DIM); idum++) {

      x[step][idum]=x[step][idum]-der_range;
      
      clear_mat(force_vir);

      do_force(log,cr,parm,nsb,force_vir,2*(step*DIM+idum),&mynrnb,
	       top,grps,x,v,f,buf,mdatoms,ener,bVerbose && !PAR(cr),
	       lambda,graph,bNS,FALSE,fr,mu_tot,FALSE);
      if (bBox)
	/* Shift back the coordinates, since we're not calling update */
	unshift_self(graph,parm->box,x);
      
      for (jdum=0; (jdum<top->atoms.nr); jdum++) {
	for (kdum=0; (kdum<DIM); kdum++) {
	  dfdx[jdum][kdum]=f[jdum][kdum];  
	}
      }
      
      x[step][idum]=x[step][idum]+2.0*der_range;
      
      clear_mat(force_vir);
      
      do_force(log,cr,parm,nsb,force_vir,2*(step*DIM+idum)+1,&mynrnb,
	       top,grps,x,v,f,buf,mdatoms,ener,bVerbose && !PAR(cr),
	       lambda,graph,bNS,FALSE,fr,mu_tot,FALSE);
      if (bBox)
	/* Shift back the coordinates, since we're not calling update */
	unshift_self(graph,parm->box,x);
      
      for (jdum=0; (jdum<top->atoms.nr); jdum++) {
	for (kdum=0; (kdum<DIM); kdum++) {
	  dfdx[jdum][kdum]=(f[jdum][kdum]-dfdx[jdum][kdum])/2.0e-6;
	}
      }

      /* store derivatives now, diagonalization later */
      xx=dfdx;
      write_traj(log,cr,ftp2fn(efMTX,nfile,fnm),
		 nsb,step,t,lambda,nrnb,nsb->natoms,xx,vv,ff,parm->box);

      /* x is restored to original */
      x[step][idum]=x[step][idum]-der_range;
      
      if (bVerbose)
	fflush(log);
      
      /*if (MASTER(cr) && bVerbose && ((step % stepout)==0))
	print_time(stderr,start_t,step,&parm->ir);*/
    }
    /* write progress */
    if (MASTER(cr) && bVerbose) {
      fprintf(stderr,"\rFinished step %d out of %d",step+1,top->atoms.nr); 
      fflush(stderr);
    }
  }
  t=t0+step*parm->ir.delta_t;
  lambda=lam0+step*parm->ir.delta_lambda;
  
  if (MASTER(cr)) {
    print_ebin(-1,FALSE,FALSE,log,step,t,eprAVER,
	       FALSE,mdebin,&(top->atoms));
    print_ebin(-1,FALSE,FALSE,log,step,t,eprRMS,
	       FALSE,mdebin,&(top->atoms));
  }
  
  /* Construct dummy particles, for last output frame */
  construct_dummies(log,x,&mynrnb,parm->ir.delta_t,v,&top->idef);
    
  /*free_nslist(log);*/
  
  return start_t;
}
