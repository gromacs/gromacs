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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
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
  tensor     force_vir,shake_vir,pme_vir;
  t_nrnb     mynrnb;
  int        i,m;
  rvec       mu_tot;
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
		     parm->ir.eDispCorr,TRICLINIC(parm->ir.compress),(parm->ir.etc==etcNOSEHOOVER),cr);

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

  /* Call do_force once to make pairlist */
  clear_mat(force_vir);
    
  bNS=TRUE;
  do_force(log,cr,parm,nsb,force_vir,pme_vir,0,&mynrnb,
	   top,grps,x,v,f,buf,mdatoms,ener,bVerbose && !PAR(cr),
	   lambda,graph,bNS,FALSE,fr,mu_tot,FALSE);
  bNS=FALSE;
  if (bBox)
    /* Shift back the coordinates, since we're not calling update */
    unshift_self(graph,parm->box,x);

  
  /* if forces not small, warn user */
  fmax=f_max(cr->left,cr->right,nsb->nnodes,&(parm->ir.opts),mdatoms,
	     0,top->atoms.nr,f);
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

      do_force(log,cr,parm,nsb,force_vir,pme_vir,2*(step*DIM+idum),&mynrnb,
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
      
      do_force(log,cr,parm,nsb,force_vir,pme_vir,2*(step*DIM+idum)+1,&mynrnb,
	       top,grps,x,v,f,buf,mdatoms,ener,bVerbose && !PAR(cr),
	       lambda,graph,bNS,FALSE,fr,mu_tot,FALSE);
      if (bBox)
	/* Shift back the coordinates, since we're not calling update */
	unshift_self(graph,parm->box,x);
      
      for (jdum=0; (jdum<top->atoms.nr); jdum++) {
	for (kdum=0; (kdum<DIM); kdum++) {
	  dfdx[jdum][kdum]=(f[jdum][kdum]-dfdx[jdum][kdum])/(2*der_range);
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
  /* NB: I have not included the communication for parallel
   * dummies here, since the rest of nm doesn't seem to
   * be parallelized. Be sure to copy the correct code from
   * e.g. md.c or steep.c if you make nm parallel!
   */
  construct_dummies(log,x,&mynrnb,parm->ir.delta_t,v,&top->idef);
    
  /*free_nslist(log);*/
  
  return start_t;
}



