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
static char *SRCID_congrad_c = "$Id$";
#include <string.h>
#include <time.h>
#include <math.h>
#include "sysstuff.h"
#include "string2.h"
#include "network.h"
#include "confio.h"
#include "copyrite.h"
#include "vcm.h"
#include "smalloc.h"
#include "nrnb.h"
#include "main.h"
#include "pbc.h"
#include "force.h"
#include "macros.h"
#include "random.h"
#include "names.h"
#include "fatal.h"
#include "txtdump.h"
#include "typedefs.h"
#include "update.h"
#include "random.h"
#include "vec.h"
#include "statutil.h"
#include "tgroup.h"
#include "mdebin.h"
#include "dummies.h"
#include "mdrun.h"

static void sp_header(FILE *out,real epot,real ftol)
{
  fprintf(out,"Conjugate gradients:\n");
  fprintf(out,"   Tolerance         = %12.5e\n",ftol);
  fprintf(out,"   Starting Energy   = %20.15e\n",epot);
}

time_t do_cg(FILE *log,int nfile,t_filenm fnm[],
	     t_parm *parm,t_topology *top,
	     t_groups *grps,t_nsborder *nsb,
	     rvec x[],rvec grad[],rvec buf[],t_mdatoms *mdatoms,
	     tensor ekin,real ener[],t_nrnb nrnb[],
	     bool bVerbose,bool bDummies,t_comm_dummies *dummycomm,
	     t_commrec *cr,t_graph *graph,
	     t_forcerec *fr,rvec box_size)
{
  static char *CG="Conjugate Gradients";
  real   step0,lambda,ftol,fmax,testf,zet,w,smin;
  rvec   *p,*f,*xprime,*xx,*ff;
  real   EpotA=0.0,EpotB=0.0,a=0.0,b,beta=0.0;
  double gpa,gpb;
  real   fnorm,pnorm,fnorm_old;
  t_vcm      *vcm;
  t_mdebin   *mdebin;
  t_nrnb mynrnb;
  bool   bNS=TRUE,bDone,bLR,bLJLR,bBHAM,b14,bRand,brerun;
  rvec   mu_tot;
  time_t start_t;
  tensor force_vir,shake_vir,pme_vir;
  int    number_steps,naccept=0,nstcg=parm->ir.nstcgsteep;
  int    fp_ene,count=0;
  int    i,m,start,end,niti,gf;
  /* not used */
  real   terminate=0;

  /* Initiate some variables */
  if (parm->ir.efep != efepNO)
    lambda       = parm->ir.init_lambda;
  else 
    lambda = 0.0;

  init_nrnb(&mynrnb);
    
  clear_rvec(mu_tot);
  calc_shifts(parm->box,box_size,fr->shift_vec,FALSE);
  
  vcm = init_vcm(stdlog,top,mdatoms,
		 START(nsb),HOMENR(nsb),parm->ir.nstcomm);
    
  /* Print to log file */
  start_t=print_date_and_time(log,cr->nodeid,"Started Conjugate Gradients");
  
  /* p is the search direction, f the force, xprime the new positions */
  snew(p,nsb->natoms);
  snew(f,nsb->natoms);
  snew(xprime,nsb->natoms);

  start = nsb->index[cr->nodeid];
  end   = nsb->homenr[cr->nodeid] + start;

  /* Set some booleans for the epot routines */
  set_pot_bools(&(parm->ir),top,&bLR,&bLJLR,&bBHAM,&b14);

  /* Open the energy file */  
  if (MASTER(cr))
    fp_ene=open_enx(ftp2fn(efENX,nfile,fnm),"w");
  else
    fp_ene=-1;
    
  /* Init bin for energy stuff */
  mdebin=init_mdebin(fp_ene,grps,&(top->atoms),&(top->idef),
		     bLR,bLJLR,bBHAM,b14,parm->ir.efep!=efepNO,parm->ir.epc,
		     parm->ir.eDispCorr,TRICLINIC(parm->ir.compress),(parm->ir.etc==etcNOSEHOOVER),cr); 

  /* Clear some matrix variables */
  clear_mat(force_vir);
  clear_mat(shake_vir);
  
  /* Set variables for stepsize (in nm). This is the largest 
   * step that we are going to make in any direction.
   */
  step0=parm->ir.em_stepsize;
  
  /* Tolerance for convergence */
  ftol=parm->ir.em_tol;
  
  /* Max number of steps */
  number_steps=parm->ir.nsteps;

  if (fr->ePBC != epbcNONE)
    /* Remove periodicity */
    do_pbc_first(log,parm,box_size,fr,graph,x);
  
  if (bDummies) {
    /* Molecules always whole, but I'm not sure whether
     * the periodicity and shift are guaranteed to be consistent
     * between different nodes when running e.g. polymers in
     * parallel. In this special case we thus unshift/shift, but
     * only when necessary. This is to make sure the coordinates
     * we move don't end up a box away...
     */
    if(dummycomm) {
      unshift_self(graph,parm->box,x);
      move_construct_x(dummycomm,x,cr);
      shift_self(graph,parm->box,x);
    }

    construct_dummies(log,x,&(nrnb[cr->nodeid]),1,NULL,&top->idef);

    if(dummycomm) {
      unshift_self(graph,parm->box,x);
      move_dummy_xv(dummycomm,x,NULL,cr);
      shift_self(graph,parm->box,x); /* maybe not necessary */
    }
  }

  /* Call the force routine and some auxiliary (neighboursearching etc.) */
  /* do_force always puts the charge groups in the box and shifts again
   * We do not unshift, so molecules are always whole in congrad.c
   */
  do_force(log,cr,parm,nsb,force_vir,pme_vir,0,&(nrnb[cr->nodeid]),top,grps,
	   x,buf,f,buf,mdatoms,ener,bVerbose && !(PAR(cr)),
	   lambda,graph,bNS,FALSE,fr,mu_tot,FALSE);
  where();

  /* Spread the force on dummy particle to the other particles... */
  if(bDummies) {
    /* We only move forces here, and they are independent of shifts */
    if(dummycomm)
      move_dummy_f(dummycomm,f,cr);

    spread_dummy_f(log,x,f,&(nrnb[cr->nodeid]),&top->idef);
  
    if(dummycomm)
      move_construct_f(dummycomm,f,cr);
  }

  /* Sum the potential energy terms from group contributions */
  sum_epot(&(parm->ir.opts),grps,ener);
  where();
  
  /* Clear var. */
  clear_mat(force_vir);
  where();
  
  /* Communicat energies etc. */
  if (PAR(cr)) 
    global_stat(log,cr,ener,force_vir,shake_vir,
		&(parm->ir.opts),grps,&mynrnb,nrnb,vcm,&terminate);
  where();
  
  if (MASTER(cr)) {
    /* Copy stuff to the energy bin for easy printing etc. */
    upd_mdebin(mdebin,NULL,mdatoms->tmass,count,(real)count,
	       ener,parm->box,shake_vir,
	       force_vir,parm->vir,parm->pres,grps,mu_tot,
	       (parm->ir.etc==etcNOSEHOOVER));
    
    print_ebin_header(log,count,count,lambda,0.0);
    print_ebin(fp_ene,TRUE,FALSE,log,count,count,eprNORMAL,
	       TRUE,mdebin,&(top->atoms));
  }
  where();
  
  /* This is the starting energy */
  EpotA=ener[F_EPOT];

  if (MASTER(cr))
    sp_header(stderr,EpotA,ftol);
  sp_header(log,EpotA,ftol);

  /* normalising step size, this saves a few hundred steps in the
   * beginning of the run.
   */
  fnorm = f_norm(cr,&(parm->ir.opts),mdatoms,start,end,f);
  fnorm_old=fnorm;

  /* Print stepsize */
  if (MASTER(cr)) {
    fprintf(stderr,"   F-Norm            = %12.5e\n",fnorm);
    fprintf(stderr,"\n");
  }
  
  /* Here starts the loop, count is the counter for the number of steps
   * bDone is a BOOLEAN variable, which will be set TRUE when
   * the minimization has converged.
   */
  for(count=1,bDone=FALSE; 
      !(bDone || ((number_steps > 0) && (count==number_steps))); count++) {

    /* start conjugate gradient, determine search interval a,b */
    gpa=0.0;
    for(i=start; i<end; i++) {
      gf = mdatoms->cFREEZE[i];
      for(m=0; m<DIM; m++) 
	if (!parm->ir.opts.nFreeze[gf][m]) {
	  p[i][m] = f[i][m] + beta*p[i][m];
	  gpa     = gpa - p[i][m]*f[i][m];
	}
    }
    gmx_sumd(1,&gpa,cr);
    pnorm=f_norm(cr,&(parm->ir.opts),mdatoms,start,end,p);

    a    = 0.0;
    b    = step0/pnorm;
    niti = 0;

    /* search a,b iteratively, if necessary */
    brerun=TRUE;
    while (brerun) {
      for (i=start; i<end; i++) {
	for(m=0; m<DIM; m++) {
	  xprime[i][m] = x[i][m] + b*p[i][m];
	}
      }
      bNS = (parm->ir.nstlist > 0);

      if (bDummies) {
	/* Molecules always whole, but I'm not sure whether
	 * the periodicity and shift are guaranteed to be consistent
	 * between different nodes when running e.g. polymers in
	 * parallel. In this special case we thus unshift/shift, but
	 * only when necessary. This is to make sure the coordinates
	 * we move don't end up a box away...
	 */
	if(dummycomm) {
	  unshift_self(graph,parm->box,xprime);
	  move_construct_x(dummycomm,xprime,cr);
	  shift_self(graph,parm->box,xprime);
	}
 
	construct_dummies(log,xprime,&(nrnb[cr->nodeid]),1,NULL,&top->idef);

	if(dummycomm) {
	  unshift_self(graph,parm->box,xprime);
	  move_dummy_xv(dummycomm,xprime,NULL,cr);
	  shift_self(graph,parm->box,xprime); /* maybe not necessary */
	}
      }
      
      /* Calc force & energy on new trial position  */
      /* do_force always puts the charge groups in the box and shifts again
       * We do not unshift, so molecules are always whole in congrad.c
       */
      do_force(log,cr,parm,nsb,force_vir,pme_vir,
	       count,&(nrnb[cr->nodeid]),top,grps,xprime,buf,f,
	       buf,mdatoms,ener,bVerbose && !(PAR(cr)),
	       lambda,graph,bNS,FALSE,fr,mu_tot,FALSE);
      
      /* Spread the force on dummy particle to the other particles... */
      if(bDummies) {
	/* We only move forces here, and they are independent of shifts */
	if(dummycomm)
	  move_dummy_f(dummycomm,f,cr);

	spread_dummy_f(log,xprime,f,&(nrnb[cr->nodeid]),&top->idef); 

	if(dummycomm)
	  move_construct_f(dummycomm,f,cr);
      }

      bNS = (parm->ir.nstlist > 0);

      gpb=0.0;
      for(i=start; i<end; i++)
	for(m=0; m<DIM; m++) 
	  gpb -= p[i][m]*f[i][m];

      gmx_sumd(1,&gpb,cr);

      /* Sum the potential energy terms from group contributions */
      sum_epot(&(parm->ir.opts),grps,ener);
      
      /* Clear stuff again */
      clear_mat(force_vir);
      clear_mat(shake_vir);
      
      /* Communicate stuff when parallel */
      if (PAR(cr)) 
	global_stat(log,cr,ener,force_vir,shake_vir,
		    &(parm->ir.opts),grps,&mynrnb,nrnb,vcm,&terminate);

      EpotB=ener[F_EPOT];
      
      if ((gpb >= 0.0) || (EpotB >= EpotA))
	brerun=FALSE;
      else {
	a     = b;
	EpotA = EpotB;
	gpa   = gpb;
	b    += b;
      }
      niti++;
    }
    /* End of while loop searching for a and b */
    
    /* find stepsize smin in interval a-b */
    zet = 3.0 * (EpotA-EpotB) / (b-a) + gpa + gpb;
    w   = zet*zet - gpa*gpb;
    if (w < 0.0) {
      fprintf(stderr,"Negative w: %20.12e\n",w);
      fprintf(stderr,"z= %20.12e\n",zet);
      fprintf(stderr,"gpa= %20.12e, gpb= %20.12e\n",gpa,gpb);
      fprintf(stderr,"a= %20.12e, b= %20.12e\n",a,b);
      fprintf(stderr,"EpotA= %20.12e, EpotB= %20.12e\n",EpotA,EpotB);
      fprintf(stderr,"Negative number for sqrt encountered (%f)\n",w);
      fprintf(stderr,"Terminating minimization\n");
      break;
    }      
    w    = sqrt(w);
    smin = b - ((gpb+w-zet)*(b-a))/((gpb-gpa)+2.0*w);

    /* new positions */
    for (i=start; i<end; i++)
      for(m=0; m<DIM; m++) 
	xprime[i][m] = x[i][m] + smin*p[i][m];

    if (bDummies) {
      /* Molecules always whole, but I'm not sure whether
       * the periodicity and shift are guaranteed to be consistent
       * between different nodes when running e.g. polymers in
       * parallel. In this special case we thus unshift/shift, but
       * only when necessary. This is to make sure the coordinates
       * we move don't end up a box away...
       */
      if(dummycomm) {
        unshift_self(graph,parm->box,xprime);
        move_construct_x(dummycomm,xprime,cr);
        shift_self(graph,parm->box,xprime);
      }

      construct_dummies(log,xprime,&(nrnb[cr->nodeid]),1,NULL,&top->idef);

      if(dummycomm) {
        unshift_self(graph,parm->box,xprime);
        move_dummy_xv(dummycomm,xprime,NULL,cr);
        shift_self(graph,parm->box,xprime); /* maybe not necessary */
      }
    }
    
    /* new energy, forces */
    /* do_force always puts the charge groups in the box and shifts again
     * We do not unshift, so molecules are always whole in congrad.c
     */
    do_force(log,cr,parm,nsb,force_vir,pme_vir,
	     count,&(nrnb[cr->nodeid]),top,grps,xprime,buf,f,
	     buf,mdatoms,ener,bVerbose && !(PAR(cr)),
	     lambda,graph,bNS,FALSE,fr,mu_tot,FALSE);
    
    /* Spread the force on dummy particle to the other particles... */
    if(bDummies) {
      /* We only move forces here, and they are independent of shifts */
      if(dummycomm)
        move_dummy_f(dummycomm,f,cr);

      spread_dummy_f(log,xprime,f,&(nrnb[cr->nodeid]),&top->idef); 

      if(dummycomm)
        move_construct_f(dummycomm,f,cr);
    }


    /* Sum the potential energy terms from group contributions */
    sum_epot(&(parm->ir.opts),grps,ener); 
    fnorm = f_norm(cr,&(parm->ir.opts),mdatoms,start,end,f);
    
    /* Clear stuff again */
    clear_mat(force_vir);
    clear_mat(shake_vir);
      
    /* Communicate stuff when parallel */
    if (PAR(cr)) 
      global_stat(log,cr,ener,force_vir,shake_vir,
		  &(parm->ir.opts),grps,&mynrnb,nrnb,vcm,&terminate);

    EpotA=ener[F_EPOT];

    /* new search direction */
    /* beta = 0 means steepest descents */
    if (nstcg && ((count % nstcg)==0)) 
      beta = 0.0;
    else
      beta=fnorm*fnorm/(fnorm_old*fnorm_old);

    /* update x, fnorm_old */
    for (i=start; i<end; i++)
      copy_rvec(xprime[i],x[i]);
    fnorm_old=fnorm;

    /* Test whether the convergence criterion is met */
    fmax=f_max(cr->left,cr->right,nsb->nnodes,&(parm->ir.opts),mdatoms,
	       start,end,f);

    /* Print it if necessary */
    if (bVerbose && MASTER(cr)) {
      fprintf(stderr,"\rStep %d, E-Pot = %16.10e, F-max = %12.5e\n",
	      count,EpotA,fmax);
      /* Store the new (lower) energies */
      upd_mdebin(mdebin,NULL,mdatoms->tmass,count,(real)count,
		 ener,parm->box,shake_vir,
		 force_vir,parm->vir,parm->pres,grps,mu_tot,(parm->ir.etc==etcNOSEHOOVER));
      /* Print the energies allways when we should be verbose */
      print_ebin_header(log,count,count,lambda,0.0);
      print_ebin(fp_ene,TRUE,FALSE,log,count,count,eprNORMAL,
		 TRUE,mdebin,&(top->atoms));
    }
    
    /* Stop when the maximum force lies below tolerance */
    bDone=(fmax < ftol);

    /*if (MASTER(cr))
      fprintf(stderr,"\n");*/
	
  } /* End of the loop */

  /* Print some shit... */
  if (MASTER(cr))
    fprintf(stderr,"\nwriting lowest energy coordinates.\n");
  xx=x;
  ff=f;
  write_traj(log,cr,ftp2fn(efTRN,nfile,fnm),
	     nsb,count,(real) count,
	     lambda,nrnb,nsb->natoms,xx,NULL,ff,parm->box);
  if (MASTER(cr))
    write_sto_conf(ftp2fn(efSTO,nfile,fnm),
		   *top->name, &(top->atoms),xx,NULL,parm->box);

  fmax=f_max(cr->left,cr->right,nsb->nnodes,&(parm->ir.opts),mdatoms,
	     start,end,f);
  if (MASTER(cr)) {
    fprintf(stderr,"Maximum force: %12.5e\n",fmax);
    if (bDone) {
      fprintf(stderr,"\n%s converged to %8.6f in %d steps\n",CG,ftol,count-1);
      fprintf(log,"%s converged to %8.6f \n",CG,ftol);
    }
    else {
      fprintf(stderr,"\n%s did not converge in %d steps\n",CG,number_steps);
      fprintf(log,"%s did not converge in %d steps\n",CG,number_steps);
    }
    fprintf(stderr,"  Function value at minimum = %12.4e\n",EpotA);
    fprintf(log,"  Function value at minimum = %12.4e\n",EpotA);

    close_enx(fp_ene);
  }
  
  /* To print the actual number of steps we needed somewhere */
  parm->ir.nsteps=count;

  return start_t;
} /* That's all folks */

