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
static char *SRCID_steep_c = "$Id$";
#include <string.h>
#include <time.h>
#include <math.h>
#include "sysstuff.h"
#include "string2.h"
#include "network.h"
#include "confio.h"
#include "copyrite.h"
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
#include "random.h"
#include "vec.h"
#include "enxio.h"
#include "tgroup.h"
#include "mdebin.h"
#include "mdrun.h"
#include "pppm.h"
#include "dummies.h"
#include "constr.h"

real f_max(int left,int right,int nnodes,
	   t_grpopts *opts,t_mdatoms *mdatoms,
	   int start,int end,rvec grad[])
{
  real fmax2,fmax2_0,fam;
  int  i,m,gf;

  /* This routine finds the largest force and returns it.
   * On parallel machines the global max is taken.
   */
  fmax2 = 0;
  
  for(i=start; i<end; i++) {
    gf = mdatoms->cFREEZE[i];
    fam = 0;
    for(m=0; m<DIM; m++)
      if (!opts->nFreeze[gf][m])
	fam += sqr(grad[i][m]);
    fmax2 = max(fmax2,fam);
  }

  if (nnodes > 1) {
    for(i=0; i<nnodes-1; i++) {
      gmx_tx(left,(void *)&fmax2,sizeof(fmax2));
      gmx_rx(right,(void *)&fmax2_0,sizeof(fmax2_0));
      gmx_wait(left,right);
      if (fmax2_0 > fmax2)
	fmax2 = fmax2_0;
    }
  }

  return sqrt(fmax2);
}

real f_norm(int left,int right,int nnodes,
	    int start,int end,rvec grad[])
{
  real fnorm;
  int  i,m;

  /* This routine finds the norm of the force
   * and returns it. Parallel machines not supported.
   */
  fnorm = 0;
  
  for(i=start; (i<end); i++) 
    for(m=0; (m<DIM); m++) { 
      fnorm=fnorm+grad[i][m]*grad[i][m]; 
    } 
  fnorm=sqrt(fnorm); 
  
  if (nnodes == 1) 
    return fnorm; 
  
  fatal_error(0,"This version of Steepest Descents cannot be run in parallel"); 
  return fnorm; 
} 

static void do_step(int start,int end,rvec x[],rvec f[], 
 		    bool bRandom,real stepsize) 
{ 
  static int seed=1993; 
  int    i,m; 
  real   r; 
  
  if (bRandom) { 
    fprintf(stderr,"\rRandom step\n"); 
    for(i=start; (i<end); i++) { 
      for(m=0; (m<DIM); m++) { 
 	r=rando(&seed); 
 	x[i][m]=x[i][m]+stepsize*r; 
      } 
    } 
  } 
  else { 
    /* New positions to try  */
    for(i=start; (i<end);i++) 
      for(m=0; (m<DIM); m++) 
 	x[i][m] = x[i][m]+stepsize*f[i][m]; 
  } 
} 

time_t do_steep(FILE *log,int nfile,t_filenm fnm[], 
 		t_parm *parm,t_topology *top, 
 		t_groups *grps,t_nsborder *nsb, 
 		rvec x[],rvec grad[],rvec buf[],t_mdatoms *mdatoms, 
 		tensor ekin,real ener[],t_nrnb nrnb[], 
 		bool bVerbose,bool bDummies, t_comm_dummies *dummycomm,
		t_commrec *cr,t_graph *graph,
		t_forcerec *fr,rvec box_size) 
{ 
  static char *SD="STEEPEST DESCENTS",sbuf[STRLEN]; 
  real   stepsize,constepsize,lambda,ftol,fmax; 
  rvec   *pos[2],*force[2],*xcf=NULL; 
  rvec   *xx,*ff; 
  real   Fmax[2]; 
  real   Epot[2]; 
  real   ustep,dvdlambda;
  t_vcm      *vcm;
  int        fp_ene; 
  t_mdebin   *mdebin; 
  t_nrnb     mynrnb; 
  t_inputrec *ir;
  bool   bDone,bAbort,bLR,bLJLR,bBHAM,b14; 
  time_t start_t; 
  tensor force_vir,shake_vir,pme_vir; 
  rvec   mu_tot;
  int    nsteps;
  int    count=0; 
  int    i,m,start,end,gf; 
  int    Min=0; 
  int    steps_accepted=0; 
  bool   bConstrain;
  /* not used */
  real   terminate=0;
#define  TRY (1-Min)
  
  /* Initiate some variables  */
  if (parm->ir.efep != efepNO)
    lambda       = parm->ir.init_lambda;
  else 
    lambda = 0.0;
  ir = &(parm->ir);

  init_nrnb(&mynrnb);

  clear_rvec(mu_tot);
  calc_shifts(parm->box,box_size,fr->shift_vec,FALSE);
   
  /* Print to log file  */
  start_t=print_date_and_time(log,cr->nodeid,"Started Steepest Descents"); 
  
  /* We need two coordinate arrays and two force arrays  */
  for(i=0; (i<2); i++) { 
    snew(pos[i],nsb->natoms); 
    snew(force[i],nsb->natoms); 
  } 

  start=nsb->index[cr->nodeid]; 
  end=nsb->homenr[cr->nodeid]+start; 
  
  /* Open the enrgy file */   
  if (MASTER(cr)) 
    fp_ene=open_enx(ftp2fn(efENX,nfile,fnm),"w"); 
  else 
    fp_ene=-1; 
  
  /* Set some booleans for the epot routines  */
  set_pot_bools(&(parm->ir),top,&bLR,&bLJLR,&bBHAM,&b14);
  
  vcm = init_vcm(stdlog,top,mdatoms,
		 START(nsb),HOMENR(nsb),parm->ir.nstcomm);
  
  /* Init bin for energy stuff  */
  mdebin=init_mdebin(fp_ene,grps,&(top->atoms),&(top->idef),bLR,bLJLR,
		     bBHAM,b14,ir->efep!=efepNO,ir->epc,
		     ir->eDispCorr,TRICLINIC(parm->ir.compress),(parm->ir.etc==etcNOSEHOOVER),cr); 
  
  /* Clear some matrix variables  */
  clear_mat(force_vir); 
  clear_mat(shake_vir); 

  /* Initiate constraint stuff */
  bConstrain=init_constraints(stdlog,top,&(parm->ir),mdatoms,
			      start,end,FALSE,cr);
  
  if (bConstrain)
    snew(xcf,nsb->natoms); 

  if (fr->ePBC != epbcNONE)
    /* Remove periodicity */
    do_pbc_first(log,parm,box_size,fr,graph,x);

  /* Copy coord vectors to our temp array  */
  for(i=0; (i<nsb->natoms); i++) { 
    copy_rvec(x[i],pos[Min][i]); 
    copy_rvec(x[i],pos[TRY][i]); 
  } 
    
  /* Set variables for stepsize (in nm). This is the largest  
   * step that we are going to make in any direction. 
   */
  ustep    = ir->em_stepsize; 
  stepsize = 0;

  /* Tolerance for conversion  */
  ftol=ir->em_tol; 
  
  /* Max number of steps  */
  nsteps=ir->nsteps; 
  
  if (MASTER(cr)) { 
    /* Print to the screen  */
    start_t=print_date_and_time(log,cr->nodeid,"Started EM"); 
    sprintf(sbuf,"%s\n   Tolerance         = %12.5g\n",SD,ftol); 
    fprintf(stderr,sbuf);
    fprintf(log,sbuf);
  } 
    
  /**** HERE STARTS THE LOOP ****
   * count is the counter for the number of steps 
   * bDone will be TRUE when the minimization has converged
   * bAbort will be TRUE when nsteps steps have been performed or when
   * the stepsize becomes smaller than is reasonable for machine precision
   */
  count=0;
  bDone=FALSE;
  bAbort=FALSE;
  while( !bDone && !bAbort ) {
    bAbort = (nsteps > 0) && (count==nsteps);

    /* set new coordinates, except for first step */
    if (count>0)
      for(i=start; i<end; i++) {
	gf = mdatoms->cFREEZE[i];
	for(m=0; m<DIM; m++) 
	  if (ir->opts.nFreeze[gf][m])
	    pos[TRY][i][m] = pos[Min][i][m];
	  else
	    pos[TRY][i][m] = pos[Min][i][m] + stepsize * force[Min][i][m]; 
      }

    if (bConstrain) {
      dvdlambda=0;
      constrain(stdlog,top,&(parm->ir),count,mdatoms,start,end,
		pos[Min],pos[TRY],NULL,parm->box,lambda,&dvdlambda,nrnb,
		TRUE);
    }

    if (bDummies) {
      /* Molecules always whole, but I'm not sure whether
       * the periodicity and shift are guaranteed to be consistent
       * between different nodes when running e.g. polymers in
       * parallel. In this special case we thus unshift/shift, but
       * only when necessary. This is to make sure the coordinates
       * we move don't end up a box away...
       */
      if(dummycomm) {
	unshift_self(graph,parm->box,pos[TRY]);
	move_construct_x(dummycomm,pos[TRY],cr);
	shift_self(graph,parm->box,pos[TRY]);
      }
      construct_dummies(log,pos[TRY],&(nrnb[cr->nodeid]),1,NULL,&top->idef);

      if(dummycomm) {
        unshift_self(graph,parm->box,pos[TRY]);
        move_dummy_xv(dummycomm,pos[TRY],NULL,cr);
	shift_self(graph,parm->box,pos[TRY]); /* maybe not necessary */
      }
    }

    /* Calc force & energy on new positions  */
    /* do_force always puts the charge groups in the box and shifts again
     * We do not unshift, so molecules are always whole in steep.c
     */
    do_force(log,cr,parm,nsb,force_vir,pme_vir,
 	     count,&(nrnb[cr->nodeid]),top,grps,pos[TRY],buf,force[TRY],buf,
	     mdatoms,ener,bVerbose && !(PAR(cr)), 
 	     lambda,graph,parm->ir.nstlist>0 || count==0,FALSE,fr,mu_tot,FALSE); 

    /* Spread the force on dummy particle to the other particles... */
    if(bDummies) {
      /* We only move forces here, and they are independent of shifts */
      if(dummycomm)
        move_dummy_f(dummycomm,force[TRY],cr);
      
      spread_dummy_f(log,pos[TRY],force[TRY],buf,&(nrnb[cr->nodeid]),
		     &top->idef);
      if(dummycomm)
        move_construct_f(dummycomm,force[TRY],cr);
    }

    /* Sum the potential energy terms from group contributions  */
    sum_epot(&(ir->opts),grps,ener); 

    if (MASTER(cr))
      print_ebin_header(log,count,count,lambda,0.0);

    if (bConstrain) {
      fmax=f_max(cr->left,cr->right,nsb->nnodes,&(ir->opts),mdatoms,start,end,
		 force[TRY]);
      constepsize=ustep/fmax;
      for(i=start; (i<end); i++)  
	for(m=0;(m<DIM);m++) 
	  xcf[i][m] = pos[TRY][i][m] + constepsize*force[TRY][i][m];
      
      dvdlambda=0;
      constrain(stdlog,top,&(parm->ir),count,mdatoms,start,end,
		pos[TRY],xcf,NULL,parm->box,lambda,&dvdlambda,nrnb,TRUE);
      
      for(i=start; (i<end); i++)  
	for(m=0;(m<DIM);m++) 
	  force[TRY][i][m] = (xcf[i][m] - pos[TRY][i][m])/constepsize;
    }

    /* Clear stuff again  */
    clear_mat(force_vir); 
    clear_mat(shake_vir); 
    
    /* Communicat stuff when parallel  */
    if (PAR(cr))  
      global_stat(log,cr,ener,force_vir,shake_vir, 
 		  &(ir->opts),grps,&mynrnb,nrnb,vcm,&terminate); 
    
    /* This is the new energy  */
    Fmax[TRY]=f_max(cr->left,cr->right,nsb->nnodes,&(ir->opts),mdatoms,
		    start,end,force[TRY]);
    Epot[TRY]=ener[F_EPOT];
    if (count == 0)
      Epot[Min] = Epot[TRY]+1;
      
    /* Print it if necessary  */
    if (MASTER(cr)) { 
      if (bVerbose) {
	fprintf(stderr,"Step = %5d, Dmax = %7.2e nm, Epot = %12.5e Fmax = %11.5e%c",
		count,ustep,Epot[TRY],Fmax[TRY],(Epot[TRY]<Epot[Min])?'\n':'\r');
      }

      if (Epot[TRY] < Epot[Min]) {
	/* Store the new (lower) energies  */
	upd_mdebin(mdebin,NULL,mdatoms->tmass,count,(real)count,
		   ener,parm->box,shake_vir, 
		   force_vir,parm->vir,parm->pres,grps,mu_tot,(parm->ir.etc==etcNOSEHOOVER)); 
	print_ebin(fp_ene,TRUE,FALSE,log,count,count,
		   eprNORMAL,TRUE,mdebin,&(top->atoms));
	fflush(log);
      }
    } 
    
    /* Now if the new energy is smaller than the previous...  
     * or if this is the first step!
     * or if we did random steps! 
     */
    
    if ( (count==0) || (Epot[TRY] < Epot[Min]) ) {
      steps_accepted++; 
      if (do_per_step(steps_accepted,ir->nstfout)) 
	ff=force[TRY];  
      else 
	ff=NULL;    
      if (do_per_step(steps_accepted,ir->nstxout)) {
	xx=pos[TRY];   
	write_traj(log,cr,ftp2fn(efTRN,nfile,fnm), 
		   nsb,count,(real) count, 
		   lambda,nrnb,nsb->natoms,xx,NULL,ff,parm->box); 
      } else 
	xx=NULL; 
      
      /* Test whether the convergence criterion is met...  */
      bDone=(Fmax[TRY] < ftol);
      
      /* Copy the arrays for force, positions and energy  */
      /* The 'Min' array always holds the coords and forces of the minimal 
	 sampled energy  */
      Min = TRY; 
      
      /* increase step size  */
      if (count>0)
	ustep *= 1.2; 
      
    } else
      /* If energy is not smaller make the step smaller...  */
      ustep *= 0.5;
    
    /* Determine new step  */
    fmax = f_max(cr->left,cr->right,nsb->nnodes,&(ir->opts),mdatoms,start,end,
		 force[Min]);
    stepsize=ustep/fmax;
    
    /* Check if stepsize is too small, with 1 nm as a characteristic length */
#ifdef DOUBLE
    if (ustep < 1e-12) {
#else
    if (ustep < 1e-6) {
#endif
      sprintf(sbuf,"\nStepsize too small "
	      "Converged to machine precision,\n"
	      "but not to the requested precision (%g)\n%s",
	      ftol,bConstrain ? "You might need to increase your constraint accuracy\n" : "");
      fprintf(stderr,sbuf);
      fprintf(log,sbuf);
      bAbort=TRUE;
    }
    
    count++;
  } /* End of the loop  */
    
  /* Print some shit...  */
  if (MASTER(cr)) 
    fprintf(stderr,"\nwriting lowest energy coordinates.\n"); 
  xx=pos[Min]; 
  ff=force[Min]; 
  write_traj(log,cr,ftp2fn(efTRN,nfile,fnm), 
	     nsb,count,(real) count, 
	     lambda,nrnb,nsb->natoms,xx,NULL,ff,parm->box); 
  if (MASTER(cr)) {
    write_sto_conf(ftp2fn(efSTO,nfile,fnm),
		   *top->name, &(top->atoms),xx,NULL,parm->box);
    
    fprintf(stderr,"Maximum force: %12.5e\n",Fmax[Min]); 
    if (bDone)
      sprintf(sbuf,"\n%s converged to %g\n",SD,ftol); 
    else 
      sprintf(sbuf,"\n%s did not converge in %d steps\n",SD,min(count,nsteps));
    fprintf(stderr,sbuf);
    fprintf(log,sbuf);
    sprintf(sbuf,"  Potential Energy  = %12.5e\n",Epot[Min]); 
    fprintf(stderr,sbuf);
    fprintf(log,sbuf);
  }
  if (MASTER(cr))
    close_enx(fp_ene);
    
  /* Put the coordinates back in the x array (otherwise the whole
   * minimization would be in vain)
   */
  for(i=0; (i<nsb->natoms); i++)
    copy_rvec(pos[Min][i],x[i]);
  
  /* To print the actual number of steps we needed somewhere */
  ir->nsteps=count;
  
  return start_t;
} /* That's all folks */

