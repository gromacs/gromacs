/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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

 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Getting the Right Output Means no Artefacts in Calculating Stuff
 */

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
#include "constr.h"
#include "vec.h"
#include "statutil.h"
#include "tgroup.h"
#include "mdebin.h"
#include "dummies.h"
#include "mdrun.h"

static void sp_header(FILE *out,const char *minimizer,real ftol,int nsteps)
{
  fprintf(out,"%s:\n",minimizer);
  fprintf(out,"   Tolerance         = %12.5e\n",ftol);
  fprintf(out,"   Number of steps   = %12d\n",nsteps);
}

static void warn_step(FILE *fp,real ustep,real ftol,bool bConstrain)
{
  fprintf(fp,"\nStepsize too small (%g nm)"
	  "Converged to machine precision,\n"
	  "but not to the requested precision (%g)\n",
	  ustep,ftol);
  if (bConstrain)
    fprintf(fp,"You might need to increase your constraint accuracy, or turn\n"
	    "off constraints alltogether (set constraints = none in mdp file)\n");
}

static void print_converged(FILE *fp,const char *alg,real ftol,int count,bool bDone,
			    int nsteps,real epot,real fmax)
{
  if (bDone)
    fprintf(fp,"\n%s converged to %g in %d steps\n",alg,ftol,count); 
  else 
    fprintf(fp,"\n%s did not converge in %d steps\n",alg,min(count,nsteps));
  fprintf(fp,"  Potential Energy  = %12.5e\n",epot); 
  fprintf(fp,"Maximum force: %12.5e\n",fmax); 
}

static real f_max(int left,int right,int nnodes,
		  t_grpopts *opts,t_mdatoms *mdatoms,
		  int start,int end,rvec grad[],int *nfmax)
{
  real fmax2,fmax2_0,fam,nfm;
  int  i,m,gf;

  /* This routine finds the largest force and returns it.
   * On parallel machines the global max is taken.
   */
  fmax2 = 0;
  nfm   = -1;
  for(i=start; i<end; i++) {
    gf = mdatoms->cFREEZE[i];
    fam = 0;
    for(m=0; m<DIM; m++)
      if (!opts->nFreeze[gf][m])
	fam += sqr(grad[i][m]);
    if (fam > fmax2) {
      fmax2 = fam;
      nfm   = i;
    }
  }

  *nfmax = nfm;
  if (nnodes > 1) {
    for(i=0; (i<nnodes-1); i++) {
      gmx_tx(left,(void *)&fmax2,sizeof(fmax2));
      gmx_rx(right,(void *)&fmax2_0,sizeof(fmax2_0));
      gmx_wait(left,right);
      gmx_tx(left,(void *)nfmax,sizeof(*nfmax));
      gmx_rx(right,(void *)&nfm,sizeof(nfm));
      gmx_wait(left,right);
      if (fmax2_0 > fmax2) {
	fmax2  = fmax2_0;
	*nfmax = nfm;
      }
    }
  }
  return sqrt(fmax2);
}

static real f_norm(t_commrec *cr,
		   t_grpopts *opts,t_mdatoms *mdatoms,
		   int start,int end,rvec grad[])
{
  double fnorm2;
  int    i,m,gf;

  /* This routine returns the norm of the force */
  fnorm2 = 0;
  
  for(i=start; i<end; i++) {
    gf = mdatoms->cFREEZE[i];
    for(m=0; m<DIM; m++)
      if (!opts->nFreeze[gf][m])
	fnorm2 += sqr(grad[i][m]); 
  } 
  
  if (PAR(cr))
    gmx_sumd(1,&fnorm2,cr);

  return sqrt(fnorm2); 
} 

static void init_em(FILE *log,const char *title,
		    t_parm *parm,real *lambda,t_nrnb *mynrnb,rvec mu_tot,rvec box_size,
		    t_forcerec *fr,t_mdatoms *mdatoms,t_topology *top,t_nsborder *nsb,
		    t_commrec *cr,t_vcm **vcm,int *start,int *end)
{
  fprintf(log,"Initiating %s\n",title);

  /* Initiate some variables */
  if (parm->ir.efep != efepNO)
    *lambda = parm->ir.init_lambda;
  else 
    *lambda = 0.0;

  init_nrnb(mynrnb);
    
  clear_rvec(mu_tot);
  calc_shifts(parm->box,box_size,fr->shift_vec);
  
  *start = nsb->index[cr->nodeid];
  *end   = nsb->homenr[cr->nodeid] + *start;

  /* Set initial values for invmass etc. */
  update_mdatoms(mdatoms,*lambda,TRUE);

  *vcm = init_vcm(log,top,cr,mdatoms,
		  *start,HOMENR(nsb),parm->ir.nstcomm);
}

time_t do_cg(FILE *log,int nfile,t_filenm fnm[],
	     t_parm *parm,t_topology *top,
	     t_groups *grps,t_nsborder *nsb,
	     rvec x[],rvec grad[],rvec buf[],t_mdatoms *mdatoms,
	     tensor ekin,real ener[],t_fcdata *fcd,t_nrnb nrnb[],
	     bool bVerbose,bool bDummies,t_comm_dummies *dummycomm,
	     t_commrec *cr,t_commrec *mcr,t_graph *graph,
	     t_forcerec *fr,rvec box_size)
{
  const char *CG="Conjugate Gradients";
  double gpa,gpb;
  double EpotA=0.0,EpotB=0.0,a=0.0,b,beta=0.0,zet,w;
  real   lambda,fmax,testf,smin;
  rvec   *p,*f,*xprime,*xx,*ff;
  real   fnorm,pnorm,fnorm_old;
  t_vcm      *vcm;
  t_mdebin   *mdebin;
  t_nrnb mynrnb;
  bool   bNS=TRUE,bDone,bLR,bLJLR,bBHAM,b14,bRand,brerun,bBeta0=FALSE;
  rvec   mu_tot;
  time_t start_t;
  tensor force_vir,shake_vir,pme_vir;
  int    number_steps,naccept=0,nstcg=parm->ir.nstcgsteep;
  int    fp_ene,count=0;
  int    i,m,nfmax,start,end,niti,gf;
  /* not used */
  real   terminate=0;

  init_em(log,CG,parm,&lambda,&mynrnb,mu_tot,box_size,fr,mdatoms,top,
	  nsb,cr,&vcm,&start,&end);
  
  /* Print to log file */
  start_t=print_date_and_time(log,cr->nodeid,"Started Conjugate Gradients");
  
  /* p is the search direction, f the force, xprime the new positions */
  snew(p,nsb->natoms);
  snew(f,nsb->natoms);
  snew(xprime,nsb->natoms);

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
		     parm->ir.eDispCorr,TRICLINIC(parm->ir.compress),
		     (parm->ir.etc==etcNOSEHOOVER),cr); 

  /* Clear some matrix variables */
  clear_mat(force_vir);
  clear_mat(shake_vir);
  
  if (fr->ePBC != epbcNONE)
    /* Remove periodicity */
    do_pbc_first(log,parm,box_size,fr,graph,x);
  
  /* Max number of steps */
  number_steps=parm->ir.nsteps;

  if (MASTER(cr))
    sp_header(stderr,CG,parm->ir.em_tol,number_steps);
  sp_header(log,CG,parm->ir.em_tol,number_steps);

  if (bDummies)
    construct_dummies(log,x,&(nrnb[cr->nodeid]),1,NULL,&top->idef,
		      graph,cr,parm->box,dummycomm);
  
  /* Call the force routine and some auxiliary (neighboursearching etc.) */
  /* do_force always puts the charge groups in the box and shifts again
   * We do not unshift, so molecules are always whole in congrad.c
   */
  do_force(log,cr,mcr,parm,nsb,force_vir,pme_vir,0,&(nrnb[cr->nodeid]),
	   top,grps,x,buf,f,buf,mdatoms,ener,fcd,bVerbose && !(PAR(cr)),
	   lambda,graph,bNS,FALSE,fr,mu_tot,FALSE,0.0);
  where();

  /* Spread the force on dummy particle to the other particles... */
  if (bDummies)
    spread_dummy_f(log,x,f,&(nrnb[cr->nodeid]),&top->idef,dummycomm,cr);

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
    
    print_ebin_header(log,count,count,lambda);
    print_ebin(fp_ene,TRUE,FALSE,FALSE,FALSE,log,count,count,eprNORMAL,
	       TRUE,mdebin,fcd,&(top->atoms),&(parm->ir.opts));
  }
  where();
  
  /* This is the starting energy */
  EpotA=ener[F_EPOT];

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
  bDone=FALSE; 
  for(count=1; (!bDone && (count < number_steps)); count++) {
    
    /* start conjugate gradient, determine search interval a,b */
    gpa = 0.0;
    for(i=start; (i<end); i++) {
      gf = mdatoms->cFREEZE[i];
      for(m=0; m<DIM; m++) 
	if (!parm->ir.opts.nFreeze[gf][m]) {
	  p[i][m] = f[i][m] + beta*p[i][m];
	  gpa     = gpa - p[i][m]*f[i][m];
	}
    }
    if (PAR(cr))
      gmx_sumd(1,&gpa,cr);
    pnorm = f_norm(cr,&(parm->ir.opts),mdatoms,start,end,p);

    /* Set variables for stepsize (in nm). This is the largest 
     * step that we are going to make in any direction.
     */
    a    = 0.0;
    b    = parm->ir.em_stepsize/pnorm;
    niti = 0;

    /* search a,b iteratively, if necessary */
    brerun=TRUE;
    while (brerun) {
      for (i=start; (i<end); i++) {
	for(m=0; (m<DIM); m++) {
	  xprime[i][m] = x[i][m] + b*p[i][m];
	}
      }
      bNS = (parm->ir.nstlist > 0);

      if (bDummies)
	construct_dummies(log,xprime,&(nrnb[cr->nodeid]),1,NULL,&top->idef,
			  graph,cr,parm->box,dummycomm);

      
      /* Calc force & energy on new trial position  */
      /* do_force always puts the charge groups in the box and shifts again
       * We do not unshift, so molecules are always whole in conjugate
       * gradients
       */
      do_force(log,cr,mcr,parm,nsb,force_vir,pme_vir,
	       count,&(nrnb[cr->nodeid]),top,grps,xprime,buf,f,
	       buf,mdatoms,ener,fcd,bVerbose && !(PAR(cr)),
	       lambda,graph,bNS,FALSE,fr,mu_tot,FALSE,0.0);
      
      /* Spread the force on dummy particle to the other particles... */
      if (bDummies) 
	spread_dummy_f(log,xprime,f,&(nrnb[cr->nodeid]),&top->idef,
		       dummycomm,cr); 
      
      /* Sum the potential energy terms from group contributions */
      sum_epot(&(parm->ir.opts),grps,ener);
      where();

      bNS = (parm->ir.nstlist > 0);

      gpb=0.0;
      for(i=start; (i<end); i++)
	gpb -= iprod(p[i],f[i]);

      if (PAR(cr))
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

      EpotB = ener[F_EPOT];
      
      if ((gpb >= 0.0) || (EpotB >= EpotA))
	brerun=FALSE;
      else {
	a     = b;
	EpotA = EpotB;
	gpa   = gpb;
	b     = 2.0*b;
      }
      niti++;
    }
    /* End of while loop searching for a and b */
    
    /* find stepsize smin in interval a-b */
    zet = 3.0 * (EpotA-EpotB) / (b-a) + gpa + gpb;
    w   = zet*zet - gpa*gpb;
    if (w < 0.0) {
      if (debug) {
	fprintf(debug,"Negative w: %20.12e\n",w);
	fprintf(debug,"z= %20.12e\n",zet);
	fprintf(debug,"gpa= %20.12e, gpb= %20.12e\n",gpa,gpb);
	fprintf(debug,"a= %20.12e, b= %20.12e\n",a,b);
	fprintf(debug,"EpotA= %20.12e, EpotB= %20.12e\n",EpotA,EpotB);
      }
      bBeta0 = TRUE;      
    } 
    else {
      w    = sqrt(w);
      smin = b - ((gpb+w-zet)*(b-a))/((gpb-gpa)+2.0*w);
      
      /* new positions */
      for (i=start; i<end; i++)
	for(m=0; m<DIM; m++) 
	  xprime[i][m] = x[i][m] + smin*p[i][m];
      
      if (bDummies) 
	construct_dummies(log,xprime,&(nrnb[cr->nodeid]),1,NULL,&top->idef,
			  graph,cr,parm->box,dummycomm);
      
      /* new energy, forces
       * do_force always puts the charge groups in the box and shifts again
       * We do not unshift, so molecules are always whole in congrad.c
       */
      do_force(log,cr,mcr,parm,nsb,force_vir,pme_vir,
	       count,&(nrnb[cr->nodeid]),top,grps,xprime,buf,f,
	       buf,mdatoms,ener,fcd,bVerbose && !(PAR(cr)),
	       lambda,graph,bNS,FALSE,fr,mu_tot,FALSE,0.0);
      
      /* Spread the force on dummy particle to the other particles... */
      if(bDummies)
	spread_dummy_f(log,xprime,f,&(nrnb[cr->nodeid]),&top->idef,dummycomm,cr); 
      
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
      
      EpotA  = ener[F_EPOT];
      
      bBeta0 = FALSE;
    }   
    /* new search direction */
    /* beta = 0 means steepest descents */
    if (bBeta0 || (nstcg && ((count % nstcg)==0))) {
      if (bVerbose)
	fprintf(stderr,"\rNew search direction\n");
      beta = 0.0;
    }
    else
      beta = sqr(fnorm/fnorm_old);
    
    /* update x, fnorm_old */
    for (i=start; i<end; i++)
      copy_rvec(xprime[i],x[i]);
    fnorm_old=fnorm;
    
    /* Test whether the convergence criterion is met */
    fmax=f_max(cr->left,cr->right,nsb->nnodes,&(parm->ir.opts),mdatoms,
	       start,end,f,&nfmax);
    
    /* Print it if necessary */
    if (bVerbose && MASTER(cr)) {
      fprintf(stderr,"\rStep %d, E-Pot = %16.10e, Fmax = %12.5e, atom = %d\n",
	      count,EpotA,fmax,nfmax);
	/* Store the new (lower) energies */
      upd_mdebin(mdebin,NULL,mdatoms->tmass,count,(real)count,
		 ener,parm->box,shake_vir,
		 force_vir,parm->vir,parm->pres,grps,mu_tot,
		 (parm->ir.etc==etcNOSEHOOVER));
      /* Print the energies allways when we should be verbose */
      print_ebin_header(log,count,count,lambda);
      print_ebin(fp_ene,TRUE,FALSE,FALSE,FALSE,log,count,count,eprNORMAL,
		 TRUE,mdebin,fcd,&(top->atoms),&(parm->ir.opts));
    }
    
    /* Stop when the maximum force lies below tolerance */
    bDone=(fmax < parm->ir.em_tol);
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
	     start,end,f,&nfmax);
  if (MASTER(cr)) {
    print_converged(stderr,CG,parm->ir.em_tol,count,bDone,
		    number_steps,EpotA,fmax);
    print_converged(log,CG,parm->ir.em_tol,count,bDone,
		    number_steps,EpotA,fmax);
    close_enx(fp_ene);
  }
  
  /* To print the actual number of steps we needed somewhere */
  parm->ir.nsteps=count;

  return start_t;
} /* That's all folks */

time_t do_steep(FILE *log,int nfile,t_filenm fnm[], 
 		t_parm *parm,t_topology *top, 
 		t_groups *grps,t_nsborder *nsb, 
 		rvec x[],rvec grad[],rvec buf[],t_mdatoms *mdatoms, 
 		tensor ekin,real ener[],t_fcdata *fcd,t_nrnb nrnb[], 
 		bool bVerbose,bool bDummies, t_comm_dummies *dummycomm,
		t_commrec *cr,t_commrec *mcr,t_graph *graph,
		t_forcerec *fr,rvec box_size) 
{ 
  const char *SD="Steepest Descents"; 
  real   stepsize,constepsize,lambda,fmax; 
  rvec   *pos[2],*force[2],*xcf=NULL; 
  rvec   *xx,*ff; 
  real   Fmax[2]; 
  real   Epot[2]; 
  real   ustep,dvdlambda;
  t_vcm      *vcm;
  int        fp_ene; 
  t_mdebin   *mdebin; 
  t_nrnb     mynrnb; 
  bool   bDone,bAbort,bLR,bLJLR,bBHAM,b14; 
  time_t start_t; 
  tensor force_vir,shake_vir,pme_vir; 
  rvec   mu_tot;
  int    nfmax,nsteps;
  int    count=0; 
  int    i,m,start,end,gf; 
  int    Min=0; 
  int    steps_accepted=0; 
  bool   bConstrain;
  /* not used */
  real   terminate=0;
#define  TRY (1-Min)

  init_em(log,SD,parm,&lambda,&mynrnb,mu_tot,box_size,fr,mdatoms,top,
	  nsb,cr,&vcm,&start,&end);
   
  /* Print to log file  */
  start_t=print_date_and_time(log,cr->nodeid,"Started Steepest Descents"); 
  
  /* We need two coordinate arrays and two force arrays  */
  for(i=0; (i<2); i++) { 
    snew(pos[i],nsb->natoms); 
    snew(force[i],nsb->natoms); 
  } 

  /* Set some booleans for the epot routines  */
  set_pot_bools(&(parm->ir),top,&bLR,&bLJLR,&bBHAM,&b14);
  
  /* Open the enrgy file */   
  if (MASTER(cr)) 
    fp_ene=open_enx(ftp2fn(efENX,nfile,fnm),"w"); 
  else 
    fp_ene=-1; 
  
  /* Init bin for energy stuff  */
  mdebin=init_mdebin(fp_ene,grps,&(top->atoms),&(top->idef),bLR,bLJLR,
		     bBHAM,b14,parm->ir.efep!=efepNO,parm->ir.epc,
		     parm->ir.eDispCorr,TRICLINIC(parm->ir.compress),
		     (parm->ir.etc==etcNOSEHOOVER),cr); 
  
  /* Clear some matrix variables  */
  clear_mat(force_vir); 
  clear_mat(shake_vir); 

  if (fr->ePBC != epbcNONE)
    /* Remove periodicity */
    do_pbc_first(log,parm,box_size,fr,graph,x);

  /* Initiate constraint stuff */
  bConstrain=init_constraints(stdlog,top,&(parm->ir),mdatoms,
			      start,end,FALSE,cr);
  
  if (bConstrain)
    snew(xcf,nsb->natoms); 

  /* Copy coord vectors to our temp array  */
  for(i=0; (i<nsb->natoms); i++) { 
    copy_rvec(x[i],pos[Min][i]); 
    copy_rvec(x[i],pos[TRY][i]); 
  } 
    
  /* Set variables for stepsize (in nm). This is the largest  
   * step that we are going to make in any direction. 
   */
  ustep    = parm->ir.em_stepsize; 
  stepsize = 0;

  /* Max number of steps  */
  nsteps=parm->ir.nsteps; 
  
  if (MASTER(cr)) 
    /* Print to the screen  */
    sp_header(stderr,SD,parm->ir.em_tol,nsteps);
  sp_header(log,SD,parm->ir.em_tol,nsteps);
    
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
	  if (parm->ir.opts.nFreeze[gf][m])
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

    if (bDummies)
      construct_dummies(log,pos[TRY],&(nrnb[cr->nodeid]),1,NULL,&top->idef,
			graph,cr,parm->box,dummycomm);

    /* Calc force & energy on new positions
     * do_force always puts the charge groups in the box and shifts again
     * We do not unshift, so molecules are always whole in steep.c
     */
    do_force(log,cr,mcr,parm,nsb,force_vir,pme_vir,
 	     count,&(nrnb[cr->nodeid]),top,grps,pos[TRY],buf,force[TRY],buf,
	     mdatoms,ener,fcd,bVerbose && !(PAR(cr)), 
 	     lambda,graph,parm->ir.nstlist>0 || count==0,FALSE,fr,mu_tot,
	     FALSE,0.0); 

    /* Spread the force on dummy particle to the other particles... */
    if (bDummies) 
      spread_dummy_f(log,pos[TRY],force[TRY],&(nrnb[cr->nodeid]),
		     &top->idef,dummycomm,cr);
    
    /* Sum the potential energy terms from group contributions  */
    sum_epot(&(parm->ir.opts),grps,ener); 

    if (MASTER(cr))
      print_ebin_header(log,count,count,lambda);

    if (bConstrain) {
      fmax=f_max(cr->left,cr->right,nsb->nnodes,&(parm->ir.opts),mdatoms,start,end,
		 force[TRY],&nfmax);
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
 		  &(parm->ir.opts),grps,&mynrnb,nrnb,vcm,&terminate); 
    
    /* This is the new energy  */
    Fmax[TRY]=f_max(cr->left,cr->right,nsb->nnodes,&(parm->ir.opts),mdatoms,
		    start,end,force[TRY],&nfmax);
    Epot[TRY]=ener[F_EPOT];
    if (count == 0)
      Epot[Min] = Epot[TRY]+1;
      
    /* Print it if necessary  */
    if (MASTER(cr)) { 
      if (bVerbose) {
	fprintf(stderr,"Step=%5d, Dmax= %6.1e nm, Epot= %12.5e Fmax= %11.5e, atom= %d%c",
		count,ustep,Epot[TRY],Fmax[TRY],nfmax+1,
		(Epot[TRY]<Epot[Min])?'\n':'\r');
      }

      if (Epot[TRY] < Epot[Min]) {
	/* Store the new (lower) energies  */
	upd_mdebin(mdebin,NULL,mdatoms->tmass,count,(real)count,
		   ener,parm->box,shake_vir, 
		   force_vir,parm->vir,parm->pres,grps,mu_tot,(parm->ir.etc==etcNOSEHOOVER)); 
	print_ebin(fp_ene,TRUE,
		   do_per_step(steps_accepted,parm->ir.nstdisreout),
		   do_per_step(steps_accepted,parm->ir.nstorireout),
		   do_per_step(steps_accepted,parm->ir.nstdihreout),
		   log,count,count,eprNORMAL,TRUE,mdebin,fcd,&(top->atoms),&(parm->ir.opts));
	fflush(log);
      }
    } 
    
    /* Now if the new energy is smaller than the previous...  
     * or if this is the first step!
     * or if we did random steps! 
     */
    
    if ( (count==0) || (Epot[TRY] < Epot[Min]) ) {
      steps_accepted++; 
      if (do_per_step(steps_accepted,parm->ir.nstfout)) 
	ff=force[TRY];  
      else 
	ff=NULL;    
      if (do_per_step(steps_accepted,parm->ir.nstxout)) {
	xx=pos[TRY];   
	write_traj(log,cr,ftp2fn(efTRN,nfile,fnm), 
		   nsb,count,(real) count, 
		   lambda,nrnb,nsb->natoms,xx,NULL,ff,parm->box); 
      } else 
	xx=NULL; 
      
      /* Test whether the convergence criterion is met...  */
      bDone=(Fmax[TRY] < parm->ir.em_tol);
      
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
    fmax = f_max(cr->left,cr->right,nsb->nnodes,&(parm->ir.opts),mdatoms,start,end,
		 force[Min],&nfmax);
    stepsize=ustep/fmax;
    
    /* Check if stepsize is too small, with 1 nm as a characteristic length */
#ifdef DOUBLE
    if (ustep < 1e-12) {
#else
    if (ustep < 1e-6) {
#endif
      warn_step(stderr,ustep,parm->ir.em_tol,bConstrain);
      warn_step(log,ustep,parm->ir.em_tol,bConstrain);
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
    
    print_converged(stderr,SD,parm->ir.em_tol,count,bDone,nsteps,Epot[Min],Fmax[Min]);
    print_converged(log,SD,parm->ir.em_tol,count,bDone,nsteps,Epot[Min],Fmax[Min]);
    close_enx(fp_ene);
  }
      
  /* Put the coordinates back in the x array (otherwise the whole
   * minimization would be in vain)
   */
  for(i=0; (i<nsb->natoms); i++)
    copy_rvec(pos[Min][i],x[i]);
  
  /* To print the actual number of steps we needed somewhere */
  parm->ir.nsteps=count;
  
  return start_t;
} /* That's all folks */

time_t do_nm(FILE *log,t_commrec *cr,int nfile,t_filenm fnm[],
	     bool bVerbose,bool bCompact,int stepout,
	     t_parm *parm,t_groups *grps,
	     t_topology *top,real ener[],t_fcdata *fcd,
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
  int        i,m,nfmax;
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

#ifndef DOUBLE
  fprintf(stderr,
	  "WARNING: This version of GROMACS has been compiled in single precision,\n"
	  "         which is usually not accurate enough for normal mode analysis.\n"
	  "         For reliable results you should compile GROMACS in double precision.\n\n");
#endif

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
    calc_shifts(parm->box,box_size,fr->shift_vec);
    fprintf(log,"Removing pbc first time\n");
    mk_mshift(log,graph,parm->box,x);
    shift_self(graph,parm->box,x);
  }

  fp_ene=-1;
  set_pot_bools(&(parm->ir),top,&bLR,&bLJLR,&bBHAM,&b14);
  mdebin=init_mdebin(fp_ene,grps,&(top->atoms),&(top->idef),bLR,bLJLR,
		     bBHAM,b14,parm->ir.efep!=efepNO,parm->ir.epc,
		     parm->ir.eDispCorr,TRICLINIC(parm->ir.compress),(parm->ir.etc==etcNOSEHOOVER),cr);

   ener[F_TEMP]=sum_ekin(&(parm->ir.opts),grps,parm->ekin,bTYZ);
  if(parm->ir.etc==etcBERENDSEN)
    berendsen_tcoupl(&(parm->ir.opts),grps,parm->ir.delta_t);
  else if(parm->ir.etc==etcNOSEHOOVER)
    nosehoover_tcoupl(&(parm->ir.opts),grps,parm->ir.delta_t);

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
  do_force(log,cr,NULL,parm,nsb,force_vir,pme_vir,0,&mynrnb,
	   top,grps,x,v,f,buf,mdatoms,ener,fcd,bVerbose && !PAR(cr),
	   lambda,graph,bNS,FALSE,fr,mu_tot,FALSE,0.0);
  bNS=FALSE;
  if (bBox)
    /* Shift back the coordinates, since we're not calling update */
    unshift_self(graph,parm->box,x);

  
  /* if forces not small, warn user */
  fmax=f_max(cr->left,cr->right,nsb->nnodes,&(parm->ir.opts),mdatoms,
	     0,top->atoms.nr,f,&nfmax);
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

      do_force(log,cr,NULL,parm,nsb,force_vir,pme_vir,2*(step*DIM+idum),
	       &mynrnb,
	       top,grps,x,v,f,buf,mdatoms,ener,fcd,bVerbose && !PAR(cr),
	       lambda,graph,bNS,FALSE,fr,mu_tot,FALSE,0.0);
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
      
      do_force(log,cr,NULL,parm,nsb,force_vir,pme_vir,2*(step*DIM+idum)+1,
	       &mynrnb,
	       top,grps,x,v,f,buf,mdatoms,ener,fcd,bVerbose && !PAR(cr),
	       lambda,graph,bNS,FALSE,fr,mu_tot,FALSE,0.0);
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
    print_ebin(-1,FALSE,FALSE,FALSE,FALSE,log,step,t,eprAVER,
	       FALSE,mdebin,fcd,&(top->atoms),&(parm->ir.opts));
    print_ebin(-1,FALSE,FALSE,FALSE,FALSE,log,step,t,eprRMS,
	       FALSE,mdebin,fcd,&(top->atoms),&(parm->ir.opts));
  }
  
  /* Construct dummy particles, for last output frame */
  /* NB: I have not included the communication for parallel
   * dummies here, since the rest of nm doesn't seem to
   * be parallelized. Be sure to copy the correct code from
   * e.g. md.c or steep.c if you make nm parallel!
   */
  construct_dummies(log,x,&mynrnb,parm->ir.delta_t,v,&top->idef,
		    graph,cr,parm->box,NULL);
    
  /*free_nslist(log);*/
  
  return start_t;
}
