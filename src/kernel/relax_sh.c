/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include "typedefs.h"
#include "smalloc.h"
#include "gmx_fatal.h"
#include "vec.h"
#include "txtdump.h"
#include "mdrun.h"
#include "partdec.h"
#include "xmdrun.h"
#include "mdatoms.h"
#include "vsite.h"
#include "network.h"
#include "names.h"
#include "constr.h"

static void do_1pos(rvec xnew,rvec xold,rvec f,real step)
{
  real xo,yo,zo;
  real dx,dy,dz;
  
  xo=xold[XX];
  yo=xold[YY];
  zo=xold[ZZ];

  dx=f[XX]*step;
  dy=f[YY]*step;
  dz=f[ZZ]*step;

  xnew[XX]=xo+dx;
  xnew[YY]=yo+dy;
  xnew[ZZ]=zo+dz;
}

static void do_1pos3(rvec xnew,rvec xold,rvec f,rvec step)
{
  real xo,yo,zo;
  real dx,dy,dz;
  
  xo=xold[XX];
  yo=xold[YY];
  zo=xold[ZZ];

  dx=f[XX]*step[XX];
  dy=f[YY]*step[YY];
  dz=f[ZZ]*step[ZZ];

  xnew[XX]=xo+dx;
  xnew[YY]=yo+dy;
  xnew[ZZ]=zo+dz;
}

static void directional_sd(FILE *log,rvec xold[],rvec xnew[],rvec acc_dir[],
			   int start,int homenr,real step)
{
  int  i;

  for(i=start; i<homenr; i++)
    do_1pos(xnew[i],xold[i],acc_dir[i],step);
}

static void shell_pos_sd(FILE *log,rvec xcur[],rvec xnew[],rvec f[],
			 int ns,t_shell s[],int count)
{
  int  i,shell,d;
  real dx,df,k_est;
#ifdef PRINT_STEP  
  real step_min,step_max;

  step_min = 1e30;
  step_max = 0;
#endif
  for(i=0; (i<ns); i++) {
    shell = s[i].shell;
    if (count == 1) {
      for(d=0; d<DIM; d++) {
	s[i].step[d] = s[i].k_1;
#ifdef PRINT_STEP
	step_min = min(step_min,s[i].step[d]);
	step_max = max(step_max,s[i].step[d]);
#endif
      }
    } else {
      for(d=0; d<DIM; d++) {
	dx = xcur[shell][d] - s[i].xold[d];
	df =    f[shell][d] - s[i].fold[d];
	if (dx != 0 && df != 0) {
	  k_est = -dx/df;
	  if (k_est >= 2*s[i].step[d]) {
	    s[i].step[d] *= 1.2;
	  } else if (k_est <= 0) {
	    s[i].step[d] *= 0.8;
	  } else {
	    s[i].step[d] = 0.8*s[i].step[d] + 0.2*k_est;
	  }
	} else if (dx != 0) {
	  s[i].step[d] *= 1.2;
	}
#ifdef PRINT_STEP
	step_min = min(step_min,s[i].step[d]);
	step_max = max(step_max,s[i].step[d]);
#endif
      }
    }
    copy_rvec(xcur[shell],s[i].xold);
    copy_rvec(f[shell],   s[i].fold);

    do_1pos3(xnew[shell],xcur[shell],f[shell],s[i].step);

    if (gmx_debug_at) {
      pr_rvec(debug,0,"fshell",f[shell],DIM,TRUE);
      pr_rvec(debug,0,"xold",xcur[shell],DIM,TRUE);
      pr_rvec(debug,0,"xnew",xnew[shell],DIM,TRUE);
    }
  }
#ifdef PRINT_STEP
  printf("step %.3e %.3e\n",step_min,step_max);
#endif
}

static void decrease_step_size(int nshell,t_shell s[])
{
  int i;
  
  for(i=0; i<nshell; i++)
    svmul(0.8,s[i].step,s[i].step);
}

static void predict_shells(FILE *log,rvec x[],rvec v[],real dt,
			   int ns,t_shell s[],
			   real mass[],bool bInit)
{
  int  i,m,s1,n1,n2,n3;
  real dt_1,dt_2,dt_3,fudge,tm,m1,m2,m3;
  rvec *ptr;
  
  /* We introduce a fudge factor for performance reasons: with this choice
   * the initial force on the shells is about a factor of two lower than 
   * without
   */
  fudge = 1.0;
    
  if (bInit) {
    fprintf(log,"RELAX: Using prediction for initial shell placement\n");
    ptr  = x;
    dt_1 = 1;
  }
  else {
    ptr  = v;
    dt_1 = fudge*dt;
  }
    
  for(i=0; (i<ns); i++) {
    s1 = s[i].shell;
    if (bInit)
      clear_rvec(x[s1]);
    switch (s[i].nnucl) {
    case 1:
      n1 = s[i].nucl1;
      for(m=0; (m<DIM); m++)
	x[s1][m]+=ptr[n1][m]*dt_1;
      break;
    case 2:
      n1 = s[i].nucl1;
      n2 = s[i].nucl2;
      m1 = mass[n1];
      m2 = mass[n2];
      tm = dt_1/(m1+m2);
      for(m=0; (m<DIM); m++)
	x[s1][m]+=(m1*ptr[n1][m]+m2*ptr[n2][m])*tm;
      break;
    case 3:
      n1 = s[i].nucl1;
      n2 = s[i].nucl2;
      n3 = s[i].nucl3;
      m1 = mass[n1];
      m2 = mass[n2];
      m3 = mass[n3];
      tm = dt_1/(m1+m2+m3);
      for(m=0; (m<DIM); m++)
	x[s1][m]+=(m1*ptr[n1][m]+m2*ptr[n2][m]+m3*ptr[n3][m])*tm;
      break;
    default:
      gmx_fatal(FARGS,"Shell %d has %d nuclei!",i,s[i].nnucl);
    }
  }
}

static void print_epot(FILE *fp,int mdstep,int count,real epot,real df,
		       int ndir,real sf_dir)
{
  fprintf(fp,"MDStep=%5d/%2d EPot: %12.8e, rmsF: %6.2e",
	  mdstep,count,epot,df);
  if (ndir)
    fprintf(fp,", dir. rmsF: %6.2e\n",sqrt(sf_dir/ndir));
  else
    fprintf(fp,"\n");
}


static real rms_force(t_commrec *cr,rvec f[],int ns,t_shell s[],
		      int ndir,real *sf_dir,real *Epot)
{
  int  i,shell,ntot;
  double buf[4];

  buf[0] = *sf_dir;
  for(i=0; i<ns; i++) {
    shell = s[i].shell;
    buf[0]  += norm2(f[shell]);
  }
  ntot = ns;

  if (PAR(cr)) {
    buf[1] = ntot;
    buf[2] = *sf_dir;
    buf[3] = *Epot;
    gmx_sumd(4,buf,cr);
    ntot = (int)(buf[1] + 0.5);
    *sf_dir = buf[2];
    *Epot   = buf[3];
  }
  ntot += ndir;

  return (ntot ? sqrt(buf[0]/ntot) : 0);
}

static void check_pbc(FILE *fp,rvec x[],int shell)
{
  int m,now;
  
  now = shell-4;
  for(m=0; (m<DIM); m++)
    if (fabs(x[shell][m]-x[now][m]) > 0.3) {
      pr_rvecs(fp,0,"SHELL-X",x+now,5);
      break;
    }
}

static void dump_shells(FILE *fp,rvec x[],rvec f[],real ftol,int ns,t_shell s[])
{
  int  i,shell;
  real ft2,ff2;
  
  ft2 = sqr(ftol);
  
  for(i=0; (i<ns); i++) {
    shell = s[i].shell;
    ff2   = iprod(f[shell],f[shell]);
    if (ff2 > ft2)
      fprintf(fp,"SHELL %5d, force %10.5f  %10.5f  %10.5f, |f| %10.5f\n",
	      shell,f[shell][XX],f[shell][YY],f[shell][ZZ],sqrt(ff2));
    check_pbc(fp,x,shell);
  }
}

static void init_adir(FILE *log,
		      gmx_constr_t constr,t_topology *top,t_inputrec *ir,
		      gmx_domdec_t *dd,
		      int step,t_mdatoms *md,int start,int end,
		      rvec *x_old,rvec *x_init,rvec *x,
		      rvec *f,rvec *acc_dir,matrix box,
		      real lambda,real *dvdlambda,t_nrnb *nrnb)
{
  static rvec *xnold=NULL,*xnew=NULL;
  double w_dt;
  int    gf,ga,gt;
  real   dt,scale;
  int    n,d; 
  unsigned short *ptype;
  rvec   p,dx;
  
  if (xnew == NULL) {
    snew(xnold,end-start);
    snew(xnew,end-start);
  }
    
  ptype = md->ptype;

  dt = ir->delta_t;

  /* Does NOT work with freeze or acceleration groups (yet) */
  for (n=start; n<end; n++) {  
    w_dt = md->invmass[n]*dt;
    
    for (d=0; d<DIM; d++) {
      if ((ptype[n] != eptVSite) && (ptype[n] != eptShell)) {
	xnold[n-start][d] = x[n][d] - (x_init[n][d] - x_old[n][d]);
	xnew[n-start][d] = 2*x[n][d] - x_old[n][d] + f[n][d]*w_dt*dt;
      } else {
	xnold[n-start][d] = x[n][d];
	xnew[n-start][d] = x[n][d];
      }
    }
  }
  constrain(log,FALSE,constr,top,ir,dd,step,md,
	    x,xnold-start,NULL,box,
	    lambda,dvdlambda,dt,NULL,NULL,nrnb,TRUE);
  constrain(log,FALSE,constr,top,ir,dd,step,md,
	    x,xnew-start,NULL,box,
	    lambda,dvdlambda,dt,NULL,NULL,nrnb,TRUE);

  /* Set xnew to minus the acceleration */
  for (n=start; n<end; n++) {
    for(d=0; d<DIM; d++)
      xnew[n-start][d] =
	-(2*x[n][d]-xnold[n-start][d]-xnew[n-start][d])/sqr(dt)
	- f[n][d]*md->invmass[n];
    clear_rvec(acc_dir[n]);
  }

  /* Project the acceleration on the old bond directions */
  constrain(log,FALSE,constr,top,ir,dd,step,md,
	    x_old,xnew-start,acc_dir,box,
	    lambda,dvdlambda,dt,NULL,NULL,nrnb,FALSE); 
}

int relax_shells(FILE *log,t_commrec *cr,bool bVerbose,
		 int mdstep,t_inputrec *inputrec,bool bDoNS,bool bStopCM,
		 t_topology *top,gmx_constr_t constr,
		 real ener[],t_fcdata *fcd,
		 t_state *state,rvec f[],
		 rvec buf[],t_mdatoms *md,
		 t_nrnb *nrnb,gmx_wallcycle_t wcycle,
		 t_graph *graph,t_groups *grps,
		 int nshell,t_shell shells[],
		 t_forcerec *fr,
		 real t,rvec mu_tot,
		 int natoms,bool *bConverged,
		 gmx_vsite_t *vsite,
		 FILE *fp_field)
{
  static bool bFirst=TRUE,bForceInit=FALSE,bNoPredict=FALSE;
  static rvec *pos[2],*force[2];
  static rvec *acc_dir=NULL,*x_old=NULL;
  real   Epot[2],df[2],Estore[F_NRE];
  tensor vir_el_recip[2];
  rvec   dx;
  real   sf_dir,invdt;
  real   ftol,xiH,xiS,dum=0;
  char   cbuf[56];
  bool   bCont,bInit;
  int    i,start=md->start,homenr=md->homenr,end=start+homenr,cg0,cg1;
  int    nflexcon,g,number_steps,d,Min=0,count=0;
#define  Try (1-Min)             /* At start Try = 1 */

  if (bFirst) {
    /* Allocate local arrays */
    for(i=0; (i<2); i++) {
      snew(pos[i],state->natoms);
      snew(force[i],state->natoms);
    }
    bNoPredict = getenv("NOPREDICT") != NULL;
    if (bNoPredict)
      fprintf(log,"Will never predict shell positions\n");
    else {
      bForceInit = getenv("FORCEINIT") != NULL;
      if (bForceInit)
	fprintf(log,"Will always initiate shell positions\n");
    }
    bFirst = FALSE;
  }
  
  bCont        = (mdstep == inputrec->init_step) && inputrec->bContinuation;
  bInit        = (mdstep == inputrec->init_step) || bForceInit;
  ftol         = inputrec->em_tol;
  number_steps = inputrec->niter;
  nflexcon     = n_flexible_constraints(constr);

  if (bDoNS) {
    /* This is the only time where the coordinates are used
     * before do_force is called, which normally puts all
     * charge groups in the box.
     */
    if (PARTDECOMP(cr)) {
      pd_cg_range(cr,&cg0,&cg1);
    } else {
      cg0 = 0;
      cg1 = top->blocks[ebCGS].nr;
    }
    put_charge_groups_in_box(log,cg0,cg1,fr->ePBC,state->box,
			     &(top->blocks[ebCGS]),state->x,fr->cg_cm);
    if (graph)
      mk_mshift(log,graph,fr->ePBC,state->box,state->x);
  }

  /* After this all coordinate arrays will contain whole molecules */
  if (graph)
    shift_self(graph,state->box,state->x);

  if (nflexcon) {
    if (acc_dir == NULL) {
      snew(acc_dir,homenr);
      snew(x_old,homenr);
    }
    for(i=0; i<homenr; i++) {
      for(d=0; d<DIM; d++)
        x_old[i][d] =
	  state->x[start+i][d] - state->v[start+i][d]*inputrec->delta_t;
    }
  }
  
  /* Do a prediction of the shell positions */
  if (!bNoPredict && !bCont)
    predict_shells(log,state->x,state->v,inputrec->delta_t,nshell,shells,
		   md->massT,bInit);

  /* do_force expected the charge groups to be in the box */
  if (graph)
    unshift_self(graph,state->box,state->x);

  /* Calculate the forces first time around */
  if (gmx_debug_at) {
    pr_rvecs(debug,0,"x b4 do_force",state->x + start,homenr);
  }
  do_force(log,cr,inputrec,mdstep,nrnb,wcycle,top,grps,
	   state->box,state->x,force[Min],buf,md,ener,fcd,
	   state->lambda,graph,
	   TRUE,bDoNS,FALSE,TRUE,fr,mu_tot,FALSE,t,fp_field,NULL);
  sum_lrforces(force[Min],fr,start,homenr);
  copy_mat(fr->vir_el_recip,vir_el_recip[Min]);

  sf_dir = 0;
  if (nflexcon) {
    init_adir(log,constr,top,inputrec,cr->dd,mdstep,md,start,end,
	      x_old-start,state->x,state->x,force[Min],acc_dir-start,
	      state->box,state->lambda,&dum,nrnb);

    for(i=start; i<end; i++)
      sf_dir += md->massT[i]*norm2(acc_dir[i-start]);
  }

  /* Sum the potential energy terms from group contributions */
  sum_epot(&(inputrec->opts),grps,ener);
  Epot[Min]=ener[F_EPOT];

  df[Min]=rms_force(cr,force[Min],nshell,shells,nflexcon,&sf_dir,&Epot[Min]);
  df[Try]=0;
  if (debug) {
    fprintf(debug,"df = %g  %g\n",df[Min],df[Try]);
  }

  if (gmx_debug_at) {
    pr_rvecs(debug,0,"force0",force[Min],md->nr);
  }

  if (nshell+nflexcon > 0) {
    /* Copy x to pos[Min] & pos[Try]: during minimization only the
     * shell positions are updated, therefore the other particles must
     * be set here.
     */
    memcpy(pos[Min],state->x,state->natoms*sizeof(state->x[0]));
    memcpy(pos[Try],state->x,state->natoms*sizeof(state->x[0]));
  }
  
  if (bVerbose && MASTER(cr))
    print_epot(stdout,mdstep,0,Epot[Min],df[Min],nflexcon,sf_dir);

  if (debug) {
    fprintf(debug,"%17s: %14.10e\n",
	    interaction_function[F_EKIN].longname, ener[F_EKIN]);
    fprintf(debug,"%17s: %14.10e\n",
	    interaction_function[F_EPOT].longname, ener[F_EPOT]);
    fprintf(debug,"%17s: %14.10e\n",
	    interaction_function[F_ETOT].longname, ener[F_ETOT]);
    fprintf(debug,"SHELLSTEP %d\n",mdstep);
  }
  
  /* First check whether we should do shells, or whether the force is 
   * low enough even without minimization.
   */
  *bConverged = (df[Min] < ftol);
  
  for(count=1; (!(*bConverged) && (count < number_steps)); count++) {
    if (vsite)
      construct_vsites(log,vsite,pos[Min],nrnb,inputrec->delta_t,state->v,
		       &top->idef,fr->ePBC,fr->bMolPBC,graph,cr,state->box);
     
    if (nflexcon) {
      init_adir(log,constr,top,inputrec,cr->dd,mdstep,md,start,end,
		x_old-start,state->x,pos[Min],force[Min],acc_dir-start,
		state->box,state->lambda,&dum,nrnb);
      
      directional_sd(log,pos[Min],pos[Try],acc_dir-start,start,end,
		     fr->fc_stepsize);
    }
    
    /* New positions, Steepest descent */
    shell_pos_sd(log,pos[Min],pos[Try],force[Min],nshell,shells,count); 

    /* do_force expected the charge groups to be in the box */
    if (graph)
      unshift_self(graph,state->box,pos[Try]);

    if (gmx_debug_at) {
      pr_rvecs(debug,0,"RELAX: pos[Min]  ",pos[Min] + start,homenr);
      pr_rvecs(debug,0,"RELAX: pos[Try]  ",pos[Try] + start,homenr);
    }
    /* Try the new positions */
    do_force(log,cr,inputrec,1,nrnb,wcycle,
	     top,grps,state->box,pos[Try],force[Try],buf,md,ener,fcd,
	     state->lambda,graph,
	     TRUE,FALSE,FALSE,TRUE,fr,mu_tot,FALSE,t,fp_field,NULL);
    if (vsite) 
      spread_vsite_f(log,vsite,pos[Try],force[Try],fr->fshift,nrnb,
		     &top->idef,fr->ePBC,fr->bMolPBC,graph,state->box,cr);
      
    /* Calculation of the virial must be done after vsites!    */
    /* Question: Is it correct to do the PME forces after this? */
    /*    calc_virial(log,START(nsb),HOMENR(nsb),pos[Try],force[Try],
		my_vir[Try],pme_vir[Try],graph,state->box,nrnb,fr,FALSE);
    */	  
    /* Spread the LR force on virtual site to the other particles... 
     * This is parallellized. MPI communication is performed
     * if the constructing atoms aren't local.
     */
    if (vsite && fr->bEwald) 
      spread_vsite_f(log,vsite,pos[Try],fr->f_el_recip,NULL,nrnb,
		     &top->idef,fr->ePBC,fr->bMolPBC,graph,state->box,cr);
    
    sum_lrforces(force[Try],fr,start,homenr);
    copy_mat(fr->vir_el_recip,vir_el_recip[Try]);
    
    if (gmx_debug_at) {
      pr_rvecs(debug,0,"RELAX: force[Min]",force[Min] + start,homenr);
      pr_rvecs(debug,0,"RELAX: force[Try]",force[Try] + start,homenr);
    }
    sf_dir = 0;
    if (nflexcon) {
      init_adir(log,constr,top,inputrec,cr->dd,mdstep,md,start,end,
		x_old-start,state->x,pos[Try],force[Try],acc_dir-start,
		state->box,state->lambda,&dum,nrnb);

      for(i=start; i<end; i++)
	sf_dir += md->massT[i]*norm2(acc_dir[i-start]);
    }

    /* Sum the potential energy terms from group contributions */
    sum_epot(&(inputrec->opts),grps,ener);
    Epot[Try]=ener[F_EPOT]; 
    
    df[Try]=rms_force(cr,force[Try],nshell,shells,nflexcon,&sf_dir,&Epot[Try]);

    if (debug)
      fprintf(debug,"df = %g  %g\n",df[Min],df[Try]);

    if (debug) {
      if (gmx_debug_at)
	pr_rvecs(debug,0,"F na do_force",force[Try] + start,homenr);
      if (gmx_debug_at) {
	fprintf(debug,"SHELL ITER %d\n",count);
	dump_shells(debug,pos[Try],force[Try],ftol,nshell,shells);
      }
    }

    if (bVerbose && MASTER(cr))
      print_epot(stdout,mdstep,count,Epot[Try],df[Try],nflexcon,sf_dir);
      
    *bConverged = (df[Try] < ftol);
    
    if ((df[Try] < df[Min])) {
      if (debug)
	fprintf(debug,"Swapping Min and Try\n");
      if (nflexcon) {
	/* Correct the velocities for the flexible constraints */
	invdt = 1/inputrec->delta_t;
	for(i=start; i<end; i++)
	  for(d=0; d<DIM; d++)
	    state->v[i][d] += (pos[Try][i][d] - pos[Min][i][d])*invdt;
      }
      Min  = Try;
    } else {
      decrease_step_size(nshell,shells);
    }
  }
  if (MASTER(cr) && !(*bConverged)) {
    fprintf(log,
	    "step %d: EM did not converge in %d iterations, RMS force %.3f\n",
	    mdstep,number_steps,df[Min]);
    fprintf(stderr,
	    "step %d: EM did not converge in %d iterations, RMS force %.3f\n",
	    mdstep,number_steps,df[Min]);
  }

  /* Parallelise this one! */
  if (EEL_FULL(fr->eeltype)) {
    for(i=start; (i<end); i++)
      rvec_dec(force[Min][i],fr->f_el_recip[i]);
  }
  memcpy(f,force[Min],state->natoms*sizeof(f[0]));

  /* CHECK VIRIAL */
  copy_mat(vir_el_recip[Min],fr->vir_el_recip);
  
  memcpy(state->x,pos[Min],state->natoms*sizeof(state->x[0]));

  return count; 
}

