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
 * Gromacs Runs One Microsecond At Cannonball Speeds
 */
static char *SRCID_relax_sh_c = "$Id$";
#include <string.h>
#include "assert.h"
#include "typedefs.h"
#include "smalloc.h"
#include "fatal.h"
#include "vec.h"
#include "txtdump.h"
#include "mdrun.h"
#include "xmdrun.h"
#include "mdatoms.h"
#include "network.h"
#include "names.h"
#include "constr.h"

static void do_1pos(rvec xnew,rvec xold,rvec f,real k_1,real step)
{
  real xo,yo,zo;
  real dx,dy,dz,dx2;
  
  xo=xold[XX];
  yo=xold[YY];
  zo=xold[ZZ];
  
  dx=f[XX]*k_1;
  dy=f[YY]*k_1;
  dz=f[ZZ]*k_1;
  
  xnew[XX]=xo+dx*step;
  xnew[YY]=yo+dy*step;
  xnew[ZZ]=zo+dz*step;
}

static void directional_sd(FILE *log,real step,rvec xold[],rvec xnew[],
			   rvec acc_dir[],int start,int homenr,real k)
{
  real invk;
  int  i;
  
  invk = 1.0/k;

  for(i=start; i<homenr; i++)
    do_1pos(xnew[i],xold[i],acc_dir[i],invk,step);
}

static void shell_pos_sd(FILE *log,real step,rvec xold[],rvec xnew[],rvec f[],
			 int ns,t_shell s[])
{
  int  i,shell;
  real k_1;
  
  for(i=0; (i<ns); i++) {
    shell = s[i].shell;
    k_1   = s[i].k_1;
    do_1pos(xnew[shell],xold[shell],f[shell],k_1,step);
    if (debug && 0) {
      pr_rvec(debug,0,"fshell",f[shell],DIM);
      pr_rvec(debug,0,"xold",xold[shell],DIM);
      pr_rvec(debug,0,"xnew",xnew[shell],DIM);
    }
  }
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
      fatal_error(0,"Shell %d has %d nuclei!",i,s[i].nnucl);
    }
  }
}

static void print_epot(FILE *fp,int mdstep,int count,real step,real epot,
		       real df,bool bLast)
{
  fprintf(fp,"MDStep=%5d/%2d lamb: %6g, EPot: %12.8e",
	  mdstep,count,step,epot);
  
  if (count != 0)
    fprintf(fp,", rmsF: %12.8e\n",df);
  else
    fprintf(fp,"\n");
}


static real rms_force(t_commrec *cr,rvec f[],int ns,t_shell s[],
		      int ndir,real sf_dir)
{
  int  i,shell,ntot;
  real df2;
  
  ntot = ns+ndir;

  if (!ntot)
    return 0;
  df2 = sf_dir;
  for(i=0; i<ns; i++) {
    shell = s[i].shell;
    df2  += iprod(f[shell],f[shell]);
  }
  if (PAR(cr)) {
    gmx_sum(1,&df2,cr);
    gmx_sumi(1,&ntot,cr);
  }
  return sqrt(df2/ntot);
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

static int count_zero_length_constraints(t_idef *idef)
{
  int nZeroLen,i;

  nZeroLen = 0;

  for(i=0; i<idef->il[F_SHAKE].nr; i+=3)
    if (idef->iparams[idef->il[F_SHAKE].iatoms[i]].shake.dA == 0)
      nZeroLen++;
  
  return nZeroLen;
}

static void init_adir(FILE *log,t_topology *top,t_inputrec *ir,int step,
		      t_mdatoms *md,int start,int end,
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
      if ((ptype[n] != eptDummy) && (ptype[n] != eptShell)) {
	xnold[start+n][d] = x[n][d] - (x_init[n][d] - x_old[n][d]);
	xnew[start+n][d] = 2*x[n][d] - x_old[n][d] + f[n][d]*w_dt*dt;
      } else {
	xnold[start+n][d] = x[n][d];
	xnew[start+n][d] = x[n][d];
      }
    }
  }
  constrain(log,top,ir,step,md,start,end,
	    x,xnold-start,NULL,box,
	    lambda,dvdlambda,nrnb,TRUE);
  constrain(log,top,ir,step,md,start,end,
	    x,xnew-start,NULL,box,
	    lambda,dvdlambda,nrnb,TRUE);

  /* Set xnew to minus the acceleration */
  for (n=start; n<end; n++) {
    for(d=0; d<DIM; d++)
      xnew[n-start][d] = -(2*x[n][d]-xnold[n][d]-xnew[n][d])/sqr(dt)
	- f[n][d]*md->invmass[n];
    clear_rvec(acc_dir[n]);
  }

  /* Project the accereration on the old bond directions */
  constrain(log,top,ir,step,md,start,end,
	    x_old,xnew,acc_dir,box,
	    lambda,dvdlambda,nrnb,FALSE); 
}

int relax_shells(FILE *log,t_commrec *cr,t_commrec *mcr,bool bVerbose,
		 int mdstep,t_parm *parm,bool bDoNS,bool bStopCM,
		 t_topology *top,real ener[],t_fcdata *fcd,
		 rvec x[],rvec vold[],rvec v[],rvec vt[],rvec f[],
		 rvec buf[],t_mdatoms *md,t_nsborder *nsb,t_nrnb *nrnb,
		 t_graph *graph,t_groups *grps,tensor vir_part,
		 tensor pme_vir_part,bool bShell,
		 int nshell,t_shell shells[],t_forcerec *fr,
		 char *traj,real t,real lambda,rvec mu_tot,
		 int natoms,matrix box,bool *bConverged)
{
  static bool bFirst=TRUE,bInit;
  static rvec *pos[2],*force[2];
  static rvec *acc_dir=NULL,*x_old=NULL;
  static int  ndir;
  real   Epot[2],df[2],Estore[F_NRE];
  tensor my_vir[2],vir_last,pme_vir[2];
  rvec   dx;
  real   sf_dir;
#define NEPOT asize(Epot)
  real   ftol,step,step0,xiH,xiS,dum=0;
  char   cbuf[56];
  bool   bDone;
  int    i,start=START(nsb),homenr=HOMENR(nsb),end=START(nsb)+HOMENR(nsb);
  int    g,number_steps,d,Min=0,count=0;
#define  Try (1-Min)             /* At start Try = 1 */

  if (bFirst) {
    /* Check for directional minimization */
    if (fr->fc_stepsize != 0)
      ndir = count_zero_length_constraints(&(top->idef));
    else
      ndir = 0;
    /* Allocate local arrays */
    if (bShell) {
      for(i=0; (i<2); i++) {
	snew(pos[i],nsb->natoms);
	snew(force[i],nsb->natoms);
      }
    }
    else {
      /* Copy pointers */
      pos[Min]   = x;
      force[Min] = f;
    }
    bInit  = (getenv("FORCEINIT") != NULL);
    bFirst = FALSE;
  }
  
  ftol         = parm->ir.em_tol;
  number_steps = parm->ir.niter;
  step0        = 1.0;

  if (ndir) {
    if (acc_dir == NULL) {
      snew(acc_dir,homenr);
      snew(x_old,homenr);
    }
    init_pbc(parm->box);
    for(i=0; i<homenr; i++) {
      for(d=0; d<DIM; d++)
        x_old[i][d] = x[start+i][d] - v[start+i][d]*parm->ir.delta_t;
    }
  }
  
  /* Do a prediction of the shell positions */
  predict_shells(log,x,v,parm->ir.delta_t,nshell,shells,
		 md->massT,bInit || (mdstep == 0));
   
  /* Calculate the forces first time around */
  clear_mat(my_vir[Min]);
  clear_mat(pme_vir[Min]);
  do_force(log,cr,mcr,parm,nsb,my_vir[Min],pme_vir[Min],mdstep,nrnb,
	   top,grps,x,v,force[Min],buf,md,ener,fcd,bVerbose && !PAR(cr),
	   lambda,graph,bDoNS,FALSE,fr,mu_tot,FALSE,t);
  sum_lrforces(force[Min],fr,start,homenr);

  sf_dir = 0;
  if (ndir) {
    init_adir(log,top,&(parm->ir),mdstep,md,start,end,
	      x_old-start,x,x,force[Min],acc_dir-start,
	      parm->box,lambda,&dum,nrnb);

    for(i=start; i<end; i++)
      sf_dir += md->massT[i]*norm2(acc_dir[i-start]);

    if (bVerbose)
      fprintf(stderr,"RMS dir. force: %g\n",sqrt(sf_dir/ndir));
  }
  
  df[Min]=rms_force(cr,force[Min],nshell,shells,ndir,sf_dir);
  df[Try]=0;
  if (debug) {
    fprintf(debug,"df = %g  %g\n",df[Min],df[Try]);
    sprintf(cbuf,"myvir step %d",0);
    pr_rvecs(debug,0,cbuf,my_vir[Min],DIM);
  }
    
#ifdef DEBUG
  if (debug) {
    pr_rvecs(debug,0,"force0",force[Min],md->nr);
  }
#endif
  if (nshell+ndir > 0) {
    /* Copy x to pos[Min] & pos[Try]: during minimization only the
     * shell positions are updated, therefore the other particles must
     * be set here.
     */
    memcpy(pos[Min],x,nsb->natoms*sizeof(x[0]));
    memcpy(pos[Try],x,nsb->natoms*sizeof(x[0]));
  }
  /* Sum the potential energy terms from group contributions */
  sum_epot(&(parm->ir.opts),grps,ener);
  Epot[Min]=ener[F_EPOT];

  if (PAR(cr))
    gmx_sum(NEPOT,Epot,cr);
  
  step=step0;
  
  if (bVerbose && MASTER(cr) && (nshell+ndir > 0))
    print_epot(stdout,mdstep,0,step,Epot[Min],df[Min],FALSE);

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
  *bConverged = bDone = (df[Min] < ftol) || (nshell+ndir == 0);
  
  for(count=1; (!bDone && (count < number_steps)); count++) {

    if (ndir) {
      init_adir(log,top,&(parm->ir),mdstep,md,start,end,
		x_old-start,x,pos[Min],force[Min],acc_dir-start,
		parm->box,lambda,&dum,nrnb);
      
      directional_sd(log,step,pos[Min],pos[Try],acc_dir-start,start,end,
		     fr->fc_stepsize);
    }
    
    /* New positions, Steepest descent */
    shell_pos_sd(log,step,pos[Min],pos[Try],force[Min],nshell,shells); 

    if (debug) {
      pr_rvecs(debug,0,"pos[Try] b4 do_force",pos[Try] + start,homenr);
      pr_rvecs(debug,0,"pos[Min] b4 do_force",pos[Min] + start,homenr);
    }
    /* Try the new positions */
    clear_mat(my_vir[Try]);
    clear_mat(pme_vir[Try]);
    do_force(log,cr,mcr,parm,nsb,my_vir[Try],pme_vir[Try],1,nrnb,
	     top,grps,pos[Try],v,force[Try],buf,md,ener,fcd,
	     bVerbose && !PAR(cr),
	     lambda,graph,FALSE,FALSE,fr,mu_tot,FALSE,t);
    sum_lrforces(force[Try],fr,start,homenr);
    sf_dir = 0;
    if (ndir) {
      init_adir(log,top,&(parm->ir),mdstep,md,start,end,
		x_old-start,x,pos[Try],force[Try],acc_dir-start,
		parm->box,lambda,&dum,nrnb);

      for(i=start; i<end; i++)
	sf_dir += md->massT[i]*norm2(acc_dir[i-start]);

      if (bVerbose)
	fprintf(stderr,"dir. rmsf %g\n",sqrt(sf_dir/ndir));
    }
    df[Try]=rms_force(cr,force[Try],nshell,shells,ndir,sf_dir);

    if (debug)
      fprintf(debug,"df = %g  %g\n",df[Min],df[Try]);

    if (debug) {
      /*pr_rvecs(debug,0,"F na do_force",force[Try] + start,homenr);*/
      sprintf(cbuf,"myvir step %d",count);
      pr_rvecs(debug,0,cbuf,my_vir[Try],DIM);
      /*fprintf(debug,"SHELL ITER %d\n",count);
	dump_shells(debug,pos[Try],force[Try],ftol,nshell,shells);*/
    }
    /* Sum the potential energy terms from group contributions */
    sum_epot(&(parm->ir.opts),grps,ener);
    Epot[Try]=ener[F_EPOT];

    if (PAR(cr)) 
      gmx_sum(1,&Epot[Try],cr);

    if (bVerbose && MASTER(cr))
      print_epot(stdout,mdstep,count,step,Epot[Try],df[Try],FALSE);
      
    *bConverged = (df[Try] < ftol);
    bDone       = *bConverged || (step < 0.01);
    
    /* if ((Epot[Try] < Epot[Min])) { */
    if ((df[Try] < df[Min])) {
      if (debug)
	fprintf(debug,"Swapping Min and Try\n");
      Min  = Try;
      step = step0;
    }
    else
      step *= 0.8;
  }
  if (MASTER(cr) && !bDone) 
    fprintf(stderr,"EM did not converge in %d steps\n",number_steps);

  /* Parallelise this one! */
  if (EEL_LR(fr->eeltype)) {
    for(i=start; (i<end); i++)
      rvec_dec(force[Min][i],fr->f_pme[i]);
  }
  memcpy(f,force[Min],nsb->natoms*sizeof(f[0]));

  /* CHECK VIRIAL */
  copy_mat(my_vir[Min],vir_part);
  copy_mat(pme_vir[Min],pme_vir_part);
  
  if (debug) {
    sprintf(cbuf,"myvir step %d",count);
    pr_rvecs(debug,0,cbuf,vir_part,DIM);
  }

  if (nshell+ndir > 0)
    memcpy(x,pos[Min],nsb->natoms*sizeof(x[0]));
  if (ndir > 0) {
    constrain(log,top,&(parm->ir),mdstep,md,start,end,
	      x-start,x_old-start,NULL,parm->box,
	      lambda,&dum,nrnb,TRUE);
    for(i=0; i<homenr; i++) {
      pbc_dx(x[start+i],x_old[i],dx);
      svmul(1/parm->ir.delta_t,dx,v[start+i]);
    }
  }

  return count; 
}

