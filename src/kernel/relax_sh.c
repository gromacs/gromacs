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
 * Gromacs Runs On Most of All Computer Systems
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
#include "init_sh.h"
#include "mdatoms.h"
#include "network.h"
#include "do_gct.h"
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
			   rvec min_f[],int start,int homenr,real k)
{
  int  i;
  
  for(i=start; i<homenr; i++)
    do_1pos(xnew[i],xold[i],min_f[i],-k,step);
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

static void zero_shell_forces(FILE *log,rvec f[],int ns,t_shell s[])
{
  int i,s1,n1,n2;
  real m1,m2,m3,tm;
  
  /* Subtract remaining forces on shells from the attaching atoms.
   * This way energy conservation may get better, and it is not unreasonable:
   * It corresponds to an ad-hoc change in the force constants, which is exactly
   * large enough to compensate for the remaining force on the shell. The
   * shell forces then becomes zero.
   */  
  for(i=0; (i<ns); i++) {
    s1 = s[i].shell;
    switch (s[i].nnucl) {
    case 1:
      n1 = s[i].nucl1;
      rvec_dec(f[n1],f[s1]);
      clear_rvec(f[s1]);
      break;
    case 2:
      n1 = s[i].nucl1;
      n2 = s[i].nucl2;
      /*m1 = mass[n1];
	m2 = mass[n2];
	tm = dt_1/(m1+m2);*/
      break;
    case 3:
      n1 = s[i].nucl1;
      n2 = s[i].nucl2;
      /*n3 = s[i].nucl3;
      m1 = mass[n1];
      m2 = mass[n2];
      m3 = mass[n3];
      tm = dt_1/(m1+m2+m3);*/
      break;
    default:
      fatal_error(0,"Shell %d has %d nuclei!",i,s[i].nnucl);
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


static real rms_force(t_commrec *cr,rvec f[],int ns,t_shell s[])
{
  int  i,shell;
  real df2;
  
  if (!ns)
    return 0;
  df2=0.0;
  for(i=0; (i<ns); i++) {
    shell = s[i].shell;
    df2  += iprod(f[shell],f[shell]);
  }
  if (PAR(cr)) {
    gmx_sum(1,&df2,cr);
    gmx_sumi(1,&ns,cr);
  }
  return sqrt(df2/ns);
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

int relax_shells(FILE *log,t_commrec *cr,bool bVerbose,
		 int mdstep,t_parm *parm,bool bDoNS,bool bStopCM,
		 t_topology *top,real ener[],
		 rvec x[],rvec vold[],rvec v[],rvec vt[],rvec f[],
		 rvec buf[],t_mdatoms *md,t_nsborder *nsb,t_nrnb *nrnb,
		 t_graph *graph,t_groups *grps,tensor vir_part,
		 tensor pme_vir_part,
		 int nshell,t_shell shells[],t_forcerec *fr,
		 char *traj,real t,real lambda,rvec mu_tot,
		 int natoms,matrix box,bool *bConverged)
{
  static bool bFirst=TRUE,bInit;
  static rvec *pos[2],*force[2];
  real   Epot[2],df[2],Estore[F_NRE];
  tensor my_vir[2],vir_last,pme_vir[2];
  rvec   *min_f_dir=NULL;
#define NEPOT asize(Epot)
  real   ftol,step,step0,xiH,xiS,dum=0;
  char   cbuf[56];
  bool   bDone;
  int    g;
  int    number_steps;
  int    count=0;
  int    i,start=START(nsb),homenr=HOMENR(nsb),end=START(nsb)+HOMENR(nsb);
  int    Min=0;
#define  Try (1-Min)             /* At start Try = 1 */

  if (bFirst) {
    /* Allocate local arrays */
    for(i=0; (i<2); i++) {
      if (nshell > 0)
	snew(pos[i],nsb->natoms);
      snew(force[i],nsb->natoms);
    }
    bInit  = (getenv("FORCEINIT") != NULL);
    bFirst = FALSE;
  }
  
  ftol         = parm->ir.em_tol;
  number_steps = parm->ir.niter;
  step0        = 1.0;

  /* Do a prediction of the shell positions */
  predict_shells(log,x,v,parm->ir.delta_t,nshell,shells,
		 md->massT,bInit || (mdstep == 0));
   
  /* Calculate the forces first time around */
  clear_mat(my_vir[Min]);
  clear_mat(pme_vir[Min]);
  do_force(log,cr,parm,nsb,my_vir[Min],pme_vir[Min],mdstep,nrnb,
	   top,grps,x,v,force[Min],buf,md,ener,bVerbose && !PAR(cr),
	   lambda,graph,bDoNS,FALSE,fr,mu_tot,FALSE);
  sum_lrforces(force[Min],fr,start,homenr);

  df[Min]=rms_force(cr,force[Min],nshell,shells);
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
  if (nshell > 0) {
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
  
  if (bVerbose && MASTER(cr) && (nshell > 0))
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
  *bConverged = bDone = (df[Min] < ftol) || (nshell == 0);
  
  for(count=1; (!bDone && (count < number_steps)); count++) {
    
    /* New positions, Steepest descent */
    shell_pos_sd(log,step,pos[Min],pos[Try],force[Min],nshell,shells); 

    if (fr->k_dirmin) {
      if (min_f_dir == NULL)
	snew(min_f_dir,homenr);
      constrain(log,top,&(parm->ir),mdstep,md,start,end,
		pos[Min],force[Min],min_f_dir,parm->box,
		lambda,&dum,nrnb,FALSE);
      directional_sd(log,step,pos[Min],pos[Try],min_f_dir,start,homenr,
		     fr->k_dirmin);
    }

    if (debug) {
      pr_rvecs(debug,0,"pos[Try] b4 do_force",pos[Try] + start,homenr);
      pr_rvecs(debug,0,"pos[Min] b4 do_force",pos[Min] + start,homenr);
    }
    /* Try the new positions */
    clear_mat(my_vir[Try]);
    clear_mat(pme_vir[Try]);
    do_force(log,cr,parm,nsb,my_vir[Try],pme_vir[Try],1,nrnb,
	     top,grps,pos[Try],v,force[Try],buf,md,ener,bVerbose && !PAR(cr),
	     lambda,graph,FALSE,FALSE,fr,mu_tot,FALSE);
    sum_lrforces(force[Try],fr,start,homenr);
    df[Try]=rms_force(cr,force[Try],nshell,shells);
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

  /*if (bDone)
    zero_shell_forces(log,force[Min],nshell,shells);
  */
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
  if (nshell > 0)
    memcpy(x,pos[Min],nsb->natoms*sizeof(x[0]));
    
  return count; 
}

