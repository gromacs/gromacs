#include <stdio.h>
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

static void shell_pos_sd(FILE *log,real step,rvec xold[],rvec xnew[],rvec f[],
			 int ns,t_shell s[])
{
  int  i,shell;
  real k_1;
  
  for(i=0; (i<ns); i++) {
    shell = s[i].shell;
    k_1   = s[i].k_1;
    do_1pos(xnew[shell],xold[shell],f[shell],k_1,step);
    if (debug) {
      pr_rvec(debug,0,"fshell",f[shell],DIM);
      pr_rvec(debug,0,"xold",xold[shell],DIM);
      pr_rvec(debug,0,"xnew",xnew[shell],DIM);
    }
  }
}

static void predict_shells(FILE *log,rvec x[],rvec v[],real dt,int ns,t_shell s[],
			   real mass[],bool bInit)
{
  int  i,m,s1,n1,n2,n3;
  real dt_1,dt_2,dt_3,fudge,tm,m1,m2,m3;
  rvec *ptr;
  
  /* We introduce a fudge factor for performance reasons: with this choice
   * the initial force on the shells is about a factor of two lower than without
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

static void print_epot(int mdstep,int count,real step,real epot,real df,
		       bool bOptim,bool bLast)
{
  if (bOptim) {
    if (bLast)
      fprintf(stderr,"MDStep=%5d    EPot: %12.8e\n",mdstep,epot);
    else if (count > 0)
      fprintf(stderr,"MDStep=%5d/%2d lamb: %6g, RMSF: %12.8e\n",
	      mdstep,count,step,df);
  }
  else {
    fprintf(stderr,"MDStep=%5d/%2d lamb: %6g, EPot: %12.8e",
	    mdstep,count,step,epot);
    
    if (count != 0)
      fprintf(stderr,", rmsF: %12.8e\n",df);
    else
      fprintf(stderr,"\n");
  }
}


static real rms_force(rvec f[],int ns,t_shell s[])
{
  int  i,shell;
  real df2;
  
  df2=0.0;
  for(i=0; (i<ns); i++) {
    shell = s[i].shell;
    df2  += iprod(f[shell],f[shell]);
  }
  return sqrt(df2/ns);
}

static void check_pbc(FILE *fp,rvec x[],int shell)
{
  int m;
  
  for(m=0; (m<DIM); m++)
    if (fabs(x[shell][m]-x[shell-3][m]) > 0.3)
      pr_rvecs(fp,0,"SHELL-X",x+(4*(shell / 4)),4);
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

void set_nbfmask(t_forcerec *fr,t_mdatoms *md,int ngrp,bool bMask)
{
  static bool *bShell=NULL;
  int    i,j,gid;

  /* THIS CODE IS OUT OF ORDER */
  return;
    
  if (!bShell) {
    snew(bShell,ngrp);
    /* Loop over atoms */
    for(i=0; (i<md->nr); i++) 
      /* If it is a shell, then set the value to TRUE */
      if (md->ptype[i] == eptShell) {
	gid = md->cENER[i];
	assert(gid < ngrp);
	assert(gid >= 0);
	bShell[gid] = TRUE;
      }
  }
  
  /* Loop over all the energy groups */
  for(i=0; (i<ngrp); i++)
    for(j=i; (j<ngrp); j++) {
      gid = GID(i,j,ngrp);
      /* If either of the groups denote a group with shells, then */
      /*
          if (bShell[i] || bShell[j])
	fr->bMask[gid] = bMask;
      else
	fr->bMask[gid] = !bMask;
      */
    }
}

int relax_shells(FILE *log,t_commrec *cr,bool bVerbose,
		 int mdstep,t_parm *parm,bool bDoNS,bool bStopCM,
		 t_topology *top,real ener[],
		 rvec x[],rvec vold[],rvec v[],rvec vt[],rvec f[],
		 rvec buf[],t_mdatoms *md,t_nsborder *nsb,t_nrnb *nrnb,
		 t_graph *graph,t_groups *grps,tensor vir_part,
		 int nshell,t_shell shells[],t_forcerec *fr,
		 char *traj,real t,real lambda,
		 int natoms,matrix box,t_mdebin *mdebin,bool *bConverged)
{
  static bool bFirst=TRUE,bOptim=FALSE;
  static rvec *pos[2],*force[2];
  real   Epot[2],df[2],Estore[F_NRE];
  tensor my_vir[2],vir_last;
#define NEPOT asize(Epot)
  real   ftol,step,step0,xiH,xiS;
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
      snew(pos[i],nsb->natoms);
      snew(force[i],nsb->natoms);
    }
    bOptim = (getenv("NO_SHELL_OPTIM") == NULL);
    fprintf(log,"bOptim = %s\n",bool_names[bOptim]);
    bFirst = FALSE;
  }
  
  ftol         = parm->ir.em_tol;
  number_steps = parm->ir.niter;
  step0        = 1.0;

  /* Do a prediction of the shell positions */
  predict_shells(log,x,v,parm->ir.delta_t,nshell,shells,md->massT,(mdstep == 0));
   
  /* Calculate the forces first time around */
  if (bOptim)
    set_nbfmask(fr,md,parm->ir.opts.ngener,TRUE);
  clear_mat(my_vir[Min]);
  do_force(log,cr,parm,nsb,my_vir[Min],mdstep,nrnb,
	   top,grps,x,v,force[Min],buf,md,ener,bVerbose && !PAR(cr),
	   lambda,graph,bDoNS,FALSE,fr);
  df[Min]=rms_force(force[Min],nshell,shells);
#ifdef DEBUG
  if (debug) {
    pr_rvecs(debug,0,"force0",force[Min],md->nr);
    calc_f_dev(md->nr,md->chargeA,x,force[Min],&top->idef,&xiH,&xiS);
    fprintf(log,"xiH = %e, xiS = %e\n",xiH,xiS);
  }
#endif
  /* Copy x to pos[Min] & pos[Try]: during minimization only the
   * shell positions are updated, therefore the other particles must
   * be set here.
   */
  memcpy(pos[Min],x,nsb->natoms*sizeof(x[0]));
  memcpy(pos[Try],x,nsb->natoms*sizeof(x[0]));

  /* Sum the potential energy terms from group contributions */
  sum_epot(&(parm->ir.opts),grps,ener);
  Epot[Min]=ener[F_EPOT];
  if (PAR(cr))
    gmx_sum(NEPOT,Epot,cr);
  
  step=step0;
  if (bVerbose && MASTER(cr) && (nshell > 0))
    print_epot(mdstep,0,step,Epot[Min],df[Min],bOptim,FALSE);

  if (debug) {
    fprintf(debug,"%17s: %14.10e\n",
	    interaction_function[F_EKIN].longname, ener[F_EKIN]);
    fprintf(debug,"%17s: %14.10e\n",
	    interaction_function[F_EPOT].longname, ener[F_EPOT]);
    fprintf(debug,"%17s: %14.10e\n",
	    interaction_function[F_ETOT].longname, ener[F_ETOT]);
  }

  if (debug)
    fprintf(debug,"SHELLSTEP %d\n",mdstep);

  /* First check whether we should do shells, or whether the force is low enough
   * even without minimization.
   */
  *bConverged = bDone = (df[Min] < ftol) || (nshell == 0);
  
  for(count=1; (!bDone && (count < number_steps)); count++) {
    
    /* New positions, Steepest descent */
    shell_pos_sd(log,step,pos[Min],pos[Try],force[Min],nshell,shells); 

#ifdef DEBUG
    if (debug) {
      pr_rvecs(debug,0,"pos[Try] b4 do_force",(rvec *)(pos[Try][start]),
	       homenr);
      pr_rvecs(debug,0,"pos[Min] b4 do_force",(rvec *)(pos[Min][start]),
	       homenr);
    }
#endif
    /* Try the new positions */
    clear_mat(my_vir[Try]);
    do_force(log,cr,parm,nsb,my_vir[Try],1,nrnb,
	     top,grps,pos[Try],v,force[Try],buf,md,ener,bVerbose && !PAR(cr),
	     lambda,graph,FALSE,FALSE,fr);
    df[Try]=rms_force(force[Try],nshell,shells);

    if (debug) 
      pr_rvecs(debug,0,"F na do_force",(rvec *)(force[Try][start]),homenr);

    if (debug) {
      fprintf(debug,"SHELL ITER %d\n",count);
      dump_shells(debug,pos[Try],force[Try],ftol,nshell,shells);
    }
    /* Sum the potential energy terms from group contributions */
    sum_epot(&(parm->ir.opts),grps,ener);
    Epot[Try]=ener[F_EPOT];
    if (PAR(cr)) 
      gmx_sum(NEPOT,Epot,cr);
    
    if (bVerbose && MASTER(cr))
      print_epot(mdstep,count,step,Epot[Try],df[Try],bOptim,FALSE);

    *bConverged = (df[Try] < ftol);
    bDone       = *bConverged || (step < 0.5);
    
    if ((Epot[Try] < Epot[Min])) {
      Min  = Try;
      step = step0;
    }
    else
      step *= 0.8;
  }
  if (MASTER(cr) && !bDone) 
    fprintf(stderr,"EM did not converge in %d steps\n",number_steps);
  
  /* Now compute the forces on the other particles (atom-atom) */
  if (bOptim) {
    /* Store the old energies */
    for(i=0; (i<F_NRE); i++)
      Estore[i] = ener[i];
    /* Set the mask to only do atom-atom interactions */
    set_nbfmask(fr,md,parm->ir.opts.ngener,FALSE);
    /* Compute remaining forces and energies now */
    do_force(log,cr,parm,nsb,vir_last,1,nrnb,
	     top,grps,pos[Try],v,force[Try],buf,md,ener,bVerbose && !PAR(cr),
	     lambda,graph,FALSE,TRUE,fr);
    /* Sum the potential energy terms from group contributions */
    sum_epot(&(parm->ir.opts),grps,ener);
    
    /* Sum over processors */
    if (PAR(cr)) 
      gmx_sum(NEPOT,Epot,cr);
      
    /* Add the stored energy contributions */
    for(i=0; (i<F_NRE); i++)
      ener[i] += Estore[i];
      
    /* Sum virial tensors */
    m_add(vir_last,my_vir[Min],vir_part);
    
    /* Last output line */
    if (bVerbose && MASTER(cr)) {
      Epot[Try]=ener[F_EPOT];
      print_epot(mdstep,count,step,Epot[Try],0.0,bOptim,TRUE);
    }
    /* Sum the forces from the previous calc. (shells only) 
     * and the current calc (atoms only)
     */
    for(i=0; (i<nsb->natoms); i++)
      rvec_add(force[Try][i],force[Min][i],f[i]);
  }
  else {
    /* Parallelise this one! */
    memcpy(f,force[Min],nsb->natoms*sizeof(f[0]));
    /* CHECK VIRIAL */
    copy_mat(my_vir[Min],vir_part);
  }
  memcpy(x,pos[Min],nsb->natoms*sizeof(x[0]));
  return count; 
}

