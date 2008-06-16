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
 * GROwing Monsters And Cloning Shrimps
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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
#include "gmx_fatal.h"
#include "txtdump.h"
#include "typedefs.h"
#include "update.h"
#include "random.h"
#include "constr.h"
#include "vec.h"
#include "statutil.h"
#include "tgroup.h"
#include "mdebin.h"
#include "vsite.h"
#include "force.h"
#include "mdrun.h"
#include "domdec.h"
#include "partdec.h"
#include "trnio.h"
#include "sparsematrix.h"
#include "mtxio.h"
#include "gmx_random.h"
#include "physics.h"
#include "xvgr.h"
#include "mdatoms.h"
#include "gbutil.h"
#include "ns.h"
#include "gmx_wallcycle.h"

typedef struct {
  t_state s;
  rvec    *f;
  real    epot;
  real    fnorm;
  real    fmax;
  int     a_fmax;
} em_state_t;

static em_state_t *init_em_state()
{
  em_state_t *ems;
  
  snew(ems,1);

  return ems;
}

static void sp_header(FILE *out,const char *minimizer,real ftol,int nsteps)
{
  fprintf(out,"%s:\n",minimizer);
  fprintf(out,"   Tolerance (Fmax)   = %12.5e\n",ftol);
  fprintf(out,"   Number of steps    = %12d\n",nsteps);
}

static void warn_step(FILE *fp,real ftol,bool bConstrain)
{
  fprintf(fp,"\nStepsize too small, or no change in energy.\n"
	  "Converged to machine precision,\n"
	  "but not to the requested precision Fmax < %g\n",
	  ftol);
  if (sizeof(real)<sizeof(double))
      fprintf(fp,"\nDouble precision normally gives you higher accuracy.\n");
	
  if (bConstrain)
    fprintf(fp,"You might need to increase your constraint accuracy, or turn\n"
	    "off constraints alltogether (set constraints = none in mdp file)\n");
}



static void print_converged(FILE *fp,const char *alg,real ftol,int count,bool bDone,
			    int nsteps,real epot,real fmax, int nfmax, real fnorm)
{
  if (bDone)
    fprintf(fp,"\n%s converged to Fmax < %g in %d steps\n",alg,ftol,count); 
  else if(count<nsteps)
    fprintf(fp,"\n%s converged to machine precision in %d steps,\n"
               "but did not reach the requested Fmax < %g.\n",alg,count,ftol);
  else 
    fprintf(fp,"\n%s did not converge to Fmax < %g in %d steps.\n",alg,ftol,count);

#ifdef GMX_DOUBLE
  fprintf(fp,"Potential Energy  = %21.14e\n",epot); 
  fprintf(fp,"Maximum force     = %21.14e on atom %d\n",fmax,nfmax+1); 
  fprintf(fp,"Norm of force     = %21.14e\n",fnorm); 
#else
  fprintf(fp,"Potential Energy  = %14.7e\n",epot); 
  fprintf(fp,"Maximum force     = %14.7e on atom %d\n",fmax,nfmax+1); 
  fprintf(fp,"Norm of force     = %14.7e\n",fnorm); 
#endif
}

static void get_f_norm_max(t_commrec *cr,
			   t_grpopts *opts,t_mdatoms *mdatoms,rvec *f,
			   real *fnorm,real *fmax,int *a_fmax)
{
  double fnorm2,*sum;
  real fmax2,fmax2_0,fam;
  int  la_max,a_max,start,end,i,m,gf;

  /* This routine finds the largest force and returns it.
   * On parallel machines the global max is taken.
   */
  fnorm2 = 0;
  fmax2 = 0;
  la_max = -1;
  gf = 0;
  start = mdatoms->start;
  end   = mdatoms->homenr + start;
  if (mdatoms->cFREEZE) {
    for(i=start; i<end; i++) {
      gf = mdatoms->cFREEZE[i];
      fam = 0;
      for(m=0; m<DIM; m++)
	if (!opts->nFreeze[gf][m])
	  fam += sqr(f[i][m]);
      fnorm2 += fam;
      if (fam > fmax2) {
	fmax2  = fam;
	la_max = i;
      }
    }
  } else {
    for(i=start; i<end; i++) {
      fam = norm2(f[i]);
      fnorm2 += fam;
      if (fam > fmax2) {
	fmax2  = fam;
	la_max = i;
      }
    }
  }

  if (la_max >= 0 && DOMAINDECOMP(cr)) {
    a_max = cr->dd->gatindex[la_max];
  } else {
    a_max = la_max;
  }
  if (PAR(cr)) {
    snew(sum,2*cr->nnodes+1);
    sum[2*cr->nodeid]   = fmax2;
    sum[2*cr->nodeid+1] = a_max;
    sum[2*cr->nnodes]   = fnorm2;
    gmx_sumd(2*cr->nnodes+1,sum,cr);
    fnorm2 = sum[2*cr->nnodes];
    /* Determine the global maximum */
    for(i=0; i<cr->nnodes; i++) {
      if (sum[2*i] > fmax2) {
	fmax2 = sum[2*i];
	a_max = (int)(sum[2*i+1] + 0.5);
      }
    }
    sfree(sum);
  }

  if (fnorm)
    *fnorm = sqrt(fnorm2);
  if (fmax)
    *fmax  = sqrt(fmax2);
  if (a_fmax)
    *a_fmax = a_max;
}

static void get_state_f_norm_max(t_commrec *cr,
			   t_grpopts *opts,t_mdatoms *mdatoms,
			   em_state_t *ems)
{
  get_f_norm_max(cr,opts,mdatoms,ems->f,&ems->fnorm,&ems->fmax,&ems->a_fmax);
}

void init_em(FILE *fplog,const char *title,
	     t_commrec *cr,t_inputrec *ir,t_groups *grps,
	     t_state *state_global,t_topology *top_global,
	     em_state_t *ems,t_topology **top,
	     rvec *f,rvec **buf,rvec **f_global,
	     t_nrnb *nrnb,rvec mu_tot,
	     t_forcerec *fr,t_graph *graph,t_mdatoms *mdatoms,
	     gmx_vsite_t *vsite,gmx_constr_t constr,
	     int nfile,t_filenm fnm[],int *fp_trn,int *fp_ene,
	     t_mdebin **mdebin)
{
  int  start,homenr,i;
  real dvdlambda;

  if (fplog)
    fprintf(fplog,"Initiating %s\n",title);

  state_global->ngtc = 0;
  if (ir->eI == eiCG) {
    state_global->flags |= (1<<estCGP);
    snew(state_global->cg_p,state_global->nalloc);
  }

  /* Initiate some variables */
  if (ir->efep != efepNO)
    state_global->lambda = ir->init_lambda;
  else 
    state_global->lambda = 0.0;

  init_nrnb(nrnb);
    
  if (DOMAINDECOMP(cr)) {
    if (vsite && DDMASTER(cr->dd)) {     
      construct_vsites(fplog,vsite,
		       state_global->x,nrnb,ir->delta_t,NULL,
		       &top_global->idef,ir->ePBC,ir->ePBC!=epbcNONE,graph,
		       DOMAINDECOMP(cr) ? NULL : cr,state_global->box);
    }

    *top = dd_init_local_top(top_global);

    dd_init_local_state(cr->dd,state_global,&ems->s);

    /* Distribute the charge groups over the nodes from the master node */
    dd_partition_system(fplog,ir->init_step,cr,TRUE,
			state_global,top_global,ir,
			&ems->s,&ems->f,buf,mdatoms,*top,
			fr,vsite,NULL,constr,
			nrnb,NULL,FALSE);
    dd_store_state(cr->dd,&ems->s);
    
    if (ir->nstfout)
      snew(*f_global,top_global->atoms.nr);
    else
      *f_global = NULL;
  } else {
    /* Just copy the state */
    ems->s = *state_global;
    snew(ems->s.x,ems->s.nalloc);
    snew(ems->f,ems->s.nalloc);
    for(i=0; i<state_global->natoms; i++)
      copy_rvec(state_global->x[i],ems->s.x[i]);
    copy_mat(state_global->box,ems->s.box);
    /* Copy the topology pointer */
    *top = top_global;
    *f_global = f;
  }

  clear_rvec(mu_tot);
  calc_shifts(ems->s.box,fr->shift_vec);

  if (PARTDECOMP(cr)) {
    pd_at_range(cr,&start,&homenr);
    homenr -= start;
  } else {
    start  = 0;
    homenr = (*top)->atoms.nr;
  }
  atoms2md(&(*top)->atoms,ir,(*top)->idef.il[F_ORIRES].nr,0,NULL,
	   start,homenr,mdatoms);
  update_mdatoms(mdatoms,state_global->lambda);

  if (constr) {
    if (ir->eConstrAlg == econtSHAKE && top_global->idef.il[F_CONSTR].nr > 0)
      gmx_fatal(FARGS,"Can not do energy minimization with %s, use %s\n",
		econstr_names[econtSHAKE],econstr_names[econtLINCS]);

    if (!DOMAINDECOMP(cr))
      set_constraints(constr,*top,ir,mdatoms,NULL);

    /* Constrain the starting coordinates */
    dvdlambda=0;
    constrain(PAR(cr) ? NULL : fplog,TRUE,TRUE,constr,*top,
	      ir,cr,-1,0,mdatoms,
	      ems->s.x,ems->s.x,NULL,ems->s.box,ems->s.lambda,&dvdlambda,
	      NULL,NULL,nrnb,econqCoord);
  }

  if (MASTER(cr)) {
    if (fp_trn)
      *fp_trn = open_trn(ftp2fn(efTRN,nfile,fnm),"w");
    if (fp_ene)
      *fp_ene = open_enx(ftp2fn(efENX,nfile,fnm),"w");
  } else {
    if (fp_trn)
      *fp_trn = -1;
    if (fp_ene)
      *fp_ene = -1;
  }

  /* Init bin for energy stuff */
  *mdebin = init_mdebin(*fp_ene,grps,
			&top_global->atoms,&top_global->idef,ir,cr); 
}

static void finish_em(FILE *fplog,t_commrec *cr,
		      int fp_traj,int fp_ene)
{
  if (MASTER(cr)) {
    close_trn(fp_traj);
    close_enx(fp_ene);
  }
}

static void swap_em_state(em_state_t *ems1,em_state_t *ems2)
{
  em_state_t tmp;

  tmp   = *ems1;
  *ems1 = *ems2;
  *ems2 = tmp;
}

static void copy_em_coords_back(em_state_t *ems,t_state *state)
{
  int i;

  for(i=0; (i<state->natoms); i++)
    copy_rvec(ems->s.x[i],state->x[i]);
}

static void do_em_step(t_commrec *cr,t_inputrec *ir,t_mdatoms *md,
		       em_state_t *ems1,real a,rvec *f,em_state_t *ems2,
		       gmx_constr_t constr,t_topology *top,t_nrnb *nrnb,
		       int count)

{
  t_state *s1,*s2;
  int  start,end,gf,i,m;
  rvec *x1,*x2;
  real dvdlambda;

  s1 = &ems1->s;
  s2 = &ems2->s;

  if (DOMAINDECOMP(cr) && s1->ddp_count != cr->dd->ddp_count)
    gmx_incons("state mismatch in do_em_step");

  s2->flags = s1->flags;

  if (s2->nalloc != s1->nalloc) {
    s2->nalloc = s1->nalloc;
    srenew(s2->x,s1->nalloc);
    srenew(ems2->f,  s1->nalloc);
    if (s2->flags & (1<<estCGP))
      srenew(s2->cg_p,  s1->nalloc);
  }
  
  s2->natoms = s1->natoms;
  s2->lambda = s1->lambda;
  copy_mat(s1->box,s2->box);

  start = md->start;
  end   = md->start + md->homenr;

  x1 = s1->x;
  x2 = s2->x;
  gf = 0;
  for(i=start; i<end; i++) {
    if (md->cFREEZE)
      gf = md->cFREEZE[i];
    for(m=0; m<DIM; m++) {
      if (ir->opts.nFreeze[gf][m])
	x2[i][m] = x1[i][m];
      else
	x2[i][m] = x1[i][m] + a*f[i][m];
    }
  }

  if (s2->flags & (1<<estCGP)) {
    /* Copy the CG p vector */
    x1 = s1->cg_p;
    x2 = s2->cg_p;
    for(i=start; i<end; i++)
      copy_rvec(x1[i],x2[i]);
  }

  if (DOMAINDECOMP(cr)) {
    s2->ddp_count = s1->ddp_count;
    if (s2->cg_gl_nalloc < s1->cg_gl_nalloc) {
      s2->cg_gl_nalloc = s1->cg_gl_nalloc;
      srenew(s2->cg_gl,s2->cg_gl_nalloc);
    }
    s2->ncg_gl = s1->ncg_gl;
    for(i=0; i<s2->ncg_gl; i++)
      s2->cg_gl[i] = s1->cg_gl[i];
    s2->ddp_count_cg_gl = s1->ddp_count_cg_gl;
  }

  if (constr) {
    dvdlambda = 0;
    constrain(NULL,TRUE,TRUE,constr,top,	
	      ir,cr,count,0,md,
	      s1->x,s2->x,NULL,s2->box,s2->lambda,
	      &dvdlambda,NULL,NULL,nrnb,econqCoord);
  }
}

static void do_x_step(t_commrec *cr,int n,rvec *x1,real a,rvec *f,rvec *x2)

{
  int  start,end,i,m;

  if (DOMAINDECOMP(cr)) {
    start = 0;
    end   = cr->dd->nat_home;
  } else if (PARTDECOMP(cr)) {
    pd_at_range(cr,&start,&end);
  } else {
    start = 0;
    end   = n;
  }

  for(i=start; i<end; i++) {
    for(m=0; m<DIM; m++) {
      x2[i][m] = x1[i][m] + a*f[i][m];
    }
  }
}

static void do_x_sub(t_commrec *cr,int n,rvec *x1,rvec *x2,real a,rvec *f)

{
  int  start,end,i,m;

  if (DOMAINDECOMP(cr)) {
    start = 0;
    end   = cr->dd->nat_home;
  } else if (PARTDECOMP(cr)) {
    pd_at_range(cr,&start,&end);
  } else {
    start = 0;
    end   = n;
  }

  for(i=start; i<end; i++) {
    for(m=0; m<DIM; m++) {
      f[i][m] = (x1[i][m] - x2[i][m])*a;
    }
  }
}

static void em_dd_partition_system(FILE *fplog,int step,t_commrec *cr,
				   t_topology *top_global,t_inputrec *ir,
				   em_state_t *ems,rvec **buf,t_topology *top,
				   t_mdatoms *mdatoms,t_forcerec *fr,
				   gmx_vsite_t *vsite,gmx_constr_t constr,
				   t_nrnb *nrnb,gmx_wallcycle_t wcycle)
{
  /* Repartition the domain decomposition */
  wallcycle_start(wcycle,ewcDOMDEC);
  dd_partition_system(fplog,step,cr,FALSE,
		      NULL,top_global,ir,
		      &ems->s,&ems->f,buf,
		      mdatoms,top,fr,vsite,NULL,constr,
		      nrnb,wcycle,FALSE);
  dd_store_state(cr->dd,&ems->s);
  wallcycle_stop(wcycle,ewcDOMDEC);
}
    
static void evaluate_energy(FILE *fplog,bool bVerbose,t_commrec *cr,
			    t_state *state_global,t_topology *top_global,
			    em_state_t *ems,rvec **buf,t_topology *top,
			    t_inputrec *inputrec,t_groups *grps,
			    t_nrnb *nrnb,gmx_wallcycle_t wcycle,
			    gmx_vsite_t *vsite,gmx_constr_t constr,
			    t_fcdata *fcd,
			    t_graph *graph,t_mdatoms *mdatoms,
			    t_forcerec *fr, rvec mu_tot,
			    real ener[],int count,bool bFirst)
{
  bool bNS;
  int  nabnsb;
  tensor force_vir,shake_vir,pres;
  real dvdl;
  real terminate=0;
  t_state *state;
  
  if (bFirst ||
      (DOMAINDECOMP(cr) && ems->s.ddp_count < cr->dd->ddp_count)) {
    /* This the first state or an old state used before the last ns */
    bNS = TRUE;
  } else {
    bNS = FALSE;
    if (inputrec->nstlist > 0) {
      bNS = TRUE;
    } else if (inputrec->nstlist == -1) {
      nabnsb = natoms_beyond_ns_buffer(inputrec,fr,&top->cgs,NULL,ems->s.x);
      if (PAR(cr))
	gmx_sumi(1,&nabnsb,cr);
      bNS = (nabnsb > 0);
    }
  }

  if (vsite)
    construct_vsites(fplog,vsite,ems->s.x,nrnb,1,NULL,&top->idef,
		     fr->ePBC,fr->bMolPBC,graph,cr,ems->s.box);

  if (DOMAINDECOMP(cr)) {
    if (bNS) {
      /* Repartition the domain decomposition */
      em_dd_partition_system(fplog,count,cr,top_global,inputrec,
			     ems,buf,top,mdatoms,fr,vsite,constr,
			     nrnb,wcycle);
    }
  }
      
  /* Calc force & energy on new trial position  */
  /* do_force always puts the charge groups in the box and shifts again
   * We do not unshift, so molecules are always whole in congrad.c
   */
  do_force(fplog,cr,inputrec,
	   count,nrnb,wcycle,top,grps,ems->s.box,ems->s.x,ems->f,
	   *buf,force_vir,mdatoms,ener,fcd,
	   ems->s.lambda,graph,fr,vsite,mu_tot,0.0,NULL,NULL,
	   GMX_FORCE_STATECHANGED | GMX_FORCE_ALLFORCES |
	   (bNS ? GMX_FORCE_NS : 0));

  /* Clear the unused shake virial and pressure */
  clear_mat(shake_vir);
  clear_mat(pres);

  /* Calculate long range corrections to pressure and energy */
  calc_dispcorr(fplog,inputrec,fr,count,mdatoms->nr,ems->s.box,ems->s.lambda,
		pres,force_vir,ener);

  /* Communicate stuff when parallel */
  if (PAR(cr)) 
    global_stat(fplog,cr,ener,force_vir,shake_vir,mu_tot,
		inputrec,grps,FALSE,NULL,NULL,NULL,NULL,&terminate);
    
  ener[F_ETOT] = ener[F_EPOT]; /* No kinetic energy */

  ems->epot = ener[F_EPOT];
  
  if (constr) {
    /* Project out the constraint components of the force */
    dvdl = 0;
    constrain(NULL,FALSE,FALSE,constr,top,
	      inputrec,cr,count,0,mdatoms,
	      ems->s.x,ems->f,ems->f,ems->s.box,ems->s.lambda,&dvdl,
	      NULL,NULL,nrnb,econqForce);
    if (fr->bSepDVDL && fplog)
      fprintf(fplog,sepdvdlformat,"Constraints",0.0,dvdl);
    ener[F_DVDL] += dvdl;
  }

  get_state_f_norm_max(cr,&(inputrec->opts),mdatoms,ems);
}

static double reorder_partsum(t_commrec *cr,t_grpopts *opts,t_mdatoms *mdatoms,
			      t_topology *top,
			      em_state_t *s_min,em_state_t *s_b)
{
  rvec *fm,*fb,*fmg;
  t_atom *atom;
  int ncg,*cg_gl,*index,c,cg,i,a0,a1,a,gf,m;
  double partsum;

  if (debug)
    fprintf(debug,"Doing reorder_partsum\n");

  fm = s_min->f;
  fb = s_b->f;
  atom = top->atoms.atom;

  index = top->cgs.index;

  /* Collect fm in a global vector fmg.
   * This conflics with the spirit of domain decomposition,
   * but to fully optimize this a much more complicated algorithm is required.
   */
  snew(fmg,top->atoms.nr);
  
  ncg   = s_min->s.ncg_gl;
  cg_gl = s_min->s.cg_gl;
  i = 0;
  for(c=0; c<ncg; c++) {
    cg = cg_gl[c];
    a0 = index[cg];
    a1 = index[cg+1];
    for(a=a0; a<a1; a++) {
      copy_rvec(fm[i],fmg[a]);
      i++;
    }
  }
  gmx_sum(top->atoms.nr*3,fmg[0],cr);

  /* Now we will determine the part of the sum for the cgs in state s_b */
  ncg   = s_b->s.ncg_gl;
  cg_gl = s_b->s.cg_gl;
  partsum = 0;
  i = 0;
  gf = 0;
  for(c=0; c<ncg; c++) {
    cg = cg_gl[c];
    a0 = index[cg];
    a1 = index[cg+1];
    for(a=a0; a<a1; a++) {
      if (mdatoms->cFREEZE)
	gf = atom[i].grpnr[egcFREEZE];
      for(m=0; m<DIM; m++) {
	if (!opts->nFreeze[gf][m]) {
	  partsum += (fb[i][m] - fmg[a][m])*fb[i][m];
	}
      }
      i++;
    }
  }
  
  sfree(fmg);

  return partsum;
}

static real pr_beta(t_commrec *cr,t_grpopts *opts,t_mdatoms *mdatoms,
		    t_topology *top,
		    em_state_t *s_min,em_state_t *s_b)
{
  rvec *fm,*fb;
  double sum;
  int  gf,i,m;

  /* This is just the classical Polak-Ribiere calculation of beta;
   * it looks a bit complicated since we take freeze groups into account,
   * and might have to sum it in parallel runs.
   */
  
  if (!DOMAINDECOMP(cr) ||
      (s_min->s.ddp_count == cr->dd->ddp_count &&
       s_b->s.ddp_count   == cr->dd->ddp_count)) {
    fm = s_min->f;
    fb = s_b->f;
    sum = 0;
    gf = 0;
    /* This part of code can be incorrect with DD,
     * since the atom ordering in s_b and s_min might differ.
     */
    for(i=mdatoms->start; i<mdatoms->start+mdatoms->homenr; i++) {
      if (mdatoms->cFREEZE)
	gf = mdatoms->cFREEZE[i];
      for(m=0; m<DIM; m++)
	if (!opts->nFreeze[gf][m]) {
	  sum += (fb[i][m] - fm[i][m])*fb[i][m];
	} 
    }
  } else {
    /* We need to reorder cgs while summing */
    sum = reorder_partsum(cr,opts,mdatoms,top,s_min,s_b);
  }
  if (PAR(cr))
    gmx_sumd(1,&sum,cr);

  return sum/sqr(s_min->fnorm);
}

time_t do_cg(FILE *fplog,t_commrec *cr,
	     int nfile,t_filenm fnm[],
	     bool bVerbose,bool bCompact,
	     gmx_vsite_t *vsite,gmx_constr_t constr,
	     int stepout,
	     t_inputrec *inputrec,t_groups *grps,
	     t_topology *top_global,real ener[],t_fcdata *fcd,
	     t_state *state_global,rvec f[],
	     rvec buf[],t_mdatoms *mdatoms,
	     t_nrnb *nrnb,gmx_wallcycle_t wcycle,
	     t_graph *graph,t_edsamyn *edyn,
	     t_forcerec *fr,
	     int repl_ex_nst,int repl_ex_seed,
	     real cpt_period,real max_hours,
	     unsigned long Flags)
{
  const char *CG="Polak-Ribiere Conjugate Gradients";

  em_state_t *s_min,*s_a,*s_b,*s_c;
  t_topology *top;
  rvec   *f_global,*p,*sf,*sfm;
  double gpa,gpb,gpc,tmp,sum[2],minstep;
  real   fnormn;
  real   stepsize;	
  real   a,b,c,beta=0.0;
  real   epot_repl=0;
  real   pnorm;
  t_mdebin   *mdebin;
  bool   converged,foundlower;
  rvec   mu_tot;
  time_t start_t;
  bool   do_log=FALSE,do_ene=FALSE,do_x,do_f;
  tensor force_vir,shake_vir,vir,pres;
  int    number_steps,neval=0,nstcg=inputrec->nstcgsteep;
  int    fp_trn,fp_ene;
  int    i,m,gf,step,nminstep;
  real   terminate=0;  

  step=0;

  s_min = init_em_state();
  s_a   = init_em_state();
  s_b   = init_em_state();
  s_c   = init_em_state();

  /* Init em and store the local state in s_min */
  init_em(fplog,CG,cr,inputrec,grps,
	  state_global,top_global,s_min,&top,f,&buf,&f_global,
	  nrnb,mu_tot,fr,graph,mdatoms,vsite,constr,
	  nfile,fnm,&fp_trn,&fp_ene,&mdebin);

  /* Print to log file */
  start_t=print_date_and_time(fplog,cr->nodeid,
			      "Started Polak-Ribiere Conjugate Gradients");
  wallcycle_start(wcycle,ewcRUN);
  
  /* Clear some matrix variables */
  clear_mat(shake_vir);
  clear_mat(vir);
  clear_mat(pres);
  
  /* Max number of steps */
  number_steps=inputrec->nsteps;

  if (MASTER(cr))
    sp_header(stderr,CG,inputrec->em_tol,number_steps);
  if (fplog)
    sp_header(fplog,CG,inputrec->em_tol,number_steps);

  /* Call the force routine and some auxiliary (neighboursearching etc.) */
  /* do_force always puts the charge groups in the box and shifts again
   * We do not unshift, so molecules are always whole in congrad.c
   */
  evaluate_energy(fplog,bVerbose,cr,
		  state_global,top_global,s_min,&buf,top,
		  inputrec,grps,nrnb,wcycle,
		  vsite,constr,fcd,graph,mdatoms,fr,mu_tot,ener,-1,TRUE);
  where();

  if (MASTER(cr)) {
    /* Copy stuff to the energy bin for easy printing etc. */
    upd_mdebin(mdebin,NULL,TRUE,mdatoms->tmass,step,(real)step,
	       ener,&s_min->s,s_min->s.box,shake_vir,
	       force_vir,vir,pres,grps,mu_tot,constr);
    
    print_ebin_header(fplog,step,step,s_min->s.lambda);
    print_ebin(fp_ene,TRUE,FALSE,FALSE,fplog,step,step,step,eprNORMAL,
	       TRUE,mdebin,fcd,&(top->atoms),&(inputrec->opts));
  }
  where();

  /* Estimate/guess the initial stepsize */
  stepsize = inputrec->em_stepsize/s_min->fnorm;
 
  if (MASTER(cr)) {
    fprintf(stderr,"   F-max             = %12.5e on atom %d\n",
	    s_min->fmax,s_min->a_fmax+1);
    fprintf(stderr,"   F-Norm            = %12.5e\n",
	    s_min->fnorm/sqrt(state_global->natoms));
    fprintf(stderr,"\n");
    /* and copy to the log file too... */
    fprintf(fplog,"   F-max             = %12.5e on atom %d\n",
	    s_min->fmax,s_min->a_fmax+1);
    fprintf(fplog,"   F-Norm            = %12.5e\n",
	    s_min->fnorm/sqrt(state_global->natoms));
    fprintf(fplog,"\n");
  }  
  /* Start the loop over CG steps.		
   * Each successful step is counted, and we continue until
   * we either converge or reach the max number of steps.
   */
  for(step=0,converged=FALSE;( step<=number_steps || number_steps==0) && !converged;step++) {
    
    /* start taking steps in a new direction 
     * First time we enter the routine, beta=0, and the direction is 
     * simply the negative gradient.
     */

    /* Calculate the new direction in p, and the gradient in this direction, gpa */
    p  = s_min->s.cg_p;
    sf = s_min->f;
    gpa = 0;
    gf = 0;
    for(i=mdatoms->start; i<mdatoms->start+mdatoms->homenr; i++) {
      if (mdatoms->cFREEZE) 
	gf = mdatoms->cFREEZE[i];
      for(m=0; m<DIM; m++) {
	if (!inputrec->opts.nFreeze[gf][m]) {
	  p[i][m] = sf[i][m] + beta*p[i][m];
	  gpa -= p[i][m]*sf[i][m];
	  /* f is negative gradient, thus the sign */
	} else {
          p[i][m] = 0;
	}
      }
    }
    
    /* Sum the gradient along the line across CPUs */
    if (PAR(cr))
      gmx_sumd(1,&gpa,cr);

    /* Calculate the norm of the search vector */
    get_f_norm_max(cr,&(inputrec->opts),mdatoms,p,&pnorm,NULL,NULL);
    
    /* Just in case stepsize reaches zero due to numerical precision... */
    if(stepsize<=0)	  
      stepsize = inputrec->em_stepsize/pnorm;
    
    /* 
     * Double check the value of the derivative in the search direction.
     * If it is positive it must be due to the old information in the
     * CG formula, so just remove that and start over with beta=0.
     * This corresponds to a steepest descent step.
     */
    if(gpa>0) {
      beta = 0;
      step--; /* Don't count this step since we are restarting */
      continue; /* Go back to the beginning of the big for-loop */
    }

    /* Calculate minimum allowed stepsize, before the average (norm)
     * relative change in coordinate is smaller than precision
     */
    minstep=0;
    for (i=mdatoms->start; i<mdatoms->start+mdatoms->homenr; i++) {
      for(m=0; m<DIM; m++) {
	tmp = fabs(s_min->s.x[i][m]);
	if(tmp < 1.0)
	  tmp = 1.0;
	tmp = p[i][m]/tmp;
	minstep += tmp*tmp;
      }
    }
    /* Add up from all CPUs */
    if(PAR(cr))
      gmx_sumd(1,&minstep,cr);

    minstep = GMX_REAL_EPS/sqrt(minstep/(3*state_global->natoms));

    if(stepsize<minstep) {
      converged=TRUE;
      break;
    }
    
    /* Write coordinates if necessary */
    do_x = do_per_step(step,inputrec->nstxout);
    do_f = do_per_step(step,inputrec->nstfout);
    
    write_traj(fplog,cr,fp_trn,do_x,FALSE,do_f,0,FALSE,0,NULL,FALSE,
	       top_global,inputrec->eI,step,(double)step,
	       &s_min->s,state_global,s_min->f,f_global);
    
    /* Take a step downhill.
     * In theory, we should minimize the function along this direction.
     * That is quite possible, but it turns out to take 5-10 function evaluations
     * for each line. However, we dont really need to find the exact minimum -
     * it is much better to start a new CG step in a modified direction as soon
     * as we are close to it. This will save a lot of energy evaluations.
     *
     * In practice, we just try to take a single step.
     * If it worked (i.e. lowered the energy), we increase the stepsize but
     * the continue straight to the next CG step without trying to find any minimum.
     * If it didn't work (higher energy), there must be a minimum somewhere between
     * the old position and the new one.
     * 
     * Due to the finite numerical accuracy, it turns out that it is a good idea
     * to even accept a SMALL increase in energy, if the derivative is still downhill.
     * This leads to lower final energies in the tests I've done. / Erik 
     */
    s_a->epot = s_min->epot;
    a = 0.0;
    c = a + stepsize; /* reference position along line is zero */
    
    /* Take a trial step (new coords in s_c) */
    do_em_step(cr,inputrec,mdatoms,s_min,c,s_min->s.cg_p,s_c,
	       constr,top,nrnb,-1);
    
    neval++;
    /* Calculate energy for the trial step */
    evaluate_energy(fplog,bVerbose,cr,
		    state_global,top_global,s_c,&buf,top,
		    inputrec,grps,nrnb,wcycle,
		    vsite,constr,fcd,graph,mdatoms,fr,mu_tot,ener,-1,FALSE);
    
    /* Calc derivative along line */
    p  = s_c->s.cg_p;
    sf = s_c->f;
    gpc=0;
    for(i=mdatoms->start; i<mdatoms->start+mdatoms->homenr; i++) {
      for(m=0; m<DIM; m++) 
	  gpc -= p[i][m]*sf[i][m];  /* f is negative gradient, thus the sign */
    }
    /* Sum the gradient along the line across CPUs */
    if (PAR(cr))
      gmx_sumd(1,&gpc,cr);

    /* This is the max amount of increase in energy we tolerate */
    tmp=sqrt(GMX_REAL_EPS)*fabs(s_a->epot);

    /* Accept the step if the energy is lower, or if it is not significantly higher
     * and the line derivative is still negative.
     */
    if (s_c->epot < s_a->epot || (gpc < 0 && s_c->epot < (s_a->epot + tmp))) {
      foundlower = TRUE;
      /* Great, we found a better energy. Increase step for next iteration
       * if we are still going down, decrease it otherwise
       */
      if(gpc<0)
	stepsize *= 1.618034;  /* The golden section */
      else
	stepsize *= 0.618034;  /* 1/golden section */
    } else {
      /* New energy is the same or higher. We will have to do some work
       * to find a smaller value in the interval. Take smaller step next time!
       */
      foundlower = FALSE;
      stepsize *= 0.618034;
    }    



    
    /* OK, if we didn't find a lower value we will have to locate one now - there must
     * be one in the interval [a=0,c].
     * The same thing is valid here, though: Don't spend dozens of iterations to find
     * the line minimum. We try to interpolate based on the derivative at the endpoints,
     * and only continue until we find a lower value. In most cases this means 1-2 iterations.
     *
     * I also have a safeguard for potentially really patological functions so we never
     * take more than 20 steps before we give up ...
     *
     * If we already found a lower value we just skip this step and continue to the update.
     */
    if (!foundlower) {
      nminstep=0;

      do {
	/* Select a new trial point.
	 * If the derivatives at points a & c have different sign we interpolate to zero,
	 * otherwise just do a bisection.
	 */
	if(gpa<0 && gpc>0)
	  b = a + gpa*(a-c)/(gpc-gpa);
	else
	  b = 0.5*(a+c);		
	
	/* safeguard if interpolation close to machine accuracy causes errors:
	 * never go outside the interval
	 */
	if(b<=a || b>=c)
	  b = 0.5*(a+c);
	
	if (DOMAINDECOMP(cr) && s_min->s.ddp_count != cr->dd->ddp_count) {
	  /* Reload the old state */
	  em_dd_partition_system(fplog,-1,cr,top_global,inputrec,
				 s_min,&buf,top,mdatoms,fr,vsite,constr,
				 nrnb,wcycle);
	}

	/* Take a trial step to this new point - new coords in s_b */
	do_em_step(cr,inputrec,mdatoms,s_min,b,s_min->s.cg_p,s_b,
		   constr,top,nrnb,-1);
	
	neval++;
	/* Calculate energy for the trial step */
	evaluate_energy(fplog,bVerbose,cr,
			state_global,top_global,s_b,&buf,top,
			inputrec,grps,nrnb,wcycle,
			vsite,constr,fcd,graph,mdatoms,fr,mu_tot,ener,
			-1,FALSE);
	
	/* p does not change within a step, but since the domain decomposition
	 * might change, we have to use cg_p of s_b here.
	 */
	p  = s_b->s.cg_p;
	sf = s_b->f;
	gpb=0;
	for(i=mdatoms->start; i<mdatoms->start+mdatoms->homenr; i++) {
	  for(m=0; m<DIM; m++)
	      gpb -= p[i][m]*sf[i][m];   /* f is negative gradient, thus the sign */
	}
	/* Sum the gradient along the line across CPUs */
	if (PAR(cr))
	  gmx_sumd(1,&gpb,cr);
	
	if (debug)
	  fprintf(debug,"CGE: EpotA %f EpotB %f EpotC %f gpb %f\n",
		  s_a->epot,s_b->epot,s_c->epot,gpb);

	epot_repl = s_b->epot;
	
	/* Keep one of the intervals based on the value of the derivative at the new point */
	if (gpb > 0) {
	  /* Replace c endpoint with b */
	  swap_em_state(s_b,s_c);
	  c = b;
	  gpc = gpb;
	} else {
	  /* Replace a endpoint with b */
	  swap_em_state(s_b,s_a);
	  a = b;
	  gpa = gpb;
	}
	
	/* 
	 * Stop search as soon as we find a value smaller than the endpoints.
	 * Never run more than 20 steps, no matter what.
	 */
	nminstep++;
      } while ((epot_repl > s_a->epot || epot_repl > s_c->epot) &&
	       (nminstep < 20));     
      
      if (fabs(epot_repl - s_min->epot) < fabs(s_min->epot)*GMX_REAL_EPS ||
	  nminstep >= 20) {
	/* OK. We couldn't find a significantly lower energy.
	 * If beta==0 this was steepest descent, and then we give up.
	 * If not, set beta=0 and restart with steepest descent before quitting.
         */
	if (beta == 0.0) {
	  /* Converged */
	  converged = TRUE;
	  break;
	} else {
	  /* Reset memory before giving up */
	  beta = 0.0;
	  continue;
	}
      }
      
      /* Select min energy state of A & C, put the best in B.
       */
      if (s_c->epot < s_a->epot) {
	if (debug)
	  fprintf(debug,"CGE: C (%f) is lower than A (%f), moving C to B\n",
		  s_c->epot,s_a->epot);
	swap_em_state(s_b,s_c);
	gpb = gpc;
	b = c;
      } else {
	if (debug)
	  fprintf(debug,"CGE: A (%f) is lower than C (%f), moving A to B\n",
		  s_a->epot,s_c->epot);
	swap_em_state(s_b,s_a);
	gpb = gpa;
	b = a;
      }
      
    } else {
      if (debug)
	fprintf(debug,"CGE: Found a lower energy %f, moving C to B\n",
		s_c->epot);
      swap_em_state(s_b,s_c);
      gpb = gpc;
      b = c;
    }
    
    /* new search direction */
    /* beta = 0 means forget all memory and restart with steepest descents. */
    if (nstcg && ((step % nstcg)==0)) 
      beta = 0.0;
    else {
      /* s_min->fnorm cannot be zero, because then we would have converged
       * and broken out.
       */

      /* Polak-Ribiere update.
       * Change to fnorm2/fnorm2_old for Fletcher-Reeves
       */
      beta = pr_beta(cr,&inputrec->opts,mdatoms,top_global,s_min,s_b);
    }
    /* Limit beta to prevent oscillations */
    if (fabs(beta) > 5.0)
      beta = 0.0;
    
    
    /* update positions */
    swap_em_state(s_min,s_b);
    gpa = gpb;
    
    /* Print it if necessary */
    if (MASTER(cr)) {
      ener[F_ETOT] = ener[F_EPOT]; /* No kinetic energy */
      if(bVerbose)
	fprintf(stderr,"\rStep %d, Epot=%12.6e, Fnorm=%9.3e, Fmax=%9.3e (atom %d)\n",
		step,s_min->epot,s_min->fnorm/sqrt(state_global->natoms),
		s_min->fmax,s_min->a_fmax+1);
      /* Store the new (lower) energies */
      upd_mdebin(mdebin,NULL,TRUE,mdatoms->tmass,step,(real)step,
		 ener,&s_min->s,s_min->s.box,shake_vir,
		 force_vir,vir,pres,grps,mu_tot,constr);
      do_log = do_per_step(step,inputrec->nstlog);
      do_ene = do_per_step(step,inputrec->nstenergy);
      if(do_log)
	print_ebin_header(fplog,step,step,s_min->s.lambda);
      print_ebin(fp_ene,do_ene,FALSE,FALSE,
		 do_log ? fplog : NULL,step,step,step,eprNORMAL,
		 TRUE,mdebin,fcd,&(top->atoms),&(inputrec->opts));
    }
    
    /* Stop when the maximum force lies below tolerance.
     * If we have reached machine precision, converged is already set to true.
     */	
    converged = converged || (s_min->fmax < inputrec->em_tol);
    
  } /* End of the loop */
  
  if (converged)	
    step--; /* we never took that last step in this case */
  
  if (s_min->fmax > inputrec->em_tol) {
    if (MASTER(cr)) {
      warn_step(stderr,inputrec->em_tol,FALSE);
      warn_step(fplog,inputrec->em_tol,FALSE);
    }
    converged = FALSE; 
  }
  
  if (MASTER(cr)) {
    /* If we printed energy and/or logfile last step (which was the last step)
     * we don't have to do it again, but otherwise print the final values.
     */
    if(!do_log) {
      /* Write final value to log since we didn't do anything the last step */
      print_ebin_header(fplog,step,step,s_min->s.lambda);
    }
    if (!do_ene || !do_log) {
      /* Write final energy file entries */
      print_ebin(fp_ene,!do_ene,FALSE,FALSE,
		 !do_log ? fplog : NULL,step,step,step,eprNORMAL,
		 TRUE,mdebin,fcd,&(top->atoms),&(inputrec->opts));
    }
  }
  
  if (!PAR(cr)) {
    copy_em_coords_back(s_min,state_global);
  }

  /* Print some stuff... */
  if (MASTER(cr))
    fprintf(stderr,"\nwriting lowest energy coordinates.\n");
  
	/* IMPORTANT!
	 * For accurate normal mode calculation it is imperative that we
	 * store the last conformation into the full precision binary trajectory
	 .
	 *
	 * However, we should only do it if we did NOT already write this step
	 * above (which we did if do_x or do_f was true).
	 */  
	do_x = !do_per_step(step,inputrec->nstxout);
    do_f = !do_per_step(step,inputrec->nstfout);

	write_traj(fplog,cr,fp_trn,do_x,FALSE,do_f,
	     0,FALSE,0,NULL,FALSE,
	     top_global,inputrec->eI,step,(real)step,
	     &s_min->s,state_global,s_min->f,f_global);
  if (MASTER(cr))
    write_sto_conf(ftp2fn(efSTO,nfile,fnm),
		   *top_global->name, &(top_global->atoms),
		   state_global->x,NULL,inputrec->ePBC,state_global->box);
  
  fnormn = s_min->fnorm/sqrt(state_global->natoms);
  
  if (MASTER(cr)) {
    print_converged(stderr,CG,inputrec->em_tol,step,converged,number_steps,
		    s_min->epot,s_min->fmax,s_min->a_fmax,fnormn);
    print_converged(fplog,CG,inputrec->em_tol,step,converged,number_steps,
		    s_min->epot,s_min->fmax,s_min->a_fmax,fnormn);
    
    fprintf(fplog,"\nPerformed %d energy evaluations in total.\n",neval);
  }
  
  finish_em(fplog,cr,fp_trn,fp_ene);
  
  /* To print the actual number of steps we needed somewhere */
  inputrec->nsteps=step;

  return start_t;
} /* That's all folks */


time_t do_lbfgs(FILE *fplog,t_commrec *cr,
		int nfile,t_filenm fnm[],
		bool bVerbose,bool bCompact,
		gmx_vsite_t *vsite,gmx_constr_t constr,
		int stepout,
		t_inputrec *inputrec,t_groups *grps,
		t_topology *top,real ener[],t_fcdata *fcd,
		t_state *state,rvec f[],
		rvec buf[],t_mdatoms *mdatoms,
		t_nrnb *nrnb,gmx_wallcycle_t wcycle,
		t_graph *graph,t_edsamyn *edyn,
		t_forcerec *fr,
		int repl_ex_nst,int repl_ex_seed,
		real cpt_period,real max_hours,
		unsigned long Flags)
{
  static char *LBFGS="Low-Memory BFGS Minimizer";
  em_state_t ems;
  int    ncorr,nmaxcorr,point,cp,neval,nminstep;
  double stepsize,gpa,gpb,gpc,tmp,minstep;
  real   *rho,*alpha,*ff,*xx,*p,*s,*lastx,*lastf,**dx,**dg;	
  real   *xa,*xb,*xc,*fa,*fb,*fc,*xtmp,*ftmp;
  real   a,b,c,maxdelta,delta;
  real   diag,Epot0,Epot,EpotA,EpotB,EpotC;
  real   dgdx,dgdg,sq,yr,beta;
  t_mdebin   *mdebin;
  bool   converged,first;
  rvec   mu_tot;
  real   fnorm,fmax;
  time_t start_t;
  bool   do_log,do_ene,do_x,do_f,foundlower,*frozen;
  tensor force_vir,shake_vir,vir,pres;
  int    fp_trn,fp_ene,start,end,number_steps;
  int    i,k,m,n,nfmax,gf,step;
  /* not used */
  real   terminate;
  
  if (PAR(cr))
    gmx_fatal(FARGS,"Cannot do parallel L-BFGS Minimization - yet.\n");
  
  n = 3*state->natoms;
  nmaxcorr = inputrec->nbfgscorr;
  
  /* Allocate memory */
  snew(f,state->natoms);
  /* Use pointers to real so we dont have to loop over both atoms and
   * dimensions all the time...
   * x/f are allocated as rvec *, so make new x0/f0 pointers-to-real
   * that point to the same memory.
   */
  snew(xa,n);
  snew(xb,n);
  snew(xc,n);
  snew(fa,n);
  snew(fb,n);
  snew(fc,n);
  snew(frozen,n);
  
  xx = (real *)state->x;
  ff = (real *)f;

	printf("x0: %20.12g %20.12g %20.12g\n",xx[0],xx[1],xx[2]);

  snew(p,n); 
  snew(lastx,n); 
  snew(lastf,n); 
  snew(rho,nmaxcorr);
  snew(alpha,nmaxcorr);
  
  snew(dx,nmaxcorr);
  for(i=0;i<nmaxcorr;i++)
    snew(dx[i],n);
  
  snew(dg,nmaxcorr);
  for(i=0;i<nmaxcorr;i++)
    snew(dg[i],n);

  step = 0;
  neval = 0; 

  /* Init em */
  init_em(fplog,LBFGS,cr,inputrec,grps,
	  state,top,&ems,&top,f,&buf,&f,
	  nrnb,mu_tot,fr,graph,mdatoms,vsite,constr,
	  nfile,fnm,&fp_trn,&fp_ene,&mdebin);
  /* Do_lbfgs is not completely updated like do_steep and do_cg,
   * so we free some memory again.
   */
  sfree(ems.s.x);
  sfree(ems.f);

  start = mdatoms->start;
  end   = mdatoms->homenr + start;
    
  /* Print to log file */
  start_t=print_date_and_time(fplog,cr->nodeid,
			      "Started Low-Memory BFGS Minimization");
  wallcycle_start(wcycle,ewcRUN);
  
  do_log = do_ene = do_x = do_f = TRUE;
  
  /* Clear some matrix variables */
  clear_mat(force_vir);
  clear_mat(shake_vir);
  clear_mat(vir);
  clear_mat(pres);
  
  /* Max number of steps */
  number_steps=inputrec->nsteps;

  /* Create a 3*natoms index to tell whether each degree of freedom is frozen */
  gf = 0;
  for(i=start; i<end; i++) {
    if (mdatoms->cFREEZE)
      gf = mdatoms->cFREEZE[i];
     for(m=0; m<DIM; m++) 
       frozen[3*i+m]=inputrec->opts.nFreeze[gf][m];  
  }
  if (MASTER(cr))
    sp_header(stderr,LBFGS,inputrec->em_tol,number_steps);
  if (fplog)
    sp_header(fplog,LBFGS,inputrec->em_tol,number_steps);
  
  if (vsite)
    construct_vsites(fplog,vsite,state->x,nrnb,1,NULL,&top->idef,
		     fr->ePBC,fr->bMolPBC,graph,cr,state->box);
  
  /* Call the force routine and some auxiliary (neighboursearching etc.) */
  /* do_force always puts the charge groups in the box and shifts again
   * We do not unshift, so molecules are always whole
   */
  neval++;
  ems.s.x = state->x;
  ems.f = f;
  evaluate_energy(fplog,bVerbose,cr,
		  state,top,&ems,&buf,top,
		  inputrec,grps,nrnb,wcycle,
		  vsite,constr,fcd,graph,mdatoms,fr,mu_tot,ener,-1,TRUE);
  where();
	
  if (MASTER(cr)) {
    /* Copy stuff to the energy bin for easy printing etc. */
    upd_mdebin(mdebin,NULL,TRUE,mdatoms->tmass,step,(real)step,
	       ener,state,state->box,shake_vir,
	       force_vir,vir,pres,grps,mu_tot,constr);
    
    print_ebin_header(fplog,step,step,state->lambda);
    print_ebin(fp_ene,TRUE,FALSE,FALSE,fplog,step,step,step,eprNORMAL,
	       TRUE,mdebin,fcd,&(top->atoms),&(inputrec->opts));
  }
  where();
  
  /* This is the starting energy */
  Epot=ener[F_EPOT];
  
  fnorm = ems.fnorm;
  fmax  = ems.fmax;
  nfmax = ems.a_fmax;
  
  /* Set the initial step.
   * since it will be multiplied by the non-normalized search direction 
   * vector (force vector the first time), we scale it by the
   * norm of the force.
   */
  
  if (MASTER(cr)) {
    fprintf(stderr,"Using %d BFGS correction steps.\n\n",nmaxcorr);
    fprintf(stderr,"   F-max             = %12.5e on atom %d\n",fmax,nfmax+1);
    fprintf(stderr,"   F-Norm            = %12.5e\n",fnorm/sqrt(state->natoms));
    fprintf(stderr,"\n");
    /* and copy to the log file too... */
    fprintf(fplog,"Using %d BFGS correction steps.\n\n",nmaxcorr);
    fprintf(fplog,"   F-max             = %12.5e on atom %d\n",fmax,nfmax+1);
    fprintf(fplog,"   F-Norm            = %12.5e\n",fnorm/sqrt(state->natoms));
    fprintf(fplog,"\n");
  }   
  
  point=0;
  for(i=0;i<n;i++)
    if(!frozen[i])
      dx[point][i] = ff[i];  /* Initial search direction */
    else
      dx[point][i] = 0;

  stepsize = 1.0/fnorm;
  converged = FALSE;
  
  /* Start the loop over BFGS steps.		
   * Each successful step is counted, and we continue until
   * we either converge or reach the max number of steps.
   */
  
  ncorr=0;

  /* Set the gradient from the force */
  for(step=0,converged=FALSE;( step<=number_steps || number_steps==0) && !converged;step++) {
    
    /* Write coordinates if necessary */
    do_x = do_per_step(step,inputrec->nstxout);
    do_f = do_per_step(step,inputrec->nstfout);
    
    write_traj(fplog,cr,fp_trn,do_x,FALSE,do_f,0,FALSE,0,NULL,FALSE,
	       top,inputrec->eI,step,(real)step,state,state,f,f);

    /* Do the linesearching in the direction dx[point][0..(n-1)] */
    
    /* pointer to current direction - point=0 first time here */
    s=dx[point];
    
    /* calculate line gradient */
    for(gpa=0,i=0;i<n;i++) 
	gpa-=s[i]*ff[i];

    /* Calculate minimum allowed stepsize, before the average (norm) 
     * relative change in coordinate is smaller than precision 
     */
    for(minstep=0,i=0;i<n;i++) {
      tmp=fabs(xx[i]);
      if(tmp<1.0)
	tmp=1.0;
      tmp = s[i]/tmp;
      minstep += tmp*tmp;
    }
    minstep = GMX_REAL_EPS/sqrt(minstep/n);
    
    if(stepsize<minstep) {
      converged=TRUE;
      break;
    }
    
    /* Store old forces and coordinates */
    for(i=0;i<n;i++) {
      lastx[i]=xx[i];
      lastf[i]=ff[i];
    }
    Epot0=Epot;
    
    first=TRUE;
    
    for(i=0;i<n;i++)
      xa[i]=xx[i];
    
    /* Take a step downhill.
     * In theory, we should minimize the function along this direction.
     * That is quite possible, but it turns out to take 5-10 function evaluations
     * for each line. However, we dont really need to find the exact minimum -
     * it is much better to start a new BFGS step in a modified direction as soon
     * as we are close to it. This will save a lot of energy evaluations.
     *
     * In practice, we just try to take a single step.
     * If it worked (i.e. lowered the energy), we increase the stepsize but
     * the continue straight to the next BFGS step without trying to find any minimum.
     * If it didn't work (higher energy), there must be a minimum somewhere between
     * the old position and the new one.
     * 
     * Due to the finite numerical accuracy, it turns out that it is a good idea
     * to even accept a SMALL increase in energy, if the derivative is still downhill.
     * This leads to lower final energies in the tests I've done. / Erik 
     */
    foundlower=FALSE;
    EpotA = Epot0;
    a = 0.0;
    c = a + stepsize; /* reference position along line is zero */

    /* Check stepsize first. We do not allow displacements 
     * larger than emstep.
     */
    do {
      c = a + stepsize;
      maxdelta=0;
      for(i=0;i<n;i++) {
	delta=c*s[i];
	if(delta>maxdelta)
	  maxdelta=delta;
      }
      if(maxdelta>inputrec->em_stepsize)
	stepsize*=0.1;
    } while(maxdelta>inputrec->em_stepsize);

    /* Take a trial step */
    for (i=0; i<n; i++)
      xc[i] = lastx[i] + c*s[i];
    
    neval++;
    /* Calculate energy for the trial step */
    ems.s.x = (rvec *)xc;
    ems.f   = (rvec *)fc;
    evaluate_energy(fplog,bVerbose,cr,
		    state,top,&ems,&buf,top,
		    inputrec,grps,nrnb,wcycle,
		    vsite,constr,fcd,graph,mdatoms,fr,mu_tot,ener,step,FALSE);
    EpotC = ems.epot;
    
    /* Calc derivative along line */
    for(gpc=0,i=0; i<n; i++) {
	gpc -= s[i]*fc[i];   /* f is negative gradient, thus the sign */
    }
    /* Sum the gradient along the line across CPUs */
    if (PAR(cr))
      gmx_sumd(1,&gpc,cr);
    
     /* This is the max amount of increase in energy we tolerate */
   tmp=sqrt(GMX_REAL_EPS)*fabs(EpotA);
    
    /* Accept the step if the energy is lower, or if it is not significantly higher
     * and the line derivative is still negative.
     */
    if(EpotC<EpotA || (gpc<0 && EpotC<(EpotA+tmp))) {
      foundlower = TRUE;
      /* Great, we found a better energy. Increase step for next iteration
       * if we are still going down, decrease it otherwise
       */
      if(gpc<0)
	stepsize *= 1.618034;  /* The golden section */
      else
	stepsize *= 0.618034;  /* 1/golden section */
    } else {
      /* New energy is the same or higher. We will have to do some work
       * to find a smaller value in the interval. Take smaller step next time!
       */
      foundlower = FALSE;
      stepsize *= 0.618034;
    }    
    
    /* OK, if we didn't find a lower value we will have to locate one now - there must
     * be one in the interval [a=0,c].
     * The same thing is valid here, though: Don't spend dozens of iterations to find
     * the line minimum. We try to interpolate based on the derivative at the endpoints,
     * and only continue until we find a lower value. In most cases this means 1-2 iterations.
     *
     * I also have a safeguard for potentially really patological functions so we never
     * take more than 20 steps before we give up ...
     *
     * If we already found a lower value we just skip this step and continue to the update.
     */

    if(!foundlower) {
     
      nminstep=0;
      do {
	/* Select a new trial point.
	 * If the derivatives at points a & c have different sign we interpolate to zero,
	 * otherwise just do a bisection.
	 */
	
	if(gpa<0 && gpc>0)
	  b = a + gpa*(a-c)/(gpc-gpa);
	else
	  b = 0.5*(a+c);		
	
	/* safeguard if interpolation close to machine accuracy causes errors:
	 * never go outside the interval
	 */
	if(b<=a || b>=c)
	  b = 0.5*(a+c);
	
	/* Take a trial step */
	for (i=0; i<n; i++) 
	  xb[i] = lastx[i] + b*s[i];
	
	neval++;
	/* Calculate energy for the trial step */
	ems.s.x = (rvec *)xb;
	ems.f   = (rvec *)fb;
	evaluate_energy(fplog,bVerbose,cr,
			state,top,&ems,&buf,top,
			inputrec,grps,nrnb,wcycle,
			vsite,constr,fcd,graph,mdatoms,fr,mu_tot,ener,
			step,FALSE);
	EpotB = ems.epot;
	
	fnorm = ems.fnorm;
	
	for(gpb=0,i=0; i<n; i++) 
	  gpb -= s[i]*fb[i];   /* f is negative gradient, thus the sign */
	
	/* Sum the gradient along the line across CPUs */
	if (PAR(cr))
	  gmx_sumd(1,&gpb,cr);
	
	/* Keep one of the intervals based on the value of the derivative at the new point */
	if(gpb>0) {
	  /* Replace c endpoint with b */
	  EpotC = EpotB;
	  c = b;
	  gpc = gpb;
	  /* swap coord pointers b/c */
	  xtmp = xb; 
	  ftmp = fb;
	  xb = xc; 
	  fb = fc;
	  xc = xtmp;
	  fc = ftmp;
	} else {
	  /* Replace a endpoint with b */
	  EpotA = EpotB;
	  a = b;
	  gpa = gpb;
	  /* swap coord pointers a/b */
	  xtmp = xb; 
	  ftmp = fb;
	  xb = xa; 
	  fb = fa;
	  xa = xtmp; 
	  fa = ftmp;
	}
	
	/* 
	 * Stop search as soon as we find a value smaller than the endpoints,
	 * or if the tolerance is below machine precision.
	 * Never run more than 20 steps, no matter what.
	 */
	nminstep++; 
      } while((EpotB>EpotA || EpotB>EpotC) && (nminstep<20));

      if(fabs(EpotB-Epot0)<GMX_REAL_EPS || nminstep>=20) {
	/* OK. We couldn't find a significantly lower energy.
	 * If ncorr==0 this was steepest descent, and then we give up.
	 * If not, reset memory to restart as steepest descent before quitting.
         */
	if(ncorr==0) {
	/* Converged */
	  converged=TRUE;
	  break;
	} else {
	  /* Reset memory */
	  ncorr=0;
	  /* Search in gradient direction */
	  for(i=0;i<n;i++)
	    dx[point][i]=ff[i];
	  /* Reset stepsize */
	  stepsize = 1.0/fnorm;
	  continue;
	}
      }
      
      /* Select min energy state of A & C, put the best in xx/ff/Epot
       */
      if(EpotC<EpotA) {
	Epot = EpotC;
	/* Use state C */
	for(i=0;i<n;i++) {
	  xx[i]=xc[i];
	  ff[i]=fc[i];
	}
	stepsize=c;
      } else {
	Epot = EpotA;
	/* Use state A */
	for(i=0;i<n;i++) {
	  xx[i]=xa[i];
	  ff[i]=fa[i];
	}
	stepsize=a;
      }
      
    } else {
      /* found lower */
      Epot = EpotC;
      /* Use state C */
      for(i=0;i<n;i++) {
	xx[i]=xc[i];
	ff[i]=fc[i];
      }
      stepsize=c;
    }

    /* Update the memory information, and calculate a new 
     * approximation of the inverse hessian 
     */
    
    /* Have new data in Epot, xx, ff */	
    if(ncorr<nmaxcorr)
      ncorr++;

    for(i=0;i<n;i++) {
      dg[point][i]=lastf[i]-ff[i];
      dx[point][i]*=stepsize;
    }
    
    dgdg=0;
    dgdx=0;	
    for(i=0;i<n;i++) {
      dgdg+=dg[point][i]*dg[point][i];
      dgdx+=dg[point][i]*dx[point][i];
    }
    
    diag=dgdx/dgdg;
    
    rho[point]=1.0/dgdx;
    point++;
    
    if(point>=nmaxcorr)
      point=0;
    
    /* Update */
    for(i=0;i<n;i++)
      p[i]=ff[i];
    
    cp=point;
    
    /* Recursive update. First go back over the memory points */
    for(k=0;k<ncorr;k++) {
      cp--;
      if(cp<0) 
	cp=ncorr-1;
      
      sq=0;
      for(i=0;i<n;i++)
	sq+=dx[cp][i]*p[i];
      
      alpha[cp]=rho[cp]*sq;
      
      for(i=0;i<n;i++)
	p[i] -= alpha[cp]*dg[cp][i];		
    }
    
    for(i=0;i<n;i++)
      p[i] *= diag;
    
    /* And then go forward again */
    for(k=0;k<ncorr;k++) {
      yr = 0;
      for(i=0;i<n;i++)
	yr += p[i]*dg[cp][i];
      
      beta = rho[cp]*yr;	    
      beta = alpha[cp]-beta;
      
      for(i=0;i<n;i++)
	p[i] += beta*dx[cp][i];
      
      cp++;	
      if(cp>=ncorr)
	cp=0;
    }
    
    for(i=0;i<n;i++)
      if(!frozen[i])
	dx[point][i] = p[i];
      else
	dx[point][i] = 0;

    stepsize=1.0;
    
    /* Test whether the convergence criterion is met */
    get_f_norm_max(cr,&(inputrec->opts),mdatoms,f,&fnorm,&fmax,&nfmax);
    
    /* Print it if necessary */
    if (MASTER(cr)) {
      ener[F_ETOT] = ener[F_EPOT]; /* No kinetic energy */
      if(bVerbose)
	fprintf(stderr,"\rStep %d, Epot=%12.6e, Fnorm=%9.3e, Fmax=%9.3e (atom %d)\n",
		step,Epot,fnorm/sqrt(state->natoms),fmax,nfmax+1);
      /* Store the new (lower) energies */
      upd_mdebin(mdebin,NULL,TRUE,mdatoms->tmass,step,(real)step,
		 ener,state,state->box,shake_vir,
		 force_vir,vir,pres,grps,mu_tot,constr);
      do_log = do_per_step(step,inputrec->nstlog);
      do_ene = do_per_step(step,inputrec->nstenergy);
      if(do_log)
	print_ebin_header(fplog,step,step,state->lambda);
      print_ebin(fp_ene,do_ene,FALSE,FALSE,
		 do_log ? fplog : NULL,step,step,step,eprNORMAL,
		 TRUE,mdebin,fcd,&(top->atoms),&(inputrec->opts));
    }
    
    /* Stop when the maximum force lies below tolerance.
     * If we have reached machine precision, converged is already set to true.
     */
    
    converged = converged || (fmax < inputrec->em_tol);
    
  } /* End of the loop */
  
  if(converged)	
    step--; /* we never took that last step in this case */
  
  if(fmax>inputrec->em_tol) {
    if (MASTER(cr)) {
      warn_step(stderr,inputrec->em_tol,FALSE);
      warn_step(fplog,inputrec->em_tol,FALSE);
    }
    converged = FALSE; 
  }
  
  /* If we printed energy and/or logfile last step (which was the last step)
   * we don't have to do it again, but otherwise print the final values.
   */
  if(!do_log) /* Write final value to log since we didn't do anythin last step */
    print_ebin_header(fplog,step,step,state->lambda);
  if(!do_ene || !do_log) /* Write final energy file entries */
    print_ebin(fp_ene,!do_ene,FALSE,FALSE,
	       !do_log ? fplog : NULL,step,step,step,eprNORMAL,
	       TRUE,mdebin,fcd,&(top->atoms),&(inputrec->opts));
  
  /* Print some stuff... */
  if (MASTER(cr))
    fprintf(stderr,"\nwriting lowest energy coordinates.\n");
  
	/* IMPORTANT!
	 * For accurate normal mode calculation it is imperative that we
	 * store the last conformation into the full precision binary trajectory
	 .
	 *
	 * However, we should only do it if we did NOT already write this step
	 * above (which we did if do_x or do_f was true).
	 */  
	do_x = !do_per_step(step,inputrec->nstxout);
    do_f = !do_per_step(step,inputrec->nstfout);
	write_traj(fplog,cr,fp_trn,do_x,FALSE,do_f,0,FALSE,0,NULL,FALSE,
	     top,inputrec->eI,step,(real)step,state,state,f,f);

  if (MASTER(cr))
    write_sto_conf(ftp2fn(efSTO,nfile,fnm),
		   *top->name, &(top->atoms),state->x,NULL,
		   inputrec->ePBC,state->box);
  
  if (MASTER(cr)) {
    print_converged(stderr,LBFGS,inputrec->em_tol,step,converged,
		    number_steps,Epot,fmax,nfmax,fnorm/sqrt(state->natoms));
    print_converged(fplog,LBFGS,inputrec->em_tol,step,converged,
		    number_steps,Epot,fmax,nfmax,fnorm/sqrt(state->natoms));
    
    fprintf(fplog,"\nPerformed %d energy evaluations in total.\n",neval);
  }
  
  finish_em(fplog,cr,fp_trn,fp_ene);

  /* To print the actual number of steps we needed somewhere */
  inputrec->nsteps=step;

  return start_t;
} /* That's all folks */


time_t do_steep(FILE *fplog,t_commrec *cr,
		int nfile,t_filenm fnm[],
		bool bVerbose,bool bCompact,
		gmx_vsite_t *vsite,gmx_constr_t constr,
		int stepout,
		t_inputrec *inputrec,t_groups *grps,
		t_topology *top_global,real ener[],t_fcdata *fcd,
		t_state *state_global,rvec f[],
		rvec buf[],t_mdatoms *mdatoms,
		t_nrnb *nrnb,gmx_wallcycle_t wcycle,
		t_graph *graph,t_edsamyn *edyn,
		t_forcerec *fr,
		int repl_ex_nst,int repl_ex_seed,
		real cpt_period,real max_hours,
		unsigned long Flags)
{ 
  const char *SD="Steepest Descents";
  em_state_t *s_min,*s_try;
  rvec       *f_global;
  t_topology *top;
  real   stepsize,constepsize;
  real   ustep,dvdlambda,fnormn;
  int        fp_trn,fp_ene; 
  t_mdebin   *mdebin; 
  bool   bDone,bAbort,do_x,do_f; 
  time_t start_t; 
  tensor force_vir,shake_vir,pres,vir; 
  rvec   mu_tot;
  int    nsteps;
  int    count=0; 
  int    steps_accepted=0; 
  /* not used */
  real   terminate=0;

  s_min = init_em_state();
  s_try = init_em_state();

  /* Init em and store the local state in s_try */
  init_em(fplog,SD,cr,inputrec,grps,
	  state_global,top_global,s_try,&top,f,&buf,&f_global,
	  nrnb,mu_tot,fr,graph,mdatoms,vsite,constr,
	  nfile,fnm,&fp_trn,&fp_ene,&mdebin);

  /* Print to log file  */
  start_t=print_date_and_time(fplog,cr->nodeid,"Started Steepest Descents");
  wallcycle_start(wcycle,ewcRUN);

  /* Clear some matrix variables  */
  clear_mat(shake_vir); 
  clear_mat(vir); 
  clear_mat(pres); 
    
  /* Set variables for stepsize (in nm). This is the largest  
   * step that we are going to make in any direction. 
   */
  ustep = inputrec->em_stepsize; 
  stepsize = 0;
  
  /* Max number of steps  */
  nsteps = inputrec->nsteps; 
  
  if (MASTER(cr)) 
    /* Print to the screen  */
    sp_header(stderr,SD,inputrec->em_tol,nsteps);
  if (fplog)
    sp_header(fplog,SD,inputrec->em_tol,nsteps);
    
  /**** HERE STARTS THE LOOP ****
   * count is the counter for the number of steps 
   * bDone will be TRUE when the minimization has converged
   * bAbort will be TRUE when nsteps steps have been performed or when
   * the stepsize becomes smaller than is reasonable for machine precision
   */
  count  = 0;
  bDone  = FALSE;
  bAbort = FALSE;
  while( !bDone && !bAbort ) {
    bAbort = (nsteps > 0) && (count==nsteps);
    
    /* set new coordinates, except for first step */
    if (count > 0) {
      do_em_step(cr,inputrec,mdatoms,s_min,stepsize,s_min->f,s_try,
		 constr,top,nrnb,count);
    }
    
    evaluate_energy(fplog,bVerbose,cr,
		    state_global,top_global,s_try,&buf,top,
		    inputrec,grps,nrnb,wcycle,
		    vsite,constr,fcd,graph,mdatoms,fr,mu_tot,ener,
		    count,count==0);
    
    if (MASTER(cr))
      print_ebin_header(fplog,count,count,s_try->s.lambda);

    if (count == 0)
      s_min->epot = s_try->epot + 1;
    
    /* Print it if necessary  */
    if (MASTER(cr)) {
      if (bVerbose) {
	fprintf(stderr,"Step=%5d, Dmax= %6.1e nm, Epot= %12.5e Fmax= %11.5e, atom= %d%c",
		count,ustep,s_try->epot,s_try->fmax,s_try->a_fmax+1,
		(s_try->epot < s_min->epot) ? '\n' : '\r');
      }
      
      if (s_try->epot < s_min->epot) {
	/* Store the new (lower) energies  */
	upd_mdebin(mdebin,NULL,TRUE,mdatoms->tmass,count,(real)count,
		   ener,&s_try->s,s_try->s.box,shake_vir, 
		   force_vir,vir,pres,grps,mu_tot,constr);
	print_ebin(fp_ene,TRUE,
		   do_per_step(steps_accepted,inputrec->nstdisreout),
		   do_per_step(steps_accepted,inputrec->nstorireout),
		   fplog,count,count,count,eprNORMAL,TRUE,
		   mdebin,fcd,&(top->atoms),&(inputrec->opts));
	fflush(fplog);
      }
    } 
    
    /* Now if the new energy is smaller than the previous...  
     * or if this is the first step!
     * or if we did random steps! 
     */
    
    if ( (count==0) || (s_try->epot < s_min->epot) ) {
      steps_accepted++; 

      /* Test whether the convergence criterion is met...  */
      bDone = (s_try->fmax < inputrec->em_tol);
      
      /* Copy the arrays for force, positions and energy  */
      /* The 'Min' array always holds the coords and forces of the minimal 
	 sampled energy  */
      swap_em_state(s_min,s_try);
      if (count > 0)
	ustep *= 1.2;

      /* Write to trn, if necessary */
      do_x = do_per_step(steps_accepted,inputrec->nstxout);
      do_f = do_per_step(steps_accepted,inputrec->nstfout);
      write_traj(fplog,cr,fp_trn,do_x,FALSE,do_f,0,FALSE,0,NULL,FALSE,
		 top_global,inputrec->eI,count,(real)count,
		 &s_min->s,state_global,s_min->f,f_global);
    } 
    else {
      /* If energy is not smaller make the step smaller...  */
      ustep *= 0.5;

      if (DOMAINDECOMP(cr) && s_min->s.ddp_count != cr->dd->ddp_count) {
	/* Reload the old state */
	em_dd_partition_system(fplog,count,cr,top_global,inputrec,
			       s_min,&buf,top,mdatoms,fr,vsite,constr,
			       nrnb,wcycle);
      }
    }
    
    /* Determine new step  */
    stepsize = ustep/s_min->fmax;
    
    /* Check if stepsize is too small, with 1 nm as a characteristic length */
#ifdef GMX_DOUBLE
    if (ustep < 1e-12)
#else
    if (ustep < 1e-6)
#endif
      {
	if (MASTER(cr)) {
	  warn_step(stderr,inputrec->em_tol,constr!=NULL);
	  warn_step(fplog,inputrec->em_tol,constr!=NULL);
	}
	bAbort=TRUE;
      }
    
    count++;
  } /* End of the loop  */
  
  /* In parallel write_traj copies back the coordinates,
   * otherwise we have to do it explicitly.
   */
  if (!PAR(cr)) {
    copy_em_coords_back(s_min,state_global);
  }

    /* Print some shit...  */
  if (MASTER(cr)) 
    fprintf(stderr,"\nwriting lowest energy coordinates.\n"); 
  write_traj(fplog,cr,fp_trn,TRUE,FALSE,inputrec->nstfout,0,FALSE,0,NULL,FALSE,
	     top_global,inputrec->eI,count,(real)count,
	     &s_min->s,state_global,s_min->f,f_global);

  fnormn = s_min->fnorm/sqrt(state_global->natoms);

  if (MASTER(cr)) {
    write_sto_conf(ftp2fn(efSTO,nfile,fnm),
		   *top_global->name,&(top_global->atoms),
		   state_global->x,NULL,inputrec->ePBC,state_global->box);
    
    print_converged(stderr,SD,inputrec->em_tol,count,bDone,nsteps,
		    s_min->epot,s_min->fmax,s_min->a_fmax,fnormn);
    print_converged(fplog,SD,inputrec->em_tol,count,bDone,nsteps,
		    s_min->epot,s_min->fmax,s_min->a_fmax,fnormn);
  }

  finish_em(fplog,cr,fp_trn,fp_ene);
  
  /* To print the actual number of steps we needed somewhere */
  inputrec->nsteps=count;
  
  return start_t;
} /* That's all folks */


time_t do_nm(FILE *fplog,t_commrec *cr,
	     int nfile,t_filenm fnm[],
	     bool bVerbose,bool bCompact,
	     gmx_vsite_t *vsite,gmx_constr_t constr,
	     int stepout,
	     t_inputrec *inputrec,t_groups *grps,
	     t_topology *top_global,real ener[],t_fcdata *fcd,
	     t_state *state_global,rvec f[],
	     rvec buf[],t_mdatoms *mdatoms,
	     t_nrnb *nrnb,gmx_wallcycle_t wcycle,
	     t_graph *graph,t_edsamyn *edyn,
	     t_forcerec *fr,
	     int repl_ex_nst,int repl_ex_seed,
	     real cpt_period,real max_hours,
	     unsigned long Flags)
{
    t_mdebin   *mdebin;
	const char *NM = "Normal Mode Analysis";
    int        fp_ene,step,i;
    time_t     start_t;
	rvec       *f_global;
    t_topology *top;
	real       t,lambda,t0,lam0;
    bool       bNS;
    tensor     force_vir;
    int        nfmax,count;
    rvec       mu_tot;
    rvec       *fneg,*fpos;
    bool       bSparse; /* use sparse matrix storage format */
    size_t     sz;
    gmx_sparsematrix_t * sparse_matrix = NULL;
    real *     full_matrix             = NULL;
	em_state_t *   state_work;
	
    /* added with respect to mdrun */
    int        idum,jdum,kdum,row,col;
    real       der_range=10.0*sqrt(GMX_REAL_EPS);
    real       fnorm,fmax;
    real       dfdx;
    
	state_work = init_em_state();

	/* Init em and store the local state in state_minimum */
	init_em(fplog,NM,cr,inputrec,grps,
			state_global,top_global,state_work,&top,f,&buf,&f_global,
			nrnb,mu_tot,fr,graph,mdatoms,vsite,constr,
			nfile,fnm,NULL,&fp_ene,&mdebin);
	
    snew(fneg,top_global->atoms.nr);
    snew(fpos,top_global->atoms.nr);
    
#ifndef GMX_DOUBLE
    fprintf(stderr,
            "NOTE: This version of Gromacs has been compiled in single precision,\n"
            "      which MIGHT not be accurate enough for normal mode analysis.\n"
            "      Gromacs now uses sparse matrix storage, so the memory requirements\n"
            "      are fairly modest even if you recompile in double precision.\n\n");
#endif
    
    /* Check if we can/should use sparse storage format.
     *
     * Sparse format is only useful when the Hessian itself is sparse, which it
      * will be when we use a cutoff.    
      * For small systems (n<1000) it is easier to always use full matrix format, though.
      */
    if(EEL_FULL(fr->eeltype) || fr->rlist==0.0)
    {
        fprintf(stderr,"Non-cutoff electrostatics used, forcing full Hessian format.\n");
        bSparse = FALSE;
    }
    else if(top_global->atoms.nr < 1000)
    {
        fprintf(stderr,"Small system size (N=%d), using full Hessian format.\n",top_global->atoms.nr);
        bSparse = FALSE;
    }
    else
    {
        fprintf(stderr,"Using compressed symmetric sparse Hessian format.\n");
        bSparse = TRUE;
    }
    
    sz = DIM*top_global->atoms.nr;
    
    fprintf(stderr,"Allocating Hessian memory...\n\n");

    if(bSparse)
    {
        sparse_matrix=gmx_sparsematrix_init(sz);
        sparse_matrix->compressed_symmetric = TRUE;
    }
    else
    {
        snew(full_matrix,sz*sz);
    }
    
    /* Initial values */
    t0           = inputrec->init_t;
    lam0         = inputrec->init_lambda;
    t            = t0;
    lambda       = lam0;
    
    init_nrnb(nrnb);
    
    where();
    
    /* Write start time and temperature */
    start_t=print_date_and_time(fplog,cr->nodeid,"Started nmrun");
    wallcycle_start(wcycle,ewcRUN);
    if (MASTER(cr)) 
    {
        fprintf(stderr,"starting normal mode calculation '%s'\n%d steps.\n\n",*(top_global->name),
                top_global->atoms.nr);
    }
    
	count = 0;
    evaluate_energy(fplog,bVerbose,cr,
					state_global,top_global,state_work,&buf,top,
					inputrec,grps,nrnb,wcycle,
					vsite,constr,fcd,graph,mdatoms,fr,mu_tot,ener,
					count,count==0);
	count++;
    /* if forces are not small, warn user */
    get_state_f_norm_max(cr,&(inputrec->opts),mdatoms,state_work);

    fprintf(stderr,"Maximum force:%12.5e\n",state_work->fmax);
    if (state_work->fmax > 1.0e-3) 
    {
        fprintf(stderr,"Maximum force probably not small enough to");
        fprintf(stderr," ensure that you are in an \nenergy well. ");
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
    
    inputrec->nsteps=top_global->atoms.nr;

    for (step=0; (step<inputrec->nsteps); step++) 
    {
        
        for (idum=0; (idum<DIM); idum++) 
        {
            row = DIM*step+idum;
			
			state_work->s.x[step][idum] -= der_range;
            
			evaluate_energy(fplog,bVerbose,cr,
							state_global,top_global,state_work,&buf,top,
							inputrec,grps,nrnb,wcycle,
							vsite,constr,fcd,graph,mdatoms,fr,mu_tot,ener,
							count,count==0);
			count++;
			
			for ( i=0 ; i < top_global->atoms.nr ; i++ )
			{
				copy_rvec ( state_work->f[i] , fneg[i] );
			}
			
			state_work->s.x[step][idum] += 2*der_range;
            
			evaluate_energy(fplog,bVerbose,cr,
							state_global,top_global,state_work,&buf,top,
							inputrec,grps,nrnb,wcycle,
							vsite,constr,fcd,graph,mdatoms,fr,mu_tot,ener,
							count,count==0);
			count++;

			for ( i=0 ; i < top_global->atoms.nr ; i++ )
			{
				copy_rvec ( state_work->f[i] , fpos[i] );
			}
						
            for (jdum=0; (jdum<top_global->atoms.nr); jdum++) 
            {
                for (kdum=0; (kdum<DIM); kdum++) 
                {
                    dfdx=-(fpos[jdum][kdum]-fneg[jdum][kdum])/(2*der_range);
                    col = DIM*jdum+kdum;
                    
                    if(bSparse)
                    {
                        if(col>=row && dfdx!=0.0)
                            gmx_sparsematrix_increment_value(sparse_matrix,row,col,dfdx);
                    }
                    else
                    {
                        full_matrix[row*sz+col] = dfdx;
                    }
                }
            }
            
            /* x is restored to original */
			state_work->s.x[step][idum] -= der_range;
            
            if (bVerbose && fplog)
                fflush(fplog);            
        }
        /* write progress */
        if (MASTER(cr) && bVerbose) 
        {
            fprintf(stderr,"\rFinished step %d out of %d",step+1,top_global->atoms.nr); 
            fflush(stderr);
        }
    }
    t=t0+step*inputrec->delta_t;
    lambda=lam0+step*inputrec->delta_lambda;
    
    if (MASTER(cr)) 
    {
        print_ebin(-1,FALSE,FALSE,FALSE,fplog,step,step,t,eprAVER,
                   FALSE,mdebin,fcd,&(top_global->atoms),&(inputrec->opts));
        print_ebin(-1,FALSE,FALSE,FALSE,fplog,step,step,t,eprRMS,
                   FALSE,mdebin,fcd,&(top_global->atoms),&(inputrec->opts));
    }
      
    fprintf(stderr,"\n\nWriting Hessian...\n");
    gmx_mtxio_write(ftp2fn(efMTX,nfile,fnm),sz,sz,full_matrix,sparse_matrix);
    
    
    return start_t;
}


static void global_max(t_commrec *cr,int *n)
{
  int *sum,i;

  snew(sum,cr->nnodes);
  sum[cr->nodeid] = *n;
  gmx_sumi(cr->nnodes,sum,cr);
  for(i=0; i<cr->nnodes; i++)
    *n = max(*n,sum[i]);
  
  sfree(sum);
}

static void realloc_bins(double **bin,int *nbin,int nbin_new)
{
  int i;

  if (nbin_new != *nbin) {
    srenew(*bin,nbin_new);
    for(i=*nbin; i<nbin_new; i++)
      (*bin)[i] = 0;
    *nbin = nbin_new;
  }
}

time_t do_tpi(FILE *fplog,t_commrec *cr,
	      int nfile,t_filenm fnm[],
	      bool bVerbose,bool bCompact,
	      gmx_vsite_t *vsite,gmx_constr_t constr,
	      int stepout,
	      t_inputrec *inputrec,t_groups *grps,
	      t_topology *top,real ener[],t_fcdata *fcd,
	      t_state *state,rvec f[],
	      rvec buf[],t_mdatoms *mdatoms,
	      t_nrnb *nrnb,gmx_wallcycle_t wcycle,
	      t_graph *graph,t_edsamyn *edyn,
	      t_forcerec *fr,
	      int repl_ex_nst,int repl_ex_seed,
	      real cpt_period,real max_hours,
	      unsigned long Flags)
{
  const char *TPI="Test Particle Insertion"; 
  real   lambda,t,temp,beta,drmax;
  double embU,sum_embU,*sum_UgembU,V,V_all,VembU_all;
  int    status;
  t_trxframe rerun_fr;
  bool   bDispCorr,bCharge,bRFExcl,bNotLastFrame,bStateChanged,bNS,bOurStep;
  time_t start_t; 
  tensor force_vir,shake_vir,vir,pres;
  int    cg_tp,a_tp0,a_tp1,ngid,gid_tp,nener,e;
  rvec   *x_mol;
  rvec   mu_tot,x_init,dx,x_tp;
  int    nnodes,frame,nsteps,step;
  int    i,start,end;
  static gmx_rng_t tpi_rand;
  FILE   *fp_tpi=NULL;
  char   *ptr,*dump_pdb,**leg,str[STRLEN],str2[STRLEN];
  double dbl,dump_ener;
  bool   bCavity;
  int    nat_cavity=0,d;
  real   *mass_cavity=NULL,mass_tot;
  int    nbin;
  double invbinw,*bin,refvolshift,logV,bUlogV;
  char   *tpid_leg[2]={"direct","reweighted"};

  /* Since numerical problems can lead to extreme negative energies
   * when atoms overlap, we need to set a lower limit for beta*U.
   */
  real bU_neg_limit = -50;

  /* Since there is no upper limit to the insertion energies,
   * we need to set an upper limit for the distribution output.
   */
  real bU_bin_limit      = 50;
  real bU_logV_bin_limit = bU_bin_limit + 10;

  nnodes = cr->nnodes;

  bCavity = (inputrec->eI == eiTPIC);
  if (bCavity) {
    ptr = getenv("GMX_TPIC_MASSES");
    if (ptr == NULL) {
      nat_cavity = 1;
    } else {
      /* Read (multiple) masses from env var GMX_TPIC_MASSES,
       * The center of mass of the last atoms is then used for TPIC.
       */
      nat_cavity = 0;
      while (sscanf(ptr,"%lf%n",&dbl,&i) > 0) {
	srenew(mass_cavity,nat_cavity+1);
	mass_cavity[nat_cavity] = dbl;
	fprintf(fplog,"mass[%d] = %f\n",
		nat_cavity+1,mass_cavity[nat_cavity]);
	nat_cavity++;
	ptr += i;
      }
      if (nat_cavity == 0)
	gmx_fatal(FARGS,"Found %d masses in GMX_TPIC_MASSES",nat_cavity);
    }
  }

  /*
  init_em(fplog,TPI,inputrec,&lambda,nrnb,mu_tot,
  state->box,fr,mdatoms,top,cr,nfile,fnm,NULL,NULL);*/
  /* We never need full pbc for TPI */
  fr->ePBC = epbcXYZ;
  /* Determine the temperature for the Boltzmann weighting */
  temp = inputrec->opts.ref_t[0];
  if (fplog) {
    for(i=1; (i<inputrec->opts.ngtc); i++) {
      if (inputrec->opts.ref_t[i] != temp) {
	fprintf(fplog,"\nWARNING: The temperatures of the different temperature coupling groups are not identical\n\n");
	fprintf(stderr,"\nWARNING: The temperatures of the different temperature coupling groups are not identical\n\n");
      }
    }
    fprintf(fplog,
	    "\n  The temperature for test particle insertion is %.3f K\n\n",
	    temp);
  }
  beta = 1.0/(BOLTZ*temp);

  /* Number of insertions per frame */
  nsteps = inputrec->nsteps; 

  /* Use the same neighborlist with more insertions points
   * in a sphere of radius drmax around the initial point
   */
  /* This should be a proper mdp parameter */
  drmax = inputrec->rtpi;

  /* An environment variable can be set to dump all configurations
   * to pdb with an insertion energy <= this value.
   */
  dump_pdb = getenv("GMX_TPI_DUMP");
  dump_ener = 0;
  if (dump_pdb)
    sscanf(dump_pdb,"%lf",&dump_ener);

  atoms2md(&top->atoms,inputrec,top->idef.il[F_ORIRES].nr,0,NULL,
	   0,top->atoms.nr,mdatoms);
  update_mdatoms(mdatoms,inputrec->init_lambda);

  /* Print to log file  */
  start_t=print_date_and_time(fplog,cr->nodeid,
			      "Started Test Particle Insertion"); 
  wallcycle_start(wcycle,ewcRUN);

  /* The last charge group is the group to be inserted */
  cg_tp = top->cgs.nr - 1;
  a_tp0 = top->cgs.index[cg_tp];
  a_tp1 = top->cgs.index[cg_tp+1];
  if (debug)
    fprintf(debug,"TPI cg %d, atoms %d-%d\n",cg_tp,a_tp0,a_tp1);
  if (a_tp1 - a_tp0 > 1 &&
      (inputrec->rlist < inputrec->rcoulomb ||
       inputrec->rlist < inputrec->rvdw))
    gmx_fatal(FARGS,"Can not do TPI for multi-atom molecule with a twin-range cut-off");
  snew(x_mol,a_tp1-a_tp0);

  bDispCorr = (inputrec->eDispCorr != edispcNO);
  bCharge = FALSE;
  for(i=a_tp0; i<a_tp1; i++) {
    /* Copy the coordinates of the molecule to be insterted */
    copy_rvec(state->x[i],x_mol[i-a_tp0]);
    /* Check if we need to print electrostatic energies */
    bCharge |= (top->atoms.atom[i].q!=0 || top->atoms.atom[i].qB!=0);
  }
  bRFExcl = (bCharge && EEL_RF(fr->eeltype) && fr->eeltype!=eelRF_NEC);

  calc_cgcm(fplog,cg_tp,cg_tp+1,&(top->cgs),state->x,fr->cg_cm);
  if (bCavity) {
    if (norm(fr->cg_cm[cg_tp]) > 0.5*inputrec->rlist && fplog) {
      fprintf(fplog, "WARNING: Your TPI molecule is not centered at 0,0,0\n");
      fprintf(stderr,"WARNING: Your TPI molecule is not centered at 0,0,0\n");
    }
  } else {
    /* Center the molecule to be inserted at zero */
     for(i=0; i<a_tp1-a_tp0; i++)
      rvec_dec(x_mol[i],fr->cg_cm[cg_tp]);
  }

  if (fplog) {
    fprintf(fplog,"\nWill insert %d atoms %s partial charges\n",
	    a_tp1-a_tp0,bCharge ? "with" : "without");
    
    fprintf(fplog,"\nWill insert %d times in each frame of %s\n",
	    nsteps,opt2fn("-rerun",nfile,fnm));
  }
  
  if (!bCavity) {
    if (inputrec->nstlist > 1) {
      if (drmax==0 && a_tp1-a_tp0==1) {
	gmx_fatal(FARGS,"Re-using the neighborlist %d times for insertions of a single atom in a sphere of radius %f does not make sense",inputrec->nstlist,drmax);
      }
      if (fplog)
	fprintf(fplog,"Will use the same neighborlist for %d insertions in a sphere of radius %f\n",inputrec->nstlist,drmax);
    }
  } else {
    if (fplog)
      fprintf(fplog,"Will insert randomly in a sphere of radius %f around the center of the cavity\n",drmax);
  }

  ngid = top->atoms.grps[egcENER].nr;
  gid_tp = GET_CGINFO_GID(fr->cginfo[cg_tp]);
  nener = 1 + ngid;
  if (bDispCorr)
    nener += 1;
  if (bCharge) {
    nener += ngid;
    if (bRFExcl)
      nener += 1;
  }
  snew(sum_UgembU,nener);

  /* Initialize random generator */
  tpi_rand = gmx_rng_init(inputrec->ld_seed);

  if (MASTER(cr)) {
    fp_tpi = xvgropen(opt2fn("-tpi",nfile,fnm),
		      "TPI energies","Time (ps)",
		      "(kJ mol\\S-1\\N) / (nm\\S3\\N)");
    xvgr_subtitle(fp_tpi,"f. are averages over one frame");
    snew(leg,4+nener);
    e = 0;
    sprintf(str,"-kT log(<Ve\\S-\\8b\\4U\\N>/<V>)");
    leg[e++] = strdup(str);
    sprintf(str,"f. -kT log<e\\S-\\8b\\4U\\N>");
    leg[e++] = strdup(str);
    sprintf(str,"f. <e\\S-\\8b\\4U\\N>");
    leg[e++] = strdup(str);
    sprintf(str,"f. V");
    leg[e++] = strdup(str);
    sprintf(str,"f. <Ue\\S-\\8b\\4U\\N>");
    leg[e++] = strdup(str);
    for(i=0; i<ngid; i++) {
      sprintf(str,"f. <U\\sVdW %s\\Ne\\S-\\8b\\4U\\N>",
	      *(top->atoms.grpname[top->atoms.grps[egcENER].nm_ind[i]]));
      leg[e++] = strdup(str);
    }
    if (bDispCorr) {
      sprintf(str,"f. <U\\sdisp c\\Ne\\S-\\8b\\4U\\N>");
      leg[e++] = strdup(str);
    }
    if (bCharge) {
      for(i=0; i<ngid; i++) {
	sprintf(str,"f. <U\\sCoul %s\\Ne\\S-\\8b\\4U\\N>",
		*(top->atoms.grpname[top->atoms.grps[egcENER].nm_ind[i]]));
	leg[e++] = strdup(str);
      }
      if (bRFExcl) {
	sprintf(str,"f. <U\\sRF excl\\Ne\\S-\\8b\\4U\\N>");
	leg[e++] = strdup(str);
      }
    }
    xvgr_legend(fp_tpi,4+nener,leg);
    for(i=0; i<4+nener; i++)
      sfree(leg[i]);
    sfree(leg);
  }
  clear_rvec(x_init);
  V_all = 0;
  VembU_all = 0;

  invbinw = 10;
  nbin = 10;
  snew(bin,nbin);

  start_time();

  bNotLastFrame = read_first_frame(&status,opt2fn("-rerun",nfile,fnm),
				   &rerun_fr,TRX_NEED_X);
  frame = 0;

  if (rerun_fr.natoms - (bCavity ? nat_cavity : 0) !=
      mdatoms->nr - (a_tp1 - a_tp0))
    gmx_fatal(FARGS,"Number of atoms in trajectory (%d)%s "
	      "is not equal the number in the run input file (%d) "
	      "minus the number of atoms to insert (%d)\n",
	      rerun_fr.natoms,bCavity ? " minus one" : "",
	      mdatoms->nr,a_tp1-a_tp0);

  refvolshift = log(det(rerun_fr.box));
  
  while (bNotLastFrame) {
    lambda = rerun_fr.lambda;
    t = rerun_fr.time;
    
    sum_embU = 0;
    for(e=0; e<nener; e++)
      sum_UgembU[e] = 0;

    /* Copy the coordinates from the input trajectory */
    for(i=0; i<rerun_fr.natoms; i++)
      copy_rvec(rerun_fr.x[i],state->x[i]);
      
    V = det(rerun_fr.box);
    logV = log(V);

    bStateChanged = TRUE;
    bNS = TRUE;
    for(step=0; step<nsteps; step++) {
      /* In parallel all nodes generate all random configurations.
       * In that way the result is identical to a single cpu tpi run.
       */
      if (!bCavity) {
	/* Random insertion in the whole volume */
	bNS = (step % inputrec->nstlist == 0);
	if (bNS) {
	  /* Generate a random position in the box */
          x_init[XX] = gmx_rng_uniform_real(tpi_rand)*state->box[XX][XX];
          x_init[YY] = gmx_rng_uniform_real(tpi_rand)*state->box[YY][YY];
          x_init[ZZ] = gmx_rng_uniform_real(tpi_rand)*state->box[ZZ][ZZ];
	}
	if (inputrec->nstlist == 1) {
	  copy_rvec(x_init,x_tp);
	} else {
	  /* Generate coordinates within |dx|=drmax of x_init */
	  do {
	    dx[XX] = (2*gmx_rng_uniform_real(tpi_rand) - 1)*drmax;
	    dx[YY] = (2*gmx_rng_uniform_real(tpi_rand) - 1)*drmax;
	    dx[ZZ] = (2*gmx_rng_uniform_real(tpi_rand) - 1)*drmax;
	  } while (norm2(dx) > drmax*drmax);
	  rvec_add(x_init,dx,x_tp);
	}
      } else {
	/* Random insertion around a cavity location
	 * given by the last coordinate of the trajectory.
	 */
	if (step == 0) {
	  if (nat_cavity == 1) {
	    /* Copy the location of the cavity */
	    copy_rvec(rerun_fr.x[rerun_fr.natoms-1],x_init);
	  } else {
	    /* Determine the center of mass of the last molecule */
	    clear_rvec(x_init);
	    mass_tot = 0;
	    for(i=0; i<nat_cavity; i++) {
	      for(d=0; d<DIM; d++)
		x_init[d] +=
		  mass_cavity[i]*rerun_fr.x[rerun_fr.natoms-nat_cavity+i][d];
	      mass_tot += mass_cavity[i];
	    }
	    for(d=0; d<DIM; d++)
	      x_init[d] /= mass_tot;
	  }
	}
	/* Generate coordinates within |dx|=drmax of x_init */
	do {
	  dx[XX] = (2*gmx_rng_uniform_real(tpi_rand) - 1)*drmax;
	  dx[YY] = (2*gmx_rng_uniform_real(tpi_rand) - 1)*drmax;
	  dx[ZZ] = (2*gmx_rng_uniform_real(tpi_rand) - 1)*drmax;
	} while (norm2(dx) > drmax*drmax);
	rvec_add(x_init,dx,x_tp);
      }

      if (a_tp1 - a_tp0 == 1) {
	/* Insert a single atom, just copy the insertion location */
	copy_rvec(x_tp,state->x[a_tp0]);
      } else {
	/* Copy the coordinates from the top file */
	for(i=a_tp0; i<a_tp1; i++)
	  copy_rvec(x_mol[i-a_tp0],state->x[i]);
	/* Rotate the molecule randomly */
	rotate_conf(a_tp1-a_tp0,state->x+a_tp0,NULL,
		    2*M_PI*gmx_rng_uniform_real(tpi_rand),
		    2*M_PI*gmx_rng_uniform_real(tpi_rand),
		    2*M_PI*gmx_rng_uniform_real(tpi_rand));
	/* Shift to the insertion location */
	for(i=a_tp0; i<a_tp1; i++)
	  rvec_inc(state->x[i],x_tp);
      }

      /* Check if this insertion belongs to this node */
      bOurStep = TRUE;
      if (PAR(cr)) {
	switch (inputrec->eI) {
	case eiTPI:
	  bOurStep = ((step / inputrec->nstlist) % nnodes == cr->nodeid);
	  break;
	case eiTPIC:
	  bOurStep = (step % nnodes == cr->nodeid);
	  break;
	default:
	  gmx_fatal(FARGS,"Unknown integrator %s",ei_names[inputrec->eI]);
	}
      }
      if (bOurStep) {
	/* Clear some matrix variables  */
	clear_mat(force_vir); 
	clear_mat(shake_vir);
	clear_mat(vir);
	clear_mat(pres);
	
	/* Set the charge group center of mass of the test particle */
	copy_rvec(x_init,fr->cg_cm[top->cgs.nr-1]);

	/* Calc energy (no forces) on new positions.
	 * Since we only need the intermolecular energy
	 * and the RF exclusion terms of the inserted molecule occur
	 * within a single charge group we can pass NULL for the graph.
	 * This also avoids shifts that would move charge groups
	 * out of the box.
	 */
	/* Make do_force do a single node fore calculation */
	cr->nnodes = 1;
	do_force(fplog,cr,inputrec,
		 step,nrnb,wcycle,top,grps,rerun_fr.box,state->x,f,
		 buf,force_vir,mdatoms,ener,fcd,
		 lambda,NULL,fr,NULL,mu_tot,t,NULL,NULL,
		 GMX_FORCE_NONBONDED |
		 (bNS ? GMX_FORCE_NS : 0) |
		 (bStateChanged ? GMX_FORCE_STATECHANGED : 0)); 
	cr->nnodes = nnodes;
	bStateChanged = FALSE;
	bNS = FALSE;

	/* Calculate long range corrections to pressure and energy */
	calc_dispcorr(fplog,inputrec,fr,step,mdatoms->nr,rerun_fr.box,lambda,
		      pres,vir,ener);
	
	/* If the compiler doesn't optimize this check away
	 * we catch the NAN energies.
	 * With tables extreme negative energies might occur close to r=0.
	 */
	if (ener[F_EPOT] != ener[F_EPOT] || ener[F_EPOT]*beta < bU_neg_limit) {
	  if (debug)
	    fprintf(debug,"\n  time %.3f, step %d: non-finite energy %f, using exp(-bU)=0\n",t,step,ener[F_EPOT]);
	  embU = 0;
	} else {
	  embU = exp(-beta*ener[F_EPOT]);
	  sum_embU += embU;
	  /* Determine the weighted energy contributions of each energy group */
	  e = 0;
	  sum_UgembU[e++] += ener[F_EPOT]*embU;
	  if (fr->bBHAM) {
	    for(i=0; i<ngid; i++)
	      sum_UgembU[e++] +=
		(grps->estat.ee[egBHAMSR][GID(i,gid_tp,ngid)] +
		 grps->estat.ee[egBHAMLR][GID(i,gid_tp,ngid)])*embU;
	  } else {
	    for(i=0; i<ngid; i++)
	      sum_UgembU[e++] +=
		(grps->estat.ee[egLJSR][GID(i,gid_tp,ngid)] +
		 grps->estat.ee[egLJLR][GID(i,gid_tp,ngid)])*embU;
	  }
	  if (bDispCorr)
	    sum_UgembU[e++] += ener[F_DISPCORR]*embU;
	  if (bCharge) {
	    for(i=0; i<ngid; i++)
	      sum_UgembU[e++] +=
		(grps->estat.ee[egCOULSR][GID(i,gid_tp,ngid)] +
		 grps->estat.ee[egCOULLR][GID(i,gid_tp,ngid)])*embU;
	    if (bRFExcl)
	      sum_UgembU[e++] += ener[F_RF_EXCL]*embU;
	  }
	}
	
	if (embU == 0 || beta*ener[F_EPOT] > bU_bin_limit) {
	  bin[0]++;
	} else {
	  i = (int)((bU_logV_bin_limit
		     - (beta*ener[F_EPOT] - logV + refvolshift))*invbinw
		    + 0.5);
	  if (i < 0)
	    i = 0;
	  if (i >= nbin)
	    realloc_bins(&bin,&nbin,i+10);
	  bin[i]++;
	}

	if (debug)
	  fprintf(debug,"TPI %7d %12.5e %12.5f %12.5f %12.5f\n",
		  step,ener[F_EPOT],x_tp[XX],x_tp[YY],x_tp[ZZ]);

	if (dump_pdb && ener[F_EPOT] <= dump_ener) {
	  sprintf(str,"t%g_step%d.pdb",t,step);
	  sprintf(str2,"t: %f step %d ener: %f",t,step,ener[F_EPOT]);
	  write_sto_conf(str,str2,&(top->atoms),state->x,state->v,
			 inputrec->ePBC,state->box);
	}
      }
    }

    if (PAR(cr)) {
      /* When running in parallel sum the energies over the processes */
      gmx_sumd(1,    &sum_embU, cr);
      gmx_sumd(nener,sum_UgembU,cr);
    }

    frame++;
    V_all += V;
    VembU_all += V*sum_embU/nsteps;

    if (fp_tpi) {
      if (bVerbose || frame%10==0 || frame<10)
	fprintf(stderr,"mu %10.3e <mu> %10.3e\n",
		-log(sum_embU/nsteps)/beta,-log(VembU_all/V_all)/beta);
      
      fprintf(fp_tpi,"%10.3f %12.5e %12.5e %12.5e %12.5e",
	      t,
	      VembU_all==0 ? 20/beta : -log(VembU_all/V_all)/beta,
	      sum_embU==0  ? 20/beta : -log(sum_embU/nsteps)/beta,
	      sum_embU/nsteps,V);
      for(e=0; e<nener; e++)
	fprintf(fp_tpi," %12.5e",sum_UgembU[e]/nsteps);
      fprintf(fp_tpi,"\n");
      fflush(fp_tpi);
    }

    bNotLastFrame = read_next_frame(status,&rerun_fr);
  } /* End of the loop  */

  update_time();

  close_trj(status);

  if (fp_tpi)
    fclose(fp_tpi);

  if (fplog) {
    fprintf(fplog,"\n");
    fprintf(fplog,"  <V>  = %12.5e nm^3\n",V_all/frame);
    fprintf(fplog,"  <mu> = %12.5e kJ/mol\n",-log(VembU_all/V_all)/beta);
  }
  
  /* Write the Boltzmann factor histogram */
  if (PAR(cr)) {
    /* When running in parallel sum the bins over the processes */
    i = nbin;
    global_max(cr,&i);
    realloc_bins(&bin,&nbin,i);
    gmx_sumd(nbin,bin,cr);
  }
  fp_tpi = xvgropen(opt2fn("-tpid",nfile,fnm),
		    "TPI energy distribution",
		    "\\8b\\4U - log(V/<V>)","count");
  sprintf(str,"number \\8b\\4U > %g: %9.3e",bU_bin_limit,bin[0]);
  xvgr_subtitle(fp_tpi,str);
  xvgr_legend(fp_tpi,2,tpid_leg);
  for(i=nbin-1; i>0; i--) {
    bUlogV = -i/invbinw + bU_logV_bin_limit - refvolshift + log(V_all/frame);
    fprintf(fp_tpi,"%6.2f %10d %12.5e\n",
	    bUlogV,
	    (int)(bin[i]+0.5),
	    bin[i]*exp(-bUlogV)*V_all/VembU_all);
  }
  fclose(fp_tpi);
  sfree(bin);

  sfree(sum_UgembU);

  return start_t;
}
