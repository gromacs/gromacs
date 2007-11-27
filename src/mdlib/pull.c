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


#include "futil.h"
#include "rdgroup.h"
#include "statutil.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "vec.h" 
#include "typedefs.h"
#include "network.h"
#include "filenm.h"
#include "string.h"
#include "smalloc.h"
#include "pull.h"
#include "xvgr.h"
#include "names.h"

static void pull_print_x(FILE *out,t_pull *pull, real t) 
{
  int g,m;
  
  fprintf(out, "%f", t);

  for (g=0; g<1+pull->ngrp; g++) {
    if (pull->grp[g].nat > 0) {
      for(m=0; m<DIM; m++)
	fprintf(out,"\t%f",pull->grp[g].x[m]);
    }
  }
  fprintf(out,"\n");
}

static void pull_print_f(FILE *out,t_pull *pull, real t) 
{
  int g,d;

  fprintf(out, "%f\t", t);

  for(g=1; g<1+pull->ngrp; g++) {
    if (pull->eGeom == epullgPOS) {
      for(d=0; d<DIM; d++)
	if (pull->dim[d])
	  fprintf(out,"\t%f",pull->grp[g].f[d]);
    } else {
      fprintf(out,"\t%f",pull->grp[g].f_scal);
    }
  }
  fprintf(out,"\n");
}

void pull_print_output(t_pull *pull, int step, real time)
{
  if ((pull->nstxout != 0) && (step % pull->nstxout == 0))
    pull_print_x(pull->out_x,pull,time);

  if ((pull->nstfout != 0) && (step % pull->nstfout == 0))
    pull_print_f(pull->out_f,pull,time);
}

static FILE *open_pull_out(char *fn,t_pull *pull,bool bCoord)
{
  FILE *fp;
  int  nsets,g,m;
  char **setname,buf[10];
  
  if (bCoord)
    fp = xvgropen(fn,"Pull COM",  "Time (ps)","Position (nm)");
  else
    fp = xvgropen(fn,"Pull force","Time (ps)","Force (kJ/mol)");
  
  snew(setname,(1+pull->ngrp)*DIM);
  nsets = 0;
  for(g=0; g<1+pull->ngrp; g++) {
    if (pull->grp[g].nat > 0 && (g > 0 || bCoord)) {
      if (bCoord || pull->eGeom == epullgPOS) {
	for(m=0; m<DIM; m++) {
	  sprintf(buf,"%d %c",g,'X'+m);
	  setname[nsets] = strdup(buf);
	  nsets++;
	}
      } else {
	sprintf(buf,"%d",g);
	setname[nsets] = strdup(buf);
	nsets++;
      }
    }
  }
  if (bCoord || nsets > 1)
    xvgr_legend(fp,nsets,setname);
  for(g=0; g<nsets; g++)
    sfree(setname[g]);
  sfree(setname);

  return fp;
}

/* Apply forces in a mass weighted fashion */
static void apply_forces_grp(t_commrec *cr,
			     t_pullgrp *pgrp, t_mdatoms * md, rvec * f,
			     dvec f_pull, int sign)
{
  int i,ii,m,start,end;
  double wmass,inv_wm;
  gmx_ga2la_t *ga2la=NULL;

  if (DOMAINDECOMP(cr))
    ga2la = cr->dd->ga2la;
  else
    ga2la = NULL;
  start = md->start;
  end   = md->homenr + start;

  inv_wm = pgrp->wscale*pgrp->invtm;
  
  for(i=0; i<pgrp->nat_loc; i++) {
    ii = pgrp->ind_loc[i];
    wmass = md->massT[ii];
    if (pgrp->weight_loc)
      wmass *= pgrp->weight_loc[i];
    
    for(m=0; m<DIM; m++)
      f[ii][m] += sign * wmass * f_pull[m] * inv_wm;
  }
}

/* Apply forces in a mass weighted fashion */
static void apply_forces(t_commrec *cr,
			 t_pull * pull, t_mdatoms * md, rvec * f)
{
  int i;
  t_pullgrp *pgrp;

  for(i=1; i<pull->ngrp+1; i++) {
    pgrp = &(pull->grp[i]);
    apply_forces_grp(cr,pgrp,md,f,pgrp->f,1);
    if (pull->grp[0].nat) {
      if (PULL_CYL(pull)) {
	apply_forces_grp(cr,&(pull->dyna[i]),md,f,pgrp->f,-1);
      } else {
	apply_forces_grp(cr,&(pull->grp[0]),md,f,pgrp->f,-1);
      }
    }
  }
}

void get_pullgrp_distance(t_pull *pull,int g,matrix box,double t,
			  dvec dr,dvec dev)
{
  static bool bWarned=FALSE;
  t_pullgrp *pref,*pgrp;
  int       ndim,m;
  dvec      ref;
  double    drs,inpr;
  
  pgrp = &pull->grp[g];
  
  if (PULL_CYL(pull))
    pref = &pull->dyna[g];
  else
    pref = &pull->grp[0];
  
  pull_d_pbc_dx(pull->npbcdim, box, pgrp->x, pref->x, dr);
  for(m=0; m<DIM; m++)
    dr[m] *= pull->dim[m];
  
  if (pull->eGeom == epullgPOS) {
    ndim = 3;
    for(m=0; m<ndim; m++)
      ref[m] = pgrp->init[m] + pgrp->rate*t*pgrp->vec[m];
  } else {
    ndim = 1;
    ref[0] = pgrp->init[0] + pgrp->rate*t;
  }
  
  switch (pull->eGeom) {
  case epullgDIST:
    /* Pull along the vector between the com's */
    if (ref[0] < 0 && !bWarned) {
      fprintf(stderr,"\nPull reference distance for group %d is negative (%f)\n",g,ref[0]);
      bWarned = TRUE;
    }
    drs = dnorm(dr);
    if (drs == 0) {
      /* With no vector we can not determine the direction for the force,
       * so we set the force to zero.
       */
      dev[0] = 0;
    } else {
      /* Determine the deviation */
      dev[0] = drs - ref[0];
    }
    break;
  case epullgDIR:
  case epullgCYL:
    /* Pull along vec */
    inpr = 0;
    for(m=0; m<DIM; m++)
      inpr += pgrp->vec[m]*dr[m];
    dev[0] = inpr - ref[0];
    break;
  case epullgPOS:
    /* Restrain to the location vec */
    pull_d_pbc_dx(pull->npbcdim, box, pgrp->x, ref, dev);
    for(m=0; m<DIM; m++)
      dev[m] *= pull->dim[m];
    break;
  }
}

/* Apply constraint using SHAKE */
static void do_constraint(t_pull *pull, rvec *x, rvec *v,
			  bool bVir, tensor vir,
			  matrix box, t_mdatoms *md, 
                          real dt, double t) 
{

  dvec *r_ij;  /* x[i] com of i in prev. step. Obeys constr. -> r_ij[i] */
  dvec unc_ij; /* xp[i] com of i this step, before constr.   -> unc_ij  */

  dvec *rinew;           /* current 'new' position of group i */
  dvec *rjnew;           /* current 'new' position of group j */
  dvec  ref,vec;
  double d0,inpr;
  double lambda, rm, mass, invdt=0;
  bool bConverged = FALSE;
  int niter=0,g,ii,j,m,max_iter=100;
  double q,a,b,c;  /* for solving the quadratic equation, 
		      see Num. Recipes in C ed 2 p. 184 */
  dvec *dr;        /* correction for group i */
  dvec ref_dr;     /* correction for group j */
  dvec tmp,tmp3;
  t_pullgrp *pdyna,*pgrp,*pref;

  snew(r_ij,pull->ngrp+1);
  if (PULL_CYL(pull)) {
    snew(rjnew,pull->ngrp+1);
  } else {
    snew(rjnew,1);
  }
  snew(dr,pull->ngrp+1);
  snew(rinew,pull->ngrp+1);

  /* copy the current unconstrained positions for use in iterations. We 
     iterate until rinew[i] and rjnew[j] obey the constraints. Then
     rinew - pull.x_unc[i] is the correction dr to group i */
  for(g=1; g<1+pull->ngrp; g++)
    copy_dvec(pull->grp[g].xp,rinew[g]);
  if (PULL_CYL(pull)) {
    for(g=1; g<1+pull->ngrp; g++)
      copy_dvec(pull->dyna[g].xp,rjnew[g]);
  } else {
    copy_dvec(pull->grp[0].xp,rjnew[0]);
  }

  /* Determine the constraint directions from the old positions */
  pref = &pull->grp[0];
  for(g=1; g<1+pull->ngrp; g++) {
    if (PULL_CYL(pull))
      pref = &pull->dyna[g];
    pull_d_pbc_dx(pull->npbcdim,box,pull->grp[g].x,pref->x,r_ij[g]);
    if (pull->eGeom == epullgDIR) {
      /* Select the component along vec */
      a = 0;
      for(m=0; m<DIM; m++)
	a += pull->grp[g].vec[m]*r_ij[g][m];
      for(m=0; m<DIM; m++)
	r_ij[g][m] = a*pull->grp[g].vec[m];
    } else {
      for(m=0; m<DIM; m++)
	r_ij[g][m] *= pull->dim[m];
    }
  }

  while (!bConverged && niter < max_iter) {
    /* loop over all constraints */
    for(g=1; g<1+pull->ngrp; g++) {
      pgrp = &pull->grp[g];
      if (PULL_CYL(pull))
	pref = &pull->dyna[g];
      else
	pref = &pull->grp[0];

      /* Get the current difference vector */
      pull_d_pbc_dx(pull->npbcdim,box,
		    rinew[g],rjnew[PULL_CYL(pull) ? g : 0],unc_ij);

      /* Select components we want */
      for(m=0; m<DIM; m++)
        unc_ij[m] *= pull->dim[m];

      if (pull->eGeom == epullgPOS) {
	for(m=0; m<DIM; m++)
	  ref[m] = pgrp->init[m] + pgrp->rate*t*pgrp->vec[m];
      } else {
	ref[0] = pgrp->init[0] + pgrp->rate*t;
      }

      if (debug)
        fprintf(debug,"Pull group %d, iteration %d\n",g,niter);

      rm = 1.0/(pull->grp[g].invtm + pref->invtm);

      switch (pull->eGeom) {
      case epullgDIST:
	if (ref[0] <= 0)
	  gmx_fatal(FARGS,"The pull constraint reference distance for group %d is <= 0 (%f)",g,ref);
	
	a = diprod(r_ij[g],r_ij[g]); 
	b = diprod(unc_ij,r_ij[g])*2;
	c = diprod(unc_ij,unc_ij) - sqr(ref[0]);

	if (b < 0) {
	  q = -0.5*(b - sqrt(b*b - 4*a*c));
	  lambda = -q/a;
	} else {
	  q = -0.5*(b + sqrt(b*b - 4*a*c));
	  lambda = -c/q;
	}

	if (debug)
	  fprintf(debug,"Pull ax^2+bx+c=0: a=%e b=%e c=%e lambda=%e\n",
		  a,b,c,lambda);
	
	/* The position corrections dr due to the constraints */
	dsvmul(-lambda*rm*pgrp->invtm, r_ij[g],  dr[g]);
	dsvmul( lambda*rm*pref->invtm, r_ij[g], ref_dr);
	break;
      case epullgDIR:
      case epullgCYL:
	/* A 1-dimensional constraint along a vector */
	a = 0;
	for(m=0; m<DIM; m++) {
	  vec[m] = pgrp->vec[m];
	  a += unc_ij[m]*vec[m];
	}
	/* Select only the component along the vector */
        dsvmul(a,vec,unc_ij);
	lambda = a - ref[0];
	if (debug)
	  fprintf(debug,"Pull inpr %e lambda: %e\n",a,lambda);

	/* The position corrections dr due to the constraints */
	dsvmul(-lambda*rm*pull->grp[g].invtm, vec, dr[g]);
	dsvmul( lambda*rm*       pref->invtm, vec,ref_dr);
	break;
      case epullgPOS:
	for(m=0; m<DIM; m++) {
	  if (pull->dim[m]) {
	    lambda = r_ij[g][m] - pgrp->vec[m];
	    /* The position corrections dr due to the constraints */
	    dr[g][m]  = -lambda*rm*pull->grp[g].invtm;
	    ref_dr[m] =  lambda*rm*pref->invtm;
	  } else {
	    dr[g][m]  = 0;
	    ref_dr[m] = 0;
	  }
	}
	break;
      }

      /* DEBUG */
      if (debug) {
        j = (PULL_CYL(pull) ? g : 0);
	pull_d_pbc_dx(pull->npbcdim, box, rinew[g], rjnew[j], tmp);
        pull_d_pbc_dx(pull->npbcdim, box,    dr[g],   ref_dr, tmp3);
        for(m=0; m<DIM; m++) {
          tmp[m]  *= pull->dim[m];
          tmp3[m] *= pull->dim[m];
        }
	fprintf(debug,
		"Pull cur %8.5f %8.5f %8.5f j:%8.5f %8.5f %8.5f d: %8.5f\n",
		rinew[g][0],rinew[g][1],rinew[g][2], 
		rjnew[j][0],rjnew[j][1],rjnew[j][2], dnorm(tmp));
	if (pull->eGeom == epullgPOS)
	  fprintf(debug,
		  "Pull ref %8.5f %8.5f %8.5f\n",
		  pgrp->vec[0],pgrp->vec[1],pgrp->vec[2]);
	else
	  fprintf(debug,
		  "Pull ref %8s %8s %8s   %8s %8s %8s d: %8.5f %8.5f %8.5f\n",
		  "","","","","","",ref[0],ref[1],ref[2]);
	fprintf(debug,
		"Pull cor %8.5f %8.5f %8.5f j:%8.5f %8.5f %8.5f d: %8.5f\n",
		dr[g][0],dr[g][1],dr[g][2],
		ref_dr[0],ref_dr[1],ref_dr[2],
		dnorm(tmp3));
	fprintf(debug,
		"Pull cor %10.7f %10.7f %10.7f\n",
		dr[g][0],dr[g][1],dr[g][2]);
      } /* END DEBUG */

      /* Update the COMs with dr */
      dvec_inc(rinew[g],                     dr[g]);
      dvec_inc(rjnew[PULL_CYL(pull) ? g : 0],ref_dr);
    }

    /* Check if all constraints are fullfilled now */
    for(g=1; g<1+pull->ngrp; g++) {
      pgrp = &pull->grp[g];

      pull_d_pbc_dx(pull->npbcdim, box,
		    rinew[g], rjnew[PULL_CYL(pull) ? g : 0],unc_ij);
      for(m=0; m<DIM; m++)
        unc_ij[m] *= pull->dim[m];

      switch (pull->eGeom) {
      case epullgDIST:
	bConverged = fabs(dnorm(unc_ij) - ref[0]) < pull->constr_tol;
	break;
      case epullgDIR:
      case epullgCYL:
	for(m=0; m<DIM; m++)
	  vec[m] = pgrp->vec[m];
	inpr = diprod(unc_ij,vec);
	dsvmul(inpr,vec,unc_ij);
	bConverged =
	  fabs(diprod(unc_ij,vec) - ref[0]) < pull->constr_tol;
	break;
      case epullgPOS:
	bConverged = TRUE;
	for(m=0; m<DIM; m++) {
	  if (pull->dim[m] && 
	      fabs(unc_ij[m] - pgrp->vec[m]) >= pull->constr_tol)
	    bConverged = FALSE;
	}
	break;
      }

      /* DEBUG */
      if (debug) {
	if (!bConverged)
	  fprintf(debug,"NOT CONVERGED YET: Group %d:"
		  "d_ref = %f %f %f, current d = %f\n",
		  g,ref[0],ref[1],ref[2],dnorm(unc_ij));
      } /* END DEBUG */
    }

    niter++;
    /* if after all constraints are dealt with and bConverged is still TRUE
       we're finished, if not we do another iteration */
  }
  if (niter > max_iter)
    gmx_fatal(FARGS,"Too many iterations for constraint run: %d",niter);

  /* DONE ITERATING, NOW UPDATE COORDINATES AND CALC. CONSTRAINT FORCES */

  if (v)
    invdt = 1/dt;

  /* update the normal groups */
  for(g=1; g<1+pull->ngrp; g++) {
    pgrp = &pull->grp[g];
    /* get the final dr and constraint force for group i */
    dvec_sub(rinew[g],pgrp->xp,dr[g]);
    /* select components of dr */
    for(m=0; m<DIM; m++)
      dr[g][m] *= pull->dim[m];
    dsvmul(1.0/(pgrp->invtm*dt*dt),dr[g],pgrp->f);
    switch (pull->eGeom) {
    case epullgDIST:
      pgrp->f_scal = 0;
      for(m=0; m<DIM; m++)
	pgrp->f_scal += r_ij[g][m]*pgrp->f[m];
      pgrp->f_scal /= dnorm(r_ij[g]);
      break;
    case epullgDIR:
    case epullgCYL:
      pgrp->f_scal = 0;
      for(m=0; m<DIM; m++)
	pgrp->f_scal += pgrp->vec[m]*pgrp->f[m];
      break;
    case epullgPOS:
      break;
    }

    if (bVir) {
      /* Add the pull contribution to the virial */
      for(j=0; j<DIM; j++)
	for(m=0; m<DIM; m++)
	  vir[j][m] += 0.5*pull->grp[g].f[j]*r_ij[g][m];
    }

    /* update the atom positions */
    copy_dvec(dr[g],tmp);
    for(j=0;j<pgrp->nat_loc;j++) {
      ii = pgrp->ind_loc[j];
      if (pgrp->weight_loc)
	dsvmul(pgrp->wscale*pgrp->weight_loc[j],dr[g],tmp); 
      for(m=0; m<DIM; m++)
	x[ii][m] += tmp[m];
      if (v)
	for(m=0; m<DIM; m++)
	  v[ii][m] += invdt*tmp[m];
    }
  }

  /* update the reference groups */
  if (PULL_CYL(pull)) {
    /* update the dynamic reference groups */
    for(g=1; g<1+pull->ngrp; g++) {
      pdyna = &pull->dyna[g];
      dvec_sub(rjnew[g],pdyna->xp,ref_dr);
      /* select components of ref_dr */
      for(m=0; m<DIM; m++)
        ref_dr[m] *= pull->dim[m];

      for(j=0;j<pdyna->nat_loc;j++) {
        /* reset the atoms with dr, weighted by w_i */
        dsvmul(pdyna->wscale*pdyna->weight_loc[j],ref_dr,tmp); 
        ii = pdyna->ind_loc[j];
	for(m=0; m<DIM; m++)
	  x[ii][m] += tmp[m];
	if (v)
	  for(m=0; m<DIM; m++)
	    v[ii][m] += invdt*tmp[m];
      }
    }
  } else {
    pgrp = &pull->grp[0];
    /* update the reference group */
    dvec_sub(rjnew[0],pgrp->xp, ref_dr); 
    /* select components of ref_dr */
    for(m=0;m<DIM;m++)
      ref_dr[m] *= pull->dim[m];
    
    copy_dvec(ref_dr,tmp);
    for(j=0; j<pgrp->nat_loc;j++) {
      ii = pgrp->ind_loc[j];
      if (pgrp->weight_loc)
	dsvmul(pgrp->wscale*pgrp->weight_loc[j],ref_dr,tmp); 
      for(m=0; m<DIM; m++)
	x[ii][m] += tmp[m];
      if (v)
	for(m=0; m<DIM; m++)
	  v[ii][m] += invdt*tmp[m];
    }
  }

  /* finished! I hope. Give back some memory */
  sfree(r_ij);
  sfree(rinew);
  sfree(rjnew);
  sfree(dr);
}

/* Pulling with a harmonic umbrella potential or constant force */
static real do_pull_pot(t_commrec *cr,int ePull,
			t_pull *pull, rvec *f, bool bVir, tensor vir, 
			matrix box, t_mdatoms *md,
			double t)
{
  int       g,j,m;
  dvec      dr,dev;
  double    ndr,invdr;
  real      V;
  t_pullgrp *pgrp;

  /* loop over the groups that are being pulled */
  V = 0;
  for(g=1; g<1+pull->ngrp; g++) {
    pgrp = &pull->grp[g];
    get_pullgrp_distance(pull,g,box,t,dr,dev);

    switch (pull->eGeom) {
    case epullgDIST:
      ndr   = dnorm(dr);
      invdr = 1/ndr;
      if (ePull == epullUMBRELLA) {
	pgrp->f_scal  = -pgrp->k*dev[0];
	V            += 0.5*pgrp->k*sqr(dev[0]);
      } else {
	pgrp->f_scal  = -pgrp->k;
	V            += pgrp->k*ndr;
      }
      for(m=0; m<DIM; m++)
	pgrp->f[m]    = pgrp->f_scal*dr[m]*invdr;
      break;
    case epullgDIR:
    case epullgCYL:
      if (ePull == epullUMBRELLA) {
	pgrp->f_scal  = -pgrp->k*dev[0];
	V            += 0.5*pgrp->k*sqr(dev[0]);
      } else {
	ndr = 0;
	for(m=0; m<DIM; m++)
	  ndr += pgrp->vec[m]*dr[m];
	pgrp->f_scal  = -pgrp->k;
	V            += pgrp->k*ndr;
      }
      for(m=0; m<DIM; m++)
	pgrp->f[m]    = pgrp->f_scal*pgrp->vec[m];
      break;
    case epullgPOS:
      for(m=0; m<DIM; m++) {
	if (ePull == epullUMBRELLA) {
	  pgrp->f[m]  = -pgrp->k*dev[m];
	  V          += 0.5*pgrp->k*sqr(dev[m]);
	} else {
      	  pgrp->f[m]  = -pgrp->k*pull->dim[m];
	  V          += pgrp->k*dr[m]*pull->dim[m];
	}
      }
      break;
    }
    
    if (bVir) {
      /* Add the pull contribution to the virial */
      for(j=0; j<DIM; j++)
	for(m=0;m<DIM;m++)
	  vir[j][m] += 0.5*pgrp->f[j]*dr[m];
    }
  }

  /* Distribute forces over pulled groups */
  apply_forces(cr, pull, md, f);

  return V;
}

real pull_potential(int ePull,t_pull *pull, rvec *x, rvec *f, tensor vir, 
		    matrix box, t_topology *top, double t,
		    t_mdatoms *md, t_commrec *cr) 
{
  pull_calc_coms(cr,pull,md,x,NULL,box);

  return do_pull_pot(cr,ePull,pull,f,MASTER(cr),vir,box,md,t);
}

void pull_constraint(t_pull *pull, rvec *x, rvec *xp, rvec *v, tensor vir, 
		     matrix box, t_topology *top, real dt, double t,
		     t_mdatoms *md, t_commrec *cr) 
{
  pull_calc_coms(cr,pull,md,x,xp,box);

  do_constraint(pull,xp,v,MASTER(cr),vir,box,md,dt,t);
}

static void dd_make_local_pull_group(gmx_domdec_t *dd,
				     t_pullgrp *pg,t_mdatoms *md)
{
  int i,ii;
  gmx_ga2la_t *ga2la=NULL;

  ga2la = dd->ga2la;
  pg->nat_loc = 0;
  for(i=0; i<pg->nat; i++) {
    if (ga2la[pg->ind[i]].cell == 0) {
      ii = ga2la[pg->ind[i]].a;
      if (ii < md->start+md->homenr) {
	/* This is a home atom, add it to the local pull group */
	pg->ind_loc[pg->nat_loc] = ii;
	if (pg->weight)
	  pg->weight_loc[pg->nat_loc] = pg->weight[i];
	pg->nat_loc++;
      }
    }
  }
}

void dd_make_local_pull_groups(gmx_domdec_t *dd,t_pull *pull,t_mdatoms *md)
{
  int g;

  if (pull->eGeom != epullgPOS)
    dd_make_local_pull_group(dd,&pull->grp[0],md);
  for(g=1; g<1+pull->ngrp; g++)
    dd_make_local_pull_group(dd,&pull->grp[g],md);
}

static void init_pull_group_index(FILE *log,t_commrec *cr,
				  int start,int end,
				  int g,t_pullgrp *pg,ivec pulldims,
				  t_atoms *atoms,ivec nFreeze[])
{
  int i,ii,d,nfreeze,ndim;
  real m,w;
  double wmass,wwmass,buf[3];
  bool bDomDec;
  gmx_ga2la_t *ga2la=NULL;

  bDomDec = DOMAINDECOMP(cr);
  if (bDomDec)
    ga2la = cr->dd->ga2la;

  if (PAR(cr)) {
    snew(pg->ind_loc,pg->nat);
    pg->nat_loc = 0;
    if (pg->weight)
      snew(pg->weight_loc,pg->nat);
  } else {
    pg->nat_loc = pg->nat;
    pg->ind_loc = pg->ind;
    pg->weight_loc = pg->weight;
  }

  nfreeze = 0;
  wmass = 0;
  wwmass = 0;
  for(i=0; i<pg->nat; i++) {
    ii = pg->ind[i];
    if (PAR(cr) && !bDomDec && ii >= start && ii < end)
      pg->ind_loc[pg->nat_loc++] = ii;
    if (nFreeze) {
      for(d=0; d<DIM; d++)
	if (pulldims[d] && nFreeze[atoms->atom[ii].grpnr[egcFREEZE]][d])
	  nfreeze++;
    }
    m = atoms->atom[ii].m;
    if (pg->weight)
      w = pg->weight[i];
    else
      w = 1;
    wmass += w*m;
    wwmass += w*w*m;
  }

  if (nfreeze == 0) {
    pg->wscale = wmass/wwmass;
    pg->invtm  = 1.0/(pg->wscale*wmass);
  } else {
    ndim = 0;
    for(d=0; d<DIM; d++)
      ndim += pulldims[d]*pg->nat;
    if (nfreeze > 0 && nfreeze < ndim)
      fprintf(log,
	      "\nWARNING: In pull group %d some, but not all of the degrees of freedom\n"
	      "         that are subject to pulling are frozen.\n"
	      "         For pulling the whole group will be frozen.\n\n",
	      g);
    pg->wscale = 1.0;
    pg->invtm  = 0.0;
  }
}

void init_pull(FILE *fplog,t_inputrec *ir,int nfile,t_filenm fnm[],
	       rvec *x,t_atoms *atoms,matrix box,
	       t_commrec *cr,int start,int end)
{
  t_pull    *pull;
  t_pullgrp *pgrp;
  int       g;

  pull = ir->pull;

  pull->ePBC = ir->ePBC;
  switch (pull->ePBC) {
  case epbcNONE: pull->npbcdim = 0; break;
  case epbcXY:   pull->npbcdim = 2; break;
  default:       pull->npbcdim = 3; break;
  }

  if (fplog) {
    fprintf(fplog,"\nWill apply %s COM pulling in geometry '%s'\n",
	    EPULLTYPE(ir->ePull),EPULLGEOM(pull->eGeom));
    if (pull->grp[0].nat > 0)
      fprintf(fplog,"between a reference group and %d group%s\n",
	      pull->ngrp,pull->ngrp==1 ? "" : "s");
    else
      fprintf(fplog,"with an absolute reference on %d group%s\n",
	      pull->ngrp,pull->ngrp==1 ? "" : "s");
  }

  for(g=0; g<pull->ngrp+1; g++) {
    pgrp = &pull->grp[g];
    if (pgrp->nat > 0) {
      /* Set the indices */
      init_pull_group_index(fplog,cr,start,end,
			    g,pgrp,pull->dim,atoms,ir->opts.nFreeze);
    }
  }      
  
  /* if we use dynamic reference groups, do some initialising for them */
  if (PULL_CYL(pull)) {
    if (pull->grp[0].nat == 0)
      gmx_fatal(FARGS, "Dynamic reference groups are not support when using absolute reference!\n");
    
    snew(pull->dyna,pull->ngrp+1);
    for(g=1; g<1+pull->ngrp; g++) {
      snew(pull->dyna[g].ind,pull->grp[0].nat);
      snew(pull->dyna[g].weight_loc,pull->grp[0].nat);
    }
  }
  
  /* Only do I/O if we are the MASTER */
  if (MASTER(cr)) {
    if (pull->nstxout > 0) {
      pull->out_x = open_pull_out(opt2fn("-px",nfile,fnm),pull,TRUE);
    }
    if (pull->nstfout > 0) {
      pull->out_f = open_pull_out(opt2fn("-pf",nfile,fnm),pull,FALSE);
    }
  }
}
