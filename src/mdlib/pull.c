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
#include "pull_internal.h"

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
    if (pgrp->weight_loc > 0)
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

/* Apply constraint using SHAKE */
static void do_constraint(t_pull *pull, rvec *x, rvec *v,
			  bool bVir, tensor vir,
			  matrix box, t_mdatoms *md, 
                          real dt, int step) 
{

  dvec *r_ij;  /* x[i] com of i in prev. step. Obeys constr. -> r_ij[i] */
  dvec unc_ij; /* xp[i] com of i this step, before constr.   -> unc_ij  */

  dvec *rinew;           /* current 'new' position of group i */
  dvec *rjnew;           /* current 'new' position of group j */
  double d0,ref,inpr;
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

  ref = 0;
  
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
      a = diprod(pull->grp[g].vec,r_ij[g]);
      dsvmul(a,pull->grp[g].vec,r_ij[g]);
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

      ref = pgrp->init + pgrp->rate*step*dt;

      if (debug)
        fprintf(debug,"Pull group %d, iteration %d\n",g,niter);

      rm = 1.0/(pull->grp[g].invtm + pref->invtm);

      switch (pull->eGeom) {
      case epullgDIST:
	if (ref <= 0)
	  gmx_fatal(FARGS,"The pull constraint reference distance for group %d is <= 0 (%f)",g,ref);
	
	a = diprod(r_ij[g],r_ij[g]); 
	b = diprod(unc_ij,r_ij[g])*2;
	c = diprod(unc_ij,unc_ij) - sqr(ref);

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
	for(m=0; m<DIM; m++)
	  a += unc_ij[m]*pgrp->vec[m];
	/* Select only the component along the vector */
        dsvmul(a,pgrp->vec,unc_ij);
	lambda = a - ref;
	if (debug)
	  fprintf(debug,"Pull inpr %e lambda: %e\n",a,lambda);

	/* The position corrections dr due to the constraints */
	dsvmul(-lambda*rm*pull->grp[g].invtm, pgrp->vec, dr[g]);
	dsvmul( lambda*rm*       pref->invtm, pgrp->vec,ref_dr);
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
		  "Pull ref %8s %8s %8s   %8s %8s %8s d: %8.5f\n",
		  "","","","","","",ref);
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
	bConverged = fabs(dnorm(unc_ij) - ref) < pull->constr_tol;
	break;
      case epullgDIR:
      case epullgCYL:
	inpr = diprod(unc_ij,pgrp->vec);
	dsvmul(inpr,pgrp->vec,unc_ij);
	bConverged = fabs(diprod(unc_ij,pgrp->vec) - ref) < pull->constr_tol;
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
	  fprintf(debug,"NOT CONVERGED YET: Group %d (%s):"
		  "d_ref = %f, current d = %f\n",
		  g,pgrp->name,ref,dnorm(unc_ij));
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

/* Pulling with a harmonic umbrella potential */
static void do_umbrella(t_commrec *cr,
			t_pull *pull, rvec *f, bool bVir, tensor vir, 
			matrix box, t_mdatoms *md,
			real dt, int step)
{
  static bool bWarned=FALSE;
  int g,j,m;
  dvec dr,ds;     /* extension of the springs */
  double drs,dss,inpr,ref;
  t_pullgrp *pgrp,*pref;

  /* loop over the groups that are being pulled */
  for(g=1; g<1+pull->ngrp; g++) {
    pgrp = &pull->grp[g];
    if (PULL_CYL(pull))
      pref = &pull->dyna[g];
    else
      pref = &pull->grp[0];

    pull_d_pbc_dx(pull->npbcdim, box, pgrp->x, pref->x, dr);
    for(m=0; m<DIM; m++)
      dr[m] *= pull->dim[m];

    switch (pull->eGeom) {
    case epullgDIST:
      ref = pgrp->init + pgrp->rate*step*dt;
      /* Pull along the vector between the com's */
      drs = dnorm(dr);
      if (drs == 0) {
	dss = 0;
      } else {
	/* Determine the extension */
	dss = drs - ref;
	if (ref < 0) {
	  if (!bWarned) {
	    fprintf(stdlog,"\nPull reference distance for group %d is negative (%f) will use zero\n",g,ref);
	    fprintf(stderr,"\nPull reference distance for group %d is negative (%f) will use zero\n",g,ref);
	    bWarned = TRUE;
	  }
	  dss = 0;
	}
      }
      pgrp->f_scal = -pgrp->k*dss;
      for(m=0; m<DIM; m++)
	pull->grp[g].f[m] = pgrp->f_scal*dr[m]/drs;
      break;
    case epullgDIR:
    case epullgCYL:
      ref = pgrp->init + pgrp->rate*step*dt;
      /* Pull along vec */
      drs = dnorm(dr);
      inpr = 0;
      for(m=0; m<DIM; m++)
	inpr += pgrp->vec[m]*dr[m];
      dss = inpr - ref;
      pgrp->f_scal = -pgrp->k*dss;
      for(m=0; m<DIM; m++)
	pgrp->f[m] = pgrp->f_scal*dr[m]/drs;
      break;
    case epullgPOS:
      /* Restrain to the location vec */
      pull_d_pbc_dx(pull->npbcdim, box, pgrp->x, pgrp->vec, ds);
      for(m=0; m<DIM; m++)
	pgrp->f[m] = -pgrp->k*pull->dim[m]*ds[m];
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
}

void pull_umbrella(t_pull *pull, rvec *x, rvec *f, tensor vir, 
		   matrix box, t_topology *top, real dt, int step,
		   t_mdatoms *md, t_commrec *cr) 
{
  pull_calc_coms(cr,pull,md,x,NULL,box);

  do_umbrella(cr,pull,f,MASTER(cr),vir,box,md,dt,step);
}

void pull_constraint(t_pull *pull, rvec *x, rvec *xp, rvec *v, tensor vir, 
		     matrix box, t_topology *top, real dt, int step,
		     t_mdatoms *md, t_commrec *cr) 
{
  pull_calc_coms(cr,pull,md,x,xp,box);

  do_constraint(pull,xp,v,MASTER(cr),vir,box,md,dt,step);
}
