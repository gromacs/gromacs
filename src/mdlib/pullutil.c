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


#include <stdlib.h>
#include "sysstuff.h"
#include "princ.h"
#include "futil.h"
#include "statutil.h"
#include "vec.h"
#include "smalloc.h"
#include "typedefs.h"
#include "names.h"
#include "gmx_fatal.h"
#include "macros.h"
#include "rdgroup.h"
#include "symtab.h"
#include "index.h"
#include "confio.h"
#include "network.h"
#include "pbc.h"
#include "pull.h"
#include "pull_internal.h"

void pull_d_pbc_dx(int npbcdim,matrix box,
		   const dvec x1, const dvec x2, dvec dx)
{
  int m,d;
  
  /* Only correct for atom pairs with a distance within
   * half of the smallest diagonal element of box.
   */
  dvec_sub(x1,x2,dx);
  for(m=npbcdim-1; m>=0; m--) {
    while (dx[m] < -0.5*box[m][m]) {
      for(d=0; d<DIM; d++)
	dx[d] += box[m][d];
    }
    while (dx[m] >=  0.5*box[m][m]) {
      for(d=0; d<DIM; d++)
	dx[d] -= box[m][d];
    }
  }
}

static void pull_pbc_dx(int npbcdim,matrix box,
			const rvec x1, const rvec x2, rvec dx)
{
  int m,d;
  
  /* Only correct for atom pairs with a distance within
   * half of the smallest diagonal element of box.
   */
  rvec_sub(x1,x2,dx);
  for(m=npbcdim-1; m>=0; m--) {
    while (dx[m] < -0.5*box[m][m]) {
      for(d=0; d<DIM; d++)
	dx[d] += box[m][d];
    }
    while (dx[m] >= 0.5*box[m][m]) {
      for(d=0; d<DIM; d++)
	dx[d] -= box[m][d];
    }
  }
}

static void pull_set_pbcatom(t_commrec *cr, t_pullgrp *pg,
			     t_mdatoms *md, rvec *x,
			     rvec x_pbc)
{
  int a,m;

  if (DOMAINDECOMP(cr)) {
    if (cr->dd->ga2la[pg->pbcatom].cell == 0)
      a = cr->dd->ga2la[pg->pbcatom].a;
    else
      a = -1;
  } else {
    a = pg->pbcatom;
  }

  if (a >= md->start && a < md->start+md->homenr) {
    copy_rvec(x[a],x_pbc);
  } else {
   clear_rvec(x_pbc);
  }
}

static void pull_set_pbcatoms(t_commrec *cr, t_pull *pull,
			      t_mdatoms *md, rvec *x,
			      rvec *x_pbc)
{
  int g,n;
  
  n = 0;
  for(g=0; g<1+pull->ngrp; g++) {
    if ((g==0 && PULL_CYL(pull)) || pull->grp[g].pbcatom == -1) {
      clear_rvec(x_pbc[g]);
    } else {
      pull_set_pbcatom(cr,&pull->grp[g],md,x,x_pbc[g]);
      n++;
    }
  }
  
  if (PAR(cr) && n > 0) {
    /* Sum over the nodes to get x_pbc from the home node of pbcatom */
    gmx_sumf((1+pull->ngrp)*DIM,x_pbc[0],cr);
  }
}

/* switch function, x between r and w */
static real get_weight(real x, real r1, real r0)
{
  real weight; 

  if (x >= r0)
    weight = 0;
  else if (x <= r1)
    weight = 1;
  else
    weight = (r0 - x)/(r0 - r1);

  return weight;
}

static double get_cylinder_distance(rvec x, dvec com, matrix box) {
  double dr, dx, dy, box_xx, box_yy, box_yx;

  box_xx = box[XX][XX];
  box_yy = box[YY][YY];
  box_yx = box[YY][YY];
  
  dx = x[XX] - com[XX];
  dy = x[YY] - com[YY];
  
  while (dy > 0.5*box_yy) {
    dy -= box_yy;
    dx -= box_yx;
  }
  while (dy < -0.5*box_yy) {
    dy += box_yy;
    dx += box_yx;
  }
  while (dx > 0.5*box_xx)
    dx -= box_xx;
  while (dx < 0.5*box_xx)
    dx += box_xx;

  dr = sqrt(dx*dx+dy*dy);
#ifdef CYLDEBUG
  fprintf(stderr,"x:%8.3f%8.3f%8.3f com:%8.3f%8.3f%8.3f dx,dy,dr:%8.3f%8.3f%8.3f\n",x[0],x[1],x[2],com[0],com[1],com[2],dx,dy,dr);
#endif
  return dr;
}

static void make_cyl_refgrps(t_commrec *cr,t_pull *pull,t_mdatoms *md,
			     rvec *x,matrix box) 
{
  int g,i,ii,m,start,end;

  double dr,mass,weight,wmass,wwmass;
  dvec test;
  t_pullgrp *pref,*pdyna;
  gmx_ga2la_t *ga2la=NULL;

  if (DOMAINDECOMP(cr))
    ga2la = cr->dd->ga2la;

  start = md->start;
  end   = md->homenr + start;

  /* loop over all groups to make a reference group for each*/
  pref = &pull->grp[0];
  for(g=1; g<1+pull->ngrp+1; g++) {
    pdyna = &pull->dyna[g];
    wmass = 0;
    wwmass = 0;
    pdyna->nat_loc = 0;

    /* loop over all atoms in the main ref group */
    for(i=0; i<pref->nat; i++) {
      ii = pull->grp[0].ind[i];
      if (ga2la) {
	if (ga2la[pref->ind[i]].cell == 0)
	  ii = ga2la[pref->ind[i]].a;
	else
	  ii = -1;
      }
      if (ii >= start && ii < end) {
	/* get_distance takes pbc into account */
	dr = get_cylinder_distance(x[ii],pull->grp[g].x,box);

	if (dr < pull->cyl_r0) {
	  /* add to index, to sum of COM, to weight array */
	  pdyna->ind_loc[pdyna->nat_loc] = ii;
	  mass = md->massT[ii];
	  weight = get_weight(dr,pull->cyl_r1,pull->cyl_r0);
	  pdyna->weight_loc[pdyna->nat_loc] = weight;
	  pdyna->x[ZZ] += mass*weight*x[ii][ZZ];
	  wmass += mass*weight;
	  wwmass += mass*sqr(weight);
	  pdyna->nat_loc++;
	}
      }
    }

    pdyna->wscale = wmass/wwmass;
    pdyna->invtm = 1.0/(pdyna->wscale*wmass);

    /* normalize the new 'x_unc' */
    dsvmul(1/wmass,pdyna->x,pdyna->x);
    if (debug)
      fprintf(debug,"Pull cylinder group %d:%8.3f%8.3f%8.3f m:%8.3f\n",
              g,pdyna->x[0],pdyna->x[1],
              pdyna->x[2],1.0/pdyna->invtm);
  }
}

/* calculates center of mass of selection index from all coordinates x */
void pull_calc_coms(t_commrec *cr,
		    t_pull *pull, t_mdatoms *md, 
		    rvec x[], rvec *xp, matrix box)
{
  static rvec *rbuf=NULL;
  static dvec *dbuf=NULL;
  int  g,i,ii,m,npbcdim;
  real wm;
  dvec com,comp;
  rvec *xx[2],x_pbc={0,0,0},dx;
  t_pullgrp *pgrp;

  if (rbuf == NULL)
    snew(rbuf,1+pull->ngrp);
  if (dbuf == NULL)
    snew(dbuf,2*(1+pull->ngrp));

  npbcdim = pull->npbcdim;

  pull_set_pbcatoms(cr,pull,md,x,rbuf);
  
  for (g=0; g<1+pull->ngrp; g++) {
    pgrp = &pull->grp[g];
    clear_dvec(com);
    clear_dvec(comp);
    if (!(g==0 && PULL_CYL(pull))) {
      if (pgrp->pbcatom >= 0) {
	/* Set the pbc atom */
	copy_rvec(rbuf[g],x_pbc);
      }
      for(i=0; i<pgrp->nat_loc; i++) {
	ii = pgrp->ind_loc[i];
	wm = md->massT[ii];
	if (pgrp->weight_loc)
	  wm *= pgrp->weight_loc[i];
	if (pgrp->pbcatom == -1) {
	  /* Sum the coordinates */
	  for(m=0; m<DIM; m++)
	    com[m] += wm*x[ii][m];
	} else {
	  /* Sum the difference with the reference atom */
	  pull_pbc_dx(npbcdim,box,x[ii],x_pbc,dx);
	  for(m=0; m<DIM; m++)
	    com[m] += wm*dx[m];
	}
	if (xp) {
	  if (pgrp->pbcatom == -1) {
	    /* Sum the coordinates */
	    for(m=0; m<DIM; m++)
	      comp[m] += wm*xp[ii][m];
	  } else {
	    /* Sum the difference with the reference atom */
	    pull_pbc_dx(npbcdim,box,xp[ii],x_pbc,dx);
	    for(m=0; m<DIM; m++)
	      comp[m] += wm*dx[m];
	  }
	}
      }
    }
    copy_dvec(com,dbuf[g]);
    if (xp)
      copy_dvec(comp,dbuf[1+pull->ngrp+g]);
  }

  if (PAR(cr)) {
    /* Sum the contributions over the nodes */
    gmx_sumd((xp ? 2 : 1)*(1+pull->ngrp)*DIM,dbuf[0],cr);
  }
  
  for (g=0; g<1+pull->ngrp; g++) {
    pgrp = &pull->grp[g];
    /* Divide by the total mass */
    for(m=0; m<DIM; m++)
      pgrp->x[m] = dbuf[g][m]*pgrp->wscale*pgrp->invtm;
    if (pgrp->pbcatom >= 0) {
      /* Add the location of the reference atom */
      for(m=0; m<DIM; m++)
	pgrp->x[m] += rbuf[g][m];
    }
    if (xp) {
      /* Divide by the total mass */
      for(m=0; m<DIM; m++)
	pgrp->xp[m] = dbuf[1+pull->ngrp+g][m]*pgrp->wscale*pgrp->invtm;
      if (pgrp->pbcatom >= 0) {
	/* Add the location of the reference atom */
	for(m=0; m<DIM; m++)
	  pgrp->xp[m] += rbuf[g][m];
      }
    }
  }

  if (PULL_CYL(pull)) {
    /* Calculate the COMs for the cyclinder reference groups */
    make_cyl_refgrps(cr,pull,md,x,box);
  }  
}
