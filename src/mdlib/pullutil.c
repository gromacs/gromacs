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
#include "fatal.h"
#include "macros.h"
#include "rdgroup.h"
#include "symtab.h"
#include "index.h"
#include "confio.h"
#include "pbc.h"
#include "pull.h"
#include "pull_internal.h"

void d_pbc_dx(matrix box,const dvec x1, const dvec x2, dvec dx)
{
  int m,d;
  
  /* Not correct for all trilinic boxes !!!
     Should be as pbc_dx and moved to pbc.c
  */
  dvec_sub(x1,x2,dx);
  for(m=DIM-1; m>=0; m--) {
    while (dx[m]   < -0.5*box[m][m])
      for(d=0; d<DIM; d++)
	dx[d] += box[m][d];
    while (dx[m]   >=  0.5*box[m][m]) 
      for(d=0; d<DIM; d++)
	dx[d] -= box[m][d];
  }
}

void put_dvec_in_box(matrix box, dvec v)
{
  int m,d;
  
  for(m=DIM-1; m>=0; m--) {
    while (v[m] < 0)
      for(d=0; d<DIM; d++)
	v[d] += box[m][d];
    while (v[m] >= box[m][m])
      for(d=0; d<DIM; d++)
	v[d] -= box[m][d];
  }
}

/* calculates center of mass of selection index from all coordinates x */
void calc_com(t_pullgrp *pg, rvec x[], t_mdatoms *md, matrix box)
{
  int  i,ii,m;
  real wm;
  dvec com;

  clear_dvec(com);
  for(i=0; i<pg->ngx; i++) {
    ii = pg->idx[i];
    wm = md->massT[ii];
    if (pg->nweight > 0)
      wm *= pg->weight[i];
    for(m=0; m<DIM; m++)
      com[m] += wm*x[ii][m];
  }
  for(m=0; m<DIM; m++)
    com[m] *= pg->wscale*pg->invtm;
  put_dvec_in_box(box,com);
  copy_dvec(com,pg->x_unc);
}

/* calculates com of all atoms in x[], *index has their index numbers
   to get the masses from atom[] */
void calc_com2(t_pullgrp *pg, rvec x[], t_mdatoms *md, matrix box)
{
  int  i,ii,m;
  real wm;
  dvec com;

  clear_dvec(com);
  for(i=0; i<pg->ngx; i++) {
    ii = pg->idx[i];
    wm = md->massT[ii];
    if (pg->nweight > 0)
      wm *= pg->weight[i];
    for(m=0; m<DIM; m++)
      com[m] += wm*x[i][m];
  }
  for(m=0; m<DIM; m++)
    com[m] *= pg->wscale*pg->invtm;
  put_dvec_in_box(box,com);
  copy_dvec(com,pg->x_unc);
}

void calc_running_com(t_pull *pull) {
  int i,j,n;
  dvec ave;
  double tm;

  /* for group i, we have nhist[i] points of history. The maximum nr of points
     is pull->reflag. The array comhist[i][0..nhist[i]-1] has the positions
     of the center of mass over the last nhist[i] steps. x_unc[i] has the
     new ones. Remove the oldest one from the list, add the new one, calculate
     the average, put that in x_unc instead and return. We just move all coms
     1 down in the list and add the latest one to the top. 
   */

  if(pull->bCyl) {
    /* act on dyna groups */
    for(i=0;i<pull->ngrp;i++) {
      clear_dvec(ave);
      for(j=0;j<(pull->reflag-1);j++) {
        copy_dvec(pull->dyna[i].comhist[j+1],pull->dyna[i].comhist[j]);
        dvec_add(ave,pull->dyna[i].comhist[j],ave);
      }
      copy_dvec(pull->dyna[i].x_unc,pull->dyna[i].comhist[j]);
      dvec_add(ave,pull->dyna[i].x_unc,ave);
      dsvmul(1.0/pull->reflag,ave,ave); 

      /* now ave has the running com for group i, copy it to x_unc[i] */
      copy_dvec(ave,pull->dyna[i].x_unc);

#ifdef DEBUG
      if(pull->bVerbose)
        for(n=0;n<pull->reflag;n++)
          fprintf(stderr,"Comhist %d, grp %d: %8.3f%8.3f%8.3f\n",n,i,
                  pull->dyna[i].comhist[n][XX],
                  pull->dyna[i].comhist[n][YY],
                  pull->dyna[i].comhist[n][ZZ]);
#endif
    }
  } else {
    /* act on ref group */
    clear_dvec(ave);

    for(j=0;j<(pull->reflag-1);j++) {
      copy_dvec(pull->ref.comhist[j+1],pull->ref.comhist[j]);
      dvec_add(ave,pull->ref.comhist[j],ave);
    }

    copy_dvec(pull->ref.x_unc,pull->ref.comhist[j]);
    dvec_add(ave,pull->ref.x_unc,ave);
    dsvmul(1.0/pull->reflag,ave,ave); 
    /* now ave has the running com for group i, copy it to x_unc */
    copy_dvec(ave,pull->ref.x_unc);

#ifdef DEBUG
    if(pull->bVerbose)
      for(i=0;i<pull->reflag;i++)
        fprintf(stderr,"Comhist %d: %8.3f%8.3f%8.3f\n",i,
                pull->ref.comhist[i][XX],
                pull->ref.comhist[i][YY],
                pull->ref.comhist[i][ZZ]);
#endif
  }
}

void correct_t0_pbc(t_pull *pull, rvec x[], t_mdatoms *md, matrix box) {
  int i,ii,j,m;
  rvec dx;

  /* loop over all atoms in index for group i. Check if they moved
     more than half a box with respect to xp. If so add/subtract a box 
     from x0. Copy x to xp.
   */
  for(i=0;i<pull->ref.ngx;i++) {
    ii = pull->ref.idx[i];

    /* correct for jumps across the box */
    pbc_dx(x[ii],pull->ref.xp[i],dx);
    for(m=0; m<DIM; m++) {
      pull->ref.x0[i][m] += dx[m];
      pull->ref.xp[i][m]  = x[ii][m];
    }
  }
}

/* switch function, x between r and w */
real get_weight(real x, real r, real w) {
  
  /* Removed static stuff. Wasn't being used anyway */
  /*  static bool bFirst = TRUE;
  static real rw, a0, a1, a2, a3; */
  
  real weight; 

  /*  if (bFirst) {
    rw = r - w;
    a0 = (3*r*w*w-w*w*w)/(rw*rw*rw);
    a1 = -6*r*w/(rw*rw*rw);
    a2 = 3*(r+w)/(rw*rw*rw);
    a3 = -2/(rw*rw*rw);
    bFirst = FALSE;
  } 
  */

  if (x <= r)
    weight = 1;
  else if (x >= w)
    weight = 0;
  else
    weight = (w - x)/(w - r);

  return weight;
}

static double get_cylinder_distance(rvec x, dvec com, matrix box) {
  /* Triclinic compatible ??? */
  double dr, dx, dy, boxx, boxy;

  boxx = box[XX][XX];
  boxy = box[YY][YY];

  dx = fabs(x[XX] - com[XX]); 
  while(dx > 0.5*boxx) dx -= boxx;

  dy = fabs(x[YY] - com[YY]);
  while(dy > 0.5*boxy) dy -= boxy;

  dr = sqrt(dx*dx+dy*dy);
#ifdef CYLDEBUG
  fprintf(stderr,"x:%8.3f%8.3f%8.3f com:%8.3f%8.3f%8.3f dx,dy,dr:%8.3f%8.3f%8.3f\n",x[0],x[1],x[2],com[0],com[1],com[2],dx,dy,dr);
#endif
  return dr;
}

void make_refgrps(t_pull *pull,matrix box,t_mdatoms *md) 
{
  int i,ii,j,k,m;

  double dr,mass,wmass,wwmass;
  dvec test;
  t_pullgrp *pdyna;

  /* loop over all groups to make a reference group for each*/
  for(i=0;i<pull->ngrp;i++) {
    pdyna = &pull->dyna[i];
    k=0;
    wmass = 0;
    wwmass = 0;
    pdyna->ngx = 0;

    /* loop over all atoms in the main ref group */
    for(j=0;j<pull->ref.ngx;j++) {
      ii = pull->ref.idx[j];

      /* get_distance takes pbc into account */
      dr = get_cylinder_distance(pull->ref.x0[j],pull->grp[i].x_unc,box);

      if(dr < pull->rc) {
        /* add to index, to sum of COM, to weight array */
        mass = md->massT[ii];
        pdyna->ngx++;
        pdyna->weight[k] = get_weight(dr,pull->r,pull->rc);
        pdyna->idx[k] = ii;
        for(m=0;m<DIM;m++)
          pdyna->x_unc[m] += mass*pdyna->weight[k]*pull->ref.x0[j][m];
        wmass += mass*pdyna->weight[k];
	wwmass += mass*sqr(pdyna->weight[k]);
        k++;
      }
    }

    pdyna->wscale = wmass/wwmass;
    pdyna->invtm = 1.0/(pdyna->wscale*wmass);

    /* normalize the new 'x_unc' */
    dsvmul(1/wmass,pdyna->x_unc,pdyna->x_unc);
    if(pull->bVerbose)
      fprintf(stderr,"Made group %d:%8.3f%8.3f%8.3f m:%8.3f\n",
              i,pdyna->x_unc[0],pdyna->x_unc[1],
              pdyna->x_unc[2],1.0/pdyna->invtm);
  }
}

