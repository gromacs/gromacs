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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
static char *SRCID_pullutil_c = "$Id$";
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
#include "pull-internal.h"

/* calculates center of mass of selection index from all coordines x */
real calc_com(rvec x[],int gnx,atom_id *index,t_mdatoms *md,
		     rvec com,matrix box)
{
  int  i,ii,m;
  real m0,tm;

  clear_rvec(com);
  tm=0;
  for(i=0; (i<gnx); i++) {
    ii=index[i];
    m0=md->massT[ii];
    tm+=m0;
    for(m=0; (m<DIM); m++)
      com[m]+=m0*x[ii][m];
  }
  svmul(1/tm,com,com);
  for(m=DIM-1; m>=0; m--) {
    if (com[m] < 0        ) rvec_inc(com,box[m]);
    if (com[m] > box[m][m]) rvec_dec(com,box[m]);
  }

  return tm;
}

/* calculates com of all atoms in x[], *index has their index numbers
   to get the masses from atom[] */
real calc_com2(rvec x[],int gnx,atom_id *index,t_mdatoms *md,rvec com,
	       matrix box)
{
  int  i,ii,m;
  real m0,tm;

  clear_rvec(com);
  tm=0;
  for(i=0; (i<gnx); i++) {
    ii=index[i];
    m0=md->massT[ii];
    tm+=m0;
    for(m=0; (m<DIM); m++)
      com[m]+=m0*x[i][m];
  }
  svmul(1/tm,com,com);
  for(m=DIM-1; m>=0; m--) {
    /* next two lines used to be commented out */
    if (com[m] < 0        ) rvec_inc(com,box[m]);
    if (com[m] > box[m][m]) rvec_dec(com,box[m]); 
  }
  return tm;
}

void calc_running_com(t_pull *pull) {
  int i,j,n;
  rvec ave;
  real tm;
  
  /* for group i, we have nhist[i] points of history. The maximum nr of points
     is pull->reflag. The array comhist[i][0..nhist[i]-1] has the positions
     of the center of mass over the last nhist[i] steps. x_unc[i] has the
     new ones. Remove the oldest one from the list, add the new one, calculate
     the average, put that in x_unc instead and return. We just move all coms
     1 down in the list and add the latest one to the top. 
   */
  
  if (pull->bCyl) {
    /* act on dyna groups */
    for (i=0;i<pull->pull.n;i++) {
      clear_rvec(ave);
      for (j=0;j<(pull->reflag-1);j++) {
	copy_rvec(pull->dyna.comhist[i][j+1],pull->dyna.comhist[i][j]);
	rvec_add(ave,pull->dyna.comhist[i][j],ave);
      }
      copy_rvec(pull->dyna.x_unc[i],pull->dyna.comhist[i][j]);
      rvec_add(ave,pull->dyna.x_unc[i],ave);
      svmul(1.0/pull->reflag,ave,ave); 

      /* now ave has the running com for group i, copy it to x_unc[i] */
      copy_rvec(ave,pull->dyna.x_unc[i]);

#ifdef DEBUG
      if (pull->bVerbose) 
	for (n=0;n<pull->reflag;n++) 
	  fprintf(stderr,"Comhist %d, grp %d: %8.3f%8.3f%8.3f\n",n,i,
		  pull->dyna.comhist[i][n][XX],
		  pull->dyna.comhist[i][n][YY],
  		  pull->dyna.comhist[i][n][ZZ]);
#endif
    }
  } else {
    /* act on ref group */
    clear_rvec(ave);
    
    for (j=0;j<(pull->reflag-1);j++) {
      copy_rvec(pull->ref.comhist[0][j+1],pull->ref.comhist[0][j]);
      rvec_add(ave,pull->ref.comhist[0][j],ave);
    }
    
    copy_rvec(pull->ref.x_unc[0],pull->ref.comhist[0][j]);
    rvec_add(ave,pull->ref.x_unc[0],ave);
    svmul(1.0/pull->reflag,ave,ave); 
    /* now ave has the running com for group i, copy it to x_unc[0] */
    copy_rvec(ave,pull->ref.x_unc[0]);

#ifdef DEBUG
    if (pull->bVerbose) 
      for (i=0;i<pull->reflag;i++) 
	fprintf(stderr,"Comhist %d: %8.3f%8.3f%8.3f\n",i,
		pull->ref.comhist[0][i][XX],
		pull->ref.comhist[0][i][YY],
		pull->ref.comhist[0][i][ZZ]);
#endif
  }
}

void correct_t0_pbc(t_pull *pull, rvec x[], t_mdatoms *md, matrix box) {
  int i,ii,j,m;
  real tm;
  rvec com,dx;

  /* loop over all atoms in index for group i. Check if they moved
     more than half a box with respect to xp. If so add/subtract a box 
     from x0. Copy x to xp.
   */
  for (i=0;i<pull->ref.ngx[0];i++) {
    ii = pull->ref.idx[0][i];
    
    /* correct for jumps across the box */
    rvec_sub(x[ii],pull->ref.xp[0][i],dx);
    for (m=DIM-1; m>=0; m--) {
      if (dx[m] < -0.5*box[m][m]) {
	rvec_inc(dx,box[m]);
	if (pull->bVerbose && pull->dims[m])
	  fprintf(stderr,"Jumped +box: nr %d dir: %d old:%8.3f\n",ii,m,
		  pull->ref.x0[0][i][m]); 
      }
      
      if (dx[m] >  0.5*box[m][m]) {
	rvec_dec(dx,box[m]);
	if (pull->bVerbose && pull->dims[m]) 
	  fprintf(stderr,"Jumped -box: nr %d dir: %d old:%8.3f\n",ii,m,
		  pull->ref.x0[0][i][m]); 
      }

      pull->ref.x0[0][i][m] += dx[m];
      pull->ref.xp[0][i][m]  = x[ii][m];
    }
  }
  tm = calc_com2(pull->ref.x0[0],pull->ref.ngx[0],pull->ref.idx[0],
		 md,com,box);
  if (pull->bVerbose) 
    fprintf(stderr,"correct_t0: Group %s: mass:%8.3f com:%8.3f%8.3f%8.3f\n",
	    pull->ref.grps[0],tm,com[0],com[1],com[2]);
}

void string2rvec(char buf[], rvec nums) {
  float a,b,c;
  if (sscanf(buf,"%f%f%f",&a,&b,&c) != 3)
    fatal_error(0,"Expected three numbers at input line %s",buf);
  nums[0]=a; nums[1]=b; nums[2]=c;
}

/* switch function, x between r and w */
real get_weight(real x, real r, real w) {
  static bool bFirst = TRUE;
  static real rw, a0, a1, a2, a3;
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
  
  if (x < r)
    weight = 1;
  else if (x > w)
    weight = 0;
  else  
    weight = -w/(r-w) + x/(r-w);
    
  return weight;
}
  
static real get_cylinder_distance(rvec x, rvec com, matrix box) {
  /* Triclinic compatible ??? */
  real dr, dx, dy, boxx, boxy;

  boxx = box[XX][XX];
  boxy = box[YY][YY];

  dx = fabs(x[XX] - com[XX]); 
  while (dx > 0.5*boxx) dx -= boxx;

  dy = fabs(x[YY] - com[YY]);
  while (dy > 0.5*boxy) dy -= boxy;
  
  dr = sqrt(dx*dx+dy*dy);
#ifdef CYLDEBUG
  fprintf(stderr,"x:%8.3f%8.3f%8.3f com:%8.3f%8.3f%8.3f dx,dy,dr:%8.3f%8.3f%8.3f\n",x[0],x[1],x[2],com[0],com[1],com[2],dx,dy,dr);
#endif
  return dr;
}

void make_refgrps(t_pull *pull,matrix box,t_mdatoms *md) 
{
  int ngrps,i,ii,j,k,m;
  static bool bFirst = TRUE;
  real dr,mass;
  real truemass;
  rvec test;

  ngrps = pull->pull.n;
  if (bFirst) {
    snew(pull->dyna.ngx,ngrps);
    snew(pull->dyna.idx,ngrps);
    snew(pull->dyna.weights,ngrps);
    for (i=0;i<ngrps;i++) {
      snew(pull->dyna.idx[i],pull->ref.ngx[0]);    /* more than nessary */
      snew(pull->dyna.weights[i],pull->ref.ngx[0]);
    }      
    bFirst = FALSE;
  }

  /* loop over all groups to make a reference group for each*/
  for (i=0;i<ngrps;i++) {
    k=0;
    truemass=0;
    pull->dyna.tmass[i] = 0;
    pull->dyna.ngx[i] = 0;

    /* loop over all atoms in the main ref group */
    for (j=0;j<pull->ref.ngx[0];j++) {
      ii = pull->ref.idx[0][j];

      /* get_distance takes pbc into account */
      dr = get_cylinder_distance(pull->ref.x0[0][j],pull->pull.x_unc[i],box);

      if (dr < pull->rc) {
	/* add to index, to sum of COM, to weight array */
	mass = md->massT[ii];
	truemass += mass;
	pull->dyna.ngx[i]++;
	pull->dyna.weights[i][k] = get_weight(dr,pull->r,pull->rc);
	pull->dyna.idx[i][k] = ii;
	for (m=0;m<DIM;m++) 
	  pull->dyna.x_unc[i][m] += mass*pull->dyna.weights[i][k]*
	    pull->ref.x0[0][j][m];
	pull->dyna.tmass[i] += mass*pull->dyna.weights[i][k];
	k++;
      }
    }

    /* normalize the new 'x_unc' */
    svmul(1/pull->dyna.tmass[i],pull->dyna.x_unc[i],pull->dyna.x_unc[i]);
    if (pull->bVerbose) 
      fprintf(stderr,"Made group %d:%8.3f%8.3f%8.3f wm:%8.3f m:%8.3f\n",
	      i,pull->dyna.x_unc[i][0],pull->dyna.x_unc[i][1],
	      pull->dyna.x_unc[i][2],pull->dyna.tmass[i],truemass);
  }
}

