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

/* update COM of pulled and reference groups */
void update_COM(t_pull * pull, rvec * x, t_mdatoms * md, matrix box, int step)
{
  int i;

  /*  Correct for particles jumping across the box relative to t=0 if any of
      the T0 reference types were selected
  */
  if (pull->AbsoluteRef) {
    for(i=0; i<DIM; i++)
      pull->ref.x_unc[i] = 0;
  } else {
    if(pull->reftype == eComT0 || pull->reftype == eDynT0) {
      correct_t0_pbc(pull, x, md, box);
    } else {
      for(i=0; i<pull->ref.ngx;i++)
	copy_rvec(x[pull->ref.idx[i]], pull->ref.x0[i]);
    }

    /* Dynamic case */
    if (pull->bCyl) {
      if (step % pull->update == 0 )   /* Should we update? */
	make_refgrps(pull,box,md);
      else {
	for (i=0; i<pull->ngrp; ++i) {
	  calc_com2(&pull->dyna[i],pull->ref.x0,md,box);
	  if (pull->bVerbose)
	    fprintf(stderr,"dynacom: %8.3f%8.3f%8.3f\n",pull->dyna[i].x_unc[0],
		    pull->dyna[i].x_unc[1],pull->dyna[i].x_unc[2]);
	}
      }
    } else {
      /* Normal case */
      calc_com2(&pull->ref,pull->ref.x0,md,box);

      /* Do running average if necessary */
      if (pull->reflag > 1) {
	calc_running_com(pull);
      }
    }
  } /* Done updating reference group */
  
  /* Update COM of pulled groups */
  for (i=0; i<pull->ngrp; ++i)
    calc_com(&pull->grp[i],x,md,box);
}

/* Apply forces in a mass weighted fashion */
void apply_forces(t_pull * pull, t_mdatoms * md, rvec * f, int start,
		  int homenr)
{
  int i, j, ii, m;
  double wmass,inv_wm;
  t_pullgrp *pgrp;
  
  for(i=0; i<pull->ngrp; i++) {
    pgrp = &pull->grp[i];
    inv_wm = pgrp->wscale*pgrp->invtm;

    for(j=0; j<pgrp->ngx; j++) {
      ii=pgrp->idx[j];
      wmass = md->massT[ii];
      if (pgrp->nweight > 0)
	wmass *= pgrp->weight[i];

      for(m=0; m<DIM; m++) {
        f[ii][m] += wmass * pgrp->f[m] * inv_wm;
      }
    }
  }
}

static void do_umbrella(t_pull *pull, rvec *x,rvec *f,matrix box, 
                        t_mdatoms *md, int start, int homenr) 
{
  int i,ii=0,j,m,g;
  dvec cm;    /* center of mass displacement of reference */
  dvec dr;    /* distance from com to potential center */

  /* loop over the groups that are being umbrella sampled */
  for(i=0;i<pull->ngrp;i++) {

    /* Fix distance stuff
       pull->ref.x_unc[0] has position of reference group
       pull->x_ref is the coordinate that we would like to be at
       pull->spring has the corrected position of the pulled group
    */

    /* Set desired position in x_ref[i] */
    if(pull->bCyl)    /* Dynamic Case */
      dvec_add(pull->dyna[i].x_unc, pull->grp[i].UmbPos, pull->grp[i].x_ref);
    else              /* Normal Case */
      dvec_add(pull->ref.x_unc,  pull->grp[i].UmbPos, pull->grp[i].x_ref);


    /* Find vector between current and desired positions */
    d_pbc_dx(box, pull->grp[i].x_ref, pull->grp[i].x_unc, dr);

    /* Select the components we want */
    for(m=DIM-1;m>=0;m--) {
      dr[m] *= pull->dims[m];
    }
    copy_dvec(dr,pull->grp[i].spring);

    /* f = um_width*dx */
    dsvmul(pull->grp[i].UmbCons,dr,pull->grp[i].f);
  }
  apply_forces(pull, md, f, start, homenr);
}

/* this implements a constraint run like SHAKE does. */
static void do_constraint(t_pull *pull, rvec *x, matrix box, t_mdatoms *md, 
                          real dt, int step, int *niter) 
{

  dvec r_ij, /* x_con[i] com of i in prev. step. Obeys constr. -> r_ij   */
  unc_ij,    /* x_unc[i] com of i this step, before constr.    -> unc_ij */
  ref_ij;    /* x_ref[i] com of i at t0, not updated           -> ref_ij */

  dvec *rinew;           /* current 'new' position of group i */
  dvec *rjnew;           /* current 'new' position of group j */
  double *direction;       /* direction of dr relative to r_ij  */
  double d0,refinc,ref,inpr;
  double lambda, rm, mass;
  bool bConverged = FALSE;
  int n=0,i,ii,j,m,max_iter=1000;
  double q,a,b,c;  /* for solving the quadratic equation, 
		      see Num. Recipes in C ed 2 p. 184 */
  dvec *dr;              /* correction for group i */
  dvec *ref_dr;          /* correction for group j */
  dvec tmp,tmp2,tmp3,sum;
  t_pullgrp *pdyna,*pgrp;

  if(pull->bCyl) {
    snew(ref_dr,pull->ngrp);
    snew(rjnew,pull->ngrp);
  } else {
    snew(ref_dr,1);
    snew(rjnew,1);
  }
  snew(dr,pull->ngrp);
  snew(rinew,pull->ngrp);
  snew(direction,pull->ngrp);

  /* copy the current unconstraint positions for use in iterations. We 
     iterate until rinew[i] and rjnew[j] obey the constraints. Then
     rinew - pull.x_unc[i] is the correction dr to group i */
  for(i=0;i<pull->ngrp;i++)
    copy_dvec(pull->grp[i].x_unc,rinew[i]);
  if(pull->bCyl)
    for(i=0;i<pull->ngrp;i++)
      copy_dvec(pull->dyna[i].x_unc,rjnew[i]);
  else
    copy_dvec(pull->ref.x_unc,rjnew[0]);

  while(!bConverged && n<max_iter) {
    /* loop over all constraints */
    for(i=0; i<pull->ngrp; i++) {

      if(pull->bVerbose)
        fprintf(stderr,"\ngroup %d, iteration %d",i,n);

      if(pull->bCyl) {
        d_pbc_dx(box,pull->grp[i].x_con,pull->dyna[i].x_con,r_ij);
        d_pbc_dx(box,rinew[i],rjnew[i],unc_ij);
        d_pbc_dx(box,pull->grp[i].x_ref,pull->dyna[i].x_ref,ref_ij);
      } else {
        d_pbc_dx(box,pull->grp[i].x_con,pull->ref.x_con,r_ij);
        d_pbc_dx(box,rinew[i],rjnew[0],unc_ij);
        d_pbc_dx(box,pull->grp[i].x_ref,pull->ref.x_ref,ref_ij);
      }

      for(m=DIM-1; m>=0; m--) {
        /* select components we want */
        r_ij[m]     *= pull->dims[m];
        unc_ij[m]   *= pull->dims[m];
        ref_ij[m]   *= pull->dims[m];
      }

      d0     = pull->grp[i].constr_d0; 
      refinc = pull->grp[i].constr_rate*step*dt;

      if (pull->bDir) {
	a = 0;
	b = 0;
	for(m=0; m<DIM; m++) {
	  a += unc_ij[m]*pull->dir[m];
	  b += ref_ij[m]*pull->dir[m];
	}
	if (d0 > 0) {
	  ref = d0 + refinc;
	} else {
	  ref = b + refinc;
	}
	copy_dvec(pull->dir,r_ij);
	dsvmul(a,pull->dir,unc_ij);
	lambda = a - ref;
	if(pull->bVerbose)
	  fprintf(stderr,"\nlambda:%e\n",lambda);
      } else {
	a = diprod(r_ij,r_ij); 
	b = diprod(unc_ij,r_ij)*2;
	if (d0 > 0) {
	  ref = d0 + refinc;
	} else {
	  ref = dnorm(ref_ij) + refinc;
	}
	c = diprod(unc_ij,unc_ij) - sqr(ref);

	if (b < 0) {
	  q = -0.5*(b - sqrt(b*b - 4*a*c));
	  lambda = -q/a;
	} else {
	  q = -0.5*(b + sqrt(b*b - 4*a*c));
	  lambda = -c/q;
	}

	if(pull->bVerbose)
	  fprintf(stderr,"ax^2+bx+c=0: a=%e b=%e c=%e lambda=%e\n",
		  a,b,c,lambda);
      }
      
      /* the position corrections dr due to the constraint are: */
      if(pull->bCyl) {
	rm = 1.0/(pull->grp[i].invtm + pull->dyna[i].invtm);
        dsvmul(-lambda*rm*pull->grp[i].invtm, r_ij,dr[i]);
        dsvmul( lambda*rm*pull->dyna[i].invtm,r_ij,ref_dr[i]);
      } else {
	rm = 1.0/(pull->grp[i].invtm + pull->ref.invtm);
        dsvmul(-lambda*rm*pull->grp[i].invtm,r_ij,dr[i]);
        dsvmul( lambda*rm*pull->ref.invtm,   r_ij,ref_dr[0]);
      }

      /* and the direction of the constraint force: */
      if (dr[i][0] != 0 || dr[i][1] != 0 || dr[i][2] != 0) {
	direction[i] = diprod(r_ij,dr[i])
	  /sqrt(diprod(r_ij,r_ij)*diprod(dr[i],dr[i]));
      } else {
	direction[i] = 0;
      }

      /* DEBUG */
      if(pull->bVerbose) {
        fprintf(stderr,"Direction: %f\n",direction[i]);
        if(pull->bCyl) {
          d_pbc_dx(box,rinew[i],rjnew[i],tmp);
        } else {
          d_pbc_dx(box,rinew[i],rjnew[0],tmp);
        }
        d_pbc_dx(box,dr[i],ref_dr[0],tmp3);
        for(m=DIM-1; m>=0; m--) {
          tmp[m]  *= pull->dims[m];
          tmp2[m] *= pull->dims[m];
          tmp3[m] *= pull->dims[m];
        }

        if(pull->bCyl)
          fprintf(stderr,
                  "cur. i:%8f %8f %8f j:%8f %8f %8f d: %8f\n"
                  "ref.   %8s %8s %8s   %8s %8s %8s d: %8f\n"
                  "cor. i:%8f %8f %8f j:%8f %8f %8f d: %8f\n",
                  rinew[i][0],rinew[i][1],rinew[i][2], 
                  rjnew[i][0],rjnew[i][1],rjnew[i][2], dnorm(tmp),
                  "","","","","","",ref,
                  dr[i][0],dr[i][1],dr[i][2],
                  ref_dr[0][0],ref_dr[0][1],ref_dr[0][2],
                  dnorm(tmp3));
        else
          fprintf(stderr,
                  "cur. i:%8f %8f %8f j:%8f %8f %8f d: %8f\n"
                  "ref.   %8s %8s %8s   %8s %8s %8s d: %8f\n"
                  "cor. i:%8f %8f %8f j:%8f %8f %8f d: %8f\n",
                  rinew[i][0],rinew[i][1],rinew[i][2], 
                  rjnew[0][0],rjnew[0][1],rjnew[0][2], dnorm(tmp),
		  "","","","","","",ref,
                  dr[i][0],dr[i][1],dr[i][2],
                  ref_dr[0][0],ref_dr[0][1],ref_dr[0][2],
                  dnorm(tmp3));
      } /* END DEBUG */

      /* update the positions with dr */
      dvec_add(rinew[i],dr[i],rinew[i]);

      if(pull->bCyl) {
        dvec_add(rjnew[i],ref_dr[i],rjnew[i]);

        /* calculate new distance between the two groups */
        d_pbc_dx(box,rinew[i],rjnew[i],unc_ij);

        /* select components */
        for(m=DIM-1; m>=0; m--) {
          unc_ij[m] *= pull->dims[m];
        }
      } else {
        dvec_add(rjnew[0],ref_dr[0],rjnew[0]);

        /* calculate new distance between the two groups */
        d_pbc_dx(box,rinew[i],rjnew[0],unc_ij);

        /* select components again */
        for(m=DIM-1; m>=0; m--) {
          unc_ij[m] *= pull->dims[m];
        }
      }
    }

    /* check if all constraints are fullfilled now */
    for(i=0;i<pull->ngrp;i++) {
      if(pull->bCyl) {
        d_pbc_dx(box,rinew[i],rjnew[i],unc_ij);
      } else {
        d_pbc_dx(box,rinew[i],rjnew[0],unc_ij);
      }

      for(m=DIM-1; m>=0; m--) {
        unc_ij[m] *= pull->dims[m];
      }
      if (pull->bDir) {
	inpr = diprod(unc_ij,pull->dir);
	dsvmul(inpr,pull->dir,unc_ij);
	bConverged =
	  fabs(diprod(unc_ij,pull->dir) - ref) < pull->constr_tol;
      } else {
	bConverged =
	  fabs(dnorm(unc_ij) - ref) < pull->constr_tol;
      }

      /* DEBUG */
      if(pull->bVerbose) {
	if(!bConverged)
	  fprintf(stderr,"NOT CONVERGED YET: Group %d (%s):"
		  "d_ref = %f, current d = %f\n",
		  i,pull->grp[i].name,ref,dnorm(unc_ij));
      } /* END DEBUG */
    }

    n++;
    /* if after all constraints are dealt with and bConverged is still TRUE
       we're finished, if not we do another iteration */
  }
  if (n > max_iter)
    gmx_fatal(FARGS,"Too many iterations for constraint run");

  /* DONE ITERATING, NOW UPDATE COORDINATES AND CALC. CONSTRAINT FORCES */

  /* update the normal groups */
  for(i=0; i<pull->ngrp; i++) {
    pgrp = &pull->grp[i];
    /* get the final dr and constraint force for group i */
    dvec_sub(rinew[i],pgrp->x_unc,dr[i]);
    /* select components of dr */
    for(m=0; m<DIM; m++)
      dr[i][m] *= pull->dims[m];
    dsvmul(1.0/(pgrp->invtm*dt*dt),dr[i],tmp);
    /* get the direction of dr */
    pgrp->f[ZZ] = -dnorm(tmp)*direction[i];

    /* copy the new x_unc to x_con */
    copy_dvec(rinew[i],pgrp->x_con);

    /* update the atom positions */
    clear_dvec(sum);
    for(j=0;j<pgrp->ngx;j++) {
      ii = pgrp->idx[j];
      if (pgrp->nweight == 0) {
	copy_dvec(dr[i],tmp);
      } else {
	dsvmul(pgrp->wscale*pgrp->weight[j],dr[i],tmp); 
      }
      for(m=0; m<DIM; m++)
	x[ii][m] += tmp[m];
      dsvmul(md->massT[ii],tmp,tmp2);
      dvec_add(tmp2,sum,sum);
    }
    if(pull->bVerbose)
      fprintf(stderr,"Group %i: correction %e %e %e\n",
              i,sum[0],sum[1],sum[2]);
  }

  /* update the reference groups */
  if(pull->bCyl) {
    /* update the dynamic reference groups */
    for(i=0;i<pull->ngrp;i++) {
      pdyna = &pull->dyna[i];
      dvec_sub(rjnew[i],pdyna->x_unc,ref_dr[i]);
      /* copy the new x_unc to x_con */
      copy_dvec(rjnew[i],pdyna->x_con);
      /* select components of ref_dr */
      for(m=0;m<DIM;m++)
        ref_dr[i][m] *= pull->dims[m];

      clear_dvec(sum);
      for(j=0;j<pdyna->ngx;j++) {
        /* reset the atoms with dr, weighted by w_i */
        dsvmul(pdyna->wscale*pdyna->weight[j],ref_dr[i],tmp); 
        ii = pdyna->idx[j];
	for(m=0; m<DIM; m++)
	  x[ii][m] += tmp[m];
        dsvmul(md->massT[ii],tmp,tmp2);
        dvec_add(tmp2,sum,sum);
      }
      if(pull->bVerbose)
        fprintf(stderr,"Dyna grp %i: correction %e %e %e\n",
                i,sum[0],sum[1],sum[2]);
    }
  } else {
    /* update the reference group */
    dvec_sub(rjnew[0],pull->ref.x_unc, ref_dr[0]); 
    /* copy the new x_unc to x_con */
    copy_dvec(rjnew[0],pull->ref.x_con);
    /* select components of ref_dr */
    for(m=0;m<DIM;m++)
      ref_dr[0][m] *= pull->dims[m];

    clear_dvec(sum);
    for(j=0;j<pull->ref.ngx;j++) {
      ii = pull->ref.idx[j];
      if (pull->ref.nweight == 0) {
	copy_dvec(ref_dr[0],tmp);
      } else {
	dsvmul(pull->ref.wscale*pull->ref.weight[j],ref_dr[0],tmp); 
      }
      for(m=0; m<DIM; m++)
	x[ii][m] += tmp[m];
      dsvmul(md->massT[ii],tmp,tmp2);
      dvec_add(tmp2,sum,sum);
    }
    if(pull->bVerbose)
      fprintf(stderr,"Reference: correction %e %e %e\n",
              sum[0],sum[1],sum[2]);

  }

  /* finished! I hope. Give back some memory */
  sfree(ref_dr);
  sfree(rinew);
  sfree(rjnew);
  sfree(dr);
  sfree(direction);
  *niter = n;
}

/* mimicks an AFM experiment, groups are pulled via a spring */
static void do_afm(t_pull *pull, rvec *f, matrix box, t_mdatoms *md,
		   int start, int homenr, real dt, int step)
{
  int i,m;
  dvec dr;     /* extension of the springs */

  /* loop over the groups that are being pulled */
  for(i=0;i<pull->ngrp;i++) {

    /* move pulling spring along dir, over pull->rate  */
    for(m=0;m<DIM;m++)
      pull->grp[i].spring[m] = pull->grp[i].AfmInit[m] + pull->ref.x_unc[m] +
	pull->grp[i].AfmVec[m]*pull->grp[i].AfmRate*step*dt;

    d_pbc_dx(box,pull->grp[i].spring, pull->grp[i].x_unc, dr);
    
    /* compute how far the springs are stretched */
    for(m=DIM-1; m>=0; m--) {
      /* dr[m]=pull->dims[m]*(pull->grp.spring[i][m]-pull->grp.x_unc[i][m]);*/
      dr[m]*=pull->dims[m];
    }

    /* calculate force from the spring on the pull group: f = - k*dr */
    for(m=0;m<DIM;m++)
      pull->grp[i].f[m] = pull->grp[i].AfmK*dr[m];
  }
  /* Distribute forces over pulled groups */
  apply_forces(pull, md, f, start, homenr);
}

void pull(t_pull *pull,rvec *x,rvec *f,matrix box, t_topology *top, 
          real dt, int step, real time, t_mdatoms *md, int start, int homenr,
          t_commrec * cr) 
{
  int i,niter;
  rvec *x_s = NULL;
  bool bShakeFirst;

  bShakeFirst = (f == NULL);

  
  snew(x_s,md->nr); /* can't rely on natoms */

  /* copy x to temp array x_s. We assume all molecules are whole already */
  for(i=0;i<md->nr;i++)
    copy_rvec(x[i],x_s[i]);  
  
  /* Update COM's */
  update_COM(pull, x_s, md, box, step);

  switch(pull->runtype) {
  case eAfm:
    if(!bShakeFirst) {
      do_afm(pull,f,box,md,start,homenr, dt, step);
      
      if(MASTER(cr))
        print_afm(pull,step,time);
    }
    break;

  case eConstraint:
    /* do the actual constraint calculation */
    do_constraint(pull,x,box,md,dt,step,&niter);
    print_constraint(pull,step,time); 
    break;

  case eUmbrella:
    if(!bShakeFirst) {
      do_umbrella(pull,x,f,box,md,start,homenr);

      if(MASTER(cr))
        print_umbrella(pull,step,time);
    }
    break;

  default:
    gmx_fatal(FARGS,"undetermined runtype");
  }

  sfree(x_s);
}
