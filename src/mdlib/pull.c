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

  /*  Correct for particles jumping across the box relative to t=0 if any of the
      T0 reference types were selected
  */
  if(pull->reftype == eComT0 || pull->reftype == eDynT0) {
    correct_t0_pbc(pull, x, md, box);
  } else {
    for(i=0; i<pull->ref.ngx[0];i++)
      copy_rvec(x[pull->ref.idx[0][i]], pull->ref.x0[0][i]);
  }

  /* Get centers of mass of pulled groups */
  for(i=0; i<pull->pull.n; i++) {
    calc_com(x, pull->pull.ngx[i], pull->pull.idx[i], md, pull->pull.x_unc[i],box);
  }

  /* Get ceneteros of mass of reference groups.  Don't do this for AFM*/
  if(!pull->runtype == eAfm) {
    /* Dynamic case */
    if(pull->bCyl) {
      if(step % pull->update == 0)      /* Do we need to make new ones? */
        make_refgrps(pull, box, md);
      else {
        for(i=0; i<pull->pull.n; i++) {
          calc_com2(pull->ref.x0[0], pull->dyna.ngx[i], pull->dyna.idx[i],
                    md, pull->dyna.x_unc[i], box);
        }
      }
    }

    /* Normal case */
    if(!pull->bCyl)
      calc_com2(pull->ref.x0[0], pull->ref.ngx[0], pull->ref.idx[0], md,
                pull->ref.x_unc[0], box);

    /* Do running average if necessary */
    if(pull->reflag > 1)
      calc_running_com(pull);
  }
}

/* Apply forces in a mass weighted fashion */
void pull_do_forces(t_pull * pull, t_mdatoms * md, rvec * f, int start,
                    int homenr)
{
  int i, j, ii, m;
  real mass;
  real total_mass;

  for(i=0; i<pull->pull.n; i++) {
    total_mass = pull->pull.tmass[i];

    for(j=0; j<pull->pull.ngx[i]; j++) {
      ii=pull->pull.idx[i][j];
      mass = md->massT[ii];

      for(m=0; m<DIM; m++) {
        f[ii][m] += mass * pull->pull.f[i][m] / total_mass;
      }
    }
  }
}
     




/* check if we're close enough to the target position */
static bool check_convergence(t_pull *pull) {
  bool bTest = TRUE;
  real tol;
  int i,m;
  rvec tmp1;
  rvec dr; /* distance between target and current position */

  tol = pull->tolerance;

  for(i=0;i<pull->pull.n;i++) {
    if(pull->bCyl)
      rvec_sub(pull->dyna.x_unc[i],pull->pull.x_unc[i],tmp1);
    else
      rvec_sub(pull->ref.x_unc[0],pull->pull.x_unc[i],tmp1);

    rvec_add(pull->pull.xtarget[i],tmp1,dr);
    /* multiply by pull->dims to pick the elements we want */
    for(m=0;m<DIM;m++)
      dr[m] = pull->dims[m]*dr[m];
    bTest = bTest && (norm(dr) < tol);
  }

  return bTest;
}

static void do_start(t_pull *pull, rvec *x, rvec *f, matrix box, t_mdatoms *md, 
                     real dt, int step, t_topology *top, int start, int homenr)
{
  int i,j,ii,m;
  rvec dr,dx,tmp;     
  bool bThereYet,bDump;
  static int nout = 0;
  rvec ds;
  real k;

  /* pull.xtarget[i] is the desired position of group i with respect
     to the reference group */

  /* bThereYet is true if all groups have reached their desired position */
  bThereYet = check_convergence(pull);

  /* If all groups haven't converged */
  if(!bThereYet) {
    for(i=0;i<pull->pull.n;i++) {
      /* Calculate vector from starting position */
      if(pull->bCyl)
        rvec_sub(pull->dyna.x_unc[i], pull->pull.x_unc[i], tmp);
      else
        rvec_sub(pull->ref.x_unc[0], pull->pull.x_unc[i], tmp);

      rvec_add(tmp, pull->pull.xtarget[i], dr);

      /* multiply by pull->dims to pick the elements we want */
      for(m=0;m<DIM;m++)
        dr[m] = pull->dims[m]*dr[m];

      /* calculate current value of k */
      k = pull->start_k0 + (pull->start_k1 - pull->start_k0)/pull->k_rate *
          pull->k_step;

      svmul(k, dr, pull->pull.f[i]);
    }

    if((pull->k_step) < (pull->k_rate)) pull->k_step++;

    /* Distribute forces */
    pull_do_forces(pull, md, f, start, homenr);

  } else {  /* All groups are converged -> Do output */

    for(i=0;i<pull->pull.n;i++) {
      for(m=0;m<DIM;m++) {
        pull->pull.xtarget[i][m] += pull->pull.dir[i][m]*pull->xlt_incr/
                                    norm(pull->pull.dir[i]);
      }

      if(pull->bVerbose)
        fprintf(stderr,"New target position: %8.3f%8.3f%8.3f\n",
                pull->pull.xtarget[i][0],pull->pull.xtarget[i][1],
                pull->pull.xtarget[i][2]);
    }
    dump_conf(pull,x,box,top,nout,step*dt/1000); 
    nout++;
    pull->k_step=0;
  }
}

static void do_umbrella(t_pull *pull, rvec *x,rvec *f,matrix box, 
                        t_mdatoms *md, int start, int homenr) 
{
  int i,ii=0,j,m,g;
  real mi;
  rvec cm;    /* center of mass displacement of reference */
  rvec dr;    /* distance from com to potential center */

  /* loop over the groups that are being umbrella sampled */
  for(i=0;i<pull->pull.n;i++) {

    /* Fix distance stuff
       pull->ref.x_unc[0] has position of reference group
       pull->x_ref is the coordinate that we would like to be at
       pull->spring has the corrected position of the pulled group
    */

    /* Set desired position in x_ref[i] */
    if(pull->bCyl)    /* Dynamic Case */
      rvec_add(pull->dyna.x_unc[i], pull->UmbPos[i], pull->pull.x_ref[i]);
    else              /* Normal Case */
      rvec_add(pull->ref.x_unc[0],  pull->UmbPos[i], pull->pull.x_ref[i]);


    /* Find vector between current and desired positions */
    rvec_sub(pull->pull.x_ref[i], pull->pull.x_unc[i], dr);

    /* Select the components we want */
    for(m=DIM-1;m>=0;m--) {
      if(dr[m] > 0.5 * box[m][m]) rvec_dec(dr,box[m]);
      if(dr[m] <-0.5 * box[m][m]) rvec_inc(dr,box[m]);
      dr[m] *= pull->dims[m];
    }
    copy_rvec(dr,pull->pull.spring[i]);

    /* f = um_width*dx */
    svmul(pull->UmbCons[i],dr,pull->pull.f[i]);
  }
  pull_do_forces(pull, md, f, start, homenr);
}

/* this implements a constraint run like SHAKE does. */
static void do_constraint(t_pull *pull, rvec *x, matrix box, t_mdatoms *md, 
                          real dt, int *niter) 
{

  rvec r_ij, /* x_con[i] com of i in prev. step. Obeys constr. -> r_ij   */
  unc_ij,    /* x_unc[i] com of i this step, before constr.    -> unc_ij */
  ref_ij;    /* x_ref[i] com of i at t0, not updated           -> ref_ij */

  rvec *rinew;           /* current 'new' position of group i */
  rvec *rjnew;           /* current 'new' position of group j */
  real *direction;       /* direction of dr relative to r_ij  */
  double lambda, rm, mass;
  bool bConverged = FALSE;
  int n=0,i,ii,j,m,max_iter=1000;
  int ref;
  double x1,x2,q,a,b,c;  /* for solving the quadratic equation, 
                            see Num. Recipes in C ed 2 p. 184 */
  rvec *dr;              /* correction for group i */
  rvec *ref_dr;          /* correction for group j */
  rvec tmp,tmp2,tmp3,sum;

  if(pull->bCyl) {
    snew(ref_dr,pull->pull.n);
    snew(rjnew,pull->pull.n);
  } else {
    snew(ref_dr,1);
    snew(rjnew,1);
  }
  snew(dr,pull->pull.n);
  snew(rinew,pull->pull.n);
  snew(direction,pull->pull.n);

  /* copy the current unconstraint positions for use in iterations. We 
     iterate until rinew[i] and rjnew[j] obey the constraints. Then
     rinew - pull.x_unc[i] is the correction dr to group i */
  for(i=0;i<pull->pull.n;i++)
    copy_rvec(pull->pull.x_unc[i],rinew[i]);
  if(pull->bCyl)
    for(i=0;i<pull->pull.n;i++)
      copy_rvec(pull->dyna.x_unc[i],rjnew[i]);
  else
    copy_rvec(pull->ref.x_unc[0],rjnew[0]);

  while(!bConverged && n<max_iter) {
    /* loop over all constraints */
    for(i=0;(i<pull->pull.n);i++) {

      if(pull->bVerbose)
        fprintf(stderr,"group %d, iteration %d\n",i,n);

      if(pull->bCyl) {
        rvec_sub(pull->dyna.x_con[i],pull->pull.x_con[i],r_ij);
        rvec_sub(rjnew[i],rinew[i],unc_ij);
        rvec_sub(pull->dyna.x_ref[i],pull->pull.x_ref[i],ref_ij);
      } else {
        rvec_sub(pull->ref.x_con[0],pull->pull.x_con[i],r_ij);
        rvec_sub(rjnew[0],rinew[i],unc_ij);
        rvec_sub(pull->ref.x_ref[0],pull->pull.x_ref[i],ref_ij);
      }

      for(m=DIM-1; m>=0; m--) {
        /* correct for PBC */
        if(r_ij[m]   < -0.5*box[m][m]) rvec_inc(r_ij,box[m]);
        if(r_ij[m]   >  0.5*box[m][m]) rvec_dec(r_ij,box[m]);
        if(unc_ij[m] < -0.5*box[m][m]) rvec_inc(unc_ij,box[m]);
        if(unc_ij[m] >  0.5*box[m][m]) rvec_dec(unc_ij,box[m]);
        if(ref_ij[m] < -0.5*box[m][m]) rvec_inc(ref_ij,box[m]);
        if(ref_ij[m] >  0.5*box[m][m]) rvec_dec(ref_ij,box[m]);
        /* select components we want */
        r_ij[m]     *= pull->dims[m];
        unc_ij[m]   *= pull->dims[m];
        ref_ij[m]   *= pull->dims[m];
      }

      if(pull->bCyl)
        rm = 1/pull->pull.tmass[i] + 1/pull->dyna.tmass[i];
      else
        rm = 1/pull->pull.tmass[i] + 1/pull->ref.tmass[0];

      a = iprod(r_ij,r_ij)*dt*dt*dt*dt*rm*rm; 
      b = iprod(unc_ij,r_ij)*2*dt*dt*rm;
      c = iprod(unc_ij,unc_ij) - norm2(ref_ij);

      if(b<0)
        q = -0.5*(b-sqrt(b*b-4*a*c));
      else
        q = -0.5*(b+sqrt(b*b-4*a*c));
      x1 = q/a; x2 = c/q;
      lambda = x1 > 0 ? x1 : x2;

      if(pull->bVerbose)
        fprintf(stderr,"\nax^2+bx+c=0: a=%e b=%e c=%e\n"
                "x1=%e x2=%e sum:%e,%e, lambda:%e\n",a,b,c,x1,x2,
                a*x1*x1+b*x1+c,a*x2*x2+b*x2+c,lambda);

      /* the position corrections dr due to the constraint are: */
      if(pull->bCyl) {
        svmul(-dt*dt*lambda/pull->pull.tmass[i],r_ij,dr[i]);
        svmul(dt*dt*lambda/pull->dyna.tmass[i],r_ij,ref_dr[i]);
      } else {
        svmul(-dt*dt*lambda/pull->pull.tmass[i],r_ij,dr[i]);
        svmul(dt*dt*lambda/pull->ref.tmass[0],r_ij,ref_dr[0]);
      }

      /* and the direction of the constraint force: */
      direction[i] = cos_angle(r_ij,dr[i]);

      /* DEBUG */
      if(pull->bVerbose) {
        fprintf(stderr,"Direction: %f\n",direction[i]);
        if(pull->bCyl) {
          rvec_sub(rinew[i],rjnew[i],tmp);
          rvec_sub(pull->pull.x_ref[i],pull->dyna.x_ref[i],tmp2);
        } else {
          rvec_sub(pull->pull.x_ref[i],pull->ref.x_ref[0],tmp2);
          rvec_sub(rinew[i],rjnew[0],tmp);
        }
        rvec_sub(dr[i],ref_dr[0],tmp3);
        for(m=DIM-1; m>=0; m--) {
          if(tmp[m]  < -0.5*box[m][m]) rvec_inc(tmp,box[m]);
          if(tmp[m]  >  0.5*box[m][m]) rvec_dec(tmp,box[m]);
          if(tmp2[m] < -0.5*box[m][m]) rvec_inc(tmp2,box[m]);
          if(tmp2[m] >  0.5*box[m][m]) rvec_dec(tmp2,box[m]);
          if(tmp3[m] < -0.5*box[m][m]) rvec_inc(tmp3,box[m]);
          if(tmp3[m] >  0.5*box[m][m]) rvec_dec(tmp3,box[m]);
          tmp[m]  *= pull->dims[m];
          tmp2[m] *= pull->dims[m];
          tmp3[m] *= pull->dims[m];
        }

        if(pull->bCyl)
          fprintf(stderr,
                  "cur. i:%f %f %f j:%f %f %f d: %f\n"
                  "ref. i:%f %f %f j:%f %f %f d: %f\n"
                  "cor. i:%f %f %f j:%f %f %f d: %f\n",
                  rinew[i][0],rinew[i][1],rinew[i][2], 
                  rjnew[i][0],rjnew[i][1],rjnew[i][2], norm(tmp),
                  pull->pull.x_ref[i][0],pull->pull.x_ref[i][1],
                  pull->pull.x_ref[i][2],pull->dyna.x_ref[i][0],
                  pull->dyna.x_ref[i][1],pull->dyna.x_ref[i][2],
                  norm(tmp2),
                  dr[i][0],dr[i][1],dr[i][2],
                  ref_dr[0][0],ref_dr[0][1],ref_dr[0][2],
                  norm(tmp3));
        else
          fprintf(stderr,
                  "cur. i:%f %f %f j:%f %f %f d: %f\n"
                  "ref. i:%f %f %f j:%f %f %f d: %f\n"
                  "cor. i:%f %f %f j:%f %f %f d: %f\n",
                  rinew[i][0],rinew[i][1],rinew[i][2], 
                  rjnew[0][0],rjnew[0][1],rjnew[0][2], norm(tmp),
                  pull->pull.x_ref[i][0],pull->pull.x_ref[i][1],
                  pull->pull.x_ref[i][2],pull->ref.x_ref[0][0],
                  pull->ref.x_ref[0][1],pull->ref.x_ref[0][2],
                  norm(tmp2),
                  dr[i][0],dr[i][1],dr[i][2],
                  ref_dr[0][0],ref_dr[0][1],ref_dr[0][2],
                  norm(tmp3));
      } /* END DEBUG */

      /* update the positions with dr */
      rvec_add(rinew[i],dr[i],rinew[i]);

      if(pull->bCyl) {
        rvec_add(rjnew[i],ref_dr[i],rjnew[i]);

        /* calculate new distance between the two groups */
        rvec_sub(rjnew[i],rinew[i],unc_ij);

        /* select components and check PBC again */
        for(m=DIM-1; m>=0; m--) {
          if(unc_ij[m] < -0.5*box[m][m]) rvec_inc(unc_ij,box[m]);
          if(unc_ij[m] >  0.5*box[m][m]) rvec_dec(unc_ij,box[m]);
          unc_ij[m] *= pull->dims[m];
        }
      } else {
        rvec_add(rjnew[0],ref_dr[0],rjnew[0]);

        /* calculate new distance between the two groups */
        rvec_sub(rjnew[0],rinew[i],unc_ij);

        /* select components again and check PBC again */
        for(m=DIM-1; m>=0; m--) {
          if(unc_ij[m] < -0.5*box[m][m]) rvec_inc(unc_ij,box[m]);
          if(unc_ij[m] >  0.5*box[m][m]) rvec_dec(unc_ij,box[m]);
          unc_ij[m] *= pull->dims[m];
        }
      }
    }

    /* check if all constraints are fullfilled now */
    bConverged = TRUE;
    for(i=0;i<pull->pull.n;i++) {
      if(pull->bCyl) {
        rvec_sub(rjnew[i],rinew[i],unc_ij);
        rvec_sub(pull->dyna.x_ref[i],pull->pull.x_ref[i],ref_ij);
      } else {
        rvec_sub(rjnew[0],rinew[i],unc_ij);
        rvec_sub(pull->ref.x_ref[0],pull->pull.x_ref[i],ref_ij);
      }

      for(m=DIM-1; m>=0; m--) {
        if(unc_ij[m] < -0.5*box[m][m]) rvec_inc(unc_ij,box[m]);
        if(unc_ij[m] >  0.5*box[m][m]) rvec_dec(unc_ij,box[m]);
        if(ref_ij[m] < -0.5*box[m][m]) rvec_inc(ref_ij,box[m]);
        if(ref_ij[m] >  0.5*box[m][m]) rvec_dec(ref_ij,box[m]);
        ref_ij[m] *= pull->dims[m];
        unc_ij[m] *= pull->dims[m];
      }

      bConverged = bConverged && (fabs(norm(unc_ij)-norm(ref_ij)) < 
                                  pull->constr_tol);
    }


    /* DEBUG */
    if(pull->bVerbose) {
      if(!bConverged)
        fprintf(stderr,"NOT CONVERGED YET: Group %d (%s):"
                "d_ref = %f, current d = %f\n",
                i,pull->pull.grps[i], norm(ref_ij),norm(unc_ij));
    } /* END DEBUG */

    n++;
    /* if after all constraints are dealt with and bConverged is still TRUE
       we're finished, if not we do another iteration */
  }
  if(n>max_iter)
    fatal_error(0,"Too many iterations for constraint run");

  /* DONE ITERATING, NOW UPDATE COORDINATES AND CALC. CONSTRAINT FORCES */

  /* update the normal groups */
  for(i=0;i<pull->pull.n;i++) {
    /* get the final dr and constraint force for group i */
    rvec_sub(rinew[i],pull->pull.x_unc[i],dr[i]);
    /* select components of dr */
    for(m=0;m<DIM;m++)
      dr[i][m] *= pull->dims[m];
    svmul(pull->pull.tmass[i]/(dt*dt),dr[i],tmp);
    /* get the direction of dr */
    pull->pull.f[i][ZZ] = -norm(tmp)*direction[i];

    /* copy the new x_unc to x_con */
    copy_rvec(rinew[i],pull->pull.x_con[i]);

    /* update the atom positions */
    clear_rvec(sum);
    for(j=0;j<pull->pull.ngx[i];j++) {
      ii = pull->pull.idx[i][j];
      rvec_add(x[ii],dr[i],x[ii]);
      svmul(md->massT[ii],dr[i],tmp);
      rvec_add(tmp,sum,sum);
    }
    if(pull->bVerbose)
      fprintf(stderr,"Group %i: correction %e %e %e\n",
              i,sum[0],sum[1],sum[2]);
  }

  /* update the reference groups */
  if(pull->bCyl) {
    /* update the dynamic reference groups */
    for(i=0;i<pull->pull.n;i++) {
      rvec_sub(rjnew[i],pull->dyna.x_unc[i],ref_dr[i]);
      /* copy the new x_unc to x_con */
      copy_rvec(rjnew[i],pull->dyna.x_con[i]);
      /* select components of ref_dr */
      for(m=0;m<DIM;m++)
        ref_dr[i][m] *= pull->dims[m];

      clear_rvec(sum);
      for(j=0;j<pull->dyna.ngx[i];j++) {
        /* reset the atoms with dr, weighted by w_i */
        svmul(pull->dyna.weights[i][j],ref_dr[i],tmp); 
        ii = pull->dyna.idx[i][j];
        rvec_add(x[ii],tmp,x[ii]);
        svmul(md->massT[ii],tmp,tmp2);
        rvec_add(tmp2,sum,sum);
      }
      if(pull->bVerbose)
        fprintf(stderr,"Dyna grp %i: correction %e %e %e\n",
                i,sum[0],sum[1],sum[2]);
    }
  } else {
    /* update the reference group */
    rvec_sub(rjnew[0],pull->ref.x_unc[0], ref_dr[0]); 
    /* copy the new x_unc to x_con */
    copy_rvec(rjnew[0],pull->ref.x_con[0]);
    /* select components of ref_dr */
    for(m=0;m<DIM;m++)
      ref_dr[0][m] *= pull->dims[m];

    clear_rvec(sum);
    for(j=0;j<pull->ref.ngx[0];j++) {
      ii = pull->ref.idx[0][j];
      rvec_add(x[ii],ref_dr[0],x[ii]);
      svmul(md->massT[ii],ref_dr[0],tmp);
      rvec_add(tmp,sum,sum);
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
static void do_afm(t_pull *pull,rvec *f,matrix box,t_mdatoms *md, int start, int homenr, int step)
{
  int i,ii,j,m,g;
  real mi; 
  rvec cm;     /* center of mass displacement of reference */
  rvec dr;     /* extension of the springs */

  /* loop over the groups that are being pulled */
  for(i=0;i<pull->pull.n;i++) {
    /* compute how far the springs are stretched */
    for(m=DIM-1; m>=0; m--) {
      while(dr[m] >  0.5*box[m][m]) rvec_dec(dr,box[m]);
      while(dr[m] < -0.5*box[m][m]) rvec_inc(dr,box[m]);
      dr[m]=pull->dims[m]*(pull->pull.spring[i][m]-pull->pull.x_unc[i][m]);
    }

    /* calculate force from the spring on the pull group: f = - k*dr */
    for(m=0;m<DIM;m++)
      pull->pull.f[i][m] = pull->k*dr[m];

    /* move pulling spring along dir, over pull->rate  */
    for(m=0;m<DIM;m++)
      pull->pull.spring[i][m] = pull->pull.dir[i][m]*pull->rate * step;
  }
  /* Distribute forces over pulled groups */
  pull_do_forces(pull, md, f, start, homenr);
}

void pull(t_pull *pull,rvec *x,rvec *f,matrix box, t_topology *top, 
          real dt, int step, int natoms, t_mdatoms *md, int start, int homenr,
          t_commrec * cr) 
{
  int i,niter;
  static rvec *x_s = NULL;
  bool bShakeFirst;

  bShakeFirst = (f == NULL);

  if(!x_s)
    snew(x_s,md->nr); /* can't rely on natoms */

  /* copy x to temp array x_s. We assume all molecules are whole already */
  for(i=0;i<md->nr;i++)
    copy_rvec(x[i],x_s[i]);  

  switch(pull->runtype) {
  case eAfm:
    if(!bShakeFirst) {
      /* Update COM's */
      update_COM(pull, x_s, md, box, step);
      
      do_afm(pull,f,box,md,start,homenr, step);
      
      if(MASTER(cr))
        print_afm(pull,step);
    }
    break;

  case eStart:
    if(!bShakeFirst) {
      update_COM(pull, x_s, md, box, step);

      do_start(pull, x, f, box, md, dt, step, top, start, homenr);
      print_start(pull,step);
    }
    break; 

  case eConstraint:
    if(PAR(cr)) fatal_error(1, "Pull: Cannont run constraint force calculations in parallel.");
    
    update_COM(pull,x_s,md,box,step);

    /* print some debug info if necessary */
    if(pull->bVerbose) {
      if(pull->bCyl) {
        for(i=0;i<pull->pull.n;i++) {
          fprintf(stderr,"I      :%9.6f %9.6f %9.6f\n",pull->pull.x_unc[i][0],
                  pull->pull.x_unc[i][1],pull->pull.x_unc[i][2]);
          fprintf(stderr,"dyna xref: unconstr. com:%9.6f %9.6f %9.6f\n",
                  pull->dyna.x_unc[i][0],
                  pull->dyna.x_unc[i][1],pull->dyna.x_unc[i][2]);
        }
      } else {
        fprintf(stderr,"xref: unconstr. com:%9.6f %9.6f %9.6f\n",
                pull->ref.x_unc[0][0], pull->ref.x_unc[0][1],
                pull->ref.x_unc[0][2]);
      }
    }

    /* do the actual constraint calculation */
    do_constraint(pull,x,box,md,dt,&niter);
    print_constraint(pull,f,step,box,niter); 
    break;

  case eUmbrella:
    if(!bShakeFirst) {
      update_COM(pull, x_s, md, box, step);

      do_umbrella(pull,x,f,box,md,start,homenr);

      if(MASTER(cr))
        print_umbrella(pull,step);
    }
    break;

  case eTest:
    if(!bShakeFirst) {
      /* code to test reference groups, without actually doing anything 
         else 
      */
      (void)calc_com(x,pull->ref.ngx[0],pull->ref.idx[0],
                     md,pull->ref.x_unc[0],box);
      fprintf(stderr,"ref: %8.3f %8.3f %8.3f\n",pull->ref.x_unc[0][XX],
              pull->ref.x_unc[0][YY],pull->ref.x_unc[0][ZZ]);
      correct_t0_pbc(pull,x,md,box);
      (void)calc_com2(pull->ref.x0[0],pull->ref.ngx[0],pull->ref.idx[0],
                      md,pull->ref.x_unc[0],box);
      fprintf(stderr,"ref_t0: %8.3f %8.3f %8.3f\n",pull->ref.x_unc[0][XX],
              pull->ref.x_unc[0][YY],pull->ref.x_unc[0][ZZ]);
    }
    break;

  default:
    fatal_error(0,"undetermined runtype");
  }
}






