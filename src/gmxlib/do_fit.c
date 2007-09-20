/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.2
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2007, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "maths.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "nrjac.h"
#include "vec.h"
#include "txtdump.h"
#include "smalloc.h"

#define EPS 1.0e-09

real calc_similar_ind(bool bRho,int nind,atom_id *index,real mass[],
		      rvec x[],rvec xp[])
{
  int i, j, d;
  real m, tm, xs, xd, rs, rd;
  
  tm=0;
  rs=0;
  rd=0;
  for(j=0; j<nind; j++) {
    if (index)
      i = index[j];
    else
      i = j;
    m = mass[i];
    tm += m;
    for(d=0 ; d<DIM; d++) {
      xd = x[i][d] - xp[i][d];
      rd += m * sqr(xd);
      if (bRho) {
	xs = x[i][d] + xp[i][d];
	rs += m * sqr(xs);
      }
    }
  }
  if (bRho)
    return 2*sqrt(rd/rs);
  else
    return sqrt(rd/tm);
}

real rmsdev_ind(int nind,atom_id index[],real mass[],rvec x[],rvec xp[])
{
  return calc_similar_ind(FALSE, nind, index, mass, x, xp);
}

real rmsdev(int natoms,real mass[],rvec x[],rvec xp[])
{
  return calc_similar_ind(FALSE, natoms, NULL, mass, x, xp);
}

real rhodev_ind(int nind,atom_id index[],real mass[],rvec x[],rvec xp[])
{
  return calc_similar_ind(TRUE, nind, index, mass, x, xp);
}

real rhodev(int natoms,real mass[],rvec x[],rvec xp[])
{
  return calc_similar_ind(TRUE, natoms, NULL, mass, x, xp);
}

void calc_fit_R(int natoms,real *w_rls,rvec *xp,rvec *x,matrix R)
{
  int    c,r,n,j,m,i,irot;
  double **omega,**om;
  double d[2*DIM],xnr,xpc;
  matrix vh,vk,u;
  real   mn;
  int    index;
  real   max_d;

  snew(omega,2*DIM);
  snew(om,2*DIM);
  for(i=0; i<2*DIM; i++) {
    snew(omega[i],2*DIM);
    snew(om[i],2*DIM);
  }
  
  for(i=0; i<2*DIM; i++) {
    d[i]=0;
    for(j=0; j<2*DIM; j++) {
      omega[i][j]=0;
      om[i][j]=0;
    }
  }
  
  /*calculate the matrix U*/
  clear_mat(u);
  for(n=0;(n<natoms);n++)
    if ((mn = w_rls[n]) != 0.0)
      for(c=0; (c<DIM); c++) {
	xpc=xp[n][c];
	for(r=0; (r<DIM); r++) {
	  xnr=x[n][r];
	  u[c][r]+=mn*xnr*xpc;
	}
      }
  
  /*construct omega*/
  /*omega is symmetric -> omega==omega' */
  for(r=0; r<2*DIM; r++)
    for(c=0; c<=r; c++)
      if (r>=DIM && c<DIM) {
        omega[r][c]=u[r-DIM][c];
        omega[c][r]=u[r-DIM][c];
      } else {
        omega[r][c]=0;
        omega[c][r]=0;
      }
  
  /*determine h and k*/
  jacobi(omega,2*DIM,d,om,&irot);
  /*real   **omega = input matrix a[0..n-1][0..n-1] must be symmetric
   *int     natoms = number of rows and columns
   *real      NULL = d[0]..d[n-1] are the eigenvalues of a[][]
   *real       **v = v[0..n-1][0..n-1] contains the vectors in columns
   *int      *irot = number of jacobi rotations
   */
  
  if (debug && irot==0) fprintf(debug,"IROT=0\n");
  
  index=0; /* For the compiler only */

  /* Copy only the first two eigenvectors */  
  for(j=0; j<2; j++) {
    max_d=-1000;
    for(i=0; i<2*DIM; i++)
      if (d[i]>max_d) {
        max_d=d[i];
        index=i;
      }
    d[index]=-10000;
    for(i=0; i<DIM; i++) {
      vh[j][i]=M_SQRT2*om[i][index];
      vk[j][i]=M_SQRT2*om[i+DIM][index];
    }
  }
  /* Calculate the last eigenvector as the outer-product of the first two.
   * This insures that the conformation is not mirrored and
   * prevents problems with completely flat reference structures.
   */  
  oprod(vh[0],vh[1],vh[2]);
  oprod(vk[0],vk[1],vk[2]);

  /*determine R*/
  for(r=0; r<DIM; r++)
    for(c=0; c<DIM; c++)
      R[r][c] = vk[0][r]*vh[0][c] +
	        vk[1][r]*vh[1][c] +
	        vk[2][r]*vh[2][c];

  for(i=0; i<2*DIM; i++) {
    sfree(omega[i]);
    sfree(om[i]);
  }
  sfree(omega);
  sfree(om);
}

void do_fit(int natoms,real *w_rls,rvec *xp,rvec *x)
{
  int    i,j,m,r,c;
  matrix R;
  rvec   x_old;

  /* Calculate the rotation matrix R */
  calc_fit_R(natoms,w_rls,xp,x,R);

  /*rotate X*/
  for(j=0; j<natoms; j++) {
    for(m=0; m<DIM; m++)
      x_old[m]=x[j][m];
    for(r=0; r<DIM; r++) {
      x[j][r]=0;
      for(c=0; c<DIM; c++)
        x[j][r]+=R[r][c]*x_old[c];
    }
  }
}

void reset_x(int ncm,atom_id ind_cm[],
	     int nreset,atom_id *ind_reset,rvec x[],real mass[])
{
  int  i,m,ai;
  rvec xcm;
  real tm,mm;
  
  tm=0.0;
  clear_rvec(xcm);
  for(i=0; i<ncm; i++) {
    ai=ind_cm[i];
    mm=mass[ai];
    for(m=0; m<DIM; m++)
      xcm[m]+=mm*x[ai][m];
    tm+=mm;
  }
  for(m=0; m<DIM; m++)
    xcm[m]/=tm;
    
  if (ind_reset)
    for(i=0; i<nreset; i++)
      rvec_dec(x[ind_reset[i]],xcm);
  else
    for(i=0; i<nreset; i++)
      rvec_dec(x[i],xcm);
}

