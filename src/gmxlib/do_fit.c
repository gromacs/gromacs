/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * GROningen MAchine for Chemical Simulation
 */
static char *SRCID_do_fit_c = "$Id$";

#include "maths.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "gstat.h"
#include "vec.h"
#include "txtdump.h"

#define EPS  1.0e-09

void do_fit(int natoms,real *w_rls,rvec *xp,rvec *x)
{
  int    c,r,n,j,m,i,irot;
  double omega[7][7],om[7][7],d[7],xnr,xpc;
  matrix vh,vk,R,u;
  /*
  matrix vh,vk,R,vh_d,vk_d,u;
  */
  real   mn;
  int    index;
  real   max_d;
  rvec   x_old;

  for(i=0;(i<7);i++) {
    d[i]=0;
    for(j=0;(j<7);j++) {
      omega[i][j]=0;
      om[i][j]=0;
    }
  }

  /*calculate the matrix U*/
  clear_mat(u);
  for(n=0;(n<natoms);n++) {
    if ((mn = w_rls[n]) != 0.0) {
      for(c=0; (c<DIM); c++) {
	xpc=xp[n][c];
	for(r=0; (r<DIM); r++) {
	  xnr=x[n][r];
	  u[c][r]+=mn*xnr*xpc;
	}
      }
    }
  }
  
  /*construct omega*/
  /*omega is symmetric -> omega==omega' */
  for(r=0;(r<6);r++)
    for(c=0;(c<=r);c++)
      if ((r>=3) && (c<3)) {
        omega[r+1][c+1]=u[r-3][c];
        omega[c+1][r+1]=u[r-3][c];
      }
      else {
        omega[r+1][c+1]=0;
        omega[c+1][r+1]=0;
      }

  /*determine h and k*/
  jacobi(omega,6,d,om,&irot);
  /*real   **omega = input matrix a[1..n][1..n] must be symmetric
   *int     natoms = number of rows and columns
   *real      NULL = d[1]..d[n] are the eigenvalues of a[][]
   *real       **v = v[1..n][1..n] contains the vectors in columns
   *int      *irot = number of jacobi rotations
   */

  if (irot==0) {
    fprintf(stderr,"IROT=0\n");
  }

  index=0; /* For the compiler only */

  /* Copy only the first two eigenvectors */  
  for(j=0;(j<2);j++) {
    max_d=-1000;
    for(i=0;(i<6);i++)
      if (d[i+1]>max_d) {
        max_d=d[i+1];
        index=i;
      }
    d[index+1]=-10000;
    for(i=0;(i<3);i++) {
      vh[j][i]=M_SQRT2*om[i+1][index+1];
      vk[j][i]=M_SQRT2*om[i+4][index+1];
    }
  }
  /* Calculate the last eigenvector as the outer-product of the first two.
   * This insures that the conformation is not mirrored and
   * prevents problems with completely flat reference structures.
   */  
  oprod(vh[0],vh[1],vh[2]);
  oprod(vk[0],vk[1],vk[2]);

  /*determine R*/
  for(c=0;(c<3);c++)
    for(r=0;(r<3);r++)
      R[c][r]=vk[0][r]*vh[0][c]+
	      vk[1][r]*vh[1][c]+
	      vk[2][r]*vh[2][c];

  /*rotate X*/
  for(j=0;(j<natoms);j++) {
    for(m=0;(m<3);m++)
      x_old[m]=x[j][m];
    for(r=0;(r<3);r++) {
      x[j][r]=0;
      for(c=0;(c<3);c++)
        x[j][r]+=R[c][r]*x_old[c];
    }
  }
}

void reset_x(int ncm,atom_id ind_cm[],
	     int nrms,atom_id ind_rms[],
	     rvec x[],real mass[])
{
  int  i,m,ai;
  rvec xcm;
  real tm,mm;
  
  tm=0.0;
  clear_rvec(xcm);
  for(i=0; (i<ncm); i++) {
    ai=ind_cm[i];
    mm=mass[ai];
    for(m=0; (m<DIM); m++)
      xcm[m]+=mm*x[ai][m];
    tm+=mm;
  }
  for(m=0; (m<DIM); m++)
    xcm[m]/=tm;
    
  for(i=0; (i<nrms); i++)
    rvec_dec(x[ind_rms[i]],xcm);
}

