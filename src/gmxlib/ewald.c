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
static char *SRCID_ewald_c = "$Id$";

#include <stdio.h>
#include <math.h>
#include "typedefs.h"
#include "vec.h"
#include "complex.h"
#include "futil.h"
#include "lrutil.h"
#include "fftgrid.h"
#include "pppm.h"
#include "fatal.h"
/*#include "ghat.h"*/

#define TOL 2e-5

static t_complex Qk(rvec k,int nj,rvec rj[],real qj[])
/* Return fourier component of discrete charge distribution for all particles
 * (Eq. 33)
 */
{
  int       i;
  real      ipr;
  t_complex c,cqi,cpk;
  
  c = cnul;
  for(i=0; (i<nj); i++) {
    if (qj[i] != 0.0) {
      ipr = iprod(k,rj[i]);
      cpk = rcexp(-ipr);
      cqi = rcmul(qj[i],cpk);
      c   = cadd(c,cqi);
    }
  }
  /*if (fabs(c.re) < TOL) 
    c.re = 0;
    if (fabs(c.im) < TOL)
    c.im = 0;*/
  
  return c;
}

static void mk_qtab(int nx,int ny,int nz,t_complex ***qtab,
		    rvec box,int nj,rvec rj[],real charge[])
{
  int       ix,iy,iz;
  int       m;
  rvec      k,lll;
  t_complex c;
  
  calc_lll(box,lll);
  
  /* Loop over lattice vectors in fourier space */    
  for(ix=0; (ix < nx); ix++) {
    for(iy=0; (iy < ny); iy++) {
      for(iz=0; (iz < nz); iz++) {
	if ((ix == 0) && (iy == 0) && (iz == 0)) 
	  c = cnul;
	else {
	  calc_k(lll,ix,iy,iz,nx,ny,nz,k);
	  
	  c = Qk(k,nj,rj,charge);
	}
	qtab[ix][iy][iz] = c;
      }
    }
  }
}

static t_complex phi_k(rvec k,rvec rj,real gk,t_complex qk)
{
  real      phikr;
  t_complex cp,cpkr;
  
  phikr  = iprod(k,rj);
  cpkr   = rcexp(phikr);
  cp     = rcmul(gk,cmul(qk,cpkr));
  
  return cp;
}

void print_qkgrid(FILE *out,char *title,t_complex ***grid,real factor,
		  char *pdb,rvec box,int nx,int ny,int nz)
{
  static char *pdbformat="%-6s%5d  %-4.4s%3.3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n";
  FILE     *fp;
  int      i,ix,iy,iz;
  real     fac=50.0;
  rvec     boxfac;
  t_complex g;
  
  if (pdb)
    fp = ffopen(pdb,"w");
  else
    fp = out;
  if (!fp)
    return;

  boxfac[XX] = fac*box[XX]/nx;
  boxfac[YY] = fac*box[YY]/ny;
  boxfac[ZZ] = fac*box[ZZ]/nz;
  
  fprintf(fp,"Printing all non-zero complex elements of %s\n",title);
  for(i=ix=0; (ix<nx); ix++)
    for(iy=0; (iy<ny); iy++)
      for(iz=0; (iz<nz); iz++,i++) {
	g = grid[ix][iy][iz];
	if (fabs(g.im) > TOL) {
	  if (pdb) {
	    fprintf(fp,pdbformat,"ATOM",i,"H","H",' ',
		    i,ix*boxfac[XX],iy*boxfac[YY],iz*boxfac[ZZ],
		    1.0,g.im);
	  } 
	  else {
	    fprintf(fp,"%s[%2d][%2d][%2d] = %12.5e + i %12.5e%s\n",
		    title,ix,iy,iz,g.re*factor,g.im*factor,
		    (g.im != 0) ? " XXX" : "");
	  }
	}
      }
  fflush(fp);
}

real do_ewald(FILE *log,       t_inputrec *ir,
	      int natoms,      rvec x[],rvec f[],
	      real charge[],   rvec box,
	      real phi[],      t_commrec *cr)
{
  static    bool bFirst = TRUE;
  static    t_complex ***qtab;
  static    real      ***gtab;
  static    int       nx,ny,nz;
  
  real      Gk,energy;
  int       m,ix,iy,iz,i;
  rvec      k,lll;
  t_complex phik,phitot,Ei[3],Qk;

  if (bFirst) {  
    if (cr != NULL) {
      if (cr->nprocs > 1)
	fatal_error(0,"No parallel Ewald (yet)");
    }
  
    nx = ir->nkx;
    ny = ir->nky;
    nz = ir->nkz;
  
    fprintf(log,"Will do Ewald summation using %d x %d x %d k-vectors\n",
	    nx,ny,nz);
    fprintf(log,"This code is not optimized, so be prepared to wait!\n"); 
       
    gtab = mk_rgrid(nx,ny,nz);
    qtab = mk_cgrid(nx,ny,nz);

    mk_ghat(NULL,nx,ny,nz,gtab,box,ir->rshort,ir->rlong,TRUE);
    
    bFirst = FALSE;
  }
  
  /* Tabulate charge spread & charge distribution in fourier space */
  mk_qtab(nx,ny,nz,qtab,box,natoms,x,charge);
  
  print_qkgrid(log,"qk",qtab,1.0,"qkfour.pdb",box,nx,ny,nz);
    
  /* Now start the actual solution! */
  calc_lll(box,lll);
  
  energy = 0;
  for(i=0; (i<natoms); i++) {
    /* Loop over atoms */
    phitot = cnul;
    
    /* Electric field vector */
    for(m=0; (m<DIM); m++)
      Ei[m] = cnul;
      
    /* Loop over lattice vectors in fourier space */    
    for(ix=0; (ix < nx); ix++) {
      for(iy=0; (iy < ny); iy++) {
	for(iz=0; (iz < nz); iz++) {
	  /* exclude k = (0,0,0) */
	  if ((ix == 0) && (iy == 0) && (iz == 0))
	    continue;	  
	  
	  calc_k(lll,ix,iy,iz,nx,ny,nz,k);

	  Qk     = qtab[ix][iy][iz];
	  Gk     = gtab[ix][iy][iz];
	  
	  /* Calculate potential */
	  phik   = phi_k(k,x[i],Gk,Qk);
	  phitot = cadd(phik,phitot);
	  
	  for(m=0; (m<DIM); m++)
	    Ei[m] = cadd(Ei[m],rcmul(k[m],phik));
	}
      }
    }
    /* Potential at atom i, and energy contribution */
#ifdef DEBUG
    fprintf(log,"phi[%3d] = %10.3e + i %10.3e\n",i,phitot.re,phitot.im);
#endif
    phi[i]  = phitot.re;
    energy += phi[i]*charge[i];
    
    /* Force on the atom */
    for(m=0; (m<DIM); m++) {
      f[i][m] = charge[i]*Ei[m].im;
#ifdef DEBUG
      if (fabs(Ei[m].re) > 1e-6)
	fprintf(log,"Ei[%1d] atom %3d = %10.3e + i%10.3e\n",
		m,i,Ei[m].re,Ei[m].im);
#endif
    }
  }
  return energy;
}

