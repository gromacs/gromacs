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

#include <stdio.h>
#include <math.h>
#include "typedefs.h"
#include "vec.h"
#include "gmxcomplex.h"
#include "smalloc.h"
#include "futil.h"
#include "ewald_util.h"
#include "shift_util.h"
#include "fftgrid.h"
#include "fatal.h"
#include "physics.h"

#define TOL 2e-5



  

/* the other routines are in complex.h */
static t_complex conjmul(t_complex a,t_complex b)
{
  t_complex c;
  
  c.re = a.re*b.re + a.im*b.im;
  c.im = a.im*b.re - a.re*b.im;
  
  return c;
}

  
	  

static void tabulate_eir(int natom,rvec x[],int kmax,cvec **eir,rvec lll)
{
  int  i,j,m;
  
  if (kmax < 1) {
      printf("Go away! kmax = %d\n",kmax);
      exit(1);
  }
  
  for(i=0; (i<natom); i++) {
    for(m=0; (m<3); m++) {
      eir[0][i][m].re = 1;
      eir[0][i][m].im = 0;
    }
    
    for(m=0; (m<3); m++) {
      eir[1][i][m].re = cos(x[i][m]*lll[m]);
      eir[1][i][m].im = sin(x[i][m]*lll[m]); 
    }
    for(j=2; (j<kmax); j++) 
      for(m=0; (m<3); m++)
	eir[j][i][m] = cmul(eir[j-1][i][m],eir[1][i][m]);
  }
}




real do_ewald(FILE *log,       bool bVerbose,
	      t_inputrec *ir,
	      rvec x[],        rvec f[],
	      real chargeA[],  real chargeB[],
	      rvec box,
	      t_commrec *cr,   t_nsborder *nsb,
	      matrix lrvir,    real ewaldcoeff,
	      real lambda,     real *dvdlambda)
{
  static    bool bFirst = TRUE;
  static    int       nx,ny,nz,kmax;
  static    cvec      **eir;
  static    t_complex  *tab_xy,*tab_qxyz;
  real factor=-1.0/(4*ewaldcoeff*ewaldcoeff);
  real *charge,energy_AB[2],energy;
  rvec lll;
  int  lowiy,lowiz,ix,iy,iz,n,q;
  real tmp,cs,ss,ak,akv,mx,my,mz,m2,scale;
  bool bFreeEnergy;

  bFreeEnergy = (ir->efep != efepNO);

  if (bFirst) {
      if (bVerbose)
    fprintf(log,"Will do ordinary reciprocal space Ewald sum.\n");

    if (cr != NULL) {
      if (cr->nnodes > 1 || cr->nthreads>1)
	fatal_error(0,"No parallel Ewald. Use PME instead.\n");
    }
    
    nx = ir->nkx+1;
    ny = ir->nky+1;
    nz = ir->nkz+1;
    kmax = max(nx,max(ny,nz));
    snew(eir,kmax);
    for(n=0;n<kmax;n++)
      snew(eir[n],HOMENR(nsb));
    snew(tab_xy,HOMENR(nsb));
    snew(tab_qxyz,HOMENR(nsb));
    bFirst = FALSE;
  }
  clear_mat(lrvir);
  
  calc_lll(box,lll);
  /* make tables for the structure factor parts */
  tabulate_eir(HOMENR(nsb),x,kmax,eir,lll);

  for(q=0; q<(bFreeEnergy ? 2 : 1); q++) {
    if (!bFreeEnergy) {
      charge = chargeA;
      scale = 1.0;
    } else if (q==0) {
      charge = chargeA;
      scale = 1.0 - lambda;
    } else {
      charge = chargeB;
      scale = lambda;
    }
    lowiy=0;
    lowiz=1;
    energy_AB[q]=0;
    for(ix=0;ix<nx;ix++) {
      mx=ix*lll[XX];
      for(iy=lowiy;iy<ny;iy++) {
	my=iy*lll[YY];
	if(iy>=0) 
	  for(n=0;n<HOMENR(nsb);n++) 
	    tab_xy[n]=cmul(eir[ix][n][XX],eir[iy][n][YY]);
	else 
	  for(n=0;n<HOMENR(nsb);n++) 
	    tab_xy[n]=conjmul(eir[ix][n][XX],eir[-iy][n][YY]); 
	for(iz=lowiz;iz<nz;iz++) {
	  mz=iz*lll[ZZ];	       
	  m2=mx*mx+my*my+mz*mz;
	  ak=exp(m2*factor)/m2;
	  akv=2.0*ak*(1.0/m2-factor);  
	  if(iz>=0) 
	    for(n=0;n<HOMENR(nsb);n++) 
	      tab_qxyz[n]=rcmul(charge[n],cmul(tab_xy[n],eir[iz][n][ZZ]));
	  else 
	    for(n=0;n<HOMENR(nsb);n++) 
	      tab_qxyz[n]=rcmul(charge[n],conjmul(tab_xy[n],eir[-iz][n][ZZ]));
	  
	  cs=ss=0;
	  for(n=0;n<HOMENR(nsb);n++) {
	    cs+=tab_qxyz[n].re;
	    ss+=tab_qxyz[n].im;
	  }
	  energy_AB[q]+=ak*(cs*cs+ss*ss);
	  tmp=scale*akv*(cs*cs+ss*ss);	       
	  lrvir[XX][XX]-=tmp*mx*mx;
	  lrvir[XX][YY]-=tmp*mx*my;
	  lrvir[XX][ZZ]-=tmp*mx*mz;
	  lrvir[YY][YY]-=tmp*my*my;
	  lrvir[YY][ZZ]-=tmp*my*mz;
	  lrvir[ZZ][ZZ]-=tmp*mz*mz;
	  for(n=0;n<HOMENR(nsb);n++) {
	    tmp=ak*(cs*tab_qxyz[n].im-ss*tab_qxyz[n].re);
	    f[n][XX]+=tmp*mx;
	    f[n][YY]+=tmp*my;
	    f[n][ZZ]+=tmp*mz;
	  }
	  lowiz=1-nz;
	}
	lowiy=1-ny;
      }
    }
  }
  
  tmp=4.0*M_PI/(box[XX]*box[YY]*box[ZZ])*ONE_4PI_EPS0;

  if (!bFreeEnergy) {
    energy = energy_AB[0];
  } else {
    energy = (1.0 - lambda)*energy_AB[0] + lambda*energy_AB[1];
    *dvdlambda += tmp*(energy_AB[1] - energy_AB[0]);
  }
  for(n=0;n<HOMENR(nsb);n++) {
    f[n][XX]*=2*tmp;
    f[n][YY]*=2*tmp;
    f[n][ZZ]*=2*tmp;
  }
  lrvir[XX][XX]=-0.5*tmp*(lrvir[XX][XX]+energy);
  lrvir[XX][YY]=-0.5*tmp*(lrvir[XX][YY]);
  lrvir[XX][ZZ]=-0.5*tmp*(lrvir[XX][ZZ]);
  lrvir[YY][YY]=-0.5*tmp*(lrvir[YY][YY]+energy);
  lrvir[YY][ZZ]=-0.5*tmp*(lrvir[YY][ZZ]);
  lrvir[ZZ][ZZ]=-0.5*tmp*(lrvir[ZZ][ZZ]+energy);
  
  lrvir[YY][XX]=lrvir[XX][YY];
  lrvir[ZZ][XX]=lrvir[XX][ZZ];
  lrvir[ZZ][YY]=lrvir[YY][ZZ];
  
  energy*=tmp;
  
  return energy;
}
