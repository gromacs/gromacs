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
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * GROningen MAchine for Chemical Simulation
 */
static char *SRCID_pme_c = "$Id$";

#include <stdio.h>
#include <math.h>
#include "typedefs.h"
#include "vec.h"
#include "complex.h"
#include "smalloc.h"
#include "futil.h"
#include "lrutil.h"
#include "fftgrid.h"
#include "fatal.h"
#include "physics.h"

#define TOL 2e-5


real calc_ewaldcoeff(real rc,real dtol)
{
  real x=5,low,high;
  int n,i=0;
  
  
  do {
    i++;
    x*=2;
  } while((erfc(x*rc)/rc)>dtol);

  n=i+60; /* search tolerance is 2^-60 */
  low=0;
  high=x;
  for(i=0;i<n;i++) {
    x=(low+high)/2;
    if((erfc(x*rc)/rc)>dtol)
      low=x;
    else 
      high=x;
  }

  return x;
}

  

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
		  real charge[],   rvec box,
	      t_commrec *cr,	      t_nsborder *nsb,
	      matrix lrvir, real ewaldcoeff)
{
  static    bool bFirst = TRUE;
  static    int       nx,ny,nz,kmax;
  static    cvec      **eir;
  static    t_complex  *tab_xy,*tab_qxyz;
  real factor=-1.0/(4*ewaldcoeff*ewaldcoeff);
  real energy;
  rvec lll;
  int  lowiy,lowiz,ix,iy,iz,n;
  real tmp,cs,ss,ak,akv,mx,my,mz,m2;
  
  if (bFirst) {
      if (bVerbose)
    fprintf(log,"Will do ordinary reciprocal space Ewald sum.\n");

    if (cr != NULL) {
      if (cr->nprocs > 1)
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
  
  lowiy=0;
  lowiz=1;
  energy=0;
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
	energy+=ak*(cs*cs+ss*ss);
	tmp=akv*(cs*cs+ss*ss);	       
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
  tmp=4.0*M_PI/(box[XX]*box[YY]*box[ZZ])*ONE_4PI_EPS0;
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
