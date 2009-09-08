/*
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

#include "typedefs.h"
#include "smalloc.h"
#include "update.h"
#include "vec.h"
#include "macros.h"
#include "physics.h"
#include "names.h"
#include "gmx_fatal.h"
#include "txtdump.h"
#include "nrnb.h"
#include "gmx_random.h"

/* 
 * This file implements temperature and pressure coupling algorithms:
 * For now only the Weak coupling and the modified weak coupling.
 *
 * Furthermore computation of pressure and temperature is done here
 *
 */

real calc_pres(int ePBC,int nwall,matrix box,tensor ekin,tensor vir,
	       tensor pres,real Elr)
{
  int  n,m;
  real fac,Plr;

  if (ePBC==epbcNONE || (ePBC==epbcXY && nwall!=2))
    clear_mat(pres);
  else {
    /* Uitzoeken welke ekin hier van toepassing is, zie Evans & Morris - E. 
     * Wrs. moet de druktensor gecorrigeerd worden voor de netto stroom in  
     * het systeem...       
     */
    
    /* Long range correction for periodic systems, see
     * Neumann et al. JCP
     * divide by 6 because it is multiplied by fac later on.
     * If Elr = 0, no correction is made.
     */

    /* This formula should not be used with Ewald or PME, 
     * where the full long-range virial is calculated. EL 990823
     */
    Plr = Elr/6.0;
    
    fac=PRESFAC*2.0/det(box);
    for(n=0; (n<DIM); n++)
      for(m=0; (m<DIM); m++)
	pres[n][m]=(ekin[n][m]-vir[n][m]+Plr)*fac;
	
    if (debug) {
      pr_rvecs(debug,0,"PC: pres",pres,DIM);
      pr_rvecs(debug,0,"PC: ekin",ekin,DIM);
      pr_rvecs(debug,0,"PC: vir ",vir, DIM);
      pr_rvecs(debug,0,"PC: box ",box, DIM);
    }
  }
  return trace(pres)/3.0;
}

real calc_temp(real ekin,real nrdf)
{
  if (nrdf > 0)
    return (2.0*ekin)/(nrdf*BOLTZ);
  else
    return 0;
}

void parrinellorahman_pcoupl(FILE *fplog,gmx_large_int_t step,
			     t_inputrec *ir,real dt,tensor pres,
			     tensor box,tensor box_rel,tensor boxv,
			     tensor M,matrix mu,bool bFirstStep)
{
  /* This doesn't do any coordinate updating. It just
   * integrates the box vector equations from the calculated
   * acceleration due to pressure difference. We also compute
   * the tensor M which is used in update to couple the particle
   * coordinates to the box vectors.
   *
   * In Nose and Klein (Mol.Phys 50 (1983) no 5., p 1055) this is
   * given as
   *            -1    .           .     -1
   * M_nk = (h')   * (h' * h + h' h) * h
   *
   * with the dots denoting time derivatives and h is the transformation from
   * the scaled frame to the real frame, i.e. the TRANSPOSE of the box. 
   * This also goes for the pressure and M tensors - they are transposed relative
   * to ours. Our equation thus becomes:
   *
   *                  -1       .    .           -1
   * M_gmx = M_nk' = b  * (b * b' + b * b') * b'
   * 
   * where b is the gromacs box matrix.                       
   * Our box accelerations are given by
   *   ..                                    ..
   *   b = vol/W inv(box') * (P-ref_P)     (=h')
   */
  
  int    d,n;
  tensor winv;
  real   vol=box[XX][XX]*box[YY][YY]*box[ZZ][ZZ];
  real   atot,arel,change,maxchange,xy_pressure;
  tensor invbox,pdiff,t1,t2;

  real maxl;

  m_inv_ur0(box,invbox);

  if (!bFirstStep) {
    /* Note that PRESFAC does not occur here.
     * The pressure and compressibility always occur as a product,
     * therefore the pressure unit drops out.
     */
    maxl=max(box[XX][XX],box[YY][YY]);
    maxl=max(maxl,box[ZZ][ZZ]);
    for(d=0;d<DIM;d++)
      for(n=0;n<DIM;n++)
	winv[d][n]=
	  (4*M_PI*M_PI*ir->compress[d][n])/(3*ir->tau_p*ir->tau_p*maxl);
    
    m_sub(pres,ir->ref_p,pdiff);
    
    if(ir->epct==epctSURFACETENSION) {
      /* Unlike Berendsen coupling it might not be trivial to include a z
       * pressure correction here? On the other hand we don't scale the
       * box momentarily, but change accelerations, so it might not be crucial.
       */
      xy_pressure=0.5*(pres[XX][XX]+pres[YY][YY]);
      for(d=0;d<ZZ;d++)
	pdiff[d][d]=(xy_pressure-(pres[ZZ][ZZ]-ir->ref_p[d][d]/box[d][d]));
    }
    
    tmmul(invbox,pdiff,t1);
    /* Move the off-diagonal elements of the 'force' to one side to ensure
     * that we obey the box constraints.
     */
    for(d=0;d<DIM;d++) {
      for(n=0;n<d;n++) {
	t1[d][n] += t1[n][d];
	t1[n][d] = 0;
      }
    }
    
    switch (ir->epct) {
    case epctANISOTROPIC:
      for(d=0;d<DIM;d++) 
	for(n=0;n<=d;n++)
	  t1[d][n] *= winv[d][n]*vol;
      break;
    case epctISOTROPIC:
      /* calculate total volume acceleration */
      atot=box[XX][XX]*box[YY][YY]*t1[ZZ][ZZ]+
	box[XX][XX]*t1[YY][YY]*box[ZZ][ZZ]+
	t1[XX][XX]*box[YY][YY]*box[ZZ][ZZ];
      arel=atot/(3*vol);
      /* set all RELATIVE box accelerations equal, and maintain total V
       * change speed */
      for(d=0;d<DIM;d++)
	for(n=0;n<=d;n++)
	  t1[d][n] = winv[0][0]*vol*arel*box[d][n];    
      break;
    case epctSEMIISOTROPIC:
    case epctSURFACETENSION:
      /* Note the correction to pdiff above for surftens. coupling  */
      
      /* calculate total XY volume acceleration */
      atot=box[XX][XX]*t1[YY][YY]+t1[XX][XX]*box[YY][YY];
      arel=atot/(2*box[XX][XX]*box[YY][YY]);
      /* set RELATIVE XY box accelerations equal, and maintain total V
       * change speed. Dont change the third box vector accelerations */
      for(d=0;d<ZZ;d++)
	for(n=0;n<=d;n++)
	  t1[d][n] = winv[d][n]*vol*arel*box[d][n];
      for(n=0;n<DIM;n++)
	t1[ZZ][n] *= winv[d][n]*vol;
      break;
    default:
      gmx_fatal(FARGS,"Parrinello-Rahman pressure coupling type %s "
		  "not supported yet\n",EPCOUPLTYPETYPE(ir->epct));
      break;
    }
    
    maxchange=0;
    for(d=0;d<DIM;d++)
      for(n=0;n<=d;n++) {
	boxv[d][n] += dt*t1[d][n];
	/* We do NOT update the box vectors themselves here, since
	 * we need them for shifting later. It is instead done last
	 * in the update() routine.
	 */
	
	/* Calculate the change relative to diagonal elements -
	 * since it's perfectly ok for the off-diagonal ones to
	 * be zero it doesn't make sense to check the change relative
	 * to its current size.
	 */
	change=fabs(dt*boxv[d][n]/box[d][d]);
	if(change>maxchange)
	  maxchange=change;
      }
    
    if (maxchange > 0.01 && fplog) {
      char buf[22];
      fprintf(fplog,"\nStep %s  Warning: Pressure scaling more than 1%%.\n",
	      gmx_step_str(step,buf));
    }
  }
  
  preserve_box_shape(ir,box_rel,boxv);

  mtmul(boxv,box,t1);       /* t1=boxv * b' */
  mmul(invbox,t1,t2);
  mtmul(t2,invbox,M);

  /* Determine the scaling matrix mu for the coordinates */
  for(d=0;d<DIM;d++)
    for(n=0;n<=d;n++)
      t1[d][n] = box[d][n] + dt*boxv[d][n];
  preserve_box_shape(ir,box_rel,t1);
  /* t1 is the box at t+dt, determine mu as the relative change */
  mmul_ur0(invbox,t1,mu);
}

void berendsen_pcoupl(FILE *fplog,gmx_large_int_t step,
		      t_inputrec *ir,real dt,tensor pres,matrix box,
		      matrix mu)
{
  int    d,n;
  real   scalar_pressure, xy_pressure, p_corr_z;
  char   *ptr,buf[STRLEN];

  /*
   *  Calculate the scaling matrix mu
   */
  scalar_pressure=0;
  xy_pressure=0;
  for(d=0; d<DIM; d++) {
    scalar_pressure += pres[d][d]/DIM;
    if (d != ZZ)
      xy_pressure += pres[d][d]/(DIM-1);
  }
  /* Pressure is now in bar, everywhere. */
#define factor(d,m) (ir->compress[d][m]*dt/ir->tau_p)
  
  /* mu has been changed from pow(1+...,1/3) to 1+.../3, since this is
   * necessary for triclinic scaling
   */
  clear_mat(mu);
  switch (ir->epct) {
  case epctISOTROPIC:
    for(d=0; d<DIM; d++)
      mu[d][d] = 1.0 - factor(d,d)*(ir->ref_p[d][d] - scalar_pressure)/DIM;
    break;
  case epctSEMIISOTROPIC:
    for(d=0; d<ZZ; d++)
      mu[d][d] = 1.0 - factor(d,d)*(ir->ref_p[d][d]-xy_pressure)/DIM;
    mu[ZZ][ZZ] = 
      1.0 - factor(ZZ,ZZ)*(ir->ref_p[ZZ][ZZ] - pres[ZZ][ZZ])/DIM;
    break;
  case epctANISOTROPIC:
    for(d=0; d<DIM; d++)
      for(n=0; n<DIM; n++)
	mu[d][n] = (d==n ? 1.0 : 0.0) 
	  -factor(d,n)*(ir->ref_p[d][n] - pres[d][n])/DIM;
    break;
  case epctSURFACETENSION:
    /* ir->ref_p[0/1] is the reference surface-tension times *
     * the number of surfaces                                */
    if (ir->compress[ZZ][ZZ])
      p_corr_z = dt/ir->tau_p*(ir->ref_p[ZZ][ZZ] - pres[ZZ][ZZ]);
    else
      /* when the compressibity is zero, set the pressure correction   *
       * in the z-direction to zero to get the correct surface tension */
      p_corr_z = 0;
    mu[ZZ][ZZ] = 1.0 - ir->compress[ZZ][ZZ]*p_corr_z;
    for(d=0; d<DIM-1; d++)
      mu[d][d] = 1.0 + factor(d,d)*(ir->ref_p[d][d]/(mu[ZZ][ZZ]*box[ZZ][ZZ])
				    - (pres[ZZ][ZZ]+p_corr_z - xy_pressure))/(DIM-1);
    break;
  default:
    gmx_fatal(FARGS,"Berendsen pressure coupling type %s not supported yet\n",
		EPCOUPLTYPETYPE(ir->epct));
    break;
  }
  /* To fullfill the orientation restrictions on triclinic boxes
   * we will set mu_yx, mu_zx and mu_zy to 0 and correct
   * the other elements of mu to first order.
   */
  mu[YY][XX] += mu[XX][YY];
  mu[ZZ][XX] += mu[XX][ZZ];
  mu[ZZ][YY] += mu[YY][ZZ];
  mu[XX][YY] = 0;
  mu[XX][ZZ] = 0;
  mu[YY][ZZ] = 0;

  if (debug) {
    pr_rvecs(debug,0,"PC: pres ",pres,3);
    pr_rvecs(debug,0,"PC: mu   ",mu,3);
  }
  
  if (mu[XX][XX]<0.99 || mu[XX][XX]>1.01 ||
      mu[YY][YY]<0.99 || mu[YY][YY]>1.01 ||
      mu[ZZ][ZZ]<0.99 || mu[ZZ][ZZ]>1.01) {
    char buf2[22];
    sprintf(buf,"\nStep %s  Warning: pressure scaling more than 1%%, "
	    "mu: %g %g %g\n",
	    gmx_step_str(step,buf2),mu[XX][XX],mu[YY][YY],mu[ZZ][ZZ]);
    if (fplog)
      fprintf(fplog,"%s",buf);
    fprintf(stderr,"%s",buf);
  }
}

void berendsen_pscale(t_inputrec *ir,matrix mu,
		      matrix box,matrix box_rel,
		      int start,int nr_atoms,
		      rvec x[],unsigned short cFREEZE[],
		      t_nrnb *nrnb)
{
  ivec   *nFreeze=ir->opts.nFreeze;
  int    n,d,g=0;
      
  /* Scale the positions */
  for (n=start; n<start+nr_atoms; n++) {
    if (cFREEZE)
      g = cFREEZE[n];
    
    if (!nFreeze[g][XX])
      x[n][XX] = mu[XX][XX]*x[n][XX]+mu[YY][XX]*x[n][YY]+mu[ZZ][XX]*x[n][ZZ];
    if (!nFreeze[g][YY])
      x[n][YY] = mu[YY][YY]*x[n][YY]+mu[ZZ][YY]*x[n][ZZ];
    if (!nFreeze[g][ZZ])
      x[n][ZZ] = mu[ZZ][ZZ]*x[n][ZZ];
  }
  /* compute final boxlengths */
  for (d=0; d<DIM; d++) {
    box[d][XX] = mu[XX][XX]*box[d][XX]+mu[YY][XX]*box[d][YY]+mu[ZZ][XX]*box[d][ZZ];
    box[d][YY] = mu[YY][YY]*box[d][YY]+mu[ZZ][YY]*box[d][ZZ];
    box[d][ZZ] = mu[ZZ][ZZ]*box[d][ZZ];
  }      

  preserve_box_shape(ir,box_rel,box);
  
  /* (un)shifting should NOT be done after this,
   * since the box vectors might have changed
   */
  inc_nrnb(nrnb,eNR_PCOUPL,nr_atoms);
}

void berendsen_tcoupl(t_grpopts *opts,gmx_ekindata_t *ekind,real dt)
{
  int    i;
  
  real   T,reft=0,lll; 

  for(i=0; (i<opts->ngtc); i++) {
    T = ekind->tcstat[i].Th;
    
    if ((opts->tau_t[i] > 0) && (T > 0.0)) {
 
      reft = max(0.0,opts->ref_t[i]);
      lll  = sqrt(1.0 + (dt/opts->tau_t[i])*(reft/T-1.0));
      ekind->tcstat[i].lambda = max(min(lll,1.25),0.8);
    }
    else {
       ekind->tcstat[i].lambda = 1.0;
    }

    if (debug)
      fprintf(debug,"TC: group %d: T: %g, Lambda: %g\n",
	      i,T,ekind->tcstat[i].lambda);
  }
}

void nosehoover_tcoupl(t_grpopts *opts,gmx_ekindata_t *ekind,real dt,
		       real xi[],double ixi[])
{
  real  Qinv;
  int   i;
  real  reft,oldxi;

  for(i=0; (i<opts->ngtc); i++) {
    if ((opts->tau_t[i] > 0) && (opts->ref_t[i] > 0))
      Qinv=1.0/(opts->tau_t[i]*opts->tau_t[i]*opts->ref_t[i]/(4*M_PI*M_PI));
    else
      Qinv=0.0;
    reft = max(0.0,opts->ref_t[i]);
    oldxi = xi[i];
    xi[i]  += dt*Qinv*(ekind->tcstat[i].Th - reft);
    ixi[i] += dt*(oldxi + xi[i])*0.5;
  }
}

real nosehoover_energy(t_grpopts *opts,gmx_ekindata_t *ekind,
		       real *xi,double *ixi)
{
  int  i,nd;
  real ener_nh;
  
  ener_nh = 0;
  for(i=0; i<opts->ngtc; i++) {
    nd = opts->nrdf[i];
    if (nd > 0)
      ener_nh += (sqr(xi[i]*opts->tau_t[i]/(2*M_PI))*0.5 +
		  ixi[i])*nd*BOLTZ*max(opts->ref_t[i],0);
  }

  return ener_nh;
}

static real vrescale_gamdev(int ia, gmx_rng_t rng)
/* Gamma distribution, adapted from numerical recipes */
{
  int j;
  real am,e,s,v1,v2,x,y;
  
  if (ia < 6) {
    x = 1.0;
    for(j=1; j<=ia; j++) {
      x *= gmx_rng_uniform_real(rng);
    }
    x = -log(x);
  } else {
    do {
      do {
        do {
          v1 = gmx_rng_uniform_real(rng);
          v2 = 2.0*gmx_rng_uniform_real(rng)-1.0;
        } while (v1*v1 + v2*v2 > 1.0);
        y = v2/v1;
        am = ia - 1;
        s = sqrt(2.0*am + 1.0);
        x = s*y + am;
      } while (x <= 0.0);
      e = (1.0 + y*y)*exp(am*log(x/am) - s*y);
    } while (gmx_rng_uniform_real(rng) > e);
  }

  return x;
}

static real vrescale_sumnoises(int nn,gmx_rng_t rng)
{
/*
 * Returns the sum of n independent gaussian noises squared
 * (i.e. equivalent to summing the square of the return values
 * of nn calls to gmx_rng_gaussian_real).xs
 */
  real rr;

  if (nn == 0) {
    return 0.0;
  } else if (nn == 1) {
    rr = gmx_rng_gaussian_real(rng);
    return rr*rr;
  } else if (nn % 2 == 0) {
    return 2.0*vrescale_gamdev(nn/2,rng);
  } else {
    rr = gmx_rng_gaussian_real(rng);
    return 2.0*vrescale_gamdev((nn-1)/2,rng) + rr*rr;
  }
}

static real vrescale_resamplekin(real kk,real sigma, int ndeg, real taut,
				 gmx_rng_t rng)
{
/*
 * Generates a new value for the kinetic energy,
 * according to Bussi et al JCP (2007), Eq. (A7)
 * kk:    present value of the kinetic energy of the atoms to be thermalized (in arbitrary units)
 * sigma: target average value of the kinetic energy (ndeg k_b T/2)  (in the same units as kk)
 * ndeg:  number of degrees of freedom of the atoms to be thermalized
 * taut:  relaxation time of the thermostat, in units of 'how often this routine is called'
 */
  real factor,rr;

  if (taut > 0.1) {
    factor = exp(-1.0/taut);
  } else {
    factor = 0.0;
  }
  rr = gmx_rng_gaussian_real(rng);
  return
    kk +
    (1.0 - factor)*(sigma*(vrescale_sumnoises(ndeg-1,rng) + rr*rr)/ndeg - kk) +
    2.0*rr*sqrt(kk*sigma/ndeg*(1.0 - factor)*factor);
}

void vrescale_tcoupl(t_grpopts *opts,gmx_ekindata_t *ekind,real dt,
		     double therm_integral[],
		     gmx_rng_t rng)
{
  int    i;
  real   Ek,Ek_ref1,Ek_ref,Ek_new; 

  for(i=0; (i<opts->ngtc); i++) {
    Ek = trace(ekind->tcstat[i].ekinh);
    
    if (opts->tau_t[i] >= 0 && opts->nrdf[i] > 0 && Ek > 0) {
      Ek_ref1   = 0.5*opts->ref_t[i]*BOLTZ;
      Ek_ref    = Ek_ref1*opts->nrdf[i];

      Ek_new =
	vrescale_resamplekin(Ek,Ek_ref,opts->nrdf[i],opts->tau_t[i]/dt,rng);

      /* Analytically Ek_new>=0, but we check for rounding errors */
      if (Ek_new <= 0) {
	ekind->tcstat[i].lambda = 0.0;
      } else {
	ekind->tcstat[i].lambda = sqrt(Ek_new/Ek);
      }
      therm_integral[i] -= Ek_new - Ek;

      if (debug) {
	fprintf(debug,"TC: group %d: Ekr %g, Ek %g, Ek_new %g, Lambda: %g\n",
		i,Ek_ref,Ek,Ek_new,ekind->tcstat[i].lambda);
      }
    } else {
       ekind->tcstat[i].lambda = 1.0;
    }
  }
}

real vrescale_energy(t_grpopts *opts,double therm_integral[])
{
  int i;
  real ener;

  ener = 0;
  for(i=0; i<opts->ngtc; i++) {
    ener += therm_integral[i];
  }
  
  return ener;
}

/* set target temperatures if we are annealing */
void update_annealing_target_temp(t_grpopts *opts,real t)
{
  int i,j,n,npoints;
  real pert,thist=0,x;

  for(i=0;i<opts->ngtc;i++) {
    npoints = opts->anneal_npoints[i];
    switch (opts->annealing[i]) {
    case eannNO:
      continue;
    case  eannPERIODIC:
      /* calculate time modulo the period */
      pert  = opts->anneal_time[i][npoints-1];
      n     = t / pert;
      thist = t - n*pert; /* modulo time */
      /* Make sure rounding didn't get us outside the interval */
      if (fabs(thist-pert) < GMX_REAL_EPS*100)
	thist=0;
      break;
    case eannSINGLE:
      thist = t;
      break;
    default:
      gmx_fatal(FARGS,"Death horror in update_annealing_target_temp (i=%d/%d npoints=%d)",i,opts->ngtc,npoints);
    }
    /* We are doing annealing for this group if we got here, 
     * and we have the (relative) time as thist.
     * calculate target temp */
    j=0;
    while ((j < npoints-1) && (thist>(opts->anneal_time[i][j+1])))
      j++;
    if (j < npoints-1) {
      /* Found our position between points j and j+1. 
       * Interpolate: x is the amount from j+1, (1-x) from point j 
       * First treat possible jumps in temperature as a special case.
       */
      if ((opts->anneal_time[i][j+1]-opts->anneal_time[i][j]) < GMX_REAL_EPS*100)
	opts->ref_t[i]=opts->anneal_temp[i][j+1];
      else {
	x = ((thist-opts->anneal_time[i][j])/
	     (opts->anneal_time[i][j+1]-opts->anneal_time[i][j]));
	opts->ref_t[i] = x*opts->anneal_temp[i][j+1]+(1-x)*opts->anneal_temp[i][j];
      }
    }
    else {
      opts->ref_t[i] = opts->anneal_temp[i][npoints-1];
    }
  }
}
