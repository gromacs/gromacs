/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_coupling_c = "$Id$";

#include "typedefs.h"
#include "smalloc.h"
#include "update.h"
#include "vec.h"
#include "macros.h"
#include "physics.h"
#include "names.h"
#include "fatal.h"
#include "txtdump.h"
#include "nrnb.h"

/* 
 * This file implements temperature and pressure coupling algorithms:
 * For now only the Weak coupling and the modified weak coupling.
 *
 * Furthermore computation of pressure and temperature is done here
 *
 */

void calc_pres(int ePBC,matrix box,tensor ekin,tensor vir,tensor pres,real Elr)
{
  int  n,m;
  real fac,Plr;

  if (ePBC == epbcNONE)
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
}

real calc_temp(real ekin,real nrdf)
{
  return (2.0*ekin)/(nrdf*BOLTZ);
}

void do_pcoupl(t_inputrec *ir,int step,tensor pres,
	       matrix box,int start,int nr_atoms,
	       rvec x[],unsigned short cFREEZE[],
	       t_nrnb *nrnb,ivec nFreeze[])
{
  int    n,d,g;
  real   scalar_pressure, xy_pressure, p_corr_z;
  matrix mu;
  char   *ptr,buf[STRLEN];

  /*
   *  PRESSURE SCALING 
   *  Step (2P)
   */

#define factor(d,m) (ir->compress[d][m]*ir->delta_t/ir->tau_p)

  scalar_pressure=0;
  xy_pressure=0;
  for(d=0; d<DIM; d++) {
    scalar_pressure += pres[d][d]/DIM;
    if (d != ZZ)
      xy_pressure += pres[d][d]/(DIM-1);
  }
  
  /* Pressure is now in bar, everywhere. */
  /* mu has been changed from pow(1+...,1/3) to 1+.../3, since this is
   * necessary for triclinic scaling
   */
  if (scalar_pressure != 0.0) {
    clear_mat(mu);
    switch (ir->epc) {
    case epcNO:
      /* do_pcoupl should not be called in this case to save some work */
      for(d=0; d<DIM; d++)
	mu[d][d] = 1.0;
      break;
    case epcISOTROPIC:
      for(d=0; d<DIM; d++)
	mu[d][d] = 1.0 - factor(d,d)*(ir->ref_p[d][d] - scalar_pressure)/DIM;
      break;
    case epcSEMIISOTROPIC:
      for(d=0; d<ZZ; d++)
	mu[d][d] = 1.0 - factor(d,d)*(ir->ref_p[d][d]-xy_pressure)/DIM;
      mu[ZZ][ZZ] = 
	1.0 - factor(d,d)*(ir->ref_p[ZZ][ZZ] - pres[ZZ][ZZ])/DIM;
      break;
    case epcANISOTROPIC:
      for(d=0; d<DIM; d++)
	for(n=0; n<DIM; n++)
	  mu[d][n] = (d==n ? 1.0 : 0.0) 
	    -factor(d,n)*(ir->ref_p[d][n] - pres[d][n])/DIM;
      break;
    case epcSURFACETENSION:
      /* ir->ref_p[0/1] is the reference surface-tension times *
       * the number of surfaces                                */
      if (ir->compress[ZZ])
	p_corr_z = ir->delta_t/ir->tau_p*(ir->ref_p[ZZ][ZZ] - pres[ZZ][ZZ]);
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
      fatal_error(0,"Pressure coupling type %s not supported yet\n",
		  EPCOUPLTYPE(ir->epc));
    }
    /* To fullfill the orientation restrictions on triclinic boxes
     * set mu_yx, mu_zx and mu_zy to 0 and correct the other elements
     * of mu to first order.
     */
    mu[XX][YY] += mu[YY][XX];
    mu[XX][ZZ] += mu[ZZ][XX];
    mu[YY][ZZ] += mu[ZZ][YY];
    
    if (debug) {
      pr_rvecs(debug,0,"PC: pres ",pres,3);
      pr_rvecs(debug,0,"PC: mu   ",mu,3);
    }
    
    if (mu[XX][XX]<0.99 || mu[XX][XX]>1.01 ||
	mu[YY][YY]<0.99 || mu[YY][YY]>1.01 ||
	mu[ZZ][ZZ]<0.99 || mu[ZZ][ZZ]>1.01) {
      sprintf(buf,"\nStep %d  Warning: pressure scaling more than 1%%, "
	      "mu: %g %g %g\n",step,mu[XX][XX],mu[YY][YY],mu[ZZ][ZZ]);
      fprintf(stdlog,"%s",buf);
      fprintf(stderr,"%s",buf);
    }

    /* Scale the positions */
    for (n=start; n<start+nr_atoms; n++) {
      g=cFREEZE[n];
      
      if (!nFreeze[g][XX])
	x[n][XX] = mu[XX][XX]*x[n][XX]+mu[XX][YY]*x[n][YY]+mu[XX][ZZ]*x[n][ZZ];
      if (!nFreeze[g][YY])
	x[n][YY] = mu[YY][YY]*x[n][YY]+mu[YY][ZZ]*x[n][ZZ];
      if (!nFreeze[g][ZZ])
	x[n][ZZ] = mu[ZZ][ZZ]*x[n][ZZ];
    }
    /* compute final boxlengths */
    for (d=0; d<DIM; d++) {
      box[d][XX] = mu[XX][XX]*box[d][XX]+mu[XX][YY]*box[d][YY]
	+mu[XX][ZZ]*box[n][ZZ];
      box[d][YY] = mu[YY][YY]*box[d][YY]+mu[YY][ZZ]*box[d][ZZ];
      box[d][ZZ] = mu[ZZ][ZZ]*box[d][ZZ];
    }

    /* check if the box still obeys the restrictions, if not, correct it */
    if (box[ZZ][YY] > BOX_MARGIN*box[YY][YY]) {
      fprintf(stdlog,0,"Correcting invalid box:\n");
      pr_rvecs(stdlog,0,"old box",box,DIM);
      rvec_dec(box[ZZ],box[YY]);
      pr_rvecs(stdlog,0,"new box",box,DIM);
    } else if (-box[ZZ][YY] > BOX_MARGIN*box[YY][YY]) {
      fprintf(stdlog,0,"Correcting invalid box:\n");
      pr_rvecs(stdlog,0,"old box",box,DIM);
      rvec_inc(box[ZZ],box[YY]);
      pr_rvecs(stdlog,0,"new box",box,DIM);
    }
    if (fabs(box[YY][XX])+fabs(box[ZZ][XX]) > BOX_MARGIN*box[XX][XX]) {
      fprintf(stdlog,0,"Correcting invalid box:\n");
      pr_rvecs(stdlog,0,"old box",box,DIM);
      if (fabs(box[YY][XX]) > fabs(box[ZZ][XX]))
	d = YY; 
      else
	d = ZZ;
      if (box[d][XX] > 0)
	rvec_dec(box[d],box[XX]);
      else
	rvec_inc(box[d],box[XX]);
      pr_rvecs(stdlog,0,"new box",box,DIM);
    }
    /* (un)shifting should NOT be done after this,
     * since the box vectors might have changed
     */

    inc_nrnb(nrnb,eNR_PCOUPL,nr_atoms);
  }
}

void tcoupl(bool bTC,t_grpopts *opts,t_groups *grps,
	    real dt,real SAfactor)
{
  int    i;
  real   T,reft,lll;

  for(i=0; (i<opts->ngtc); i++) {
    reft=opts->ref_t[i]*SAfactor;
    if (reft < 0)
      reft=0;
    
    T = grps->tcstat[i].T;
    
    if ((bTC) && (T != 0.0)) {
      lll=sqrt(1.0 + (dt/opts->tau_t[i])*(reft/T-1.0));
      grps->tcstat[i].lambda=max(min(lll,1.25),0.8);
    }
    else
      grps->tcstat[i].lambda=1.0;
#ifdef DEBUGTC
    fprintf(stdlog,"group %d: T: %g, Lambda: %g\n",
	    i,T,grps->tcstat[i].lambda);
#endif
  }
}


