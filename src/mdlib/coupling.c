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
  if (nrdf > 0)
    return (2.0*ekin)/(nrdf*BOLTZ);
  else
    return 0;
}

void parinellorahman_pcoupl(t_inputrec *ir,int step,tensor pres,
			    tensor box,tensor boxv,tensor M)
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
  static tensor winv;
  static bool bFirst=TRUE;
  real   vol=box[XX][XX]*box[YY][YY]*box[ZZ][ZZ];
  real   fac=vol/PRESFAC;
  real   atot,arel,change,maxchange,xy_pressure;
  tensor invbox,pdiff,t1,t2;

  if(bFirst) {
    real maxl;
    maxl=max(box[XX][XX],box[YY][YY]);
    maxl=max(maxl,box[ZZ][ZZ]);
    for(d=0;d<DIM;d++)
      for(n=0;n<DIM;n++)
	winv[d][n]=(4*M_PI*M_PI*ir->compress[d][n])/(3*ir->tau_p*ir->tau_p*maxl);
    bFirst=FALSE;
  }
  
  m_inv(box,invbox);
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

  switch (ir->epct) {
  case epctANISOTROPIC:
    for(d=0;d<DIM;d++) 
      for(n=0;n<=d;n++)
	t1[d][n]*=winv[d][n]*fac;
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
	t1[d][n]=winv[d][n]*fac*arel*box[d][n];    
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
	t1[d][n]=winv[d][n]*fac*arel*box[d][n];
    for(n=0;n<DIM;n++)
      t1[ZZ][n]*=winv[d][n]*fac;
    break;
  default:
    fatal_error(0,"Parinello-Rahman pressure coupling type %s "
		"not supported yet\n",EPCOUPLTYPETYPE(ir->epct));
    break;
  }
  
  maxchange=0;
  for(d=0;d<DIM;d++)
    for(n=0;n<=d;n++) {
      boxv[d][n] += ir->delta_t*t1[d][n];
      /* We do NOT update the box vectors themselves here, since
       * we need them for shifting later. It is instead done last
       * in the update() routine.
       */
      if(box[d][n]!=0)
	change=fabs(ir->delta_t*boxv[d][n]/box[d][n]);
      else
	change=fabs(ir->delta_t*boxv[d][n]/box[d][d]);
      if(change>maxchange)
	maxchange=change;
    }
  
  if(maxchange>1.01) 
    fprintf(stdlog,"\nStep %d  Warning: Pressure scaling more than 1%%.\n",step);
 
  mtmul(boxv,box,t1);       /* t1=boxv * b' */
  for(d=0;d<DIM;d++)
    for(n=0;n<DIM;n++)
      t1[d][n] += t1[n][d]; /* t1=t1+t1' */
  mmul(invbox,t1,t2);
  mtmul(t2,invbox,M);
}


static void correct_box_elem(tensor box,int v,int d)
{
  /* correct elem d of vector v with vector d */
  while (box[v][d] > BOX_MARGIN*box[d][d]) {
    fprintf(stdlog,"Correcting invalid box:\n");
    pr_rvecs(stdlog,0,"old box",box,DIM);
    rvec_dec(box[v],box[d]);
    pr_rvecs(stdlog,0,"new box",box,DIM);
  } 
  while (-box[v][d] > BOX_MARGIN*box[d][d]) {
    fprintf(stdlog,"Correcting invalid box:\n");
    pr_rvecs(stdlog,0,"old box",box,DIM);
    rvec_inc(box[v],box[d]);
    pr_rvecs(stdlog,0,"new box",box,DIM);
  }
}

void correct_box(tensor box)
{
  /* check if the box still obeys the restrictions, if not, correct it */
  correct_box_elem(box,ZZ,YY);
  correct_box_elem(box,ZZ,XX);
  correct_box_elem(box,YY,XX);
}


void berendsen_pcoupl(t_inputrec *ir,int step,tensor pres,
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
    scalar_pressure=0;
    xy_pressure=0;
    for(d=0; d<DIM; d++) {
      scalar_pressure += pres[d][d]/DIM;
      if (d != ZZ)
	xy_pressure += pres[d][d]/(DIM-1);
    }
    /* Pressure is now in bar, everywhere. */
#define factor(d,m) (ir->compress[d][m]*ir->delta_t/ir->tau_p)
    
    /* mu has been changed from pow(1+...,1/3) to 1+.../3, since this is
     * necessary for triclinic scaling
     */
    if (scalar_pressure != 0.0) {
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
	fatal_error(0,"Berendsen pressure coupling type %s not supported yet\n",
		    EPCOUPLTYPETYPE(ir->epct));
	break;
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
	  +mu[XX][ZZ]*box[d][ZZ];
	box[d][YY] = mu[YY][YY]*box[d][YY]+mu[YY][ZZ]*box[d][ZZ];
	box[d][ZZ] = mu[ZZ][ZZ]*box[d][ZZ];
      }      

      /* (un)shifting should NOT be done after this,
       * since the box vectors might have changed
       */
      inc_nrnb(nrnb,eNR_PCOUPL,nr_atoms);
    }
}

void berendsen_tcoupl(t_grpopts *opts,t_groups *grps,
		      real dt,real SAfactor)
{
  int    i;
  real   T,reft=0,lll; 

  for(i=0; (i<opts->ngtc); i++) {
    reft=opts->ref_t[i]*SAfactor;
    if (reft < 0)
      reft=0;

    T = grps->tcstat[i].T;
    
    if (T != 0.0) {
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

void nosehoover_tcoupl(t_grpopts *opts,t_groups *grps,
		      real dt,real SAfactor)
{
  int i;
  real reft=0,xit,oldxi;
  real *Q=NULL,*Qinv=NULL;

  if(Qinv==NULL) {
    snew(Q,opts->ngtc);
    snew(Qinv,opts->ngtc);
    for(i=0;i<opts->ngtc;i++) {
      Q[i]=(opts->tau_t[i]*opts->tau_t[i]*opts->ref_t[i])/(4*M_PI*M_PI);
      Qinv[i]=1.0/Q[i];
    /* Use inputrec ref_t - Q shouldnt change during the run. */
    }
  }

  for(i=0; (i<opts->ngtc); i++) {
    reft=opts->ref_t[i]*SAfactor;
    if (reft < 0)
      reft=0;
    grps->tcstat[i].xi += dt*Qinv[i]*(grps->tcstat[i].T-reft);
  }
}
  
