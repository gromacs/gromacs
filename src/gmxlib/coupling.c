#include <stdio.h>
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

void calc_pres(matrix box,tensor ekin,tensor vir,tensor pres,real Elr)
{
  int  n,m;
  real fac,Plr;

  /* Uitzoeken welke ekin hier van toepassing is, zie Evans & Morris - E. */ 
  /* Wrs. moet de druktensor gecorrigeerd worden voor de netto stroom in het */
  /* systeem...       */
  
  /* Long range correctie for periodic systems, see
   * Neumann et al. JCP
   * divide by 6 because it is multiplied by fac later on.
   * If Elr = 0, no correction is made.
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

real calc_temp(real ekin,int nrdf)
{
  return (2.0*ekin)/(nrdf*BOLTZ);
}

real run_aver(real old,real cur,int step,int nmem)
{
  nmem   = max(1,nmem);
  
  return ((nmem-1)*old+cur)/nmem;
}


void do_pcoupl(t_inputrec *ir,int step,tensor pres,
	       matrix box,int start,int nr_atoms,
	       rvec x[],ushort cFREEZE[],
	       t_nrnb *nrnb,rvec freezefac[])
{
  static bool bFirst=TRUE;
  static rvec PPP;
  int    n,d,m,g,ncoupl=0;
  real   scalar_pressure, xy_pressure, p_corr_z;
  real   X,Y,Z,dx,dy,dz;
  rvec   factor;
  tensor mu;
  real   muxx,muxy,muxz,muyx,muyy,muyz,muzx,muzy,muzz;
  real   fgx,fgy,fgz;
  
  /*
   *  PRESSURE SCALING 
   *  Step (2P)
   */
  if (bFirst) {
    /* Initiate the pressure to the reference one */
    for(m=0; m<DIM; m++)
      PPP[m] = ir->ref_p[m];
    bFirst=FALSE;
  }
  scalar_pressure=0;
  xy_pressure=0;
  for(m=0; m<DIM; m++) {
    PPP[m]           = run_aver(PPP[m],pres[m][m],step,ir->npcmemory);
    scalar_pressure += PPP[m]/DIM;
    if (m != ZZ)
      xy_pressure += PPP[m]/(DIM-1);
  }
  
  /* Pressure is now in bar, everywhere. */
  if ((ir->epc != epcNO) && (scalar_pressure != 0.0)) {
    for(m=0; m<DIM; m++)
      factor[m] = ir->compress[m]*ir->delta_t/ir->tau_p;
    clear_mat(mu);
    switch (ir->epc) {
    case epcISOTROPIC:
      for(m=0; m<DIM; m++)
	mu[m][m] = pow(1.0-factor[m]*(ir->ref_p[m]-scalar_pressure),1.0/DIM);
      break;
    case epcSEMIISOTROPIC:
      for(m=0; m<ZZ; m++)
	mu[m][m] = pow(1.0-factor[m]*(ir->ref_p[m]-xy_pressure),1.0/DIM);
      mu[ZZ][ZZ] = pow(1.0-factor[ZZ]*(ir->ref_p[ZZ] - PPP[ZZ]),1.0/DIM);
      break;
    case epcANISOTROPIC:
      for (m=0; m<DIM; m++)
	mu[m][m] = pow(1.0-factor[m]*(ir->ref_p[m] - PPP[m]),1.0/DIM);
      break;
    case epcSURFACETENSION:
      /* ir->ref_p[0/1] is the reference surface-tension times *
       * the number of surfaces                                */
      if (ir->compress[ZZ])
	p_corr_z = ir->delta_t/ir->tau_p*(ir->ref_p[ZZ] - PPP[ZZ]);
      else
	/* when the compressibity is zero, set the pressure correction   *
	 * in the z-direction to zero to get the correct surface tension */
	p_corr_z = 0;
      mu[ZZ][ZZ] = 1.0 - ir->compress[ZZ]*p_corr_z;
      for(m=0; m<ZZ; m++)
	mu[m][m] = sqrt(1.0+factor[m]*(ir->ref_p[m]/(mu[ZZ][ZZ]*box[ZZ][ZZ]) - 
	(PPP[ZZ]+p_corr_z - xy_pressure)));
      break;
    case epcTRICLINIC:
    default:
      fatal_error(0,"Pressure coupling type %s not supported yet\n",
		  EPCOUPLTYPE(ir->epc));
    }
    if (debug) {
      pr_rvecs(debug,0,"PC: PPP ",&PPP,1);
      pr_rvecs(debug,0,"PC: fac ",&factor,1);
      pr_rvecs(debug,0,"PC: mu  ",mu,DIM);
    }
    /* Scale the positions using matrix operation */
    nr_atoms+=start;
    muxx=mu[XX][XX],muxy=mu[XX][YY],muxz=mu[XX][ZZ];
    muyx=mu[YY][XX],muyy=mu[YY][YY],muyz=mu[YY][ZZ];
    muzx=mu[ZZ][XX],muzy=mu[ZZ][YY],muzz=mu[ZZ][ZZ];
    for (n=start; n<nr_atoms; n++) {
      g=cFREEZE[n];
      fgx=freezefac[g][XX];
      fgy=freezefac[g][YY];
      fgz=freezefac[g][ZZ];
      
      X=x[n][XX];
      Y=x[n][YY];
      Z=x[n][ZZ];
      dx=muxx*X+muxy*Y+muxz*Z;
      dy=muyx*X+muyy*Y+muyz*Z;
      dz=muzx*X+muzy*Y+muzz*Z;
      x[n][XX]=X+fgx*(dx-X);
      x[n][YY]=Y+fgy*(dy-Y);
      x[n][ZZ]=Z+fgz*(dz-Z);
      
      ncoupl++;
    }
    /* compute final boxlengths */
    for (d=0; d<DIM; d++)
      for (m=0; m<DIM; m++)
	box[d][m] *= mu[d][m];
  }
  inc_nrnb(nrnb,eNR_PCOUPL,ncoupl);
}

void tcoupl(bool bTC,t_grpopts *opts,t_groups *grps,
	    real dt,real SAfactor,int step,int nmem)
{
  static real *Told=NULL;
  int    i;
  real   T,reft,lll;

  if (!Told) {
    snew(Told,opts->ngtc);
    for(i=0; (i<opts->ngtc); i++) 
      Told[i]=opts->ref_t[i]*SAfactor;
  }
  
  for(i=0; (i<opts->ngtc); i++) {
    reft=opts->ref_t[i]*SAfactor;
    if (reft < 0)
      reft=0;
    
    Told[i] = run_aver(Told[i],grps->tcstat[i].T,step,nmem);
    T       = Told[i];
    
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


