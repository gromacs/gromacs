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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
static char *SRCID_ewald_util_c = "$Id$";

#include <stdio.h>
#include <math.h>
#include "assert.h"
#include "typedefs.h"
#include "vec.h"
#include "ewald_util.h"
#include "smalloc.h"
#include "physics.h"
#include "txtdump.h"
#include "futil.h"
#include "names.h"
#include "fftgrid.h"
#include "writeps.h"
#include "macros.h"
#include "xvgr.h"

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



real ewald_LRcorrection(FILE *fp,t_nsborder *nsb,t_commrec *cr,t_forcerec *fr,
			real charge[],t_block *excl,rvec x[],
			rvec box_size,matrix lr_vir)
{
  static  bool bFirst=TRUE;
  static  real Vself;

  int     i,i1,i2,j,k,m,iv,jv;
  atom_id *AA;
  double  qq; /* Necessary for precision */
  real    vc,qi,dr,ddd,dr2,rinv,fscal,Vexcl,rinv2,ewc=fr->ewaldcoeff;
  rvec    df,dx;
  rvec    *flr=fr->flr;
  /*#define TABLES*/
#ifdef TABLES
  real    tabscale=fr->tabscale;
  real    eps,eps2,VV,FF,F,Y,Geps,Heps2,Fp,fijC,r1t;
  real    *VFtab=fr->VFtab;
  int     n0,n1,nnn;
#else
  double  isp=0.564189583547756;
#endif
  int     start = START(nsb);
  int     end   = start+HOMENR(nsb);
  
  if (bFirst) {
    qq =0;  
    for(i=start; (i<end); i++) 
      qq  += charge[i]*charge[i];
    
    /* Charge and dipole correction remain to be implemented... */

    Vself=ewc*ONE_4PI_EPS0*qq/sqrt(M_PI);
    if (debug) {
      fprintf(fp,"Long Range corrections for Ewald interactions:\n");
      fprintf(fp,"start=%d,natoms=%d\n",start,end-start);
      fprintf(fp,"qq = %g, Vself=%g\n",qq,Vself);
    }
  }
  
  AA    = excl->a;
  Vexcl = 0;
 
  for(i=start; (i<end); i++) {
      /* Initiate local variables (for this i-particle) to 0 */
      qi  = charge[i]*ONE_4PI_EPS0;
      i1  = excl->index[i];
      i2  = excl->index[i+1];
      /* Loop over excluded neighbours */
      for(j=i1; (j<i2); j++) {
	k = AA[j];
	  /* 
	   * First we must test whether k <> i, and then, because the
	   * exclusions are all listed twice i->k and k->i we must select
	   * just one of the two.
       * As a minor optimization we only compute forces when the charges
       * are non-zero.
       */
      if (k > i) {
	qq = qi*charge[k];
	if (qq != 0.0) {
	  /* Compute distance vector, no PBC check! */
	  dr2 = 0;
	  for(m=0; (m<DIM); m++) {
	    ddd = x[i][m] - x[k][m];
	      if(ddd>box_size[m]/2) {  /* ugly hack,   */
		ddd-=box_size[m];      /* to fix pbc.. */
		/* Can this be done better? */
	      } else if (ddd<-box_size[m]/2)
		  ddd+=box_size[m];
	      
	    dx[m] = ddd;
	    dr2  += ddd*ddd;
	  }

	  /* It might be possible to optimize this slightly 
	   * when using spc or similar water models:
	   * Potential could be stored once, at the beginning,
	   * and forces could use that bonds are constant 
	   */
	  /* The Ewald interactions are tabulated. If you 
	   * remove the tabulation you must replace this
	   * part too
	   */

	  rinv              = invsqrt(dr2);
	  rinv2             = rinv*rinv;
	  dr                = 1.0/rinv;      
#ifdef TABLES
	  r1t               = tabscale*dr;   
	  n0                = r1t;
	  assert(n0 >= 3);
	  n1                = 12*n0;
	  eps               = r1t-n0;
	  eps2              = eps*eps;
	  nnn               = n1;
	  Y                 = VFtab[nnn];
	  F                 = VFtab[nnn+1];
	  Geps              = eps*VFtab[nnn+2];
	  Heps2             = eps2*VFtab[nnn+3];
	  Fp                = F+Geps+Heps2;
	  VV                = Y+eps*Fp;
	  FF                = Fp+Geps+2.0*Heps2;
	  vc                = qq*(rinv-VV);
	  fijC              = qq*FF;
	  Vexcl            += vc;
	  
	  fscal             = vc*rinv2+fijC*tabscale*rinv;
	  /* End of tabulated interaction part */
#else

	  /* This is the code you would want instead if not using
	   * tables: 
	   */
	  vc      = qq*erf(ewc*dr)*rinv;
	  Vexcl  += vc;
	  fscal   = rinv2*(vc-2.0*qq*ewc*isp*exp(-ewc*ewc*dr2));
#endif
	  
	  /* The force vector is obtained by multiplication with the 
	   * distance vector 
	   */
	  svmul(fscal,dx,df);
	  if (debug)
	    fprintf(debug,"dr=%8.4f, fscal=%8.0f, df=%10.0f,%10.0f,%10.0f\n",
		    dr,fscal,df[XX],df[YY],df[ZZ]);
	  rvec_inc(flr[k],df);
	  rvec_dec(flr[i],df);
	  for(iv=0; (iv<DIM); iv++)
	    for(jv=0; (jv<DIM); jv++)
	      lr_vir[iv][jv]+=0.5*dx[iv]*df[jv];
	}
      }
    }
  }
  if (bFirst && debug)
    fprintf(fp,"Long Range correction: Vexcl=%g\n",Vexcl);

  bFirst = FALSE;
  
  return (Vself+Vexcl);
}
  
