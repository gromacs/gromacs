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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <math.h>
#include "typedefs.h"
#include "smalloc.h"
#include "maths.h"
#include "glaasje.h"
#include "macros.h"
	
void do_glas(FILE *log,int start,int homenr,rvec x[],rvec f[],
	     t_forcerec *fr,t_mdatoms *md,int atnr,t_inputrec *ir,
	     real ener[])
{
  static bool   bFirst=TRUE,bGlas;
  static real   d[2],pi6,pi12,rc9,rc4,rc10,rc3,rc;
  static real   *c6,*c12;
  real wd,wdd,zi,fz,dd,d10,d4,d9,d3,r9,r3,sign,cc6,cc12;
  int  *type;
  int  i,j,ti;
  
  type=md->typeA;
  if (bFirst) {
    pi6  = ir->userreal1;
    pi12 = ir->userreal2;
    d[0] = ir->userreal3;
    d[1] = ir->userreal4;
    
    /* Check whether these constants have been set. */
    bGlas = (pi6 != 0) && (pi12 != 0) && (d[0] != 0) && (d[1] != 0);
    
    if (bGlas) {
      if (ir->eDispCorr != edispcNO) {
	fatal_error(0,"Can not have Long Range C6 corrections and GLASMD");
      }
      rc   = max(fr->rvdw,fr->rlist);
      rc3  = rc*rc*rc;
      rc4  = rc3*rc;
      rc9  = rc3*rc3*rc3;
      rc10 = rc9*rc;
    
      fprintf(log,
	      "Constants for GLASMD: pi6 = %10g, pi12 = %10g\n"
	      "                      d1  = %10g, d2   = %10g\n"
	      "                      rc3 = %10g, rc4  = %10g\n"
	      "                      rc9 = %10g, rc10 = %10g\n",
	      pi6,pi12,d[0],d[1],rc3,rc4,rc9,rc10);
      if (d[0] > d[1])
	fatal_error(0,"d1 > d2 for GLASMD (check log file)");
    
      snew(c6,atnr);
      snew(c12,atnr);
    
      for(i=0; (i<atnr); i++) {
	c6[i]  = C6 (fr->nbfp,atnr,i,i);
	c12[i] = C12(fr->nbfp,atnr,i,i);
      }
    }
    else
      fprintf(stderr,"No glasmd!\n");
    bFirst = FALSE;
  }
  
  if (bGlas) {
    wd=0;
    for(i=start; (i<start+homenr); i++) {
      ti   = type[i];
      if ((c6[ti] != 0) || (c12[ti] != 0)) {
	zi   = x[i][ZZ];
	cc6  = M_PI*sqrt(c6[ti]*pi6);
	cc12 = M_PI*sqrt(c12[ti]*pi12);
	
	/* Use a factor for the sign, this saves taking absolute values */
	sign = 1;
	for(j=0; (j<2); j++) {
	  dd = sign*(zi-d[j]);
	  if (dd >= rc) {
	    d3  = dd*dd*dd;
	    d9  = d3*d3*d3;
	    wdd = cc12/(45.0*d9) - cc6/(6.0*d3);
	    d4  = d3*dd;
	    d10 = d9*dd;
	    fz  = sign*(cc12/(5.0*d10) - cc6/(2.0*d4));
	  }
	  else {
	    wdd = cc12*(2.0/(9.0*rc9) - dd/(5.0*rc10)) -
	      cc6*(2.0/(3.0*rc3) - dd/(2.0*rc4));
	    fz  = sign*(cc12/(5.0*rc10)-cc6/(2.0*rc4));
	  }
	  wd       += wdd;
	  f[i][ZZ] += fz;
	  sign      = -sign;
	}
      }
    }
    ener[F_LJ_LR] = wd;
  }
}

