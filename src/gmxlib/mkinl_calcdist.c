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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
/* This file is NOT threadsafe, but it is only used to create
 * the innerloops during the build process, so it will never be
 * executed by multiple threads.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "mkinl.h"
#include <string.h>


void call_vectorized_routines()
{
  /* m was updated a last time after the assignments. Reduce
   * it if we are using fortran, so we avoid referencing m-1
   * repeatedly (and some functions need an adress to a variable
   * rather than a value)
   */
  if(!bC)
    decrement("m","1");
  
  /* always start with vectorized invsqrt if we should do it
   */
  if(loop.vectorize_invsqrt) {
    /* if we need the sqrt (non-reciprocal) later
     * we have to save the square values in buf1 and write to buf2,
     * otherwise we can gain some speed by writing back to the
     * same array buf1
     */
#ifdef HAVE_LIBMASSV_ANY
    /* IBM provide special fast single and double routines in their vectorized MASS libraries */
#ifdef DOUBLE
    code( bC ? "vrsqrt(%s,buf1,&m);" : "call vrsqrt(%s,buf1,m)",
	  OVERWRITE_RSQ ? "buf1" : "buf2" );	  
#else
    code( bC ? "vsrsqrt(%s,buf1,&m);" : "call vsrsqrt(%s,buf1,m)",
	  OVERWRITE_RSQ ? "buf1" : "buf2" );
#endif /* End of IBM routines */
#elif defined(USE_AXP_ASM)
    /* Compaq provide a math library too, but not an inverse square root or
     * a fast reciprocal. Our own assembly routines are much faster! */
#ifdef DOUBLE
    code( bC ? "sqrtiv_(buf1,%s,&m);" : "call sqrtiv(buf1,%s,m)",
	  OVERWRITE_RSQ ? "buf1" : "buf2" );
#else
    code( bC ? "ssqrtiv_(buf1,%s,&m);" : "call ssqrtiv(buf1,%s,m)",
	  OVERWRITE_RSQ ? "buf1" : "buf2" );
#endif /* End of code for compaq alpha */
#else 
    /* this is a vanilla GMX routine, suitable for most chips w/o hardware sqrt.
     * It is quite fast though - about 18 clocks on intel! */

    code( bC ? "vecinvsqrt(buf1,%s,m);" : "call f77vecinvsqrt(buf1,%s,m)",
	  OVERWRITE_RSQ ? "buf1" : "buf2" );
 
#endif
  } else if(loop.vectorize_recip) {
      /* Only do this if we didnt do invsqrt. If we did,
       * it is probably faster to square the value when
       * we use it.
       */
#ifdef HAVE_LIBMASSV_ANY
    /* IBM stuff again. They provide a nice vectorized reciprocal, too */
#ifdef DOUBLE
    code( bC ? "vrec(%s,buf1,&m);" : "call vrec(%s,buf1,m)",
	  OVERWRITE_RSQ ? "buf1" : "buf2" );
#else
    code( bC ? "vsrec(%s,buf1,&m);" : "call vsrec(%s,buf1,m)",
	  OVERWRITE_RSQ ? "buf1" : "buf2" );
#endif
#else
	/* on some architectures, though, the compiler can 
	 * link with appropriate libs if the code is easily recognized as
	 * vectorizable.
	 */
    vector_pragma();
    code( bC ? "vecrecip(buf1,%s,m);" : "call f77vecrecip(buf1,%s,m)",
	  OVERWRITE_RSQ ? "buf1" : "buf2" );
   
#endif
  }
}


int calc_rsq(int i, int j)
{
  char lm[12],lm2[32],sqbuf[64];
  int  m;

  sqbuf[0]=0;
  for(m='x'; (m<='z'); m++) {
    sprintf(lm,"d%c%d%d",m,i,j);
    sprintf(lm2,"%s*%s",lm,lm);
    add_to_buffer(sqbuf,lm2);
  }
  if(DO_VECTORIZE) {
    assign(ARRAY(buf1,m),sqbuf);
    increment("m","1");
  } else 
    assign("rsq%d%d",sqbuf,i,j);
  return 5;
}


int calc_dist()
{
  int nflop=0;
  int i,j,m;
  char buf[10];
  /* First calculate dr, and possibly store it */

  comment("Calculate distance vector"); 

  for(i=1;i<=loop.ni;i++)
      for(j=1;j<=loop.nj;j++) {
          /* first calculate dr */
          for(m='x'; (m<='z'); m++) { 
              subtract("d%c%d%d","i%c%d","j%c%d",m,i,j,m,i,m,j);
              nflop ++;
              if(DO_VECTORIZE && DO_FORCE) {
                  assign(ARRAY(drbuf,m3),"d%c%d%d",m,i,j);
                  increment("m3","1");
              }
          }
          nflop += calc_rsq(i,j);
      }

   return nflop;
}


int calc_rinv_and_rinvsq()
{
  int nflop=0;
  int i,j;

  if(!DO_VECTORIZE && !opt.delay_invsqrt) {
    if(loop.invsqrt)
      nflop += calc_invsqrt();
    else if(loop.recip) {
      for(j=1;j<=loop.nj;j++)
	for(i=1;i<=loop.ni;i++) {
	  assign("rinvsq%d", "1.0/rsq%d",i*10+j,i*10+j);
#ifdef DOUBLE
	  nflop += 6;
#else
	  nflop += 3;
#endif
	}
    }
  }
  return nflop;
}
