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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <typedefs.h>
#include <vec.h>
#include <detectcpu.h>

/* Fast vectorized versions of 1/sqrt(x) and 1/x.
 * Assembly routines are used for Alpha processors, AMD processors
 * with 3DNow. On IBM we use the MASS libraries if they are present.
 * 
 * We cannot call the SSE/3DNow/SSE2/Altivec from this general loop,
 * since the CPU detection is costly and we dont want to maintain
 * a separate variable for the cpu detected. Call detectcpu()
 * yourself, and use the SSE/3DNow/SSE2/Altivec versions directly
 * if the flags match. 
 */
 
void vecinvsqrt(real in[],real out[],int n)
{
  /* Always use assembly on Alpha CPU. */
#ifdef USE_AXP_ASM
#  ifdef DOUBLE
   sqrtiv_(in,out,&n);
#  else /* SINGLE */
   ssqrtiv_(in,out,&n);
#  endif /* SINGLE/DOUBLE */   
#elif defined HAVE_LIBMASSV_ANY
 /* On IBM we should definitely use vectorized MASS if present. */ 
#  ifdef DOUBLE
   vrsqrt(out,in,&n);
#  else /* SINGLE */
   vsrsqrt(out,in,&n);
#  endif  /* SINGLE/DOUBLE */   
#elif defined SOFTWARE_INVSQRT /* gromacs software 1/sqrt*/
  const real  half=0.5;
  const real  three=3.0;
  t_convert   result,bit_pattern;
  unsigned int exp,fract;
  float       lu;
  real        x;
#ifdef DOUBLE
  real        y;
#endif /* DOUBLE */
#else /* NO INVSQRT */
  int i; 
    for(i=0;i<n;i++)
      out[i]=1.0f/sqrt(in[i]);
#endif /* SOFTWARE_INVSQRT */
}



void vecrecip(real in[],real out[],int n)
{
/* No vectorized 1/x on alpha chips */

/* On IBM we should definitely use vectorized MASS if present. */ 
#ifdef HAVE_LIBMASSV_ANY
#  ifdef DOUBLE
   vrec(out,in,&n);
#  else /* SINGLE */
   vsrec(out,in,&n);
#  endif 
#else /* not IBM with MASS */

  int i;

    for(i=0;i<n;i++)
      out[i]=1.0f/(in[i]);
#endif
}






