/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Giving Russians Opium May Alter Current Situation
 */

#include <typedefs.h>
#include <vec.h>
#include <detectcpu.h>

/* Fast vectorized versions of 1/sqrt(x) and 1/x.
 * Assembly routines are used for Alpha processors, AMD processors
 * with 3DNow. On IBM we use the MASS libraries if they are present.
 * 
 * WE ALSO HAVE AN ASSEMBLY LOOP FOR PENTIUM PROCESSORS, BUT IT CAN
 * ONLY BE USED IF YOUR INPUT AND OUTPUT ARRAYS ARE ALIGNED TO
 * THE CACHE LINE!
 * This is not complicated to do, but you must do it when you allocate
 * your memory. Start by allocating 31 bytes more than you need, put
 * this in a temp variable (e.g. _buf, so you can free it later), and
 * create your aligned array buf with
 * 
 *  buf=(real *) ( ( (unsigned long int)_buf + 31 ) & (~0x1f) );
 *
 * And, of course, simliar for your output buffer. 
 * If you have an SSE-enabled CPU (pentium III and later) and OS
 * (Linux 2.4 and later) you will now be able to perform 1/sqrt(x)
 * in 3-4 clocks/element for long lists!
 */
 
void vecinvsqrt(real in[],real out[],int n)
{
  /* Always use assembly on Alpha CPU. */
#ifdef USE_AXP_ASM
#  ifdef DOUBLE
   sqrtiv_(in,out,&n);
#  else /* SINGLE */
   ssqrtiv_(in,out,&n);
#  endif    
#elif defined HAVE_LIBMASSV_ANY
 /* On IBM we should definitely use vectorized MASS if present. */ 
#  ifdef DOUBLE
   vrsqrt(out,in,&n);
#  else /* SINGLE */
   vsrsqrt(out,in,&n);
#  endif 
#else /* not alpha, and not IBM with MASS */
   /* Software routines and calls to x86 assembly. */
#ifdef SOFTWARE_INVSQRT
  const real  half=0.5;
  const real  three=3.0;
  t_convert   result,bit_pattern;
  unsigned int exp,fract;
  float       lu,x;
#ifdef DOUBLE
  real        y;
#endif
#endif /* VARIABLES FOR SOFTWARE_INVSQRT */
  int i;
 
#if (defined USE_X86_ASM && !defined DOUBLE)
  static bool bFirst=TRUE;
  static int cpu_capabilities;
  
  if(bFirst) {
    cpu_capabilities=check_x86cpu(NULL);
    bFirst=FALSE;
  }

  if((cpu_capabilities & X86_SSE_SUPPORT) && !((unsigned long int)in & 0x1f) && !((unsigned long int)out & 0x1f)) /* SSE data must be cache aligned */
    vecinvsqrt_sse(in,out,n);
  else if(cpu_capabilities & X86_3DNOW_SUPPORT)
    vecinvsqrt_3dnow(in,out,n);
  else
#endif /* no x86 optimizations */
#ifdef SOFTWARE_INVSQRT
    for(i=0;i<n;i++) {
      x=in[i];
      bit_pattern.fval=x;
      exp   = EXP_ADDR(bit_pattern.bval);
      fract = FRACT_ADDR(bit_pattern.bval);
      result.bval=cinvsqrtexptab[exp] | cinvsqrtfracttab[fract];
      lu    = result.fval;
      
#ifdef DOUBLE
      y=(half*lu*(three-((x*lu)*lu)));
      out[i]=(half*y*(three-((x*y)*y)));
#else
      out[i]=(half*lu*(three-((x*lu)*lu)));
#endif
    }
#else  /* no gmx invsqrt */
    for(i=0;i<n;i++)
      out[i]=1.0f/sqrt(in[i]);
#endif /* SOFTWARE_SQRT */
#endif
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


#ifdef SOFTWARE_RECIP
  const real  two=2.0;
  t_convert   result,bit_pattern;
  unsigned int exp,fract;
  float       lu,x;
#ifdef DOUBLE
  real        y;
#endif
#endif /* SOFTWARE_RECIP */
  int i;

#if (defined USE_X86_ASM && !defined DOUBLE)
  static bool bFirst=TRUE;
  static int cpu_capabilities;

  if(bFirst) {
    cpu_capabilities=check_x86cpu(NULL);
    bFirst=FALSE;
  }

  if((cpu_capabilities & X86_SSE_SUPPORT) && !((unsigned long int)in & 0x1f) && !((unsigned long int)out & 0x1f)) /* SSE data must be cache aligned */
    vecrecip_sse(in,out,n);
  else if(cpu_capabilities & X86_3DNOW_SUPPORT)
    vecrecip_3dnow(in,out,n);
  else
#endif /* no x86 optimizations */
#ifdef SOFTWARE_RECIP
    for(i=0;i<n;i++) {
      x=in[i];
      bit_pattern.fval=x;
      exp   = EXP_ADDR(bit_pattern.bval);
      fract = FRACT_ADDR(bit_pattern.bval);
      result.bval=crecipexptab[exp] | crecipfracttab[fract];
      lu    = result.fval;
      
#ifdef DOUBLE
      y=lu*(two-x*lu);
      out[i]=y*(two-x*y);
#else
      out[i]=lu*(two-x*lu);
#endif
    }
#else /* No gmx recip */ 
    for(i=0;i<n;i++)
      out[i]=1.0f/(in[i]);
#endif /* SOFTWARE_RECIP */
#endif
}






