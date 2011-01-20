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
 * Gromacs Runs On Most of All Computer Systems
 */
#ifndef _complex_h
#define _complex_h

#include <math.h>
#include "typedefs.h"

typedef struct {
  real re,im;
} t_complex;

typedef t_complex cvec[DIM];

static t_complex cnul = { 0.0, 0.0 };

static t_complex rcmul(real r,t_complex c)
{
  t_complex d;
  
  d.re = r*c.re;
  d.im = r*c.im;
  
  return d;
}

static t_complex rcexp(real r)
{
  t_complex c;
  
  c.re = (real)cos(r);
  c.im = (real)sin(r);
  
  return c;
}


static t_complex cadd(t_complex a,t_complex b)
{
  t_complex c;
  
  c.re = a.re+b.re;
  c.im = a.im+b.im;
  
  return c;
}

static t_complex csub(t_complex a,t_complex b)
{
  t_complex c;
  
  c.re = a.re-b.re;
  c.im = a.im-b.im;
  
  return c;
}

static t_complex cmul(t_complex a,t_complex b)
{
  t_complex c;
  
  c.re = a.re*b.re - a.im*b.im;
  c.im = a.re*b.im + a.im*b.re;
  
  return c;
}

static t_complex conjugate(t_complex c)
{
  t_complex d;
  
  d.re =  c.re;
  d.im = -c.im;
  
  return d;
}

static t_complex cdiv(t_complex teller,t_complex noemer)
{
  t_complex res,anoemer;
  
  anoemer = cmul(conjugate(noemer),noemer);
  res = cmul(teller,conjugate(noemer));
  
  return rcmul(1.0/anoemer.re,res);
}
#endif
