/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

#ifndef _complex_h
#define _complex_h

static char *SRCID_complex_h = "$Id$";

#include "typedefs.h"

typedef struct {
  real re,im;
} t_complex;

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
  
  c.re = cos(r);
  c.im = sin(r);
  
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
