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
 * S  C  A  M  O  R  G
 */
static char *SRCID_xdrd_c = "$Id$";
#include "typedefs.h"
#include "xdrf.h"
#include "fatal.h"
#include "smalloc.h"

int xdr_real(XDR *xdrs,real *r)
{
#ifdef DOUBLE
  float f;
  int   ret;
  
  f=*r;
  ret=xdr_float(xdrs,&f);
  *r=f;

  return ret;
#else
  return xdr_float(xdrs,(float *)r);
#endif
}

int xdr3drcoord(XDR *xdrs, real *fp, int *size, real *precision)
{
#ifdef DOUBLE
  static float *ffp=NULL;
  float  fprec;
  int    i,ret,isize;
  
  isize=*size*DIM;
  if (ffp == NULL)  {
    if (isize > 0) {
      snew(ffp,isize);
    }
    else
      fatal_error(0,"Don't know what to malloc for ffp (file %s, line %d)",
		  __FILE__,__LINE__);
  }
  for(i=0; (i<isize); i++)
    ffp[i]=fp[i];
  fprec=*precision;
  ret=xdr3dfcoord(xdrs,ffp,size,&fprec);
  
  *precision=fprec;
  for(i=0; (i<isize); i++)
    fp[i]=ffp[i];
  
  return ret;
#else
  return xdr3dfcoord(xdrs,(float *)fp,size,(float *)precision);
#endif
}
