/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * Gromacs Runs On Most of All Computer Systems
 */
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
  return xdr_float(xdrs,r);
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
  return xdr3dfcoord(xdrs,fp,size,precision);
#endif
}
