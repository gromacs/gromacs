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
#include "fatal.h"
#include "xdrf.h"
	
int xdropen(XDR *xdrs, const char *filename, const char *type)
{
  fatal_error(0,"xdropen called");
  
  return 0;
}

int xdrclose(XDR *xdrs) 
{
  fatal_error(0,"xdrclose called");
  
  return 0;
}

int xdr3dfcoord(XDR *xdrs, float *fp, int *size, float *precision) 
{
  fatal_error(0,"xdr3dfcoord called");
  
  return 0;
}

bool_t	xdr_int(XDR *xdr, int *i)
{
  fatal_error(0,"xdr_int called");
  
  return (bool_t) 0;
}

bool_t	xdr_float(XDR *xdr, float *f)
{
  fatal_error(0,"xdr_float called");
  
  return (bool_t) 0;
}
   
bool_t	xdr_double(XDR *xdr, double *d)
{
  fatal_error(0,"xdr_double called");
  
  return (bool_t) 0;
}

bool_t xdr_string(XDR *xdr,char **s,int size)
{
  fatal_error(0,"xdr_string called");
  
  return (bool_t) 0;
}

