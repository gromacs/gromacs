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
 * Great Red Owns Many ACres of Sand 
 */
static char *SRCID_dumxdrf_c = "$Id$";
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

