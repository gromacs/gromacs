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
 * GRoups of Organic Molecules in ACtion for Science
 */

#ifndef _xdrf_h
#define _xdrf_h

static char *SRCID_xdrf_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef USE_XDR
#ifdef __GNUC__  
#ifdef FALSE
#undef FALSE
#endif
#ifdef TRUE
#undef TRUE
#endif
#endif

#ifdef __PGI    /*Portland group compiler*/
#define int64_t long long
#endif

#include <rpc/rpc.h>
#include <rpc/xdr.h>
#else
typedef int XDR;
typedef int bool_t;
#endif

#include "typedefs.h"

extern int xdropen(XDR *xdrs, const char *filename, const char *type);

extern int xdrclose(XDR *xdrs);

extern int xdr3dfcoord(XDR *xdrs, float *fp, int *size, float *precision);
/* Read or write reduced precision *float* coordinates */

extern int xdr_real(XDR *xdrs,real *r); 
/* Read or write a *real* value (stored as float) */

extern int xdr3drcoord(XDR *xdrs,real *fp,int *size,real *precision);
/* Read or write reduced precision *real* coordinates */

#endif
