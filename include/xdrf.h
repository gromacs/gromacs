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
 * GROningen MAchine for Chemical Simulation
 */

#ifndef _xdrf_h
#define _xdrf_h

#ifdef USE_XDR
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
