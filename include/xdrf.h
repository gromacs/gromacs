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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _xdrf_h
#define _xdrf_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"

#ifdef __PGI    /*Portland group compiler*/
#define int64_t long long
#endif

#if (defined WIN32 || defined _WIN32 || defined WIN64  || defined _WIN64 || defined __CYGWIN__ || defined __CYGWIN32__ || defined GMX_INTERNAL_XDR)
#include "gmx_system_xdr.h"
#else
#include <rpc/rpc.h>
#include <rpc/xdr.h>
#endif

int 
xdropen(XDR *xdrs, const char *filename, const char *type);


int 
xdrclose(XDR *xdrs);


/* Read or write reduced precision *float* coordinates */
int 
xdr3dfcoord(XDR *xdrs, float *fp, int *size, float *precision);


/* Read or write a *real* value (stored as float) */
int 
xdr_real(XDR *xdrs,real *r); 


/* Read or write reduced precision *real* coordinates */
int 
xdr3drcoord(XDR *xdrs,real *fp,int *size,real *precision);


int 
xtc_seek_time(real time, int fp, int natoms);


int 
xtc_seek_frame(int frame, int fp, int natoms);


float 
xtc_get_last_frame_time(int fp, int natoms, bool * bOK);


int 
xtc_get_last_frame_number(int fp, int natoms, bool * bOK);

#endif







