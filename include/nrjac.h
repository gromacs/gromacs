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
 * Good ROcking Metal Altar for Chronical Sinners
 */

#ifndef _nrjac_h
#define _nrjac_h

static char *SRCID_nrjac_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

extern void jacobi(double **a,int n,double d[],double **v,int *nrot);
/* 
 * real   **omega = input matrix a[0..n-1][0..n-1] must be symmetric
 * int     natoms = number of rows and columns
 * real      NULL = d[0]..d[n-1] are the eigenvalues of a[][]
 * real       **v = v[0..n-1][0..n-1] contains the vectors in columns
 * int      *irot = number of jacobi rotations
 */
#endif
