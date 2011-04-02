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

#ifndef _nrjac_h
#define _nrjac_h

#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

void jacobi(double **a,int n,double d[],double **v,int *nrot);
/* 
 * real   **omega = input matrix a[0..n-1][0..n-1] must be symmetric
 * int     natoms = number of rows and columns
 * real      NULL = d[0]..d[n-1] are the eigenvalues of a[][]
 * real       **v = v[0..n-1][0..n-1] the eigenvectors:
 *                                    v[i][j] is component i of vector j
 * int      *irot = number of jacobi rotations
 */

int m_inv_gen(real **m,int n,real **minv);
/* Produces minv, a generalized inverse of m.
 * Inversion is done via diagonalization,
 * eigenvalues smaller than 1e-6 times the average diagonal element
 * are assumed to be zero.
 * For zero eigenvalues 1/eigenvalue is set to zero for the inverse matrix.
 * Returns the number of zero eigenvalues.
 */

#ifdef __cplusplus
}
#endif

#endif
