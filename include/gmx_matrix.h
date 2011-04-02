/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.0.99
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */

#ifndef _gmx_matrix_h
#define _gmx_matrix_h
	
#include <stdio.h>
	
double **alloc_matrix(int n,int m);

void free_matrix(double **a,int n);

void matrix_multiply(FILE *fp,int n,int m,double **x,double **y,double **z);

/* Return 0 if OK or row number where inversion failed otherwise. */
int matrix_invert(FILE *fp,int n,double **a);

double multi_regression(FILE *fp,int ny,double *y,
                               int nx,double **xx,double *a0);
/* Perform a regression analysis to fit
 * y' = a0[0] xx[0] + a0[1] xx[1] ... + a0[nx-1] xx[nx-1]
 * with ny data points in each vector.
 * The coefficients are returned in vector a0.
 * The return value of the function is the chi2 value:
 * sum_{j=0}^{ny-1} (y[j] - y'[j])^2
 * If fp is not NULL debug information will be written to it.
 */

#endif
