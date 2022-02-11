/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#ifndef GMX_LINEARALGEBRA_MATRIX_H
#define GMX_LINEARALGEBRA_MATRIX_H

#include <stdio.h>

double** alloc_matrix(int n, int m);

void free_matrix(double** a);

void matrix_multiply(FILE* fp, int n, int m, double** x, double** y, double** z);

/* Return 0 if OK or row number where inversion failed otherwise. */
int matrix_invert(FILE* fp, int n, double** a);

double multi_regression(FILE* fp, int ny, double* y, int nx, double** xx, double* a0);
/* Perform a regression analysis to fit
 * y' = a0[0] xx[0] + a0[1] xx[1] ... + a0[nx-1] xx[nx-1]
 * with ny data points in each vector.
 * The coefficients are returned in vector a0.
 * The return value of the function is the chi2 value:
 * sum_{j=0}^{ny-1} (y[j] - y'[j])^2
 * If fp is not NULL debug information will be written to it.
 */

#endif
