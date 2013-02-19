/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include "gmx_blas.h"
#include "gmx_lapack.h"

void
F77_FUNC(dgetrs,DGETRS)(const char *trans, 
	int *n, 
	int *nrhs, 
	double *a, 
	int *lda, 
	int *ipiv,
	double *b, 
	int *ldb, 
	int *info)
{
    int a_dim1, a_offset, b_dim1, b_offset;
    int notran;
    int c__1 = 1;
    int c_n1 = -1;
    double one = 1.0;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    *info = 0;
    notran = (*trans=='N' || *trans=='n');

    if (*n <= 0 || *nrhs <= 0) 
	return;

    if (notran) {
	F77_FUNC(dlaswp,DLASWP)(nrhs, &b[b_offset], ldb, &c__1, n, &ipiv[1], &c__1);
	F77_FUNC(dtrsm,DTRSM)("Left", "Lower", "No transpose", "Unit", n, nrhs, &one, 
	       &a[a_offset], lda, &b[b_offset], ldb);

	F77_FUNC(dtrsm,DTRSM)("Left", "Upper", "No transpose", "Non-unit", n, nrhs, &one, 
	       &a[a_offset], lda, &b[b_offset], ldb);
    } else {
	F77_FUNC(dtrsm,DTRSM)("Left", "Upper", "Transpose", "Non-unit", n, nrhs, &one, 
	       &a[a_offset], lda, &b[b_offset], ldb);
	F77_FUNC(dtrsm,DTRSM)("Left", "Lower", "Transpose", "Unit", n, nrhs, &one, 
	       &a[a_offset], lda, &b[b_offset], ldb);

	F77_FUNC(dlaswp,DLASWP)(nrhs, &b[b_offset], ldb, &c__1, n, &ipiv[1], &c_n1);
    }

    return;

} 
