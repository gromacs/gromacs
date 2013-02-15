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
#include "lapack_limits.h"

void
F77_FUNC(sgetri,SGETRI)(int *n, 
	float *a, 
	int *lda, 
	int *ipiv, 
	float *work, 
	int *lwork, 
	int *info)
{
    int a_dim1, a_offset, i__1, i__2, i__3;

    int i__, j, jb, nb, jj, jp, nn, iws;
    int nbmin;
    int ldwork;
    int lwkopt;
    int c__1 = 1;
    float c_b20 = -1.;
    float c_b22 = 1.;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    --work;

    *info = 0;
    nb = DGETRI_BLOCKSIZE;
    lwkopt = *n * nb;
    work[1] = (float) lwkopt;

    if (*n < 0) {
	*info = -1;
    } else if (*lda < (*n)) {
	*info = -3;
    } else if (*lwork < (*n) && *lwork!=-1) {
	*info = -6;
    }
    if (*info != 0) {
	i__1 = -(*info);
	return;
    } else if (*lwork == -1) {
	return;
    }

    if (*n == 0) {
	return;
    }

    F77_FUNC(strtri,STRTRI)("Upper", "Non-unit", n, &a[a_offset], lda, info);
    if (*info > 0) {
	return;
    }

    nbmin = 2;
    ldwork = *n;
    if (nb > 1 && nb < *n) {
	i__1 = ldwork * nb;
	iws = (i__1>1) ? i__1 : 1;
	if (*lwork < iws) {
	    nb = *lwork / ldwork;
	    nbmin = DGETRI_MINBLOCKSIZE;
	}
    } else {
	iws = *n;
    }

    if (nb < nbmin || nb >= *n) {

	for (j = *n; j >= 1; --j) {

	    i__1 = *n;
	    for (i__ = j + 1; i__ <= i__1; ++i__) {
		work[i__] = a[i__ + j * a_dim1];
		a[i__ + j * a_dim1] = 0.;
	    }

	    if (j < *n) {
		i__1 = *n - j;
		F77_FUNC(sgemv,SGEMV)("No transpose", n, &i__1, &c_b20, &a[(j + 1) * a_dim1 
			+ 1], lda, &work[j + 1], &c__1, &c_b22, &a[j * a_dim1 
			+ 1], &c__1);
	    }
	}
    } else {

	nn = (*n - 1) / nb * nb + 1;
	i__1 = -nb;
	for (j = nn; i__1 < 0 ? j >= 1 : j <= 1; j += i__1) {
	    i__2 = nb, i__3 = *n - j + 1;
	    jb = (i__2<i__3) ? i__2 : i__3;

	    i__2 = j + jb - 1;
	    for (jj = j; jj <= i__2; ++jj) {
		i__3 = *n;
		for (i__ = jj + 1; i__ <= i__3; ++i__) {
		    work[i__ + (jj - j) * ldwork] = a[i__ + jj * a_dim1];
		    a[i__ + jj * a_dim1] = 0.;
		}
	    }

	    if (j + jb <= *n) {
		i__2 = *n - j - jb + 1;
		F77_FUNC(sgemm,SGEMM)("No transpose", "No transpose", n, &jb, &i__2, &c_b20, 
			&a[(j + jb) * a_dim1 + 1], lda, &work[j + jb], &
			ldwork, &c_b22, &a[j * a_dim1 + 1], lda);
	    }
	    F77_FUNC(strsm,STRSM)("Right", "Lower", "No transpose", "Unit", n, &jb, &c_b22, &
		    work[j], &ldwork, &a[j * a_dim1 + 1], lda);
	}
    }

    for (j = *n - 1; j >= 1; --j) {
	jp = ipiv[j];
	if (jp != j) {
	    F77_FUNC(sswap,SSWAP)(n, &a[j * a_dim1 + 1], &c__1, &a[jp * a_dim1 + 1], &c__1);
	}
    }

    work[1] = (float) iws;
    return;

}


