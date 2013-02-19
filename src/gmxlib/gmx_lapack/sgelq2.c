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
#include "gmx_lapack.h"

void 
F77_FUNC(sgelq2,SGELQ2)(int *m, 
                        int *n, 
                        float *a,
                        int *lda, 
                        float *tau, 
                        float *work, 
                        int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    int i__, k;
    float aii;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    *info = 0;
    
    i__4 = (*m > 1) ? *m : 1;
    
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < i__4) {
	*info = -4;
    }
    if (*info != 0) {
	return;
    }

    
    k = (*m < *n ) ? *m : *n;
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n - i__ + 1;
	i__3 = i__ + 1;
    i__4 = (i__3 < *n) ? i__3 : *n;
	F77_FUNC(slarfg,SLARFG)(&i__2, &a[i__ + i__ * a_dim1], &a[i__ + i__4 * a_dim1],
            lda, &tau[i__]);
	if (i__ < *m) {
	    aii = a[i__ + i__ * a_dim1];
	    a[i__ + i__ * a_dim1] = 1.f;
	    i__2 = *m - i__;
	    i__3 = *n - i__ + 1;
	    F77_FUNC(slarf,SLARF)("R", &i__2, &i__3, &a[i__ + i__ * a_dim1], lda, 
               &tau[i__], &a[i__ + 1 + i__ * a_dim1], lda, &work[1]);
	    a[i__ + i__ * a_dim1] = aii;
	}
    }
    return;
}


