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
#include <math.h>


#include "gmx_lapack.h"

double 
F77_FUNC(dlansy,DLANSY)(const char *norm, const char *uplo, int *n, double *a, int 
	*lda, double *work)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2;
    double ret_val, d__1, d__2, d__3;
    int c__1 = 1;

    /* Local variables */
    int i__, j;
    double sum, absa, scale;
    double value =0.0;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --work;

    if (*n == 0) {
	value = 0.;
    } else if (*norm=='M' || *norm=='m') {

	value = 0.;
	if (*uplo=='U' || *uplo=='u') {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = j;
		for (i__ = 1; i__ <= i__2; ++i__) {
		  d__2 = value;
		  d__3 = fabs(a[i__ + j * a_dim1]);
		  value = (d__2>d__3) ? d__2 : d__3;
		}
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = j; i__ <= i__2; ++i__) {
		  d__2 = value;
		  d__3 = fabs(a[i__ + j * a_dim1]);
		    value =  (d__2>d__3) ? d__2 : d__3;
		}
	    }
	}
    } else if (*norm=='I' || *norm=='i' || *norm=='O' || *norm=='o' || *norm=='1') {

	value = 0.;
	if (*uplo=='U' || *uplo=='u') {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		sum = 0.;
		i__2 = j - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    absa = fabs(a[i__ + j * a_dim1]);
		    sum += absa;
		    work[i__] += absa;
		}
		work[j] = sum + fabs(a[j + j * a_dim1]);
	    }
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		d__1 = value, d__2 = work[i__];
		value =  (d__1>d__2) ? d__1 : d__2;
	    }
	} else {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		work[i__] = 0.;
	    }
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		sum = work[j] + fabs(a[j + j * a_dim1]);
		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    absa = fabs(a[i__ + j * a_dim1]);
		    sum += absa;
		    work[i__] += absa;
		}
		if(sum>value)
		  value = sum;
	    }
	}
    } else if (*norm=='F' || *norm=='f' || *norm=='E' || *norm=='e') {

	scale = 0.;
	sum = 1.;
	if (*uplo=='U' || *uplo=='u') {
	    i__1 = *n;
	    for (j = 2; j <= i__1; ++j) {
		i__2 = j - 1;
		F77_FUNC(dlassq,DLASSQ)(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
	    }
	} else {
	    i__1 = *n - 1;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n - j;
		F77_FUNC(dlassq,DLASSQ)(&i__2, &a[j + 1 + j * a_dim1], &c__1, &scale, &sum);
	    }
	}
	sum *= 2;
	i__1 = *lda + 1;
	F77_FUNC(dlassq,DLASSQ)(n, &a[a_offset], &i__1, &scale, &sum);
	value = scale * sqrt(sum);
    }

    ret_val = value;
    return ret_val;
}


