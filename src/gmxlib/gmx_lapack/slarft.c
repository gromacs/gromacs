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
#include <types/simple.h>

#include "gmx_blas.h"
#include "gmx_lapack.h"

void 
F77_FUNC(slarft,SLARFT)(const char *direct, 
	const char *storev, 
	int *n, 
	int *k, 
	float *v, 
	int *ldv, 
	float *tau, 
	float *t, 
	int *ldt)
{
    /* System generated locals */
    int t_dim1, t_offset, v_dim1, v_offset, i__1, i__2, i__3;
    float d__1;

    /* Local variables */
    int i__, j;
    float vii;
    int c__1 = 1;
    float zero = 0.0;

    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --tau;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;

    if (*n == 0) {
	return;
    }

    if (*direct=='F' || *direct=='f') {
	i__1 = *k;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (fabs(tau[i__])<GMX_FLOAT_MIN) {

		i__2 = i__;
		for (j = 1; j <= i__2; ++j) {
		    t[j + i__ * t_dim1] = 0.;
		}
	    } else {

		vii = v[i__ + i__ * v_dim1];
		v[i__ + i__ * v_dim1] = 1.;
		if (*storev=='C' || *storev=='c') {

		    i__2 = *n - i__ + 1;
		    i__3 = i__ - 1;
		    d__1 = -tau[i__];
		    F77_FUNC(sgemv,SGEMV)("Transpose", &i__2, &i__3, &d__1, &v[i__ + v_dim1],
			     ldv, &v[i__ + i__ * v_dim1], &c__1, &zero, &t[
			    i__ * t_dim1 + 1], &c__1);
		} else {

		    i__2 = i__ - 1;
		    i__3 = *n - i__ + 1;
		    d__1 = -tau[i__];
		    F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__3, &d__1, &v[i__ * 
			    v_dim1 + 1], ldv, &v[i__ + i__ * v_dim1], ldv, &
			    zero, &t[i__ * t_dim1 + 1], &c__1);
		}
		v[i__ + i__ * v_dim1] = vii;


		i__2 = i__ - 1;
		F77_FUNC(strmv,STRMV)("Upper", "No transpose", "Non-unit", &i__2, &t[
			t_offset], ldt, &t[i__ * t_dim1 + 1], &c__1);
		t[i__ + i__ * t_dim1] = tau[i__];
	    }
	}
    } else {
	for (i__ = *k; i__ >= 1; --i__) {
	    if (fabs(tau[i__])<GMX_FLOAT_MIN) {

		i__1 = *k;
		for (j = i__; j <= i__1; ++j) {
		    t[j + i__ * t_dim1] = 0.;
		}
	    } else {

		if (i__ < *k) {
		    if (*storev=='C' || *storev=='c') {
			vii = v[*n - *k + i__ + i__ * v_dim1];
			v[*n - *k + i__ + i__ * v_dim1] = 1.;

			i__1 = *n - *k + i__;
			i__2 = *k - i__;
			d__1 = -tau[i__];
			F77_FUNC(sgemv,SGEMV)("Transpose", &i__1, &i__2, &d__1, &v[(i__ + 1) 
				* v_dim1 + 1], ldv, &v[i__ * v_dim1 + 1], &
				c__1, &zero, &t[i__ + 1 + i__ * t_dim1], &
				c__1);
			v[*n - *k + i__ + i__ * v_dim1] = vii;
		    } else {
			vii = v[i__ + (*n - *k + i__) * v_dim1];
			v[i__ + (*n - *k + i__) * v_dim1] = 1.;

			i__1 = *k - i__;
			i__2 = *n - *k + i__;
			d__1 = -tau[i__];
			F77_FUNC(sgemv,SGEMV)("No transpose", &i__1, &i__2, &d__1, &v[i__ + 
				1 + v_dim1], ldv, &v[i__ + v_dim1], ldv, &
				zero, &t[i__ + 1 + i__ * t_dim1], &c__1);
			v[i__ + (*n - *k + i__) * v_dim1] = vii;
		    }

		    i__1 = *k - i__;
		    F77_FUNC(strmv,STRMV)("Lower", "No transpose", "Non-unit", &i__1, &t[i__ 
			    + 1 + (i__ + 1) * t_dim1], ldt, &t[i__ + 1 + i__ *
			     t_dim1], &c__1)
			    ;
		}
		t[i__ + i__ * t_dim1] = tau[i__];
	    }
	}
    }
    return;


}
