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

#include "gmx_lapack.h"
#include "lapack_limits.h"



void
F77_FUNC(slagtf,SLAGTF)(int *n, 
	float *a, 
	float *lambda, 
	float *b, 
	float *c__, 
	float *tol, 
	float *d__, 
	int *in, 
	int *info)
{
    int i__1;

    int k;
    float tl, eps, piv1, piv2, temp, mult, scale1, scale2;

    --in;
    --d__;
    --c__;
    --b;
    --a;

    *info = 0;
    if (*n < 0) {
	*info = -1;
	i__1 = -(*info);
	return;
    }

    if (*n == 0) 
	return;
    
    a[1] -= *lambda;
    in[*n] = 0;
    if (*n == 1) {
	if (fabs(a[1])<GMX_FLOAT_MIN) {
	    in[1] = 1;
	}
	return;
    }

    eps = GMX_FLOAT_EPS;

    tl = (*tol>eps) ? *tol : eps;
    scale1 = fabs(a[1]) + fabs(b[1]);
    i__1 = *n - 1;
    for (k = 1; k <= i__1; ++k) {
	a[k + 1] -= *lambda;
	scale2 = fabs(c__[k]) + fabs(a[k + 1]);
	if (k < *n - 1) {
	    scale2 += fabs(b[k + 1]);
	}
	if (fabs(a[k])<GMX_FLOAT_MIN) {
	    piv1 = 0.;
	} else {
	    piv1 = fabs(a[k]) / scale1;
	}
	if (fabs(c__[k])<GMX_FLOAT_MIN) {
	    in[k] = 0;
	    piv2 = 0.;
	    scale1 = scale2;
	    if (k < *n - 1) {
		d__[k] = 0.;
	    }
	} else {
	    piv2 = fabs(c__[k]) / scale2;
	    if (piv2 <= piv1) {
		in[k] = 0;
		scale1 = scale2;
		c__[k] /= a[k];
		a[k + 1] -= c__[k] * b[k];
		if (k < *n - 1) {
		    d__[k] = 0.;
		}
	    } else {
		in[k] = 1;
		mult = a[k] / c__[k];
		a[k] = c__[k];
		temp = a[k + 1];
		a[k + 1] = b[k] - mult * temp;
		if (k < *n - 1) {
		    d__[k] = b[k + 1];
		    b[k + 1] = -mult * d__[k];
		}
		b[k] = temp;
		c__[k] = mult;
	    }
	}
	if (((piv1>piv2) ? piv1 : piv2) <= tl && in[*n] == 0) {
	    in[*n] = k;
	}
    }
    if (fabs(a[*n]) <= scale1 * tl && in[*n] == 0) {
	in[*n] = *n;
    }

    return;

}


