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
F77_FUNC(slasq3,SLASQ3)(int *i0, 
                        int *n0, 
                        float *z__, 
                        int *pp, 
                        float *dmin__, 
                        float *sigma,
                        float *desig,
                        float *qmax, 
                        int *nfail, 
                        int *iter, 
                        int *ndiv, 
	int *ieee)
{

    int ttype = 0;
    float dmin1 = 0.;
    float dmin2 = 0.;
    float dn = 0.;
    float dn1 = 0.;
    float dn2 = 0.;
    float tau = 0.;

    int i__1;
    float d__1, d__2;
    float s, t;
    int j4, nn;
    float eps, tol;
    int n0in, ipn4;
    float tol2, temp;
    --z__;

    n0in = *n0;
    eps = GMX_FLOAT_EPS;
    tol = eps * 100.;
    d__1 = tol;
    tol2 = d__1 * d__1;


L10:

    if (*n0 < *i0) {
	return;
    }
    if (*n0 == *i0) {
	goto L20;
    }
    nn = (*n0 << 2) + *pp;
    if (*n0 == *i0 + 1) {
	goto L40;
    }

    if (z__[nn - 5] > tol2 * (*sigma + z__[nn - 3]) && z__[nn - (*pp << 1) - 
	    4] > tol2 * z__[nn - 7]) {
	goto L30;
    }

L20:

    z__[(*n0 << 2) - 3] = z__[(*n0 << 2) + *pp - 3] + *sigma;
    --(*n0);
    goto L10;

L30:

    if (z__[nn - 9] > tol2 * *sigma && z__[nn - (*pp << 1) - 8] > tol2 * z__[
	    nn - 11]) {
	goto L50;
    }

L40:

    if (z__[nn - 3] > z__[nn - 7]) {
	s = z__[nn - 3];
	z__[nn - 3] = z__[nn - 7];
	z__[nn - 7] = s;
    }
    if (z__[nn - 5] > z__[nn - 3] * tol2) {
	t = (z__[nn - 7] - z__[nn - 3] + z__[nn - 5]) * .5;
	s = z__[nn - 3] * (z__[nn - 5] / t);
	if (s <= t) {
	    s = z__[nn - 3] * (z__[nn - 5] / (t * (sqrt(s / t + 1.) + 1.)));
	} else {
	    s = z__[nn - 3] * (z__[nn - 5] / (t + sqrt(t) * sqrt(t + s)));
	}
	t = z__[nn - 7] + (s + z__[nn - 5]);
	z__[nn - 3] *= z__[nn - 7] / t;
	z__[nn - 7] = t;
    }
    z__[(*n0 << 2) - 7] = z__[nn - 7] + *sigma;
    z__[(*n0 << 2) - 3] = z__[nn - 3] + *sigma;
    *n0 += -2;
    goto L10;

L50:
    if (*pp == 2) {
	*pp = 0;
    }

    if (*dmin__ <= 0. || *n0 < n0in) {
	if (z__[(*i0 << 2) + *pp - 3] * 1.5 < z__[(*n0 << 2) + *pp - 3]) {
	    ipn4 = 4*(*i0 + *n0);
	    i__1 = 2*(*i0 + *n0 - 1);
	    for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
		temp = z__[j4 - 3];
		z__[j4 - 3] = z__[ipn4 - j4 - 3];
		z__[ipn4 - j4 - 3] = temp;
		temp = z__[j4 - 2];
		z__[j4 - 2] = z__[ipn4 - j4 - 2];
		z__[ipn4 - j4 - 2] = temp;
		temp = z__[j4 - 1];
		z__[j4 - 1] = z__[ipn4 - j4 - 5];
		z__[ipn4 - j4 - 5] = temp;
		temp = z__[j4];
		z__[j4] = z__[ipn4 - j4 - 4];
		z__[ipn4 - j4 - 4] = temp;
	    }
	    if (*n0 - *i0 <= 4) {
		z__[(*n0 << 2) + *pp - 1] = z__[(*i0 << 2) + *pp - 1];
		z__[(*n0 << 2) - *pp] = z__[(*i0 << 2) - *pp];
	    }
	    d__1 = dmin2, d__2 = z__[(*n0 << 2) + *pp - 1];
	    dmin2 = ((d__1<d__2) ? d__1 : d__2);
	    d__1 = z__[(*n0 << 2) + *pp - 1], d__2 = z__[(*i0 << 2) + *pp - 1]
		    , d__1 = ((d__1<d__2) ? d__1 : d__2), d__2 = z__[(*i0 << 2) + *pp + 3];
	    z__[(*n0 << 2) + *pp - 1] = ((d__1<d__2) ? d__1 : d__2);
	    d__1 = z__[(*n0 << 2) - *pp], d__2 = z__[(*i0 << 2) - *pp], d__1 =
		     ((d__1<d__2) ? d__1 : d__2), d__2 = z__[(*i0 << 2) - *pp + 4];
	    z__[(*n0 << 2) - *pp] = ((d__1<d__2) ? d__1 : d__2);
	    d__1 = *qmax;
	    d__2 = z__[(*i0 << 2) + *pp - 3];
	    d__1 = (d__1>d__2) ? d__1 : d__2;
	    d__2 = z__[(*i0 << 2) + *pp + 1];
	    *qmax = ((d__1>d__2) ? d__1 : d__2);
	    *dmin__ = -0.;
	}
    }


    F77_FUNC(slasq4,SLASQ4)(i0, n0, &z__[1], pp, &n0in, dmin__, &dmin1, &dmin2, &dn, &dn1, &
	    dn2, &tau, &ttype);

L70:

    F77_FUNC(slasq5,SLASQ5)(i0, n0, &z__[1], pp, &tau, dmin__, &dmin1, &dmin2, &dn, &dn1, &
	    dn2, ieee);

    *ndiv += *n0 - *i0 + 2;
    ++(*iter);

    if (*dmin__ >= 0. && dmin1 > 0.) {

	goto L90;

    } else if (*dmin__ < 0. && dmin1 > 0. && z__[4*(*n0 - 1) - *pp] < tol *
	     (*sigma + dn1) && fabs(dn) < tol * *sigma) {

	z__[4*(*n0 - 1) - *pp + 2] = 0.;
	*dmin__ = 0.;
	goto L90;
    } else if (*dmin__ < 0.) {

	++(*nfail);
	if (ttype < -22) {

	    tau = 0.;
	} else if (dmin1 > 0.) {

	    tau = (tau + *dmin__) * (1. - eps * 2.);
	    ttype += -11;
	} else {

	    tau *= .25;
	    ttype += -12;
	}
	goto L70;
    }
    else {
        
        goto L80;
    }

L80:
    F77_FUNC(slasq6,SLASQ6)(i0, n0, &z__[1], pp, dmin__, &dmin1, &dmin2, &dn, &dn1, &dn2);
    *ndiv += *n0 - *i0 + 2;
    ++(*iter);
    tau = 0.;

L90:
    if (tau < *sigma) {
	*desig += tau;
	t = *sigma + *desig;
	*desig -= t - *sigma;
    } else {
	t = *sigma + tau;
	*desig = *sigma - (t - tau) + *desig;
    }
    *sigma = t;

    return;
}
