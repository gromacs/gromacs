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

void
F77_FUNC(slasq5,SLASQ5)(int *i0, 
	int *n0,
	float *z__, 
	int *pp, 
	float *tau,
	float *dmin__, 
	float *dmin1, 
	float *dmin2, 
	float *dn,
	float *dnm1, 
	float *dnm2,
	int *ieee)
{
    int i__1;
    float d__1, d__2;

    static float d__;
    static int j4, j4p2;
    static float emin, temp;

    --z__;

    if (*n0 - *i0 - 1 <= 0) {
	return;
    }

    j4 = (*i0 << 2) + *pp - 3;
    emin = z__[j4 + 4];
    d__ = z__[j4] - *tau;
    *dmin__ = d__;
    *dmin1 = -z__[j4];

    if (*ieee) {

	if (*pp == 0) {
	    i__1 = 4*(*n0 - 3);
	    for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
		z__[j4 - 2] = d__ + z__[j4 - 1];
		temp = z__[j4 + 1] / z__[j4 - 2];
		d__ = d__ * temp - *tau;
                if(d__<*dmin__)
                  *dmin__ = d__;
		z__[j4] = z__[j4 - 1] * temp;
		d__1 = z__[j4];
                if(d__1<emin)
                  emin = d__1;
	    }
	} else {
	    i__1 = 4*(*n0 - 3);
	    for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
		z__[j4 - 3] = d__ + z__[j4];
		temp = z__[j4 + 2] / z__[j4 - 3];
		d__ = d__ * temp - *tau;
                if(d__<*dmin__)
                  *dmin__ = d__;
		z__[j4 - 1] = z__[j4] * temp;
		d__1 = z__[j4 - 1];
                if(d__1<emin)
                  emin = d__1;
	    }
	}

	*dnm2 = d__;
	*dmin2 = *dmin__;
	j4 = 4*(*n0 - 2) - *pp;
	j4p2 = j4 + (*pp << 1) - 1;
	z__[j4 - 2] = *dnm2 + z__[j4p2];
	z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
	*dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]) - *tau;
        if(*dnm1<*dmin__)
          *dmin__ = *dnm1;

	*dmin1 = *dmin__;
	j4 += 4;
	j4p2 = j4 + (*pp << 1) - 1;
	z__[j4 - 2] = *dnm1 + z__[j4p2];
	z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
	*dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]) - *tau;
        if(*dn<*dmin__)
          *dmin__ = *dn;

    } else {

	if (*pp == 0) {
	    i__1 = 4*(*n0 - 3);
	    for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
		z__[j4 - 2] = d__ + z__[j4 - 1];
		if (d__ < 0.) {
		    return;
		} else {
		    z__[j4] = z__[j4 + 1] * (z__[j4 - 1] / z__[j4 - 2]);
		    d__ = z__[j4 + 1] * (d__ / z__[j4 - 2]) - *tau;
		}
                if(d__<*dmin__)
                  *dmin__ = d__;
		d__1 = emin, d__2 = z__[j4];
		emin = (d__1<d__2) ? d__1 : d__2;
	    }
	} else {
	    i__1 = 4*(*n0 - 3);
	    for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
		z__[j4 - 3] = d__ + z__[j4];
		if (d__ < 0.) {
		    return;
		} else {
		    z__[j4 - 1] = z__[j4 + 2] * (z__[j4] / z__[j4 - 3]);
		    d__ = z__[j4 + 2] * (d__ / z__[j4 - 3]) - *tau;
		}
                if(d__<*dmin__)
                  *dmin__ = d__;
		d__1 = emin, d__2 = z__[j4 - 1];
		emin = (d__1<d__2) ? d__1 : d__2;
	    }
	}

	*dnm2 = d__;
	*dmin2 = *dmin__;
	j4 = 4*(*n0 - 2) - *pp;
	j4p2 = j4 + (*pp << 1) - 1;
	z__[j4 - 2] = *dnm2 + z__[j4p2];
	if (*dnm2 < 0.) {
	    return;
	} else {
	    z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
	    *dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]) - *tau;
	}
        if(*dnm1<*dmin__)
          *dmin__ = *dnm1;

	*dmin1 = *dmin__;
	j4 += 4;
	j4p2 = j4 + (*pp << 1) - 1;
	z__[j4 - 2] = *dnm1 + z__[j4p2];
	if (*dnm1 < 0.) {
	    return;
	} else {
	    z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
	    *dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]) - *tau;
	}
        if(*dn<*dmin__)
          *dmin__ = *dn;

    }

    z__[j4 + 2] = *dn;
    z__[(*n0 << 2) - *pp] = emin;
    return;

}

