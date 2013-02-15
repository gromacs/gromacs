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


void
F77_FUNC(dlaev2,DLAEV2)(double *   a, 
	double *   b, 
	double *   c__, 
	double *   rt1, 
	double *   rt2, 
	double *   cs1, 
	double *   sn1)
{
    double d__1;

    double ab, df, cs, ct, tb, sm, tn, rt, adf, acs;
    int sgn1, sgn2;
    double acmn, acmx;

    sm = *a + *c__;
    df = *a - *c__;
    adf = fabs(df);
    tb = *b + *b;
    ab = fabs(tb);
    if (fabs(*a) > fabs(*c__)) {
	acmx = *a;
	acmn = *c__;
    } else {
	acmx = *c__;
	acmn = *a;
    }
    if (adf > ab) {
	d__1 = ab / adf;
	rt = adf * sqrt(d__1 * d__1 + 1.);
    } else if (adf < ab) {
	d__1 = adf / ab;
	rt = ab * sqrt(d__1 * d__1 + 1.);
    } else {

	rt = ab * sqrt(2.);
    }
    if (sm < 0.) {
	*rt1 = (sm - rt) * .5;
	sgn1 = -1;

	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else if (sm > 0.) {
	*rt1 = (sm + rt) * .5;
	sgn1 = 1;
	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else {
	*rt1 = rt * .5;
	*rt2 = rt * -.5;
	sgn1 = 1;
    }
    if (df >= 0.) {
	cs = df + rt;
	sgn2 = 1;
    } else {
	cs = df - rt;
	sgn2 = -1;
    }
    acs = fabs(cs);
    if (acs > ab) {
	ct = -tb / cs;
	*sn1 = 1. / sqrt(ct * ct + 1.);
	*cs1 = ct * *sn1;
    } else {
	if (fabs(ab)<GMX_DOUBLE_MIN) {
	    *cs1 = 1.;
	    *sn1 = 0.;
	} else {
	    tn = -cs / tb;
	    *cs1 = 1. / sqrt(tn * tn + 1.);
	    *sn1 = tn * *cs1;
	}
    }
    if (sgn1 == sgn2) {
	tn = *cs1;
	*cs1 = -(*sn1);
	*sn1 = tn;
    }
    return;

}


