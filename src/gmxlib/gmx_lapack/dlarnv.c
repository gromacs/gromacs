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
F77_FUNC(dlarnv,DLARNV)(int *idist, 
	int *iseed, 
	int *n, 
	double *x)
{
    int i__1, i__2, i__3;

    int i__;
    double u[128];
    int il, iv, il2;

    --x;
    --iseed;

    i__1 = *n;
    for (iv = 1; iv <= i__1; iv += 64) {
	i__2 = 64, i__3 = *n - iv + 1;
	il = (i__2<i__3) ? i__2 : i__3;
	if (*idist == 3) {
	    il2 = il << 1;
	} else {
	    il2 = il;
	}

	F77_FUNC(dlaruv,DLARUV)(&iseed[1], &il2, u);

	if (*idist == 1) {

	    i__2 = il;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		x[iv + i__ - 1] = u[i__ - 1];
	    }
	} else if (*idist == 2) {

	    i__2 = il;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		x[iv + i__ - 1] = u[i__ - 1] * 2. - 1.;
	    }
	} else if (*idist == 3) {

	    i__2 = il;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		x[iv + i__ - 1] = sqrt(log(u[(i__ << 1) - 2]) * -2.) * 
		  cos(u[(i__ << 1) - 1] * 
		      (double)6.2831853071795864769252867663);
	    }
	}
    }
    return;

}
