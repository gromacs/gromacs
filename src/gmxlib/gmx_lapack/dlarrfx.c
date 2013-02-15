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
#include "lapack_limits.h"


void
F77_FUNC(dlarrfx,DLARRFX)(int *n, 
	double *d__, 
	double *l, 
	double *ld, 
	double *lld, 
	int *ifirst, 
	int *ilast, 
	double *w, 
	double *sigma, 
	double *dplus, 
	double *lplus, 
	double *work,
	int *info)
{
    int i1 = 1;
    int i__1;
    double d__2, d__3;

    int i__;
    double s, eps, tmp, dmax1, dmax2, delta;
    --work;
    --lplus;
    --dplus;
    --w;
    --lld;
    --ld;
    --l;
    --d__;
    *info = 0;
    eps = GMX_DOUBLE_EPS;
    *sigma = w[*ifirst];
    delta = eps * 2.;

L10:
    s = -(*sigma);
    dplus[1] = d__[1] + s;
    dmax1 = fabs(dplus[1]);
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	lplus[i__] = ld[i__] / dplus[i__];
	s = s * lplus[i__] * l[i__] - *sigma;
	dplus[i__ + 1] = d__[i__ + 1] + s;
	d__2 = dmax1, d__3 = fabs(dplus[i__ + 1]);
	dmax1 = (d__2>d__3) ? d__2 : d__3;
    }
    if (! (dmax1 > 0. || dmax1 < 1.)) {
	*sigma -= fabs(*sigma) * delta;
	delta *= 2.;
	goto L10;
    }

    tmp = w[*ilast];
    delta = eps * 2.;
L30:
    s = -tmp;
    work[1] = d__[1] + s;
    dmax2 = fabs(work[1]);
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	work[*n + i__] = ld[i__] / work[i__];
	s = s * work[*n + i__] * l[i__] - tmp;
	work[i__ + 1] = d__[i__ + 1] + s;
	d__2 = dmax2, d__3 = fabs(work[i__ + 1]);
	dmax2 = (d__2>d__3) ? d__2 : d__3;
    }
    if (! (dmax2 > 0. || dmax2 < 1.)) {
	tmp += fabs(tmp) * delta;
	delta *= 2.;
	goto L30;
    }
    if (dmax2 < dmax1) {
	*sigma = tmp;
	F77_FUNC(dcopy,DCOPY)(n, &work[1], &i1, &dplus[1], &i1);
	i__1 = *n - 1;
	F77_FUNC(dcopy,DCOPY)(&i__1, &work[*n + 1], &i1, &lplus[1], &i1);
    }

    return;
}
