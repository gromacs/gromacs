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
#include "gmx_blas.h"
#include "gmx_lapack.h"
#include "lapack_limits.h"

#include <types/simple.h>

void
F77_FUNC(dstegr,DSTEGR)(const char *jobz, 
	const char *range, 
	int *n, 
	double *d__, 
	double *e, 
	double *vl, 
	double *vu, 
	int *il, 
	int *iu, 
	double *abstol, 
	int *m, 
	double *w, 
	double *z__, 
	int *ldz, 
	int *isuppz,
	double *work, 
	int *lwork, 
	int *iwork, 
	int *liwork, 
	int *info)
{
    int z_dim1, z_offset, i__1, i__2;
    double d__1, d__2;
    int c__1 = 1;

    int i__, j;
    int jj;
    double eps, tol, tmp, rmin, rmax;
    int itmp;
    double tnrm;
    double scale;
    int iinfo, iindw;
    int lwmin;
    int wantz;
    int iindbl;
    int valeig,alleig,indeig;
    double safmin,minval;
    double bignum;
    int iindwk, indgrs;
    double thresh;
    int iinspl, indwrk, liwmin, nsplit;
    double smlnum;
    int lquery;


    --d__;
    --e;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --isuppz;
    --work;
    --iwork;

    wantz = (*jobz=='V' || *jobz=='v');
    alleig = (*range=='A' || *range=='a');
    valeig = (*range=='V' || *range=='v');
    indeig = (*range=='I' || *range=='i');

    lquery = *lwork == -1 || *liwork == -1;
    lwmin = *n * 17;
    liwmin = *n * 10;

    *info = 0;
    if (! (wantz || (*jobz=='N' || *jobz=='n'))) {
	*info = -1;
    } else if (! (alleig || valeig || indeig)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (valeig && *n > 0 && *vu <= *vl) {
	*info = -7;
    } else if (indeig && (*il < 1 || *il > *n)) {
	*info = -8;
    } else if (indeig && (*iu < *il || *iu > *n)) {
	*info = -9;
    } else if (*ldz < 1 || (wantz && *ldz < *n)) {
	*info = -14;
    } else if (*lwork < lwmin && ! lquery) {
	*info = -17;
    } else if (*liwork < liwmin && ! lquery) {
	*info = -19;
    }
    if (*info == 0) {
	work[1] = (double) lwmin;
	iwork[1] = liwmin;
    }

    if (*info != 0) {
	i__1 = -(*info);
	return;
    } else if (lquery) {
	return;
    }

    *m = 0;
    if (*n == 0) {
	return;
    }

    if (*n == 1) {
	if (alleig || indeig) {
	    *m = 1;
	    w[1] = d__[1];
	} else {
	    if (*vl < d__[1] && *vu >= d__[1]) {
		*m = 1;
		w[1] = d__[1];
	    }
	}
	if (wantz) {
	    z__[z_dim1 + 1] = 1.;
	}
	return;
    }

    minval = GMX_DOUBLE_MIN;
    safmin = minval*(1.0+GMX_DOUBLE_EPS);
    eps = GMX_DOUBLE_EPS;
    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin = sqrt(smlnum);
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
    rmax = (d__1<d__2) ? d__1 : d__2;
    scale = 1.;
    tnrm = F77_FUNC(dlanst,DLANST)("M", n, &d__[1], &e[1]);
    if (tnrm > 0. && tnrm < rmin) {
	scale = rmin / tnrm;
    } else if (tnrm > rmax) {
	scale = rmax / tnrm;
    }
    if ( fabs(scale-1.0)>GMX_DOUBLE_EPS) {
	F77_FUNC(dscal,DSCAL)(n, &scale, &d__[1], &c__1);
	i__1 = *n - 1;
	F77_FUNC(dscal,DSCAL)(&i__1, &scale, &e[1], &c__1);
	tnrm *= scale;
    }
    indgrs = 1;
    indwrk = (*n << 1) + 1;

    iinspl = 1;
    iindbl = *n + 1;
    iindw = (*n << 1) + 1;
    iindwk = *n * 3 + 1;

    thresh = eps * tnrm;
    F77_FUNC(dlarrex,DLARREX)(range, n, vl, vu, il, iu, &d__[1], &e[1], &thresh, &nsplit, &
	    iwork[iinspl], m, &w[1], &iwork[iindbl], &iwork[iindw], &work[
	    indgrs], &work[indwrk], &iwork[iindwk], &iinfo);
    
    if (iinfo != 0) {
	*info = 1;
	return;
    }

    if (wantz) {
	d__1 = *abstol, d__2 = (double) (*n) * eps;
	tol = (d__1>d__2) ? d__1 : d__2;
	F77_FUNC(dlarrvx,DLARRVX)(n, &d__[1], &e[1], &iwork[iinspl], m, &w[1], &iwork[iindbl], &
		iwork[iindw], &work[indgrs], &tol, &z__[z_offset], ldz, &
		isuppz[1], &work[indwrk], &iwork[iindwk], &iinfo);
	if (iinfo != 0) {
	    *info = 2;
	    return;
	}
    }

    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	itmp = iwork[iindbl + j - 1];
	w[j] += e[iwork[iinspl + itmp - 1]];
    } 

    if (fabs(scale-1.0)>GMX_DOUBLE_EPS) {
	d__1 = 1. / scale;
	F77_FUNC(dscal,DSCAL)(m, &d__1, &w[1], &c__1);
    }
    if (nsplit > 1) {
	i__1 = *m - 1;
	for (j = 1; j <= i__1; ++j) {
	    i__ = 0;
	    tmp = w[j];
	    i__2 = *m;
	    for (jj = j + 1; jj <= i__2; ++jj) {
		if (w[jj] < tmp) {
		    i__ = jj;
		    tmp = w[jj];
		}
	    }
	    if (i__ != 0) {
		w[i__] = w[j];
		w[j] = tmp;
		if (wantz) {
		    F77_FUNC(dswap,DSWAP)(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 
			    + 1], &c__1);
		    itmp = isuppz[(i__ << 1) - 1];
		    isuppz[(i__ << 1) - 1] = isuppz[(j << 1) - 1];
		    isuppz[(j << 1) - 1] = itmp;
		    itmp = isuppz[i__ * 2];
		    isuppz[i__ * 2] = isuppz[j * 2];
		    isuppz[j * 2] = itmp;
		}
	    }
	}
    }

    work[1] = (double) lwmin;
    iwork[1] = liwmin;
    return;

} 
