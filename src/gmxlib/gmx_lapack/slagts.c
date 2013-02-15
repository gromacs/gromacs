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
F77_FUNC(slagts,SLAGTS)(int *job, 
	int *n, 
	float *a, 
	float *b, 
	float *c__, 
	float *d__, 
	int *in, 
	float *y, 
	float *tol, 
	int *info)
{
    int i__1;
    float d__1, d__2, d__4, d__5;

    int k;
    float ak, eps, temp, pert, absak, sfmin;
    float bignum,minval;
    --y;
    --in;
    --d__;
    --c__;
    --b;
    --a;

    *info = 0;
    if (fabs(*job) > 2 || *job == 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    }
    if (*info != 0) {
	i__1 = -(*info);
	return;
    }

    if (*n == 0) {
	return;
    }
    eps = GMX_FLOAT_EPS;
    minval = GMX_FLOAT_MIN;
    sfmin = minval / eps;

    bignum = 1. / sfmin;

    if (*job < 0) {
	if (*tol <= 0.) {
	    *tol = fabs(a[1]);
	    if (*n > 1) {
		d__1 = *tol;
		d__2 = fabs(a[2]);
		d__1 = (d__1>d__2) ? d__1 : d__2;
		d__2 = fabs(b[1]);
		*tol = (d__1>d__2) ? d__1 : d__2;
	    }
	    i__1 = *n;
	    for (k = 3; k <= i__1; ++k) {
	      d__4 = *tol;
	      d__5 = fabs(a[k]);
	      d__4 = (d__4>d__5) ? d__4 : d__5;
	      d__5 = fabs(b[k - 1]);
	      d__4 = (d__4>d__5) ? d__4 : d__5;
	      d__5 = fabs(d__[k - 2]);
	      *tol = (d__4>d__5) ? d__4 : d__5;
	    }
	    *tol *= eps;
	    if (fabs(*tol)<GMX_FLOAT_MIN) {
		*tol = eps;
	    }
	}
    }

    if (fabs(fabs(*job)-1.0)<GMX_FLOAT_MIN) {
	i__1 = *n;
	for (k = 2; k <= i__1; ++k) {
	    if (in[k - 1] == 0) {
		y[k] -= c__[k - 1] * y[k - 1];
	    } else {
		temp = y[k - 1];
		y[k - 1] = y[k];
		y[k] = temp - c__[k - 1] * y[k];
	    }
	}
	if (*job == 1) {
	    for (k = *n; k >= 1; --k) {
		if (k <= *n - 2) {
		    temp = y[k] - b[k] * y[k + 1] - d__[k] * y[k + 2];
		} else if (k == *n - 1) {
		    temp = y[k] - b[k] * y[k + 1];
		} else {
		    temp = y[k];
		}
		ak = a[k];
		absak = fabs(ak);
		if (absak < 1.) {
		    if (absak < sfmin) {
			if (fabs(absak)<GMX_FLOAT_MIN || fabs(temp) * sfmin > absak) {
			    *info = k;
			    return;
			} else {
			    temp *= bignum;
			    ak *= bignum;
			}
		    } else if (fabs(temp) > absak * bignum) {
			*info = k;
			return;
		    }
		}
		y[k] = temp / ak;
	    }
	} else {
	    for (k = *n; k >= 1; --k) {
		if (k <= *n - 2) {
		    temp = y[k] - b[k] * y[k + 1] - d__[k] * y[k + 2];
		} else if (k == *n - 1) {
		    temp = y[k] - b[k] * y[k + 1];
		} else {
		    temp = y[k];
		}
		ak = a[k];

		pert = *tol;
		if(ak<0)
		  pert *= -1.0;
L40:
		absak = fabs(ak);
		if (absak < 1.) {
		    if (absak < sfmin) {
			if (fabs(absak)<GMX_FLOAT_MIN || fabs(temp) * sfmin > absak) {
			    ak += pert;
			    pert *= 2;
			    goto L40;
			} else {
			    temp *= bignum;
			    ak *= bignum;
			}
		    } else if (fabs(temp) > absak * bignum) {
			ak += pert;
			pert *= 2;
			goto L40;
		    }
		}
		y[k] = temp / ak;
	    }
	}
    } else {

	if (*job == 2) {
	    i__1 = *n;
	    for (k = 1; k <= i__1; ++k) {
		if (k >= 3) {
		    temp = y[k] - b[k - 1] * y[k - 1] - d__[k - 2] * y[k - 2];
		} else if (k == 2) {
		    temp = y[k] - b[k - 1] * y[k - 1];
		} else {
		    temp = y[k];
		}
		ak = a[k];
		absak = fabs(ak);
		if (absak < 1.) {
		    if (absak < sfmin) {
			if (fabs(absak)<GMX_FLOAT_MIN || fabs(temp) * sfmin > absak) {
			    *info = k;
			    return;
			} else {
			    temp *= bignum;
			    ak *= bignum;
			}
		    } else if (fabs(temp) > absak * bignum) {
			*info = k;
			return;
		    }
		}
		y[k] = temp / ak;
	    }
	} else {
	    i__1 = *n;
	    for (k = 1; k <= i__1; ++k) {
		if (k >= 3) {
		    temp = y[k] - b[k - 1] * y[k - 1] - d__[k - 2] * y[k - 2];
		} else if (k == 2) {
		    temp = y[k] - b[k - 1] * y[k - 1];
		} else {
		    temp = y[k];
		}
		ak = a[k];

		pert = *tol;
		if(ak<0)
		  pert *= -1.0;

L70:
		absak = fabs(ak);
		if (absak < 1.) {
		    if (absak < sfmin) {
			if (fabs(absak)<GMX_FLOAT_MIN || fabs(temp) * sfmin > absak) {
			    ak += pert;
			    pert *= 2;
			    goto L70;
			} else {
			    temp *= bignum;
			    ak *= bignum;
			}
		    } else if (fabs(temp) > absak * bignum) {
			ak += pert;
			pert *= 2;
			goto L70;
		    }
		}
		y[k] = temp / ak;
	    }
	}

	for (k = *n; k >= 2; --k) {
	    if (in[k - 1] == 0) {
		y[k - 1] -= c__[k - 1] * y[k];
	    } else {
		temp = y[k - 1];
		y[k - 1] = y[k];
		y[k] = temp - c__[k - 1] * y[k];
	    }
	}
    }

    return;
}


