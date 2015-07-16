#include <cmath>
#include "gromacs/utility/real.h"

#include "../gmx_lapack.h"
#include "lapack_limits.h"



void
F77_FUNC(dlagtf,DLAGTF)(int *n, 
	double *a, 
	double *lambda, 
	double *b, 
	double *c__, 
	double *tol, 
	double *d__, 
	int *in, 
	int *info)
{
    int i__1;

    int k;
    double tl, eps, piv1, piv2, temp, mult, scale1, scale2;

    --in;
    --d__;
    --c__;
    --b;
    --a;

    *info = 0;
    if (*n < 0) {
	*info = -1;
	return;
    }

    if (*n == 0) 
	return;
    
    a[1] -= *lambda;
    in[*n] = 0;
    if (*n == 1) {
	if (std::abs(a[1])<GMX_DOUBLE_MIN) {
	    in[1] = 1;
	}
	return;
    }

    eps = GMX_DOUBLE_EPS;

    tl = (*tol>eps) ? *tol : eps;
    scale1 = std::abs(a[1]) + std::abs(b[1]);
    i__1 = *n - 1;
    for (k = 1; k <= i__1; ++k) {
	a[k + 1] -= *lambda;
	scale2 = std::abs(c__[k]) + std::abs(a[k + 1]);
	if (k < *n - 1) {
	    scale2 += std::abs(b[k + 1]);
	}
	if (std::abs(a[k])<GMX_DOUBLE_MIN) {
	    piv1 = 0.;
	} else {
	    piv1 = std::abs(a[k]) / scale1;
	}
	if (std::abs(c__[k])<GMX_DOUBLE_MIN) {
	    in[k] = 0;
	    piv2 = 0.;
	    scale1 = scale2;
	    if (k < *n - 1) {
		d__[k] = 0.;
	    }
	} else {
	    piv2 = std::abs(c__[k]) / scale2;
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
    if (std::abs(a[*n]) <= scale1 * tl && in[*n] == 0) {
	in[*n] = *n;
    }

    return;

}


