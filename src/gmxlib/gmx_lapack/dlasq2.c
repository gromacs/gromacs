#include <math.h>
#include "gmx_lapack.h"
#include "lapack_limits.h"

#include <types/simple.h>

void 
F77_FUNC(dlasq2,DLASQ2)(int *n, 
                        double *z__, 
                        int *info)
{
    int i__1, i__2, i__3;
    double d__1, d__2;

    double d__, e;
    int k;
    double s, t;
    int i0, i4, n0, pp;
    double dee, eps, tol;
    int ipn4;
    double tol2;
    int ieee;
    int nbig;
    double dmin__, emin, emax;
    int kmin, ndiv, iter;
    double qmin, temp, qmax, zmax;
    int splt, nfail;
    double desig, trace, sigma;
    int iinfo;
    double deemin;
    int iwhila, iwhilb;
    double oldemn, safmin, minval;
    double posinf,neginf,negzro,newzro;
    double zero = 0.0;
    double one = 1.0;

    --z__;

    *info = 0;
    eps = GMX_DOUBLE_EPS;
    minval = GMX_DOUBLE_MIN;
    safmin = minval*(1.0+eps);

    tol = eps * 100.;

    d__1 = tol;
    tol2 = d__1 * d__1;

    if (*n < 0) {
	*info = -1;
	return;
    } else if (*n == 0) {
	return;
    } else if (*n == 1) {

	if (z__[1] < 0.) {
	    *info = -201;
	}
	return;
    } else if (*n == 2) {

	if (z__[2] < 0. || z__[3] < 0.) {
	    *info = -2;
	    return;
	} else if (z__[3] > z__[1]) {
	    d__ = z__[3];
	    z__[3] = z__[1];
	    z__[1] = d__;
	}
	z__[5] = z__[1] + z__[2] + z__[3];
	if (z__[2] > z__[3] * tol2) {
	    t = (z__[1] - z__[3] + z__[2]) * .5;
	    s = z__[3] * (z__[2] / t);
	    if (s <= t) {
		s = z__[3] * (z__[2] / (t * (sqrt(s / t + 1.) + 1.)));
	    } else {
		s = z__[3] * (z__[2] / (t + sqrt(t) * sqrt(t + s)));
	    }
	    t = z__[1] + (s + z__[2]);
	    z__[3] *= z__[1] / t;
	    z__[1] = t;
	}
	z__[2] = z__[3];
	z__[6] = z__[2] + z__[1];
	return;
    }
    z__[*n * 2] = 0.;
    emin = z__[2];
    qmax = 0.;
    zmax = 0.;
    d__ = 0.;
    e = 0.;

    i__1 = 2*(*n - 1);
    for (k = 1; k <= i__1; k += 2) {
	if (z__[k] < 0.) {
	    *info = -(k + 200);
	    return;
	} else if (z__[k + 1] < 0.) {
	    *info = -(k + 201);
	    return;
	}
	d__ += z__[k];
	e += z__[k + 1];
	d__1 = qmax, d__2 = z__[k];
	qmax = (d__1>d__2) ? d__1 : d__2;
	d__1 = emin, d__2 = z__[k + 1];
	emin = (d__1<d__2) ? d__1 : d__2;
	d__1 = (qmax>zmax) ? qmax : zmax;
	d__2 = z__[k + 1];
	zmax = (d__1>d__2) ? d__1 : d__2;
    }
    if (z__[(*n << 1) - 1] < 0.) {
	*info = -((*n << 1) + 199);
	return;
    }
    d__ += z__[(*n << 1) - 1];
    d__1 = qmax, d__2 = z__[(*n << 1) - 1];
    qmax = (d__1>d__2) ? d__1 : d__2;
    zmax = (qmax>zmax) ? qmax : zmax;

    if (fabs(e)<GMX_DOUBLE_MIN) {
	i__1 = *n;
	for (k = 2; k <= i__1; ++k) {
	    z__[k] = z__[(k << 1) - 1];
	}
	F77_FUNC(dlasrt,DLASRT)("D", n, &z__[1], &iinfo);
	z__[(*n << 1) - 1] = d__;
	return;
    }

    trace = d__ + e;

    if (fabs(trace)<GMX_DOUBLE_MIN) {
	z__[(*n << 1) - 1] = 0.;
	return;
    }

    ieee = 1;
    posinf = one/zero;
    if(posinf<=1.0)
      ieee = 0;
    neginf = -one/zero;
    if(neginf>=0.0)
      ieee = 0;
    negzro = one/(neginf+one);
    if(fabs(negzro)>GMX_DOUBLE_MIN)
      ieee = 0;
    neginf = one/negzro;
    if(neginf>=0)
      ieee = 0;
    newzro = negzro + zero;
    if(fabs(newzro-zero)>GMX_DOUBLE_MIN)
      ieee = 0;
    posinf = one /newzro;
    if(posinf<=one)
      ieee = 0;
    neginf = neginf*posinf;
    if(neginf>=zero)
      ieee = 0;
    posinf = posinf*posinf;
    if(posinf<=1.0)
      ieee = 0;

    for (k = *n << 1; k >= 2; k += -2) {
	z__[k * 2] = 0.;
	z__[(k << 1) - 1] = z__[k];
	z__[(k << 1) - 2] = 0.;
	z__[(k << 1) - 3] = z__[k - 1];
    }

    i0 = 1;
    n0 = *n;

    if (z__[(i0 << 2) - 3] * 1.5 < z__[(n0 << 2) - 3]) {
	ipn4 = 4*(i0 + n0);
	i__1 = 2*(i0 + n0 - 1);
	for (i4 = i0 << 2; i4 <= i__1; i4 += 4) {
	    temp = z__[i4 - 3];
	    z__[i4 - 3] = z__[ipn4 - i4 - 3];
	    z__[ipn4 - i4 - 3] = temp;
	    temp = z__[i4 - 1];
	    z__[i4 - 1] = z__[ipn4 - i4 - 5];
	    z__[ipn4 - i4 - 5] = temp;
	}
    }

    pp = 0;

    for (k = 1; k <= 2; ++k) {

	d__ = z__[(n0 << 2) + pp - 3];
	i__1 = (i0 << 2) + pp;
	for (i4 = 4*(n0 - 1) + pp; i4 >= i__1; i4 += -4) {
	    if (z__[i4 - 1] <= tol2 * d__) {
		z__[i4 - 1] = -0.;
		d__ = z__[i4 - 3];
	    } else {
		d__ = z__[i4 - 3] * (d__ / (d__ + z__[i4 - 1]));
	    }
	}

	emin = z__[(i0 << 2) + pp + 1];
	d__ = z__[(i0 << 2) + pp - 3];
	i__1 = 4*(n0 - 1) + pp;
	for (i4 = (i0 << 2) + pp; i4 <= i__1; i4 += 4) {
	    z__[i4 - (pp << 1) - 2] = d__ + z__[i4 - 1];
	    if (z__[i4 - 1] <= tol2 * d__) {
		z__[i4 - 1] = -0.;
		z__[i4 - (pp << 1) - 2] = d__;
		z__[i4 - (pp << 1)] = 0.;
		d__ = z__[i4 + 1];
	    } else if (safmin * z__[i4 + 1] < z__[i4 - (pp << 1) - 2] && 
		    safmin * z__[i4 - (pp << 1) - 2] < z__[i4 + 1]) {
		temp = z__[i4 + 1] / z__[i4 - (pp << 1) - 2];
		z__[i4 - (pp << 1)] = z__[i4 - 1] * temp;
		d__ *= temp;
	    } else {
		z__[i4 - (pp << 1)] = z__[i4 + 1] * (z__[i4 - 1] / z__[i4 - (
			pp << 1) - 2]);
		d__ = z__[i4 + 1] * (d__ / z__[i4 - (pp << 1) - 2]);
	    }
	    d__1 = emin, d__2 = z__[i4 - (pp << 1)];
	    emin = (d__1<d__2) ? d__1 : d__2;
	}
	z__[(n0 << 2) - pp - 2] = d__;


	qmax = z__[(i0 << 2) - pp - 2];
	i__1 = (n0 << 2) - pp - 2;
	for (i4 = (i0 << 2) - pp + 2; i4 <= i__1; i4 += 4) {
	    d__1 = qmax, d__2 = z__[i4];
	    qmax = (d__1>d__2) ? d__1 : d__2;
	}

	pp = 1 - pp;
    }

    iter = 2;
    nfail = 0;
    ndiv = 2*(n0 - i0);

    i__1 = *n + 1;
    for (iwhila = 1; iwhila <= i__1; ++iwhila) {
	if (n0 < 1) {
	    goto L170;
	}

	desig = 0.;
	if (n0 == *n) {
	    sigma = 0.;
	} else {
	    sigma = -z__[(n0 << 2) - 1];
	}
	if (sigma < 0.) {
	    *info = 1;
	    return;
	}

	emax = 0.;
	if (n0 > i0) {
	    emin = fabs(z__[(n0 << 2) - 5]);
	} else {
	    emin = 0.;
	}
	qmin = z__[(n0 << 2) - 3];
	qmax = qmin;
	for (i4 = n0 << 2; i4 >= 8; i4 += -4) {
	    if (z__[i4 - 5] <= 0.) {
		goto L100;
	    }
	    if (qmin >= emax * 4.) {
		d__1 = qmin, d__2 = z__[i4 - 3];
		qmin = (d__1<d__2) ? d__1 : d__2;
		d__1 = emax, d__2 = z__[i4 - 5];
		emax = (d__1>d__2) ? d__1 : d__2;
	    }
	    d__1 = qmax, d__2 = z__[i4 - 7] + z__[i4 - 5];
	    qmax = (d__1>d__2) ? d__1 : d__2;
	    d__1 = emin, d__2 = z__[i4 - 5];
	    emin = (d__1<d__2) ? d__1 : d__2;
	}
	i4 = 4;

L100:
	i0 = i4 / 4;
	pp = 0;

	if (n0 - i0 > 1) {
	    dee = z__[(i0 << 2) - 3];
	    deemin = dee;
	    kmin = i0;
	    i__2 = (n0 << 2) - 3;
	    for (i4 = (i0 << 2) - 3; i4 <= i__2; i4 += 4) {
		dee = z__[i4] * (dee / (dee + z__[i4 - 2]));
		if (dee <= deemin) {
		    deemin = dee;
		    kmin = (i4 + 3) / 4;
		}
	    }
	    if (2*(kmin - i0) < n0 - kmin && deemin <= z__[(n0 << 2) - 3] * 
		    .5) {
		ipn4 = 4*(i0 + n0);
		pp = 2;
		i__2 = 2*(i0 + n0 - 1);
		for (i4 = i0 << 2; i4 <= i__2; i4 += 4) {
		    temp = z__[i4 - 3];
		    z__[i4 - 3] = z__[ipn4 - i4 - 3];
		    z__[ipn4 - i4 - 3] = temp;
		    temp = z__[i4 - 2];
		    z__[i4 - 2] = z__[ipn4 - i4 - 2];
		    z__[ipn4 - i4 - 2] = temp;
		    temp = z__[i4 - 1];
		    z__[i4 - 1] = z__[ipn4 - i4 - 5];
		    z__[ipn4 - i4 - 5] = temp;
		    temp = z__[i4];
		    z__[i4] = z__[ipn4 - i4 - 4];
		    z__[ipn4 - i4 - 4] = temp;
		}
	    }
	}


	d__1 = 0., d__2 = qmin - sqrt(qmin) * 2. * sqrt(emax);
	dmin__ = -((d__1>d__2) ? d__1 : d__2);

	nbig = (n0 - i0 + 1) * 30;
	i__2 = nbig;
	for (iwhilb = 1; iwhilb <= i__2; ++iwhilb) {
	    if (i0 > n0) {
		goto L150;
	    }

	    F77_FUNC(dlasq3,DLASQ3)(&i0, &n0, &z__[1], &pp, &dmin__, &sigma, &desig, &qmax, &
		    nfail, &iter, &ndiv, &ieee);

	    pp = 1 - pp;

	    if (pp == 0 && n0 - i0 >= 3) {
		if (z__[n0 * 4] <= tol2 * qmax || z__[(n0 << 2) - 1] <= tol2 *
			 sigma) {
		    splt = i0 - 1;
		    qmax = z__[(i0 << 2) - 3];
		    emin = z__[(i0 << 2) - 1];
		    oldemn = z__[i0 * 4];
		    i__3 = 4*(n0 - 3);
		    for (i4 = i0 << 2; i4 <= i__3; i4 += 4) {
			if (z__[i4] <= tol2 * z__[i4 - 3] || z__[i4 - 1] <= 
				tol2 * sigma) {
			    z__[i4 - 1] = -sigma;
			    splt = i4 / 4;
			    qmax = 0.;
			    emin = z__[i4 + 3];
			    oldemn = z__[i4 + 4];
			} else {
			    d__1 = qmax, d__2 = z__[i4 + 1];
			    qmax = (d__1>d__2) ? d__1 : d__2;
			    d__1 = emin, d__2 = z__[i4 - 1];
			    emin = (d__1<d__2) ? d__1 : d__2;
			    d__1 = oldemn, d__2 = z__[i4];
			    oldemn = (d__1<d__2) ? d__1 : d__2;
			}
		    }
		    z__[(n0 << 2) - 1] = emin;
		    z__[n0 * 4] = oldemn;
		    i0 = splt + 1;
		}
	    }
	}

	*info = 2;
	return;

L150:
	;
    }

    *info = 3;
    return;


L170:

    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
	z__[k] = z__[(k << 2) - 3];
    }

    F77_FUNC(dlasrt,DLASRT)("D", n, &z__[1], &iinfo);

    e = 0.;
    for (k = *n; k >= 1; --k) {
	e += z__[k];
    }


    z__[(*n << 1) + 1] = trace;
    z__[(*n << 1) + 2] = e;
    z__[(*n << 1) + 3] = (double) iter;
    i__1 = *n;
    z__[(*n << 1) + 4] = (double) ndiv / (double) (i__1 * i__1);
    z__[(*n << 1) + 5] = nfail * 100. / (double) iter;

    return;

}



