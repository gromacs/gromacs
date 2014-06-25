#include <math.h>
#include "../gmx_lapack.h"
#include "lapack_limits.h"

#include "gromacs/utility/real.h"

void
F77_FUNC(ssterf,SSTERF)(int *n, 
	float *d__, 
	float *e, 
	int *info)
{
    int i__1;
    float d__1;

    float c__;
    int i__, l, m;
    float p, r__, s;
    int l1;
    float bb, rt1, rt2, eps, rte;
    int lsv;
    float eps2, oldc;
    int lend, jtot;
    float gamma, alpha, sigma, anorm;
      int iscale;
    float oldgam;
    float safmax;
    int lendsv;
    float ssfmin;
    int nmaxit;
    float ssfmax;
    int c__0 = 0;
    int c__1 = 1;
    float c_b32 = 1.;
    const float safmin = GMX_FLOAT_MIN*(1.0+GMX_FLOAT_EPS);

    --e;
    --d__;

    *info = 0;

    if (*n < 0) {
	*info = -1;
	i__1 = -(*info);
	return;
    }
    if (*n <= 1) {
	return;
    }

    eps = GMX_FLOAT_EPS;
    d__1 = eps;
    eps2 = d__1 * d__1;
    safmax = 1. / safmin;
    ssfmax = sqrt(safmax) / 3.;
    ssfmin = sqrt(safmin) / eps2;

    nmaxit = *n * 30;
    sigma = 0.;
    jtot = 0;

    l1 = 1;

L10:
    if (l1 > *n) {
      F77_FUNC(slasrt,SLASRT)("I", n, &d__[1], info);
      return;
    }
    if (l1 > 1) {
	e[l1 - 1] = 0.;
    }
    i__1 = *n - 1;
    for (m = l1; m <= i__1; ++m) {
	if (fabs(e[m]) <= sqrt(fabs(d__[m])) * 
		sqrt(fabs(d__[m + 1])) * eps) {
	    e[m] = 0.;
	    goto L30;
	}
    }
    m = *n;

L30:
    l = l1;
    lsv = l;
    lend = m;
    lendsv = lend;
    l1 = m + 1;
    if (lend == l) {
	goto L10;
    }

    i__1 = lend - l + 1;
    anorm = F77_FUNC(slanst,SLANST)("I", &i__1, &d__[l], &e[l]);
    iscale = 0;
    if (anorm > ssfmax) {
	iscale = 1;
	i__1 = lend - l + 1;
	F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &d__[l], n, 
		info);
	i__1 = lend - l;
	F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &e[l], n, 
		info);
    } else if (anorm < ssfmin) {
	iscale = 2;
	i__1 = lend - l + 1;
	F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &d__[l], n, 
		info);
	i__1 = lend - l;
	F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &e[l], n, 
		info);
    }

    i__1 = lend - 1;
    for (i__ = l; i__ <= i__1; ++i__) {
	d__1 = e[i__];
	e[i__] = d__1 * d__1;
    }

    if (fabs(d__[lend]) < fabs(d__[l])) {
	lend = lsv;
	l = lendsv;
    }

    if (lend >= l) {

L50:
	if (l != lend) {
	    i__1 = lend - 1;
	    for (m = l; m <= i__1; ++m) {
		if (fabs(e[m]) <= eps2 * fabs(d__[m] * d__[m + 1])) {
		    goto L70;
		}
	    }
	}
	m = lend;

L70:
	if (m < lend) {
	    e[m] = 0.;
	}
	p = d__[l];
	if (m == l) {
	    goto L90;
	}
	if (m == l + 1) {
	    rte = sqrt(e[l]);
	    F77_FUNC(slae2,SLAE2)(&d__[l], &rte, &d__[l + 1], &rt1, &rt2);
	    d__[l] = rt1;
	    d__[l + 1] = rt2;
	    e[l] = 0.;
	    l += 2;
	    if (l <= lend) {
		goto L50;
	    }
	    goto L150;
	}

	if (jtot == nmaxit) {
	    goto L150;
	}
	++jtot;

	rte = sqrt(e[l]);
	sigma = (d__[l + 1] - p) / (rte * 2.);
	r__ = F77_FUNC(slapy2,SLAPY2)(&sigma, &c_b32);
	sigma = p - rte / (sigma + ( (sigma>0) ? r__ : -r__));

	c__ = 1.;
	s = 0.;
	gamma = d__[m] - sigma;
	p = gamma * gamma;

	i__1 = l;
	for (i__ = m - 1; i__ >= i__1; --i__) {
	    bb = e[i__];
	    r__ = p + bb;
	    if (i__ != m - 1) {
		e[i__ + 1] = s * r__;
	    }
	    oldc = c__;
	    c__ = p / r__;
	    s = bb / r__;
	    oldgam = gamma;
	    alpha = d__[i__];
	    gamma = c__ * (alpha - sigma) - s * oldgam;
	    d__[i__ + 1] = oldgam + (alpha - gamma);
	    if (fabs(c__)>GMX_FLOAT_MIN) {
		p = gamma * gamma / c__;
	    } else {
		p = oldc * bb;
	    }
	}

	e[l] = s * p;
	d__[l] = sigma + gamma;
	goto L50;

L90:
	d__[l] = p;

	++l;
	if (l <= lend) {
	    goto L50;
	}
	goto L150;

    } else {

L100:
	i__1 = lend + 1;
	for (m = l; m >= i__1; --m) {
	    if (fabs(e[m - 1]) <= eps2 * fabs(d__[m] * d__[m - 1])) {
		goto L120;
	    }
	}
	m = lend;

L120:
	if (m > lend) {
	    e[m - 1] = 0.;
	}
	p = d__[l];
	if (m == l) {
	    goto L140;
	}

	if (m == l - 1) {
	    rte = sqrt(e[l - 1]);
	    F77_FUNC(slae2,SLAE2)(&d__[l], &rte, &d__[l - 1], &rt1, &rt2);
	    d__[l] = rt1;
	    d__[l - 1] = rt2;
	    e[l - 1] = 0.;
	    l += -2;
	    if (l >= lend) {
		goto L100;
	    }
	    goto L150;
	}

	if (jtot == nmaxit) {
	    goto L150;
	}
	++jtot;

	rte = sqrt(e[l - 1]);
	sigma = (d__[l - 1] - p) / (rte * 2.);
	r__ = F77_FUNC(slapy2,SLAPY2)(&sigma, &c_b32);
	sigma = p - rte / (sigma + ( (sigma>0) ? r__ : -r__));

	c__ = 1.;
	s = 0.;
	gamma = d__[m] - sigma;
	p = gamma * gamma;

	i__1 = l - 1;
	for (i__ = m; i__ <= i__1; ++i__) {
	    bb = e[i__];
	    r__ = p + bb;
	    if (i__ != m) {
		e[i__ - 1] = s * r__;
	    }
	    oldc = c__;
	    c__ = p / r__;
	    s = bb / r__;
	    oldgam = gamma;
	    alpha = d__[i__ + 1];
	    gamma = c__ * (alpha - sigma) - s * oldgam;
	    d__[i__] = oldgam + (alpha - gamma);
	    if (fabs(c__)>GMX_FLOAT_MIN) {
		p = gamma * gamma / c__;
	    } else {
		p = oldc * bb;
	    }
	}

	e[l - 1] = s * p;
	d__[l] = sigma + gamma;
	goto L100;

L140:
	d__[l] = p;

	--l;
	if (l >= lend) {
	    goto L100;
	}
	goto L150;

    }

L150:
    if (iscale == 1) {
	i__1 = lendsv - lsv + 1;
	F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info);
    }
    if (iscale == 2) {
	i__1 = lendsv - lsv + 1;
	F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info);
    }

    if (jtot < nmaxit) {
	goto L10;
    }
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (fabs(e[i__])>GMX_FLOAT_MIN) {
	    ++(*info);
	}
    }
    return;
}


