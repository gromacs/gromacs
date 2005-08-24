#include <math.h>
#include <types/simple.h>

#include "gmx_blas.h"
#include "gmx_lapack.h"
#include "lapack_limits.h"

void
F77_FUNC(dsteqr,DSTEQR)(char *    compz, 
	int *     n, 
	double *  d__, 
	double *  e, 
	double *  z__, 
	int *     ldz, 
	double *  work, 
	int *     info)
{
    double c_b9 = 0.;
    double c_b10 = 1.;
    int c__0 = 0;
    int c__1 = 1;
    int c__2 = 2;
    int z_dim1, z_offset, i__1, i__2;
    double d__1, d__2;

    double b, c__, f, g;
    int i__, j, k, l, m;
    double p, r__, s;
    int l1, ii, mm, lm1, mm1, nm1;
    double rt1, rt2, eps;
    int lsv;
    double tst, eps2;
    int lend, jtot;
    double anorm;
    int lendm1, lendp1;
    int iscale;
    double safmin,minval;
    double safmax;
    int lendsv;
    double ssfmin;
    int nmaxit, icompz;
    double ssfmax;


    --d__;
    --e;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;

    *info = 0;

    if (*compz=='N' || *compz=='n') {
	icompz = 0;
    } else if (*compz=='V' || *compz=='v') {
	icompz = 1;
    } else if (*compz=='I' || *compz=='i') {
	icompz = 2;
    } else {
	icompz = -1;
    }
    if (icompz < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*ldz < 1 || (icompz > 0 && *ldz < ((*n>1) ? *n : 1))) {
	*info = -6;
    }
    if (*info != 0) {
	return;
    }


    if (*n == 0) {
	return;
    }

    if (*n == 1) {
	if (icompz == 2) {
	    z__[z_dim1 + 1] = 1.;
	}
	return;
    }

    eps = GMX_DOUBLE_EPS;
    d__1 = eps;
    eps2 = d__1 * d__1;
    minval = GMX_DOUBLE_MIN;
    safmin = minval*(1.0+GMX_DOUBLE_EPS);

    safmax = 1. / safmin;
    ssfmax = sqrt(safmax) / 3.;
    ssfmin = sqrt(safmin) / eps2;

    if (icompz == 2) {
	F77_FUNC(dlaset,DLASET)("Full", n, n, &c_b9, &c_b10, &z__[z_offset], ldz);
    }

    nmaxit = *n * 30;
    jtot = 0;

    l1 = 1;
    nm1 = *n - 1;

L10:
    if (l1 > *n) {
	goto L160;
    }
    if (l1 > 1) {
	e[l1 - 1] = 0.;
    }
    if (l1 <= nm1) {
	i__1 = nm1;
	for (m = l1; m <= i__1; ++m) {
	    tst = fabs(e[m]);
	    if (fabs(tst)<GMX_DOUBLE_MIN) {
		goto L30;
	    }
	    if (tst <= sqrt(fabs(d__[m])) * sqrt(fabs(d__[m + 1])) * eps) {
		e[m] = 0.;
		goto L30;
	    }
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
    anorm = F77_FUNC(dlanst,DLANST)("I", &i__1, &d__[l], &e[l]);
    iscale = 0;
    if (fabs(anorm)<GMX_DOUBLE_MIN) {
	goto L10;
    }
    if (anorm > ssfmax) {
	iscale = 1;
	i__1 = lend - l + 1;
	F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &d__[l], n, 
		info);
	i__1 = lend - l;
	F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &e[l], n, 
		info);
    } else if (anorm < ssfmin) {
	iscale = 2;
	i__1 = lend - l + 1;
	F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &d__[l], n, 
		info);
	i__1 = lend - l;
	F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &e[l], n, 
		info);
    }

    if (fabs(d__[lend]) < fabs(d__[l])) {
	lend = lsv;
	l = lendsv;
    }

    if (lend > l) {

L40:
	if (l != lend) {
	    lendm1 = lend - 1;
	    i__1 = lendm1;
	    for (m = l; m <= i__1; ++m) {
  	        d__2 = fabs(e[m]);
		tst = d__2 * d__2;
		if (tst <= eps2 * fabs(d__[m]) * fabs(d__[m+ 1]) + safmin) {
		    goto L60;
		}
	    }
	}

	m = lend;

L60:
	if (m < lend) {
	    e[m] = 0.;
	}
	p = d__[l];
	if (m == l) {
	    goto L80;
	}

	if (m == l + 1) {
	    if (icompz > 0) {
		F77_FUNC(dlaev2,DLAEV2)(&d__[l], &e[l], &d__[l + 1], &rt1, &rt2, &c__, &s);
		work[l] = c__;
		work[*n - 1 + l] = s;
		F77_FUNC(dlasr,DLASR)("R", "V", "B", n, &c__2, &work[l], &work[*n - 1 + l], &
			z__[l * z_dim1 + 1], ldz);
	    } else {
		F77_FUNC(dlae2,DLAE2)(&d__[l], &e[l], &d__[l + 1], &rt1, &rt2);
	    }
	    d__[l] = rt1;
	    d__[l + 1] = rt2;
	    e[l] = 0.;
	    l += 2;
	    if (l <= lend) {
		goto L40;
	    }
	    goto L140;
	}

	if (jtot == nmaxit) {
	    goto L140;
	}
	++jtot;

	g = (d__[l + 1] - p) / (e[l] * 2.);
	r__ = F77_FUNC(dlapy2,DLAPY2)(&g, &c_b10);
	g = d__[m] - p + e[l] / (g + ( (g>0) ? r__ : -r__ ) );

	s = 1.;
	c__ = 1.;
	p = 0.;

	mm1 = m - 1;
	i__1 = l;
	for (i__ = mm1; i__ >= i__1; --i__) {
	    f = s * e[i__];
	    b = c__ * e[i__];
	    F77_FUNC(dlartg,DLARTG)(&g, &f, &c__, &s, &r__);
	    if (i__ != m - 1) {
		e[i__ + 1] = r__;
	    }
	    g = d__[i__ + 1] - p;
	    r__ = (d__[i__] - g) * s + c__ * 2. * b;
	    p = s * r__;
	    d__[i__ + 1] = g + p;
	    g = c__ * r__ - b;

	    if (icompz > 0) {
		work[i__] = c__;
		work[*n - 1 + i__] = -s;
	    }
	}

	if (icompz > 0) {
	    mm = m - l + 1;
	    F77_FUNC(dlasr,DLASR)("R", "V", "B", n, &mm, &work[l], &work[*n - 1 + l], &z__[l 
		    * z_dim1 + 1], ldz);
	}

	d__[l] -= p;
	e[l] = g;
	goto L40;

L80:
	d__[l] = p;

	++l;
	if (l <= lend) {
	    goto L40;
	}
	goto L140;

    } else {

L90:
	if (l != lend) {
	    lendp1 = lend + 1;
	    i__1 = lendp1;
	    for (m = l; m >= i__1; --m) {
		d__2 = fabs(e[m - 1]);
		tst = d__2 * d__2;
		if (tst <= eps2 * fabs(d__[m]) * fabs(d__[m- 1]) + safmin) {
		    goto L110;
		}
	    }
	}

	m = lend;

L110:
	if (m > lend) {
	    e[m - 1] = 0.;
	}
	p = d__[l];
	if (m == l) {
	    goto L130;
	}
	if (m == l - 1) {
	    if (icompz > 0) {
		F77_FUNC(dlaev2,DLAEV2)(&d__[l - 1], &e[l - 1], &d__[l], &rt1, &rt2, &c__, &s)
			;
		work[m] = c__;
		work[*n - 1 + m] = s;
		F77_FUNC(dlasr,DLASR)("R", "V", "F", n, &c__2, &work[m], &work[*n - 1 + m], &
			z__[(l - 1) * z_dim1 + 1], ldz);
	    } else {
		F77_FUNC(dlae2,DLAE2)(&d__[l - 1], &e[l - 1], &d__[l], &rt1, &rt2);
	    }
	    d__[l - 1] = rt1;
	    d__[l] = rt2;
	    e[l - 1] = 0.;
	    l += -2;
	    if (l >= lend) {
		goto L90;
	    }
	    goto L140;
	}

	if (jtot == nmaxit) {
	    goto L140;
	}
	++jtot;

	g = (d__[l - 1] - p) / (e[l - 1] * 2.);
	r__ = F77_FUNC(dlapy2,DLAPY2)(&g, &c_b10);
	g = d__[m] - p + e[l - 1] / (g + ( (g>0) ? r__ : -r__ ));

	s = 1.;
	c__ = 1.;
	p = 0.;

	lm1 = l - 1;
	i__1 = lm1;
	for (i__ = m; i__ <= i__1; ++i__) {
	    f = s * e[i__];
	    b = c__ * e[i__];
	    F77_FUNC(dlartg,DLARTG)(&g, &f, &c__, &s, &r__);
	    if (i__ != m) {
		e[i__ - 1] = r__;
	    }
	    g = d__[i__] - p;
	    r__ = (d__[i__ + 1] - g) * s + c__ * 2. * b;
	    p = s * r__;
	    d__[i__] = g + p;
	    g = c__ * r__ - b;

	    if (icompz > 0) {
		work[i__] = c__;
		work[*n - 1 + i__] = s;
	    }
	}

	if (icompz > 0) {
	    mm = l - m + 1;
	    F77_FUNC(dlasr,DLASR)("R", "V", "F", n, &mm, &work[m], &work[*n - 1 + m], &z__[m 
		    * z_dim1 + 1], ldz);
	}

	d__[l] -= p;
	e[lm1] = g;
	goto L90;

L130:
	d__[l] = p;

	--l;
	if (l >= lend) {
	    goto L90;
	}
	goto L140;

    }

L140:
    if (iscale == 1) {
	i__1 = lendsv - lsv + 1;
	F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info);
	i__1 = lendsv - lsv;
	F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &e[lsv], n, 
		info);
    } else if (iscale == 2) {
	i__1 = lendsv - lsv + 1;
	F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &d__[lsv], 
		n, info);
	i__1 = lendsv - lsv;
	F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &e[lsv], n, 
		info);
    }

    if (jtot < nmaxit) {
	goto L10;
    }
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (fabs(e[i__])>GMX_DOUBLE_MIN) {
	    ++(*info);
	}
    }
    goto L190;

L160:
    if (icompz == 0) {

	F77_FUNC(dlasrt,DLASRT)("I", n, &d__[1], info);

    } else {

	i__1 = *n;
	for (ii = 2; ii <= i__1; ++ii) {
	    i__ = ii - 1;
	    k = i__;
	    p = d__[i__];
	    i__2 = *n;
	    for (j = ii; j <= i__2; ++j) {
		if (d__[j] < p) {
		    k = j;
		    p = d__[j];
		}
	    }
	    if (k != i__) {
		d__[k] = d__[i__];
		d__[i__] = p;
		F77_FUNC(dswap,DSWAP)(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[k * z_dim1 + 1],
			 &c__1);
	    }
	}
    }

L190:
    return;
}


