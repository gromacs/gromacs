#include <ctype.h>
#include <math.h>

#include "gmx_blas.h"
#include "gmx_lapack.h"

#include <types/simple.h>

void 
F77_FUNC(dbdsqr,DBDSQR)(char *uplo,
                        int *n,
                        int *ncvt,
                        int *nru, 
                        int *ncc, 
                        double *d__,
                        double *e,
                        double *vt, 
                        int *ldvt,
                        double *u, 
                        int *ldu,
                        double *c__, 
                        int *ldc,
                        double *work,
                        int *info)
{
    char xuplo = toupper(*uplo);
    int c_dim1, c_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1, 
	    i__2;
    double r__1, r__2, r__3, r__4;
    double c_b15 = -.125;

    int c__1 = 1;
    double c_b49 = 1.f;
    double c_b72 = -1.f;

    double f, g, h__;
    int i__, j, m;
    double r__, cs;
    int ll;
    double sn, mu;
    int nm1, nm12, nm13, lll;
    double eps, sll, tol, abse;
    int idir;
    double abss;
    int oldm;
    double cosl;
    int isub, iter;
    double unfl, sinl, cosr, smin, smax, sinr;
    double oldcs;
    int oldll;
    double shift, sigmn, oldsn;
    int maxit;
    double sminl;
    double sigmx;
    int lower;
    double sminoa;
    double thresh;
    int rotate;
    double sminlo, tolmul;
    int itmp1,itmp2;
    double ftmp;
    
    --d__;
    --e;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;

    *info = 0;
    
    itmp1 = (*n > 1) ? *n : 1;
    itmp2 = (*nru > 1) ? *nru : 1;
    
    lower = (xuplo == 'L');
    if ( (xuplo!='U') && !lower) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*ncvt < 0) {
	*info = -3;
    } else if (*nru < 0) {
	*info = -4;
    } else if (*ncc < 0) {
	*info = -5;
    } else if ( ((*ncvt == 0) && (*ldvt < 1)) || ((*ncvt > 0) && (*ldvt < itmp1)) ) {
	*info = -9;
    } else if (*ldu < itmp2) {
	*info = -11;
    } else if ( ((*ncc == 0) && (*ldc < 1)) || ((*ncc > 0) && (*ldc < itmp1))) {
	*info = -13;
    }
    if (*info != 0) {
	return;
    }
    if (*n == 0) {
	return;
    }
    if (*n == 1) {
	goto L160;
    }

    rotate = *ncvt > 0 || *nru > 0 || *ncc > 0;

    if (! rotate) {
	F77_FUNC(dlasq1,DLASQ1)(n, &d__[1], &e[1], &work[1], info);
	return;
    }

    nm1 = *n - 1;
    nm12 = nm1 + nm1;
    nm13 = nm12 + nm1;
    idir = 0;

    eps = GMX_DOUBLE_EPS;
    unfl = GMX_DOUBLE_MIN/GMX_DOUBLE_EPS;

    if (lower) {
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    F77_FUNC(dlartg,DLARTG)(&d__[i__], &e[i__], &cs, &sn, &r__);
	    d__[i__] = r__;
	    e[i__] = sn * d__[i__ + 1];
	    d__[i__ + 1] = cs * d__[i__ + 1];
	    work[i__] = cs;
	    work[nm1 + i__] = sn;
	}

	if (*nru > 0) {
	    F77_FUNC(dlasr,DLASR)("R", "V", "F", nru, n, &work[1], &work[*n], &u[u_offset], 
		    ldu);
	}
	if (*ncc > 0) {
	    F77_FUNC(dlasr,DLASR)("L", "V", "F", n, ncc, &work[1], &work[*n], &c__[c_offset],
		     ldc);
	}
    }

    r__3 = 100.f, r__4 = pow(GMX_DOUBLE_EPS,c_b15);
    r__1 = 10.f, r__2 = (r__3<r__4) ? r__3 : r__4;
    tolmul = (r__1>r__2) ? r__1 : r__2;
    tol = tolmul * eps;
    smax = 0.f;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r__2 = smax, r__3 = (r__1 = d__[i__], fabs(r__1));
	smax = (r__2>r__3) ? r__2 : r__3;
    }
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r__2 = smax, r__3 = (r__1 = e[i__], fabs(r__1));
	smax = (r__2>r__3) ? r__2 : r__3;
    }
    sminl = 0.f;
    if (tol >= 0.f) {
	sminoa = fabs(d__[1]);
	if (sminoa == 0.f) {
	    goto L50;
	}
	mu = sminoa;
	i__1 = *n;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    mu = (r__2 = d__[i__], fabs(r__2)) * (mu / (mu + (r__1 = e[i__ - 
		    1], fabs(r__1))));
	    sminoa = (sminoa<mu) ? sminoa : mu;
	    if (sminoa == 0.f) {
		goto L50;
	    }
	}
L50:
	sminoa /= sqrt((double) (*n));
	r__1 = tol * sminoa, r__2 = *n * 6 * *n * unfl;
	thresh = (r__1>r__2) ? r__1 : r__2;
    } else {
	r__1 = fabs(tol) * smax, r__2 = *n * 6 * *n * unfl;
	thresh = (r__1>r__2) ? r__1 : r__2;
    }
    maxit = *n * 6 * *n;
    iter = 0;
    oldll = -1;
    oldm = -1;
    m = *n;

L60:

    if (m <= 1) {
	goto L160;
    }
    if (iter > maxit) {
	goto L200;
    }

    if (tol < 0.f && (r__1 = d__[m], fabs(r__1)) <= thresh) {
	d__[m] = 0.f;
    }
    smax = (r__1 = d__[m], fabs(r__1));
    smin = smax;
    i__1 = m - 1;
    for (lll = 1; lll <= i__1; ++lll) {
	ll = m - lll;
	abss = (r__1 = d__[ll], fabs(r__1));
	abse = (r__1 = e[ll], fabs(r__1));
	if (tol < 0.f && abss <= thresh) {
	    d__[ll] = 0.f;
	}
	if (abse <= thresh) {
	    goto L80;
	}
	smin = (smin<abss) ? smin : abss;
	r__1 = (smax>abss) ? smax : abss;
	smax = (r__1>abse) ? r__1 : abse;
    }
    ll = 0;
    goto L90;
L80:
    e[ll] = 0.f;
    if (ll == m - 1) {
	--m;
	goto L60;
    }
L90:
    ++ll;
    if (ll == m - 1) {
	F77_FUNC(dlasv2,DLASV2)(&d__[m - 1], &e[m - 1], &d__[m], &sigmn, &sigmx, &sinr, &cosr,
		 &sinl, &cosl);
	d__[m - 1] = sigmx;
	e[m - 1] = 0.f;
	d__[m] = sigmn;
	if (*ncvt > 0) {
	    F77_FUNC(drot,DROT)(ncvt, &vt[m - 1 + vt_dim1], ldvt, &vt[m + vt_dim1], ldvt, &
		    cosr, &sinr);
	}
	if (*nru > 0) {
	    F77_FUNC(drot,DROT)(nru, &u[(m - 1) * u_dim1 + 1], &c__1, &u[m * u_dim1 + 1], &
		    c__1, &cosl, &sinl);
	}
	if (*ncc > 0) {
	    F77_FUNC(drot,DROT)(ncc, &c__[m - 1 + c_dim1], ldc, &c__[m + c_dim1], ldc, &
		    cosl, &sinl);
	}
	m += -2;
	goto L60;
    }
    if (ll > oldm || m < oldll) {
	if ((r__1 = d__[ll], fabs(r__1)) >= (r__2 = d__[m], fabs(r__2))) {
	    idir = 1;
	} else {
	    idir = 2;
	}
    }
    if (idir == 1) {

        if( (fabs(e[m-1]) <= fabs(tol) * fabs(d__[m])) ||
            (tol<0.0 && fabs(e[m-1])<=thresh)) {
	    e[m - 1] = 0.f;
	    goto L60;
	}
	if (tol >= 0.f) {
	    mu = (r__1 = d__[ll], fabs(r__1));
	    sminl = mu;
	    i__1 = m - 1;
	    for (lll = ll; lll <= i__1; ++lll) {
		if ((r__1 = e[lll], fabs(r__1)) <= tol * mu) {
		    e[lll] = 0.f;
		    goto L60;
		}
		sminlo = sminl;
		mu = (r__2 = d__[lll + 1], fabs(r__2)) * (mu / (mu + (r__1 = 
			e[lll], fabs(r__1))));
		sminl = (sminl<mu) ? sminl : mu;
	    }
	}
    } else {
        if( (fabs(e[ll]) <= fabs(tol)*fabs(d__[ll])) ||
            (tol<0.0 && fabs(e[ll])<=thresh)) {
	    e[ll] = 0.f;
	    goto L60;
	}
	if (tol >= 0.f) {
	    mu = (r__1 = d__[m], fabs(r__1));
	    sminl = mu;
	    i__1 = ll;
	    for (lll = m - 1; lll >= i__1; --lll) {
		if ((r__1 = e[lll], fabs(r__1)) <= tol * mu) {
		    e[lll] = 0.f;
		    goto L60;
		}
		sminlo = sminl;
		mu = (r__2 = d__[lll], fabs(r__2)) * (mu / (mu + (r__1 = e[
			lll], fabs(r__1))));
		sminl = (sminl<mu) ? sminl : mu;
	    }
	}
    }
    oldll = ll;
    oldm = m;

    r__1 = eps, r__2 = tol * .01f;
    if (tol >= 0.f && *n * tol * (sminl / smax) <= ((r__1>r__2) ? r__1 : r__2)) {
	shift = 0.f;
    } else {
	if (idir == 1) {
	    sll = (r__1 = d__[ll], fabs(r__1));
	    F77_FUNC(dlas2,DLAS2)(&d__[m - 1], &e[m - 1], &d__[m], &shift, &r__);
	} else {
	    sll = (r__1 = d__[m], fabs(r__1));
	    F77_FUNC(dlas2,DLAS2)(&d__[ll], &e[ll], &d__[ll + 1], &shift, &r__);
	}
	if (sll > 0.f) {
	    r__1 = shift / sll;
	    if (r__1 * r__1 < eps) {
		shift = 0.f;
	    }
	}
    }
    iter = iter + m - ll;
    if (shift == 0.f) {
	if (idir == 1) {
	    cs = 1.f;
	    oldcs = 1.f;
	    i__1 = m - 1;
	    for (i__ = ll; i__ <= i__1; ++i__) {
		r__1 = d__[i__] * cs;
		F77_FUNC(dlartg,DLARTG)(&r__1, &e[i__], &cs, &sn, &r__);
		if (i__ > ll) {
		    e[i__ - 1] = oldsn * r__;
		}
		r__1 = oldcs * r__;
		r__2 = d__[i__ + 1] * sn;
		F77_FUNC(dlartg,DLARTG)(&r__1, &r__2, &oldcs, &oldsn, &d__[i__]);
		work[i__ - ll + 1] = cs;
		work[i__ - ll + 1 + nm1] = sn;
		work[i__ - ll + 1 + nm12] = oldcs;
		work[i__ - ll + 1 + nm13] = oldsn;
	    }
	    h__ = d__[m] * cs;
	    d__[m] = h__ * oldcs;
	    e[m - 1] = h__ * oldsn;
	    if (*ncvt > 0) {
		i__1 = m - ll + 1;
		F77_FUNC(dlasr,DLASR)("L", "V", "F", &i__1, ncvt, &work[1], &work[*n], &vt[
			ll + vt_dim1], ldvt);
	    }
	    if (*nru > 0) {
		i__1 = m - ll + 1;
		F77_FUNC(dlasr,DLASR)("R", "V", "F", nru, &i__1, &work[nm12 + 1], &work[nm13 
			+ 1], &u[ll * u_dim1 + 1], ldu);
	    }
	    if (*ncc > 0) {
		i__1 = m - ll + 1;
		F77_FUNC(dlasr,DLASR)("L", "V", "F", &i__1, ncc, &work[nm12 + 1], &work[nm13 
			+ 1], &c__[ll + c_dim1], ldc);
	    }
	    if ((r__1 = e[m - 1], fabs(r__1)) <= thresh) {
		e[m - 1] = 0.f;
	    }
	} else {
	    cs = 1.f;
	    oldcs = 1.f;
	    i__1 = ll + 1;
	    for (i__ = m; i__ >= i__1; --i__) {
		r__1 = d__[i__] * cs;
		F77_FUNC(dlartg,DLARTG)(&r__1, &e[i__ - 1], &cs, &sn, &r__);
		if (i__ < m) {
		    e[i__] = oldsn * r__;
		}
		r__1 = oldcs * r__;
		r__2 = d__[i__ - 1] * sn;
		F77_FUNC(dlartg,DLARTG)(&r__1, &r__2, &oldcs, &oldsn, &d__[i__]);
		work[i__ - ll] = cs;
		work[i__ - ll + nm1] = -sn;
		work[i__ - ll + nm12] = oldcs;
		work[i__ - ll + nm13] = -oldsn;
	    }
	    h__ = d__[ll] * cs;
	    d__[ll] = h__ * oldcs;
	    e[ll] = h__ * oldsn;
	    if (*ncvt > 0) {
		i__1 = m - ll + 1;
		F77_FUNC(dlasr,DLASR)("L", "V", "B", &i__1, ncvt, &work[nm12 + 1], &work[
			nm13 + 1], &vt[ll + vt_dim1], ldvt);
	    }
	    if (*nru > 0) {
		i__1 = m - ll + 1;
		F77_FUNC(dlasr,DLASR)("R", "V", "B", nru, &i__1, &work[1], &work[*n], &u[ll *
			 u_dim1 + 1], ldu);
	    }
	    if (*ncc > 0) {
		i__1 = m - ll + 1;
		F77_FUNC(dlasr,DLASR)("L", "V", "B", &i__1, ncc, &work[1], &work[*n], &c__[
			ll + c_dim1], ldc);
	    }
	    if ((r__1 = e[ll], fabs(r__1)) <= thresh) {
		e[ll] = 0.f;
	    }
	}
    } else {

	if (idir == 1) {
	    f = ((r__1 = d__[ll], fabs(r__1)) - shift) * ( ((d__[ll] > 0) ? c_b49 : -c_b49) + shift / d__[ll]);
	    g = e[ll];
	    i__1 = m - 1;
	    for (i__ = ll; i__ <= i__1; ++i__) {
		F77_FUNC(dlartg,DLARTG)(&f, &g, &cosr, &sinr, &r__);
		if (i__ > ll) {
		    e[i__ - 1] = r__;
		}
		f = cosr * d__[i__] + sinr * e[i__];
		e[i__] = cosr * e[i__] - sinr * d__[i__];
		g = sinr * d__[i__ + 1];
		d__[i__ + 1] = cosr * d__[i__ + 1];
		F77_FUNC(dlartg,DLARTG)(&f, &g, &cosl, &sinl, &r__);
		d__[i__] = r__;
		f = cosl * e[i__] + sinl * d__[i__ + 1];
		d__[i__ + 1] = cosl * d__[i__ + 1] - sinl * e[i__];
		if (i__ < m - 1) {
		    g = sinl * e[i__ + 1];
		    e[i__ + 1] = cosl * e[i__ + 1];
		}
		work[i__ - ll + 1] = cosr;
		work[i__ - ll + 1 + nm1] = sinr;
		work[i__ - ll + 1 + nm12] = cosl;
		work[i__ - ll + 1 + nm13] = sinl;
	    }
	    e[m - 1] = f;

	    if (*ncvt > 0) {
		i__1 = m - ll + 1;
		F77_FUNC(dlasr,DLASR)("L", "V", "F", &i__1, ncvt, &work[1], &work[*n], &vt[
			ll + vt_dim1], ldvt);
	    }
	    if (*nru > 0) {
		i__1 = m - ll + 1;
		F77_FUNC(dlasr,DLASR)("R", "V", "F", nru, &i__1, &work[nm12 + 1], &work[nm13 
			+ 1], &u[ll * u_dim1 + 1], ldu);
	    }
	    if (*ncc > 0) {
		i__1 = m - ll + 1;
		F77_FUNC(dlasr,DLASR)("L", "V", "F", &i__1, ncc, &work[nm12 + 1], &work[nm13 
			+ 1], &c__[ll + c_dim1], ldc);
	    }
	    if ((r__1 = e[m - 1], fabs(r__1)) <= thresh) {
		e[m - 1] = 0.f;
	    }
	} else {

	    f = ((r__1 = d__[m], fabs(r__1)) - shift) * ( ((d__[m] > 0) ? c_b49 : -c_b49) + shift / d__[m]);
	    g = e[m - 1];
	    i__1 = ll + 1;
	    for (i__ = m; i__ >= i__1; --i__) {
		F77_FUNC(dlartg,DLARTG)(&f, &g, &cosr, &sinr, &r__);
		if (i__ < m) {
		    e[i__] = r__;
		}
		f = cosr * d__[i__] + sinr * e[i__ - 1];
		e[i__ - 1] = cosr * e[i__ - 1] - sinr * d__[i__];
		g = sinr * d__[i__ - 1];
		d__[i__ - 1] = cosr * d__[i__ - 1];
		F77_FUNC(dlartg,DLARTG)(&f, &g, &cosl, &sinl, &r__);
		d__[i__] = r__;
		f = cosl * e[i__ - 1] + sinl * d__[i__ - 1];
		d__[i__ - 1] = cosl * d__[i__ - 1] - sinl * e[i__ - 1];
		if (i__ > ll + 1) {
		    g = sinl * e[i__ - 2];
		    e[i__ - 2] = cosl * e[i__ - 2];
		}
		work[i__ - ll] = cosr;
		work[i__ - ll + nm1] = -sinr;
		work[i__ - ll + nm12] = cosl;
		work[i__ - ll + nm13] = -sinl;
	    }
	    e[ll] = f;

	    if ((r__1 = e[ll], fabs(r__1)) <= thresh) {
		e[ll] = 0.f;
	    }
	    if (*ncvt > 0) {
		i__1 = m - ll + 1;
		F77_FUNC(dlasr,DLASR)("L", "V", "B", &i__1, ncvt, &work[nm12 + 1], &work[
			nm13 + 1], &vt[ll + vt_dim1], ldvt);
	    }
	    if (*nru > 0) {
		i__1 = m - ll + 1;
		F77_FUNC(dlasr,DLASR)("R", "V", "B", nru, &i__1, &work[1], &work[*n], &u[ll *
			 u_dim1 + 1], ldu);
	    }
	    if (*ncc > 0) {
		i__1 = m - ll + 1;
		F77_FUNC(dlasr,DLASR)("L", "V", "B", &i__1, ncc, &work[1], &work[*n], &c__[
			ll + c_dim1], ldc);
	    }
	}
    }

    goto L60;

L160:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (d__[i__] < 0.f) {
	    d__[i__] = -d__[i__];

	    if (*ncvt > 0) {
		F77_FUNC(dscal,DSCAL)(ncvt, &c_b72, &vt[i__ + vt_dim1], ldvt);
	    }
	}
    }

    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {

	isub = 1;
	smin = d__[1];
	i__2 = *n + 1 - i__;
	for (j = 2; j <= i__2; ++j) {
	    if (d__[j] <= smin) {
		isub = j;
		smin = d__[j];
	    }
	}
	if (isub != *n + 1 - i__) {
	    d__[isub] = d__[*n + 1 - i__];
	    d__[*n + 1 - i__] = smin;
	    if (*ncvt > 0) {
		F77_FUNC(dswap,DSWAP)(ncvt, &vt[isub + vt_dim1], ldvt, &vt[*n + 1 - i__ + 
			vt_dim1], ldvt);
	    }
	    if (*nru > 0) {
		F77_FUNC(dswap,DSWAP)(nru, &u[isub * u_dim1 + 1], &c__1, &u[(*n + 1 - i__) * 
			u_dim1 + 1], &c__1);
	    }
	    if (*ncc > 0) {
		F77_FUNC(dswap,DSWAP)(ncc, &c__[isub + c_dim1], ldc, &c__[*n + 1 - i__ + 
			c_dim1], ldc);
	    }
	}
    }
    goto L220;

L200:
    *info = 0;
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (e[i__] != 0.f) {
	    ++(*info);
	}
    }
L220:
    return;

}


