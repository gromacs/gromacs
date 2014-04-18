#include <math.h>

#include "gromacs/utility/real.h"

#include "../gmx_lapack.h"
#include "lapack_limits.h"

void
F77_FUNC(slarrbx,SLARRBX)(int *n, 
	 float *d__, 
	 float *l, 
	 float *ld, 
	 float *lld, 
	 int *ifirst, 
	 int *ilast, 
	 float *rtol1, 
	 float *rtol2, 
	 int *offset, 
	 float *w, 
	 float *wgap, 
	 float *werr, 
	 float *work,
	 int *iwork, 
	 int *info)
{
    int i__1, i__2, i__3;
    float d__1, d__2;

    int i__, j, k, p;
    float s;
    int i1, i2, ii, kk;
    float fac, gap, mid;
    int cnt;
    float eps, tmp, left;
    int nint, prev, next, nleft;
    float right, width, dplus, error;
    int nright, olnint;
    k = 0;
    right = 0.0;

    --iwork;
    --work;
    --werr;
    --wgap;
    --w;
    --lld;
    --ld;
    --l;
    --d__;

    *info = 0;
    eps = GMX_FLOAT_EPS;
    i__1 = *n << 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iwork[i__] = 0;
    }
    i1 = *ifirst;
    i2 = *ifirst;
    prev = 0;
    i__1 = *ilast;
    for (i__ = *ifirst; i__ <= i__1; ++i__) {
	ii = i__ - *offset;
	if (i__ == *ifirst) {
	    gap = wgap[ii];
	} else if (i__ == *ilast) {
	    gap = wgap[ii - 1];
	} else {
	    d__1 = wgap[ii - 1], d__2 = wgap[ii];
	    gap = (d__1<d__2) ? d__1 : d__2;
	}
	error = werr[ii];
	k = i__ << 1;
	iwork[k - 1] = 1;
	i2 = i__;
    }

    i__ = i1;
    nint = 0;
L30:
    if (i__ <= i2) {
	ii = i__ - *offset;
	if (iwork[(i__ << 1) - 1] == 1) {
	    fac = 1.;
	    left = w[ii] - werr[ii];


L40:
	    if (i__ > i1 && left <= right) {
		left = right;
		cnt = i__ - 1;
	    } else {
		s = -left;
		cnt = 0;
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
		    dplus = d__[j] + s;
		    s = s * lld[j] / dplus - left;
		    if (dplus < 0.) {
			++cnt;
		    }
		}
		dplus = d__[*n] + s;
		if (dplus < 0.) {
		    ++cnt;
		}
		if (! (s > 0. || s < 1.)) {

		    cnt = 0;
		    s = -left;
		    i__1 = *n - 1;
		    for (j = 1; j <= i__1; ++j) {
			dplus = d__[j] + s;
			if (dplus < 0.) {
			    ++cnt;
			}
			tmp = lld[j] / dplus;
			if (fabs(tmp)<GMX_FLOAT_MIN) {
			    s = lld[j] - left;
			} else {
			    s = s * tmp - left;
			}
		    }
		    dplus = d__[*n] + s;
		    if (dplus < 0.) {
			++cnt;
		    }
		}
		if (cnt > i__ - 1) {
		    left -= werr[ii] * fac;
		    fac *= 2.;
		    goto L40;
		}
	    }
	    nleft = cnt + 1;
	    i1 = (i1<nleft) ? i1 : nleft;
	    fac = 1.;
	    right = w[ii] + werr[ii];
L60:
	    s = -right;
	    cnt = 0;
	    i__1 = *n - 1;
	    for (j = 1; j <= i__1; ++j) {
		dplus = d__[j] + s;
		s = s * lld[j] / dplus - right;
		if (dplus < 0.) {
		    ++cnt;
		}
	    }
	    dplus = d__[*n] + s;
	    if (dplus < 0.) {
		++cnt;
	    }
	    if (! (s > 0. || s < 1.)) {

		cnt = 0;
		s = -right;
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
		    dplus = d__[j] + s;
		    if (dplus < 0.) {
			++cnt;
		    }
		    tmp = lld[j] / dplus;
		    if (fabs(tmp)<GMX_FLOAT_MIN) {
			s = lld[j] - right;
		    } else {
			s = s * tmp - right;
		    }
		}
		dplus = d__[*n] + s;
		if (dplus < 0.) {
		    ++cnt;
		}
	    }
	    if (cnt < i__) {
		right += werr[ii] * fac;
		fac *= 2.;
		goto L60;
	    }
	    cnt = (cnt<i2) ? cnt : i2;
	    ++nint;
	    k = nleft << 1;
	    work[k - 1] = left;
	    work[k] = right;
	    i__ = cnt + 1;
	    iwork[k - 1] = i__;
	    iwork[k] = cnt;
	    if (prev != nleft - 1) {
		work[k - 2] = left;
	    }
	    prev = nleft;
	} else {
	    right = work[i__ * 2];

	    ++iwork[k - 1];
	    prev = i__;
	    ++i__;
	}
	goto L30;
    }
    if (i__ <= *n && iwork[(i__ << 1) - 1] != -1) {
	work[(i__ << 1) - 1] = work[prev * 2];
    }

L80:
    prev = i1 - 1;
    olnint = nint;
    i__ = i1;
    i__1 = olnint;
    for (p = 1; p <= i__1; ++p) {
	k = i__ << 1;
	left = work[k - 1];
	right = work[k];
	next = iwork[k - 1];
	nright = iwork[k];
	mid = (left + right) * .5;
	width = right - mid;
	d__1 = fabs(left);
	d__2 = fabs(right);
	tmp = (d__1>d__2) ? d__1 : d__2;

	gap = 0.;
	if (i__ == nright) {
	    if (prev > 0 && next <= *n) {
		d__1 = left - work[k - 2], d__2 = work[k + 1] - right;
		gap = (d__1<d__2) ? d__1 : d__2;
	    } else if (prev > 0) {
		gap = left - work[k - 2];
	    } else if (next <= *n) {
		gap = work[k + 1] - right;
	    }
	}
	d__1 = *rtol1 * gap, d__2 = *rtol2 * tmp;
	if (width < ((d__1>d__2) ? d__1 : d__2)) {
	    --nint;
	    iwork[k - 1] = 0;
	    kk = k;
	    i__2 = nright;
	    for (j = i__ + 1; j <= i__2; ++j) {
		kk += 2;
		iwork[kk - 1] = 0;
		work[kk - 1] = left;
		work[kk] = right;
		wgap[j - 1 - *offset] = 0.;
	    }
	    if (i1 == i__) {
		i1 = next;
	    } else {
		iwork[(prev << 1) - 1] = next;
	    }
	    i__ = next;
	    continue;
	}
	prev = i__;

	s = -mid;
	cnt = 0;
	i__2 = *n - 1;
	for (j = 1; j <= i__2; ++j) {
	    dplus = d__[j] + s;
	    s = s * lld[j] / dplus - mid;
	    if (dplus < 0.) {
		++cnt;
	    }
	}
	dplus = d__[*n] + s;
	if (dplus < 0.) {
	    ++cnt;
	}
	if (! (s > 0. || s < 1.)) {
	    cnt = 0;
	    s = -mid;
	    i__2 = *n - 1;
	    for (j = 1; j <= i__2; ++j) {
		dplus = d__[j] + s;
		if (dplus < 0.) {
		    ++cnt;
		}
		tmp = lld[j] / dplus;
		if (fabs(tmp)<GMX_FLOAT_MIN) {
		    s = lld[j] - mid;
		} else {
		    s = s * tmp - mid;
		}
	    }
	    dplus = d__[*n] + s;
	    if (dplus < 0.) {
		++cnt;
	    }
	}
	i__2 = i__ - 1, i__3 = (nright<cnt) ? nright : cnt;
	cnt = (i__2>i__3) ? i__2 : i__3;
	if (cnt == i__ - 1) {
	    work[k - 1] = mid;
	} else if (cnt == nright) {
	    work[k] = mid;
	} else {
	    iwork[k] = cnt;
	    ++cnt;
	    iwork[k - 1] = cnt;
	    kk = cnt << 1;
	    iwork[kk - 1] = next;
	    iwork[kk] = nright;
	    work[k] = mid;
	    work[kk - 1] = mid;
	    work[kk] = right;
	    prev = cnt;
	    if (cnt - 1 > i__) {
		work[kk - 2] = mid;
	    }
	    if (cnt > *ifirst && cnt <= *ilast) {
		++nint;
	    } else if (cnt <= *ifirst) {
		i1 = cnt;
	    }
	}
	i__ = next;
    }
    if (nint > 0) {
	goto L80;
    }
    i__1 = *ilast;
    for (i__ = *ifirst; i__ <= i__1; ++i__) {
	k = i__ << 1;
	ii = i__ - *offset;
	if (iwork[k - 1] != -1) {
	    w[ii] = (work[k - 1] + work[k]) * .5;
	    werr[ii] = work[k] - w[ii];
	    if (i__ != *ilast) {
		wgap[ii] = work[k + 1] - work[k];
	    }
	}
    }

    return;

} 
