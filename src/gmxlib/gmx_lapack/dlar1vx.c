#include <math.h>

#include <types/simple.h>

#include "gmx_lapack.h"
#include "lapack_limits.h"

void F77_FUNC(dlar1vx,DLAR1VX)(int *n, 
	      int *b1, 
	      int *bn,
	      double *sigma, 
	      double *d__, 
	      double *l, 
	      double *ld, 
	      double *lld, 
	      double *eval, 
	      double *gersch, 
	      double *z__, 
	      double *ztz, 
	      double *mingma, 
	      int *r__, 
	      int *isuppz, 
	      double *work)
{
    int i__1;

    int i__, j;
    double s;
    int r1, r2;
    int to;
    double eps, tmp;
    int indp, inds, from;
    double dplus;
    int sawnan;
    int indumn;
    double dminus;

    --work;
    --isuppz;
    --z__;
    --gersch;
    --lld;
    --ld;
    --l;
    --d__;

    /* Function Body */
    eps = GMX_DOUBLE_EPS;
    if (*r__ == 0) {

	r1 = *b1;
	r2 = *bn;
	i__1 = *bn;
	for (i__ = *b1; i__ <= i__1; ++i__) {
	    if (*eval >= gersch[(i__ << 1) - 1] && *eval <= gersch[i__ * 2]) {
		r1 = i__;
		goto L20;
	    }
	}
	goto L40;
L20:
	i__1 = *b1;
	for (i__ = *bn; i__ >= i__1; --i__) {
	    if (*eval >= gersch[(i__ << 1) - 1] && *eval <= gersch[i__ * 2]) {
		r2 = i__;
		goto L40;
	    }
	}
    } else {
	r1 = *r__;
	r2 = *r__;
    }

L40:
    indumn = *n;
    inds = (*n << 1) + 1;
    indp = *n * 3 + 1;
    sawnan = 0;

    if (*b1 == 1) {
	work[inds] = 0.;
    } else {
	work[inds] = lld[*b1 - 1];
    }
    s = work[inds] - *sigma;
    i__1 = r2 - 1;
    for (i__ = *b1; i__ <= i__1; ++i__) {
	dplus = d__[i__] + s;
	work[i__] = ld[i__] / dplus;
	work[inds + i__] = s * work[i__] * l[i__];
	s = work[inds + i__] - *sigma;
    }

    if (! (s > 0. || s < 1.)) {

	sawnan = 1;
	j = *b1 + 1;
L60:
	if (work[inds + j] > 0. || work[inds + j] < 1.) {
	    ++j;
	    goto L60;
	}
	work[inds + j] = lld[j];
	s = work[inds + j] - *sigma;
	i__1 = r2 - 1;
	for (i__ = j + 1; i__ <= i__1; ++i__) {
	    dplus = d__[i__] + s;
	    work[i__] = ld[i__] / dplus;
	    if (fabs(work[i__])<GMX_DOUBLE_MIN) {
		work[inds + i__] = lld[i__];
	    } else {
		work[inds + i__] = s * work[i__] * l[i__];
	    }
	    s = work[inds + i__] - *sigma;
	}
    }

    work[indp + *bn - 1] = d__[*bn] - *sigma;
    i__1 = r1;
    for (i__ = *bn - 1; i__ >= i__1; --i__) {
	dminus = lld[i__] + work[indp + i__];
	tmp = d__[i__] / dminus;
	work[indumn + i__] = l[i__] * tmp;
	work[indp + i__ - 1] = work[indp + i__] * tmp - *sigma;
    }
    tmp = work[indp + r1 - 1];
    if (! (tmp > 0. || tmp < 1.)) {

	sawnan = 1;
	j = *bn - 3;
L90:
	if (work[indp + j] > 0. || work[indp + j] < 1.) {
	    --j;
	    goto L90;
	}
	work[indp + j] = d__[j + 1] - *sigma;
	i__1 = r1;
	for (i__ = j; i__ >= i__1; --i__) {
	    dminus = lld[i__] + work[indp + i__];
	    tmp = d__[i__] / dminus;
	    work[indumn + i__] = l[i__] * tmp;
	    if (fabs(tmp)<GMX_DOUBLE_MIN) {
		work[indp + i__ - 1] = d__[i__] - *sigma;
	    } else {
		work[indp + i__ - 1] = work[indp + i__] * tmp - *sigma;
	    }
	}
    }

    *mingma = work[inds + r1 - 1] + work[indp + r1 - 1];
    if (fabs(*mingma)<GMX_DOUBLE_MIN) {
	*mingma = eps * work[inds + r1 - 1];
    }
    *r__ = r1;
    i__1 = r2 - 1;
    for (i__ = r1; i__ <= i__1; ++i__) {
	tmp = work[inds + i__] + work[indp + i__];
	if (fabs(tmp)<GMX_DOUBLE_MIN) {
	    tmp = eps * work[inds + i__];
	}
	if (fabs(tmp) < fabs(*mingma)) {
	    *mingma = tmp;
	    *r__ = i__ + 1;
	}
    }

    isuppz[1] = *b1;
    isuppz[2] = *bn;
    z__[*r__] = 1.;
    *ztz = 1.;
    if (! sawnan) {
	from = *r__ - 1;
	i__1 = *r__ - 32;
	to = (i__1>(*b1)) ? i__1 : (*b1);
L120:
	if (from >= *b1) {
	    i__1 = to;
	    for (i__ = from; i__ >= i__1; --i__) {
		z__[i__] = -(work[i__] * z__[i__ + 1]);
		*ztz += z__[i__] * z__[i__];
	    }
	    if (fabs(z__[to]) <= eps && fabs(z__[to + 1]) <= eps) {
		isuppz[1] = to + 2;
	    } else {
		from = to - 1;
		i__1 = to - 32;
		to = (i__1>*b1) ? i__1 : *b1;
		goto L120;
	    }
	}
	from = *r__ + 1;
	i__1 = *r__ + 32;
	to = (i__1<*bn) ? i__1 : *bn;
L140:
	if (from <= *bn) {
	    i__1 = to;
	    for (i__ = from; i__ <= i__1; ++i__) {
		z__[i__] = -(work[indumn + i__ - 1] * z__[i__ - 1]);
		*ztz += z__[i__] * z__[i__];
	    }
	    if (fabs(z__[to]) <= eps && fabs(z__[to - 1]) <= eps) {
		isuppz[2] = to - 2;
	    } else {
		from = to + 1;
		i__1 = to + 32;
		to = (i__1<*bn) ? i__1 : *bn;
		goto L140;
	    }
	}
    } else {
	i__1 = *b1;
	for (i__ = *r__ - 1; i__ >= i__1; --i__) {
	    if (fabs(z__[i__ + 1])<GMX_DOUBLE_MIN) {
		z__[i__] = -(ld[i__ + 1] / ld[i__]) * z__[i__ + 2];
	    } else {
		z__[i__] = -(work[i__] * z__[i__ + 1]);
	    }
	    if (fabs(z__[i__]) <= eps && fabs(z__[i__ + 1]) <= eps) {
		isuppz[1] = i__ + 2;
		goto L170;
	    }
	    *ztz += z__[i__] * z__[i__];
	}
L170:
	i__1 = *bn - 1;
	for (i__ = *r__; i__ <= i__1; ++i__) {
	    if (fabs(z__[i__])<GMX_DOUBLE_MIN) {
		z__[i__ + 1] = -(ld[i__ - 1] / ld[i__]) * z__[i__ - 1];
	    } else {
		z__[i__ + 1] = -(work[indumn + i__] * z__[i__]);
	    }
	    if (fabs(z__[i__]) <= eps && fabs(z__[i__ + 1]) <= eps) {
		isuppz[2] = i__ - 1;
		break;
	    }
	    *ztz += z__[i__ + 1] * z__[i__ + 1];
	}
    }

    return;

}


