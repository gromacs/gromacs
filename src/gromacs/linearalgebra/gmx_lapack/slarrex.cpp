#include <cctype>
#include <cmath>

#include "gromacs/utility/real.h"

#include "../gmx_blas.h"
#include "../gmx_lapack.h"
#include "lapack_limits.h"



void
F77_FUNC(slarrex,SLARREX)(const char *range,
	 int *n, 
	 float *vl, 
	 float *vu, 
	 int *il, 
	 int *iu, 
	 float *d__, 
	 float *e, 
	 float *tol, 
	 int *nsplit, 
	 int *isplit, 
	 int *m, 
	 float *w, 
	 int *iblock, 
	 int *indexw, 
	 float *gersch, 
	 float *work,
	 int *iwork, 
	 int *info)
{
    int i__1, i__2, i__3;
    float d__1, d__2;
    int c__1 = 1;
    int c__0 = 0;

    int i__, j, k;
    float s, gl;
    int in;
    float gu;
    int cnt;
    float eps, tau, nrm, tmp, vvl, vvu, offd;
    int iend, jblk, till, itmp;
    float rtol, delta, sigma;
    int iinfo;
    float width;
    int ibegin;
    int irange;
    float sgndef;
    int maxcnt;
    --iwork;
    --work;
    --gersch;
    --indexw;
    --iblock;
    --w;
    --isplit;
    --e;
    --d__;

    sigma = 0;
    irange = 0;
    sgndef = 0;
    maxcnt = 0;

    *info = 0;

    if (*range=='A' || *range=='a')
	irange = 1;
    else if (*range=='V' || *range=='v')
	irange = 2;
    else if (*range=='I' || *range=='i')
	irange = 3;
    

    *m = 0;
    eps = GMX_FLOAT_EPS;

    *nsplit = 1;
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (std::abs(e[i__]) <= *tol) {
	    isplit[*nsplit] = i__;
	    ++(*nsplit);
	}
    }
    isplit[*nsplit] = *n;

    ibegin = 1;
    i__1 = *nsplit;
    for (jblk = 1; jblk <= i__1; ++jblk) {
	iend = isplit[jblk];
	if (ibegin == iend) {
	    ++(*m);
	    w[*m] = d__[ibegin];
	    iblock[*m] = jblk;
	    indexw[*m] = 1;
	    e[iend] = 0.;
	    ibegin = iend + 1;
	    goto L170;
	}
	in = iend - ibegin + 1;

	gl = d__[ibegin] - std::abs(e[ibegin]);
	gu = d__[ibegin] + std::abs(e[ibegin]);
	gersch[(ibegin << 1) - 1] = gl;
	gersch[ibegin * 2] = gu;
	gersch[(iend << 1) - 1] = d__[iend] - std::abs(e[iend - 1]);
	gersch[iend * 2] = d__[iend] + std::abs(e[iend - 1]);
	d__1 = gersch[(iend << 1) - 1];
	gl = (d__1<gl) ? d__1 : gl;
	d__1 = gersch[iend * 2];
	gu = (d__1>gu) ? d__1 : gu;
	i__2 = iend - 1;
	for (i__ = ibegin + 1; i__ <= i__2; ++i__) {
	    offd = std::abs(e[i__ - 1]) + std::abs(e[i__]);
	    gersch[(i__ << 1) - 1] = d__[i__] - offd;
	    d__1 = gersch[(i__ << 1) - 1];
	    gl = (d__1<gl) ? d__1 : gl;
	    gersch[i__ * 2] = d__[i__] + offd;
	    d__1 = gersch[i__ * 2];
	    gu = (d__1>gu) ? d__1 : gu;
	}
	d__1 = std::abs(gl), d__2 = std::abs(gu);
	nrm = (d__1>d__2) ? d__1 : d__2;

	width = gu - gl;
	i__2 = iend - 1;
	for (i__ = ibegin; i__ <= i__2; ++i__) {
	    work[i__] = e[i__] * e[i__];
	}
	for (j = 1; j <= 2; ++j) {
	    if (j == 1) {
		tau = gl + width * .25;
	    } else {
		tau = gu - width * .25;
	    }
	    tmp = d__[ibegin] - tau;
	    if (tmp < 0.) {
		cnt = 1;
	    } else {
		cnt = 0;
	    }
	    i__2 = iend;
	    for (i__ = ibegin + 1; i__ <= i__2; ++i__) {
		tmp = d__[i__] - tau - work[i__ - 1] / tmp;
		if (tmp < 0.) {
		    ++cnt;
		}
	    }
	    if (cnt == 0) {
		gl = tau;
	    } else if (cnt == in) {
		gu = tau;
	    }
	    if (j == 1) {
		maxcnt = cnt;
		sigma = gl;
		sgndef = 1.;
	    } else {
		if (in - cnt > maxcnt) {
		    sigma = gu;
		    sgndef = -1.;
		}
	    }
	}

	work[in * 3] = 1.;
	delta = eps;
	tau = sgndef * nrm;
L60:
	sigma -= delta * tau;
	work[1] = d__[ibegin] - sigma;
	j = ibegin;
	i__2 = in - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work[(in << 1) + i__] = 1. / work[i__];
	    tmp = e[j] * work[(in << 1) + i__];
	    work[i__ + 1] = d__[j + 1] - sigma - tmp * e[j];
	    work[in + i__] = tmp;
	    ++j;
	}
	for (i__ = in; i__ >= 1; --i__) {
	    tmp = sgndef * work[i__];
	    if (tmp < 0. || std::abs(work[(in << 1) + i__])<GMX_FLOAT_MIN || std::isnan(tmp)) {
		delta *= 2.;
		goto L60;
	    }
	}

	F77_FUNC(scopy,SCOPY)(&in, &work[1], &c__1, &d__[ibegin], &c__1);
	i__2 = in - 1;
	F77_FUNC(scopy,SCOPY)(&i__2, &work[in + 1], &c__1, &e[ibegin], &c__1);
	i__2 = in - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work[in * 3 + i__] = work[i__] * work[in + i__];
	    work[(in << 2) + i__] = work[in * 3 + i__] * work[in + i__];
	}
	if (sgndef > 0.) {
	    cnt = 1;
	    work[1] = (gl + gu) / 2. - sigma;
	    work[in + 1] = 0.;
	    work[(in << 1) + 1] = (gu - gl) / 2.;
	} else {
	    cnt = in;
	    work[in] = (gl + gu) / 2. - sigma;
	    work[in * 2] = 0.;
	    work[in * 3] = (gu - gl) / 2.;
	}
	rtol = eps * 4.;
	F77_FUNC(slarrbx,SLARRBX)(&in, &d__[ibegin], &e[ibegin], &work[in * 3 + 1], &work[(in <<
		 2) + 1], &cnt, &cnt, &rtol, &rtol, &c__0, &work[1], &work[in 
		+ 1], &work[(in << 1) + 1], &work[in * 5 + 1], &iwork[1], &
		iinfo);
	if (sgndef > 0.) {
	    tau = work[1] - work[(in << 1) + 1];
	} else {
	    tau = work[in] + work[in * 3];
	}

	work[in * 3] = 1.;
	delta = eps * 2.;
L100:
	tau *= 1. - delta;

	s = -tau;
	j = ibegin;
	i__2 = in - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work[i__] = d__[j] + s;
	    work[(in << 1) + i__] = 1. / work[i__];
	    work[in + i__] = e[j] * d__[j] * work[(in << 1) + i__];
	    s = s * work[in + i__] * e[j] - tau;
	    ++j;
	}
	work[in] = d__[iend] + s;

	for (i__ = in; i__ >= 1; --i__) {
	    tmp = sgndef * work[i__];
	    if (tmp < 0. || std::abs(work[(in << 1) + i__])<GMX_FLOAT_MIN || std::isnan(tmp)) {
		delta *= 2.;
		goto L100;
	    }
	}

	sigma += tau;
	F77_FUNC(scopy,SCOPY)(&in, &work[1], &c__1, &d__[ibegin], &c__1);
	i__2 = in - 1;
	F77_FUNC(scopy,SCOPY)(&i__2, &work[in + 1], &c__1, &e[ibegin], &c__1);
	e[iend] = sigma;
	tmp = (float) in * 4. * eps * (std::abs(sigma) + std::abs(tau));
	i__2 = iend;
	for (i__ = ibegin; i__ <= i__2; ++i__) {
	    gersch[(i__ << 1) - 1] = gersch[(i__ << 1) - 1] - sigma - tmp;
	    gersch[i__ * 2] = gersch[i__ * 2] - sigma + tmp;
	}

	j = ibegin;
	i__2 = in - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work[(i__ << 1) - 1] = std::abs(d__[j]);
	    work[i__ * 2] = e[j] * e[j] * work[(i__ << 1) - 1];
	    ++j;
	}
	work[(in << 1) - 1] = std::abs(d__[iend]);

	F77_FUNC(slasq2,SLASQ2)(&in, &work[1], info);
	if (*info != 0) {
	    return;
	}

	if (sgndef > 0.) {
	    i__2 = in;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		++(*m);
		w[*m] = work[in - i__ + 1];
		iblock[*m] = jblk;
		indexw[*m] = i__;
	    }
	} else {
	    i__2 = in;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		++(*m);
		w[*m] = -work[i__];
		iblock[*m] = jblk;
		indexw[*m] = i__;
	    }
	}
	ibegin = iend + 1;
L170:
	;
    }
    if (irange == 2) {
	*m = 0;
	ibegin = 1;
	i__1 = *nsplit;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    iend = isplit[i__];
	    vvl = *vl - e[iend];
	    vvu = *vu - e[iend];
	    i__2 = iend;
	    for (j = ibegin; j <= i__2; ++j) {
		if (vvl <= w[j] && w[j] <= vvu) {
		    ++(*m);
		    w[*m] = w[j];
		    iblock[*m] = i__;
		    indexw[*m] = j - ibegin + 1;
		}
	    }
	    ibegin = iend + 1;
	}
    } else if (irange == 3) {
	*m = *iu - *il + 1;
	if (*nsplit == 1) {
	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		w[i__] = w[*il + i__ - 1];
		indexw[i__] = *il + i__ - 1;
	    }
	} else {
	    ibegin = 1;
	    i__1 = *nsplit;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		iend = isplit[i__];
		i__2 = iend;
		for (j = ibegin; j <= i__2; ++j) {
		    work[j] = w[j] + e[iend];
		}
		ibegin = iend + 1;
	    }
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		iwork[i__] = i__;
		iwork[*n + i__] = iblock[i__];
	    }
	    F77_FUNC(slasrt2,SLASRT2)("I", n, &work[1], &iwork[1], &iinfo);
	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		itmp = iwork[*il + i__ - 1];
		work[i__] = w[itmp];
		iblock[i__] = iwork[*n + itmp];
	    }
	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		iwork[*n + i__] = iwork[*il + i__ - 1];
		iwork[i__] = i__;
	    }
	    F77_FUNC(ilasrt2,ILASRT2)("I", m, &iblock[1], &iwork[1], &iinfo);
	    j = 1;
	    itmp = iblock[j];
	    cnt = iwork[*n + iwork[j]];
	    if (itmp == 1) {
		ibegin = 1;
	    } else {
		ibegin = isplit[itmp - 1] + 1;
	    }
	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		w[i__] = work[iwork[i__]];
		if (iblock[i__] != itmp || i__ == *m) {
		    if (iblock[i__] == itmp) {
			till = *m;
		    } else {
			till = i__ - 1;
		    }
		    i__2 = till - j + 1;
		    F77_FUNC(slasrt,SLASRT)("I", &i__2, &w[j], &iinfo);
		    cnt = cnt - ibegin + 1;
		    i__2 = till;
		    for (k = j; k <= i__2; ++k) {
			indexw[k] = cnt + k - j;
		    }
		    j = i__;
		    itmp = iblock[j];
		    cnt = iwork[*n + iwork[j]];
		    ibegin = isplit[itmp - 1] + 1;
		    if (i__ == *m && till < *m) {
			indexw[*m] = cnt - ibegin + 1;
		    }
		} else {
		    i__2 = cnt, i__3 = iwork[*n + iwork[i__]];
		    cnt = (i__2<i__3) ? i__2 : i__3;
		}
	    }
	}
    }

    return;

}


