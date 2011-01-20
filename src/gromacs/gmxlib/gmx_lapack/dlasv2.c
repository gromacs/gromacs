#include <math.h>
#include "gmx_lapack.h"
#include "lapack_limits.h"

#include <types/simple.h>

void 
F77_FUNC(dlasv2,DLASV2)(double *f, 
                        double *g, 
                        double *h__, 
                        double *ssmin, 
                        double *ssmax, 
                        double *snr, 
                        double *csr, 
                        double *snl, 
                        double *csl)
{
    double d__1;

    double a, d__, l, m, r__, s, t, fa, ga, ha, ft, gt, ht, mm, tt,
	     clt, crt, slt, srt;
    int pmax;
    double temp;
    int swap;
    double tsign=1.0;
    int gasmal;

    ft = *f;
    fa = fabs(ft);
    ht = *h__;
    ha = fabs(*h__);

    pmax = 1;
    swap = ha > fa;
    if (swap) {
	pmax = 3;
	temp = ft;
	ft = ht;
	ht = temp;
	temp = fa;
	fa = ha;
	ha = temp;

    }
    gt = *g;
    ga = fabs(gt);
    if (fabs(ga)<GMX_DOUBLE_MIN) {

	*ssmin = ha;
	*ssmax = fa;
	clt = 1.;
	crt = 1.;
	slt = 0.;
	srt = 0.;
    } else {
	gasmal = 1;
	if (ga > fa) {
	    pmax = 2;
	    if (fa / ga < GMX_DOUBLE_EPS) {

		gasmal = 0;
		*ssmax = ga;
		if (ha > 1.) {
		    *ssmin = fa / (ga / ha);
		} else {
		    *ssmin = fa / ga * ha;
		}
		clt = 1.;
		slt = ht / gt;
		srt = 1.;
		crt = ft / gt;
	    }
	}
	if (gasmal) {

	    d__ = fa - ha;
	    if ( fabs( fa - d__ )<GMX_DOUBLE_EPS*fabs( fa + d__ )) {
		l = 1.;
	    } else {
		l = d__ / fa;
	    }

	    m = gt / ft;
	    t = 2. - l;

	    mm = m * m;
	    tt = t * t;
	    s = sqrt(tt + mm);

	    if ( fabs(l)<GMX_DOUBLE_MIN) {
		r__ = fabs(m);
	    } else {
		r__ = sqrt(l * l + mm);
	    }
	    a = (s + r__) * .5;

	    *ssmin = ha / a;
	    *ssmax = fa * a;
	    if ( fabs(mm)<GMX_DOUBLE_MIN) {

		if (fabs(l)<GMX_DOUBLE_MIN) {
		    t = ( (ft>0) ? 2.0 : -2.0) * ( (gt>0) ? 1.0 : -1.0);
		} else {
		    t = gt / ( (ft>0) ? d__ : d__) + m / t;
		}
	    } else {
		t = (m / (s + t) + m / (r__ + l)) * (a + 1.);
	    }
	    l = sqrt(t * t + 4.);
	    crt = 2. / l;
	    srt = t / l;
	    clt = (crt + srt * m) / a;
	    slt = ht / ft * srt / a;
	}
    }
    if (swap) {
	*csl = srt;
	*snl = crt;
	*csr = slt;
	*snr = clt;
    } else {
	*csl = clt;
	*snl = slt;
	*csr = crt;
	*snr = srt;
    }

    if (pmax == 1) {
	tsign = ( (*csr>0) ? 1.0 : -1.0) * ( (*csl>0) ? 1.0 : -1.0) * ( (*f>0) ? 1.0 : -1.0);
    }
    if (pmax == 2) {
	tsign = ( (*snr>0) ? 1.0 : -1.0) * ( (*csl>0) ? 1.0 : -1.0) * ( (*g>0) ? 1.0 : -1.0);
    }
    if (pmax == 3) {
	tsign = ( (*snr>0) ? 1.0 : -1.0) * ( (*snl>0) ? 1.0 : -1.0) * ( (*h__>0) ? 1.0 : -1.0);
    }
    if(tsign<0)
      *ssmax *= -1.0;
    d__1 = tsign * ( (*f>0) ? 1.0 : -1.0) * ( (*h__>0) ? 1.0 : -1.0);
    if(d__1<0)
      *ssmin *= -1.0;
    return;

}
