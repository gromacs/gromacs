#include <math.h>

#include "gmx_lapack.h"

#include <types/simple.h>

void
F77_FUNC(dlaed6,DLAED6)(int *kniter, 
                        int *orgati, 
                        double *rho, 
                        double *d__,
                        double *z__, 
                        double *finit, 
                        double *tau, 
                        int *info)
{
    int i__1;
    double r__1, r__2, r__3, r__4;

    double a, b, c__, f;
    int i__;
    double fc, df, ddf, eta, eps, base;
    int iter;
    double temp, temp1, temp2, temp3, temp4;
    int scale;
    int niter;
    double small1, small2, sminv1, sminv2, dscale[3], sclfac;
    double zscale[3], erretm;
    double safemin;
    double sclinv = 0;
    
    --z__;
    --d__;

    *info = 0;

    niter = 1;
    *tau = 0.f;
    if (*kniter == 2) {
	if (*orgati) {
	    temp = (d__[3] - d__[2]) / 2.f;
	    c__ = *rho + z__[1] / (d__[1] - d__[2] - temp);
	    a = c__ * (d__[2] + d__[3]) + z__[2] + z__[3];
	    b = c__ * d__[2] * d__[3] + z__[2] * d__[3] + z__[3] * d__[2];
	} else {
	    temp = (d__[1] - d__[2]) / 2.f;
	    c__ = *rho + z__[3] / (d__[3] - d__[2] - temp);
	    a = c__ * (d__[1] + d__[2]) + z__[1] + z__[2];
	    b = c__ * d__[1] * d__[2] + z__[1] * d__[2] + z__[2] * d__[1];
	}
        r__1 = fabs(a), r__2 = fabs(b), r__1 = ((r__1>r__2)? r__1:r__2), r__2 = fabs(c__);
        temp = (r__1>r__2) ? r__1 : r__2;
	a /= temp;
	b /= temp;
	c__ /= temp;
	if (c__ == 0.f) {
	    *tau = b / a;
	} else if (a <= 0.f) {
	    *tau = (a - sqrt((r__1 = a * a - b * 4.f * c__, fabs(r__1)))) / (
		    c__ * 2.f);
	} else {
	    *tau = b * 2.f / (a + sqrt((r__1 = a * a - b * 4.f * c__, fabs(r__1))));
	}

	temp = *rho + z__[1] / (d__[1] - *tau) + z__[2] / (d__[2] - *tau) + 
		z__[3] / (d__[3] - *tau);
	if (fabs(*finit) <= fabs(temp)) {
	    *tau = 0.f;
	}
    }

    eps = GMX_DOUBLE_EPS;
    base = 2;
    safemin = GMX_DOUBLE_MIN*(1.0+GMX_DOUBLE_EPS);
    i__1 = (int) (log(safemin) / log(base) / 3.f);
    small1 = pow(base, i__1);
    sminv1 = 1.f / small1;
    small2 = small1 * small1;
    sminv2 = sminv1 * sminv1;

    if (*orgati) {
	r__3 = (r__1 = d__[2] - *tau, fabs(r__1)), r__4 = (r__2 = d__[3] - *
		tau, fabs(r__2));
        temp = (r__3<r__4) ? r__3 : r__4;
    } else {
	r__3 = (r__1 = d__[1] - *tau, fabs(r__1)), r__4 = (r__2 = d__[2] - *
		tau, fabs(r__2));
	temp = (r__3<r__4) ? r__3 : r__4;
    }
    scale = 0;
    if (temp <= small1) {
	scale = 1;
	if (temp <= small2) {

	    sclfac = sminv2;
	    sclinv = small2;
	} else {

	    sclfac = sminv1;
	    sclinv = small1;

	}

	for (i__ = 1; i__ <= 3; ++i__) {
	    dscale[i__ - 1] = d__[i__] * sclfac;
	    zscale[i__ - 1] = z__[i__] * sclfac;
	}
	*tau *= sclfac;
    } else {

	for (i__ = 1; i__ <= 3; ++i__) {
	    dscale[i__ - 1] = d__[i__];
	    zscale[i__ - 1] = z__[i__];
	}
    }
    fc = 0.f;
    df = 0.f;
    ddf = 0.f;
    for (i__ = 1; i__ <= 3; ++i__) {
	temp = 1.f / (dscale[i__ - 1] - *tau);
	temp1 = zscale[i__ - 1] * temp;
	temp2 = temp1 * temp;
	temp3 = temp2 * temp;
	fc += temp1 / dscale[i__ - 1];
	df += temp2;
	ddf += temp3;
    }
    f = *finit + *tau * fc;

    if (fabs(f) <= 0.f) {
	goto L60;
    }
    iter = niter + 1;
    for (niter = iter; niter <= 20; ++niter) {
	if (*orgati) {
	    temp1 = dscale[1] - *tau;
	    temp2 = dscale[2] - *tau;
	} else {
	    temp1 = dscale[0] - *tau;
	    temp2 = dscale[1] - *tau;
	}
	a = (temp1 + temp2) * f - temp1 * temp2 * df;
	b = temp1 * temp2 * f;
	c__ = f - (temp1 + temp2) * df + temp1 * temp2 * ddf;
	r__1 = fabs(a), r__2 = fabs(b), r__1 = ((r__1>r__2)? r__1:r__2), r__2 = fabs(c__);
	temp = (r__1>r__2) ? r__1 : r__2;
	a /= temp;
	b /= temp;
	c__ /= temp;
	if (c__ == 0.f) {
	    eta = b / a;
	} else if (a <= 0.f) {
	    eta = (a - sqrt((r__1 = a * a - b * 4.f * c__, fabs(r__1)))) / ( c__ * 2.f);
	} else {
	    eta = b * 2.f / (a + sqrt((r__1 = a * a - b * 4.f * c__, fabs( r__1))));
	}
	if (f * eta >= 0.f) {
	    eta = -f / df;
	}
	temp = eta + *tau;
	if (*orgati) {
	    if (eta > 0.f && temp >= dscale[2]) {
		eta = (dscale[2] - *tau) / 2.f;
	    }

	    if (eta < 0.f && temp <= dscale[1]) {
		eta = (dscale[1] - *tau) / 2.f;
	    }
	} else {
	    if (eta > 0.f && temp >= dscale[1]) {
		eta = (dscale[1] - *tau) / 2.f;
	    }
	    if (eta < 0.f && temp <= dscale[0]) {
		eta = (dscale[0] - *tau) / 2.f;
	    }
	}
	*tau += eta;
	fc = 0.f;
	erretm = 0.f;
	df = 0.f;
	ddf = 0.f;
	for (i__ = 1; i__ <= 3; ++i__) {
	    temp = 1.f / (dscale[i__ - 1] - *tau);
	    temp1 = zscale[i__ - 1] * temp;
	    temp2 = temp1 * temp;
	    temp3 = temp2 * temp;
	    temp4 = temp1 / dscale[i__ - 1];
	    fc += temp4;
	    erretm += fabs(temp4);
	    df += temp2;
	    ddf += temp3;
	}
	f = *finit + *tau * fc;
	erretm = (fabs(*finit) + fabs(*tau) * erretm) * 8.f + fabs(*tau) * df;
	if (fabs(f) <= eps * erretm) {
	    goto L60;
	}
    }
    *info = 1;
L60:
    if (scale) {
	*tau *= sclinv;
    }
    return;
} 


