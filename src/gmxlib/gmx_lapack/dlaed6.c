#include <math.h>
#include "gmx_lapack.h"
#include "lapack_limits.h"

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
    double d__1, d__2, d__3, d__4;

    double a, b, c__, f;
    int i__;
    double fc, df, ddf, eta;
    int iter;
    double temp, temp1, temp2, temp3, temp4;
    int scale;
    int niter;
    double dscale[3], sclfac, zscale[3], erretm, sclinv;
    const double safemin = 
      (1.0+LAPACK_EPS_DOUBLE)/LAPACK_MAX_DOUBLE/LAPACK_EPS_DOUBLE;
    const double small1 = ( 1 << ((int)(log(safemin)/log(2.0)/3.0)));
    const double sminv1 = 1.0/small1;
    const double small2 = small1*small1;
    const double sminv2 = sminv1*sminv1;

    --z__;
    --d__;


    *info = 0;

    sclinv = 1.0;
    niter = 1;
    *tau = 0.;
    if (*kniter == 2) {
	if (*orgati) {
	    temp = (d__[3] - d__[2]) / 2.;
	    c__ = *rho + z__[1] / (d__[1] - d__[2] - temp);
	    a = c__ * (d__[2] + d__[3]) + z__[2] + z__[3];
	    b = c__ * d__[2] * d__[3] + z__[2] * d__[3] + z__[3] * d__[2];
	} else {
	    temp = (d__[1] - d__[2]) / 2.;
	    c__ = *rho + z__[3] / (d__[3] - d__[2] - temp);
	    a = c__ * (d__[1] + d__[2]) + z__[1] + z__[2];
	    b = c__ * d__[1] * d__[2] + z__[1] * d__[2] + z__[2] * d__[1];
	}
	d__1 = fabs(a);
	d__2 = fabs(b);
	if(d__2>d__1)
	   d__1 = d__2;
	d__2 = fabs(c__);
	temp = (d__1>d__2) ? d__1 : d__2;
	a /= temp;
	b /= temp;
	c__ /= temp;
	if (c__ == 0.) {
	    *tau = b / a;
	} else if (a <= 0.) {
	    *tau = (a - sqrt(fabs(a * a - b * 4. * c__))) / ( c__ * 2.);
	} else {
	  *tau = b * 2. / (a + sqrt(fabs(a * a - b * 4. * c__)));
	}
	temp = *rho + z__[1] / (d__[1] - *tau) + z__[2] / (d__[2] - *tau) + 
		z__[3] / (d__[3] - *tau);
	if (fabs(*finit) <= fabs(temp)) {
	    *tau = 0.;
	}
    }

    if (*orgati) {
	d__3 = fabs(d__[2] - *tau);
	d__4 = fabs(d__[3] - *tau);
	temp = (d__3<d__4) ? d__3 : d__4;
    } else {
	d__3 = fabs(d__[1] - *tau);
	d__4 = fabs(d__[2] - *tau);
	temp = (d__3<d__4) ? d__3 : d__4;
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

    fc = 0.;
    df = 0.;
    ddf = 0.;
    for (i__ = 1; i__ <= 3; ++i__) {
	temp = 1. / (dscale[i__ - 1] - *tau);
	temp1 = zscale[i__ - 1] * temp;
	temp2 = temp1 * temp;
	temp3 = temp2 * temp;
	fc += temp1 / dscale[i__ - 1];
	df += temp2;
	ddf += temp3;
    }
    f = *finit + *tau * fc;

    if (fabs(f) <= 0.) {
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
	d__1 = fabs(a);
	d__2 = fabs(b);
	if(d__2>d__1)
	  d__1 = d__2;
	d__2 = fabs(c__);
	temp = (d__1>d__2) ? d__1 : d__2;
	a /= temp;
	b /= temp;
	c__ /= temp;
	if (c__ == 0.) {
	    eta = b / a;
	} else if (a <= 0.) {
	    eta = (a - sqrt(fabs(a * a - b * 4. * c__))) / (c__ * 2.);
	} else {
	  eta = b * 2. / (a + sqrt(fabs(a * a - b * 4. * c__)));
	}
	if (f * eta >= 0.) {
	    eta = -f / df;
	}

	temp = eta + *tau;
	if (*orgati) {
	    if (eta > 0. && temp >= dscale[2]) {
		eta = (dscale[2] - *tau) / 2.;
	    }
	    if (eta < 0. && temp <= dscale[1]) {
		eta = (dscale[1] - *tau) / 2.;
	    }
	} else {
	    if (eta > 0. && temp >= dscale[1]) {
		eta = (dscale[1] - *tau) / 2.;
	    }
	    if (eta < 0. && temp <= dscale[0]) {
		eta = (dscale[0] - *tau) / 2.;
	    }
	}
	*tau += eta;

	fc = 0.;
	erretm = 0.;
	df = 0.;
	ddf = 0.;
	for (i__ = 1; i__ <= 3; ++i__) {
	    temp = 1. / (dscale[i__ - 1] - *tau);
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
	erretm = (fabs(*finit) + fabs(*tau) * erretm) * 8. + fabs(*tau) * df;
	if (fabs(f) <= LAPACK_EPS_DOUBLE * erretm) {
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


