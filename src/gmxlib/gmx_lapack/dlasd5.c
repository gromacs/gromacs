#include <math.h>
#include "gmx_lapack.h"

void 
F77_FUNC(dlasd5,DLASD5)(int *i__, 
	double *d__, 
	double *z__, 
	double *delta, 
	double *rho, 
	double *dsigma, 
	double *work)
{
    double b, c__, w, del, tau, delsq;

    --work;
    --delta;
    --z__;
    --d__;

    del = d__[2] - d__[1];
    delsq = del * (d__[2] + d__[1]);
    if (*i__ == 1) {
	w = *rho * 4. * (z__[2] * z__[2] / (d__[1] + d__[2] * 3.) - z__[1] * 
		z__[1] / (d__[1] * 3. + d__[2])) / del + 1.;
	if (w > 0.) {
	    b = delsq + *rho * (z__[1] * z__[1] + z__[2] * z__[2]);
	    c__ = *rho * z__[1] * z__[1] * delsq;

	    tau = c__ * 2. / (b + sqrt(fabs(b * b - c__ * 4.)));

	    tau /= d__[1] + sqrt(d__[1] * d__[1] + tau);
	    *dsigma = d__[1] + tau;
	    delta[1] = -tau;
	    delta[2] = del - tau;
	    work[1] = d__[1] * 2. + tau;
	    work[2] = d__[1] + tau + d__[2];
	} else {
	    b = -delsq + *rho * (z__[1] * z__[1] + z__[2] * z__[2]);
	    c__ = *rho * z__[2] * z__[2] * delsq;

	    if (b > 0.) {
		tau = c__ * -2. / (b + sqrt(b * b + c__ * 4.));
	    } else {
		tau = (b - sqrt(b * b + c__ * 4.)) / 2.;
	    }

	    tau /= d__[2] + sqrt(fabs(d__[2] * d__[2] + tau));
	    *dsigma = d__[2] + tau;
	    delta[1] = -(del + tau);
	    delta[2] = -tau;
	    work[1] = d__[1] + tau + d__[2];
	    work[2] = d__[2] * 2. + tau;
	}
    } else {

	b = -delsq + *rho * (z__[1] * z__[1] + z__[2] * z__[2]);
	c__ = *rho * z__[2] * z__[2] * delsq;

	if (b > 0.) {
	    tau = (b + sqrt(b * b + c__ * 4.)) / 2.;
	} else {
	    tau = c__ * 2. / (-b + sqrt(b * b + c__ * 4.));
	}
	tau /= d__[2] + sqrt(d__[2] * d__[2] + tau);
	*dsigma = d__[2] + tau;
	delta[1] = -(del + tau);
	delta[2] = -tau;
	work[1] = d__[1] + tau + d__[2];
	work[2] = d__[2] * 2. + tau;
    }
    return;

} 
