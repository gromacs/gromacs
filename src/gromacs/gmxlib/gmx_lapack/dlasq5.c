#include <math.h>
#include "gmx_lapack.h"

void
F77_FUNC(dlasq5,DLASQ5)(int *i0, 
	int *n0,
	double *z__, 
	int *pp, 
	double *tau,
	double *dmin__, 
	double *dmin1, 
	double *dmin2, 
	double *dn,
	double *dnm1, 
	double *dnm2,
	int *ieee)
{
    int i__1;
    double d__1, d__2;

    static double d__;
    static int j4, j4p2;
    static double emin, temp;

    --z__;

    if (*n0 - *i0 - 1 <= 0) {
	return;
    }

    j4 = (*i0 << 2) + *pp - 3;
    emin = z__[j4 + 4];
    d__ = z__[j4] - *tau;
    *dmin__ = d__;
    *dmin1 = -z__[j4];

    if (*ieee) {

	if (*pp == 0) {
	    i__1 = 4*(*n0 - 3);
	    for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
		z__[j4 - 2] = d__ + z__[j4 - 1];
		temp = z__[j4 + 1] / z__[j4 - 2];
		d__ = d__ * temp - *tau;
                if(d__<*dmin__)
                  *dmin__ = d__;
		z__[j4] = z__[j4 - 1] * temp;
		d__1 = z__[j4];
                if(d__1<emin)
                  emin = d__1;
	    }
	} else {
	    i__1 = 4*(*n0 - 3);
	    for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
		z__[j4 - 3] = d__ + z__[j4];
		temp = z__[j4 + 2] / z__[j4 - 3];
		d__ = d__ * temp - *tau;
                if(d__<*dmin__)
                  *dmin__ = d__;
		z__[j4 - 1] = z__[j4] * temp;
		d__1 = z__[j4 - 1];
                if(d__1<emin)
                  emin = d__1;
	    }
	}

	*dnm2 = d__;
	*dmin2 = *dmin__;
	j4 = 4*(*n0 - 2) - *pp;
	j4p2 = j4 + (*pp << 1) - 1;
	z__[j4 - 2] = *dnm2 + z__[j4p2];
	z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
	*dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]) - *tau;
        if(*dnm1<*dmin__)
          *dmin__ = *dnm1;

	*dmin1 = *dmin__;
	j4 += 4;
	j4p2 = j4 + (*pp << 1) - 1;
	z__[j4 - 2] = *dnm1 + z__[j4p2];
	z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
	*dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]) - *tau;
        if(*dn<*dmin__)
          *dmin__ = *dn;

    } else {

	if (*pp == 0) {
	    i__1 = 4*(*n0 - 3);
	    for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
		z__[j4 - 2] = d__ + z__[j4 - 1];
		if (d__ < 0.) {
		    return;
		} else {
		    z__[j4] = z__[j4 + 1] * (z__[j4 - 1] / z__[j4 - 2]);
		    d__ = z__[j4 + 1] * (d__ / z__[j4 - 2]) - *tau;
		}
                if(d__<*dmin__)
                  *dmin__ = d__;
		d__1 = emin, d__2 = z__[j4];
		emin = (d__1<d__2) ? d__1 : d__2;
	    }
	} else {
	    i__1 = 4*(*n0 - 3);
	    for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
		z__[j4 - 3] = d__ + z__[j4];
		if (d__ < 0.) {
		    return;
		} else {
		    z__[j4 - 1] = z__[j4 + 2] * (z__[j4] / z__[j4 - 3]);
		    d__ = z__[j4 + 2] * (d__ / z__[j4 - 3]) - *tau;
		}
                if(d__<*dmin__)
                  *dmin__ = d__;
		d__1 = emin, d__2 = z__[j4 - 1];
		emin = (d__1<d__2) ? d__1 : d__2;
	    }
	}

	*dnm2 = d__;
	*dmin2 = *dmin__;
	j4 = 4*(*n0 - 2) - *pp;
	j4p2 = j4 + (*pp << 1) - 1;
	z__[j4 - 2] = *dnm2 + z__[j4p2];
	if (*dnm2 < 0.) {
	    return;
	} else {
	    z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
	    *dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]) - *tau;
	}
        if(*dnm1<*dmin__)
          *dmin__ = *dnm1;

	*dmin1 = *dmin__;
	j4 += 4;
	j4p2 = j4 + (*pp << 1) - 1;
	z__[j4 - 2] = *dnm1 + z__[j4p2];
	if (*dnm1 < 0.) {
	    return;
	} else {
	    z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
	    *dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]) - *tau;
	}
        if(*dn<*dmin__)
          *dmin__ = *dn;

    }

    z__[j4 + 2] = *dn;
    z__[(*n0 << 2) - *pp] = emin;
    return;

}

