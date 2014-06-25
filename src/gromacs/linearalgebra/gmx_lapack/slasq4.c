#include <math.h>
#include "gromacs/utility/real.h"

#include "../gmx_lapack.h"

void 
F77_FUNC(slasq4,SLASQ4)(int *i0, 
	int *n0, 
	float *z__, 
	int *pp, 
	int *n0in, 
	float *dmin__, 
	float *dmin1, 
	float *dmin2, 
	float *dn, 
	float *dn1, 
	float *dn2, 
	float *tau, 
	int *ttype)
{
    float g = 0.;
    int i__1;
    float d__1, d__2;

    float s, a2, b1, b2;
    int i4, nn, np;
    float gam, gap1, gap2;


    if (*dmin__ <= 0.) {
	*tau = -(*dmin__);
	*ttype = -1;
	return;
    }

    s = 0.0;

    nn = (*n0 << 2) + *pp;
    if (*n0in == *n0) {

	if ( fabs(*dmin__ - *dn)<GMX_FLOAT_EPS*fabs(*dmin__ + *dn) ||
         fabs(*dmin__ - *dn1)<GMX_FLOAT_EPS*fabs(*dmin__ + *dn1)) {

	    b1 = sqrt(z__[nn - 3]) * sqrt(z__[nn - 5]);
	    b2 = sqrt(z__[nn - 7]) * sqrt(z__[nn - 9]);
	    a2 = z__[nn - 7] + z__[nn - 5];

        if ( fabs(*dmin__ - *dn)<GMX_FLOAT_EPS*fabs(*dmin__ + *dn) &&
             fabs(*dmin1 - *dn1)<GMX_FLOAT_EPS*fabs(*dmin1 + *dn1)) {

            gap2 = *dmin2 - a2 - *dmin2 * .25;
		if (gap2 > 0. && gap2 > b2) {
		    gap1 = a2 - *dn - b2 / gap2 * b2;
		} else {
		    gap1 = a2 - *dn - (b1 + b2);
		}
		if (gap1 > 0. && gap1 > b1) {
		    d__1 = *dn - b1 / gap1 * b1, d__2 = *dmin__ * .5;
		    s = (d__1>d__2) ? d__1 : d__2;
		    *ttype = -2;
		} else {
		    s = 0.;
		    if (*dn > b1) {
			s = *dn - b1;
		    }
		    if (a2 > b1 + b2) {
			d__1 = s, d__2 = a2 - (b1 + b2);
			s = (d__1<d__2) ? d__1 : d__2;
		    }
		    d__1 = s, d__2 = *dmin__ * .333;
		    s = (d__1>d__2) ? d__1 : d__2;
		    *ttype = -3;
		}
	    } else {


		*ttype = -4;
		s = *dmin__ * .25;
		if (fabs(*dmin__ - *dn)<GMX_FLOAT_EPS*fabs(*dmin__ + *dn)) {
		    gam = *dn;
		    a2 = 0.;
		    if (z__[nn - 5] > z__[nn - 7]) {
			return;
		    }
		    b2 = z__[nn - 5] / z__[nn - 7];
		    np = nn - 9;
		} else {
		    np = nn - (*pp << 1);
		    b2 = z__[np - 2];
		    gam = *dn1;
		    if (z__[np - 4] > z__[np - 2]) {
			return;
		    }
		    a2 = z__[np - 4] / z__[np - 2];
		    if (z__[nn - 9] > z__[nn - 11]) {
			return;
		    }
		    b2 = z__[nn - 9] / z__[nn - 11];
		    np = nn - 13;
		}


		a2 += b2;
		i__1 = (*i0 << 2) - 1 + *pp;
		for (i4 = np; i4 >= i__1; i4 += -4) {
		    if (fabs(b2)<GMX_FLOAT_MIN) {
			goto L20;
		    }
		    b1 = b2;
		    if (z__[i4] > z__[i4 - 2]) {
			return;
		    }
		    b2 *= z__[i4] / z__[i4 - 2];
		    a2 += b2;
		    if (((b2>b1) ? b2 : b1) * 100. < a2 || .563 < a2) {
			goto L20;
		    }
		}
L20:
		a2 *= 1.05;


		if (a2 < .563) {
		    s = gam * (1. - sqrt(a2)) / (a2 + 1.);
		}
	    }
	} else if (fabs(*dmin__ - *dn2)<GMX_FLOAT_EPS*fabs(*dmin__ + *dn2)) {

	    *ttype = -5;
	    s = *dmin__ * .25;

	    np = nn - (*pp << 1);
	    b1 = z__[np - 2];
	    b2 = z__[np - 6];
	    gam = *dn2;
	    if (z__[np - 8] > b2 || z__[np - 4] > b1) {
		return;
	    }
	    a2 = z__[np - 8] / b2 * (z__[np - 4] / b1 + 1.);


	    if (*n0 - *i0 > 2) {
		b2 = z__[nn - 13] / z__[nn - 15];
		a2 += b2;
		i__1 = (*i0 << 2) - 1 + *pp;
		for (i4 = nn - 17; i4 >= i__1; i4 += -4) {
		    if (fabs(b2)<GMX_FLOAT_MIN) {
			goto L40;
		    }
		    b1 = b2;
		    if (z__[i4] > z__[i4 - 2]) {
			return;
		    }
		    b2 *= z__[i4] / z__[i4 - 2];
		    a2 += b2;
		    if (((b2>b1) ? b2 : b1) * 100. < a2 || .563 < a2) {
			goto L40;
		    }
		}
L40:
		a2 *= 1.05;
	    }

	    if (a2 < .563) {
		s = gam * (1. - sqrt(a2)) / (a2 + 1.);
	    }
	} else {

	    if (*ttype == -6) {
		g += (1. - g) * .333;
	    } else if (*ttype == -18) {
		g = .083250000000000005;
	    } else {
		g = .25;
	    }
	    s = g * *dmin__;
	    *ttype = -6;
	}

    } else if (*n0in == *n0 + 1) {

        if ( fabs(*dmin1 - *dn1)<GMX_FLOAT_EPS*fabs(*dmin1 + *dn1) &&
             fabs(*dmin2 - *dn2)<GMX_FLOAT_EPS*fabs(*dmin2 + *dn2)) {

	    *ttype = -7;
	    s = *dmin1 * .333;
	    if (z__[nn - 5] > z__[nn - 7]) {
		return;
	    }
	    b1 = z__[nn - 5] / z__[nn - 7];
	    b2 = b1;
	    if (fabs(b2)<GMX_FLOAT_MIN) {
		goto L60;
	    }
	    i__1 = (*i0 << 2) - 1 + *pp;
	    for (i4 = (*n0 << 2) - 9 + *pp; i4 >= i__1; i4 += -4) {
		a2 = b1;
		if (z__[i4] > z__[i4 - 2]) {
		    return;
		}
		b1 *= z__[i4] / z__[i4 - 2];
		b2 += b1;
		if (((a2>b1) ? a2 : b1) * 100. < b2) {
		    goto L60;
		}
	    }
L60:
	    b2 = sqrt(b2 * 1.05);
	    d__1 = b2;
	    a2 = *dmin1 / (d__1 * d__1 + 1.);
	    gap2 = *dmin2 * .5 - a2;
	    if (gap2 > 0. && gap2 > b2 * a2) {
		d__1 = s, d__2 = a2 * (1. - a2 * 1.01 * (b2 / gap2) * b2);
		s = (d__1>d__2) ? d__1 : d__2;
	    } else {
		d__1 = s, d__2 = a2 * (1. - b2 * 1.01);
		s = (d__1>d__2) ? d__1 : d__2;
		*ttype = -8;
	    }
	} else {

	    s = *dmin1 * .25;
	    if (fabs(*dmin1 - *dn1)<GMX_FLOAT_EPS*fabs(*dmin1 + *dn1)) {
		s = *dmin1 * .5;
	    }
	    *ttype = -9;
	}

    } else if (*n0in == *n0 + 2) {

	if (fabs(*dmin2 - *dn2)<GMX_FLOAT_EPS*fabs(*dmin2 + *dn2) &&
        z__[nn - 5] * 2. < z__[nn - 7]) {
	    *ttype = -10;
	    s = *dmin2 * .333;
	    if (z__[nn - 5] > z__[nn - 7]) {
		return;
	    }
	    b1 = z__[nn - 5] / z__[nn - 7];
	    b2 = b1;
	    if (fabs(b2)<GMX_FLOAT_MIN) {
		goto L80;
	    }
	    i__1 = (*i0 << 2) - 1 + *pp;
	    for (i4 = (*n0 << 2) - 9 + *pp; i4 >= i__1; i4 += -4) {
		if (z__[i4] > z__[i4 - 2]) {
		    return;
		}
		b1 *= z__[i4] / z__[i4 - 2];
		b2 += b1;
		if (b1 * 100. < b2) {
		    goto L80;
		}
	    }
L80:
	    b2 = sqrt(b2 * 1.05);
	    d__1 = b2;
	    a2 = *dmin2 / (d__1 * d__1 + 1.);
	    gap2 = z__[nn - 7] + z__[nn - 9] - sqrt(z__[nn - 11]) * sqrt(z__[
		    nn - 9]) - a2;
	    if (gap2 > 0. && gap2 > b2 * a2) {
		d__1 = s, d__2 = a2 * (1. - a2 * 1.01 * (b2 / gap2) * b2);
		s = (d__1>d__2) ? d__1 : d__2;
	    } else {
		d__1 = s, d__2 = a2 * (1. - b2 * 1.01);
		s = (d__1>d__2) ? d__1 : d__2;
	    }
	} else {
	    s = *dmin2 * .25;
	    *ttype = -11;
	}
    } else if (*n0in > *n0 + 2) {

	s = 0.;
	*ttype = -12;
    }

    *tau = s;
    return;

}


