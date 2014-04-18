#include <math.h>
#include "../gmx_lapack.h"
#include "lapack_limits.h"

#include "gromacs/utility/real.h"

void 
F77_FUNC(slasq6,SLASQ6)(int *i0, 
	int *n0, 
	float *z__, 
	int *pp, 
	float *dmin__, 
	float *dmin1, 
	float *dmin2,
	float *dn, 
	float *dnm1, 
	float *dnm2)
{
    int i__1;
    float d__1, d__2;

    /* Local variables */
    float d__;
    int j4, j4p2;
    float emin, temp;
    const float safemin = GMX_FLOAT_MIN*(1.0+GMX_FLOAT_EPS);

    --z__;

    if (*n0 - *i0 - 1 <= 0) {
	return;
    }

    j4 = (*i0 << 2) + *pp - 3;
    emin = z__[j4 + 4];
    d__ = z__[j4];
    *dmin__ = d__;

    if (*pp == 0) {
	i__1 = 4*(*n0 - 3);
	for (j4 = *i0*4; j4 <= i__1; j4 += 4) {
	    z__[j4 - 2] = d__ + z__[j4 - 1];
	    if (fabs(z__[j4 - 2])<GMX_FLOAT_MIN) {
		z__[j4] = 0.;
		d__ = z__[j4 + 1];
		*dmin__ = d__;
		emin = 0.;
	    } else if (safemin * z__[j4 + 1] < z__[j4 - 2] && safemin * z__[j4 
		    - 2] < z__[j4 + 1]) {
		temp = z__[j4 + 1] / z__[j4 - 2];
		z__[j4] = z__[j4 - 1] * temp;
		d__ *= temp;
	    } else {
		z__[j4] = z__[j4 + 1] * (z__[j4 - 1] / z__[j4 - 2]);
		d__ = z__[j4 + 1] * (d__ / z__[j4 - 2]);
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
	    if (fabs(z__[j4 - 3])<GMX_FLOAT_MIN) {
		z__[j4 - 1] = 0.;
		d__ = z__[j4 + 2];
		*dmin__ = d__;
		emin = 0.;
	    } else if (safemin * z__[j4 + 2] < z__[j4 - 3] && safemin * z__[j4 
		    - 3] < z__[j4 + 2]) {
		temp = z__[j4 + 2] / z__[j4 - 3];
		z__[j4 - 1] = z__[j4] * temp;
		d__ *= temp;
	    } else {
		z__[j4 - 1] = z__[j4 + 2] * (z__[j4] / z__[j4 - 3]);
		d__ = z__[j4 + 2] * (d__ / z__[j4 - 3]);
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
    if (fabs(z__[j4 - 2])<GMX_FLOAT_MIN) {
	z__[j4] = 0.;
	*dnm1 = z__[j4p2 + 2];
	*dmin__ = *dnm1;
	emin = 0.;
    } else if (safemin * z__[j4p2 + 2] < z__[j4 - 2] && safemin * z__[j4 - 2] < 
	    z__[j4p2 + 2]) {
	temp = z__[j4p2 + 2] / z__[j4 - 2];
	z__[j4] = z__[j4p2] * temp;
	*dnm1 = *dnm2 * temp;
    } else {
	z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
	*dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]);
    }
    if(*dnm1<*dmin__)
      *dmin__ = *dnm1;

    *dmin1 = *dmin__;
    j4 += 4;
    j4p2 = j4 + (*pp << 1) - 1;
    z__[j4 - 2] = *dnm1 + z__[j4p2];
    if (fabs(z__[j4 - 2])<GMX_FLOAT_MIN) {
	z__[j4] = 0.;
	*dn = z__[j4p2 + 2];
	*dmin__ = *dn;
	emin = 0.;
    } else if (safemin * z__[j4p2 + 2] < z__[j4 - 2] && safemin * z__[j4 - 2] < 
	    z__[j4p2 + 2]) {
	temp = z__[j4p2 + 2] / z__[j4 - 2];
	z__[j4] = z__[j4p2] * temp;
	*dn = *dnm1 * temp;
    } else {
	z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
	*dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]);
    }
    if(*dn<*dmin__)
      *dmin__ = *dn;

    z__[j4 + 2] = *dn;
    z__[(*n0 << 2) - *pp] = emin;
    return;


} 
