#include <cmath>
#include "gromacs/utility/real.h"

#include "../gmx_lapack.h"


void
F77_FUNC(slaev2,SLAEV2)(float *   a, 
	float *   b, 
	float *   c__, 
	float *   rt1, 
	float *   rt2, 
	float *   cs1, 
	float *   sn1)
{
    float d__1;

    float ab, df, cs, ct, tb, sm, tn, rt, adf, acs;
    int sgn1, sgn2;
    float acmn, acmx;

    sm = *a + *c__;
    df = *a - *c__;
    adf = std::abs(df);
    tb = *b + *b;
    ab = std::abs(tb);
    if (std::abs(*a) > std::abs(*c__)) {
	acmx = *a;
	acmn = *c__;
    } else {
	acmx = *c__;
	acmn = *a;
    }
    if (adf > ab) {
	d__1 = ab / adf;
	rt = adf *  std::sqrt(d__1 * d__1 + 1.);
    } else if (adf < ab) {
	d__1 = adf / ab;
	rt = ab *  std::sqrt(d__1 * d__1 + 1.);
    } else {

	rt = ab *  std::sqrt(2.);
    }
    if (sm < 0.) {
	*rt1 = (sm - rt) * .5;
	sgn1 = -1;

	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else if (sm > 0.) {
	*rt1 = (sm + rt) * .5;
	sgn1 = 1;
	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else {
	*rt1 = rt * .5;
	*rt2 = rt * -.5;
	sgn1 = 1;
    }
    if (df >= 0.) {
	cs = df + rt;
	sgn2 = 1;
    } else {
	cs = df - rt;
	sgn2 = -1;
    }
    acs = std::abs(cs);
    if (acs > ab) {
	ct = -tb / cs;
	*sn1 = 1. /  std::sqrt(ct * ct + 1.);
	*cs1 = ct * *sn1;
    } else {
	if (std::abs(ab)<GMX_FLOAT_MIN) {
	    *cs1 = 1.;
	    *sn1 = 0.;
	} else {
	    tn = -cs / tb;
	    *cs1 = 1. /  std::sqrt(tn * tn + 1.);
	    *sn1 = tn * *cs1;
	}
    }
    if (sgn1 == sgn2) {
	tn = *cs1;
	*cs1 = -(*sn1);
	*sn1 = tn;
    }
    return;

}


