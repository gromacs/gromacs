#include <math.h>
#include "gmx_lapack.h"


void
F77_FUNC(slae2,SLAE2)(float *a, 
       float *b,
       float *c__, 
       float *rt1, 
       float *rt2)
{
    float d__1;
    float ab, df, tb, sm, rt, adf, acmn, acmx;


    sm = *a + *c__;
    df = *a - *c__;
    adf = fabs(df);
    tb = *b + *b;
    ab = fabs(tb);
    if (fabs(*a) > fabs(*c__)) {
	acmx = *a;
	acmn = *c__;
    } else {
	acmx = *c__;
	acmn = *a;
    }
    if (adf > ab) {
	d__1 = ab / adf;
	rt = adf * sqrt(d__1 * d__1 + 1.);
    } else if (adf < ab) {
	d__1 = adf / ab;
	rt = ab * sqrt(d__1 * d__1 + 1.);
    } else {

	rt = ab * sqrt(2.);
    }
    if (sm < 0.) {
	*rt1 = (sm - rt) * .5;
	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else if (sm > 0.) {
	*rt1 = (sm + rt) * .5;
	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else {
	*rt1 = rt * .5;
	*rt2 = rt * -.5;
    }
    return;

}


