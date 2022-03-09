#include <cmath>

#include "../gmx_lapack.h"


void F77_FUNC(slae2, SLAE2)(float* a, float* b, float* c__, float* rt1, float* rt2)
{
    float d__1;
    float ab, df, tb, sm, rt, adf, acmn, acmx;


    sm  = *a + *c__;
    df  = *a - *c__;
    adf = std::abs(df);
    tb  = *b + *b;
    ab  = std::abs(tb);
    if (std::abs(*a) > std::abs(*c__))
    {
        acmx = *a;
        acmn = *c__;
    }
    else
    {
        acmx = *c__;
        acmn = *a;
    }
    if (adf > ab)
    {
        d__1 = ab / adf;
        rt   = adf * std::sqrt(d__1 * d__1 + 1.);
    }
    else if (adf < ab)
    {
        d__1 = adf / ab;
        rt   = ab * std::sqrt(d__1 * d__1 + 1.);
    }
    else
    {

        rt = ab * std::sqrt(2.);
    }
    if (sm < 0.)
    {
        *rt1 = (sm - rt) * .5;
        *rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    }
    else if (sm > 0.)
    {
        *rt1 = (sm + rt) * .5;
        *rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    }
    else
    {
        *rt1 = rt * .5;
        *rt2 = rt * -.5;
    }
    return;
}
