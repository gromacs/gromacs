#include <cmath>

#include "../gmx_lapack.h"

#include "gromacs/utility/real.h"

void F77_FUNC(dlas2, DLAS2)(double* f, double* g, double* h, double* ssmin, double* ssmax)
{
    double fa = std::abs(*f);
    double ga = std::abs(*g);
    double ha = std::abs(*h);
    double fhmin, fhmax, tmax, tmin, tmp1, tmp2;
    double as, at, au, c;

    fhmin = (fa < ha) ? fa : ha;
    fhmax = (fa > ha) ? fa : ha;

    if (std::abs(fhmin) < GMX_DOUBLE_MIN)
    {
        *ssmin = 0.0;
        if (std::abs(fhmax) < GMX_DOUBLE_MIN)
            *ssmax = ga;
        else
        {
            tmax   = (fhmax > ga) ? fhmax : ga;
            tmin   = (fhmax < ga) ? fhmax : ga;
            tmp1   = tmin / tmax;
            tmp1   = tmp1 * tmp1;
            *ssmax = tmax * std::sqrt(1.0 + tmp1);
        }
    }
    else
    {
        if (ga < fhmax)
        {
            as     = 1.0 + fhmin / fhmax;
            at     = (fhmax - fhmin) / fhmax;
            au     = (ga / fhmax);
            au     = au * au;
            c      = 2.0 / (std::sqrt(as * as + au) + std::sqrt(at * at + au));
            *ssmin = fhmin * c;
            *ssmax = fhmax / c;
        }
        else
        {
            au = fhmax / ga;
            if (std::abs(au) < GMX_DOUBLE_MIN)
            {
                *ssmin = (fhmin * fhmax) / ga;
                *ssmax = ga;
            }
            else
            {
                as     = 1.0 + fhmin / fhmax;
                at     = (fhmax - fhmin) / fhmax;
                tmp1   = as * au;
                tmp2   = at * au;
                c      = 1.0 / (std::sqrt(1.0 + tmp1 * tmp1) + std::sqrt(1.0 + tmp2 * tmp2));
                *ssmin = (fhmin * c) * au;
                *ssmin = *ssmin + *ssmin;
                *ssmax = ga / (c + c);
            }
        }
    }
    return;
}
