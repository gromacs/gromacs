#include <cmath>

#include "gromacs/utility/real.h"

#include "../gmx_blas.h"
#include "../gmx_lapack.h"
#include "lapack_limits.h"


void
F77_FUNC(slarrfx,SLARRFX)(int *n, 
	float *d__, 
	float *l, 
	float *ld, 
	float *lld, 
	int *ifirst, 
	int *ilast, 
	float *w, 
	float *sigma, 
	float *dplus, 
	float *lplus, 
	float *work,
	int *info)
{
    int i1 = 1;
    int i__1;
    float d__2, d__3;

    int i__;
    float s, eps, tmp, dmax1, dmax2, delta;
    --work;
    --lplus;
    --dplus;
    --w;
    --lld;
    --ld;
    --l;
    --d__;
    *info = 0;
    eps = GMX_FLOAT_EPS;
    *sigma = w[*ifirst];
    delta = eps * 2.;

L10:
    s = -(*sigma);
    dplus[1] = d__[1] + s;
    dmax1 = std::abs(dplus[1]);
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	lplus[i__] = ld[i__] / dplus[i__];
	s = s * lplus[i__] * l[i__] - *sigma;
	dplus[i__ + 1] = d__[i__ + 1] + s;
	d__2 = dmax1, d__3 = std::abs(dplus[i__ + 1]);
	dmax1 = (d__2>d__3) ? d__2 : d__3;
    }
    if (std::isnan(dmax1)) {
	*sigma -= std::abs(*sigma) * delta;
	delta *= 2.;
	goto L10;
    }

    tmp = w[*ilast];
    delta = eps * 2.;
L30:
    s = -tmp;
    work[1] = d__[1] + s;
    dmax2 = std::abs(work[1]);
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	work[*n + i__] = ld[i__] / work[i__];
	s = s * work[*n + i__] * l[i__] - tmp;
	work[i__ + 1] = d__[i__ + 1] + s;
	d__2 = dmax2, d__3 = std::abs(work[i__ + 1]);
	dmax2 = (d__2>d__3) ? d__2 : d__3;
    }
    if (std::isnan(dmax2)) {
	tmp += std::abs(tmp) * delta;
	delta *= 2.;
	goto L30;
    }
    if (dmax2 < dmax1) {
	*sigma = tmp;
	F77_FUNC(scopy,SCOPY)(n, &work[1], &i1, &dplus[1], &i1);
	i__1 = *n - 1;
	F77_FUNC(scopy,SCOPY)(&i__1, &work[*n + 1], &i1, &lplus[1], &i1);
    }

    return;
}
