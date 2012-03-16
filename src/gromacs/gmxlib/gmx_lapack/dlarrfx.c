#include <math.h>

#include <types/simple.h>

#include "gmx_blas.h"
#include "gmx_lapack.h"
#include "lapack_limits.h"


void
F77_FUNC(dlarrfx,DLARRFX)(int *n, 
	double *d__, 
	double *l, 
	double *ld, 
	double *lld, 
	int *ifirst, 
	int *ilast, 
	double *w, 
	double *sigma, 
	double *dplus, 
	double *lplus, 
	double *work,
	int *info)
{
    int i1 = 1;
    int i__1;
    double d__2, d__3;

    int i__;
    double s, eps, tmp, dmax1, dmax2, delta;
    --work;
    --lplus;
    --dplus;
    --w;
    --lld;
    --ld;
    --l;
    --d__;
    *info = 0;
    eps = GMX_DOUBLE_EPS;
    *sigma = w[*ifirst];
    delta = eps * 2.;

L10:
    s = -(*sigma);
    dplus[1] = d__[1] + s;
    dmax1 = fabs(dplus[1]);
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	lplus[i__] = ld[i__] / dplus[i__];
	s = s * lplus[i__] * l[i__] - *sigma;
	dplus[i__ + 1] = d__[i__ + 1] + s;
	d__2 = dmax1, d__3 = fabs(dplus[i__ + 1]);
	dmax1 = (d__2>d__3) ? d__2 : d__3;
    }
    if (! (dmax1 > 0. || dmax1 < 1.)) {
	*sigma -= fabs(*sigma) * delta;
	delta *= 2.;
	goto L10;
    }

    tmp = w[*ilast];
    delta = eps * 2.;
L30:
    s = -tmp;
    work[1] = d__[1] + s;
    dmax2 = fabs(work[1]);
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	work[*n + i__] = ld[i__] / work[i__];
	s = s * work[*n + i__] * l[i__] - tmp;
	work[i__ + 1] = d__[i__ + 1] + s;
	d__2 = dmax2, d__3 = fabs(work[i__ + 1]);
	dmax2 = (d__2>d__3) ? d__2 : d__3;
    }
    if (! (dmax2 > 0. || dmax2 < 1.)) {
	tmp += fabs(tmp) * delta;
	delta *= 2.;
	goto L30;
    }
    if (dmax2 < dmax1) {
	*sigma = tmp;
	F77_FUNC(dcopy,DCOPY)(n, &work[1], &i1, &dplus[1], &i1);
	i__1 = *n - 1;
	F77_FUNC(dcopy,DCOPY)(&i__1, &work[*n + 1], &i1, &lplus[1], &i1);
    }

    return;
}
