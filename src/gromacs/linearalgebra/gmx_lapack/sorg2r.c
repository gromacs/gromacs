#include <math.h>

#include "../gmx_blas.h"
#include "../gmx_lapack.h"

void 
F77_FUNC(sorg2r,SORG2R)(int *m, 
                        int *n,
                        int *k, 
                        float *a, 
                        int *lda,
                        float *tau,
                        float *work,
                        int *info)
{
    int a_dim1, a_offset, i__1, i__2;
    float r__1;
    int c__1 = 1;

    int i__, j, l;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    *info = 0;

    if (*n <= 0) {
        return;
    }

    i__1 = *n;
    for (j = *k + 1; j <= i__1; ++j) {
	i__2 = *m;
	for (l = 1; l <= i__2; ++l) {
	    a[l + j * a_dim1] = 0.0;
	}
	a[j + j * a_dim1] = 1.0;
    }
    for (i__ = *k; i__ >= 1; --i__) {
	if (i__ < *n) {
	    a[i__ + i__ * a_dim1] = 1.0;
	    i__1 = *m - i__ + 1;
	    i__2 = *n - i__;
	    F77_FUNC(slarf,SLARF)("L", &i__1, &i__2, &a[i__ + i__ * a_dim1], &c__1, 
                              &tau[i__], &a[i__ + (i__ + 1) * a_dim1], lda, &work[1]);
	}
	if (i__ < *m) {
	    i__1 = *m - i__;
	    r__1 = -tau[i__];
	    F77_FUNC(sscal,SSCAL)(&i__1, &r__1, &a[i__ + 1 + i__ * a_dim1], &c__1);
	}
	a[i__ + i__ * a_dim1] = 1.0 - tau[i__];
	i__1 = i__ - 1;
	for (l = 1; l <= i__1; ++l) {
	    a[l + i__ * a_dim1] = 0.0;
	}
    }
    return;

}


