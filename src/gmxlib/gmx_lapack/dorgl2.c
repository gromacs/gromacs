#include "gmx_blas.h"
#include "gmx_lapack.h"

void
F77_FUNC(dorgl2,DORGL2)(int *m,
                        int *n, 
                        int *k, 
                        double *a, 
                        int *lda, 
                        double *tau, 
                        double *work, 
                        int *info)
{
    int a_dim1, a_offset, i__1, i__2;
    double r__1;

    int i__, j, l;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    i__ = (*m > 1) ? *m : 1;
    
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < *m) {
	*info = -2;
    } else if (*k < 0 || *k > *m) {
	*info = -3;
    } else if (*lda < i__) {
	*info = -5;
    }
    if (*info != 0) {
	return;
    }
    if (*m <= 0) {
	return;
    }

    if (*k < *m) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (l = *k + 1; l <= i__2; ++l) {
		a[l + j * a_dim1] = 0.0;
	    }
	    if (j > *k && j <= *m) {
		a[j + j * a_dim1] = 1.0;
	    }
	}
    }

    for (i__ = *k; i__ >= 1; --i__) {
	if (i__ < *n) {
	    if (i__ < *m) {
		a[i__ + i__ * a_dim1] = 1.0;
		i__1 = *m - i__;
		i__2 = *n - i__ + 1;
		F77_FUNC(dlarf,DLARF)("R", &i__1, &i__2, &a[i__ + i__ * a_dim1], lda, 
               &tau[i__], &a[i__ + 1 + i__ * a_dim1], lda, &work[1]);
	    }
	    i__1 = *n - i__;
	    r__1 = -tau[i__];
	    F77_FUNC(dscal,DSCAL)(&i__1, &r__1, &a[i__ + (i__ + 1) * a_dim1], lda);
	}
	a[i__ + i__ * a_dim1] = 1.0 - tau[i__];
	i__1 = i__ - 1;
	for (l = 1; l <= i__1; ++l) {
	    a[i__ + l * a_dim1] = 0.0;
	}
    }
    return;

}



