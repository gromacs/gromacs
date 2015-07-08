#include "../gmx_lapack.h"

void 
F77_FUNC(dgelq2,DGELQ2)(int *m, 
                        int *n, 
                        double *a,
                        int *lda, 
                        double *tau, 
                        double *work, 
                        int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    int i__, k;
    double aii;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    *info = 0;
    
    i__4 = (*m > 1) ? *m : 1;
    
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < i__4) {
	*info = -4;
    }
    if (*info != 0) {
	return;
    }

    
    k = (*m < *n ) ? *m : *n;
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n - i__ + 1;
	i__3 = i__ + 1;
    i__4 = (i__3 < *n) ? i__3 : *n;
	F77_FUNC(dlarfg,DLARFG)(&i__2, &a[i__ + i__ * a_dim1], &a[i__ + i__4 * a_dim1],
                            lda, &tau[i__]);
	if (i__ < *m) {
	    aii = a[i__ + i__ * a_dim1];
	    a[i__ + i__ * a_dim1] = 1.f;
	    i__2 = *m - i__;
	    i__3 = *n - i__ + 1;
	    F77_FUNC(dlarf,DLARF)("R", &i__2, &i__3, &a[i__ + i__ * a_dim1], lda, 
                              &tau[i__], &a[i__ + 1 + i__ * a_dim1], lda, &work[1]);
	    a[i__ + i__ * a_dim1] = aii;
	}
    }
    return;
}


