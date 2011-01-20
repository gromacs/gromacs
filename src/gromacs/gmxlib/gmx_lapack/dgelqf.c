#include <math.h>
#include "gmx_lapack.h"
#include "lapack_limits.h"



void
F77_FUNC(dgelqf,DGELQF)(int *m,
	int *n, 
	double *a, 
	int *lda, 
	double *tau,
	double *work, 
	int *lwork, 
	int *info)
{
    int a_dim1, a_offset, i__1, i__2, i__3, i__4;

    int i__, k, ib, nb, nx, iws, nbmin, iinfo;
    int ldwork, lwkopt;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    *info = 0;
    nb = DGELQF_BLOCKSIZE;
    lwkopt = *m * nb;
    work[1] = (double) lwkopt;

    if (*lwork==-1) {
	return;
    }

    k =(*m < *n) ? *m : *n;
    if (k == 0) {
	work[1] = 1.;
	return;
    }

    nbmin = 2;
    nx = 0;
    iws = *m;
    if (nb > 1 && nb < k) {
	nx = DGELQF_CROSSOVER;
	if (nx < k) {
	    ldwork = *m;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

		nb = *lwork / ldwork;
		nbmin = DGELQF_MINBLOCKSIZE;
	    }
	}
    }

    if (nb >= nbmin && nb < k && nx < k) {

	i__1 = k - nx;
	i__2 = nb;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	    i__3 = k - i__ + 1;
	    ib = (i__3 < nb) ? i__3 : nb;

	    i__3 = *n - i__ + 1;
	    F77_FUNC(dgelq2,DGELQ2)(&ib, &i__3, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[
		    1], &iinfo);
	    if (i__ + ib <= *m) {

		i__3 = *n - i__ + 1;
		F77_FUNC(dlarft,DLARFT)("Forward", "Rowwise", &i__3, &ib, &a[i__ + i__ * 
			a_dim1], lda, &tau[i__], &work[1], &ldwork);

		i__3 = *m - i__ - ib + 1;
		i__4 = *n - i__ + 1;
		F77_FUNC(dlarfb,DLARFB)("Right", "No transpose", "Forward", "Rowwise", &i__3, 
			&i__4, &ib, &a[i__ + i__ * a_dim1], lda, &work[1], &
			ldwork, &a[i__ + ib + i__ * a_dim1], lda, &work[ib + 
			1], &ldwork);
	    }
	}
    } else {
	i__ = 1;
    }

    if (i__ <= k) {
	i__2 = *m - i__ + 1;
	i__1 = *n - i__ + 1;
	F77_FUNC(dgelq2,DGELQ2)(&i__2, &i__1, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[1]
		, &iinfo);
    }

    work[1] = (double) iws;
    return;

}
