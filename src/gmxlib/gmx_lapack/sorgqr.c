#include "gmx_lapack.h"
#include "lapack_limits.h"


void 
F77_FUNC(sorgqr,SORGQR)(int *m, 
	int *n, 
	int *k, 
	float *a, 
	int *lda, 
	float *tau, 
	float *work, 
	int *lwork, 
	int *info)
{
    int a_dim1, a_offset, i__1, i__2, i__3;

    int i__, j, l, ib, nb, ki, kk, nx, iws, nbmin, iinfo;
    int ldwork, lwkopt;
    int lquery;
 
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    ki = 0;
    *info = 0;
    nb = DORGQR_BLOCKSIZE;
    lwkopt = (*n) * nb;
    work[1] = (float) lwkopt;
    lquery = *lwork == -1;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0 || *n > *m) {
	*info = -2;
    } else if (*k < 0 || *k > *n) {
	*info = -3;
    } else if (*lda < (*m)) {
	*info = -5;
    } else if (*lwork < (*n) && ! lquery) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	return;
    } else if (lquery) {
	return;
    }

    if (*n <= 0) {
	work[1] = 1.;
	return;
    }

    nbmin = 2;
    nx = 0;
    iws = *n;
    if (nb > 1 && nb < *k) {

	nx = DORGQR_CROSSOVER;
	if (nx < *k) {

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

		nb = *lwork / ldwork;
		nbmin = DORGQR_MINBLOCKSIZE;
	    }
	}
    }

    if (nb >= nbmin && nb < *k && nx < *k) {

	ki = (*k - nx - 1) / nb * nb;
	i__1 = *k, i__2 = ki + nb;
	kk = (i__1<i__2) ? i__1 : i__2;

	i__1 = *n;
	for (j = kk + 1; j <= i__1; ++j) {
	    i__2 = kk;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		a[i__ + j * a_dim1] = 0.;
	    }
	}
    } else {
	kk = 0;
    }

    if (kk < *n) {
	i__1 = *m - kk;
	i__2 = *n - kk;
	i__3 = *k - kk;
	F77_FUNC(sorg2r,SORG2R)(&i__1, &i__2, &i__3, &a[kk + 1 + (kk + 1) * a_dim1], lda, &
		tau[kk + 1], &work[1], &iinfo);
    }

    if (kk > 0) {

	i__1 = -nb;
	for (i__ = ki + 1; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
	    i__2 = nb, i__3 = *k - i__ + 1;
	    ib = (i__2<i__3) ? i__2 : i__3;
	    if (i__ + ib <= *n) {

		i__2 = *m - i__ + 1;
		F77_FUNC(slarft,SLARFT)("Forward", "Columnwise", &i__2, &ib, &a[i__ + i__ * 
			a_dim1], lda, &tau[i__], &work[1], &ldwork);

		i__2 = *m - i__ + 1;
		i__3 = *n - i__ - ib + 1;
		F77_FUNC(slarfb,SLARFB)("Left", "No transpose", "Forward", "Columnwise", &
			i__2, &i__3, &ib, &a[i__ + i__ * a_dim1], lda, &work[
			1], &ldwork, &a[i__ + (i__ + ib) * a_dim1], lda, &
			work[ib + 1], &ldwork);
	    }

	    i__2 = *m - i__ + 1;
	    F77_FUNC(sorg2r,SORG2R)(&i__2, &ib, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &
		    work[1], &iinfo);

	    i__2 = i__ + ib - 1;
	    for (j = i__; j <= i__2; ++j) {
		i__3 = i__ - 1;
		for (l = 1; l <= i__3; ++l) {
		    a[l + j * a_dim1] = 0.;
		}
	    }
	}
    }

    work[1] = (float) iws;
    return;

} 
