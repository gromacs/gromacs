#include "gmx_lapack.h"

#define SORGLQ_BLOCKSIZE    32
#define SORGLQ_MINBLOCKSIZE 2
#define SORGLQ_CROSSOVER    128


void 
F77_FUNC(sorglq,SORGLQ)(int *m, 
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

    *info = 0;
    ki = 0;
    nb = SORGLQ_BLOCKSIZE;
    lwkopt = (*m) * nb;
    work[1] = (float) lwkopt;
    lquery = *lwork == -1;
    if (*m < 0) {
	*info = -1;
    } else if (*n < *m) {
	*info = -2;
    } else if (*k < 0 || *k > *m) {
	*info = -3;
    } else if (*lda < (*m)) {
	*info = -5;
    } else if (*lwork < (*m) && ! lquery) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	return;
    } else if (lquery) {
	return;
    }

    if (*m <= 0) {
	work[1] = 1.;
	return;
    }

    nbmin = 2;
    nx = 0;
    iws = *m;
    if (nb > 1 && nb < *k) {

	nx = SORGLQ_CROSSOVER;
	if (nx < *k) {

	    ldwork = *m;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

		nb = *lwork / ldwork;
		nbmin = SORGLQ_MINBLOCKSIZE;
	    }
	}
    }

    if (nb >= nbmin && nb < *k && nx < *k) {

	ki = (*k - nx - 1) / nb * nb;
	i__1 = *k, i__2 = ki + nb;
	kk = (i__1<i__2) ? i__1 : i__2;

	i__1 = kk;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = kk + 1; i__ <= i__2; ++i__) {
		a[i__ + j * a_dim1] = 0.;
	    }
	}
    } else {
	kk = 0;
    }
    if (kk < *m) {
	i__1 = *m - kk;
	i__2 = *n - kk;
	i__3 = *k - kk;
	F77_FUNC(sorgl2,SORGL2)(&i__1, &i__2, &i__3, &a[kk + 1 + (kk + 1) * a_dim1], lda, &
		tau[kk + 1], &work[1], &iinfo);
    }

    if (kk > 0) {

	i__1 = -nb;
	for (i__ = ki + 1; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
	    i__2 = nb, i__3 = *k - i__ + 1;
	    ib = (i__2<i__3) ? i__2 : i__3;
	    if (i__ + ib <= *m) {

		i__2 = *n - i__ + 1;
		F77_FUNC(slarft,SLARFT)("Forward", "Rowwise", &i__2, &ib, &a[i__ + i__ * 
			a_dim1], lda, &tau[i__], &work[1], &ldwork);

		i__2 = *m - i__ - ib + 1;
		i__3 = *n - i__ + 1;
		F77_FUNC(slarfb,SLARFB)("Right", "Transpose", "Forward", "Rowwise", &i__2, &
			i__3, &ib, &a[i__ + i__ * a_dim1], lda, &work[1], &
			ldwork, &a[i__ + ib + i__ * a_dim1], lda, &work[ib + 
			1], &ldwork);
	    }

	    i__2 = *n - i__ + 1;
	    F77_FUNC(sorgl2,SORGL2)(&ib, &i__2, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &
		    work[1], &iinfo);

	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = i__ + ib - 1;
		for (l = i__; l <= i__3; ++l) {
		    a[l + j * a_dim1] = 0.;
		}
	    }
	}
    }

    work[1] = (float) iws;
    return;

}


