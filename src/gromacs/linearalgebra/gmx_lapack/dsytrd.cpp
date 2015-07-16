#include "../gmx_blas.h"
#include "../gmx_lapack.h"
#include "lapack_limits.h"

void
F77_FUNC(dsytrd,DSYTRD)(const char *uplo, int *n, double *a, int *
	lda, double *d__, double *e, double *tau, double *
	work, int *lwork, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    int i__, j, nb, kk, nx, iws;
    int nbmin, iinfo;
    int upper;
    int ldwork, lwkopt;
    int lquery;
    double c_b22 = -1.;
    double c_b23 = 1.;


    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    upper = (*uplo=='U' || *uplo=='u');
    lquery = (*lwork == -1);

    if (! upper && ! (*uplo=='L' || *uplo=='l')) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < ((1>*n) ? 1 : *n)) {
	*info = -4;
    } else if (*lwork < 1 && ! lquery) {
	*info = -9;
    }

    if (*info == 0) {

      nb = DSYTRD_BLOCKSIZE;
      lwkopt = *n * nb;
      work[1] = (double) lwkopt;
    } else
      return;

    if (lquery) 
      return;
  
    if (*n == 0) {
	work[1] = 1.;
	return;
    }

    nx = *n;
    if (nb > 1 && nb < *n) {

	nx = DSYTRD_CROSSOVER;
	if (nx < *n) {

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

		i__1 = *lwork / ldwork;
		nb = (i__1>1) ? i__1 : 1;
		nbmin = DSYTRD_MINBLOCKSIZE;
		if (nb < nbmin) {
		    nx = *n;
		}
	    }
	} else {
	    nx = *n;
	}
    } else {
	nb = 1;
    }

    if (upper) {

	kk = *n - (*n - nx + nb - 1) / nb * nb;
	i__1 = kk + 1;
	i__2 = -nb;
	for (i__ = *n - nb + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
		i__2) {

	    i__3 = i__ + nb - 1;
	    F77_FUNC(dlatrd,DLATRD)(uplo, &i__3, &nb, &a[a_offset], lda, &e[1], &tau[1], &
		    work[1], &ldwork);

	    i__3 = i__ - 1;
	    F77_FUNC(dsyr2k,DSYR2K)(uplo, "No transpose", &i__3, &nb, &c_b22, &a[i__ * a_dim1 
		    + 1], lda, &work[1], &ldwork, &c_b23, &a[a_offset], lda);

	    i__3 = i__ + nb - 1;
	    for (j = i__; j <= i__3; ++j) {
		a[j - 1 + j * a_dim1] = e[j - 1];
		d__[j] = a[j + j * a_dim1];

	    }

	}

	F77_FUNC(dsytd2,DSYTD2)(uplo, &kk, &a[a_offset], lda, &d__[1], &e[1], &tau[1], &iinfo);
    } else {

	i__2 = *n - nx;
	i__1 = nb;
	for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {


	    i__3 = *n - i__ + 1;
	    F77_FUNC(dlatrd,DLATRD)(uplo, &i__3, &nb, &a[i__ + i__ * a_dim1], lda, &e[i__], &
		    tau[i__], &work[1], &ldwork);

	    i__3 = *n - i__ - nb + 1;
	    F77_FUNC(dsyr2k,DSYR2K)(uplo, "No transpose", &i__3, &nb, &c_b22, &a[i__ + nb + 
		    i__ * a_dim1], lda, &work[nb + 1], &ldwork, &c_b23, &a[
		    i__ + nb + (i__ + nb) * a_dim1], lda);


	    i__3 = i__ + nb - 1;
	    for (j = i__; j <= i__3; ++j) {
		a[j + 1 + j * a_dim1] = e[j];
		d__[j] = a[j + j * a_dim1];

	    }

	}


	i__1 = *n - i__ + 1;
	F77_FUNC(dsytd2,DSYTD2)(uplo, &i__1, &a[i__ + i__ * a_dim1], lda, &d__[i__], &e[i__], 
		&tau[i__], &iinfo);
    }

    work[1] = (double) lwkopt;
    return;

}


