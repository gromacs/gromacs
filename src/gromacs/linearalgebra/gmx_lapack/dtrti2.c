#include "../gmx_blas.h"
#include "../gmx_lapack.h"
#include "lapack_limits.h"

void
F77_FUNC(dtrti2,DTRTI2)(const char *uplo,
	const char *diag, 
	int *n, 
	double *a,
	int *lda,
	int *info)
{
    int a_dim1, a_offset, i__1, i__2;

    int j;
    double ajj;
    int upper, nounit;
    int c__1 = 1;


    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    *info = 0;
    upper = (*uplo=='U' || *uplo=='u');
    nounit = (*diag=='N' || *diag=='n');

    if (*info != 0) {
	i__1 = -(*info);
	return;
    }

    if (upper) {

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (nounit) {
		a[j + j * a_dim1] = 1. / a[j + j * a_dim1];
		ajj = -a[j + j * a_dim1];
	    } else {
		ajj = -1.;
	    }

	    i__2 = j - 1;
	    F77_FUNC(dtrmv,DTRMV)("Upper", "No transpose", diag, &i__2, &a[a_offset], lda, &
		    a[j * a_dim1 + 1], &c__1);
	    i__2 = j - 1;
	    F77_FUNC(dscal,DSCAL)(&i__2, &ajj, &a[j * a_dim1 + 1], &c__1);
	}
    } else {

	for (j = *n; j >= 1; --j) {
	    if (nounit) {
		a[j + j * a_dim1] = 1. / a[j + j * a_dim1];
		ajj = -a[j + j * a_dim1];
	    } else {
		ajj = -1.;
	    }
	    if (j < *n) {

		i__1 = *n - j;
		F77_FUNC(dtrmv,DTRMV)("Lower", "No transpose", diag, &i__1, &a[j + 1 + (j + 
			1) * a_dim1], lda, &a[j + 1 + j * a_dim1], &c__1);
		i__1 = *n - j;
		F77_FUNC(dscal,DSCAL)(&i__1, &ajj, &a[j + 1 + j * a_dim1], &c__1);
	    }
	}
    }
    return;
}
