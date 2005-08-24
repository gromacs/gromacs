#include <math.h>
#include "gmx_blas.h"
#include "gmx_lapack.h"
#include "lapack_limits.h"

#include <types/simple.h>

void
F77_FUNC(dtrtri,DTRTRI)(char *uplo,
	char *diag, 
	int *n,
	double *a, 
	int *lda,
	int *info)
{
    int a_dim1, a_offset, i__1, i__3, i__4, i__5;
    int j, jb, nb, nn;
    double c_b18 = 1.;
    double c_b22 = -1.;

    int upper;
    int nounit;

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

    if (*n == 0) {
	return;
    }

    if (nounit) {
	i__1 = *n;
	for (*info = 1; *info <= i__1; ++(*info)) {
	    if (fabs(a[*info + *info * a_dim1])<GMX_DOUBLE_MIN) {
		return;
	    }
	}
	*info = 0;
    }

    nb = DTRTRI_BLOCKSIZE;
    if (nb <= 1 || nb >= *n) {

	F77_FUNC(dtrti2,DTRTI2)(uplo, diag, n, &a[a_offset], lda, info);
    } else {

	if (upper) {

	    i__1 = *n;
	    i__3 = nb;
	    for (j = 1; i__3 < 0 ? j >= i__1 : j <= i__1; j += i__3) {
		i__4 = nb, i__5 = *n - j + 1;
		jb = (i__4<i__5) ? i__4 : i__5;

		i__4 = j - 1;
		F77_FUNC(dtrmm,DTRMM)("Left", "Upper", "No transpose", diag, &i__4, &jb, &
			c_b18, &a[a_offset], lda, &a[j * a_dim1 + 1], lda);
		i__4 = j - 1;
		F77_FUNC(dtrsm,DTRSM)("Right", "Upper", "No transpose", diag, &i__4, &jb, &
			c_b22, &a[j + j * a_dim1], lda, &a[j * a_dim1 + 1], 
			lda);

		F77_FUNC(dtrti2,DTRTI2)("Upper", diag, &jb, &a[j + j * a_dim1], lda, info);
	    }
	} else {

	    nn = (*n - 1) / nb * nb + 1;
	    i__3 = -nb;
	    for (j = nn; i__3 < 0 ? j >= 1 : j <= 1; j += i__3) {
		i__1 = nb, i__4 = *n - j + 1;
		jb = (i__1<i__4) ? i__1 : i__4;
		if (j + jb <= *n) {

		    i__1 = *n - j - jb + 1;
		    F77_FUNC(dtrmm,DTRMM)("Left", "Lower", "No transpose", diag, &i__1, &jb, 
			    &c_b18, &a[j + jb + (j + jb) * a_dim1], lda, &a[j 
			    + jb + j * a_dim1], lda);
		    i__1 = *n - j - jb + 1;
		    F77_FUNC(dtrsm,DTRSM)("Right", "Lower", "No transpose", diag, &i__1, &jb,
			     &c_b22, &a[j + j * a_dim1], lda, &a[j + jb + j * 
			    a_dim1], lda);
		}

		F77_FUNC(dtrti2,DTRTI2)("Lower", diag, &jb, &a[j + j * a_dim1], lda, info);
	    }
	}
    }
    return;
}


