#include "gmx_blas.h"
#include "gmx_lapack.h"
#include "lapack_limits.h"

void
F77_FUNC(sgetri,SGETRI)(int *n, 
	float *a, 
	int *lda, 
	int *ipiv, 
	float *work, 
	int *lwork, 
	int *info)
{
    int a_dim1, a_offset, i__1, i__2, i__3;

    int i__, j, jb, nb, jj, jp, nn, iws;
    int nbmin;
    int ldwork;
    int lwkopt;
    int c__1 = 1;
    float c_b20 = -1.;
    float c_b22 = 1.;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    --work;

    *info = 0;
    nb = DGETRI_BLOCKSIZE;
    lwkopt = *n * nb;
    work[1] = (float) lwkopt;

    if (*n < 0) {
	*info = -1;
    } else if (*lda < (*n)) {
	*info = -3;
    } else if (*lwork < (*n) && *lwork!=-1) {
	*info = -6;
    }
    if (*info != 0) {
	i__1 = -(*info);
	return;
    } else if (*lwork == -1) {
	return;
    }

    if (*n == 0) {
	return;
    }

    F77_FUNC(strtri,STRTRI)("Upper", "Non-unit", n, &a[a_offset], lda, info);
    if (*info > 0) {
	return;
    }

    nbmin = 2;
    ldwork = *n;
    if (nb > 1 && nb < *n) {
	i__1 = ldwork * nb;
	iws = (i__1>1) ? i__1 : 1;
	if (*lwork < iws) {
	    nb = *lwork / ldwork;
	    nbmin = DGETRI_MINBLOCKSIZE;
	}
    } else {
	iws = *n;
    }

    if (nb < nbmin || nb >= *n) {

	for (j = *n; j >= 1; --j) {

	    i__1 = *n;
	    for (i__ = j + 1; i__ <= i__1; ++i__) {
		work[i__] = a[i__ + j * a_dim1];
		a[i__ + j * a_dim1] = 0.;
	    }

	    if (j < *n) {
		i__1 = *n - j;
		F77_FUNC(sgemv,SGEMV)("No transpose", n, &i__1, &c_b20, &a[(j + 1) * a_dim1 
			+ 1], lda, &work[j + 1], &c__1, &c_b22, &a[j * a_dim1 
			+ 1], &c__1);
	    }
	}
    } else {

	nn = (*n - 1) / nb * nb + 1;
	i__1 = -nb;
	for (j = nn; i__1 < 0 ? j >= 1 : j <= 1; j += i__1) {
	    i__2 = nb, i__3 = *n - j + 1;
	    jb = (i__2<i__3) ? i__2 : i__3;

	    i__2 = j + jb - 1;
	    for (jj = j; jj <= i__2; ++jj) {
		i__3 = *n;
		for (i__ = jj + 1; i__ <= i__3; ++i__) {
		    work[i__ + (jj - j) * ldwork] = a[i__ + jj * a_dim1];
		    a[i__ + jj * a_dim1] = 0.;
		}
	    }

	    if (j + jb <= *n) {
		i__2 = *n - j - jb + 1;
		F77_FUNC(sgemm,SGEMM)("No transpose", "No transpose", n, &jb, &i__2, &c_b20, 
			&a[(j + jb) * a_dim1 + 1], lda, &work[j + jb], &
			ldwork, &c_b22, &a[j * a_dim1 + 1], lda);
	    }
	    F77_FUNC(strsm,STRSM)("Right", "Lower", "No transpose", "Unit", n, &jb, &c_b22, &
		    work[j], &ldwork, &a[j * a_dim1 + 1], lda);
	}
    }

    for (j = *n - 1; j >= 1; --j) {
	jp = ipiv[j];
	if (jp != j) {
	    F77_FUNC(sswap,SSWAP)(n, &a[j * a_dim1 + 1], &c__1, &a[jp * a_dim1 + 1], &c__1);
	}
    }

    work[1] = (float) iws;
    return;

}


