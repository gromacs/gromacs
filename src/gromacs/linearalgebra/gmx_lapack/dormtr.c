#include "../gmx_lapack.h"
#include "lapack_limits.h"


void
F77_FUNC(dormtr,DORMTR)(const char *side, 
	const char *uplo,
	const char *trans, 
	int *m, 
	int *n,
	double *a, 
	int *lda, 
	double *tau, 
	double *c__, 
	int *ldc,
	double *work, 
	int *lwork, 
	int *info)
{
    int a_dim1, a_offset, c_dim1, c_offset, i__2;

    int i1, i2, nb, mi, ni, nq, nw;
    int left;
    int iinfo;
    int upper;
    int lwkopt;
    int lquery;


    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;

    *info = 0;
    left = (*side=='L' || *side=='l');
    upper = (*uplo=='U' || *uplo=='u');
    lquery = *lwork == -1;

    if (left) {
	nq = *m;
	nw = *n;
    } else {
	nq = *n;
	nw = *m;
    }


    nb = DORMQL_BLOCKSIZE;
    lwkopt = nw * nb;
    work[1] = (double) lwkopt;
    
    if (*info != 0) {
	i__2 = -(*info);
	return;
    } else if (lquery) {
	return;
    }

    if (*m == 0 || *n == 0 || nq == 1) {
	work[1] = 1.;
	return;
    }

    if (left) {
	mi = *m - 1;
	ni = *n;
    } else {
	mi = *m;
	ni = *n - 1;
    }

    if (upper) {
	i__2 = nq - 1;
	F77_FUNC(dormql,DORMQL)(side, trans, &mi, &ni, &i__2, &a[(a_dim1 << 1) + 1], lda, &
		tau[1], &c__[c_offset], ldc, &work[1], lwork, &iinfo);
    } else {
	if (left) {
	    i1 = 2;
	    i2 = 1;
	} else {
	    i1 = 1;
	    i2 = 2;
	}
	i__2 = nq - 1;
	F77_FUNC(dormqr,DORMQR)(side, trans, &mi, &ni, &i__2, &a[a_dim1 + 2], lda, &tau[1], &
		c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, &iinfo);
    }
    work[1] = (double) lwkopt;
    return;

}


