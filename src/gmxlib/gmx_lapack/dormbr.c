#include "gmx_lapack.h"
#include "lapack_limits.h"

void 
F77_FUNC(dormbr,DORMBR)(char *vect, 
	char *side, 
	char *trans, 
	int *m, 
	int *n, 
	int *k, 
	double *a, 
	int *lda, 
	double *tau, 
	double *c__, 
	int *ldc, 
	double *work, 
	int *lwork, 
	int *info)
{
    int a_dim1, a_offset, c_dim1, c_offset, i__1;
 

    int i1, i2, nb, mi, ni, nq, nw;
    int left;
    int iinfo;
    int notran;
    int applyq;
    char transt[1];
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
    applyq = (*vect=='Q' || *vect=='q');
    left = (*side=='L' || *side=='l');
    notran = (*trans=='N' || *trans=='n');
    lquery = *lwork == -1;

    if (left) {
	nq = *m;
	nw = *n;
    } else {
	nq = *n;
	nw = *m;
    }

    nb = DORMQR_BLOCKSIZE;
    lwkopt = nw * nb;
    work[1] = (double) lwkopt;
    
    if (*info != 0) {
	i__1 = -(*info);
	return;
    } else if (lquery) {
	return;
    }

    work[1] = 1.;
    if (*m == 0 || *n == 0) {
	return;
    }

    if (applyq) {

	if (nq >= *k) {

	    F77_FUNC(dormqr,DORMQR)(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		    c_offset], ldc, &work[1], lwork, &iinfo);
	} else if (nq > 1) {

	    if (left) {
		mi = *m - 1;
		ni = *n;
		i1 = 2;
		i2 = 1;
	    } else {
		mi = *m;
		ni = *n - 1;
		i1 = 1;
		i2 = 2;
	    }
	    i__1 = nq - 1;
	    F77_FUNC(dormqr,DORMQR)(side, trans, &mi, &ni, &i__1, &a[a_dim1 + 2], lda, &tau[1]
		    , &c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, &iinfo);
	}
    } else {

	if (notran) {
	    *(unsigned char *)transt = 'T';
	} else {
	    *(unsigned char *)transt = 'N';
	}
	if (nq > *k) {

	    F77_FUNC(dormlq,DORMLQ)(side, transt, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		    c_offset], ldc, &work[1], lwork, &iinfo);
	} else if (nq > 1) {

	    if (left) {
		mi = *m - 1;
		ni = *n;
		i1 = 2;
		i2 = 1;
	    } else {
		mi = *m;
		ni = *n - 1;
		i1 = 1;
		i2 = 2;
	    }
	    i__1 = nq - 1;
	    F77_FUNC(dormlq,DORMLQ)(side, transt, &mi, &ni, &i__1, &a[(a_dim1 << 1) + 1], lda,
		     &tau[1], &c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, &
		    iinfo);
	}
    }
    work[1] = (double) lwkopt;
    return;


}


