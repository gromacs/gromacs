#include "../gmx_lapack.h"
#include "lapack_limits.h"

void 
F77_FUNC(sormqr,SORMQR)(const char *side, 
	const char *trans, 
	int *m, 
	int *n, 
	int *k, 
	float *a, 
	int *lda, 
	float *tau, 
	float *c__, 
	int *ldc, 
	float *work, 
	int *lwork, 
	int *info)
{
   int a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__4, i__5;

    int i__;
    float t[4160];
    int i1, i2, i3, ib, ic, jc, nb, mi, ni, nq, nw, iws;
    int left;
    int nbmin, iinfo;
    int notran;
    int ldwork, lwkopt;
    int lquery;
    int ldt = 65;

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
    notran = (*trans=='N' || *trans=='n');
    lquery = *lwork == -1;

    if (left) {
	nq = *m;
	nw = *n;
    } else {
	nq = *n;
	nw = *m;
    }

     ic = jc = 0;
     nb = DORMQR_BLOCKSIZE;
     lwkopt = nw * nb;
     work[1] = (float) lwkopt;

    if (*info != 0) {
	i__1 = -(*info);
	return;
    } else if (lquery) {
      return;
    }

    if (*m == 0 || *n == 0 || *k == 0) {
	work[1] = 1.;
	return;
    }

    nbmin = 2;
    ldwork = nw;
    if (nb > 1 && nb < *k) {
	iws = nw * nb;
	if (*lwork < iws) {
	    nb = *lwork / ldwork;
	    nbmin = DORMQR_MINBLOCKSIZE;
	}
    } else {
	iws = nw;
    }

    if (nb < nbmin || nb >= *k) {

	F77_FUNC(sorm2r,SORM2R)(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo);
    } else {

	if ((left && !notran) || (!left && notran)) {
	    i1 = 1;
	    i2 = *k;
	    i3 = nb;
	} else {
	    i1 = (*k - 1) / nb * nb + 1;
	    i2 = 1;
	    i3 = -nb;
	}

	if (left) {
	    ni = *n;
	    jc = 1;
	} else {
	    mi = *m;
	    ic = 1;
	}

	i__1 = i2;
	i__2 = i3;
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	    i__4 = nb, i__5 = *k - i__ + 1;
	    ib = (i__4<i__5) ? i__4 : i__5;

	    i__4 = nq - i__ + 1;
	    F77_FUNC(slarft,SLARFT)("Forward", "Columnwise", &i__4, &ib, &a[i__ + i__ * 
		    a_dim1], lda, &tau[i__], t, &ldt);
	    if (left) {

		mi = *m - i__ + 1;
		ic = i__;
	    } else {
		ni = *n - i__ + 1;
		jc = i__;
	    }

	    F77_FUNC(slarfb,SLARFB)(side, trans, "Forward", "Columnwise", &mi, &ni, &ib, &a[
		    i__ + i__ * a_dim1], lda, t, &ldt, &c__[ic + jc * 
		    c_dim1], ldc, &work[1], &ldwork);
	}
    }
    work[1] = (float) lwkopt;
    return;


}


