#include "gmx_lapack.h"

void
F77_FUNC(dorm2l,DORM2L)(char *side, 
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
	int *info)
{
    int a_dim1, a_offset, c_dim1, c_offset, i__1, i__2;
    int c__1 = 1;

    int i__, i1, i2, i3, mi, ni, nq;
    double aii;
    int left;
    int notran;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;

    /* Function Body */
    *info = 0;
    left = (*side=='L' || *side=='l');
    notran = (*trans=='N' || *trans=='n');

    if (left) {
	nq = *m;
    } else {
	nq = *n;
    }
    if (*info != 0) {
	i__1 = -(*info);
	return;
    }

    if (*m == 0 || *n == 0 || *k == 0) {
	return;
    }

    if ((left && notran) || (! left && ! notran)) {
	i1 = 1;
	i2 = *k;
	i3 = 1;
    } else {
	i1 = *k;
	i2 = 1;
	i3 = -1;
    }

    if (left) {
	ni = *n;
    } else {
	mi = *m;
    }

    i__1 = i2;
    i__2 = i3;
    for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	if (left) {

	    mi = *m - *k + i__;
	} else {

	    ni = *n - *k + i__;
	}

	aii = a[nq - *k + i__ + i__ * a_dim1];
	a[nq - *k + i__ + i__ * a_dim1] = 1.;
	F77_FUNC(dlarf,DLARF)(side, &mi, &ni, &a[i__ * a_dim1 + 1], &c__1, &tau[i__], &c__[
		c_offset], ldc, &work[1]);
	a[nq - *k + i__ + i__ * a_dim1] = aii;
    }
    return;
}
