#include "../gmx_lapack.h"

void 
F77_FUNC(sorm2r,SORM2R)(const char *side, 
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
	int *info)
{
    int a_dim1, a_offset, c_dim1, c_offset, i__1, i__2;

    int i__, i1, i2, i3, ic, jc, mi, ni, nq;
    float aii;
    int left;
    int notran;
    int c__1 = 1;

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

    ic = jc = 0;

    if (left) {
	nq = *m;
    } else {
	nq = *n;
    }

    if (*m <= 0 || *n <= 0 || *k <= 0) {
	return;
    }

    if ((left && !notran) || (!left && notran)) {
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
	jc = 1;
    } else {
	mi = *m;
	ic = 1;
    }

    i__1 = i2;
    i__2 = i3;
    for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	if (left) {

	    mi = *m - i__ + 1;
	    ic = i__;
	} else {

	    ni = *n - i__ + 1;
	    jc = i__;
	}


	aii = a[i__ + i__ * a_dim1];
	a[i__ + i__ * a_dim1] = 1.;
	F77_FUNC(slarf,SLARF)(side, &mi, &ni, &a[i__ + i__ * a_dim1], &c__1, &tau[i__], &c__[
		ic + jc * c_dim1], ldc, &work[1]);
	a[i__ + i__ * a_dim1] = aii;
    }
    return;

} 
