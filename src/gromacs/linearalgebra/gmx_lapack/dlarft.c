#include <math.h>
#include "gromacs/utility/real.h"

#include "../gmx_blas.h"
#include "../gmx_lapack.h"

void 
F77_FUNC(dlarft,DLARFT)(const char *direct, 
	const char *storev, 
	int *n, 
	int *k, 
	double *v, 
	int *ldv, 
	double *tau, 
	double *t, 
	int *ldt)
{
    /* System generated locals */
    int t_dim1, t_offset, v_dim1, v_offset, i__1, i__2, i__3;
    double d__1;

    /* Local variables */
    int i__, j;
    double vii;
    int c__1 = 1;
    double zero = 0.0;

    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --tau;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;

    if (*n == 0) {
	return;
    }

    if (*direct=='F' || *direct=='f') {
	i__1 = *k;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (fabs(tau[i__])<GMX_DOUBLE_MIN) {

		i__2 = i__;
		for (j = 1; j <= i__2; ++j) {
		    t[j + i__ * t_dim1] = 0.;
		}
	    } else {

		vii = v[i__ + i__ * v_dim1];
		v[i__ + i__ * v_dim1] = 1.;
		if (*storev=='C' || *storev=='c') {

		    i__2 = *n - i__ + 1;
		    i__3 = i__ - 1;
		    d__1 = -tau[i__];
		    F77_FUNC(dgemv,DGEMV)("Transpose", &i__2, &i__3, &d__1, &v[i__ + v_dim1],
			     ldv, &v[i__ + i__ * v_dim1], &c__1, &zero, &t[
			    i__ * t_dim1 + 1], &c__1);
		} else {

		    i__2 = i__ - 1;
		    i__3 = *n - i__ + 1;
		    d__1 = -tau[i__];
		    F77_FUNC(dgemv,DGEMV)("No transpose", &i__2, &i__3, &d__1, &v[i__ * 
			    v_dim1 + 1], ldv, &v[i__ + i__ * v_dim1], ldv, &
			    zero, &t[i__ * t_dim1 + 1], &c__1);
		}
		v[i__ + i__ * v_dim1] = vii;


		i__2 = i__ - 1;
		F77_FUNC(dtrmv,DTRMV)("Upper", "No transpose", "Non-unit", &i__2, &t[
			t_offset], ldt, &t[i__ * t_dim1 + 1], &c__1);
		t[i__ + i__ * t_dim1] = tau[i__];
	    }
	}
    } else {
	for (i__ = *k; i__ >= 1; --i__) {
	    if (fabs(tau[i__])<GMX_DOUBLE_MIN) {

		i__1 = *k;
		for (j = i__; j <= i__1; ++j) {
		    t[j + i__ * t_dim1] = 0.;
		}
	    } else {

		if (i__ < *k) {
		    if (*storev=='C' || *storev=='c') {
			vii = v[*n - *k + i__ + i__ * v_dim1];
			v[*n - *k + i__ + i__ * v_dim1] = 1.;

			i__1 = *n - *k + i__;
			i__2 = *k - i__;
			d__1 = -tau[i__];
			F77_FUNC(dgemv,DGEMV)("Transpose", &i__1, &i__2, &d__1, &v[(i__ + 1) 
				* v_dim1 + 1], ldv, &v[i__ * v_dim1 + 1], &
				c__1, &zero, &t[i__ + 1 + i__ * t_dim1], &
				c__1);
			v[*n - *k + i__ + i__ * v_dim1] = vii;
		    } else {
			vii = v[i__ + (*n - *k + i__) * v_dim1];
			v[i__ + (*n - *k + i__) * v_dim1] = 1.;

			i__1 = *k - i__;
			i__2 = *n - *k + i__;
			d__1 = -tau[i__];
			F77_FUNC(dgemv,DGEMV)("No transpose", &i__1, &i__2, &d__1, &v[i__ + 
				1 + v_dim1], ldv, &v[i__ + v_dim1], ldv, &
				zero, &t[i__ + 1 + i__ * t_dim1], &c__1);
			v[i__ + (*n - *k + i__) * v_dim1] = vii;
		    }

		    i__1 = *k - i__;
		    F77_FUNC(dtrmv,DTRMV)("Lower", "No transpose", "Non-unit", &i__1, &t[i__ 
			    + 1 + (i__ + 1) * t_dim1], ldt, &t[i__ + 1 + i__ *
			     t_dim1], &c__1)
			    ;
		}
		t[i__ + i__ * t_dim1] = tau[i__];
	    }
	}
    }
    return;


}
