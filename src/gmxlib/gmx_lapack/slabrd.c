#include <math.h>
#include "gmx_blas.h"
#include "gmx_lapack.h"


void 
F77_FUNC(slabrd,SLABRD)(int *m, 
	int *n, 
	int *nb,
	float *a, 
	int *lda, 
	float *d__,
	float *e,
	float *tauq, 
	float *taup,
	float *x,
	int *ldx,
	float *y,
	int *ldy)
{
    int a_dim1, a_offset, x_dim1, x_offset, y_dim1, y_offset;
    int i__1, i__2, i__3;
    float one = 1.0;
    float minusone = -1.0;
    float zero = 0.0;
    int c__1 = 1;
    int i__;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --tauq;
    --taup;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    y_dim1 = *ldy;
    y_offset = 1 + y_dim1;
    y -= y_offset;

    if (*m <= 0 || *n <= 0) {
	return;
    }

    if (*m >= *n) {

	i__1 = *nb;
	for (i__ = 1; i__ <= i__1; ++i__) {

	    i__2 = *m - i__ + 1;
	    i__3 = i__ - 1;
	    F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__3, &minusone, &a[i__ + a_dim1], lda,
		     &y[i__ + y_dim1], ldy, &one, &a[i__ + i__ * a_dim1], &c__1);
	    i__2 = *m - i__ + 1;
	    i__3 = i__ - 1;
	    F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__3, &minusone, &x[i__ + x_dim1], ldx,
		   &a[i__*a_dim1+1],&c__1,&one,&a[i__+i__*a_dim1],&c__1);

	    i__2 = *m - i__ + 1;
	    i__3 = i__ + 1;
	    if(*m<i__3)
	      i__3 = *m;
	    F77_FUNC(slarfg,SLARFG)(&i__2, &a[i__ + i__ * a_dim1], &a[i__3 + i__ * a_dim1], 
		    &c__1, &tauq[i__]);
	    d__[i__] = a[i__ + i__ * a_dim1];
	    if (i__ < *n) {
		a[i__ + i__ * a_dim1] = 1.;

		i__2 = *m - i__ + 1;
		i__3 = *n - i__;
		F77_FUNC(sgemv,SGEMV)("Transpose", &i__2, &i__3, &one, &a[i__ + (i__ + 1) * 
			a_dim1], lda, &a[i__ + i__ * a_dim1], &c__1, &zero, &
			y[i__ + 1 + i__ * y_dim1], &c__1);
		i__2 = *m - i__ + 1;
		i__3 = i__ - 1;
		F77_FUNC(sgemv,SGEMV)("Transpose", &i__2, &i__3, &one, &a[i__ + a_dim1], 
			lda, &a[i__ + i__ * a_dim1], &c__1, &zero, &y[i__ * 
			y_dim1 + 1], &c__1);
		i__2 = *n - i__;
		i__3 = i__ - 1;
		F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__3, &minusone, &y[i__ + 1 + 
			y_dim1], ldy, &y[i__ * y_dim1 + 1], &c__1, &one, &y[
			i__ + 1 + i__ * y_dim1], &c__1);
		i__2 = *m - i__ + 1;
		i__3 = i__ - 1;
		F77_FUNC(sgemv,SGEMV)("Transpose", &i__2, &i__3, &one, &x[i__ + x_dim1], 
			ldx, &a[i__ + i__ * a_dim1], &c__1, &zero, &y[i__ * 
			y_dim1 + 1], &c__1);
		i__2 = i__ - 1;
		i__3 = *n - i__;
		F77_FUNC(sgemv,SGEMV)("Transpose", &i__2, &i__3, &minusone, &a[(i__ + 1) * 
			a_dim1 + 1], lda, &y[i__ * y_dim1 + 1], &c__1, &one, 
			&y[i__ + 1 + i__ * y_dim1], &c__1);
		i__2 = *n - i__;
		F77_FUNC(sscal,SSCAL)(&i__2, &tauq[i__], &y[i__ + 1 + i__ * y_dim1], &c__1);

		i__2 = *n - i__;
		F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__, &minusone, &y[i__ + 1 + 
			y_dim1], ldy, &a[i__ + a_dim1], lda, &one, &a[i__ + (
			i__ + 1) * a_dim1], lda);
		i__2 = i__ - 1;
		i__3 = *n - i__;
		F77_FUNC(sgemv,SGEMV)("Transpose", &i__2, &i__3, &minusone, &a[(i__ + 1) * 
			a_dim1 + 1], lda, &x[i__ + x_dim1], ldx, &one, &a[
			i__ + (i__ + 1) * a_dim1], lda);

		i__2 = *n - i__;
		i__3 = i__ + 2;
		if(*n<i__3)
		  i__3 = *n;
		F77_FUNC(slarfg,SLARFG)(&i__2, &a[i__ + (i__ + 1) * a_dim1], 
			&a[i__ + i__3 * a_dim1], lda, &taup[i__]);
		e[i__] = a[i__ + (i__ + 1) * a_dim1];
		a[i__ + (i__ + 1) * a_dim1] = 1.;

		i__2 = *m - i__;
		i__3 = *n - i__;
		F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__3, &one, &a[i__ + 1 + (i__ 
			+ 1) * a_dim1], lda, &a[i__ + (i__ + 1) * a_dim1], 
			lda, &zero, &x[i__ + 1 + i__ * x_dim1], &c__1);
		i__2 = *n - i__;
		F77_FUNC(sgemv,SGEMV)("Transpose", &i__2, &i__, &one, &y[i__ + 1 + y_dim1], 
			ldy, &a[i__ + (i__ + 1) * a_dim1], lda, &zero, &x[
			i__ * x_dim1 + 1], &c__1);
		i__2 = *m - i__;
		F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__, &minusone, &a[i__ + 1 + 
			a_dim1], lda, &x[i__ * x_dim1 + 1], &c__1, &one, &x[
			i__ + 1 + i__ * x_dim1], &c__1);
		i__2 = i__ - 1;
		i__3 = *n - i__;
		F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__3, &one, &a[(i__ + 1) * 
			a_dim1 + 1], lda, &a[i__ + (i__ + 1) * a_dim1], lda, &
			zero, &x[i__ * x_dim1 + 1], &c__1);
		i__2 = *m - i__;
		i__3 = i__ - 1;
		F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__3, &minusone, &x[i__ + 1 + 
			x_dim1], ldx, &x[i__ * x_dim1 + 1], &c__1, &one, &x[
			i__ + 1 + i__ * x_dim1], &c__1);
		i__2 = *m - i__;
		F77_FUNC(sscal,SSCAL)(&i__2, &taup[i__], &x[i__ + 1 + i__ * x_dim1], &c__1);
	    }
	}
    } else {

	i__1 = *nb;
	for (i__ = 1; i__ <= i__1; ++i__) {

	    i__2 = *n - i__ + 1;
	    i__3 = i__ - 1;
	    F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__3, &minusone, &y[i__ + y_dim1], ldy,
		     &a[i__ + a_dim1], lda, &one, &a[i__ + i__ * a_dim1],lda);
	    i__2 = i__ - 1;
	    i__3 = *n - i__ + 1;
	    F77_FUNC(sgemv,SGEMV)("Transpose", &i__2, &i__3, &minusone, &a[i__ * a_dim1 + 1], 
		    lda, &x[i__ + x_dim1], ldx, &one,&a[i__+i__*a_dim1],lda);

	    i__2 = *n - i__ + 1;
	    i__3 = i__ + 1;
	    if(*n<i__3)
	      i__3 = *n;
	    F77_FUNC(slarfg,SLARFG)(&i__2, &a[i__ + i__ * a_dim1], 
		    &a[i__ + i__3 * a_dim1], lda, &taup[i__]);
	    d__[i__] = a[i__ + i__ * a_dim1];
	    if (i__ < *m) {
		a[i__ + i__ * a_dim1] = 1.;

		i__2 = *m - i__;
		i__3 = *n - i__ + 1;
		F77_FUNC(sgemv,SGEMV)("No transpose",&i__2,&i__3,&one,&a[i__+1+i__*a_dim1], 
		       lda, &a[i__ + i__ * a_dim1], lda, &zero, 
		       &x[i__ + 1 + i__ * x_dim1], &c__1);
		i__2 = *n - i__ + 1;
		i__3 = i__ - 1;
		F77_FUNC(sgemv,SGEMV)("Transpose", &i__2, &i__3, &one, &y[i__ + y_dim1], 
			ldy, &a[i__ + i__ * a_dim1], lda, &zero, &x[i__ * 
			x_dim1 + 1], &c__1);
		i__2 = *m - i__;
		i__3 = i__ - 1;
		F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__3, &minusone, &a[i__ + 1 + 
			a_dim1], lda, &x[i__ * x_dim1 + 1], &c__1, &one, &x[
			i__ + 1 + i__ * x_dim1], &c__1);
		i__2 = i__ - 1;
		i__3 = *n - i__ + 1;
		F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__3, &one, &a[i__ * a_dim1 + 
			1], lda, &a[i__ + i__ * a_dim1], lda, &zero, &x[i__ *
			 x_dim1 + 1], &c__1);
		i__2 = *m - i__;
		i__3 = i__ - 1;
		F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__3, &minusone, &x[i__ + 1 + 
			x_dim1], ldx, &x[i__ * x_dim1 + 1], &c__1, &one, &x[
			i__ + 1 + i__ * x_dim1], &c__1);
		i__2 = *m - i__;
		F77_FUNC(sscal,SSCAL)(&i__2, &taup[i__], &x[i__ + 1 + i__ * x_dim1], &c__1);

		i__2 = *m - i__;
		i__3 = i__ - 1;
		F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__3, &minusone, &a[i__ + 1 + 
			a_dim1], lda, &y[i__ + y_dim1], ldy, &one, &a[i__ + 
			1 + i__ * a_dim1], &c__1);
		i__2 = *m - i__;
		F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__, &minusone, &x[i__ + 1 + 
			x_dim1], ldx, &a[i__ * a_dim1 + 1], &c__1, &one, &a[
			i__ + 1 + i__ * a_dim1], &c__1);

		i__2 = *m - i__;
		i__3 = i__ + 2;
		if(*m<i__3)
		  i__3 = *m;
		F77_FUNC(slarfg,SLARFG)(&i__2, &a[i__ + 1 + i__ * a_dim1], 
			&a[i__3 + i__ * a_dim1], &c__1, &tauq[i__]);
		e[i__] = a[i__ + 1 + i__ * a_dim1];
		a[i__ + 1 + i__ * a_dim1] = 1.;

		i__2 = *m - i__;
		i__3 = *n - i__;
		F77_FUNC(sgemv,SGEMV)("Transpose", &i__2, &i__3, &one, &a[i__ + 1 + (i__ + 
			1) * a_dim1], lda, &a[i__ + 1 + i__ * a_dim1], &c__1, 
			&zero, &y[i__ + 1 + i__ * y_dim1], &c__1);
		i__2 = *m - i__;
		i__3 = i__ - 1;
		F77_FUNC(sgemv,SGEMV)("Transpose", &i__2, &i__3, &one, &a[i__ + 1 + a_dim1],
			 lda, &a[i__ + 1 + i__ * a_dim1], &c__1, &zero, &y[
			i__ * y_dim1 + 1], &c__1);
		i__2 = *n - i__;
		i__3 = i__ - 1;
		F77_FUNC(sgemv,SGEMV)("No transpose", &i__2, &i__3, &minusone, &y[i__ + 1 + 
			y_dim1], ldy, &y[i__ * y_dim1 + 1], &c__1, &one, &y[
			i__ + 1 + i__ * y_dim1], &c__1);
		i__2 = *m - i__;
		F77_FUNC(sgemv,SGEMV)("Transpose", &i__2, &i__, &one, &x[i__ + 1 + x_dim1], 
			ldx, &a[i__ + 1 + i__ * a_dim1], &c__1, &zero, &y[
			i__ * y_dim1 + 1], &c__1);
		i__2 = *n - i__;
		F77_FUNC(sgemv,SGEMV)("Transpose", &i__, &i__2, &minusone, &a[(i__ + 1) * a_dim1 
			+ 1], lda, &y[i__ * y_dim1 + 1], &c__1, &one, &y[i__ 
			+ 1 + i__ * y_dim1], &c__1);
		i__2 = *n - i__;
		F77_FUNC(sscal,SSCAL)(&i__2, &tauq[i__], &y[i__ + 1 + i__ * y_dim1], &c__1);
	    }
	}
    }
    return;
} 

