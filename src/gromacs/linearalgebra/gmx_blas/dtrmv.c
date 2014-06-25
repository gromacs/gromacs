#include <math.h>

#include "gromacs/utility/real.h"
#include "../gmx_blas.h"

void 
F77_FUNC(dtrmv,DTRMV)(const char *uplo, 
       const char *trans,
       const char *diag, 
       int *n__, 
       double *a, 
       int *lda__, 
       double *x, 
       int *incx__)
{
    int a_dim1, a_offset, i__1, i__2;

    int i__, j, ix, jx, kx, info;
    double temp;
    int nounit;
    
    int n = *n__;
    int lda = *lda__;
    int incx = *incx__;
    
    a_dim1 = lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;

    info = 0;

    if (n == 0) {
	return;
    }

    nounit = (*diag=='n' || *diag=='N');

    if (incx <= 0) {
	kx = 1 - (n - 1) * incx;
    } else {
	kx = 1;
    }

    if (*trans=='N' || *trans=='n') {

	if (*uplo=='U' || *uplo=='u') {
	    if (incx == 1) {
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    if (fabs(x[j])>GMX_DOUBLE_MIN) {
			temp = x[j];
			i__2 = j - 1;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    x[i__] += temp * a[i__ + j * a_dim1];
			}
			if (nounit) {
			    x[j] *= a[j + j * a_dim1];
			}
		    }
		}
	    } else {
		jx = kx;
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    if (fabs(x[jx])>GMX_DOUBLE_MIN) {
			temp = x[jx];
			ix = kx;
			i__2 = j - 1;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    x[ix] += temp * a[i__ + j * a_dim1];
			    ix += incx;
			}
			if (nounit) {
			    x[jx] *= a[j + j * a_dim1];
			}
		    }
		    jx += incx;
		}
	    }
	} else {
	    if (incx == 1) {
		for (j = n; j >= 1; --j) {
		    if (fabs(x[j])>GMX_DOUBLE_MIN) {
			temp = x[j];
			i__1 = j + 1;
			for (i__ = n; i__ >= i__1; --i__) {
			    x[i__] += temp * a[i__ + j * a_dim1];
			}
			if (nounit) {
			    x[j] *= a[j + j * a_dim1];
			}
		    }
		}
	    } else {
		kx += (n - 1) * incx;
		jx = kx;
		for (j = n; j >= 1; --j) {
		    if (fabs(x[jx])>GMX_DOUBLE_MIN) {
			temp = x[jx];
			ix = kx;
			i__1 = j + 1;
			for (i__ = n; i__ >= i__1; --i__) {
			    x[ix] += temp * a[i__ + j * a_dim1];
			    ix -= incx;
			}
			if (nounit) {
			    x[jx] *= a[j + j * a_dim1];
			}
		    }
		    jx -= incx;
		}
	    }
	}
    } else {

	if (*uplo=='U' || *uplo=='u') {
	    if (incx == 1) {
		for (j = n; j >= 1; --j) {
		    temp = x[j];
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    for (i__ = j - 1; i__ >= 1; --i__) {
			temp += a[i__ + j * a_dim1] * x[i__];
		    }
		    x[j] = temp;
		}
	    } else {
		jx = kx + (n - 1) * incx;
		for (j = n; j >= 1; --j) {
		    temp = x[jx];
		    ix = jx;
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    for (i__ = j - 1; i__ >= 1; --i__) {
			ix -= incx;
			temp += a[i__ + j * a_dim1] * x[ix];
		    }
		    x[jx] = temp;
		    jx -= incx;
		}
	    }
	} else {
	    if (incx == 1) {
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[j];
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    i__2 = n;
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
			temp += a[i__ + j * a_dim1] * x[i__];
		    }
		    x[j] = temp;
		}
	    } else {
		jx = kx;
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[jx];
		    ix = jx;
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    i__2 = n;
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
			ix += incx;
			temp += a[i__ + j * a_dim1] * x[ix];
		    }
		    x[jx] = temp;
		    jx += incx;
		}
	    }
	}
    }

    return;

}


