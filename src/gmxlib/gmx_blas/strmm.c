#include <math.h>

#include <types/simple.h>
#include "gmx_blas.h"

void 
F77_FUNC(strmm,STRMM)(char *side, 
                      char *uplo, 
                      char *transa, 
                      char *diag, 
                      int *m__, 
                      int *n__, 
                      float *alpha__, 
                      float *a, 
                      int *lda__, 
                      float *b, 
                      int *ldb__)
{
    int a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    
    int m = *m__;
    int n = *n__;
    int lda = *lda__;
    int ldb = *ldb__;
    float alpha = *alpha__;
    
    /* Local variables */
    int i__, j, k, info;
    float temp;
    int lside;
    int nrowa;
    int upper;
    int nounit;
    a_dim1 = lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    lside = (*side=='L' || *side=='l');
    if (lside) {
	nrowa = m;
    } else {
	nrowa = n;
    }
    nounit = (*diag=='N' || *diag=='n');
    upper = (*uplo=='U' || *uplo=='u');

    info = 0;

    if (n == 0) {
	return;
    }
    if (fabs(alpha)<GMX_FLOAT_MIN) {
	i__1 = n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		b[i__ + j * b_dim1] = 0.;
	    }
	}
	return;
    }
    if (lside) {
	if (*transa=='N' || *transa=='n') {
	    if (upper) {
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = m;
		    for (k = 1; k <= i__2; ++k) {
			if ( fabs(b[k + j * b_dim1])>GMX_FLOAT_MIN) {
			    temp = alpha * b[k + j * b_dim1];
			    i__3 = k - 1;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] += temp * a[i__ + k * 
					a_dim1];
			    }
			    if (nounit) {
				temp *= a[k + k * a_dim1];
			    }
			    b[k + j * b_dim1] = temp;
			}
		    }
		}
	    } else {
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    for (k = m; k >= 1; --k) {
			if (fabs(b[k + j * b_dim1])>GMX_FLOAT_MIN) {
			    temp = alpha * b[k + j * b_dim1];
			    b[k + j * b_dim1] = temp;
			    if (nounit) {
				b[k + j * b_dim1] *= a[k + k * a_dim1];
			    }
			    i__2 = m;
			    for (i__ = k + 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] += temp * a[i__ + k * 
					a_dim1];
			    }
			}
		    }
		}
	    }
	} else {

	    if (upper) {
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    for (i__ = m; i__ >= 1; --i__) {
			temp = b[i__ + j * b_dim1];
			if (nounit) {
			    temp *= a[i__ + i__ * a_dim1];
			}
			i__2 = i__ - 1;
			for (k = 1; k <= i__2; ++k) {
			    temp += a[k + i__ * a_dim1] * b[k + j * b_dim1];
			}
			b[i__ + j * b_dim1] = alpha * temp;
		    }
		}
	    } else {
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			temp = b[i__ + j * b_dim1];
			if (nounit) {
			    temp *= a[i__ + i__ * a_dim1];
			}
			i__3 = m;
			for (k = i__ + 1; k <= i__3; ++k) {
			    temp += a[k + i__ * a_dim1] * b[k + j * b_dim1];
			}
			b[i__ + j * b_dim1] = alpha * temp;
		    }
		}
	    }
	}
    } else {
	if (*transa=='N' || *transa=='n') {

	    if (upper) {
		for (j = n; j >= 1; --j) {
		    temp = alpha;
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    i__1 = m;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
		    }
		    i__1 = j - 1;
		    for (k = 1; k <= i__1; ++k) {
			if ( fabs(a[k + j * a_dim1])>GMX_FLOAT_MIN) {
			    temp = alpha * a[k + j * a_dim1];
			    i__2 = m;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
			    }
			}
		    }
		}
	    } else {
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
		    temp = alpha;
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    i__2 = m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
		    }
		    i__2 = n;
		    for (k = j + 1; k <= i__2; ++k) {
			if ( fabs(a[k + j * a_dim1])>GMX_FLOAT_MIN) {
			    temp = alpha * a[k + j * a_dim1];
			    i__3 = m;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
			    }
			}
		    }
		}
	    }
	} else {

	    if (upper) {
		i__1 = n;
		for (k = 1; k <= i__1; ++k) {
		    i__2 = k - 1;
		    for (j = 1; j <= i__2; ++j) {
			if ( fabs(a[j + k * a_dim1])>GMX_FLOAT_MIN ) {
			    temp = alpha * a[j + k * a_dim1];
			    i__3 = m;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
			    }
			}
		    }
		    temp = alpha;
		    if (nounit) {
			temp *= a[k + k * a_dim1];
		    }
		    if ( fabs(temp-1.0)>GMX_FLOAT_EPS) {
			i__2 = m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
			}
		    }
		}
	    } else {
		for (k = n; k >= 1; --k) {
		    i__1 = n;
		    for (j = k + 1; j <= i__1; ++j) {
			if ( fabs(a[j + k * a_dim1])>GMX_FLOAT_MIN) {
			    temp = alpha * a[j + k * a_dim1];
			    i__2 = m;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
			    }
			}
		    }
		    temp = alpha;
		    if (nounit) {
			temp *= a[k + k * a_dim1];
		    }
		    if ( fabs(temp-1.0)>GMX_FLOAT_EPS) {
			i__1 = m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
			}
		    }
		}
	    }
	}
    }

    return;

}


