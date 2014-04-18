#include <math.h>

#include "gromacs/utility/real.h"
#include "../gmx_lapack.h"

void 
F77_FUNC(dlasr,DLASR)(const char *side, 
       const char *pivot, 
       const char *direct, 
       int *m,
       int *n, 
       double *c__, 
       double *s, 
       double *a, 
       int *lda)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    int i__, j, info;
    double temp;
    double ctemp, stemp;

    --c__;
    --s;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    info = 0;

    if (*m == 0 || *n == 0) {
	return;
    }
    if (*side=='L' || *side=='l') {

	if (*pivot=='V' || *pivot=='v') {
	    if (*direct=='F' || *direct=='f') {
		i__1 = *m - 1;
		for (j = 1; j <= i__1; ++j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (fabs(ctemp-1.0)>GMX_DOUBLE_EPS || fabs(stemp)>GMX_DOUBLE_MIN) {
			i__2 = *n;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a[j + 1 + i__ * a_dim1];
			    a[j + 1 + i__ * a_dim1] = ctemp * temp - stemp * 
				    a[j + i__ * a_dim1];
			    a[j + i__ * a_dim1] = stemp * temp + ctemp * a[j 
				    + i__ * a_dim1];
			}
		    }
		}
	    } else if (*direct=='B' || *direct=='b') {
		for (j = *m - 1; j >= 1; --j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (fabs(ctemp-1.0)>GMX_DOUBLE_EPS || fabs(stemp)>GMX_DOUBLE_MIN) {
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a[j + 1 + i__ * a_dim1];
			    a[j + 1 + i__ * a_dim1] = ctemp * temp - stemp * 
				    a[j + i__ * a_dim1];
			    a[j + i__ * a_dim1] = stemp * temp + ctemp * a[j 
				    + i__ * a_dim1];
			}
		    }
		}
	    }
	} else if (*pivot=='T' || *pivot=='t') {
	    if (*direct=='F' || *direct=='f') {
		i__1 = *m;
		for (j = 2; j <= i__1; ++j) {
		    ctemp = c__[j - 1];
		    stemp = s[j - 1];
		    if (fabs(ctemp-1.0)>GMX_DOUBLE_EPS || fabs(stemp)>GMX_DOUBLE_MIN) {
			i__2 = *n;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a[j + i__ * a_dim1];
			    a[j + i__ * a_dim1] = ctemp * temp - stemp * a[
				    i__ * a_dim1 + 1];
			    a[i__ * a_dim1 + 1] = stemp * temp + ctemp * a[
				    i__ * a_dim1 + 1];
			}
		    }
		}
	    } else if (*direct=='B' || *direct=='b') {
		for (j = *m; j >= 2; --j) {
		    ctemp = c__[j - 1];
		    stemp = s[j - 1];
		    if (fabs(ctemp-1.0)>GMX_DOUBLE_EPS || fabs(stemp)>GMX_DOUBLE_MIN) {
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a[j + i__ * a_dim1];
			    a[j + i__ * a_dim1] = ctemp * temp - stemp * a[
				    i__ * a_dim1 + 1];
			    a[i__ * a_dim1 + 1] = stemp * temp + ctemp * a[
				    i__ * a_dim1 + 1];
			}
		    }
		}
	    }
	} else if (*pivot=='B' || *pivot=='b') {
	    if (*direct=='F' || *direct=='f') {
		i__1 = *m - 1;
		for (j = 1; j <= i__1; ++j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (fabs(ctemp-1.0)>GMX_DOUBLE_EPS || fabs(stemp)>GMX_DOUBLE_MIN) {
			i__2 = *n;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a[j + i__ * a_dim1];
			    a[j + i__ * a_dim1] = stemp * a[*m + i__ * a_dim1]
				     + ctemp * temp;
			    a[*m + i__ * a_dim1] = ctemp * a[*m + i__ * 
				    a_dim1] - stemp * temp;
			}
		    }
		}
	    } else if (*direct=='B' || *direct=='b') {
		for (j = *m - 1; j >= 1; --j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (fabs(ctemp-1.0)>GMX_DOUBLE_EPS || fabs(stemp)>GMX_DOUBLE_MIN) {
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a[j + i__ * a_dim1];
			    a[j + i__ * a_dim1] = stemp * a[*m + i__ * a_dim1]
				     + ctemp * temp;
			    a[*m + i__ * a_dim1] = ctemp * a[*m + i__ * 
				    a_dim1] - stemp * temp;
			}
		    }
		}
	    }
	}
    } else if (*side=='R' || *side=='r') {

	if (*pivot=='V' || *pivot=='v') {
	    if (*direct=='F' || *direct=='f') {
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (fabs(ctemp-1.0)>GMX_DOUBLE_EPS || fabs(stemp)>GMX_DOUBLE_MIN) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a[i__ + (j + 1) * a_dim1];
			    a[i__ + (j + 1) * a_dim1] = ctemp * temp - stemp *
				     a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = stemp * temp + ctemp * a[
				    i__ + j * a_dim1];
			}
		    }
		}
	    } else if (*direct=='B' || *direct=='b') {
		for (j = *n - 1; j >= 1; --j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (fabs(ctemp-1.0)>GMX_DOUBLE_EPS || fabs(stemp)>GMX_DOUBLE_MIN) {
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a[i__ + (j + 1) * a_dim1];
			    a[i__ + (j + 1) * a_dim1] = ctemp * temp - stemp *
				     a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = stemp * temp + ctemp * a[
				    i__ + j * a_dim1];
			}
		    }
		}
	    }
	} else if (*pivot=='T' || *pivot=='t') {
	    if (*direct=='F' || *direct=='f') {
		i__1 = *n;
		for (j = 2; j <= i__1; ++j) {
		    ctemp = c__[j - 1];
		    stemp = s[j - 1];
		    if (fabs(ctemp-1.0)>GMX_DOUBLE_EPS || fabs(stemp)>GMX_DOUBLE_MIN) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = ctemp * temp - stemp * a[
				    i__ + a_dim1];
			    a[i__ + a_dim1] = stemp * temp + ctemp * a[i__ + 
				    a_dim1];
			}
		    }
		}
	    } else if (*direct=='B' || *direct=='b') {
		for (j = *n; j >= 2; --j) {
		    ctemp = c__[j - 1];
		    stemp = s[j - 1];
		    if (fabs(ctemp-1.0)>GMX_DOUBLE_EPS || fabs(stemp)>GMX_DOUBLE_MIN) {
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = ctemp * temp - stemp * a[
				    i__ + a_dim1];
			    a[i__ + a_dim1] = stemp * temp + ctemp * a[i__ + 
				    a_dim1];
			}
		    }
		}
	    }
	} else if (*pivot=='B' || *pivot=='b') {
	    if (*direct=='F' || *direct=='f') {
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (fabs(ctemp-1.0)>GMX_DOUBLE_EPS || fabs(stemp)>GMX_DOUBLE_MIN) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    temp = a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = stemp * a[i__ + *n * a_dim1]
				     + ctemp * temp;
			    a[i__ + *n * a_dim1] = ctemp * a[i__ + *n * 
				    a_dim1] - stemp * temp;
			}
		    }
		}
	    } else if (*direct=='B' || *direct=='b') {
		for (j = *n - 1; j >= 1; --j) {
		    ctemp = c__[j];
		    stemp = s[j];
		    if (fabs(ctemp-1.0)>GMX_DOUBLE_EPS || fabs(stemp)>GMX_DOUBLE_MIN) {
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    temp = a[i__ + j * a_dim1];
			    a[i__ + j * a_dim1] = stemp * a[i__ + *n * a_dim1]
				     + ctemp * temp;
			    a[i__ + *n * a_dim1] = ctemp * a[i__ + *n * 
				    a_dim1] - stemp * temp;
			}
		    }
		}
	    }
	}
    }

    return;

}


