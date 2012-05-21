#include "../gmx_blas.h"
#include "../gmx_lapack.h"


void 
F77_FUNC(dlarfb,DLARFB)(const char *side, 
	const char *trans, 
	const char *direct, 
	const char *storev, 
	int *m, 
	int *n, 
	int *k, 
	double *v, 
	int *ldv, 
	double *t, 
	int *ldt, 
	double *c__,
	int *ldc, 
	double *work, 
	int *ldwork)
{
    int c_dim1, c_offset, t_dim1, t_offset, v_dim1, v_offset, work_dim1, 
	    work_offset, i__1, i__2;

    int i__, j;
    char transt[1];
    int c__1 = 1;
    double one = 1.0;
    double minusone = -1.0;

    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    work_dim1 = *ldwork;
    work_offset = 1 + work_dim1;
    work -= work_offset;

    if (*m <= 0 || *n <= 0) {
	return;
    }
    if (*trans=='N' || *trans=='n') {
      *(unsigned char *)transt = 'T';
    } else {
	*(unsigned char *)transt = 'N';
    }
    
    if (*storev=='C' || *storev=='c') {

	if (*direct=='F' || *direct=='f') {
	  if (*side=='l' || *side=='L') {

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    F77_FUNC(dcopy,DCOPY)(n, &c__[j + c_dim1], ldc, &work[j * work_dim1 + 1],
			     &c__1);
		}

		F77_FUNC(dtrmm,DTRMM)("Right", "Lower", "No transpose", "Unit", n, k, &one,
			 &v[v_offset], ldv, &work[work_offset], ldwork);
		if (*m > *k) {

		    i__1 = *m - *k;
		    F77_FUNC(dgemm,DGEMM)("Transpose", "No transpose", n, k, &i__1, &one, &
			    c__[*k + 1 + c_dim1], ldc, &v[*k + 1 + v_dim1], 
			    ldv, &one, &work[work_offset], ldwork);
		}

		F77_FUNC(dtrmm,DTRMM)("Right", "Upper", transt, "Non-unit", n, k, &one, &t[
			t_offset], ldt, &work[work_offset], ldwork);

		if (*m > *k) {
		    i__1 = *m - *k;
		    F77_FUNC(dgemm,DGEMM)("No transpose", "Transpose", &i__1, n, k, &minusone, &
			    v[*k + 1 + v_dim1], ldv, &work[work_offset], 
			    ldwork, &one, &c__[*k + 1 + c_dim1], ldc);
		}

		F77_FUNC(dtrmm,DTRMM)("Right", "Lower", "Transpose", "Unit", n, k, &one, &
			v[v_offset], ldv, &work[work_offset], ldwork);

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[j + i__ * c_dim1] -= work[i__ + j * work_dim1];
		    }
		}

	    } else if (*side=='r' || *side=='R') {

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    F77_FUNC(dcopy,DCOPY)(m, &c__[j * c_dim1 + 1], &c__1, &work[j * 
			    work_dim1 + 1], &c__1);
		}

		F77_FUNC(dtrmm,DTRMM)("Right", "Lower", "No transpose", "Unit", m, k, &one,
			 &v[v_offset], ldv, &work[work_offset], ldwork);
		if (*n > *k) {

		    i__1 = *n - *k;
		    F77_FUNC(dgemm,DGEMM)("No transpose", "No transpose", m, k, &i__1, &
			    one, &c__[(*k + 1) * c_dim1 + 1], ldc, &v[*k + 
			    1 + v_dim1], ldv, &one, &work[work_offset], 
			    ldwork);
		}

		F77_FUNC(dtrmm,DTRMM)("Right", "Upper", trans, "Non-unit", m, k, &one, &t[
			t_offset], ldt, &work[work_offset], ldwork);

		if (*n > *k) {
		    i__1 = *n - *k;
		    F77_FUNC(dgemm,DGEMM)("No transpose", "Transpose", m, &i__1, k, &minusone, &
			    work[work_offset], ldwork, &v[*k + 1 + v_dim1], 
			    ldv, &one, &c__[(*k + 1) * c_dim1 + 1], ldc);
		}

		F77_FUNC(dtrmm,DTRMM)("Right", "Lower", "Transpose", "Unit", m, k, &one, &
			v[v_offset], ldv, &work[work_offset], ldwork);

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] -= work[i__ + j * work_dim1];
		    }
		}
	    }

	} else {

	  if (*side=='l' || *side=='L') {
		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    F77_FUNC(dcopy,DCOPY)(n, &c__[*m - *k + j + c_dim1], ldc, &work[j * 
			    work_dim1 + 1], &c__1);
		}

		F77_FUNC(dtrmm,DTRMM)("Right", "Upper", "No transpose", "Unit", n, k, &one,
			 &v[*m - *k + 1 + v_dim1], ldv, &work[work_offset], 
			ldwork);
		if (*m > *k) {
		    i__1 = *m - *k;
		    F77_FUNC(dgemm,DGEMM)("Transpose", "No transpose", n, k, &i__1, &one, &
			    c__[c_offset], ldc, &v[v_offset], ldv, &one, &
			    work[work_offset], ldwork);
		}

		F77_FUNC(dtrmm,DTRMM)("Right", "Lower", transt, "Non-unit", n, k, &one, &t[
			t_offset], ldt, &work[work_offset], ldwork);

		if (*m > *k) {

		    i__1 = *m - *k;
		    F77_FUNC(dgemm,DGEMM)("No transpose", "Transpose", &i__1, n, k, &minusone, &
			    v[v_offset], ldv, &work[work_offset], ldwork, &
			    one, &c__[c_offset], ldc)
			    ;
		}

		F77_FUNC(dtrmm,DTRMM)("Right", "Upper", "Transpose", "Unit", n, k, &one, &
			v[*m - *k + 1 + v_dim1], ldv, &work[work_offset], 
			ldwork);

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[*m - *k + j + i__ * c_dim1] -= work[i__ + j * 
				work_dim1];
		    }
		}

	    } else if (*side=='r' || *side=='R') {
		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    F77_FUNC(dcopy,DCOPY)(m, &c__[(*n - *k + j) * c_dim1 + 1], &c__1, &work[
			    j * work_dim1 + 1], &c__1);
		}

		F77_FUNC(dtrmm,DTRMM)("Right", "Upper", "No transpose", "Unit", m, k, &one,
			 &v[*n - *k + 1 + v_dim1], ldv, &work[work_offset], 
			ldwork);
		if (*n > *k) {
		    i__1 = *n - *k;
		    F77_FUNC(dgemm,DGEMM)("No transpose", "No transpose", m, k, &i__1, &
			    one, &c__[c_offset], ldc, &v[v_offset], ldv, &
			    one, &work[work_offset], ldwork);
		}
		F77_FUNC(dtrmm,DTRMM)("Right", "Lower", trans, "Non-unit", m, k, &one, &t[
			t_offset], ldt, &work[work_offset], ldwork);
		if (*n > *k) {
		    i__1 = *n - *k;
		    F77_FUNC(dgemm,DGEMM)("No transpose", "Transpose", m, &i__1, k, &minusone, &
			    work[work_offset], ldwork, &v[v_offset], ldv, &
			    one, &c__[c_offset], ldc)
			    ;
		}
		F77_FUNC(dtrmm,DTRMM)("Right", "Upper", "Transpose", "Unit", m, k, &one, &
			v[*n - *k + 1 + v_dim1], ldv, &work[work_offset], 
			ldwork);
		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + (*n - *k + j) * c_dim1] -= work[i__ + j * 
				work_dim1];
		    }
		}
	    }
	}

    } else  if (*storev=='r' || *storev=='R') {
      if (*direct=='F' || *direct=='f') {
	  if (*side=='l' || *side=='L') {
		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    F77_FUNC(dcopy,DCOPY)(n, &c__[j + c_dim1], ldc, &work[j * work_dim1 + 1],
			     &c__1);
		}
		F77_FUNC(dtrmm,DTRMM)("Right", "Upper", "Transpose", "Unit", n, k, &one, &
			v[v_offset], ldv, &work[work_offset], ldwork);
		if (*m > *k) {
		    i__1 = *m - *k;
		    F77_FUNC(dgemm,DGEMM)("Transpose", "Transpose", n, k, &i__1, &one, &
			    c__[*k + 1 + c_dim1], ldc, &v[(*k + 1) * v_dim1 + 
			    1], ldv, &one, &work[work_offset], ldwork);
		}

		F77_FUNC(dtrmm,DTRMM)("Right", "Upper", transt, "Non-unit", n, k, &one, &t[
			t_offset], ldt, &work[work_offset], ldwork);
		if (*m > *k) {

		    i__1 = *m - *k;
		    F77_FUNC(dgemm,DGEMM)("Transpose", "Transpose", &i__1, n, k, &minusone, &v[(
			    *k + 1) * v_dim1 + 1], ldv, &work[work_offset], 
			    ldwork, &one, &c__[*k + 1 + c_dim1], ldc);
		}

		F77_FUNC(dtrmm,DTRMM)("Right", "Upper", "No transpose", "Unit", n, k, &one,
			 &v[v_offset], ldv, &work[work_offset], ldwork);

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[j + i__ * c_dim1] -= work[i__ + j * work_dim1];
		    }
		}

	    } else if (*side=='r' || *side=='R') {

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    F77_FUNC(dcopy,DCOPY)(m, &c__[j * c_dim1 + 1], &c__1, &work[j * 
			    work_dim1 + 1], &c__1);
		}

		F77_FUNC(dtrmm,DTRMM)("Right", "Upper", "Transpose", "Unit", m, k, &one, &
			v[v_offset], ldv, &work[work_offset], ldwork);
		if (*n > *k) {

		    i__1 = *n - *k;
		    F77_FUNC(dgemm,DGEMM)("No transpose", "Transpose", m, k, &i__1, &one, &
			    c__[(*k + 1) * c_dim1 + 1], ldc, &v[(*k + 1) * 
			    v_dim1 + 1], ldv, &one, &work[work_offset], 
			    ldwork);
		}

		F77_FUNC(dtrmm,DTRMM)("Right", "Upper", trans, "Non-unit", m, k, &one, &t[
			t_offset], ldt, &work[work_offset], ldwork);

		if (*n > *k) {

		    i__1 = *n - *k;
		    F77_FUNC(dgemm,DGEMM)("No transpose", "No transpose", m, &i__1, k, &
			    minusone, &work[work_offset], ldwork, &v[(*k + 1) * 
			    v_dim1 + 1], ldv, &one, &c__[(*k + 1) * c_dim1 
			    + 1], ldc);
		}
		F77_FUNC(dtrmm,DTRMM)("Right", "Upper", "No transpose", "Unit", m, k, &one,
			 &v[v_offset], ldv, &work[work_offset], ldwork);
		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] -= work[i__ + j * work_dim1];
		    }
		}

	    }

	} else {

	    if (*side=='l' || *side=='L') {

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    F77_FUNC(dcopy,DCOPY)(n, &c__[*m - *k + j + c_dim1], ldc, &work[j * 
			    work_dim1 + 1], &c__1);
		}

		F77_FUNC(dtrmm,DTRMM)("Right", "Lower", "Transpose", "Unit", n, k, &one, &
			v[(*m - *k + 1) * v_dim1 + 1], ldv, &work[work_offset]
			, ldwork);
		if (*m > *k) {

		    i__1 = *m - *k;
		    F77_FUNC(dgemm,DGEMM)("Transpose", "Transpose", n, k, &i__1, &one, &
			    c__[c_offset], ldc, &v[v_offset], ldv, &one, &
			    work[work_offset], ldwork);
		}

		F77_FUNC(dtrmm,DTRMM)("Right", "Lower", transt, "Non-unit", n, k, &one, &t[
			t_offset], ldt, &work[work_offset], ldwork);

		if (*m > *k) {

		    i__1 = *m - *k;
		    F77_FUNC(dgemm,DGEMM)("Transpose", "Transpose", &i__1, n, k, &minusone, &v[
			    v_offset], ldv, &work[work_offset], ldwork, &
			    one, &c__[c_offset], ldc);
		}

		F77_FUNC(dtrmm,DTRMM)("Right", "Lower", "No transpose", "Unit", n, k, &one,
			 &v[(*m - *k + 1) * v_dim1 + 1], ldv, &work[
			work_offset], ldwork);

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[*m - *k + j + i__ * c_dim1] -= work[i__ + j * 
				work_dim1];
		    }
		}

	    } else if (*side=='r' || *side=='R') {

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    F77_FUNC(dcopy,DCOPY)(m, &c__[(*n - *k + j) * c_dim1 + 1], &c__1, &work[
			    j * work_dim1 + 1], &c__1);
		}

		F77_FUNC(dtrmm,DTRMM)("Right", "Lower", "Transpose", "Unit", m, k, &one, &
			v[(*n - *k + 1) * v_dim1 + 1], ldv, &work[work_offset]
			, ldwork);
		if (*n > *k) {

		    i__1 = *n - *k;
		    F77_FUNC(dgemm,DGEMM)("No transpose", "Transpose", m, k, &i__1, &one, &
			    c__[c_offset], ldc, &v[v_offset], ldv, &one, &
			    work[work_offset], ldwork);
		}

		F77_FUNC(dtrmm,DTRMM)("Right", "Lower", trans, "Non-unit", m, k, &one, &t[
			t_offset], ldt, &work[work_offset], ldwork);

		if (*n > *k) {

		    i__1 = *n - *k;
		    F77_FUNC(dgemm,DGEMM)("No transpose", "No transpose", m, &i__1, k, &
			    minusone, &work[work_offset], ldwork, &v[v_offset], 
			    ldv, &one, &c__[c_offset], ldc);
		}

		F77_FUNC(dtrmm,DTRMM)("Right", "Lower", "No transpose", "Unit", m, k, &one,
			 &v[(*n - *k + 1) * v_dim1 + 1], ldv, &work[
			work_offset], ldwork);

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + (*n - *k + j) * c_dim1] -= work[i__ + j * 
				work_dim1];
		    }
		}

	    }

	}
    }

    return;


}

