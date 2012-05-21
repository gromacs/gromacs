#include <math.h>
#include "../gmx_blas.h"
#include "../gmx_lapack.h"

void 
F77_FUNC(dlasd3,DLASD3)(int *nl, 
	int *nr,
	int *sqre, 
	int *k, 
	double *d__, 
	double *q, 
	int *ldq, 
	double *dsigma, 
	double *u, 
	int *ldu, 
	double *u2, 
	int *ldu2, 
	double *vt, 
	int *ldvt, 
	double *vt2, 
	int *ldvt2, 
	int *idxc, 
	int *ctot, 
	double *z__, 
	int *info)
{
    int q_dim1, q_offset, u_dim1, u_offset, u2_dim1, u2_offset, vt_dim1, 
	    vt_offset, vt2_dim1, vt2_offset, i__1, i__2;
    double d__2;

    int i__, j, m, n, jc;
    double rho;
    int nlp1, nlp2, nrp1;
    double temp;
    int ctemp;
    int ktemp;
    int c__1 = 1;
    int c__0 = 0;
    double zero = 0.0;
    double one = 1.0;
    double *p1,*p2,t1,t2;

    p1 = &t1;
    p2 = &t2;

    --d__;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --dsigma;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    u2_dim1 = *ldu2;
    u2_offset = 1 + u2_dim1;
    u2 -= u2_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    vt2_dim1 = *ldvt2;
    vt2_offset = 1 + vt2_dim1;
    vt2 -= vt2_offset;
    --idxc;
    --ctot;
    --z__;

    /* Function Body */
    *info = 0;

    if (*nl < 1) {
	*info = -1;
    } else if (*nr < 1) {
	*info = -2;
    } else if (*sqre != 1 && *sqre != 0) {
	*info = -3;
    }

    n = *nl + *nr + 1;
    m = n + *sqre;
    nlp1 = *nl + 1;
    nlp2 = *nl + 2;

    if (*k == 1) {
	d__[1] = fabs(z__[1]);
	F77_FUNC(dcopy,DCOPY)(&m, &vt2[vt2_dim1 + 1], ldvt2, &vt[vt_dim1 + 1], ldvt);
	if (z__[1] > 0.) {
	    F77_FUNC(dcopy,DCOPY)(&n, &u2[u2_dim1 + 1], &c__1, &u[u_dim1 + 1], &c__1);
	} else {
	    i__1 = n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		u[i__ + u_dim1] = -u2[i__ + u2_dim1];
	    }
	}
	return;
    }

    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
      t1 = dsigma[i__];
      t2 = dsigma[i__];
      /* force store and reload from memory */
      t1 = (*p1) + (*p2) - dsigma[i__];

      dsigma[i__] = t1;
    }

    F77_FUNC(dcopy,DCOPY)(k, &z__[1], &c__1, &q[q_offset], &c__1);

    rho = F77_FUNC(dnrm2,DNRM2)(k, &z__[1], &c__1);
    F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &rho, &one, k, &c__1, &z__[1], k, info);
    rho *= rho;


    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	F77_FUNC(dlasd4,DLASD4)(k, &j, &dsigma[1], &z__[1], &u[j * u_dim1 + 1], &rho, &d__[j],
		 &vt[j * vt_dim1 + 1], info);

	if (*info != 0) {
	    return;
	}
    }

    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z__[i__] = u[i__ + *k * u_dim1] * vt[i__ + *k * vt_dim1];
	i__2 = i__ - 1;
	for (j = 1; j <= i__2; ++j) {
	    z__[i__] *= u[i__ + j * u_dim1] * vt[i__ + j * vt_dim1] / (dsigma[
		    i__] - dsigma[j]) / (dsigma[i__] + dsigma[j]);
	}
	i__2 = *k - 1;
	for (j = i__; j <= i__2; ++j) {
	    z__[i__] *= u[i__ + j * u_dim1] * vt[i__ + j * vt_dim1] / (dsigma[
		    i__] - dsigma[j + 1]) / (dsigma[i__] + dsigma[j + 1]);
	}
	d__2 = sqrt(fabs(z__[i__]));
	z__[i__] = (q[i__ + q_dim1] > 0) ? d__2 : -d__2;
    }

    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	vt[i__ * vt_dim1 + 1] = z__[1] / u[i__ * u_dim1 + 1] / vt[i__ * 
		vt_dim1 + 1];
	u[i__ * u_dim1 + 1] = -1.;
	i__2 = *k;
	for (j = 2; j <= i__2; ++j) {
	    vt[j + i__ * vt_dim1] = z__[j] / u[j + i__ * u_dim1] / vt[j + i__ 
		    * vt_dim1];
	    u[j + i__ * u_dim1] = dsigma[j] * vt[j + i__ * vt_dim1];
	}
	temp = F77_FUNC(dnrm2,DNRM2)(k, &u[i__ * u_dim1 + 1], &c__1);
	q[i__ * q_dim1 + 1] = u[i__ * u_dim1 + 1] / temp;
	i__2 = *k;
	for (j = 2; j <= i__2; ++j) {
	    jc = idxc[j];
	    q[j + i__ * q_dim1] = u[jc + i__ * u_dim1] / temp;
	}
    }

    if (*k == 2) {
	F77_FUNC(dgemm,DGEMM)("N", "N", &n, k, k, &one, &u2[u2_offset], ldu2, &q[q_offset],
		 ldq, &zero, &u[u_offset], ldu);
	goto L100;
    }
    if (ctot[1] > 0) {
	F77_FUNC(dgemm,DGEMM)("N", "N", nl, k, &ctot[1], &one, &u2[(u2_dim1 << 1) + 1], 
		ldu2, &q[q_dim1 + 2], ldq, &zero, &u[u_dim1 + 1], ldu);
	if (ctot[3] > 0) {
	    ktemp = ctot[1] + 2 + ctot[2];
	    F77_FUNC(dgemm,DGEMM)("N", "N", nl, k, &ctot[3], &one, &u2[ktemp * u2_dim1 + 1]
		    , ldu2, &q[ktemp + q_dim1], ldq, &one, &u[u_dim1 + 1], 
		    ldu);
	}
    } else if (ctot[3] > 0) {
	ktemp = ctot[1] + 2 + ctot[2];
	F77_FUNC(dgemm,DGEMM)("N", "N", nl, k, &ctot[3], &one, &u2[ktemp * u2_dim1 + 1], 
		ldu2, &q[ktemp + q_dim1], ldq, &zero, &u[u_dim1 + 1], ldu);
    } else {
	F77_FUNC(dlacpy,DLACPY)("F", nl, k, &u2[u2_offset], ldu2, &u[u_offset], ldu);
    }
    F77_FUNC(dcopy,DCOPY)(k, &q[q_dim1 + 1], ldq, &u[nlp1 + u_dim1], ldu);
    ktemp = ctot[1] + 2;
    ctemp = ctot[2] + ctot[3];
    F77_FUNC(dgemm,DGEMM)("N", "N", nr, k, &ctemp, &one, &u2[nlp2 + ktemp * u2_dim1], ldu2,
	     &q[ktemp + q_dim1], ldq, &zero, &u[nlp2 + u_dim1], ldu);

L100:
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	temp = F77_FUNC(dnrm2,DNRM2)(k, &vt[i__ * vt_dim1 + 1], &c__1);
	q[i__ + q_dim1] = vt[i__ * vt_dim1 + 1] / temp;
	i__2 = *k;
	for (j = 2; j <= i__2; ++j) {
	    jc = idxc[j];
	    q[i__ + j * q_dim1] = vt[jc + i__ * vt_dim1] / temp;
	}
    }

    if (*k == 2) {
	F77_FUNC(dgemm,DGEMM)("N", "N", k, &m, k, &one, &q[q_offset], ldq, &vt2[vt2_offset]
		, ldvt2, &zero, &vt[vt_offset], ldvt);
	return;
    }
    ktemp = ctot[1] + 1;
    F77_FUNC(dgemm,DGEMM)("N", "N", k, &nlp1, &ktemp, &one, &q[q_dim1 + 1], ldq, &vt2[
	    vt2_dim1 + 1], ldvt2, &zero, &vt[vt_dim1 + 1], ldvt);
    ktemp = ctot[1] + 2 + ctot[2];
    if (ktemp <= *ldvt2) {
	F77_FUNC(dgemm,DGEMM)("N", "N", k, &nlp1, &ctot[3], &one, &q[ktemp * q_dim1 + 1], 
		ldq, &vt2[ktemp + vt2_dim1], ldvt2, &one, &vt[vt_dim1 + 1], 
		ldvt);
    }

    ktemp = ctot[1] + 1;
    nrp1 = *nr + *sqre;
    if (ktemp > 1) {
	i__1 = *k;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    q[i__ + ktemp * q_dim1] = q[i__ + q_dim1];
	}
	i__1 = m;
	for (i__ = nlp2; i__ <= i__1; ++i__) {
	    vt2[ktemp + i__ * vt2_dim1] = vt2[i__ * vt2_dim1 + 1];
	}
    }
    ctemp = ctot[2] + 1 + ctot[3];
    F77_FUNC(dgemm,DGEMM)("N", "N", k, &nrp1, &ctemp, &one, &q[ktemp * q_dim1 + 1], ldq, &
	    vt2[ktemp + nlp2 * vt2_dim1], ldvt2, &zero, &vt[nlp2 * vt_dim1 + 
	    1], ldvt);

    return;


}


