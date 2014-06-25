#include <math.h>
#include "../gmx_blas.h"
#include "../gmx_lapack.h"
#include "lapack_limits.h"

#include "gromacs/utility/real.h"

void 
F77_FUNC(slasd2,SLASD2)(int *nl, 
                        int *nr, 
                        int *sqre, 
                        int *k, 
                        float *d__, 
                        float *z__, 
                        float *alpha, 
                        float *beta, 
                        float *u, 
                        int *ldu, 
                        float *vt, 
                        int *ldvt, 
                        float *dsigma, 
                        float *u2, 
                        int *ldu2, 
                        float *vt2, 
                        int *ldvt2, 
                        int *idxp, 
                        int *idx, 
                        int *idxc, 
                        int *idxq, 
                        int *coltyp, 
                        int *info)
{
    int u_dim1, u_offset, u2_dim1, u2_offset, vt_dim1, vt_offset;
    int vt2_dim1, vt2_offset, i__1;
    float d__1, d__2;

    float c__;
    int i__, j, m, n;
    float s;
    int k2;
    float z1;
    int ct, jp;
    float eps, tau, tol;
    int psm[4], nlp1, nlp2, idxi, idxj;
    int ctot[4], idxjp;
    int jprev = 0;
    float hlftol;
    float zero = 0.0;
    int c__1 = 1;


    --d__;
    --z__;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    --dsigma;
    u2_dim1 = *ldu2;
    u2_offset = 1 + u2_dim1;
    u2 -= u2_offset;
    vt2_dim1 = *ldvt2;
    vt2_offset = 1 + vt2_dim1;
    vt2 -= vt2_offset;
    --idxp;
    --idx;
    --idxc;
    --idxq;
    --coltyp;

    *info = 0;

    n = *nl + *nr + 1;
    m = n + *sqre;

    nlp1 = *nl + 1;
    nlp2 = *nl + 2;

    z1 = *alpha * vt[nlp1 + nlp1 * vt_dim1];
    z__[1] = z1;
    for (i__ = *nl; i__ >= 1; --i__) {
	z__[i__ + 1] = *alpha * vt[i__ + nlp1 * vt_dim1];
	d__[i__ + 1] = d__[i__];
	idxq[i__ + 1] = idxq[i__] + 1;
    }

    i__1 = m;
    for (i__ = nlp2; i__ <= i__1; ++i__) {
	z__[i__] = *beta * vt[i__ + nlp2 * vt_dim1];
    }

    i__1 = nlp1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	coltyp[i__] = 1;
    }
    i__1 = n;
    for (i__ = nlp2; i__ <= i__1; ++i__) {
	coltyp[i__] = 2;
    }

    i__1 = n;
    for (i__ = nlp2; i__ <= i__1; ++i__) {
	idxq[i__] += nlp1;
    }

    i__1 = n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	dsigma[i__] = d__[idxq[i__]];
	u2[i__ + u2_dim1] = z__[idxq[i__]];
	idxc[i__] = coltyp[idxq[i__]];
    }

    F77_FUNC(slamrg,SLAMRG)(nl, nr, &dsigma[2], &c__1, &c__1, &idx[2]);

    i__1 = n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	idxi = idx[i__] + 1;
	d__[i__] = dsigma[idxi];
	z__[i__] = u2[idxi + u2_dim1];
	coltyp[i__] = idxc[idxi];
    }

    eps = GMX_FLOAT_EPS;
    d__1 = fabs(*alpha), d__2 = fabs(*beta);
    tol = (d__1 > d__2) ? d__1 : d__2;
    d__2 = fabs(d__[n]);
    tol = eps * 8. * ((d__2 > tol) ? d__2 : tol);

    *k = 1;
    k2 = n + 1;
    i__1 = n;
    for (j = 2; j <= i__1; ++j) {
	if (fabs(z__[j]) <= tol) {

	    --k2;
	    idxp[k2] = j;
	    coltyp[j] = 4;
	    if (j == n) {
		goto L120;
	    }
	} else {
	    jprev = j;
	    goto L90;
	}
    }
L90:
    j = jprev;
L100:
    ++j;
    if (j > n) {
	goto L110;
    }
    if (fabs(z__[j]) <= tol) {

	--k2;
	idxp[k2] = j;
	coltyp[j] = 4;
    } else {

	if (fabs(d__[j] - d__[jprev]) <= tol) {

            s = z__[jprev];
	    c__ = z__[j];

	    tau = F77_FUNC(slapy2,SLAPY2)(&c__, &s);
	    c__ /= tau;
	    s = -s / tau;
	    z__[j] = tau;
	    z__[jprev] = 0.;

	    idxjp = idxq[idx[jprev] + 1];
	    idxj = idxq[idx[j] + 1];
	    if (idxjp <= nlp1) {
		--idxjp;
	    }
	    if (idxj <= nlp1) {
		--idxj;
	    }
	    F77_FUNC(srot,SROT)(&n, &u[idxjp * u_dim1 + 1], &c__1, &u[idxj * u_dim1 + 1], &
		    c__1, &c__, &s);
	    F77_FUNC(srot,SROT)(&m, &vt[idxjp + vt_dim1], ldvt, &vt[idxj + vt_dim1], ldvt, &
		    c__, &s);
	    if (coltyp[j] != coltyp[jprev]) {
		coltyp[j] = 3;
	    }
	    coltyp[jprev] = 4;
	    --k2;
	    idxp[k2] = jprev;
	    jprev = j;
	} else {
	    ++(*k);
	    u2[*k + u2_dim1] = z__[jprev];
	    dsigma[*k] = d__[jprev];
	    idxp[*k] = jprev;
	    jprev = j;
	}
    }
    goto L100;
L110:

    ++(*k);
    u2[*k + u2_dim1] = z__[jprev];
    dsigma[*k] = d__[jprev];
    idxp[*k] = jprev;

L120:

    for (j = 1; j <= 4; ++j) {
	ctot[j - 1] = 0;
    }
    i__1 = n;
    for (j = 2; j <= i__1; ++j) {
	ct = coltyp[j];
	++ctot[ct - 1];
    }

    psm[0] = 2;
    psm[1] = ctot[0] + 2;
    psm[2] = psm[1] + ctot[1];
    psm[3] = psm[2] + ctot[2];

    i__1 = n;
    for (j = 2; j <= i__1; ++j) {
	jp = idxp[j];
	ct = coltyp[jp];
	idxc[psm[ct - 1]] = j;
	++psm[ct - 1];
    }

    i__1 = n;
    for (j = 2; j <= i__1; ++j) {
	jp = idxp[j];
	dsigma[j] = d__[jp];
	idxj = idxq[idx[idxp[idxc[j]]] + 1];
	if (idxj <= nlp1) {
	    --idxj;
	}
	F77_FUNC(scopy,SCOPY)(&n, &u[idxj * u_dim1 + 1], &c__1, &u2[j * u2_dim1 + 1], &c__1);
	F77_FUNC(scopy,SCOPY)(&m, &vt[idxj + vt_dim1], ldvt, &vt2[j + vt2_dim1], ldvt2);
    }

    dsigma[1] = 0.;
    hlftol = tol / 2.;
    if (fabs(dsigma[2]) <= hlftol) {
	dsigma[2] = hlftol;
    }
    if (m > n) {
	z__[1] = F77_FUNC(slapy2,SLAPY2)(&z1, &z__[m]);
	if (z__[1] <= tol) {
	    c__ = 1.;
	    s = 0.;
	    z__[1] = tol;
	} else {
	    c__ = z1 / z__[1];
	    s = z__[m] / z__[1];
	}
    } else {
	if (fabs(z1) <= tol) {
	    z__[1] = tol;
	} else {
	    z__[1] = z1;
	}
    }

    i__1 = *k - 1;
    F77_FUNC(scopy,SCOPY)(&i__1, &u2[u2_dim1 + 2], &c__1, &z__[2], &c__1);

    F77_FUNC(slaset,SLASET)("A", &n, &c__1, &zero, &zero, &u2[u2_offset], ldu2);
    u2[nlp1 + u2_dim1] = 1.;
    if (m > n) {
	i__1 = nlp1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    vt[m + i__ * vt_dim1] = -s * vt[nlp1 + i__ * vt_dim1];
	    vt2[i__ * vt2_dim1 + 1] = c__ * vt[nlp1 + i__ * vt_dim1];
	}
	i__1 = m;
	for (i__ = nlp2; i__ <= i__1; ++i__) {
	    vt2[i__ * vt2_dim1 + 1] = s * vt[m + i__ * vt_dim1];
	    vt[m + i__ * vt_dim1] = c__ * vt[m + i__ * vt_dim1];
	}
    } else {
	F77_FUNC(scopy,SCOPY)(&m, &vt[nlp1 + vt_dim1], ldvt, &vt2[vt2_dim1 + 1], ldvt2);
    }
    if (m > n) {
	F77_FUNC(scopy,SCOPY)(&m, &vt[m + vt_dim1], ldvt, &vt2[m + vt2_dim1], ldvt2);
    }

    if (n > *k) {
	i__1 = n - *k;
	F77_FUNC(scopy,SCOPY)(&i__1, &dsigma[*k + 1], &c__1, &d__[*k + 1], &c__1);
	i__1 = n - *k;
	F77_FUNC(slacpy,SLACPY)("A", &n, &i__1, &u2[(*k + 1) * u2_dim1 + 1], ldu2, &u[(*k + 1)
		 * u_dim1 + 1], ldu);
	i__1 = n - *k;
	F77_FUNC(slacpy,SLACPY)("A", &i__1, &m, &vt2[*k + 1 + vt2_dim1], ldvt2, &vt[*k + 1 + 
		vt_dim1], ldvt);
    }
    for (j = 1; j <= 4; ++j) {
	coltyp[j] = ctot[j - 1];
    }

    return;

}


