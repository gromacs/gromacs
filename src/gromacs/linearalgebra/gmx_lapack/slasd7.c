#include <math.h>
#include "gromacs/utility/real.h"

#include "../gmx_blas.h"
#include "../gmx_lapack.h"
#include "lapack_limits.h"

void 
F77_FUNC(slasd7,SLASD7)(int *icompq, 
	int *nl, 
	int *nr, 
	int *sqre, 
	int *k, 
	float *d__, 
	float *z__, 
	float *zw, 
	float *vf, 
	float *vfw,
	float *vl, 
	float *vlw,
	float *alpha, 
	float *beta,
	float *dsigma, 
	int *idx, 
	int *idxp,
	int *idxq, 
	int *perm, 
	int *givptr,
	int *givcol, 
	int *ldgcol, 
	float *givnum,
	int *ldgnum, 
	float *c__, 
	float *s, 
	int *info)
{
    int givcol_dim1, givcol_offset, givnum_dim1, givnum_offset, i__1;
    float d__1, d__2;

    int i__, j, m, n, k2;
    float z1;
    int jp;
    float eps, tau, tol;
    int nlp1, nlp2, idxi, idxj;
    int idxjp;
    int jprev = 0;
    float hlftol;
    int c__1 = 1;

    --d__;
    --z__;
    --zw;
    --vf;
    --vfw;
    --vl;
    --vlw;
    --dsigma;
    --idx;
    --idxp;
    --idxq;
    --perm;
    givcol_dim1 = *ldgcol;
    givcol_offset = 1 + givcol_dim1;
    givcol -= givcol_offset;
    givnum_dim1 = *ldgnum;
    givnum_offset = 1 + givnum_dim1;
    givnum -= givnum_offset;

    *info = 0;
    n = *nl + *nr + 1;
    m = n + *sqre;

    nlp1 = *nl + 1;
    nlp2 = *nl + 2;
    if (*icompq == 1) {
	*givptr = 0;
    }

    z1 = *alpha * vl[nlp1];
    vl[nlp1] = 0.;
    tau = vf[nlp1];
    for (i__ = *nl; i__ >= 1; --i__) {
	z__[i__ + 1] = *alpha * vl[i__];
	vl[i__] = 0.;
	vf[i__ + 1] = vf[i__];
	d__[i__ + 1] = d__[i__];
	idxq[i__ + 1] = idxq[i__] + 1;
    }
    vf[1] = tau;

    i__1 = m;
    for (i__ = nlp2; i__ <= i__1; ++i__) {
	z__[i__] = *beta * vf[i__];
	vf[i__] = 0.;
    }
    i__1 = n;
    for (i__ = nlp2; i__ <= i__1; ++i__) {
	idxq[i__] += nlp1;
    }

    i__1 = n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	dsigma[i__] = d__[idxq[i__]];
	zw[i__] = z__[idxq[i__]];
	vfw[i__] = vf[idxq[i__]];
	vlw[i__] = vl[idxq[i__]];
    }

    F77_FUNC(slamrg,SLAMRG)(nl, nr, &dsigma[2], &c__1, &c__1, &idx[2]);

    i__1 = n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	idxi = idx[i__] + 1;
	d__[i__] = dsigma[idxi];
	z__[i__] = zw[idxi];
	vf[i__] = vfw[idxi];
	vl[i__] = vlw[idxi];
    }

    eps = GMX_FLOAT_EPS;

    d__1 = fabs(*alpha);
    d__2 = fabs(*beta);
    tol = (d__1>d__2) ? d__1 : d__2;
    d__2 = fabs(d__[n]);
    tol = eps * 64. * ((d__2>tol) ? d__2 : tol);

    *k = 1;
    k2 = n + 1;
    i__1 = n;
    for (j = 2; j <= i__1; ++j) {
	if (fabs(z__[j]) <= tol) {

	    --k2;
	    idxp[k2] = j;
	    if (j == n) {
		goto L100;
	    }
	} else {
	    jprev = j;
	    goto L70;
	}
    }
L70:
    j = jprev;
L80:
    ++j;
    if (j > n) {
	goto L90;
    }
    if (fabs(z__[j]) <= tol) {

	--k2;
	idxp[k2] = j;
    } else {

	if (fabs(d__[j] - d__[jprev]) <= tol) {

	    *s = z__[jprev];
	    *c__ = z__[j];

	    tau = F77_FUNC(slapy2,SLAPY2)(c__, s);
	    z__[j] = tau;
	    z__[jprev] = 0.;
	    *c__ /= tau;
	    *s = -(*s) / tau;


	    if (*icompq == 1) {
		++(*givptr);
		idxjp = idxq[idx[jprev] + 1];
		idxj = idxq[idx[j] + 1];
		if (idxjp <= nlp1) {
		    --idxjp;
		}
		if (idxj <= nlp1) {
		    --idxj;
		}
		givcol[*givptr + (givcol_dim1 << 1)] = idxjp;
		givcol[*givptr + givcol_dim1] = idxj;
		givnum[*givptr + (givnum_dim1 << 1)] = *c__;
		givnum[*givptr + givnum_dim1] = *s;
	    }
	    F77_FUNC(srot,SROT)(&c__1, &vf[jprev], &c__1, &vf[j], &c__1, c__, s);
	    F77_FUNC(srot,SROT)(&c__1, &vl[jprev], &c__1, &vl[j], &c__1, c__, s);
	    --k2;
	    idxp[k2] = jprev;
	    jprev = j;
	} else {
	    ++(*k);
	    zw[*k] = z__[jprev];
	    dsigma[*k] = d__[jprev];
	    idxp[*k] = jprev;
	    jprev = j;
	}
    }
    goto L80;
L90:

    ++(*k);
    zw[*k] = z__[jprev];
    dsigma[*k] = d__[jprev];
    idxp[*k] = jprev;

L100:

    i__1 = n;
    for (j = 2; j <= i__1; ++j) {
	jp = idxp[j];
	dsigma[j] = d__[jp];
	vfw[j] = vf[jp];
	vlw[j] = vl[jp];
    }
    if (*icompq == 1) {
	i__1 = n;
	for (j = 2; j <= i__1; ++j) {
	    jp = idxp[j];
	    perm[j] = idxq[idx[jp] + 1];
	    if (perm[j] <= nlp1) {
		--perm[j];
	    }
	}
    }
    i__1 = n - *k;
    F77_FUNC(scopy,SCOPY)(&i__1, &dsigma[*k + 1], &c__1, &d__[*k + 1], &c__1);

    dsigma[1] = 0.;
    hlftol = tol / 2.;
    if (fabs(dsigma[2]) <= hlftol) {
	dsigma[2] = hlftol;
    }
    if (m > n) {
	z__[1] = F77_FUNC(slapy2,SLAPY2)(&z1, &z__[m]);
	if (z__[1] <= tol) {
	    *c__ = 1.;
	    *s = 0.;
	    z__[1] = tol;
	} else {
	    *c__ = z1 / z__[1];
	    *s = -z__[m] / z__[1];
	}
	F77_FUNC(srot,SROT)(&c__1, &vf[m], &c__1, &vf[1], &c__1, c__, s);
	F77_FUNC(srot,SROT)(&c__1, &vl[m], &c__1, &vl[1], &c__1, c__, s);
    } else {
	if (fabs(z1) <= tol) {
	    z__[1] = tol;
	} else {
	    z__[1] = z1;
	}
    }

    i__1 = *k - 1;
    F77_FUNC(scopy,SCOPY)(&i__1, &zw[2], &c__1, &z__[2], &c__1);
    i__1 = n - 1;
    F77_FUNC(scopy,SCOPY)(&i__1, &vfw[2], &c__1, &vf[2], &c__1);
    i__1 = n - 1;
    F77_FUNC(scopy,SCOPY)(&i__1, &vlw[2], &c__1, &vl[2], &c__1);

    return;

}


