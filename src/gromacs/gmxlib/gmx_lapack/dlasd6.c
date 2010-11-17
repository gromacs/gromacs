#include <math.h>
#include "gmx_blas.h"
#include "gmx_lapack.h"

void 
F77_FUNC(dlasd6,DLASD6)(int *icompq, 
	int *nl, 
	int *nr, 
	int *sqre, 
	double *d__, 
	double *vf, 
	double *vl, 
	double *alpha, 
	double *beta, 
	int *idxq, 
	int *perm, 
	int *givptr, 
	int *givcol, 
	int *ldgcol, 
	double *givnum,
	int *ldgnum, 
	double *poles, 
	double *difl, 
	double *difr, 
	double *z__, 
	int *k, 
	double *c__, 
	double *s, 
	double *work, 
	int *iwork, 
	int *info)
{
    int givcol_dim1, givcol_offset, givnum_dim1, givnum_offset, 
	    poles_dim1, poles_offset, i__1;
    double d__1, d__2;

    int i__, m, n, n1, n2, iw, idx, idxc, idxp, ivfw, ivlw;
    int isigma;
    double orgnrm;
    int c__0 = 0;
    double one = 1.0;
    int c__1 = 1;
    int c_n1 = -1;

    --d__;
    --vf;
    --vl;
    --idxq;
    --perm;
    givcol_dim1 = *ldgcol;
    givcol_offset = 1 + givcol_dim1;
    givcol -= givcol_offset;
    poles_dim1 = *ldgnum;
    poles_offset = 1 + poles_dim1;
    poles -= poles_offset;
    givnum_dim1 = *ldgnum;
    givnum_offset = 1 + givnum_dim1;
    givnum -= givnum_offset;
    --difl;
    --difr;
    --z__;
    --work;
    --iwork;

    *info = 0;
    n = *nl + *nr + 1;
    m = n + *sqre;

    isigma = 1;
    iw = isigma + n;
    ivfw = iw + m;
    ivlw = ivfw + m;

    idx = 1;
    idxc = idx + n;
    idxp = idxc + n;

    d__1 = fabs(*alpha); 
    d__2 = fabs(*beta);
    orgnrm = (d__1 > d__2) ? d__1 : d__2;
    d__[*nl + 1] = 0.;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      d__1 = fabs(d__[i__]);
	if (d__1 > orgnrm)
	    orgnrm = d__1;
    }
    F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &orgnrm, &one, &n, &c__1, &d__[1], &n, info);
    *alpha /= orgnrm;
    *beta /= orgnrm;

    F77_FUNC(dlasd7,DLASD7)(icompq, nl, nr, sqre, k, &d__[1], &z__[1], &work[iw], &vf[1], &
	    work[ivfw], &vl[1], &work[ivlw], alpha, beta, &work[isigma], &
	    iwork[idx], &iwork[idxp], &idxq[1], &perm[1], givptr, &givcol[
	    givcol_offset], ldgcol, &givnum[givnum_offset], ldgnum, c__, s, 
	    info);

    F77_FUNC(dlasd8,DLASD8)(icompq, k, &d__[1], &z__[1], &vf[1], &vl[1], &difl[1], &difr[1], 
	    ldgnum, &work[isigma], &work[iw], info);

    if (*icompq == 1) {
	F77_FUNC(dcopy,DCOPY)(k, &d__[1], &c__1, &poles[poles_dim1 + 1], &c__1);
	F77_FUNC(dcopy,DCOPY)(k, &work[isigma], &c__1, &poles[(poles_dim1 << 1) + 1], &c__1);
    }

    F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &one, &orgnrm, &n, &c__1, &d__[1], &n, info);

    n1 = *k;
    n2 = n - *k;
    F77_FUNC(dlamrg,DLAMRG)(&n1, &n2, &d__[1], &c__1, &c_n1, &idxq[1]);

    return;

}


