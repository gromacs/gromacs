#include <math.h>
#include "../gmx_lapack.h"

void 
F77_FUNC(slasd1,SLASD1)(int *nl, 
	int *nr, 
	int *sqre, 
	float *d__, 
	float *alpha, 
	float *beta, 
	float *u, 
	int *ldu, 
	float *vt, 
	int *ldvt, 
	int *idxq, 
	int *iwork, 
	float *work, 
	int *info)
{
    int u_dim1, u_offset, vt_dim1, vt_offset, i__1;
    float d__1, d__2;

    int i__, k, m, n, n1, n2, iq, iz, iu2, ldq, idx, ldu2, ivt2, 
	    idxc, idxp, ldvt2;
    int isigma;
    float orgnrm;
    int coltyp;
    int c__0 = 0;
    float one = 1.0;
    int c__1 = 1;
    int c_n1 = -1;

    --d__;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    --idxq;
    --iwork;
    --work;

    *info = 0;

    if (*nl < 1) {
	*info = -1;
    } else if (*nr < 1) {
	*info = -2;
    } else if (*sqre < 0 || *sqre > 1) {
	*info = -3;
    }
    if (*info != 0) {
	i__1 = -(*info);
	return;
    }

    n = *nl + *nr + 1;
    m = n + *sqre;


    ldu2 = n;
    ldvt2 = m;

    iz = 1;
    isigma = iz + m;
    iu2 = isigma + n;
    ivt2 = iu2 + ldu2 * n;
    iq = ivt2 + ldvt2 * m;

    idx = 1;
    idxc = idx + n;
    coltyp = idxc + n;
    idxp = coltyp + n;

    d__1 = fabs(*alpha);
    d__2 = fabs(*beta);
    orgnrm = (d__1>d__2) ? d__1 : d__2;
    d__[*nl + 1] = 0.;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (fabs(d__[i__]) > orgnrm) {
	    orgnrm = fabs(d__[i__]);
	}
    }
    F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &orgnrm, &one, &n, &c__1, &d__[1], &n, info);
    *alpha /= orgnrm;
    *beta /= orgnrm;

    F77_FUNC(slasd2,SLASD2)(nl, nr, sqre, &k, &d__[1], &work[iz], alpha, beta, &u[u_offset], 
	    ldu, &vt[vt_offset], ldvt, &work[isigma], &work[iu2], &ldu2, &
	    work[ivt2], &ldvt2, &iwork[idxp], &iwork[idx], &iwork[idxc], &
	    idxq[1], &iwork[coltyp], info);

    ldq = k;
    F77_FUNC(slasd3,SLASD3)(nl, nr, sqre, &k, &d__[1], &work[iq], &ldq, &work[isigma], &u[
	    u_offset], ldu, &work[iu2], &ldu2, &vt[vt_offset], ldvt, &work[
	    ivt2], &ldvt2, &iwork[idxc], &iwork[coltyp], &work[iz], info);
    if (*info != 0) {
	return;
    }
    F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &one, &orgnrm, &n, &c__1, &d__[1], &n, info);

    n1 = k;
    n2 = n - k;
    F77_FUNC(slamrg,SLAMRG)(&n1, &n2, &d__[1], &c__1, &c_n1, &idxq[1]);

    return;

}
