#include <ctype.h>

#include "../gmx_blas.h"
#include "../gmx_lapack.h"


void 
F77_FUNC(dlasdq,DLASDQ)(const char *uplo,
                        int *sqre,
                        int *n,
                        int *ncvt,
                        int *nru,
                        int *ncc,
                        double *d__,
                        double *e, 
                        double *vt, 
                        int *ldvt, 
                        double *u,
                        int *ldu, 
                        double *c__,
                        int *ldc,
                        double *work, 
                        int *info)
{
    const char xuplo=toupper(*uplo);
    int c_dim1, c_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1, 
	    i__2;
    int c__1 = 1;
    int itmp1,itmp2;
    int i__, j;
    double r__, cs, sn;
    int np1, isub;
    double smin;
    int sqre1;
    int iuplo;
    int rotate;

    --d__;
    --e;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;

    *info = 0;
    iuplo = 0;
    if (xuplo == 'U') {
	iuplo = 1;
    }
    if (xuplo == 'L') {
	iuplo = 2;
    }
    
    itmp1 = (*n > 1) ? *n : 1;
    itmp2 = (*nru > 1) ? *nru : 1;
    if (iuplo == 0) {
	*info = -1;
    } else if (*sqre < 0 || *sqre > 1) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*ncvt < 0) {
	*info = -4;
    } else if (*nru < 0) {
	*info = -5;
    } else if (*ncc < 0) {
	*info = -6;
    } else if ((*ncvt == 0 && *ldvt < 1) || (*ncvt > 0 && *ldvt < itmp1)) {
	*info = -10;
    } else if (*ldu < itmp2) {
	*info = -12;
    } else if ((*ncc == 0 && *ldc < 1) || (*ncc > 0 && *ldc < itmp1)) {
	*info = -14;
    }
    if (*info != 0) {
	return;
    }
    if (*n == 0) {
	return;
    }

    rotate = *ncvt > 0 || *nru > 0 || *ncc > 0;
    np1 = *n + 1;
    sqre1 = *sqre;

    if (iuplo == 1 && sqre1 == 1) {
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    F77_FUNC(dlartg,DLARTG)(&d__[i__], &e[i__], &cs, &sn, &r__);
	    d__[i__] = r__;
	    e[i__] = sn * d__[i__ + 1];
	    d__[i__ + 1] = cs * d__[i__ + 1];
	    if (rotate) {
		work[i__] = cs;
		work[*n + i__] = sn;
	    }
	}
	F77_FUNC(dlartg,DLARTG)(&d__[*n], &e[*n], &cs, &sn, &r__);
	d__[*n] = r__;
	e[*n] = 0.f;
	if (rotate) {
	    work[*n] = cs;
	    work[*n + *n] = sn;
	}
	iuplo = 2;
	sqre1 = 0;

	if (*ncvt > 0) {
	    F77_FUNC(dlasr,DLASR)("L", "V", "F", &np1, ncvt, &work[1], &work[np1], &vt[
		    vt_offset], ldvt);
	}
    }
    if (iuplo == 2) {
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    F77_FUNC(dlartg,DLARTG)(&d__[i__], &e[i__], &cs, &sn, &r__);
	    d__[i__] = r__;
	    e[i__] = sn * d__[i__ + 1];
	    d__[i__ + 1] = cs * d__[i__ + 1];
	    if (rotate) {
		work[i__] = cs;
		work[*n + i__] = sn;
	    }
	}

	if (sqre1 == 1) {
	    F77_FUNC(dlartg,DLARTG)(&d__[*n], &e[*n], &cs, &sn, &r__);
	    d__[*n] = r__;
	    if (rotate) {
		work[*n] = cs;
		work[*n + *n] = sn;
	    }
	}
	if (*nru > 0) {
	    if (sqre1 == 0) {
		F77_FUNC(dlasr,DLASR)("R", "V", "F", nru, n, &work[1], &work[np1], &u[
			u_offset], ldu);
	    } else {
		F77_FUNC(dlasr,DLASR)("R", "V", "F", nru, &np1, &work[1], &work[np1], &u[
			u_offset], ldu);
	    }
	}
	if (*ncc > 0) {
	    if (sqre1 == 0) {
		F77_FUNC(dlasr,DLASR)("L", "V", "F", n, ncc, &work[1], &work[np1], &c__[
			c_offset], ldc);
	    } else {
		F77_FUNC(dlasr,DLASR)("L", "V", "F", &np1, ncc, &work[1], &work[np1], &c__[
			c_offset], ldc);
	    }
	}
    }

    F77_FUNC(dbdsqr,DBDSQR)("U", n, ncvt, nru, ncc, &d__[1], &e[1], &vt[vt_offset], ldvt, &u[
	    u_offset], ldu, &c__[c_offset], ldc, &work[1], info);

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {

	isub = i__;
	smin = d__[i__];
	i__2 = *n;
	for (j = i__ + 1; j <= i__2; ++j) {
	    if (d__[j] < smin) {
		isub = j;
		smin = d__[j];
	    }
	}
	if (isub != i__) {
	    d__[isub] = d__[i__];
	    d__[i__] = smin;
	    if (*ncvt > 0) {
		F77_FUNC(dswap,DSWAP)(ncvt, &vt[isub + vt_dim1], ldvt, &vt[i__ + vt_dim1], 
			ldvt);
	    }
	    if (*nru > 0) {
		F77_FUNC(dswap,DSWAP)(nru, &u[isub * u_dim1 + 1], &c__1, &u[i__ * u_dim1 + 1]
			, &c__1);
	    }
	    if (*ncc > 0) {
		F77_FUNC(dswap,DSWAP)(ncc, &c__[isub + c_dim1], ldc, &c__[i__ + c_dim1], ldc)
			;
	    }
	}
    }

    return;
}


