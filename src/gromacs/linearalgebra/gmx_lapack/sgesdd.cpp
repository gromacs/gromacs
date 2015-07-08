#include <cmath>
#include "gromacs/utility/real.h"


#include "../gmx_blas.h"
#include "../gmx_lapack.h"
#include "lapack_limits.h"

void
F77_FUNC(sgesdd,SGESDD)(const char *jobz, 
                        int *m, 
                        int *n, 
                        float *a, 
                        int *lda, 
                        float *s,
                        float *u, 
                        int *ldu, 
                        float *vt, 
                        int *ldvt, 
                        float *work,
                        int *lwork, 
                        int *iwork, 
                        int *info)
{
    int a_dim1, a_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1, i__2;

    int ie, iu;
    float dum[1], eps;
    int ivt, iscl;
    float anrm;
    int idum[1], ierr, itau;
    int minmn, wrkbl, itaup, itauq, mnthr;
    int nwork;
    int wntqn;
    int bdspac;
    float bignum;
    int ldwrku, maxwrk, ldwkvt;
    float smlnum,minval, safemin;
    int lquery;
    int c__0 = 0;
    int c__1 = 1;
    float zero = 0.0;
    float one = 1.0;


    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --s;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    --work;
    --iwork;

    *info = 0;
    minmn = (*m < *n) ? *m : *n;
    mnthr = (int) (minmn * 11. / 6.);
    wntqn = (*jobz=='o' || *jobz=='O');

    maxwrk = 1;
    lquery = *lwork == -1;

    if (*info == 0 && *m > 0 && *n > 0) {
	if (*m >= *n) {

	    if (wntqn) {
		bdspac = *n * 7;
	    } else {
		bdspac = *n * 3 * *n + (*n << 2);
	    }
	    if (*m >= mnthr) {
		if (wntqn) {

		    wrkbl = *n * 67;
		    i__1 = wrkbl, i__2 = bdspac + *n;
		    maxwrk = (i__1 > i__2) ? i__1 : i__2;
		} else {

		    wrkbl = *n * 67;
		    i__1 = wrkbl, i__2 = *n + (*m << 5);
		    wrkbl = (i__1 > i__2) ? i__1 : i__2;
		    i__1 = wrkbl, i__2 = bdspac + *n * 3;
		    wrkbl = (i__1 > i__2) ? i__1 : i__2;
		    maxwrk = wrkbl + *n * *n;
		}
	    } else {

		wrkbl = *n * 3 + (*m + *n*32);
		if (wntqn) {
		    i__1 = wrkbl, i__2 = bdspac + *n * 3;
		    maxwrk = (i__1 > i__2) ? i__1 : i__2;
		} else {
		    i__1 = maxwrk, i__2 = bdspac + *n * 3;
		    maxwrk = (i__1 > i__2) ? i__1 : i__2;
		}
	    }
	} else {

	    if (wntqn) {
		bdspac = *m * 7;
	    } else {
		bdspac = *m * 3 * *m + (*m*4);
	    }
	    if (*n >= mnthr) {
		if (wntqn) {

		    wrkbl = *m * 67;
		    i__1 = wrkbl, i__2 = bdspac + *m;
		    maxwrk = (i__1 > i__2) ? i__1 : i__2;
		} else {

		    wrkbl = *m * 67;
		    i__1 = wrkbl, i__2 = *m + (*n*32);
		    wrkbl = (i__1 > i__2) ? i__1 : i__2;

		    i__1 = wrkbl, i__2 = bdspac + *m * 3;
		    wrkbl = (i__1 > i__2) ? i__1 : i__2;
		    maxwrk = wrkbl + *m * *m;
		}
	    } else {
		wrkbl = *m * 3 + (*m + *n*32);
		if (wntqn) {
		    i__1 = wrkbl, i__2 = bdspac + *m * 3;
		    maxwrk = (i__1 > i__2) ? i__1 : i__2;
		} else {
		    i__1 = wrkbl, i__2 = bdspac + *m * 3;
		    maxwrk = (i__1 > i__2) ? i__1 : i__2;
		}
	    }
	}
	work[1] = (float) maxwrk;
    }
    
    if( lquery != 0)
    {
        return;
    }
    
    if (*m == 0 || *n == 0) {
	if (*lwork >= 1) {
	    work[1] = 1.;
	}
	return;
    }
    eps = GMX_FLOAT_EPS;
    minval = GMX_FLOAT_MIN;
    safemin = minval / eps;
    smlnum =  std::sqrt(safemin) / eps;


    bignum = 1. / smlnum;


    anrm = F77_FUNC(slange,SLANGE)("M", m, n, &a[a_offset], lda, dum);
    iscl = 0;
    if (anrm > 0. && anrm < smlnum) {
	iscl = 1;
	F77_FUNC(slascl,SLASCL)("G",&c__0,&c__0,&anrm,&smlnum,m,n,&a[a_offset],lda,&ierr);
    } else if (anrm > bignum) {
	iscl = 1;
	F77_FUNC(slascl,SLASCL)("G",&c__0,&c__0,&anrm,&bignum,m,n,&a[a_offset],lda,&ierr);
    }

    if (*m >= *n) {
	if (*m >= mnthr) {

	    if (wntqn) {

		itau = 1;
		nwork = itau + *n;

		i__1 = *lwork - nwork + 1;
		F77_FUNC(sgeqrf,SGEQRF)(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

		i__1 = *n - 1;
		i__2 = *n - 1;
		F77_FUNC(slaset,SLASET)("L", &i__1, &i__2, &zero, &zero, &a[a_dim1 + 2], 
			lda);
		ie = 1;
		itauq = ie + *n;
		itaup = itauq + *n;
		nwork = itaup + *n;

		i__1 = *lwork - nwork + 1;
		F77_FUNC(sgebrd,SGEBRD)(n, n, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);
		nwork = ie + *n;

		F77_FUNC(sbdsdc,SBDSDC)("U", "N", n, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info);

	    } else {
		iu = 1;

		ldwrku = *n;
		itau = iu + ldwrku * *n;
		nwork = itau + *n;

		i__1 = *lwork - nwork + 1;
		F77_FUNC(sgeqrf,SGEQRF)(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);
		F77_FUNC(slacpy,SLACPY)("L", m, n, &a[a_offset], lda, &u[u_offset], ldu);

		i__1 = *lwork - nwork + 1;
		F77_FUNC(sorgqr,SORGQR)(m, m, n, &u[u_offset], ldu, &work[itau], &work[nwork],
			 &i__1, &ierr);

		i__1 = *n - 1;
		i__2 = *n - 1;
		F77_FUNC(slaset,SLASET)("L", &i__1, &i__2, &zero, &zero, &a[a_dim1 + 2], 
			lda);
		ie = itau;
		itauq = ie + *n;
		itaup = itauq + *n;
		nwork = itaup + *n;

		i__1 = *lwork - nwork + 1;
		F77_FUNC(sgebrd,SGEBRD)(n, n, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);

		F77_FUNC(sbdsdc,SBDSDC)("U", "I", n, &s[1], &work[ie], &work[iu], n, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info);

		i__1 = *lwork - nwork + 1;
		F77_FUNC(sormbr,SORMBR)("Q", "L", "N", n, n, n, &a[a_offset], lda, &work[
			itauq], &work[iu], &ldwrku, &work[nwork], &i__1, &
			ierr);
		i__1 = *lwork - nwork + 1;
		F77_FUNC(sormbr,SORMBR)("P", "R", "T", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr);

		F77_FUNC(sgemm,SGEMM)("N", "N", m, n, n, &one, &u[u_offset], ldu, &work[iu]
			, &ldwrku, &zero, &a[a_offset], lda);

		F77_FUNC(slacpy,SLACPY)("F", m, n, &a[a_offset], lda, &u[u_offset], ldu);

	    }

	} else {
	    ie = 1;
	    itauq = ie + *n;
	    itaup = itauq + *n;
	    nwork = itaup + *n;

	    i__1 = *lwork - nwork + 1;
	    F77_FUNC(sgebrd,SGEBRD)(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		    work[itaup], &work[nwork], &i__1, &ierr);
	    if (wntqn) {

		F77_FUNC(sbdsdc,SBDSDC)("U", "N", n, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info);
	    } else {

		F77_FUNC(slaset,SLASET)("F", m, m, &zero, &zero, &u[u_offset], ldu);
		F77_FUNC(sbdsdc,SBDSDC)("U", "I", n, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info);

		i__1 = *m - *n;
		i__2 = *m - *n;
		F77_FUNC(slaset,SLASET)("F", &i__1, &i__2, &zero, &one, &u[*n + 1 + (*n + 
			1) * u_dim1], ldu);

		i__1 = *lwork - nwork + 1;
		F77_FUNC(sormbr,SORMBR)("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr);
		i__1 = *lwork - nwork + 1;
		F77_FUNC(sormbr,SORMBR)("P", "R", "T", n, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset],ldvt,&work[nwork],&i__1,&ierr);
	    }

	}

    } else {

	if (*n >= mnthr) {

	    if (wntqn) {

		itau = 1;
		nwork = itau + *m;

		i__1 = *lwork - nwork + 1;
		F77_FUNC(sgelqf,SGELQF)(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

		i__1 = *m - 1;
		i__2 = *m - 1;
		F77_FUNC(slaset,SLASET)("U", &i__1, &i__2, &zero, &zero, &a[(a_dim1*2) + 
			1], lda);
		ie = 1;
		itauq = ie + *m;
		itaup = itauq + *m;
		nwork = itaup + *m;

		i__1 = *lwork - nwork + 1;
		F77_FUNC(sgebrd,SGEBRD)(m, m, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);
		nwork = ie + *m;

		F77_FUNC(sbdsdc,SBDSDC)("U", "N", m, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info);

	    } else {

		ivt = 1;

		ldwkvt = *m;
		itau = ivt + ldwkvt * *m;
		nwork = itau + *m;

		i__1 = *lwork - nwork + 1;
		F77_FUNC(sgelqf,SGELQF)(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);
		F77_FUNC(slacpy,SLACPY)("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt);

		i__1 = *lwork - nwork + 1;
		F77_FUNC(sorglq,SORGLQ)(n, n, m, &vt[vt_offset], ldvt, &work[itau], &work[
			nwork], &i__1, &ierr);

		i__1 = *m - 1;
		i__2 = *m - 1;
		F77_FUNC(slaset,SLASET)("U", &i__1, &i__2, &zero, &zero, &a[(a_dim1*2) + 
			1], lda);
		ie = itau;
		itauq = ie + *m;
		itaup = itauq + *m;
		nwork = itaup + *m;

		i__1 = *lwork - nwork + 1;
		F77_FUNC(sgebrd,SGEBRD)(m, m, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);

		F77_FUNC(sbdsdc,SBDSDC)("U", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &
			work[ivt], &ldwkvt, dum, idum, &work[nwork], &iwork[1]
			, info);

		i__1 = *lwork - nwork + 1;
		F77_FUNC(sormbr,SORMBR)("Q", "L", "N", m, m, m, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr);
		i__1 = *lwork - nwork + 1;
		F77_FUNC(sormbr,SORMBR)("P", "R", "T", m, m, m, &a[a_offset], lda, &work[
			itaup], &work[ivt], &ldwkvt, &work[nwork], &i__1, &
			ierr);

		F77_FUNC(sgemm,SGEMM)("N", "N", m, n, m, &one, &work[ivt], &ldwkvt, &vt[
			vt_offset], ldvt, &zero, &a[a_offset], lda);

		F77_FUNC(slacpy,SLACPY)("F", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt);

	    }

	} else {

	    ie = 1;
	    itauq = ie + *m;
	    itaup = itauq + *m;
	    nwork = itaup + *m;

	    i__1 = *lwork - nwork + 1;
	    F77_FUNC(sgebrd,SGEBRD)(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		    work[itaup], &work[nwork], &i__1, &ierr);
	    if (wntqn) {

		F77_FUNC(sbdsdc,SBDSDC)("L", "N", m, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info);
	    } else {
		F77_FUNC(slaset,SLASET)("F", n, n, &zero, &zero, &vt[vt_offset], ldvt);
		F77_FUNC(sbdsdc,SBDSDC)("L", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info);

		i__1 = *n - *m;
		i__2 = *n - *m;
		F77_FUNC(slaset,SLASET)("F", &i__1, &i__2, &zero, &one, &vt[*m + 1 + (*m + 
			1) * vt_dim1], ldvt);

		i__1 = *lwork - nwork + 1;
		F77_FUNC(sormbr,SORMBR)("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr);
		i__1 = *lwork - nwork + 1;
		F77_FUNC(sormbr,SORMBR)("P", "R", "T", n, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr);
	    }

	}

    }

    if (iscl == 1) {
	if (anrm > bignum) {
	    F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr);
	}
	if (anrm < smlnum) {
	    F77_FUNC(slascl,SLASCL)("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr);
	}
    }

    work[1] = (float) maxwrk;

    return;

}


