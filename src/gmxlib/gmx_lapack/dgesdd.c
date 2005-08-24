#include <math.h>
#include <types/simple.h>

#include "gmx_blas.h"
#include "gmx_lapack.h"
#include "lapack_limits.h"


void 
F77_FUNC(dgesdd,DGESDD)(char *jobz, 
	int *m, 
	int *n, 
	double *a, 
	int *lda, 
	double *s,
	double *u, 
	int *ldu, 
	double *vt, 
	int *ldvt, 
	double *work,
	int *lwork, 
	int *iwork, 
	int *info)
{
    int a_dim1, a_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1, i__2;

    int ie, iu;
    double dum[1], eps;
    int ivt, iscl;
    double anrm;
    int idum[1], ierr, itau;
    int minmn, wrkbl, itaup, itauq, mnthr;
    int wntqa;
    int nwork;
    int wntqn, wntqo, wntqs;
    int bdspac;
    double bignum;
    int minwrk, ldwrku, maxwrk, ldwkvt;
    double smlnum,minval, safemin;
    int wntqas, lquery;
    int c__0 = 0;
    int c__1 = 1;
    double zero = 0.0;
    double one = 1.0;


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
    wntqa  = (*jobz=='a' || *jobz=='A');
    wntqs  = (*jobz=='s' || *jobz=='S');
    wntqas = wntqa || wntqs;
    wntqn = (*jobz=='o' || *jobz=='O');
    wntqo = (*jobz=='n' || *jobz=='N');

    minwrk = 1;
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
		    minwrk = bdspac + *n;
		} else {

		    wrkbl = *n * 67;
		    i__1 = wrkbl, i__2 = *n + (*m << 5);
		    wrkbl = (i__1 > i__2) ? i__1 : i__2;
		    i__1 = wrkbl, i__2 = bdspac + *n * 3;
		    wrkbl = (i__1 > i__2) ? i__1 : i__2;
		    maxwrk = wrkbl + *n * *n;
		    minwrk = bdspac + *n * *n + *n * 3;
		}
	    } else {

		wrkbl = *n * 3 + (*m + *n*32);
		if (wntqn) {
		    i__1 = wrkbl, i__2 = bdspac + *n * 3;
		    maxwrk = (i__1 > i__2) ? i__1 : i__2;
		    minwrk = *n * 3 + ((*m > bdspac) ? *m : bdspac);
		} else {
		    i__1 = maxwrk, i__2 = bdspac + *n * 3;
		    maxwrk = (i__1 > i__2) ? i__1 : i__2;
		    minwrk = *n * 3 + ((*m > bdspac) ? *m : bdspac);
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
		    minwrk = bdspac + *m;
		} else {

		    wrkbl = *m * 67;
		    i__1 = wrkbl, i__2 = *m + (*n*32);
		    wrkbl = (i__1 > i__2) ? i__1 : i__2;

		    i__1 = wrkbl, i__2 = bdspac + *m * 3;
		    wrkbl = (i__1 > i__2) ? i__1 : i__2;
		    maxwrk = wrkbl + *m * *m;
		    minwrk = bdspac + *m * *m + *m * 3;
		}
	    } else {
		wrkbl = *m * 3 + (*m + *n*32);
		if (wntqn) {
		    i__1 = wrkbl, i__2 = bdspac + *m * 3;
		    maxwrk = (i__1 > i__2) ? i__1 : i__2;
		    minwrk = *m * 3 + ((*m > bdspac) ? *m : bdspac);
		} else {
		    i__1 = wrkbl, i__2 = bdspac + *m * 3;
		    maxwrk = (i__1 > i__2) ? i__1 : i__2;
		    minwrk = *m * 3 + ((*m > bdspac) ? *m : bdspac);
		}
	    }
	}
	work[1] = (double) maxwrk;
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
    eps = GMX_DOUBLE_EPS;
    minval = GMX_DOUBLE_MIN;
    safemin = minval / eps;
    smlnum = sqrt(safemin) / eps;


    bignum = 1. / smlnum;


    anrm = F77_FUNC(dlange,DLANGE)("M", m, n, &a[a_offset], lda, dum);
    iscl = 0;
    if (anrm > 0. && anrm < smlnum) {
	iscl = 1;
	F77_FUNC(dlascl,DLASCL)("G",&c__0,&c__0,&anrm,&smlnum,m,n,&a[a_offset],lda,&ierr);
    } else if (anrm > bignum) {
	iscl = 1;
	F77_FUNC(dlascl,DLASCL)("G",&c__0,&c__0,&anrm,&bignum,m,n,&a[a_offset],lda,&ierr);
    }

    if (*m >= *n) {
	if (*m >= mnthr) {

	    if (wntqn) {

		itau = 1;
		nwork = itau + *n;

		i__1 = *lwork - nwork + 1;
		F77_FUNC(dgeqrf,DGEQRF)(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

		i__1 = *n - 1;
		i__2 = *n - 1;
		F77_FUNC(dlaset,DLASET)("L", &i__1, &i__2, &zero, &zero, &a[a_dim1 + 2], 
			lda);
		ie = 1;
		itauq = ie + *n;
		itaup = itauq + *n;
		nwork = itaup + *n;

		i__1 = *lwork - nwork + 1;
		F77_FUNC(dgebrd,DGEBRD)(n, n, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);
		nwork = ie + *n;

		F77_FUNC(dbdsdc,DBDSDC)("U", "N", n, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info);

	    } else {
		iu = 1;

		ldwrku = *n;
		itau = iu + ldwrku * *n;
		nwork = itau + *n;

		i__1 = *lwork - nwork + 1;
		F77_FUNC(dgeqrf,DGEQRF)(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);
		F77_FUNC(dlacpy,DLACPY)("L", m, n, &a[a_offset], lda, &u[u_offset], ldu);

		i__1 = *lwork - nwork + 1;
		F77_FUNC(dorgqr,DORGQR)(m, m, n, &u[u_offset], ldu, &work[itau], &work[nwork],
			 &i__1, &ierr);

		i__1 = *n - 1;
		i__2 = *n - 1;
		F77_FUNC(dlaset,DLASET)("L", &i__1, &i__2, &zero, &zero, &a[a_dim1 + 2], 
			lda);
		ie = itau;
		itauq = ie + *n;
		itaup = itauq + *n;
		nwork = itaup + *n;

		i__1 = *lwork - nwork + 1;
		F77_FUNC(dgebrd,DGEBRD)(n, n, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);

		F77_FUNC(dbdsdc,DBDSDC)("U", "I", n, &s[1], &work[ie], &work[iu], n, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info);

		i__1 = *lwork - nwork + 1;
		F77_FUNC(dormbr,DORMBR)("Q", "L", "N", n, n, n, &a[a_offset], lda, &work[
			itauq], &work[iu], &ldwrku, &work[nwork], &i__1, &
			ierr);
		i__1 = *lwork - nwork + 1;
		F77_FUNC(dormbr,DORMBR)("P", "R", "T", n, n, n, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr);

		F77_FUNC(dgemm,DGEMM)("N", "N", m, n, n, &one, &u[u_offset], ldu, &work[iu]
			, &ldwrku, &zero, &a[a_offset], lda);

		F77_FUNC(dlacpy,DLACPY)("F", m, n, &a[a_offset], lda, &u[u_offset], ldu);

	    }

	} else {
	    ie = 1;
	    itauq = ie + *n;
	    itaup = itauq + *n;
	    nwork = itaup + *n;

	    i__1 = *lwork - nwork + 1;
	    F77_FUNC(dgebrd,DGEBRD)(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		    work[itaup], &work[nwork], &i__1, &ierr);
	    if (wntqn) {

		F77_FUNC(dbdsdc,DBDSDC)("U", "N", n, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info);
	    } else {

		F77_FUNC(dlaset,DLASET)("F", m, m, &zero, &zero, &u[u_offset], ldu);
		F77_FUNC(dbdsdc,DBDSDC)("U", "I", n, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info);

		i__1 = *m - *n;
		i__2 = *m - *n;
		F77_FUNC(dlaset,DLASET)("F", &i__1, &i__2, &zero, &one, &u[*n + 1 + (*n + 
			1) * u_dim1], ldu);

		i__1 = *lwork - nwork + 1;
		F77_FUNC(dormbr,DORMBR)("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr);
		i__1 = *lwork - nwork + 1;
		F77_FUNC(dormbr,DORMBR)("P", "R", "T", n, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset],ldvt,&work[nwork],&i__1,&ierr);
	    }

	}

    } else {

	if (*n >= mnthr) {

	    if (wntqn) {

		itau = 1;
		nwork = itau + *m;

		i__1 = *lwork - nwork + 1;
		F77_FUNC(dgelqf,DGELQF)(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);

		i__1 = *m - 1;
		i__2 = *m - 1;
		F77_FUNC(dlaset,DLASET)("U", &i__1, &i__2, &zero, &zero, &a[(a_dim1*2) + 
			1], lda);
		ie = 1;
		itauq = ie + *m;
		itaup = itauq + *m;
		nwork = itaup + *m;

		i__1 = *lwork - nwork + 1;
		F77_FUNC(dgebrd,DGEBRD)(m, m, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);
		nwork = ie + *m;

		F77_FUNC(dbdsdc,DBDSDC)("U", "N", m, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info);

	    } else {

		ivt = 1;

		ldwkvt = *m;
		itau = ivt + ldwkvt * *m;
		nwork = itau + *m;

		i__1 = *lwork - nwork + 1;
		F77_FUNC(dgelqf,DGELQF)(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &
			i__1, &ierr);
		F77_FUNC(dlacpy,DLACPY)("U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt);

		i__1 = *lwork - nwork + 1;
		F77_FUNC(dorglq,DORGLQ)(n, n, m, &vt[vt_offset], ldvt, &work[itau], &work[
			nwork], &i__1, &ierr);

		i__1 = *m - 1;
		i__2 = *m - 1;
		F77_FUNC(dlaset,DLASET)("U", &i__1, &i__2, &zero, &zero, &a[(a_dim1*2) + 
			1], lda);
		ie = itau;
		itauq = ie + *m;
		itaup = itauq + *m;
		nwork = itaup + *m;

		i__1 = *lwork - nwork + 1;
		F77_FUNC(dgebrd,DGEBRD)(m, m, &a[a_offset], lda, &s[1], &work[ie], &work[
			itauq], &work[itaup], &work[nwork], &i__1, &ierr);

		F77_FUNC(dbdsdc,DBDSDC)("U", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &
			work[ivt], &ldwkvt, dum, idum, &work[nwork], &iwork[1]
			, info);

		i__1 = *lwork - nwork + 1;
		F77_FUNC(dormbr,DORMBR)("Q", "L", "N", m, m, m, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr);
		i__1 = *lwork - nwork + 1;
		F77_FUNC(dormbr,DORMBR)("P", "R", "T", m, m, m, &a[a_offset], lda, &work[
			itaup], &work[ivt], &ldwkvt, &work[nwork], &i__1, &
			ierr);

		F77_FUNC(dgemm,DGEMM)("N", "N", m, n, m, &one, &work[ivt], &ldwkvt, &vt[
			vt_offset], ldvt, &zero, &a[a_offset], lda);

		F77_FUNC(dlacpy,DLACPY)("F", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt);

	    }

	} else {

	    ie = 1;
	    itauq = ie + *m;
	    itaup = itauq + *m;
	    nwork = itaup + *m;

	    i__1 = *lwork - nwork + 1;
	    F77_FUNC(dgebrd,DGEBRD)(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
		    work[itaup], &work[nwork], &i__1, &ierr);
	    if (wntqn) {

		F77_FUNC(dbdsdc,DBDSDC)("L", "N", m, &s[1], &work[ie], dum, &c__1, dum, &c__1,
			 dum, idum, &work[nwork], &iwork[1], info);
	    } else {
		F77_FUNC(dlaset,DLASET)("F", n, n, &zero, &zero, &vt[vt_offset], ldvt);
		F77_FUNC(dbdsdc,DBDSDC)("L", "I", m, &s[1], &work[ie], &u[u_offset], ldu, &vt[
			vt_offset], ldvt, dum, idum, &work[nwork], &iwork[1], 
			info);

		i__1 = *n - *m;
		i__2 = *n - *m;
		F77_FUNC(dlaset,DLASET)("F", &i__1, &i__2, &zero, &one, &vt[*m + 1 + (*m + 
			1) * vt_dim1], ldvt);

		i__1 = *lwork - nwork + 1;
		F77_FUNC(dormbr,DORMBR)("Q", "L", "N", m, m, n, &a[a_offset], lda, &work[
			itauq], &u[u_offset], ldu, &work[nwork], &i__1, &ierr);
		i__1 = *lwork - nwork + 1;
		F77_FUNC(dormbr,DORMBR)("P", "R", "T", n, n, m, &a[a_offset], lda, &work[
			itaup], &vt[vt_offset], ldvt, &work[nwork], &i__1, &
			ierr);
	    }

	}

    }

    if (iscl == 1) {
	if (anrm > bignum) {
	    F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr);
	}
	if (anrm < smlnum) {
	    F77_FUNC(dlascl,DLASCL)("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
		    minmn, &ierr);
	}
    }

    work[1] = (double) maxwrk;

    return;

}


