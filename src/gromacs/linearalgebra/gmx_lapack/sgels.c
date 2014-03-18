#include <ctype.h>
#include <math.h>
#include "../gmx_blas.h"
#include "../gmx_lapack.h"

#include <types/simple.h>

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))


/* Table of constant values */

static int c__1 = 1;
static int c_n1 = -1;
static float c_b33 = 0.f;
static int c__0 = 0;
static char U = 'U';
static char N = 'N';
static char T = 'T';

void
F77_FUNC(sgels, SGELS) (char *trans, int *m, int *n, int *nrhs, float *a, int *lda, float *b,
                        int *ldb,	float *work, int *lwork, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    int i__, j, nb, mn;
    float anrm, bnrm;
    int brow;
    int tpsd=0;
    int iascl, ibscl;
    int wsize=0;
    float rwork[1];
    int scllen;
    float bignum;
    float smlnum;
    int lquery;

/*  -- LAPACK driver routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SGELS solves overdetermined or underdetermined float linear systems */
/*  involving an M-by-N matrix A, or its transpose, using a QR or LQ */
/*  factorization of A.  It is assumed that A has full rank. */

/*  The following options are provided: */

/*  1. If TRANS = 'N' and m >= n:  find the least squares solution of */
/*     an overdetermined system, i.e., solve the least squares problem */
/*                  minimize || B - A*X ||. */

/*  2. If TRANS = 'N' and m < n:  find the minimum norm solution of */
/*     an underdetermined system A * X = B. */

/*  3. If TRANS = 'T' and m >= n:  find the minimum norm solution of */
/*     an undetermined system A**T * X = B. */

/*  4. If TRANS = 'T' and m < n:  find the least squares solution of */
/*     an overdetermined system, i.e., solve the least squares problem */
/*                  minimize || B - A**T * X ||. */

/*  Several right hand side vectors b and solution vectors x can be */
/*  handled in a single call; they are stored as the columns of the */
/*  M-by-NRHS right hand side matrix B and the N-by-NRHS solution */
/*  matrix X. */

/*  Arguments */
/*  ========= */

/*  TRANS   (input) CHARACTER*1 */
/*          = 'N': the linear system involves A; */
/*          = 'T': the linear system involves A**T. */

/*  M       (input) int */
/*          The number of rows of the matrix A.  M >= 0. */

/*  N       (input) int */
/*          The number of columns of the matrix A.  N >= 0. */

/*  NRHS    (input) int */
/*          The number of right hand sides, i.e., the number of */
/*          columns of the matrices B and X. NRHS >=0. */

/*  A       (input/output) float array, dimension (LDA,N) */
/*          On entry, the M-by-N matrix A. */
/*          On exit, */
/*            if M >= N, A is overwritten by details of its QR */
/*                       factorization as returned by SGEQRF; */
/*            if M <  N, A is overwritten by details of its LQ */
/*                       factorization as returned by SGELQF. */

/*  LDA     (input) int */
/*          The leading dimension of the array A.  LDA >= max(1,M). */

/*  B       (input/output) float array, dimension (LDB,NRHS) */
/*          On entry, the matrix B of right hand side vectors, stored */
/*          columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS */
/*          if TRANS = 'T'. */
/*          On exit, if INFO = 0, B is overwritten by the solution */
/*          vectors, stored columnwise: */
/*          if TRANS = 'N' and m >= n, rows 1 to n of B contain the least */
/*          squares solution vectors; the residual sum of squares for the */
/*          solution in each column is given by the sum of squares of */
/*          elements N+1 to M in that column; */
/*          if TRANS = 'N' and m < n, rows 1 to N of B contain the */
/*          minimum norm solution vectors; */
/*          if TRANS = 'T' and m >= n, rows 1 to M of B contain the */
/*          minimum norm solution vectors; */
/*          if TRANS = 'T' and m < n, rows 1 to M of B contain the */
/*          least squares solution vectors; the residual sum of squares */
/*          for the solution in each column is given by the sum of */
/*          squares of elements M+1 to N in that column. */

/*  LDB     (input) int */
/*          The leading dimension of the array B. LDB >= MAX(1,M,N). */

/*  WORK    (workspace/output) float array, dimension (MAX(1,LWORK)) */
/*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */

/*  LWORK   (input) int */
/*          The dimension of the array WORK. */
/*          LWORK >= max( 1, MN + max( MN, NRHS ) ). */
/*          For optimal performance, */
/*          LWORK >= max( 1, MN + max( MN, NRHS )*NB ). */
/*          where MN = min(M,N) and NB is the optimum block size. */

/*          If LWORK = -1, then a workspace query is assumed; the routine */
/*          only calculates the optimal size of the WORK array, returns */
/*          this value as the first entry of the WORK array, and no error */
/*          message related to LWORK is issued by XERBLA. */

/*  INFO    (output) int */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value */
/*          > 0:  if INFO =  i, the i-th diagonal element of the */
/*                triangular factor of A is zero, so that A does not have */
/*                full rank; the least squares solution could not be */
/*                computed. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments. */

    /* Parameter adjustments */
    const char ch=toupper(*trans);

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --work;

    /* Function Body */
    *info = 0;
    mn = min(*m,*n);
    lquery = *lwork == -1;
    if (! ((ch == 'N') || (ch == 'T'))) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*nrhs < 0) {
	*info = -4;
    } else if (*lda < max(1,*m)) {
	*info = -6;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = max(1,*m);
	if (*ldb < max(i__1,*n)) {
	    *info = -8;
	} else /* if(complicated condition) */ {
/* Computing MAX */
	    i__1 = 1, i__2 = mn + max(mn,*nrhs);
	    if (*lwork < max(i__1,i__2) && ! lquery) {
		*info = -10;
	    }
	}
    }

/*     Figure out optimal block size */

    if (*info == 0 || *info == -10) {

	tpsd = 1;
	if (ch == 'N') {
	    tpsd = 0;
	}

	if (*m >= *n) {
	    nb = F77_FUNC(ilaenv,ILAENV)(&c__1, "SGEQRF", " ", m, n, &c_n1, &c_n1);
	    if (tpsd) {
/* Computing MAX */
		i__1 = nb, i__2 = F77_FUNC(ilaenv,ILAENV)(&c__1, "SORMQR", "LN", m, nrhs, n, &
			c_n1);
		nb = max(i__1,i__2);
	    } else {
/* Computing MAX */
		i__1 = nb, i__2 = F77_FUNC(ilaenv,ILAENV)(&c__1, "SORMQR", "LT", m, nrhs, n, &
			c_n1);
		nb = max(i__1,i__2);
	    }
	} else {
	    nb = F77_FUNC(ilaenv,ILAENV)(&c__1, "SGELQF", " ", m, n, &c_n1, &c_n1);
	    if (tpsd) {
/* Computing MAX */
		i__1 = nb, i__2 = F77_FUNC(ilaenv,ILAENV)(&c__1, "SORMLQ", "LT", n, nrhs, m, &
			c_n1);
		nb = max(i__1,i__2);
	    } else {
/* Computing MAX */
		i__1 = nb, i__2 = F77_FUNC(ilaenv,ILAENV)(&c__1, "SORMLQ", "LN", n, nrhs, m, &
			c_n1);
		nb = max(i__1,i__2);
	    }
	}

/* Computing MAX */
	i__1 = 1, i__2 = mn + max(mn,*nrhs) * nb;
	wsize = max(i__1,i__2);
	work[1] = (float) wsize;

    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SGELS ", &i__1);
	//return;
    } else if (lquery) {
	//return;
    }

/*     Quick return if possible */

/* Computing MIN */
    i__1 = min(*m,*n);
    if (min(i__1,*nrhs) == 0) {
	i__1 = max(*m,*n);
	F77_FUNC(slaset, SLASET) ("Full", &i__1, nrhs, &c_b33, &c_b33, &b[b_offset], ldb);
	//return;
    }

/*     Get machine parameters */

    smlnum = F77_FUNC(slamch, SLAMCH) ("S") / F77_FUNC(slamch, SLAMCH) ("P");
    bignum = 1.f / smlnum;
    F77_FUNC(slabad, SLABAD)(&smlnum, &bignum);

/*     Scale A, B if max element outside range [SMLNUM,BIGNUM] */

    anrm = F77_FUNC(slange, SLANGE)("M", m, n, &a[a_offset], lda, rwork);
    iascl = 0;
    if (anrm > 0.f && anrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

	F77_FUNC(slascl, SLASCL)("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda,
		info);
	iascl = 1;
    } else if (anrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

    F77_FUNC(slascl, SLASCL)("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda,
		info);
	iascl = 2;
    } else if (anrm == 0.f) {

/*        Matrix all zero. Return zero solution. */

	i__1 = max(*m,*n);
	F77_FUNC(slaset, SLASET)("F", &i__1, nrhs, &c_b33, &c_b33, &b[b_offset], ldb);
	goto L50;
    }

    brow = *m;
    if (tpsd) {
	brow = *n;
    }
    bnrm = slange_("M", &brow, nrhs, &b[b_offset], ldb, rwork);
    ibscl = 0;
    if (bnrm > 0.f && bnrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

	F77_FUNC(slascl, SLASCL)("G", &c__0, &c__0, &bnrm, &smlnum, &brow, nrhs, &b[b_offset],
		ldb, info);
	ibscl = 1;
    } else if (bnrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

    F77_FUNC(slascl, SLASCL)("G", &c__0, &c__0, &bnrm, &bignum, &brow, nrhs, &b[b_offset],
		ldb, info);
	ibscl = 2;
    }

    if (*m >= *n) {

/*        compute QR factorization of A */

	i__1 = *lwork - mn;
	F77_FUNC(sgeqrf, SGEQRF)(m, n, &a[a_offset], lda, &work[1], &work[mn + 1], &i__1, info)
		;

/*        workspace at least N, optimally N*NB */

	if (! tpsd) {

/*           Least-Squares Problem min || A * X - B || */

/*           B(1:M,1:NRHS) := Q' * B(1:M,1:NRHS) */

	    i__1 = *lwork - mn;
	    F77_FUNC(sormqr, SORMQR)("Left", "Transpose", m, nrhs, n, &a[a_offset], lda, &work[
		    1], &b[b_offset], ldb, &work[mn + 1], &i__1, info);

/*           workspace at least NRHS, optimally NRHS*NB */

/*           B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS) */

	    F77_FUNC(strtrs, STRTRS)("Upper", "No transpose", "Non-unit", n, nrhs, &a[a_offset]
, lda, &b[b_offset], ldb, info);

	    //if (*info > 0) {
		//return;
	    //}

	    scllen = *n;

	} else {

/*           Overdetermined system of equations A' * X = B */

/*           B(1:N,1:NRHS) := inv(R') * B(1:N,1:NRHS) */

	    F77_FUNC(strtrs, STRTRS)("Upper", "Transpose", "Non-unit", n, nrhs, &a[a_offset],
		    lda, &b[b_offset], ldb, info);

	    //if (*info > 0) {
		//return;
	    //}

/*           B(N+1:M,1:NRHS) = ZERO */

	    i__1 = *nrhs;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = *n + 1; i__ <= i__2; ++i__) {
		    b[i__ + j * b_dim1] = 0.f;
/* L10: */
		}
/* L20: */
	    }

/*           B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS) */

	    i__1 = *lwork - mn;
	    F77_FUNC(sormqr, SORMQR)("Left", "No transpose", m, nrhs, n, &a[a_offset], lda, &
		    work[1], &b[b_offset], ldb, &work[mn + 1], &i__1, info);

/*           workspace at least NRHS, optimally NRHS*NB */

	    scllen = *m;

	}

    } else {

/*        Compute LQ factorization of A */

	i__1 = *lwork - mn;
	F77_FUNC(sgelqf, SGELQF)(m, n, &a[a_offset], lda, &work[1], &work[mn + 1], &i__1, info)
		;

/*        workspace at least M, optimally M*NB. */

	if (! tpsd) {

/*           underdetermined system of equations A * X = B */

/*           B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS) */

		F77_FUNC(strtrs, STRTRS)("Lower", "No transpose", "Non-unit", m, nrhs, &a[a_offset]
, lda, &b[b_offset], ldb, info);

	    //if (*info > 0) {
		//return;
	    //}

/*           B(M+1:N,1:NRHS) = 0 */

	    i__1 = *nrhs;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = *m + 1; i__ <= i__2; ++i__) {
		    b[i__ + j * b_dim1] = 0.f;
/* L30: */
		}
/* L40: */
	    }

/*           B(1:N,1:NRHS) := Q(1:N,:)' * B(1:M,1:NRHS) */

	    i__1 = *lwork - mn;
	    F77_FUNC(sormlq, SORMLQ)("Left", "Transpose", n, nrhs, m, &a[a_offset], lda, &work[
		    1], &b[b_offset], ldb, &work[mn + 1], &i__1, info);

/*           workspace at least NRHS, optimally NRHS*NB */

	    scllen = *n;

	} else {

/*           overdetermined system min || A' * X - B || */

/*           B(1:N,1:NRHS) := Q * B(1:N,1:NRHS) */

	    i__1 = *lwork - mn;
	    F77_FUNC(sormlq, SORMLQ)("Left", "No transpose", n, nrhs, m, &a[a_offset], lda, &
		    work[1], &b[b_offset], ldb, &work[mn + 1], &i__1, info);

/*           workspace at least NRHS, optimally NRHS*NB */

/*           B(1:M,1:NRHS) := inv(L') * B(1:M,1:NRHS) */

	    F77_FUNC(strtrs, STRTRS)("Lower", "Transpose", "Non-unit", m, nrhs, &a[a_offset],
		    lda, &b[b_offset], ldb, info);

	    //if (*info > 0) {
		//return;
	    //}

	    scllen = *m;

	}

    }

/*     Undo scaling */

    if (iascl == 1) {
    F77_FUNC(slascl, SLASCL)("G", &c__0, &c__0, &anrm, &smlnum, &scllen, nrhs, &b[b_offset]
, ldb, info);
    } else if (iascl == 2) {
    F77_FUNC(slascl, SLASCL)("G", &c__0, &c__0, &anrm, &bignum, &scllen, nrhs, &b[b_offset]
, ldb, info);
    }
    if (ibscl == 1) {
    F77_FUNC(slascl, SLASCL)("G", &c__0, &c__0, &smlnum, &bnrm, &scllen, nrhs, &b[b_offset]
, ldb, info);
    } else if (ibscl == 2) {
	F77_FUNC(slascl, SLASCL)("G", &c__0, &c__0, &bignum, &bnrm, &scllen, nrhs, &b[b_offset]
, ldb, info);
    }

L50:
    work[1] = (float) wsize;

    return;

/*     End of SGELS */

} /* sgels_ */
