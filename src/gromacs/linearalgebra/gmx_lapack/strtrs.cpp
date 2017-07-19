/* strtrs.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "../gmx_blas.h"
#include "gmxpre.h"

#include "../gmx_lapack.h"
#include <string.h>

#define min_def(a,b) ((a) <= (b) ? (a) : (b))
#define max_def(a,b) ((a) >= (b) ? (a) : (b))

/* Table of constant values */

static float c_b12 = 1.f;

int
F77_FUNC(strtrs, STRTRS) (const char *uplo, const char *trans, const char *diag, int *n,
							int *nrhs, float *a, int *lda, float *b, int *
							ldb, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, b_dim1, b_offset, i__1;

    /* Local variables */
    int nounit;


/*  -- LAPACK routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  STRTRS solves a triangular system of the form */

/*     A * X = B  or  A**T * X = B, */

/*  where A is a triangular matrix of order N, and B is an N-by-NRHS */
/*  matrix.  A check is made to verify that A is nonsingular. */

/*  Arguments */
/*  ========= */

/*  UPLO    (input) CHARACTER*1 */
/*          = 'U':  A is upper triangular; */
/*          = 'L':  A is lower triangular. */

/*  TRANS   (input) CHARACTER*1 */
/*          Specifies the form of the system of equations: */
/*          = 'N':  A * X = B  (No transpose) */
/*          = 'T':  A**T * X = B  (Transpose) */
/*          = 'C':  A**H * X = B  (Conjugate transpose = Transpose) */

/*  DIAG    (input) CHARACTER*1 */
/*          = 'N':  A is non-unit triangular; */
/*          = 'U':  A is unit triangular. */

/*  N       (input) INTEGER */
/*          The order of the matrix A.  N >= 0. */

/*  NRHS    (input) INTEGER */
/*          The number of right hand sides, i.e., the number of columns */
/*          of the matrix B.  NRHS >= 0. */

/*  A       (input) REAL array, dimension (LDA,N) */
/*          The triangular matrix A.  If UPLO = 'U', the leading N-by-N */
/*          upper triangular part of the array A contains the upper */
/*          triangular matrix, and the strictly lower triangular part of */
/*          A is not referenced.  If UPLO = 'L', the leading N-by-N lower */
/*          triangular part of the array A contains the lower triangular */
/*          matrix, and the strictly upper triangular part of A is not */
/*          referenced.  If DIAG = 'U', the diagonal elements of A are */
/*          also not referenced and are assumed to be 1. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= max_def(1,N). */

/*  B       (input/output) REAL array, dimension (LDB,NRHS) */
/*          On entry, the right hand side matrix B. */
/*          On exit, if INFO = 0, the solution matrix X. */

/*  LDB     (input) INTEGER */
/*          The leading dimension of the array B.  LDB >= max_def(1,N). */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0: if INFO = -i, the i-th argument had an illegal value */
/*          > 0: if INFO = i, the i-th diagonal element of A is zero, */
/*               indicating that the matrix is singular and the solutions */
/*               X have not been computed. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    *info = 0;
    nounit = strcmp(diag, "N");
    if (! strcmp(uplo, "U") && ! strcmp(uplo, "L")) {
	*info = -1;
    } else if (! strcmp(trans, "N") && ! strcmp(trans,
	    "T") && ! strcmp(trans, "C")) {
	*info = -2;
    } else if (! nounit && ! strcmp(diag, "U")) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*nrhs < 0) {
	*info = -5;
    } else if (*lda < max_def(1,*n)) {
	*info = -7;
    } else if (*ldb < max_def(1,*n)) {
	*info = -9;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("STRTRS", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Check for singularity. */

    if (nounit) {
	i__1 = *n;
	for (*info = 1; *info <= i__1; ++(*info)) {
	    if (a[*info + *info * a_dim1] == 0.f) {
		return 0;
	    }
/* L10: */
	}
    }
    *info = 0;

/*     Solve A * x = b  or  A' * x = b. */

    strsm_("Left", uplo, trans, diag, n, nrhs, &c_b12, &a[a_offset], lda, &b[
	    b_offset], ldb);

    return 0;

/*     End of STRTRS */

} /* strtrs_ */
