/* slabad.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "../gmx_lapack.h"

#define log10e 0.43429448190325182765

#ifdef KR_headers
float log();
float r_lg10(x) float *x;
#else
#undef abs
#include <math.h>
#ifdef __cplusplus
extern "C" {
#endif
float r_lg10(float *x)
#endif
{
return( log10e * log(*x) );
}
#ifdef __cplusplus
}
#endif

void
	F77_FUNC(slabad, SLABAD)(float *small, float *large)
{
    /* Builtin functions */
    float r_lg10(float *);
    double sqrt(double);


/*  -- LAPACK auxiliary routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SLABAD takes as input the values computed by SLAMCH for underflow and */
/*  overflow, and returns the square root of each of these values if the */
/*  log of LARGE is sufficiently large.  This subroutine is intended to */
/*  identify machines with a large exponent range, such as the Crays, and */
/*  redefine the underflow and overflow limits to be the square roots of */
/*  the values computed by SLAMCH.  This subroutine is needed because */
/*  SLAMCH does not compensate for poor arithmetic in the upper half of */
/*  the exponent range, as is found on a Cray. */

/*  Arguments */
/*  ========= */

/*  SMALL   (input/output) REAL */
/*          On entry, the underflow threshold as computed by SLAMCH. */
/*          On exit, if LOG10(LARGE) is sufficiently large, the square */
/*          root of SMALL, otherwise unchanged. */

/*  LARGE   (input/output) REAL */
/*          On entry, the overflow threshold as computed by SLAMCH. */
/*          On exit, if LOG10(LARGE) is sufficiently large, the square */
/*          root of LARGE, otherwise unchanged. */

/*  ===================================================================== */

/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     If it looks like we're on a Cray, take the square root of */
/*     SMALL and LARGE to avoid overflow and underflow problems. */

    if (r_lg10(large) > 2e3f) {
	*small = sqrt(*small);
	*large = sqrt(*large);
    }

/*     End of SLABAD */

} /* slabad_ */
