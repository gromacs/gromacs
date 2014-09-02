#include "../gmx_lapack.h"
#define abs(x) ((x) >= 0 ? (x) : -(x))

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

#ifdef __cplusplus
extern "C" {
#endif

#ifdef KR_headers
float pow_ri(ap, bp) float *ap; int *bp;
#else
float pow_ri(float *ap, int *bp)
#endif
{
float pow, x;
int n;
unsigned long u;

pow = 1;
x = *ap;
n = *bp;

if(n != 0)
	{
	if(n < 0)
		{
		n = -n;
		x = 1/x;
		}
	for(u = n; ; )
		{
		if(u & 01)
			pow *= x;
		if(u >>= 1)
			x *= x;
		else
			break;
		}
	}
return(pow);
}
#ifdef __cplusplus
}
#endif


/* Table of constant values */

static int c__1 = 1;
static float c_b32 = 0.f;

float F77_FUNC(slamch, SLAMCH) (char *cmach)
{
    /* Initialized data */

    static int first = 1;

    /* System generated locals */
    int i__1;
    float ret_val;

    /* Builtin functions */
    float pow_ri(float *, int *);

    /* Local variables */
    static float t;
    int it;
    static float rnd, eps, base;
    int beta;
    static float emin, prec, emax;
    int imin, imax;
    int lrnd;
    static float rmin, rmax;
    float rmach = 0.f;
    float small;
    static float sfmin;
    extern /* Subroutine */ int slamc2_(int *, int *, int *, float
	    *, int *, float *, int *, float *);


/*  -- LAPACK auxiliary routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SLAMCH determines single precision machine parameters. */

/*  Arguments */
/*  ========= */

/*  CMACH   (input) CHARACTER*1 */
/*          Specifies the value to be returned by SLAMCH: */
/*          = 'E' or 'e',   SLAMCH := eps */
/*          = 'S' or 's ,   SLAMCH := sfmin */
/*          = 'B' or 'b',   SLAMCH := base */
/*          = 'P' or 'p',   SLAMCH := eps*base */
/*          = 'N' or 'n',   SLAMCH := t */
/*          = 'R' or 'r',   SLAMCH := rnd */
/*          = 'M' or 'm',   SLAMCH := emin */
/*          = 'U' or 'u',   SLAMCH := rmin */
/*          = 'L' or 'l',   SLAMCH := emax */
/*          = 'O' or 'o',   SLAMCH := rmax */

/*          where */

/*          eps   = relative machine precision */
/*          sfmin = safe minimum, such that 1/sfmin does not overflow */
/*          base  = base of the machine */
/*          prec  = eps*base */
/*          t     = number of (base) digits in the mantissa */
/*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise */
/*          emin  = minimum exponent before (gradual) underflow */
/*          rmin  = underflow threshold - base**(emin-1) */
/*          emax  = largest exponent before overflow */
/*          rmax  = overflow threshold  - (base**emax)*(1-eps) */

/* ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Save statement .. */
/*     .. */
/*     .. Data statements .. */
/*     .. */
/*     .. Executable Statements .. */

    if (first) {
	slamc2_(&beta, &it, &lrnd, &eps, &imin, &rmin, &imax, &rmax);
	base = (float) beta;
	t = (float) it;
	if (lrnd) {
	    rnd = 1.f;
	    i__1 = 1 - it;
	    eps = pow_ri(&base, &i__1) / 2;
	} else {
	    rnd = 0.f;
	    i__1 = 1 - it;
	    eps = pow_ri(&base, &i__1);
	}
	prec = eps * base;
	emin = (float) imin;
	emax = (float) imax;
	sfmin = rmin;
	small = 1.f / rmax;
	if (small >= sfmin) {

/*           Use SMALL plus a bit, to avoid the possibility of rounding */
/*           causing overflow when computing  1/sfmin. */

	    sfmin = small * (eps + 1.f);
	}
    }

    if (*cmach=='E') {
	rmach = eps;
    } else if (*cmach=='S') {
	rmach = sfmin;
    } else if (*cmach=='B') {
	rmach = base;
    } else if (*cmach=='P') {
	rmach = prec;
    } else if (*cmach=='N') {
	rmach = t;
    } else if (*cmach=='R') {
	rmach = rnd;
    } else if (*cmach=='M') {
	rmach = emin;
    } else if (*cmach=='U') {
	rmach = rmin;
    } else if (*cmach=='L') {
	rmach = emax;
    } else if (*cmach=='O') {
	rmach = rmax;
    }

    ret_val = rmach;
    first = 0;
    return ret_val;

/*     End of SLAMCH */

} /* slamch_ */


/* *********************************************************************** */

/* Subroutine */ int slamc1_(int *beta, int *t, int *rnd, int
	*ieee1)
{
    /* Initialized data */

    static int first = 1;

    /* System generated locals */
    float r__1, r__2;

    /* Local variables */
    float a, b, c__, f, t1, t2;
    static int lt;
    float one, qtr;
    static int lrnd;
    static int lbeta;
    float savec;
    static int lieee1;
    extern double slamc3_(float *, float *);


/*  -- LAPACK auxiliary routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SLAMC1 determines the machine parameters given by BETA, T, RND, and */
/*  IEEE1. */

/*  Arguments */
/*  ========= */

/*  BETA    (output) int */
/*          The base of the machine. */

/*  T       (output) int */
/*          The number of ( BETA ) digits in the mantissa. */

/*  RND     (output) int */
/*          Specifies whether proper rounding  ( RND = .TRUE. )  or */
/*          chopping  ( RND = .FALSE. )  occurs in addition. This may not */
/*          be a reliable guide to the way in which the machine performs */
/*          its arithmetic. */

/*  IEEE1   (output) int */
/*          Specifies whether rounding appears to be done in the IEEE */
/*          'round to nearest' style. */

/*  Further Details */
/*  =============== */

/*  The routine is based on the routine  ENVRON  by Malcolm and */
/*  incorporates suggestions by Gentleman and Marovich. See */

/*     Malcolm M. A. (1972) Algorithms to reveal properties of */
/*        floating-point arithmetic. Comms. of the ACM, 15, 949-951. */

/*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms */
/*        that reveal properties of floating point arithmetic units. */
/*        Comms. of the ACM, 17, 276-277. */

/* ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Save statement .. */
/*     .. */
/*     .. Data statements .. */
/*     .. */
/*     .. Executable Statements .. */

    if (first) {
	one = 1.f;

/*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA, */
/*        IEEE1, T and RND. */

/*        Throughout this routine  we use the function  SLAMC3  to ensure */
/*        that relevant values are  stored and not held in registers,  or */
/*        are not affected by optimizers. */

/*        Compute  a = 2.0**m  with the  smallest positive int m such */
/*        that */

/*           fl( a + 1.0 ) = a. */

	a = 1.f;
	c__ = 1.f;

/* +       WHILE( C.EQ.ONE )LOOP */
L10:
	if (c__ == one) {
	    a *= 2;
	    c__ = slamc3_(&a, &one);
	    r__1 = -a;
	    c__ = slamc3_(&c__, &r__1);
	    goto L10;
	}
/* +       END WHILE */

/*        Now compute  b = 2.0**m  with the smallest positive int m */
/*        such that */

/*           fl( a + b ) .gt. a. */

	b = 1.f;
	c__ = slamc3_(&a, &b);

/* +       WHILE( C.EQ.A )LOOP */
L20:
	if (c__ == a) {
	    b *= 2;
	    c__ = slamc3_(&a, &b);
	    goto L20;
	}
/* +       END WHILE */

/*        Now compute the base.  a and c  are neighbouring floating point */
/*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so */
/*        their difference is beta. Adding 0.25 to c is to ensure that it */
/*        is truncated to beta and not ( beta - 1 ). */

	qtr = one / 4;
	savec = c__;
	r__1 = -a;
	c__ = slamc3_(&c__, &r__1);
	lbeta = c__ + qtr;

/*        Now determine whether rounding or chopping occurs,  by adding a */
/*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a. */

	b = (float) lbeta;
	r__1 = b / 2;
	r__2 = -b / 100;
	f = slamc3_(&r__1, &r__2);
	c__ = slamc3_(&f, &a);
	if (c__ == a) {
	    lrnd = 1;
	} else {
	    lrnd = 0;
	}
	r__1 = b / 2;
	r__2 = b / 100;
	f = slamc3_(&r__1, &r__2);
	c__ = slamc3_(&f, &a);
	if (lrnd && c__ == a) {
	    lrnd = 0;
	}

/*        Try and decide whether rounding is done in the  IEEE  'round to */
/*        nearest' style. B/2 is half a unit in the last place of the two */
/*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit */
/*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change */
/*        A, but adding B/2 to SAVEC should change SAVEC. */

	r__1 = b / 2;
	t1 = slamc3_(&r__1, &a);
	r__1 = b / 2;
	t2 = slamc3_(&r__1, &savec);
	lieee1 = t1 == a && t2 > savec && lrnd;

/*        Now find  the  mantissa, t.  It should  be the  int part of */
/*        log to the base beta of a,  however it is safer to determine  t */
/*        by powering.  So we find t as the smallest positive int for */
/*        which */

/*           fl( beta**t + 1.0 ) = 1.0. */

	lt = 0;
	a = 1.f;
	c__ = 1.f;

/* +       WHILE( C.EQ.ONE )LOOP */
L30:
	if (c__ == one) {
	    ++lt;
	    a *= lbeta;
	    c__ = slamc3_(&a, &one);
	    r__1 = -a;
	    c__ = slamc3_(&c__, &r__1);
	    goto L30;
	}
/* +       END WHILE */

    }

    *beta = lbeta;
    *t = lt;
    *rnd = lrnd;
    *ieee1 = lieee1;
    first = 0;
    return 0;

/*     End of SLAMC1 */

} /* slamc1_ */
typedef struct
{	int cierr;
	int ciunit;
	int ciend;
	char *cifmt;
	int cirec;
} cilist;


/* *********************************************************************** */

/* Subroutine */ int slamc2_(int *beta, int *t, int *rnd, float *
	eps, int *emin, float *rmin, int *emax, float *rmax)
{
    /* Initialized data */

    static int first = 1;
    static int iwarn = 0;

    /* Format strings */
    static char fmt_9999[] = "(//\002 WARNING. The value EMIN may be incorre"
	    "ct:-\002,\002  EMIN = \002,i8,/\002 If, after inspection, the va"
	    "lue EMIN looks\002,\002 acceptable please comment out \002,/\002"
	    " the IF block as marked within the code of routine\002,\002 SLAM"
	    "C2,\002,/\002 otherwise supply EMIN explicitly.\002,/)";

    /* System generated locals */
    int i__1;
    float r__1, r__2, r__3, r__4, r__5;

    /* Builtin functions */
    float pow_ri(float *, int *);
    int do_fio(int *, char *, int), e_wsfe(void);

    /* Local variables */
    float a, b, c__;
    int i__;
    static int lt;
    float one, two;
    int ieee;
    float half;
    int lrnd;
    static float leps;
    float zero;
    static int lbeta;
    float rbase;
    static int lemin, lemax;
    int gnmin;
    float small;
    int gpmin;
    float third;
    static float lrmin, lrmax;
    float sixth;
    int lieee1;
    extern /* Subroutine */ int slamc1_(int *, int *, int *,
	    int *);
    extern double slamc3_(float *, float *);
    extern /* Subroutine */ int slamc4_(int *, float *, int *),
	    slamc5_(int *, int *, int *, int *, int *,
	    float *);
    int ngnmin, ngpmin;

    /* Fortran I/O blocks */
    static cilist io___58 = { 0, 6, 0, fmt_9999, 0 };



/*  -- LAPACK auxiliary routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SLAMC2 determines the machine parameters specified in its argument */
/*  list. */

/*  Arguments */
/*  ========= */

/*  BETA    (output) int */
/*          The base of the machine. */

/*  T       (output) int */
/*          The number of ( BETA ) digits in the mantissa. */

/*  RND     (output) int */
/*          Specifies whether proper rounding  ( RND = .TRUE. )  or */
/*          chopping  ( RND = .FALSE. )  occurs in addition. This may not */
/*          be a reliable guide to the way in which the machine performs */
/*          its arithmetic. */

/*  EPS     (output) float */
/*          The smallest positive number such that */

/*             fl( 1.0 - EPS ) .LT. 1.0, */

/*          where fl denotes the computed value. */

/*  EMIN    (output) int */
/*          The minimum exponent before (gradual) underflow occurs. */

/*  RMIN    (output) float */
/*          The smallest normalized number for the machine, given by */
/*          BASE**( EMIN - 1 ), where  BASE  is the floating point value */
/*          of BETA. */

/*  EMAX    (output) int */
/*          The maximum exponent before overflow occurs. */

/*  RMAX    (output) float */
/*          The largest positive number for the machine, given by */
/*          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point */
/*          value of BETA. */

/*  Further Details */
/*  =============== */

/*  The computation of  EPS  is based on a routine PARANOIA by */
/*  W. Kahan of the University of California at Berkeley. */

/* ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Save statement .. */
/*     .. */
/*     .. Data statements .. */
/*     .. */
/*     .. Executable Statements .. */

    if (first) {
	zero = 0.f;
	one = 1.f;
	two = 2.f;

/*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of */
/*        BETA, T, RND, EPS, EMIN and RMIN. */

/*        Throughout this routine  we use the function  SLAMC3  to ensure */
/*        that relevant values are stored  and not held in registers,  or */
/*        are not affected by optimizers. */

/*        SLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1. */

	slamc1_(&lbeta, &lt, &lrnd, &lieee1);

/*        Start to find EPS. */

	b = (float) lbeta;
	i__1 = -lt;
	a = pow_ri(&b, &i__1);
	leps = a;

/*        Try some tricks to see whether or not this is the correct  EPS. */

	b = two / 3;
	half = one / 2;
	r__1 = -half;
	sixth = slamc3_(&b, &r__1);
	third = slamc3_(&sixth, &sixth);
	r__1 = -half;
	b = slamc3_(&third, &r__1);
	b = slamc3_(&b, &sixth);
	b = abs(b);
	if (b < leps) {
	    b = leps;
	}

	leps = 1.f;

/* +       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP */
L10:
	if (leps > b && b > zero) {
	    leps = b;
	    r__1 = half * leps;
/* Computing 5th power */
	    r__3 = two, r__4 = r__3, r__3 *= r__3;
/* Computing 2nd power */
	    r__5 = leps;
	    r__2 = r__4 * (r__3 * r__3) * (r__5 * r__5);
	    c__ = slamc3_(&r__1, &r__2);
	    r__1 = -c__;
	    c__ = slamc3_(&half, &r__1);
	    b = slamc3_(&half, &c__);
	    r__1 = -b;
	    c__ = slamc3_(&half, &r__1);
	    b = slamc3_(&half, &c__);
	    goto L10;
	}
/* +       END WHILE */

	if (a < leps) {
	    leps = a;
	}

/*        Computation of EPS complete. */

/*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)). */
/*        Keep dividing  A by BETA until (gradual) underflow occurs. This */
/*        is detected when we cannot recover the previous A. */

	rbase = one / lbeta;
	small = one;
	for (i__ = 1; i__ <= 3; ++i__) {
	    r__1 = small * rbase;
	    small = slamc3_(&r__1, &zero);
/* L20: */
	}
	a = slamc3_(&one, &small);
	slamc4_(&ngpmin, &one, &lbeta);
	r__1 = -one;
	slamc4_(&ngnmin, &r__1, &lbeta);
	slamc4_(&gpmin, &a, &lbeta);
	r__1 = -a;
	slamc4_(&gnmin, &r__1, &lbeta);
	ieee = 0;

	if (ngpmin == ngnmin && gpmin == gnmin) {
	    if (ngpmin == gpmin) {
		lemin = ngpmin;
/*            ( Non twos-complement machines, no gradual underflow; */
/*              e.g.,  VAX ) */
	    } else if (gpmin - ngpmin == 3) {
		lemin = ngpmin - 1 + lt;
		ieee = 1;
/*            ( Non twos-complement machines, with gradual underflow; */
/*              e.g., IEEE standard followers ) */
	    } else {
		lemin = min(ngpmin,gpmin);
/*            ( A guess; no known machine ) */
		iwarn = 1;
	    }

	} else if (ngpmin == gpmin && ngnmin == gnmin) {
	    if ((i__1 = ngpmin - ngnmin, abs(i__1)) == 1) {
		lemin = max(ngpmin,ngnmin);
/*            ( Twos-complement machines, no gradual underflow; */
/*              e.g., CYBER 205 ) */
	    } else {
		lemin = min(ngpmin,ngnmin);
/*            ( A guess; no known machine ) */
		iwarn = 1;
	    }

	} else if ((i__1 = ngpmin - ngnmin, abs(i__1)) == 1 && gpmin == gnmin)
		 {
	    if (gpmin - min(ngpmin,ngnmin) == 3) {
		lemin = max(ngpmin,ngnmin) - 1 + lt;
/*            ( Twos-complement machines with gradual underflow; */
/*              no known machine ) */
	    } else {
		lemin = min(ngpmin,ngnmin);
/*            ( A guess; no known machine ) */
		iwarn = 1;
	    }

	} else {
/* Computing MIN */
	    i__1 = min(ngpmin,ngnmin), i__1 = min(i__1,gpmin);
	    lemin = min(i__1,gnmin);
/*         ( A guess; no known machine ) */
	    iwarn = 1;
	}
	first = 0;
/* ** */
/* Comment out this if block if EMIN is ok
	if (iwarn) {
	    first = 1;
	    s_wsfe(&io___58);
	    do_fio(&c__1, (char *)&lemin, (int)sizeof(int));
	    e_wsfe();
	}
 ** */

/*        Assume IEEE arithmetic if we found denormalised  numbers above, */
/*        or if arithmetic seems to round in the  IEEE style,  determined */
/*        in routine SLAMC1. A true IEEE machine should have both  things */
/*        true; however, faulty machines may have one or the other. */

	ieee = ieee || lieee1;

/*        Compute  RMIN by successive division by  BETA. We could compute */
/*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during */
/*        this computation. */

	lrmin = 1.f;
	i__1 = 1 - lemin;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    r__1 = lrmin * rbase;
	    lrmin = slamc3_(&r__1, &zero);
/* L30: */
	}

/*        Finally, call SLAMC5 to compute EMAX and RMAX. */

	slamc5_(&lbeta, &lt, &lemin, &ieee, &lemax, &lrmax);
    }

    *beta = lbeta;
    *t = lt;
    *rnd = lrnd;
    *eps = leps;
    *emin = lemin;
    *rmin = lrmin;
    *emax = lemax;
    *rmax = lrmax;

    return 0;


/*     End of SLAMC2 */

} /* slamc2_ */


/* *********************************************************************** */

double slamc3_(float *a, float *b)
{
    /* System generated locals */
    float ret_val;


/*  -- LAPACK auxiliary routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SLAMC3  is intended to force  A  and  B  to be stored prior to doing */
/*  the addition of  A  and  B ,  for use in situations where optimizers */
/*  might hold one of these in a register. */

/*  Arguments */
/*  ========= */

/*  A       (input) float */
/*  B       (input) float */
/*          The values A and B. */

/* ===================================================================== */

/*     .. Executable Statements .. */

    ret_val = *a + *b;

    return ret_val;

/*     End of SLAMC3 */

} /* slamc3_ */


/* *********************************************************************** */

/* Subroutine */ int slamc4_(int *emin, float *start, int *base)
{
    /* System generated locals */
    int i__1;
    float r__1;

    /* Local variables */
    float a;
    int i__;
    float b1, b2, c1, c2, d1, d2, one, zero, rbase;
    extern double slamc3_(float *, float *);


/*  -- LAPACK auxiliary routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SLAMC4 is a service routine for SLAMC2. */

/*  Arguments */
/*  ========= */

/*  EMIN    (output) int */
/*          The minimum exponent before (gradual) underflow, computed by */
/*          setting A = START and dividing by BASE until the previous A */
/*          can not be recovered. */

/*  START   (input) float */
/*          The starting point for determining EMIN. */

/*  BASE    (input) int */
/*          The base of the machine. */

/* ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    a = *start;
    one = 1.f;
    rbase = one / *base;
    zero = 0.f;
    *emin = 1;
    r__1 = a * rbase;
    b1 = slamc3_(&r__1, &zero);
    c1 = a;
    c2 = a;
    d1 = a;
    d2 = a;
/* +    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND. */
/*    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP */
L10:
    if (c1 == a && c2 == a && d1 == a && d2 == a) {
	--(*emin);
	a = b1;
	r__1 = a / *base;
	b1 = slamc3_(&r__1, &zero);
	r__1 = b1 * *base;
	c1 = slamc3_(&r__1, &zero);
	d1 = zero;
	i__1 = *base;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    d1 += b1;
/* L20: */
	}
	r__1 = a * rbase;
	b2 = slamc3_(&r__1, &zero);
	r__1 = b2 / rbase;
	c2 = slamc3_(&r__1, &zero);
	d2 = zero;
	i__1 = *base;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    d2 += b2;
/* L30: */
	}
	goto L10;
    }
/* +    END WHILE */

    return 0;

/*     End of SLAMC4 */

} /* slamc4_ */


/* *********************************************************************** */

/* Subroutine */ int slamc5_(int *beta, int *p, int *emin,
	int *ieee, int *emax, float *rmax)
{
    /* System generated locals */
    int i__1;
    float r__1;

    /* Local variables */
    int i__;
    float y, z__;
    int try__, lexp;
    float oldy = 0.f;
    int uexp, nbits;
    extern double slamc3_(float *, float *);
    float recbas;
    int exbits, expsum;


/*  -- LAPACK auxiliary routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SLAMC5 attempts to compute RMAX, the largest machine floating-point */
/*  number, without overflow.  It assumes that EMAX + abs(EMIN) sum */
/*  approximately to a power of 2.  It will fail on machines where this */
/*  assumption does not hold, for example, the Cyber 205 (EMIN = -28625, */
/*  EMAX = 28718).  It will also fail if the value supplied for EMIN is */
/*  too large (i.e. too close to zero), probably with overflow. */

/*  Arguments */
/*  ========= */

/*  BETA    (input) int */
/*          The base of floating-point arithmetic. */

/*  P       (input) int */
/*          The number of base BETA digits in the mantissa of a */
/*          floating-point value. */

/*  EMIN    (input) int */
/*          The minimum exponent before (gradual) underflow. */

/*  IEEE    (input) int */
/*          A int flag specifying whether or not the arithmetic */
/*          system is thought to comply with the IEEE standard. */

/*  EMAX    (output) int */
/*          The largest exponent before overflow */

/*  RMAX    (output) float */
/*          The largest machine floating-point number. */

/* ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     First compute LEXP and UEXP, two powers of 2 that bound */
/*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum */
/*     approximately to the bound that is closest to abs(EMIN). */
/*     (EMAX is the exponent of the required number RMAX). */

    lexp = 1;
    exbits = 1;
L10:
    try__ = lexp << 1;
    if (try__ <= -(*emin)) {
	lexp = try__;
	++exbits;
	goto L10;
    }
    if (lexp == -(*emin)) {
	uexp = lexp;
    } else {
	uexp = try__;
	++exbits;
    }

/*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater */
/*     than or equal to EMIN. EXBITS is the number of bits needed to */
/*     store the exponent. */

    if (uexp + *emin > -lexp - *emin) {
	expsum = lexp << 1;
    } else {
	expsum = uexp << 1;
    }

/*     EXPSUM is the exponent range, approximately equal to */
/*     EMAX - EMIN + 1 . */

    *emax = expsum + *emin - 1;
    nbits = exbits + 1 + *p;

/*     NBITS is the total number of bits needed to store a */
/*     floating-point number. */

    if (nbits % 2 == 1 && *beta == 2) {

/*        Either there are an odd number of bits used to store a */
/*        floating-point number, which is unlikely, or some bits are */
/*        not used in the representation of numbers, which is possible, */
/*        (e.g. Cray machines) or the mantissa has an implicit bit, */
/*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the */
/*        most likely. We have to assume the last alternative. */
/*        If this is true, then we need to reduce EMAX by one because */
/*        there must be some way of representing zero in an implicit-bit */
/*        system. On machines like Cray, we are reducing EMAX by one */
/*        unnecessarily. */

	--(*emax);
    }

    if (*ieee) {

/*        Assume we are on an IEEE machine which reserves one exponent */
/*        for infinity and NaN. */

	--(*emax);
    }

/*     Now create RMAX, the largest machine number, which should */
/*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX . */

/*     First compute 1.0 - BETA**(-P), being careful that the */
/*     result is less than 1.0 . */

    recbas = 1.f / *beta;
    z__ = *beta - 1.f;
    y = 0.f;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z__ *= recbas;
	if (y < 1.f) {
	    oldy = y;
	}
	y = slamc3_(&y, &z__);
/* L20: */
    }
    if (y >= 1.f) {
	y = oldy;
    }

/*     Now multiply by BETA**EMAX to get RMAX. */

    i__1 = *emax;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r__1 = y * *beta;
	y = slamc3_(&r__1, &c_b32);
/* L30: */
    }

    *rmax = y;
    return 0;

/*     End of SLAMC5 */

} /* slamc5_ */
