/* hybrj.f -- translated by f2c (version 20020621).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "minpack.h"
#include <math.h>
#define real __minpack_real__

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE_ (1)
#define FALSE_ (0)

__minpack_attr__
void __minpack_func__(hybrj)(__minpack_decl_fcnder_nn__  const int *n, real *x, real *
	fvec, real *fjac, const int *ldfjac, const real *xtol, const int *
	maxfev, real *diag, const int *mode, const real *factor, const int *
	nprint, int *info, int *nfev, int *njev, real *r__, 
	const int *lr, real *qtf, real *wa1, real *wa2, 
	real *wa3, real *wa4)
{
    /* Initialized data */

#define p1 .1
#define p5 .5
#define p001 .001
#define p0001 1e-4
    const int c_false = FALSE_;
    const int c__1 = 1;

    /* System generated locals */
    int fjac_dim1, fjac_offset, i__1, i__2;
    real d__1, d__2;

    /* Local variables */
    int i__, j, l, jm1, iwa[1];
    real sum;
    int sing;
    int iter;
    real temp;
    int iflag;
    real delta;
    int jeval;
    int ncsuc;
    real ratio;
    real fnorm;
    real pnorm, xnorm, fnorm1;
    int nslow1, nslow2;
    int ncfail;
    real actred, epsmch, prered;

/*     ********** */

/*     subroutine hybrj */

/*     the purpose of hybrj is to find a zero of a system of */
/*     n nonlinear functions in n variables by a modification */
/*     of the powell hybrid method. the user must provide a */
/*     subroutine which calculates the functions and the jacobian. */

/*     the subroutine statement is */

/*       subroutine hybrj(fcn,n,x,fvec,fjac,ldfjac,xtol,maxfev,diag, */
/*                        mode,factor,nprint,info,nfev,njev,r,lr,qtf, */
/*                        wa1,wa2,wa3,wa4) */

/*     where */

/*       fcn is the name of the user-supplied subroutine which */
/*         calculates the functions and the jacobian. fcn must */
/*         be declared in an external statement in the user */
/*         calling program, and should be written as follows. */

/*         subroutine fcn(n,x,fvec,fjac,ldfjac,iflag) */
/*         integer n,ldfjac,iflag */
/*         double precision x(n),fvec(n),fjac(ldfjac,n) */
/*         ---------- */
/*         if iflag = 1 calculate the functions at x and */
/*         return this vector in fvec. do not alter fjac. */
/*         if iflag = 2 calculate the jacobian at x and */
/*         return this matrix in fjac. do not alter fvec. */
/*         --------- */
/*         return */
/*         end */

/*         the value of iflag should not be changed by fcn unless */
/*         the user wants to terminate execution of hybrj. */
/*         in this case set iflag to a negative integer. */

/*       n is a positive integer input variable set to the number */
/*         of functions and variables. */

/*       x is an array of length n. on input x must contain */
/*         an initial estimate of the solution vector. on output x */
/*         contains the final estimate of the solution vector. */

/*       fvec is an output array of length n which contains */
/*         the functions evaluated at the output x. */

/*       fjac is an output n by n array which contains the */
/*         orthogonal matrix q produced by the qr factorization */
/*         of the final approximate jacobian. */

/*       ldfjac is a positive integer input variable not less than n */
/*         which specifies the leading dimension of the array fjac. */

/*       xtol is a nonnegative input variable. termination */
/*         occurs when the relative error between two consecutive */
/*         iterates is at most xtol. */

/*       maxfev is a positive integer input variable. termination */
/*         occurs when the number of calls to fcn with iflag = 1 */
/*         has reached maxfev. */

/*       diag is an array of length n. if mode = 1 (see */
/*         below), diag is internally set. if mode = 2, diag */
/*         must contain positive entries that serve as */
/*         multiplicative scale factors for the variables. */

/*       mode is an integer input variable. if mode = 1, the */
/*         variables will be scaled internally. if mode = 2, */
/*         the scaling is specified by the input diag. other */
/*         values of mode are equivalent to mode = 1. */

/*       factor is a positive input variable used in determining the */
/*         initial step bound. this bound is set to the product of */
/*         factor and the euclidean norm of diag*x if nonzero, or else */
/*         to factor itself. in most cases factor should lie in the */
/*         interval (.1,100.). 100. is a generally recommended value. */

/*       nprint is an integer input variable that enables controlled */
/*         printing of iterates if it is positive. in this case, */
/*         fcn is called with iflag = 0 at the beginning of the first */
/*         iteration and every nprint iterations thereafter and */
/*         immediately prior to return, with x and fvec available */
/*         for printing. fvec and fjac should not be altered. */
/*         if nprint is not positive, no special calls of fcn */
/*         with iflag = 0 are made. */

/*       info is an integer output variable. if the user has */
/*         terminated execution, info is set to the (negative) */
/*         value of iflag. see description of fcn. otherwise, */
/*         info is set as follows. */

/*         info = 0   improper input parameters. */

/*         info = 1   relative error between two consecutive iterates */
/*                    is at most xtol. */

/*         info = 2   number of calls to fcn with iflag = 1 has */
/*                    reached maxfev. */

/*         info = 3   xtol is too small. no further improvement in */
/*                    the approximate solution x is possible. */

/*         info = 4   iteration is not making good progress, as */
/*                    measured by the improvement from the last */
/*                    five jacobian evaluations. */

/*         info = 5   iteration is not making good progress, as */
/*                    measured by the improvement from the last */
/*                    ten iterations. */

/*       nfev is an integer output variable set to the number of */
/*         calls to fcn with iflag = 1. */

/*       njev is an integer output variable set to the number of */
/*         calls to fcn with iflag = 2. */

/*       r is an output array of length lr which contains the */
/*         upper triangular matrix produced by the qr factorization */
/*         of the final approximate jacobian, stored rowwise. */

/*       lr is a positive integer input variable not less than */
/*         (n*(n+1))/2. */

/*       qtf is an output array of length n which contains */
/*         the vector (q transpose)*fvec. */

/*       wa1, wa2, wa3, and wa4 are work arrays of length n. */

/*     subprograms called */

/*       user-supplied ...... fcn */

/*       minpack-supplied ... dogleg,dpmpar,enorm, */
/*                            qform,qrfac,r1mpyq,r1updt */

/*       fortran-supplied ... dabs,dmax1,dmin1,mod */

/*     argonne national laboratory. minpack project. march 1980. */
/*     burton s. garbow, kenneth e. hillstrom, jorge j. more */

/*     ********** */
    /* Parameter adjustments */
    --wa4;
    --wa3;
    --wa2;
    --wa1;
    --qtf;
    --diag;
    --fvec;
    --x;
    fjac_dim1 = *ldfjac;
    fjac_offset = 1 + fjac_dim1 * 1;
    fjac -= fjac_offset;
    --r__;

    /* Function Body */

/*     epsmch is the machine precision. */

    epsmch = __minpack_func__(dpmpar)(&c__1);

    *info = 0;
    iflag = 0;
    *nfev = 0;
    *njev = 0;

/*     check the input parameters for errors. */

    if (*n <= 0 || *ldfjac < *n || *xtol < 0. || *maxfev <= 0 || *factor <= 
	    0. || *lr < *n * (*n + 1) / 2) {
	goto L300;
    }
    if (*mode != 2) {
	goto L20;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (diag[j] <= 0.) {
	    goto L300;
	}
/* L10: */
    }
L20:

/*     evaluate the function at the starting point */
/*     and calculate its norm. */

    iflag = 1;
    fcnder_nn(n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &iflag);
    *nfev = 1;
    if (iflag < 0) {
	goto L300;
    }
    fnorm = __minpack_func__(enorm)(n, &fvec[1]);

/*     initialize iteration counter and monitors. */

    iter = 1;
    ncsuc = 0;
    ncfail = 0;
    nslow1 = 0;
    nslow2 = 0;

/*     beginning of the outer loop. */

L30:
    jeval = TRUE_;

/*        calculate the jacobian matrix. */

    iflag = 2;
    fcnder_nn(n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &iflag);
    ++(*njev);
    if (iflag < 0) {
	goto L300;
    }

/*        compute the qr factorization of the jacobian. */

    __minpack_func__(qrfac)(n, n, &fjac[fjac_offset], ldfjac, &c_false, iwa, &c__1, &wa1[1], &
	    wa2[1], &wa3[1]);

/*        on the first iteration and if mode is 1, scale according */
/*        to the norms of the columns of the initial jacobian. */

    if (iter != 1) {
	goto L70;
    }
    if (*mode == 2) {
	goto L50;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	diag[j] = wa2[j];
	if (wa2[j] == 0.) {
	    diag[j] = 1.;
	}
/* L40: */
    }
L50:

/*        on the first iteration, calculate the norm of the scaled x */
/*        and initialize the step bound delta. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	wa3[j] = diag[j] * x[j];
/* L60: */
    }
    xnorm = __minpack_func__(enorm)(n, &wa3[1]);
    delta = *factor * xnorm;
    if (delta == 0.) {
	delta = *factor;
    }
L70:

/*        form (q transpose)*fvec and store in qtf. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	qtf[i__] = fvec[i__];
/* L80: */
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (fjac[j + j * fjac_dim1] == 0.) {
	    goto L110;
	}
	sum = 0.;
	i__2 = *n;
	for (i__ = j; i__ <= i__2; ++i__) {
	    sum += fjac[i__ + j * fjac_dim1] * qtf[i__];
/* L90: */
	}
	temp = -sum / fjac[j + j * fjac_dim1];
	i__2 = *n;
	for (i__ = j; i__ <= i__2; ++i__) {
	    qtf[i__] += fjac[i__ + j * fjac_dim1] * temp;
/* L100: */
	}
L110:
/* L120: */
	;
    }

/*        copy the triangular factor of the qr factorization into r. */

    sing = FALSE_;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	l = j;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L140;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    r__[l] = fjac[i__ + j * fjac_dim1];
	    l = l + *n - i__;
/* L130: */
	}
L140:
	r__[l] = wa1[j];
	if (wa1[j] == 0.) {
	    sing = TRUE_;
	}
/* L150: */
    }

/*        accumulate the orthogonal factor in fjac. */

    __minpack_func__(qform)(n, n, &fjac[fjac_offset], ldfjac, &wa1[1]);

/*        rescale if necessary. */

    if (*mode == 2) {
	goto L170;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__1 = diag[j], d__2 = wa2[j];
	diag[j] = max(d__1,d__2);
/* L160: */
    }
L170:

/*        beginning of the inner loop. */

L180:

/*           if requested, call fcn to enable printing of iterates. */

    if (*nprint <= 0) {
	goto L190;
    }
    iflag = 0;
    if ((iter - 1) % *nprint == 0) {
	fcnder_nn(n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &iflag);
    }
    if (iflag < 0) {
	goto L300;
    }
L190:

/*           determine the direction p. */

    __minpack_func__(dogleg)(n, &r__[1], lr, &diag[1], &qtf[1], &delta, &wa1[1], &wa2[1], &wa3[
	    1]);

/*           store the direction p and x + p. calculate the norm of p. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	wa1[j] = -wa1[j];
	wa2[j] = x[j] + wa1[j];
	wa3[j] = diag[j] * wa1[j];
/* L200: */
    }
    pnorm = __minpack_func__(enorm)(n, &wa3[1]);

/*           on the first iteration, adjust the initial step bound. */

    if (iter == 1) {
	delta = min(delta,pnorm);
    }

/*           evaluate the function at x + p and calculate its norm. */

    iflag = 1;
    fcnder_nn(n, &wa2[1], &wa4[1], &fjac[fjac_offset], ldfjac, &iflag);
    ++(*nfev);
    if (iflag < 0) {
	goto L300;
    }
    fnorm1 = __minpack_func__(enorm)(n, &wa4[1]);

/*           compute the scaled actual reduction. */

    actred = -1.;
    if (fnorm1 < fnorm) {
/* Computing 2nd power */
	d__1 = fnorm1 / fnorm;
	actred = 1. - d__1 * d__1;
    }

/*           compute the scaled predicted reduction. */

    l = 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
	i__2 = *n;
	for (j = i__; j <= i__2; ++j) {
	    sum += r__[l] * wa1[j];
	    ++l;
/* L210: */
	}
	wa3[i__] = qtf[i__] + sum;
/* L220: */
    }
    temp = __minpack_func__(enorm)(n, &wa3[1]);
    prered = 0.;
    if (temp < fnorm) {
/* Computing 2nd power */
	d__1 = temp / fnorm;
	prered = 1. - d__1 * d__1;
    }

/*           compute the ratio of the actual to the predicted */
/*           reduction. */

    ratio = 0.;
    if (prered > 0.) {
	ratio = actred / prered;
    }

/*           update the step bound. */

    if (ratio >= p1) {
	goto L230;
    }
    ncsuc = 0;
    ++ncfail;
    delta = p5 * delta;
    goto L240;
L230:
    ncfail = 0;
    ++ncsuc;
    if (ratio >= p5 || ncsuc > 1) {
/* Computing MAX */
	d__1 = delta, d__2 = pnorm / p5;
	delta = max(d__1,d__2);
    }
    if (fabs(ratio - 1.) <= p1) {
	delta = pnorm / p5;
    }
L240:

/*           test for successful iteration. */

    if (ratio < p0001) {
	goto L260;
    }

/*           successful iteration. update x, fvec, and their norms. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	x[j] = wa2[j];
	wa2[j] = diag[j] * x[j];
	fvec[j] = wa4[j];
/* L250: */
    }
    xnorm = __minpack_func__(enorm)(n, &wa2[1]);
    fnorm = fnorm1;
    ++iter;
L260:

/*           determine the progress of the iteration. */

    ++nslow1;
    if (actred >= p001) {
	nslow1 = 0;
    }
    if (jeval) {
	++nslow2;
    }
    if (actred >= p1) {
	nslow2 = 0;
    }

/*           test for convergence. */

    if (delta <= *xtol * xnorm || fnorm == 0.) {
	*info = 1;
    }
    if (*info != 0) {
	goto L300;
    }

/*           tests for termination and stringent tolerances. */

    if (*nfev >= *maxfev) {
	*info = 2;
    }
/* Computing MAX */
    d__1 = p1 * delta;
    if (p1 * max(d__1,pnorm) <= epsmch * xnorm) {
	*info = 3;
    }
    if (nslow2 == 5) {
	*info = 4;
    }
    if (nslow1 == 10) {
	*info = 5;
    }
    if (*info != 0) {
	goto L300;
    }

/*           criterion for recalculating jacobian. */

    if (ncfail == 2) {
	goto L290;
    }

/*           calculate the rank one modification to the jacobian */
/*           and update qtf if necessary. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	sum = 0.;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sum += fjac[i__ + j * fjac_dim1] * wa4[i__];
/* L270: */
	}
	wa2[j] = (sum - wa3[j]) / pnorm;
	wa1[j] = diag[j] * (diag[j] * wa1[j] / pnorm);
	if (ratio >= p0001) {
	    qtf[j] = sum;
	}
/* L280: */
    }

/*           compute the qr factorization of the updated jacobian. */

    __minpack_func__(r1updt)(n, n, &r__[1], lr, &wa1[1], &wa2[1], &wa3[1], &sing);
    __minpack_func__(r1mpyq)(n, n, &fjac[fjac_offset], ldfjac, &wa2[1], &wa3[1]);
    __minpack_func__(r1mpyq)(&c__1, n, &qtf[1], &c__1, &wa2[1], &wa3[1]);

/*           end of the inner loop. */

    jeval = FALSE_;
    goto L180;
L290:

/*        end of the outer loop. */

    goto L30;
L300:

/*     termination, either normal or user imposed. */

    if (iflag < 0) {
	*info = iflag;
    }
    iflag = 0;
    if (*nprint > 0) {
	fcnder_nn(n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &iflag);
    }
    return;

/*     last card of subroutine hybrj. */

} /* hybrj_ */

