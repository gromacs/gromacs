/* hybrd1.f -- translated by f2c (version 20020621).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "cminpack.h"
#include "cminpackP.h"

__cminpack_attr__
int __cminpack_func__(hybrd1)(__cminpack_decl_fcn_nn__ void *p, int n, real *x, real *
	fvec, real tol, real *wa, int lwa)
{
    /* Initialized data */

    const real factor = 100.;

    /* Local variables */
    int j, ml, lr, mu, mode, nfev;
    real xtol;
    int index;
    real epsfcn;
    int maxfev, nprint;
    int info;

/*     ********** */

/*     subroutine hybrd1 */

/*     the purpose of hybrd1 is to find a zero of a system of */
/*     n nonlinear functions in n variables by a modification */
/*     of the powell hybrid method. this is done by using the */
/*     more general nonlinear equation solver hybrd. the user */
/*     must provide a subroutine which calculates the functions. */
/*     the jacobian is then calculated by a forward-difference */
/*     approximation. */

/*     the subroutine statement is */

/*       subroutine hybrd1(fcn,n,x,fvec,tol,info,wa,lwa) */

/*     where */

/*       fcn is the name of the user-supplied subroutine which */
/*         calculates the functions. fcn must be declared */
/*         in an external statement in the user calling */
/*         program, and should be written as follows. */

/*         subroutine fcn(n,x,fvec,iflag) */
/*         integer n,iflag */
/*         double precision x(n),fvec(n) */
/*         ---------- */
/*         calculate the functions at x and */
/*         return this vector in fvec. */
/*         --------- */
/*         return */
/*         end */

/*         the value of iflag should not be changed by fcn unless */
/*         the user wants to terminate execution of hybrd1. */
/*         in this case set iflag to a negative integer. */

/*       n is a positive integer input variable set to the number */
/*         of functions and variables. */

/*       x is an array of length n. on input x must contain */
/*         an initial estimate of the solution vector. on output x */
/*         contains the final estimate of the solution vector. */

/*       fvec is an output array of length n which contains */
/*         the functions evaluated at the output x. */

/*       tol is a nonnegative input variable. termination occurs */
/*         when the algorithm estimates that the relative error */
/*         between x and the solution is at most tol. */

/*       info is an integer output variable. if the user has */
/*         terminated execution, info is set to the (negative) */
/*         value of iflag. see description of fcn. otherwise, */
/*         info is set as follows. */

/*         info = 0   improper input parameters. */

/*         info = 1   algorithm estimates that the relative error */
/*                    between x and the solution is at most tol. */

/*         info = 2   number of calls to fcn has reached or exceeded */
/*                    200*(n+1). */

/*         info = 3   tol is too small. no further improvement in */
/*                    the approximate solution x is possible. */

/*         info = 4   iteration is not making good progress. */

/*       wa is a work array of length lwa. */

/*       lwa is a positive integer input variable not less than */
/*         (n*(3*n+13))/2. */

/*     subprograms called */

/*       user-supplied ...... fcn */

/*       minpack-supplied ... hybrd */

/*     argonne national laboratory. minpack project. march 1980. */
/*     burton s. garbow, kenneth e. hillstrom, jorge j. more */

/*     ********** */
    /* Parameter adjustments */
    --fvec;
    --x;
    --wa;

    /* Function Body */

/*     check the input parameters for errors. */

    if (n <= 0 || tol < 0. || lwa < n * (n * 3 + 13) / 2) {
        return 0;
    }

/*     call hybrd. */

    maxfev = (n + 1) * 200;
    xtol = tol;
    ml = n - 1;
    mu = n - 1;
    epsfcn = 0.;
    mode = 2;
    for (j = 1; j <= n; ++j) {
	wa[j] = 1.;
    }
    nprint = 0;
    lr = n * (n + 1) / 2;
    index = n * 6 + lr;
    info = __cminpack_func__(hybrd)(__cminpack_param_fcn_nn__ p, n, &x[1], &fvec[1], xtol, maxfev, ml, mu, epsfcn, &
	    wa[1], mode, factor, nprint, &nfev, &wa[index + 1], n, &
	    wa[n * 6 + 1], lr, &wa[n + 1], &wa[(n << 1) + 1], &wa[n * 3 
	    + 1], &wa[(n << 2) + 1], &wa[n * 5 + 1]);
    if (info == 5) {
	info = 4;
    }
    return info;

/*     last card of subroutine hybrd1. */

} /* hybrd1_ */

