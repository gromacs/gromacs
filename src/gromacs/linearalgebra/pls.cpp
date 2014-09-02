/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements gmx::linearalgebra::Pls, a wrapper class for partial
 * least squares regression
 *
 * \author Jan Henning Peters <JanHPeters@gmx.net>
 * \ingroup linearalgebra
 */

#include <vector>

#include "pls.h"
#include "gmx_blas.h"
#include "gmx_lapack.h"
#include "gromacs/utility/smalloc.h"

/*! \brief Partial least squares (PLS) regression
 *
 * Implemented as described in:
 *
 * Denham, M. C. "Implementing Partial Least Squares."
 * Statistics and Computing 5, no. 3 (September 1995): 191–202. doi:10.1007/BF00142661.
 *
 * The algoritm aims to find a vector t to minimize
 *
 * abs(y-X*t)
 *
 * where the vector t=w*q is a linear combination of vectors.
 *
 * \param x  contains the centered (ln,k)-matrix X.
 * \param y  contains the centered (n)-vector y.
 * \param n  the number of "relevant" rows of the matrix X.
 * \param ln the "real" number of rows in the matrix x, or the "leading dimension"
 *           in FORTRAN terms. If all structure-value pairs are used for training, this number
 *           is equal to n, otherwise bigger.
 * \param k  the number of columns (to be used) of the matrix x.
 * \param a  < MIN(n-1,k) - the number of PLS factors to include in the regression
 *           of y on the matrix X.
 * \param w  (output) a (k,a)-matrix containing the coefficient vectors stored by column of
 *           the a PLS regression factors obtained from the matrix X on .
 * \param q  a (a)-vector containing the least squares regression coefficients of
 *           the ordinary least squares regression of y on the PLS factor matrix t.
 *
 * Arrays entered into this function should be in FORTRAN-compatible
 * column major order. This is
 */
int
pls_denham(FMatrix *x, RVector *y, int n, int k, int ln, int a, FMatrix *w, RVector *q)
{

    /*! \param x  contains the centered (ln,k)-matrix X.
        \param y  contains the centered (n)-vector y.
        \param n  the number of "relevant" rows of the matrix X.
        \param ln the "real" number of rows in the matrix x, or the "leading dimension"
                     in FORTRAN terms. If all structure-value pairs are used for training, this number
                  is equal to n, otherwise bigger.
        \param k  the number of columns (to be used) of the matrix x.
        \param a  < MIN(n-1,k) - the number of PLS factors to include in the regression
                  of y on the matrix X.
        \param w  (output) a (k,a)-matrix containing the coefficient vectors stored by column of
                  the a PLS regression factors obtained from the matrix X on .
        \param q  a (a)-vector containing the least squares regression coefficients of
                  the ordinary least squares regression of y on the PLS factor matrix t.

       Arrays entered into this function should be in FORTRAN-compatible
       column major order.
     */

	FMatrix t(n,a);
	FMatrix qrt(n,a);
	RVector b(k);
	RVector dum(n);
	RVector rsd(n);

    int i;
    /* As FORTRAN functions in C require call by reference, all constants have
       to be defined explicitly*/
    char T         = 'T';
    char N         = 'N';
    real c_done    = 1.0;
    real c_dnegone = -1.0;
    real c_dzero   = 0.0;
    int  c_ione    = 1;

    int  idum;
    int  info  = 0;
    int  lwork = -1;

    /* Note: as BLAS and LAPACK functions are less than self-explaining, this function
       contains comments that explain what the calls do. The algorithm itself is described
       in the paper cited above. */
    /*   w_1 = x' * y*/
    GMX_LAPACK(gemv, GEMV) (    &T, &n, &k, &c_done, PTOF(x), &ln, PTOF(y), &c_ione, &c_dzero, PTOF(w), &c_ione);

    /* t_1 = x * w_1 */
    GMX_LAPACK(gemv, GEMV) (&N, &n, &k, &c_done, PTOF(x), &ln, PTOF(w), &c_ione, &c_dzero, TOF(t), &c_ione);

    GMX_LAPACK(copy, COPY) (&n, TOF(t), &c_ione, TOF(qrt), &c_ione);
    GMX_LAPACK(copy, COPY) (&n, PTOF(y), &c_ione, TOF(rsd), &c_ione);

    /* perform regression
       min(rsd) abs(y - qrt * rsd) */
    GMX_LAPACK(gels, GELS) (    &N,      &n, &c_ione, &c_ione, TOF(qrt),  &n,  TOF(rsd),  &n,  TOF(dum), &lwork, &info );

    (*q)[0] = rsd[0];
    GMX_LAPACK(copy, COPY) (&n, PTOF(y), &c_ione, TOF(rsd), &c_ione);

    /* calculate residuals rsd := y - t_1*q_1
       using DGEMV          (y := alpha*A*x + beta*y)*/
    GMX_LAPACK(gemv, GEMV) (    &N,      &n, &c_ione, &c_dnegone, TOF(t),      &n,     PTOF(q), &c_ione, &c_done, TOF(rsd), &c_ione);


    for (i = 1; i < a; i++)
    {
        /* w_i := X' * rsd */
        GMX_LAPACK(gemv, GEMV) (&T, &n, &k, &c_done, PTOF(x), &ln, TOF(rsd), &c_ione, &c_dzero, &PTOF(w)[i*k], &c_ione);
        /* t_i := X * w_i */
        GMX_LAPACK(gemv, GEMV) (&N, &n, &k, &c_done, PTOF(x), &ln, &PTOF(w)[i*k], &c_ione, &c_dzero, &t[n*i], &c_ione);

        idum = n*a;
        GMX_LAPACK(copy, COPY) (&idum, TOF(t), &c_ione, TOF(qrt), &c_ione);
        GMX_LAPACK(copy, COPY) (&n, PTOF(y), &c_ione, TOF(rsd), &c_ione);

        idum = i + 1;
        /* run regression */
        GMX_LAPACK(gels, GELS) ( &N, &n, &idum, &c_ione, TOF(qrt),  &n, TOF(rsd),  &n,  TOF(dum), &lwork, &info );

        GMX_LAPACK(copy, COPY) (&a, TOF(rsd), &c_ione,  PTOF(q), &c_ione);

        GMX_LAPACK(copy, COPY) (&n, PTOF(y), &c_ione, TOF(rsd), &c_ione);

        /* calculate residuals rsd := y - t_1*q_1
           using DGEMV          (y := alpha*A*x + beta*y) */
        idum = i + 1;
        GMX_LAPACK(gemv, GEMV) (    &N,      &n, &idum, &c_dnegone, TOF(t),      &n, PTOF(q), &c_ione, &c_done, TOF(rsd), &c_ione);

    }

    GMX_LAPACK(gemv, GEMV) (  &N, &k, &a, &c_done, PTOF(w), &k, PTOF(q), &c_ione, &c_dzero, TOF(b), &c_ione);
    return 0;
}

