/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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
 * Implements partial * least squares regression
 *
 * \author Jan Henning Peters <JanHPeters@gmx.net>
 */
#include "gmxpre.h"

#include "partial_least_squares.h"

#include <vector.h>

#include "gmx_blas.h"
#include "gmx_lapack.h"


int
pls_denham(FMatrix *x, FVector *y, int n, int a, FMatrix *w, FVector *q)
{

    int k  = x->NCols();
    int ln = x->NRows();


    FMatrix t(n, a);
    FMatrix qrt(n, a);
    FVector b(k);
    FVector dum(n);
    FVector rsd(n);

    int     i;
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
    GMX_LAPACK(gemv, GEMV) (    &T, &n, &k, &c_done, x->toF(), &ln, y->toF(), &c_ione, &c_dzero, w->toF(), &c_ione);

    /* t_1 = x * w_1 */
    GMX_LAPACK(gemv, GEMV) (&N, &n, &k, &c_done, x->toF(), &ln, w->toF(), &c_ione, &c_dzero, t.toF(), &c_ione);

    GMX_LAPACK(copy, COPY) (&n, t.toF(), &c_ione, qrt.toF(), &c_ione);
    GMX_LAPACK(copy, COPY) (&n, y->toF(), &c_ione, rsd.toF(), &c_ione);

    /* perform regression
       min(rsd) abs(y - qrt * rsd) */
    GMX_LAPACK(gels, GELS) (    &N,      &n, &c_ione, &c_ione, qrt.toF(),  &n,  rsd.toF(),  &n,  dum.toF(), &lwork, &info );

    (*q)[0] = rsd[0];
    GMX_LAPACK(copy, COPY) (&n, y->toF(), &c_ione, rsd.toF(), &c_ione);

    /* calculate residuals rsd := y - t_1*q_1
       using DGEMV          (y := alpha*A*x + beta*y)*/
    GMX_LAPACK(gemv, GEMV) (    &N,      &n, &c_ione, &c_dnegone, t.toF(),      &n,     q->toF(), &c_ione, &c_done, rsd.toF(), &c_ione);


    for (i = 1; i < a; i++)
    {
        /* w_i := X' * rsd */
        GMX_LAPACK(gemv, GEMV) (&T, &n, &k, &c_done, x->toF(), &ln, rsd.toF(), &c_ione, &c_dzero, &(w->toF())[i*k], &c_ione);
        /* t_i := X * w_i */
        GMX_LAPACK(gemv, GEMV) (&N, &n, &k, &c_done, x->toF(), &ln, &(w->toF())[i*k], &c_ione, &c_dzero, &(t.toF())[n*i], &c_ione);

        idum = n*a;
        GMX_LAPACK(copy, COPY) (&idum, t.toF(), &c_ione, qrt.toF(), &c_ione);
        GMX_LAPACK(copy, COPY) (&n, y->toF(), &c_ione, rsd.toF(), &c_ione);

        idum = i + 1;
        /* run regression */
        GMX_LAPACK(gels, GELS) ( &N, &n, &idum, &c_ione, qrt.toF(),  &n, rsd.toF(),  &n,  dum.toF(), &lwork, &info );

        GMX_LAPACK(copy, COPY) (&a, rsd.toF(), &c_ione,  q->toF(), &c_ione);

        GMX_LAPACK(copy, COPY) (&n, y->toF(), &c_ione, rsd.toF(), &c_ione);

        /* calculate residuals rsd := y - t_1*q_1
           using DGEMV          (y := alpha*A*x + beta*y) */
        idum = i + 1;
        GMX_LAPACK(gemv, GEMV) (    &N,      &n, &idum, &c_dnegone, t.toF(),      &n, q->toF(), &c_ione, &c_done, rsd.toF(), &c_ione);

    }

    GMX_LAPACK(gemv, GEMV) (  &N, &k, &a, &c_done, w->toF(), &k, q->toF(), &c_ione, &c_dzero, b.toF(), &c_ione);
    return 0;
}
