/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2017, by the GROMACS development team, led by
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
#include "gromacs/math/irreg_array_math/irreg_array_2d_matrix_products.h"
#include "gromacs/math/irreg_array_math/irreg_array_2d_lu_decomposition.h"


template <class real_type>
real_type scalar_product(const std::vector<real_type> &inVec1,
                         const std::vector<real_type> &inVec2)
{
    if (inVec1.size() != inVec2.size())
    {
        std::fprintf(stderr, "Error in pls_denham::scalar_product(): size of the vectors does not match.");
        throw;
    }

    real_type product = 0;
    for (size_t idx = 0; idx < inVec1.size(); ++idx)
    {
        product += inVec1[idx] * inVec2[idx];
    }
    return product;
}


int pls_denham(const gmx::IrregArray2D<real> &X, const std::vector<real> &y, const size_t n, const int k,
               gmx::IrregArray2D<real>* W, std::vector<real>* q)
{
    if (q == nullptr || q->size() != n)
    {
        std::fprintf(stderr, "Error in pls_denham invalid q.");
        throw;
    }
    if (W == nullptr || W->getLength1() != n || int(W->getLength2()) != k)
    {
        std::fprintf(stderr, "Error in pls_denham invalid W.");
        throw;
    }
    if (y.size() != n)
    {
        std::fprintf(stderr, "Error in pls_denham invalid y.");
        throw;
    }
    if (X.getLength1() != n || int(X.getLength2()) != k)
    {
        std::fprintf(stderr, "Error in pls_denham invalid y.");
        throw;
    }


    // With : T_i = (t_1, ... , t_i)
    gmx::IrregArray2D<real> T_i(int(0), n, int(0), k);
    // With : W_i = (W_1, ... , W_i)
    gmx::IrregArray2D<real> W_i(int(0), n, int(0), k);

    // w_i = X^T × y
    std::vector<real> w_i;
    gmx::matrixProductMatrixTVector(X, y, w_i);

    // t_1 = X × w_1
    std::vector<real> t_i;
    gmx::matrixProductMatrixVector(X, w_i, t_i);

    // Keep t_i and w_i in T_i and W_i
    T_i.appendVector(t_i.data(), n);
    W_i.appendVector(w_i.data(), n);

    // r = y - t_1 × (t_1^T × t_1)^-1 × t_1^T × y
    std::vector<real> r(n);
    {
        const real t_1T_t_1            = scalar_product(t_i, t_i);
        const real t_1T_t_1_INV        = 1 / t_1T_t_1;
        const real t_1T_y              = scalar_product(t_i, y);
        const real t_1T_t_1_INV_t_1T_y = t_1T_t_1_INV * t_1T_y;

        for (size_t idxVal = 0; idxVal < n; ++idxVal)
        {
            r[idxVal] = y[idxVal] - t_i[idxVal] * t_1T_t_1_INV_t_1T_y;
        }
    }

    // for i = 2 ... k
    for (int i = 2; i <= k; i++)
    {
        // w_i = X^T × r
        gmx::matrixProductMatrixTVector(X, r, w_i);

        // t_i = X × w_i
        gmx::matrixProductMatrixVector(X, w_i, t_i);

        // Keep t_i and w_i in T_i and W_i
        T_i.appendVector(t_i.data(), n);
        W_i.appendVector(w_i.data(), n);

        // r = y - T_i × (T_i^T × T_i)^-1 × T_i^T × y
        if (i != k)  // r is no more used from i == k
        {
            gmx::IrregArray2D<real> T_iT_T_i(int(0), i, int(0), i);
            gmx::matrixProductMatrixTMatrix(T_i, T_i, T_iT_T_i);

            gmx::IrregArray2D<real> T_iT_T_i_INV(int(0), i, int(0), i);
            gmx::matrixInverseLUP(T_iT_T_i, T_iT_T_i_INV);

            std::vector<real> T_iT_y(i);
            gmx::matrixProductMatrixTVector(T_i, y, T_iT_y);

            std::vector<real> T_iT_T_i_INV_T_iT_y(i);
            gmx::matrixProductMatrixVector(T_iT_T_i_INV, T_iT_y, T_iT_T_i_INV_T_iT_y);

            std::vector<real> T_i_T_iT_T_i_INV_T_iT_y(n);
            gmx::matrixProductMatrixVector(T_i, T_iT_T_i_INV_T_iT_y, T_i_T_iT_T_i_INV_T_iT_y);

            for (size_t idxVal = 0; idxVal < n; ++idxVal)
            {
                r[idxVal] = y[idxVal] - T_i_T_iT_T_i_INV_T_iT_y[idxVal];
            }
        }
    }

    // q = (T_k^T × T_k)^-1 × T_k^T × y
    {
        const int               i = k;
        gmx::IrregArray2D<real> T_iT_T_i(int(0), i, int(0), i);
        gmx::matrixProductMatrixTMatrix(T_i, T_i, T_iT_T_i);

        gmx::IrregArray2D<real> T_iT_T_i_INV(int(0), i, int(0), i);
        gmx::matrixInverseLUP(T_iT_T_i, T_iT_T_i_INV);

        std::vector<real> T_iT_y(i);
        gmx::matrixProductMatrixTVector(T_i, y, T_iT_y);

        std::vector<real> T_iT_T_i_INV_T_iT_y(i);
        gmx::matrixProductMatrixVector(T_iT_T_i_INV, T_iT_y, T_iT_T_i_INV_T_iT_y);

        (*q) = T_iT_T_i_INV_T_iT_y;
    }

    // b = W_k × q
    std::vector<real> b;
    gmx::matrixProductMatrixVector(W_i, (*q), b);

    // Keep the results
    (*W) = std::move(W_i);

    /*ssize_t zero = 0;
       ssize_t k  = x->getLength2();// number of columns
       ssize_t ln = x->getLength1();// number of rows


       gmx::IrregArray2D<real> t(zero, n-1, zero, a-1);
       gmx::IrregArray2D<real> qrt(zero, n, zero, a-1);
       gmx::IrregArray1D<real> b(zero, k-1);
       gmx::IrregArray1D<real> dum(zero, n-1);
       gmx::IrregArray1D<real> rsd(zero, n-1);

       int     i;
       // As FORTRAN functions in C require call by reference, all constants have
       //   to be defined explicitly
       char T         = 'T';
       char N         = 'N';
       real c_done    = 1.0;
       real c_dnegone = -1.0;
       real c_dzero   = 0.0;
       int  c_ione    = 1;

       int  idum;
       int  info  = 0;
       int  lwork = -1;

       // Note: as BLAS and LAPACK functions are less than self-explaining, this function
       // contains comments that explain what the calls do. The algorithm itself is described
       // in the paper cited above.
       //   w_1 = x' * y
       GMX_LAPACK(gemv, GEMV) (    &T, &n, &k, &c_done, x->toF(), &ln, y->toF(), &c_ione, &c_dzero, w->toF(), &c_ione);

       // t_1 = x * w_1
       GMX_LAPACK(gemv, GEMV) (&N, &n, &k, &c_done, x->toF(), &ln, w->toF(), &c_ione, &c_dzero, t.toF(), &c_ione);

       GMX_LAPACK(copy, COPY) (&n, t.toF(), &c_ione, qrt.toF(), &c_ione);
       GMX_LAPACK(copy, COPY) (&n, y->toF(), &c_ione, rsd.toF(), &c_ione);

       // perform regression
       //   min_def(rsd) abs(y - qrt * rsd)
       GMX_LAPACK(gels, GELS) (    &N,      &n, &c_ione, &c_ione, qrt.toF(),  &n,  rsd.toF(),  &n,  dum.toF(), &lwork, &info );

       (*q)[0] = rsd[0];
       GMX_LAPACK(copy, COPY) (&n, y->toF(), &c_ione, rsd.toF(), &c_ione);

       // calculate residuals rsd := y - t_1*q_1
       // using DGEMV          (y := alpha*A*x + beta*y)
       GMX_LAPACK(gemv, GEMV) (    &N,      &n, &c_ione, &c_dnegone, t.toF(),      &n,     q->toF(), &c_ione, &c_done, rsd.toF(), &c_ione);


       for (i = 1; i < a; i++)
       {
        // w_i := X' * rsd
        GMX_LAPACK(gemv, GEMV) (&T, &n, &k, &c_done, x->toF(), &ln, rsd.toF(), &c_ione, &c_dzero, &(w->toF())[i*k], &c_ione);
        // t_i := X * w_i
        GMX_LAPACK(gemv, GEMV) (&N, &n, &k, &c_done, x->toF(), &ln, &(w->toF())[i*k], &c_ione, &c_dzero, &(t.toF())[n*i], &c_ione);

        idum = n*a;
        GMX_LAPACK(copy, COPY) (&idum, t.toF(), &c_ione, qrt.toF(), &c_ione);
        GMX_LAPACK(copy, COPY) (&n, y->toF(), &c_ione, rsd.toF(), &c_ione);

        idum = i + 1;
        // run regression
        GMX_LAPACK(gels, GELS) ( &N, &n, &idum, &c_ione, qrt.toF(),  &n, rsd.toF(),  &n,  dum.toF(), &lwork, &info );

        GMX_LAPACK(copy, COPY) (&a, rsd.toF(), &c_ione,  q->toF(), &c_ione);

        GMX_LAPACK(copy, COPY) (&n, y->toF(), &c_ione, rsd.toF(), &c_ione);

        // calculate residuals rsd := y - t_1*q_1
        // using DGEMV          (y := alpha*A*x + beta*y)
        idum = i + 1;
        GMX_LAPACK(gemv, GEMV) (    &N,      &n, &idum, &c_dnegone, t.toF(),      &n, q->toF(), &c_ione, &c_done, rsd.toF(), &c_ione);

       }

       GMX_LAPACK(gemv, GEMV) (  &N, &k, &a, &c_done, w->toF(), &k, q->toF(), &c_ione, &c_dzero, b.toF(), &c_ione);
     */
    return 0;
}
