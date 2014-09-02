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

#include "gromacs/math/irreg_array_math/irreg_array_2d_lu_decomposition.h"
#include "gromacs/math/irreg_array_math/irreg_array_2d_matrix_products.h"


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
    gmx::IrregArray2D<real> T_i(int(0), n, int(0), 0 /*k*/);
    // With : W_i = (W_1, ... , W_i)
    gmx::IrregArray2D<real> W_i(int(0), n, int(0), 0 /*k*/);

    // w_i = X^T × y
    std::vector<real> w_i(n);
    gmx::matrixProductMatrixTVector(X, y, w_i);

    // t_1 = X × w_1
    std::vector<real> t_i(n);
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
    for (int i = 2; i < k; i++)
    {
        // w_i = X^T × r
        gmx::matrixProductMatrixTVector(X, r, w_i);

        // t_i = X × w_i
        gmx::matrixProductMatrixVector(X, w_i, t_i);

        // Keep t_i and w_i in T_i and W_i
        T_i.appendVector(t_i.data(), n);
        W_i.appendVector(w_i.data(), n);

        // r = y - T_i × (T_i^T × T_i)^-1 × T_i^T × y
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
    {
        // For i == k
        // w_i = X^T × r
        gmx::matrixProductMatrixTVector(X, r, w_i);

        // t_i = X × w_i
        gmx::matrixProductMatrixVector(X, w_i, t_i);

        // Keep t_i and w_i in T_i and W_i
        T_i.appendVector(t_i.data(), n);
        W_i.appendVector(w_i.data(), n);
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

        (*q) = std::move(T_iT_T_i_INV_T_iT_y);
    }

    // b = W_k × q
    std::vector<real> b;
    gmx::matrixProductMatrixVector(W_i, (*q), b);

    // Keep the results
    (*W) = std::move(W_i);

    return 0;
}
