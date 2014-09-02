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

template <class MatrixType, class ColumnType>
void setColumn(MatrixType &mt, const int columnIdx, const ColumnType &cl, const int length)
{
    for (int idxRow = 0; idxRow < length; ++idxRow)
    {
        mt(idxRow)[columnIdx] = cl[idxRow];
    }
}


int pls_denham(const gmx::IrregArray2D<real> &X, const std::vector<real> &y,
               const int n, const int p, const int k,
               gmx::IrregArray2D<real>* W, std::vector<real>* q)
{
    if (p < k)
    {
        std::fprintf(stderr, "Error in pls_denham invalid k.");
        throw;
    }
    if (q == nullptr || int(q->size()) != k)
    {
        std::fprintf(stderr, "Error in pls_denham invalid q.");
        throw;
    }
    if (W == nullptr || int(W->getLength1()) != p || int(W->getLength2()) != k)
    {
        std::fprintf(stderr, "Error in pls_denham invalid W.");
        throw;
    }
    if (int(y.size()) != n)
    {
        std::fprintf(stderr, "Error in pls_denham invalid y.");
        throw;
    }
    if (int(X.getLength1()) != n || int(X.getLength2()) != p)
    {
        std::fprintf(stderr, "Error in pls_denham invalid y.");
        throw;
    }


    // With : T_k = (t_1, ... , t_k)
    gmx::IrregArray2D<real> T_k(int(0), n-1, int(0), k-1);
    // With : W_k = (W_1, ... , W_k)
    gmx::IrregArray2D<real> W_k(int(0), p-1, int(0), k-1);

    // w_i = X^T × y
    std::vector<real> w_i(p);
    gmx::matrixProductMatrixTVector(X, y, w_i);

    // t_1 = X × w_1
    std::vector<real> t_i(n);
    gmx::matrixProductMatrixVector(X, w_i, t_i);

    // Keep t_i and w_i in T_k and W_k
    setColumn(T_k, 0, t_i.begin(), n);
    setColumn(W_k, 0, w_i.begin(), p);

    // r = y - t_1 × (t_1^T × t_1)^-1 × t_1^T × y
    std::vector<real> r(n);
    {
        const real t_1T_t_1            = std::inner_product(t_i.begin(), t_i.end(), t_i.begin(), 0);
        const real t_1T_t_1_INV        = 1 / t_1T_t_1;
        const real t_1T_y              = std::inner_product(t_i.begin(), t_i.end(), y.begin(), 0);
        const real t_1T_t_1_INV_t_1T_y = t_1T_t_1_INV * t_1T_y;

        for (int idxVal = 0; idxVal < n; ++idxVal)
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

        // Keep t_i and w_i in T_k and W_k
        setColumn(T_k, i-1, t_i.begin(), n);
        setColumn(W_k, i-1, w_i.begin(), p);

        // With : T_k = (t_1, ... , t_i)
        gmx::IrregArray2D<real> T_i(int(0), n-1, int(0), i-1);
        for (int i_copy = 0; i_copy < i; ++i_copy)
        {
            setColumn(T_i, i_copy, T_k(i_copy), n);
        }


        // r = y - T_i × (T_i^T × T_i)^-1 × T_i^T × y
        {
            gmx::IrregArray2D<real> T_iT_T_i(int(0), i-1, int(0), i-1);
            gmx::matrixProductMatrixTMatrix(T_i, T_i, T_iT_T_i);

            gmx::IrregArray2D<real> T_iT_T_i_INV(int(0), i-1, int(0), i-1);
            gmx::matrixInverseLUP(T_iT_T_i, T_iT_T_i_INV);

            std::vector<real> T_iT_y(i);
            gmx::matrixProductMatrixTVector(T_i, y, T_iT_y);

            std::vector<real> T_iT_T_i_INV_T_iT_y(i);
            gmx::matrixProductMatrixVector(T_iT_T_i_INV, T_iT_y, T_iT_T_i_INV_T_iT_y);

            std::vector<real> T_i_T_iT_T_i_INV_T_iT_y(n);
            gmx::matrixProductMatrixVector(T_i, T_iT_T_i_INV_T_iT_y, T_i_T_iT_T_i_INV_T_iT_y);

            for (int idxVal = 0; idxVal < n; ++idxVal)
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

        // Keep t_i and w_i in T_k and W_k
        setColumn(T_k, k-1, t_i.begin(), n);
        setColumn(W_k, k-1, w_i.begin(), k);
    }

    // q = (T_k^T × T_k)^-1 × T_k^T × y
    {
        const int               i = k;
        gmx::IrregArray2D<real> T_kT_T_k(int(0), i-1, int(0), i-1);
        gmx::matrixProductMatrixTMatrix(T_k, T_k, T_kT_T_k);

        gmx::IrregArray2D<real> T_kT_T_k_INV(int(0), i-1, int(0), i-1);
        gmx::matrixInverseLUP(T_kT_T_k, T_kT_T_k_INV);

        std::vector<real> T_kT_y(i);
        gmx::matrixProductMatrixTVector(T_k, y, T_kT_y);

        std::vector<real> T_kT_T_k_INV_T_kT_y(i);
        gmx::matrixProductMatrixVector(T_kT_T_k_INV, T_kT_y, T_kT_T_k_INV_T_kT_y);

        (*q) = std::move(T_kT_T_k_INV_T_kT_y);
    }

    // currently not needed: b = W_k × q
    // std::vector<real> b(p);
    // gmx::matrixProductMatrixVector(W_k, (*q), b);

    // Keep the results
    for (int i = 0; i < p; ++i)
    {
        std::copy_n(W_k(i), k, (*W)(i));
    }

    return 0;
}
