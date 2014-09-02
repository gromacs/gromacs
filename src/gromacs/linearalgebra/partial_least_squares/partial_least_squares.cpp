/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2017,2018, by the GROMACS development team, led by
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

#include "gromacs/math/irreg_array_math.h"

/*!
 * \brief
 * This function copy the column cl into the column of index columnIdx of mt
 */
template <class MatrixType, class ColumnType>
void setColumn(MatrixType &mt, const int columnIdx, const ColumnType &cl, const int nbRows)
{
    if (int(mt.getLength1()) != nbRows)
    {
        std::fprintf(stderr, "Error in pls_denham setColumn invalid nbRows, mt.getLength1() %zu != nbRows %d",
                     mt.getLength1(), nbRows);
        throw;
    }
    if (int(mt.getLength2()) <= columnIdx)
    {
        std::fprintf(stderr, "Error in pls_denham setColumn invalid column, mt.getLength2() %zu != columnIdx %d",
                     mt.getLength2(), columnIdx);
        throw;
    }
    if (int(cl.size()) != nbRows)
    {
        std::fprintf(stderr, "Error in pls_denham setColumn invalid size, cl.size() %zu != nbRows %d",
                     cl.size(), nbRows);
        throw;
    }

    for (int idxRow = 0; idxRow < nbRows; ++idxRow)
    {
        mt(idxRow)[columnIdx] = cl[idxRow];
    }
}


class IrregArray2DColumnView
{
    const Pls2dArray &matrix;
    const int         columnIdx;
    public:
        IrregArray2DColumnView(const Pls2dArray              &inMatrix,
                               const int                      inColumnIdx)
            : matrix(inMatrix), columnIdx(inColumnIdx)
        {
            if (int(matrix.getLength2()) <= columnIdx)
            {
                std::fprintf(stderr, "Error in pls_denham IrregArray2DColumnView invalid columnIdx, matrix.getLength2() %zu != columnIdx %d",
                             matrix.getLength2(), columnIdx);
                throw;
            }
        }

        const real &operator[](const int inRowIdx) const
        {
            if (int(matrix.getLength1()) <= inRowIdx)
            {
                std::fprintf(stderr, "Error in pls_denham IrregArray2DColumnView invalid inRowIdx, matrix.getLength1() %zu != inRowIdx %d",
                             matrix.getLength1(), inRowIdx);
                throw;
            }
            return matrix(inRowIdx)[columnIdx];
        }

        size_t size() const
        {
            return matrix.getLength1();
        }
};


int pls_denham(const Pls2dArray &X,
               const Pls1dArray &y,
               const int n, const int p, const int k,
               Pls2dArray &W,
               Pls1dArray &q)
{
    using allocator_type = gmx::SimdSetup<real>::allocator_type;
    using index_type     = Pls2dArray::index_type;

    if (p < k)
    {
        throw std::invalid_argument("Error in pls_denham invalid k.");
    }
    if (int(q.size()) != k)
    {
        throw std::invalid_argument("Error in pls_denham invalid q.");
    }
    if (int(W.getLength1()) != p || int(W.getLength2()) != k)
    {
        throw std::invalid_argument("Error in pls_denham invalid W.");
    }
    if (int(y.size()) != n)
    {
        throw std::invalid_argument("Error in pls_denham invalid y.");
    }
    if (int(X.getLength1()) != n || int(X.getLength2()) != p)
    {
        throw std::invalid_argument("Error in pls_denham invalid y.");
    }


    // With : T_k = (t_1, ... , t_k)
    gmx::IrregArray2D<real, allocator_type> T_k(index_type(0), index_type(n-1), index_type(0), index_type(k-1));
    // With : W_k = (W_1, ... , W_k)
    gmx::IrregArray2D<real, allocator_type> W_k(index_type(0), index_type(p-1), index_type(0), index_type(k-1));

    // w_i = X^T × y
    std::vector<real, allocator_type> w_i(p);
    gmx::matrixProductMatrixTVector(X, y, w_i);

    // t_1 = X × w_1
    std::vector<real, allocator_type> t_i(n);
    gmx::matrixProductMatrixVector(X, w_i, t_i);

    // Keep t_i and w_i in T_k and W_k
    setColumn(T_k, 0, t_i, n);
    setColumn(W_k, 0, w_i, p);

    // r = y - t_1 × (t_1^T × t_1)^-1 × t_1^T × y
    std::vector<real, allocator_type> r(n);
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
        setColumn(T_k, i-1, t_i, n);
        setColumn(W_k, i-1, w_i, p);

        // With : T_k = (t_1, ... , t_i)
        gmx::IrregArray2D<real, allocator_type> T_i(index_type(0), index_type(n-1), index_type(0), index_type(i-1));
        for (int i_copy = 0; i_copy < i; ++i_copy)
        {
            setColumn(T_i, i_copy, IrregArray2DColumnView(T_k, i_copy), n);
        }


        // r = y - T_i × (T_i^T × T_i)^-1 × T_i^T × y
        {
            gmx::IrregArray2D<real, allocator_type> T_iT_T_i(index_type(0), index_type(i-1), index_type(0), index_type(i-1));
            gmx::matrixProductMatrixTMatrix(T_i, T_i, T_iT_T_i);

            gmx::IrregArray2D<real, allocator_type> T_iT_T_i_INV(index_type(0), index_type(i-1), index_type(0), index_type(i-1));
            gmx::matrixInverseLUP(T_iT_T_i, T_iT_T_i_INV);

            std::vector<real, allocator_type> T_iT_y(i);
            gmx::matrixProductMatrixTVector(T_i, y, T_iT_y);

            std::vector<real, allocator_type> T_iT_T_i_INV_T_iT_y(i);
            gmx::matrixProductMatrixVector(T_iT_T_i_INV, T_iT_y, T_iT_T_i_INV_T_iT_y);

            std::vector<real, allocator_type> T_i_T_iT_T_i_INV_T_iT_y(n);
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
        setColumn(T_k, k-1, t_i, n);
        setColumn(W_k, k-1, w_i, p);
    }

    // q = (T_k^T × T_k)^-1 × T_k^T × y
    {
        const int                               i = k;
        gmx::IrregArray2D<real, allocator_type> T_kT_T_k(index_type(0), index_type(i-1), index_type(0), index_type(i-1));
        gmx::matrixProductMatrixTMatrix(T_k, T_k, T_kT_T_k);

        gmx::IrregArray2D<real, allocator_type> T_kT_T_k_INV(index_type(0), index_type(i-1), index_type(0), index_type(i-1));
        gmx::matrixInverseLUP(T_kT_T_k, T_kT_T_k_INV);

        std::vector<real, allocator_type> T_kT_y(i);
        gmx::matrixProductMatrixTVector(T_k, y, T_kT_y);

        std::vector<real, allocator_type> T_kT_T_k_INV_T_kT_y(i);
        gmx::matrixProductMatrixVector(T_kT_T_k_INV, T_kT_y, T_kT_T_k_INV_T_kT_y);

        q = std::move(T_kT_T_k_INV_T_kT_y);
    }

    // currently not needed: b = W_k × q
    // std::vector<real, allocator_type> b(p);
    // gmx::matrixProductMatrixVector(W_k, q, b);

    // Keep the results
    for (int i = 0; i < p; ++i)
    {
        std::copy_n(W_k(i), k, W(i));
    }

    return 0;
}
