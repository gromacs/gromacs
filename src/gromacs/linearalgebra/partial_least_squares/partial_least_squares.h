/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Implements partial least squares regression
 *
 * \author Jan Henning Peters <JanHPeters@gmx.net>
 * \inlibraryapi
 */

#ifndef GMX_LINEARALGEBRA_PLS_H
#define GMX_LINEARALGEBRA_PLS_H

#include "gmxpre.h"

#include <type_traits>

#include "gromacs/math/irreg_array_math.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/data_structures/irreg_array_2d.h"


/*! \libinternal
 * \brief Partial least squares (PLS) regression
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
 * The template types must provide different methods:
 * - MatrixType must provide a size_type attribute type,
 * a index_type attribute type, a getLength1 method,
 * a getLength2 method, an operator (i,j) to access value
 * at row i and column j.
 * - VectorType must provide a size method, a [i]
 * method to access value at position i.
 *
 * \inlibraryapi
 */
template <typename MatrixType, typename VectorType>
class PartialLeastSquares
{
    /** \brief The type to represent the size of matrices and vectors */
    using size_type = typename MatrixType::size_type;
    static_assert(std::is_same<typename MatrixType::size_type, typename VectorType::size_type>::value, "size_type must be the same for a matrix and a vector");
    /** \brief The type to represent the index to iterate on matrices and vectors */
    using index_type = typename MatrixType::index_type;

    /*!
     * \brief
     * This static method copy the column cl view into the column of index columnIdx of mt
     */
    template <class ColumnType>
    static void setColumn(MatrixType* mt, const size_type columnIdx, const ColumnType &cl, const size_type nbRows)
    {
        for (size_type idxRow = 0; idxRow < nbRows; ++idxRow)
        {
            (*mt)(idxRow, columnIdx) = cl[idxRow];
        }
    }

    /*! \libinternal
     * \brief
     * The class IrregArray2DColumnView provides a view over a column of a matrix,
     * that is it allow to access a column as if it was a vector.
     *
     * \inlibraryapi
     */
    class IrregArray2DColumnView
    {
        //! \brief the matrix to cover
        const MatrixType &matrix;
        //! \brief the column index to access
        const index_type  columnIdx;
        public:
            IrregArray2DColumnView(const MatrixType                     &inMatrix,
                                   const index_type                      inColumnIdx)
                : matrix(inMatrix), columnIdx(inColumnIdx)
            {
            }

            /** \brief returns the valus at position
             * (inRowIdx, columnIdx) in matrix.
             * \param inRowIdx the index of the value to access
             */
            const real &operator[](const index_type inRowIdx) const
            {
                return matrix(inRowIdx, columnIdx);
            }

            size_type size() const
            {
                return matrix.getLength1();
            }
    };

    public:
        /**
         * @brief pls_denham computes the pls regression
         *
         * The method will throw if there is invalid arguments
         * (matrices of invalid size, or output variables equal to nullptr).
         *
         * \param [in]  X  contains the centered (n,p)-matrix X.
         * \param [in]  y  contains the centered (n)-vector y.
         * \param [in]  n  the number of rows of the matrix X to use in analysis.
         * \param [in]  p  the number of columns of the matrix X to use in analysis.
         * \param [in]  k  the number of PLS factors to include in the regression of y on the matrix X.
         * \param [out] W  a (k,p)-matrix containing the coefficient vectors stored by column of
         *                 the a PLS regression factors obtained from the matrix X.
         * \param [out] q  a (k)-vector containing the least squares regression coefficients of
         *                 the ordinary least squares regression of y on the PLS factor matrix t.
         */
        void pls_denham(const MatrixType &X,
                        const VectorType &y,
                        const size_type n, const size_type p, const size_type k,
                        MatrixType* W,
                        VectorType* q) const
        {
            if (W == nullptr)
            {
                throw std::invalid_argument("Error in pls_denham W cannot be null.");
            }
            if (q == nullptr)
            {
                throw std::invalid_argument("Error in pls_denham q cannot be null.");
            }
            if (p < k)
            {
                throw std::invalid_argument("Error in pls_denham invalid k.");
            }
            if (q->size() != k)
            {
                throw std::invalid_argument("Error in pls_denham invalid q.");
            }
            if (W->getLength1() != p || W->getLength2() != k)
            {
                throw std::invalid_argument("Error in pls_denham invalid W.");
            }
            if (y.size() != n)
            {
                throw std::invalid_argument("Error in pls_denham invalid y.");
            }
            if (X.getLength1() != n || X.getLength2() != p)
            {
                throw std::invalid_argument("Error in pls_denham invalid y.");
            }

            // With : T_k = (t_1, ... , t_k)
            MatrixType T_k(index_type(0), index_type(n-1), index_type(0), index_type(k-1));
            // With : W_k = (W_1, ... , W_k)
            MatrixType W_k(index_type(0), index_type(p-1), index_type(0), index_type(k-1));

            // The comments in this method are a copy/past from
            // the algorithm in the paper.

            // w_i = X^T × y
            VectorType w_i(p);
            gmx::matrixProductMatrixTVector(X, y, w_i);

            // t_1 = X × w_1
            VectorType t_i(n);
            gmx::matrixProductMatrixVector(X, w_i, t_i);

            // Keep t_i and w_i in T_k and W_k
            setColumn(&T_k, 0, t_i, n);
            setColumn(&W_k, 0, w_i, p);

            // r = y - t_1 × (t_1^T × t_1)^-1 × t_1^T × y
            VectorType r(n);
            {
                const real t_1T_t_1            = std::inner_product(t_i.begin(), t_i.end(), t_i.begin(), 0);

                if (t_1T_t_1 == 0)
                {
                    throw std::runtime_error("Error in pls_denham t_1T_t_1 is equal to 0.");
                }

                const real t_1T_t_1_INV        = 1 / t_1T_t_1;
                const real t_1T_y              = std::inner_product(t_i.begin(), t_i.end(), y.begin(), 0);
                const real t_1T_t_1_INV_t_1T_y = t_1T_t_1_INV * t_1T_y;

                for (index_type idxVal = 0; idxVal < index_type(n); ++idxVal)
                {
                    r[idxVal] = y[idxVal] - t_i[idxVal] * t_1T_t_1_INV_t_1T_y;
                }
            }

            // for i = 2 ... k
            for (index_type i = 2; i < index_type(k); i++)
            {
                // w_i = X^T × r
                gmx::matrixProductMatrixTVector(X, r, w_i);

                // t_i = X × w_i
                gmx::matrixProductMatrixVector(X, w_i, t_i);

                // Keep t_i and w_i in T_k and W_k
                setColumn(&T_k, i-1, t_i, n);
                setColumn(&W_k, i-1, w_i, p);

                // With : T_k = (t_1, ... , t_i)
                MatrixType T_i(index_type(0), index_type(n-1), index_type(0), index_type(i-1));
                for (index_type i_copy = 0; i_copy < i; ++i_copy)
                {
                    setColumn(&T_i, i_copy, IrregArray2DColumnView(T_k, i_copy), n);
                }


                // r = y - T_i × (T_i^T × T_i)^-1 × T_i^T × y
                {
                    MatrixType T_iT_T_i(index_type(0), index_type(i-1), index_type(0), index_type(i-1));
                    gmx::matrixProductMatrixTMatrix(T_i, T_i, T_iT_T_i);

                    MatrixType T_iT_T_i_INV(index_type(0), index_type(i-1), index_type(0), index_type(i-1));
                    const bool inversionWorked = gmx::matrixInverseLUP(T_iT_T_i, T_iT_T_i_INV);

                    if (inversionWorked == false)
                    {
                        throw std::runtime_error("Error in pls_denham t_1T_t_1 encountered singular matrix.");
                    }

                    VectorType T_iT_y(i);
                    gmx::matrixProductMatrixTVector(T_i, y, T_iT_y);

                    VectorType T_iT_T_i_INV_T_iT_y(i);
                    gmx::matrixProductMatrixVector(T_iT_T_i_INV, T_iT_y, T_iT_T_i_INV_T_iT_y);

                    VectorType T_i_T_iT_T_i_INV_T_iT_y(n);
                    gmx::matrixProductMatrixVector(T_i, T_iT_T_i_INV_T_iT_y, T_i_T_iT_T_i_INV_T_iT_y);

                    for (index_type idxVal = 0; idxVal < index_type(n); ++idxVal)
                    {
                        r[idxVal] = y[idxVal] - T_i_T_iT_T_i_INV_T_iT_y[idxVal];
                    }
                }
            }
            if (k > 1)
            {
                // For i == k
                // w_i = X^T × r
                gmx::matrixProductMatrixTVector(X, r, w_i);

                // t_i = X × w_i
                gmx::matrixProductMatrixVector(X, w_i, t_i);

                // Keep t_i and w_i in T_k and W_k
                setColumn(&T_k, k-1, t_i, n);
                setColumn(&W_k, k-1, w_i, p);
            }

            // q = (T_k^T × T_k)^-1 × T_k^T × y
            {
                const index_type                               i = k;
                MatrixType       T_kT_T_k(index_type(0), index_type(i-1), index_type(0), index_type(i-1));
                gmx::matrixProductMatrixTMatrix(T_k, T_k, T_kT_T_k);

                MatrixType T_kT_T_k_INV(index_type(0), index_type(i-1), index_type(0), index_type(i-1));
                const bool inversionWorked = gmx::matrixInverseLUP(T_kT_T_k, T_kT_T_k_INV);

                if (inversionWorked == false)
                {
                    throw std::runtime_error("Error in pls_denham t_1T_t_1 encountered singular matrix.");
                }

                VectorType T_kT_y(i);
                gmx::matrixProductMatrixTVector(T_k, y, T_kT_y);

                VectorType T_kT_T_k_INV_T_kT_y(i);
                gmx::matrixProductMatrixVector(T_kT_T_k_INV, T_kT_y, T_kT_T_k_INV_T_kT_y);

                (*q) = std::move(T_kT_T_k_INV_T_kT_y);
            }

            // currently not needed: b = W_k × q
            // VectorType b(p);
            // gmx::matrixProductMatrixVector(W_k, q, b);

            // Keep the results
            for (index_type i = 0; i < index_type(p); ++i)
            {
                std::copy_n(W_k(i), k, (*W)(i));
            }
        }
};

#endif
