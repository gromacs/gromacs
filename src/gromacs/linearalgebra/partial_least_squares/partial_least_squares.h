/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
 * \author Berenger Bramas <berenger.bramas@mpcdf.mpg.de>
 * \author R. Thomas Ullmann <thomas.ullmann@mpibpc.mpg.de>
 * \author Luka Stanisic <luka.stanisic@mpcdf.mpg.de>
 *
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
 * Implemented as described in Sec 2.3 "Helland's algorithm" of:
 *
 * Denham, M. C., 1995, "Implementing Partial Least Squares.",
 * Statistics and Computing, 5, 191-202, doi:10.1007/BF00142661
 *
 * See the section "Functional Mode Analysis" of the reference manual for
 * a brief introduction to partial least-squares regression.
 *
 * \tparam   MatrixType   data type to store the input and output matrices
 * \tparam   VectorType   data type to store the input and output vectors
 *
 * The template types must provide different methods:
 * + MatrixType must provide
 *       - a size_type attribute type
 *       - a index_type attribute type
 *       - a length1() function to determine the array length in dimension 1
 *       - a length2() method to determine the array length in dimension 2
 *       - an operator (i,j) to access the value at row i and column j.
 * + VectorType must provide
 *       - a size() method to determine the vector length
 *       - operator[i] to access vector element i
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
     * This static method copy the column cl view into the column of index columnIdx of mt.
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
             * \param inRowIdx the index of the value to access.
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
         * \param [out] W  a (k,p)-matrix containing the coefficient vectors
         *                 of the PLS regression factors obtained for matrix X
         *                 Row k contains the coefficients of regressor k for the p variables x.
         * \param [out] q  a (k)-vector containing the least squares regression coefficients of
         *                 the ordinary least squares regression of y on the PLS factor matrix T.
         */
        void pls_denham(const MatrixType &X,
                        const VectorType &y,
                        MatrixType      * W,
                        VectorType      * q) const
        {
            // the number of rows of the matrix X, equal to the number of samples analyzed
            // (i.e., the number of structure-observable pairs to be explained by the PLS model)
            const size_type n = X.length1();
            // the number of columns of the matrix X equal to the number of input variables
            // (e.g., the number of coordinates of a protein considered in the analysis)
            const size_type p = X.length2();
            // the number of PLS factors or regressors to include in the regression of y on the matrix X
            const size_type k = q->size();
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
                throw std::invalid_argument("Error in pls_denham: need at least as many analyzed variables p as regressors k, i.e., p >= k.");
            }
            if (y.size() != n)
            {
                throw std::invalid_argument("Error in pls_denham: size of the vector y does not match the number of samples (rows of matrix x).");
            }
            if (W->length1() != p)
            {
                throw std::invalid_argument("Error in pls_denham: number of rows of matrix W must match the number of analyzed variables (columns of matrix X).");
            }
            if (W->length2() != k)
            {
                throw std::invalid_argument("Error in pls_denham: number of columns of matrix W must match the number of regressors (size of vector q).");
            }

            // With: T_k = (t_1, ... , t_k)
            MatrixType T_k(n, k);
            // With: W_k = (W_1, ... , W_k)
            MatrixType W_k(p, k);

            // The comments in this method are a copy/paste from
            // the algorithm in the paper.

            // w_i = X^T * y
            VectorType w_i(p);
            gmx::matrixProductMatrixTVector(X, y, w_i);

            // t_1 = X * w_1
            VectorType t_i(n);
            gmx::matrixProductMatrixVector(X, w_i, t_i);

            // Keep t_i and w_i in T_k and W_k
            setColumn(&T_k, 0, t_i, n);
            setColumn(&W_k, 0, w_i, p);

            // r = y - t_1 * (t_1^T * t_1)^-1 * t_1^T * y
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
                // w_i = X^T * r
                gmx::matrixProductMatrixTVector(X, r, w_i);

                // t_i = X * w_i
                gmx::matrixProductMatrixVector(X, w_i, t_i);

                // Keep t_i and w_i in T_k and W_k
                setColumn(&T_k, i-1, t_i, n);
                setColumn(&W_k, i-1, w_i, p);

                // With: T_k = (t_1, ... , t_i)
                MatrixType T_i(n, i);
                for (index_type i_copy = 0; i_copy < i; ++i_copy)
                {
                    setColumn(&T_i, i_copy, IrregArray2DColumnView(T_k, i_copy), n);
                }

                // r = y - T_i * (T_i^T * T_i)^-1 * T_i^T * y
                {
                    MatrixType T_iT_T_i(i, i);
                    gmx::matrixProductMatrixTMatrix(T_i, T_i, T_iT_T_i);

                    MatrixType T_iT_T_i_INV(i, i);
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
                // w_i = X^T * r
                gmx::matrixProductMatrixTVector(X, r, w_i);

                // t_i = X * w_i
                gmx::matrixProductMatrixVector(X, w_i, t_i);

                // Keep t_i and w_i in T_k and W_k
                setColumn(&T_k, k-1, t_i, n);
                setColumn(&W_k, k-1, w_i, p);
            }

            // q = (T_k^T * T_k)^-1 * T_k^T * y
            {
                const index_type i = k;
                MatrixType       T_kT_T_k(i, i);
                gmx::matrixProductMatrixTMatrix(T_k, T_k, T_kT_T_k);

                MatrixType T_kT_T_k_INV(i, i);
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

            // currently not needed: b = W_k * q
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
