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
/* \internal \file
 * \brief
 * Implements partial least squares regression
 *
 * \author Jan Henning Peters <JanHPeters@gmx.net>
 */

#ifndef GMX_LINEARALGEBRA_PLS_H
#define GMX_LINEARALGEBRA_PLS_H

#include <vector>

#include "gromacs/utility/real.h"

/*! \brief
 * Implementation of a n*m matrix based on the vector class whose
 * internal format is consistent with FORTRAN ordering
 */
class FMatrix
{
    std::vector<real> data;
    int               m;
    int               n;

    public:
        /*! \brief
         * initialization of a 0*0 matrix
         */
        FMatrix() : data(), m(0), n(0)
        {
        }
        /*! \brief
         * initialization of an n*m matrix (all fields are set to zero)
         *
         * \param m_ number of rows of the matrix
         * \param n_ number of columns of the matrix
         */
        FMatrix(int m_, int n_) : data(m_*n_), m(m_), n(n_)
        {
        }
        /*! \brief
         * initialization of an n*m matrix with a standard vaule
         *
         * \param m_ number of rows of the matrix
         * \param n_ number of columns of the matrix
         * \param val initial value of the matrix fields
         */
        FMatrix(int m_, int n_, real val) : data(m_*n_, val), m(m_), n(n_)
        {
        }
        /*! \brief
         * resizing the matrix
         *
         * \param m_ new number of rows
         * \param n_ new number of columns
         */
        void resize(int m_, int n_)
        {
            m = m_;
            n = n_;
            data.resize(m_*n_);
        }
        /*! \brief
         * Requests that the matrix capacity is at least enough to contain
         * m_ rows and n_ columns
         *
         * \param m_ requested number of rows
         * \param n_ requested number of columns
         */
        void reserve(int m_, int n_)
        {
            data.reserve(m_*n_);
        }
        /*! \brief
         * Removes all elements from the matrix and sets the size to 0x0
         */
        void clear()
        {
            m = n = 0;
            data.clear();
        }

        /*! \brief
         * retrieve the number of rows of the matrix
         *
         * \return number of rows
         */
        int NRows() const
        {
            return m;
        }
        /*! \brief
         * retrieve the number of columns of the matrix
         *
         * \return number of columns
         */
        int NCols() const
        {
            return n;
        }
        /*! \brief
         * operator to access an element of the matrix
         *
         * \param i row of the element to retrieve
         * \param j column of the element to retrieve
         *
         * \return element (i,j)
         */
        real &operator()(int i, int j)
        {

            return data[i + j*m];
        }
        /*! \brief
         * operator to access an element of the matrix
         *
         * \param i row of the element to retrieve
         * \param j column of the element to retrieve
         *
         * \return element (i,j)
         */
        const real &operator()(int i, int j) const
        {
            return data[i + j*m];
        }
        /*! \brief
         * return a pointer to the "FORTRAN-style" matrix representation
         */
        real * toF()
        {
            return &(data.front());
        }
        /*! \brief
         * Exchange the contents of the matrix with another
         *
         * \param y Matrix to swap values with
         */
        void swap(FMatrix &y)
        {
            data.swap(y.data);
            std::swap(n, y.n);
            std::swap(m, y.m);
        }
};


/*! \brief
 * vector of real values to supply to a FORTRAN function
 */
class FVector : private FMatrix
{
    public:
        /*! \brief
         * initialization of an empty vector
         */
        FVector() : FMatrix()
        {
        }
        /*! \brief
         * initialization of an n-vector (all fields are set to zero)
         *
         * \param n_ number of elements of the vector
         */
        FVector(int n_) : FMatrix(1, n_)
        {
        }
        /*! \brief
         * initialization of an n-vector with a standard vaule
         *
         * \param n_ number of elements of the vector
         * \param val initial value of the vector fields
         */
        FVector(int n_, real val) : FMatrix(1, n_, val)
        {
        }
        /*! \brief
         * resizing the vector
         *
         * \param n_ new number of elements
         */
        void resize(int n_)
        {
            FMatrix::resize(1, n_);
        }
        /*! \brief
         * Requests that the vector capacity is at least enough to contain
         * n_ elements
         *
         * \param n_ requested number of elements
         */
        void reserve(int n_)
        {
            FMatrix::reserve(1, n_);
        }

        /*! \brief
         * retrieve the number of elements of the vector
         *
         * \return number of elements
         */
        int Length() const
        {
            return FMatrix::NCols();
        }
        /*! \brief
         * operator to access an element of the vector
         *
         * \param i element to retrieve
         *
         * \return element i
         */
        real &operator[](int i)
        {
            return FMatrix::operator()(0, i);
        }
        /*! \brief
         * operator to access an element of the vector
         *
         * \param i row of the element to retrieve
         *
         * \return element i
         */
        const real &operator[](int i) const
        {
            return FMatrix::operator()(0, i);
        }

        using FMatrix::clear;
        using FMatrix::toF;
        using FMatrix::swap;
};

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
 * \param [in]  x  contains the centered (ln,k)-matrix X.
 * \param [in]  y  contains the centered (n)-vector y.
 * \param [in]  n  the number of rows of the matrix X to use in analysis. Any row beyond that
 *                 will be ignored
 * \param [in]  a  < MIN(n-1,k) - the number of PLS factors to include in the regression
 *                 of y on the matrix X.
 * \param [out] w  a (k,a)-matrix containing the coefficient vectors stored by column of
 *                 the a PLS regression factors obtained from the matrix X on .
 * \param [out] q  a (a)-vector containing the least squares regression coefficients of
 *                 the ordinary least squares regression of y on the PLS factor matrix t.
 *
 */
int
pls_denham(FMatrix *x, FVector *y, int n, int a, FMatrix *w, FVector *q);
#endif
