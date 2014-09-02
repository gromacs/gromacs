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

#ifndef GMX_LINEARALGEBRA_PLS_H
#define GMX_LINEARALGEBRA_PLS_H

#include <vector>

#include "gmx_blas.h"
#include "gmx_lapack.h"
#include "gromacs/legacyheaders/types/simple.h"

/*! \brief
 * This macro is used to call the LAPACK and BLAS functions required for
 * the preset GROMACS precision (i.e. real or single precision)
 */
#ifdef GMX_DOUBLE
#define GMX_LAPACK(A, B)  F77_FUNC(d ## A, D ## B)
#else
#define GMX_LAPACK(A, B)  F77_FUNC(s ## A, S ## B)
#endif

/*! \brief
 *  macro to supply a Vector/Matrix pointer to a FORTRAN function
 *  */
#define PTOF(m) (&m->front())
/*! \brief
 *  macro to supply a Vector/Matrix to a FORTRAN function
 */
#define TOF(m) (&m.front())

/*! \brief
 * vector of real values
 */
typedef std::vector<real> RVector;

/*! \brief Vector-based implementation of FORTRAN Column-Major Matrices
 *
 * as implemented by Evgenii B. Rudnyi, (http://MatrixProgramming.com)
 *
 */
class FMatrix : public std::vector<real>
{
/*! \brief
 * Implementation of a n*m matrix based on the vector class whose
 * internal format is consistent with FORTRAN odering
 */
  int m;
  int n;

public:
  FMatrix() : m(0), n(0)
  {
  /*! \brief
   * initialization of a 0*0 matrix
   */
  }
  FMatrix(int m_, int n_) : std::vector<real>(m_*n_), m(m_), n(n_)
  {
  /*! \brief
   * initialization of an n*m matrix (all fields are set to zero)
   *
   * \param m_ number of rows of the matrix
   * \param n_ number of columns of the matrix
   */
  }
  FMatrix(int m_, int n_, real val) : std::vector<real>(m_*n_, val), m(m_), n(n_)
  {
  /*! \brief
   * initialization of an n*m matrix with a standard vaule
   *
   * \param m_ number of rows of the matrix
   * \param n_ number of columns of the matrix
   * \param val initial value of the matrix fields
   */
  }
  void resize(int m_, int n_)
  {
  /*! \brief
   * resizing the matrix
   *
   * \param m_ new number of rows
   * \param n_ new number of columns
   */
	  m = m_; n = n_; std::vector<real>::resize(m_*n_);
   }
  void reserve(int m_, int n_)
  {
  /*! \brief
   * Requests that the matrix capacity is at least enough to contain
   * m_ rows and n_ columns
   *
   * \param m_ requested number of rows
   * \param n_ requested number of columns
   */
	  std::vector<real>::reserve(m_*n_);
   }
  void clear()
    {
  /*! \brief
   * Removes all elements from the matrix and sets the size to 0x0
   */
	  m = n = 0; std::vector<real>::clear();
    }

  int NRows() const
  {
  /*! \brief
   * retrieve the number of rows of the matrix
   *
   * \return number of rows
   */
	  return m;
  }
  int NCols() const
  {
  /*! \brief
   * retrieve the number of columns of the matrix
   *
   * \return number of columns
   */
   return n;
  }
  real& operator()(int i, int j)
  {
  /*! \brief
   * operator to access an element of the matrix
   *
   * \param i row of the element to retrieve
   * \param j column of the element to retrieve
   *
   * \return element (i,j)
   */

	  return operator[](i + j*m);
  }
  const real& operator()(int i, int j) const
  {
  /*! \brief
   * operator to access an element of the matrix
   *
   * \param i row of the element to retrieve
   * \param j column of the element to retrieve
   *
   * \return element (i,j)
   */
    return operator[](i + j*m);
  }
  void swap(FMatrix &y)
  {
  /*! \brief
   * Exchange the contents of the matrix with another
   *
   * \param y Matrix to swap values with
   */
	std::vector<real>::swap(y);
	std::swap(n, y.n);
	std::swap(m, y.m);
  }
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
pls_denham(FMatrix *x, RVector *y, int n, int k, int ln, int a, FMatrix *w, RVector *q);
#endif
