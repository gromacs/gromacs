/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
/*! \libinternal\file

   \brief   matrix and matrix vector products utilizing IrregArray2D and std::vector

   The functions are intentionally declared as free functions to allow
   storage of non-arithmetic data types in IrregArray2D, too.

 + vector matrix product \f$ \mathbf{w} = \mathbf{v}^{\mathrm{T}}\,\mathbf{A} \f$
       - input matrix \f$ \mathbf{A} \f$ stored as is                                       : \ref matrixProductVectorMatrix()
       - input matrix \f$ \mathbf{A} \f$ stored as transpose \f$ \mathbf{A}^{\mathrm{T}} \f$: \ref matrixProductVectorMatrixT()
 + matrix vector product \f$ \mathbf{w} = \mathbf{A}\,\mathbf{v} \f$
       - input matrix \f$ \mathbf{A} \f$ stored as is                                       : \ref matrixProductMatrixVector()
       - input matrix \f$ \mathbf{A} \f$ stored as transpose \f$ \mathbf{A}^{\mathrm{T}} \f$: \ref matrixProductMatrixTVector()
 + matrix product \f$ \mathbf{C} = \mathbf{A} \,\mathbf{B} \f$
       - input matrices stored as \f$ \mathbf{A},             \mathbf{B} \f$             : \ref matrixProductMatrixMatrix() \f$ \equiv \f$ \ref matrixProduct()
       - input matrices stored as \f$ \mathbf{A},             \mathbf{B}^{\mathrm{T}} \f$: \ref matrixProductMatrixMatrixT()
       - input matrices stored as \f$ \mathbf{A}^{\mathrm{T}},\mathbf{B} \f$             : \ref matrixProductMatrixTMatrix()
       - input matrices stored as \f$ \mathbf{A}^{\mathrm{T}},\mathbf{B}^{\mathrm{T}} \f$: \ref matrixProductMatrixTMatrixT()

   \note organization of the functions:
         For each mathematical function, there is a master function without any suffix
         that is intended to be called by user code (listed above) and suffixed
         implementation functions. The master function calls a scalar reference
         implementation (suffix Ref). A subsequent patch introduces SIMD implementations
         that will be called by the master function if there is SIMD support.

   \author R. Thomas Ullmann <tullman@gwdg.de>
   \inlibraryapi
   \ingroup module_math
 */
#ifndef GMX_MATH_IRREG_ARRAY_MATH_IRREG_ARRAY_2D_MATRIX_PRODUCTS_H
#define GMX_MATH_IRREG_ARRAY_MATH_IRREG_ARRAY_2D_MATRIX_PRODUCTS_H
#include <numeric>

#include "gromacs/utility/data_structures/irreg_array_2d.h"

namespace gmx
{

// --------------------------------------------------------------------------------------------------------------------------------------------------------

/*! \brief vector-matrix product \f$ \mathbf{w} = \mathbf{v}^{\mathrm{T}}\,\mathbf{A} \f$
           reference implementation

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    allocator type for T used by the function arguments

    \param[in]   a        matrix \f$ \mathbf{A} \f$
    \param[in]   v        vector \f$ \mathbf{v} \f$
    \param[out]  w        the result vector \f$ \mathbf{w} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline void matrixProductMatrixVectorRef(const IrregArray2D<T, Alloc, false> &a,
                                         const std::vector<T, Alloc> &v,
                                         std::vector<T, Alloc> &w)
{
    typedef typename IrregArray2D<T, Alloc, false>::value_type value_type;
    typedef typename IrregArray2D<T, Alloc, false>::index_type index_type;

    // number of rows    of matrix A
    const index_type m = a.length1();
    // number of columns of matrix A
    const index_type n = a.length2(0);

    std::fill(w.begin(), w.end(), static_cast<value_type>(0));

    for (index_type i = 0; i < m; ++i)
    {
        // w[i][1] = a[i][:] v[:][1]
        // scalar product of row vector i of matrix a and column vector v (col 1 of the n x 1 matrix)
        for (index_type j = 0; j < n; ++j)
        {
            w[i] += a(i, j) * v[j];
        }
    }
}

/*! \brief vector-matrix product \f$ \mathbf{w} = \mathbf{v}^{\mathrm{T}}\,\mathbf{A} \f$

    Formally, the function computes the matrix product of the \f$ 1 \times m \f$ vector
    \f$ \mathbf{v} \f$ and the \f$ m \times n \f$ matrix \f$ \mathbf{A} \f$ resulting in
    the \f$ n \times  1 \f$ vector \f$ \mathbf{w} \f$.
    Matrix \f$ \mathbf{A} \f$ is expected in transpose form.

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   a        matrix \f$ \mathbf{A} \f$
    \param[in]   v        vector \f$ \mathbf{v} \f$
    \param[out]  w        the result vector \f$ \mathbf{w} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline void matrixProductMatrixVector(const IrregArray2D<T, Alloc, false> &a,
                                      const std::vector<T, Alloc> &v,
                                      std::vector<T, Alloc> &w)
{
    return matrixProductMatrixVectorRef(a, v, w);
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------

/*! \brief vector-matrix product \f$ \mathbf{w} = \mathbf{v}^{\mathrm{T}}\,\mathbf{A} \f$
           reference implementation

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    allocator type for T used by the function arguments

    \param[in]   v        vector \f$ \mathbf{v} \f$
    \param[in]   aT       matrix \f$ \mathbf{A}^{\mathrm{T}} \f$
    \param[out]  w        the result vector \f$ \mathbf{w} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline void matrixProductVectorMatrixTRef(const std::vector<T, Alloc> &v,
                                          const IrregArray2D<T, Alloc, false> &aT,
                                          std::vector<T, Alloc> &w)
{
    typedef typename IrregArray2D<T, Alloc, false>::value_type value_type;
    typedef typename IrregArray2D<T, Alloc, false>::index_type index_type;

    // number of rows    of matrix A, columns of A^T
    const index_type m = aT.length2(0);
    // number of columns of matrix A, rows    of A^T
    const index_type n = aT.length1();

    std::fill(w.begin(), w.end(), static_cast<value_type>(0));

    for (index_type i = 0; i < n; ++i)
    {
        // w[1][i] = v[1][:] b[i][:], since b = a^T
        // w[1][i] = v[1][:] a[:][i]
        // scalar product of vector v (row 1 of the 1 x m matrix) and column vector j of matrix a
        for (index_type j = 0; j < m; ++j)
        {
            w[i] += v[j] * aT(i, j);
        }
    }
}

/*! \brief vector-matrix product \f$ \mathbf{w} = \mathbf{v}^{\mathrm{T}}\,\mathbf{A} \f$

    Formally, the function computes the matrix product of the \f$ 1 \times m \f$ vector
    \f$ \mathbf{v} \f$ and the \f$ m \times n \f$ matrix \f$ \mathbf{A} \f$ resulting in
    the \f$ n \times  1 \f$ vector \f$ w \f$.
    Matrix \f$ \mathbf{A} \f$ is expected in transpose form.

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   v        vector \f$ \mathbf{v} \f$
    \param[in]   aT       matrix \f$ \mathbf{A}^{\mathrm{T}} \f$
    \param[out]  w        the result matrix \f$ \mathbf{w} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline void matrixProductVectorMatrixT(const std::vector<T, Alloc> &v,
                                       const IrregArray2D<T, Alloc, false> &aT,
                                       std::vector<T, Alloc> &w)
{
    matrixProductVectorMatrixTRef(v, aT, w);
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------

/*! \brief vector-matrix product \f$ \mathbf{w} = \mathbf{v}^{\mathrm{T}}\,\mathbf{A} \f$,
           reference implementation

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   v        vector \f$ \mathbf{v} \f$
    \param[in]   a        matrix \f$ \mathbf{A} \f$
    \param[out]  w        the result vector \f$ \mathbf{w} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline void matrixProductVectorMatrixRef(const std::vector<T, Alloc> &v,
                                         const IrregArray2D<T, Alloc, false> &a,
                                         std::vector<T, Alloc> &w)
{
    typedef typename IrregArray2D<T, Alloc, false>::value_type value_type;
    typedef typename IrregArray2D<T, Alloc, false>::index_type index_type;

    // number of rows    of matrix A
    const index_type m = a.length1();
    // number of columns of matrix A
    const index_type n = a.length2(0);

    std::fill(w.begin(), w.end(), static_cast<value_type>(0));

    for (index_type i = 0; i < m; ++i)
    {
        // w[1][j] = v[1][:] a[:][j]
        // scalar product of vector v (row 1 of the 1 x m matrix) and column vector j of matrix a
        for (index_type j = 0; j < n; ++j)
        {
            w[j] += v[i] * a(i, j);
        }
    }
}

/*! \brief vector-matrix product \f$ \mathbf{w} = \mathbf{v}^{\mathrm{T}}\,\mathbf{A} \f$

    Formally, the function computes the matrix product of the \f$ 1 \times m \f$ vector
    \f$ \mathbf{v} \f$ and the \f$ m \times n \f$ matrix \f$ \mathbf{A} \f$ resulting in
    the \f$ n \times  1 \f$ vector \f$ w \f$.

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   v        vector \f$ \mathbf{v} \f$
    \param[in]   a        matrix \f$ \mathbf{A} \f$
    \param[out]  w        the result vector \f$ \mathbf{w} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline void matrixProductVectorMatrix(const std::vector<T, Alloc> &v,
                                      const IrregArray2D<T, Alloc, false> &a,
                                      std::vector<T, Alloc> &w)
{
    return matrixProductVectorMatrixRef(v, a, w);
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------

/*! \brief vector-matrix product \f$ \mathbf{w} = \mathbf{A}\,\mathbf{v} \f$,
           reference implementation

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   aT       matrix \f$ \mathbf{A}^{\mathrm{T}} \f$
    \param[in]   v        vector \f$ \mathbf{v} \f$
    \param[out]  w        the result vector \f$ \mathbf{w} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline void matrixProductMatrixTVectorRef(const IrregArray2D<T, Alloc, false> &aT,
                                          const std::vector<T, Alloc> &v,
                                          std::vector<T, Alloc> &w)
{
    typedef typename IrregArray2D<T, Alloc, false>::value_type value_type;
    typedef typename IrregArray2D<T, Alloc, false>::index_type index_type;

    // number of rows    of matrix A, columns of A^T
    const index_type m = aT.length2(0);
    // number of columns of matrix A, rows    of A^T
    const index_type n = aT.length1();

    std::fill(w.begin(), w.end(), static_cast<value_type>(0));

    for (index_type i = 0; i < m; ++i)
    {
        // w[1][j] = v[1][:] a[:][j]
        // scalar product of column vector j of matrix aT and column vector v (col 1 of the n x 1 matrix)
        for (index_type j = 0; j < n; ++j)
        {
            w[i] += aT(j, i) * v[j];
        }
    }
}

/*! \brief vector-matrix product \f$ \mathbf{w} = \mathbf{A}\,\mathbf{v} \f$

    Formally, the function computes the matrix product of the \f$ 1 \times m \f$ vector
    \f$ \mathbf{v} \f$ and the \f$ m \times n \f$ matrix \f$ \mathbf{A} \f$ resulting in
    the \f$ n \times  1 \f$ vector \f$ w \f$.

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   aT       matrix \f$ \mathbf{A}^{\mathrm{T}} \f$
    \param[in]   v        vector \f$ \mathbf{v} \f$
    \param[out]  w        the result vector \f$ \mathbf{w} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline void matrixProductMatrixTVector(const IrregArray2D<T, Alloc, false> &aT,
                                       const std::vector<T, Alloc> &v,
                                       std::vector<T, Alloc> &w)
{
    return matrixProductMatrixTVectorRef(aT, v, w);
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------

/*! \brief matrix product \f$ \mathbf{C} = \mathbf{A}\,\mathbf{B} \f$,
           reference implementation

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   a        \f$ m \times l \f$ matrix \f$ \mathbf{A} \f$
    \param[in]   b        \f$ l \times n \f$ matrix \f$ \mathbf{B} \f$
    \param[out]  c        the result: \f$ m \times  n \f$ matrix \f$ \mathbf{C} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline void matrixProductRef(const IrregArray2D<T, Alloc, false> &a,
                             const IrregArray2D<T, Alloc, false> &b,
                             IrregArray2D<T, Alloc, false> &c)
{
    typedef typename IrregArray2D<T, Alloc, false>::value_type value_type;
    typedef typename IrregArray2D<T, Alloc, false>::index_type index_type;

    // number of rows    of matrix A
    const index_type m = a.length1();
    // number of columns of matrix A
    const index_type l = a.length2(0);
    // number of columns of matrix B
    const index_type n = b.length2(0);

    // initialize matrix c one row vector at a time
    for (index_type i = 0; i < m; ++i)
    {
        std::fill(c(i), c(i) + static_cast<ptrdiff_t>(n), static_cast<value_type>(0));
    }

    for (index_type i = 0; i < m; ++i)
    {
        // c[i][j] = a[i][:] b[:][j]
        // scalar product of row vector i of a and column vector j of b
        for (index_type k = 0; k < l; ++k)
        {
            for (index_type j = 0; j < n; ++j)
            {
                c(i, j) += a(i, k) * b(k, j);
            }
        }
    }
}

/*! \brief matrix product \f$ \mathbf{C} = \mathbf{A}\,\mathbf{B} \f$

    Formally, the function computes the matrix product of the \f$ m \times l \f$ matrix
    \f$ \mathbf{A} \f$ and the \f$ l \times n \f$ matrix \f$ \mathbf{B} \f$ resulting in
    the \f$ m \times  n \f$ matrix \f$ \mathbf{C} \f$.

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   a        \f$ m \times l \f$ matrix \f$ \mathbf{A} \f$
    \param[in]   b        \f$ l \times n \f$ matrix \f$ \mathbf{B} \f$
    \param[out]  c        the result: \f$ m \times  n \f$ matrix \f$ \mathbf{C} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline void matrixProduct(const IrregArray2D<T, Alloc, false> &a,
                          const IrregArray2D<T, Alloc, false> &b,
                          IrregArray2D<T, Alloc, false> &c)
{
    return matrixProductRef(a, b, c);
}


/*! \brief matrix product \f$ \mathbf{C} = \mathbf{A}\,\mathbf{B} \f$

    Formally, the function computes the matrix product of the \f$ m \times l \f$ matrix
    \f$ \mathbf{A} \f$ and the \f$ l \times n \f$ matrix \f$ \mathbf{B} \f$ resulting in
    the \f$ m \times  n \f$ matrix \f$ \mathbf{C} \f$.

    Returns \ref matrixProduct(), here for accessibility/consistent naming.

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   a        \f$ m \times l \f$ matrix \f$ \mathbf{A} \f$
    \param[in]   b        \f$ l \times n \f$ matrix \f$ \mathbf{B} \f$
    \param[out]  c        the result: \f$ m \times  n \f$ matrix \f$ \mathbf{C} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline void matrixProductMatrixMatrix(const IrregArray2D<T, Alloc, false> &a,
                                      const IrregArray2D<T, Alloc, false> &b,
                                      IrregArray2D<T, Alloc, false> &c)
{
    return matrixProduct(a, b, c);
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------

/*! \brief matrix product \f$ \mathbf{C} = \mathbf{A}\,\mathbf{B} \f$,
           reference implementation

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   a        \f$ m \times l \f$ matrix \f$ \mathbf{A} \f$
    \param[in]   bT       \f$ n \times l \f$ matrix \f$ \mathbf{B}^{\mathrm{T}} \f$
    \param[out]  c        the result: \f$ m \times  n \f$ matrix \f$ \mathbf{C} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline void matrixProductMatrixMatrixTRef(const IrregArray2D<T, Alloc, false> &a,
                                          const IrregArray2D<T, Alloc, false> &bT,
                                          IrregArray2D<T, Alloc, false> &c)
{
    typedef typename IrregArray2D<T, Alloc, false>::value_type value_type;
    typedef typename IrregArray2D<T, Alloc, false>::index_type index_type;

    // number of rows    of matrix A
    const index_type m = a.length1();
    // number of columns of matrix A
    const index_type l = a.length2(0);
    // number of columns of matrix A
    const index_type n = bT.length1();

    // initialize matrix c one row vector at a time
    for (index_type i = 0; i < m; ++i)
    {
        std::fill(c(i), c(i) + static_cast<ptrdiff_t>(n), static_cast<value_type>(0));
    }

    for (index_type i = 0; i < m; ++i)
    {
        // c[i][j] = a[i][:] bT[:][j]
        // scalar product of row vector i of a and column vector j of b
        for (index_type j = 0; j < n; ++j)
        {
            T tmpf = 0;
            for (index_type k = 0; k < l; ++k)
            {
                tmpf += a(i, k) * bT(j, k);
            }
            c(i, j) = tmpf;
        }
    }
}

/*! \brief matrix product \f$ \mathbf{C} = \mathbf{A}\,\mathbf{B} \f$

    Formally, the function computes the matrix product of the \f$ m \times l \f$ matrix
    \f$ \mathbf{A} \f$ and the \f$ l \times n \f$ matrix \f$ \mathbf{B} \f$ resulting in
    the \f$ m \times  n \f$ matrix \f$ \mathbf{C} \f$.

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   a        \f$ m \times l \f$ matrix \f$ \mathbf{A} \f$
    \param[in]   bT       \f$ n \times l \f$ matrix \f$ \mathbf{B}^{\mathrm{T}} \f$
    \param[out]  c        the result: \f$ m \times  n \f$ matrix \f$ \mathbf{C} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline void matrixProductMatrixMatrixT(const IrregArray2D<T, Alloc, false> &a,
                                       const IrregArray2D<T, Alloc, false> &bT,
                                       IrregArray2D<T, Alloc, false> &c)
{
    return matrixProductMatrixMatrixTRef(a, bT, c);
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------

/*! \brief matrix product \f$ \mathbf{C} = \mathbf{A}\,\mathbf{B} \f$,
           reference implementation

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   aT       \f$ l \times m \f$ matrix \f$ \mathbf{A}^{\mathrm{T}} \f$
    \param[in]   bT       \f$ n \times l \f$ matrix \f$ \mathbf{B}^{\mathrm{T}} \f$
    \param[out]  c        the result: \f$ m \times  n \f$ matrix \f$ \mathbf{C} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline void matrixProductMatrixTMatrixTRef(const IrregArray2D<T, Alloc, false> &aT,
                                           const IrregArray2D<T, Alloc, false> &bT,
                                           IrregArray2D<T, Alloc, false> &c)
{
    typedef typename IrregArray2D<T, Alloc, false>::value_type value_type;
    typedef typename IrregArray2D<T, Alloc, false>::index_type index_type;

    // number of rows    of matrix A^T
    const index_type l = aT.length1();
    // number of columns of matrix A^T
    const index_type m = aT.length2(0);
    // number of rows    of matrix B^T
    const index_type n = bT.length1();

    // initialize matrix c one row vector at a time
    for (index_type i = 0; i < m; ++i)
    {
        std::fill(c(i), c(i) + static_cast<ptrdiff_t>(n), static_cast<value_type>(0));
    }

    for (index_type j = 0; j < n; ++j)
    {
        for (index_type i = 0; i < m; ++i)
        {
            // c[i][j] = aT[:][i] bT[j][:]
            // scalar product of row vector i of a and column vector j of b
            T tmpf = 0;
            for (index_type k = 0; k < l; ++k)
            {
                tmpf += aT(k, i) * bT(j, k);
            }
            c(i, j) = tmpf;
        }
    }
}

/*! \brief matrix product \f$ \mathbf{C} = \mathbf{A}\,\mathbf{B} \f$

    Formally, the function computes the matrix product of the \f$ m \times l \f$ matrix
    \f$ \mathbf{A} \f$ and the \f$ l \times n \f$ matrix \f$ \mathbf{B} \f$ resulting in
    the \f$ m \times  n \f$ matrix \f$ \mathbf{C} \f$.

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   aT       \f$ l \times m \f$ matrix \f$ \mathbf{A}^{\mathrm{T}} \f$
    \param[in]   bT       \f$ n \times l \f$ matrix \f$ \mathbf{B}^{\mathrm{T}} \f$
    \param[out]  c        the result: \f$ m \times  n \f$ matrix \f$ \mathbf{C} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline void matrixProductMatrixTMatrixT(const IrregArray2D<T, Alloc, false> &aT,
                                        const IrregArray2D<T, Alloc, false> &bT,
                                        IrregArray2D<T, Alloc, false> &c)
{
    return matrixProductMatrixTMatrixTRef(aT, bT, c);
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------

/*! \brief matrix product \f$ \mathbf{C} = \mathbf{A}\,\mathbf{B} \f$,
           reference implementation

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   aT       \f$ l \times m \f$ matrix \f$ \mathbf{A}^{\mathrm{T}} \f$
    \param[in]   b        \f$ l \times n \f$ matrix \f$ \mathbf{B} \f$
    \param[out]  c        the result: \f$ m \times  n \f$ matrix \f$ \mathbf{C} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline void matrixProductMatrixTMatrixRef(const IrregArray2D<T, Alloc, false> &aT,
                                          const IrregArray2D<T, Alloc, false> &b,
                                          IrregArray2D<T, Alloc, false> &c)
{
    typedef typename IrregArray2D<T, Alloc, false>::value_type value_type;
    typedef typename IrregArray2D<T, Alloc, false>::index_type index_type;

    // number of rows    of matrix A
    const index_type m = aT.length2(0);
    // number of columns of matrix A
    const index_type l = aT.length1();
    // number of columns of matrix A
    const index_type n = b.length2(0);

    // initialize matrix c one row vector at a time
    for (index_type i = 0; i < m; ++i)
    {
        std::fill(c(i), c(i) + static_cast<ptrdiff_t>(n), static_cast<value_type>(0));
    }

    for (index_type k = 0; k < l; ++k)
    {
        // c[i][j] = a[i][:] b[:][j]
        // scalar product of row vector i of a and column vector j of b
        for (index_type i = 0; i < m; ++i)
        {
            for (index_type j = 0; j < n; ++j)
            {
                c(i, j) += aT(k, i) * b(k, j);
            }
        }
    }
}

/*! \brief matrix product \f$ \mathbf{C} = \mathbf{A}\,\mathbf{B} \f$

    Formally, the function computes the matrix product of the \f$ m \times l \f$ matrix
    \f$ \mathbf{A} \f$ and the \f$ l \times n \f$ matrix \f$ \mathbf{B} \f$ resulting in
    the \f$ m \times  n \f$ matrix \f$ \mathbf{C} \f$.

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   aT       \f$ l \times m \f$ matrix \f$ \mathbf{A}^{\mathrm{T}} \f$
    \param[in]   b        \f$ l \times n \f$ matrix \f$ \mathbf{B} \f$
    \param[out]  c        the result: \f$ m \times  n \f$ matrix \f$ \mathbf{C} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline void matrixProductMatrixTMatrix(const IrregArray2D<T, Alloc, false> &aT,
                                       const IrregArray2D<T, Alloc, false> &b,
                                       IrregArray2D<T, Alloc, false> &c)
{
    return matrixProductMatrixTMatrixRef(aT, b, c);
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------

} // end namespace gmx

#endif
