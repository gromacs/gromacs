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
         implementation (suffix Ref) if SIMD is disabled at compile time or a SIMD
         vectorized function (suffix Simd) if SIMD is enabled. There are three variants
         of the SIMD function only one of which is enabled at compile time:
         1) Typically, SIMD operations are only supported for float and double.
            If there is no hardware SIMD support for the employed floating point
            data type, the selected function XyzSimd simply returns the scalar version.
         2) If there is SIMD support and support for vector operations on (shorter)
            vectors of length 4 (SIMD4), a SIMD function is enabled that processes
            longer stretches of data in full SIMD width, the rest as far as possible
            with SIMD4 operations and the remainder with scalar math.
         3) If there is SIMD support, but no SIMD4 support, the enabled function uses
            the full SIMD width as far as possible and processes the remainder of the
            data with scalar math.

   \author R. Thomas Ullmann <tullman@gwdg.de>
   \inlibraryapi
   \ingroup module_math
 */
#ifndef GMX_MATH_IRREG_ARRAY_MATH_IRREG_ARRAY_2D_MATRIX_PRODUCTS_H
#define GMX_MATH_IRREG_ARRAY_MATH_IRREG_ARRAY_2D_MATRIX_PRODUCTS_H
#include <numeric>
#include <type_traits>

#include "gromacs/simd/simd_math.h"
#include "gromacs/simd/simd_setup.h"
#include "gromacs/utility/data_structures/irreg_array_2d.h"

namespace gmx
{

// --------------------------------------------------------------------------------------------------------------------------------------------------------

/*! \brief vector-matrix product \f$ \mathbf{w} = \mathbf{v}^{\mathrm{T}}\,\mathbf{A} \f$
           Currently, SIMD support is limited to float and/or double depending on the architecture.
           This template function is only active for unsupported data types and returns the reference
           implementation.

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   a        matrix \f$ \mathbf{A} \f$
    \param[in]   v        vector \f$ \mathbf{v} \f$
    \param[out]  w        the result vector \f$ \mathbf{w} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline
typename std::enable_if<!(gmx::SimdSetup<T>::doSimd_), void>::type
matrixProductMatrixVectorSimd(const IrregArray2D<T, Alloc, false> &a,
                              const std::vector<T, Alloc> &v,
                              std::vector<T, Alloc> &w)
{
    return matrixProductMatrixVectorRef(a, v, w);
}

#if GMX_SIMD // can't get rid of these guards because the functions load store etc are undefined if !GMX_SIMD

/*! \brief vector-matrix product \f$ \mathbf{w} = \mathbf{v}^{\mathrm{T}}\,\mathbf{A} \f$
           SIMD implementation, only enabled for float or double, if the respective data
           type is supported by the GROMACS simd module on the current architecture.
           This variant also requires reduced-SIMD-width SIMD4 support.

           The data is expected to be aligned. In principle, the alignment requirements
           could be matched against values provided by the allocator at compile time, but,
           since there is no standardized way for allocators to provide this information,
           the template function would not work for allocators that provide this information
           differently or not at all.

           \verbatim
            m row vectors    column vector      column vector
            of length n      of length n        of length m
            - - - -
            a a a a          v |                w |
            : : ...          v |                :
            : : ...          v |                :
            : : ...      x   v |             =  :

            : : ...          :                  :
            : : ...          :                  :
           \endverbatim

           use SIMD vector operations for the marked parts that match the SIMD width

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    allocator type for T used by the function arguments

    \param[in]   a        matrix \f$ \mathbf{A} \f$
    \param[in]   v        vector \f$ \mathbf{v} \f$
    \param[out]  w        the result vector \f$ \mathbf{w} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline
typename std::enable_if<gmx::SimdSetup<T>::doSimd_ &&gmx::SimdSetup<T>::doSimd4_, void>::type
matrixProductMatrixVectorSimd(const IrregArray2D<T, Alloc, false> &a,
                              const std::vector<T, Alloc> &v,
                              std::vector<T, Alloc> &w)
{
    typedef typename IrregArray2D<T, Alloc, false>::value_type value_type;
    typedef typename IrregArray2D<T, Alloc, false>::index_type index_type;

    // set up at compile time
    static constexpr index_type simdWidth  = SimdSetup<T>::simdWidth_;
    static constexpr index_type simd4Width = SimdSetup<T>::simd4Width_;
    typedef typename   SimdSetup<T>::simd_type     simd_float_type;
    typedef typename   SimdSetup<T>::simd4_type   simd4_float_type;

    // number of rows    of matrix A
    const index_type m = a.length1();
    std::fill(w.begin(), w.end(), static_cast<value_type>(0));
    // number of columns of matrix A
    const index_type n = a.length2(0);

    // stores stripes of the row vector v
    simd_float_type  vSimd(0);
    // stores stripes of column vectors of a
    simd_float_type  aSimd(0);
    // stores stripes of the result column vectors w
    simd_float_type  wSimd(0);

    const index_type simdBlockWidth = n - n % simdWidth;
    if (simdBlockWidth != 0)
    {
        for (index_type i = 0; i < m; ++i)
        {
            // scalar product of vector v (row 1 of the 1 x m matrix) and column vector j of matrix a
            for (index_type j = 0; j < simdBlockWidth; j += simdWidth)
            {
                vSimd = gmx::load<simd_float_type>(v.data() + j);
                aSimd = gmx::load<simd_float_type>(a(i) + j);
                wSimd = gmx::fma(vSimd, aSimd, wSimd);
            }
            w[i] += gmx::reduce(wSimd);
            wSimd = gmx::setZero();
        }
    }

    const index_type rest            = n - simdBlockWidth;
    const index_type simd4BlockWidth = rest - rest % simd4Width;
    // the Simd4 block should be completely optimized away if not active
    if (simd4BlockWidth != 0)
    {
        // stores stripes of the row vector v
        simd4_float_type vSimd4(0);
        // stores stripes of column vectors of a
        simd4_float_type aSimd4(0);
        // stores stripes of the result column vectors w
        simd4_float_type wSimd4(0);

        for (index_type i = 0; i < m; ++i)
        {
            // scalar product of vector v (row 1 of the 1 x m matrix) and column vector j of matrix a
            for (index_type j = simdBlockWidth; j < (simdBlockWidth + simd4BlockWidth); j += simd4Width)
            {
                vSimd4 = gmx::load4(v.data() + j);
                aSimd4 = gmx::load4(a(i) + j);
                wSimd4 = gmx::fma(vSimd4, aSimd4, wSimd4);
            }
            w[i]  += gmx::reduce(wSimd4);
            wSimd4 = gmx::setZero();
        }
    }

    // compute the remaining parts with scalar math
    for (index_type i = 0; i < m; ++i)
    {
        // c[1][j] = v[1][:] b[j][:], since b = a^T
        // c[1][j] = v[1][:] a[:][j]
        // scalar product of vector v (row 1 of the 1 x m matrix) and column vector j of matrix a
        for (index_type j = (simdBlockWidth + simd4BlockWidth); j < n; ++j)
        {
            w[i] += a(i, j) * v[j];
        }
    }
}

/*! \brief vector-matrix product \f$ \mathbf{w} = \mathbf{v}^{\mathrm{T}}\,\mathbf{A} \f$
           SIMD implementation, only enabled for float or double, if the respective data
           type is supported by the GROMACS simd module on the current architecture.

           \verbatim
            m row vectors    column vector      column vector
            of length n      of length n        of length m
            - - - -
            a a a a          v |                w |
            : : ...          v |                :
            : : ...          v |                :
            : : ...      x   v |             =  :

            : : ...     :    :                  :
            : : ...     :    :                  :
           \endverbatim

           use SIMD vector operations for the marked parts that match the SIMD width

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    allocator type for T used by the function arguments

    \param[in]   a        matrix \f$ \mathbf{A} \f$
    \param[in]   v        vector \f$ \mathbf{v} \f$
    \param[out]  w        the result vector \f$ \mathbf{w} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline
typename std::enable_if<gmx::SimdSetup<T>::doSimd_ && !(gmx::SimdSetup<T>::doSimd4_), void>::type
matrixProductMatrixVectorSimd(const IrregArray2D<T, Alloc, false> &a,
                              const std::vector<T, Alloc> &v,
                              std::vector<T, Alloc> &w)
{
    typedef typename IrregArray2D<T, Alloc, false>::value_type value_type;
    typedef typename IrregArray2D<T, Alloc, false>::index_type index_type;

    // set up at compile time
    static constexpr index_type simdWidth  = SimdSetup<T>::simdWidth_;
    typedef typename   SimdSetup<T>::simd_type     simd_float_type;

    // number of rows    of matrix A
    const index_type m = a.length1();
    std::fill(w.begin(), w.end(), static_cast<value_type>(0));
    // number of columns of matrix A
    const index_type n = a.length2(0);

    // stores stripes of the row vector v
    simd_float_type  vSimd(0);
    // stores stripes of column vectors of a
    simd_float_type  aSimd(0);
    // stores stripes of the result column vectors w
    simd_float_type  wSimd(0);

    const index_type simdBlockWidth = n - n % simdWidth;
    if (simdBlockWidth != 0)
    {
        for (index_type i = 0; i < m; ++i)
        {
            // scalar product of vector v (row 1 of the 1 x m matrix) and column vector j of matrix a
            for (index_type j = 0; j < simdBlockWidth; j += simdWidth)
            {
                vSimd = gmx::load<simd_float_type>(v.data() + j);
                aSimd = gmx::load<simd_float_type>(a(i) + j);
                wSimd = gmx::fma(vSimd, aSimd, wSimd);
            }
            w[i] += gmx::reduce(wSimd);
            wSimd = gmx::setZero();
        }
    }

    // compute the remaining parts with scalar math
    for (index_type i = 0; i < m; ++i)
    {
        // c[1][j] = v[1][:] b[j][:], since b = a^T
        // c[1][j] = v[1][:] a[:][j]
        // scalar product of vector v (row 1 of the 1 x m matrix) and column vector j of matrix a
        for (index_type j = simdBlockWidth; j < n; ++j)
        {
            w[i] += a(i, j) * v[j];
        }
    }
}

#endif // GMX_SIMD

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
#if GMX_SIMD
    return matrixProductMatrixVectorSimd(a, v, w);
#else
    return matrixProductMatrixVectorRef(a, v, w);
#endif
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------

/*! \brief vector-matrix product \f$ \mathbf{w} = \mathbf{v}^{\mathrm{T}}\,\mathbf{A} \f$
           Currently, SIMD support is limited to float and/or double depending on the architecture.
           This template function is only active for unsupported data types and returns the reference
           implementation.

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   v        vector \f$ \mathbf{v} \f$
    \param[in]   aT       matrix \f$ \mathbf{A}^{\mathrm{T}} \f$
    \param[out]  w        the result vector \f$ \mathbf{w} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline
typename std::enable_if<!(gmx::SimdSetup<T>::doSimd_), void>::type
matrixProductVectorMatrixTSimd(const std::vector<T, Alloc> &v,
                               const IrregArray2D<T, Alloc, false> &aT,
                               std::vector<T, Alloc> &w)
{
    return matrixProductVectorMatrixTRef(v, aT, w);
}

#if GMX_SIMD // can't get rid of these guards because the functions load store etc are undefined if !GMX_SIMD

/*! \brief vector-matrix product \f$ \mathbf{w} = \mathbf{v}^{\mathrm{T}}\,\mathbf{A} \f$
           SIMD implementation, only enabled for float or double, if the respective data
           type is supported by the GROMACS simd module on the current architecture.
           This variant also requires reduced-SIMD-width SIMD4 support.

           \verbatim
            row vector    n column vectors    row vector
            of length n   of of length m      of length m
            ----                              -
            vvvv ... vv       a | ...         wwww ... ww
                              a | ...
                              a | ...
                         x    a | ...        =
                              :   ...
                              a   ...
                              a   ...
           \endverbatim

           use SIMD vector operations for the marked parts that match the SIMD width

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    allocator type for T used by the function arguments

    \param[in]   v        vector \f$ \mathbf{v} \f$
    \param[in]   aT       matrix \f$ \mathbf{A}^{\mathrm{T}} \f$
    \param[out]  w        the result vector \f$ \mathbf{w} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline
typename std::enable_if<gmx::SimdSetup<T>::doSimd_ &&gmx::SimdSetup<T>::doSimd4_, void>::type
matrixProductVectorMatrixTSimd(const std::vector<T, Alloc> &v,
                               const IrregArray2D<T, Alloc, false> &aT,
                               std::vector<T, Alloc> &w)
{
    typedef typename IrregArray2D<T, Alloc, false>::value_type value_type;
    typedef typename IrregArray2D<T, Alloc, false>::index_type index_type;

    // set up at compile time
    static constexpr index_type simdWidth  = SimdSetup<T>::simdWidth_;
    static constexpr index_type simd4Width = SimdSetup<T>::simd4Width_;
    typedef typename   SimdSetup<T>::simd_type     simd_float_type;
    typedef typename   SimdSetup<T>::simd4_type   simd4_float_type;

    // number of rows    of matrix A, columns of A^T
    const index_type m = aT.length2(0);
    // number of columns of matrix A, rows    of A^T
    const index_type n = aT.length1();
    std::fill(w.begin(), w.end(), static_cast<value_type>(0));

    // stores stripes of the row vector v
    simd_float_type  vSimd(0);
    // stores stripes of column vectors of a
    simd_float_type  aSimd(0);
    // stores stripes of the result column vectors w
    simd_float_type  wSimd(0);

    const index_type simdBlockWidth = m - m % simdWidth;
    if (simdBlockWidth != 0)
    {
        for (index_type i = 0; i < n; ++i)
        {
            // scalar product of vector v (row 1 of the 1 x m matrix) and column vector j of matrix a
            for (index_type j = 0; j < simdBlockWidth; j += simdWidth)
            {
                vSimd = gmx::load<simd_float_type>(v.data() + j);
                aSimd = gmx::load<simd_float_type>(aT(i) + j);
                wSimd = gmx::fma(vSimd, aSimd, wSimd);
            }
            w[i] += gmx::reduce(wSimd);
            wSimd = gmx::setZero();
        }
    }

    const index_type rest            = m - simdBlockWidth;
    const index_type simd4BlockWidth = rest - rest % simd4Width;
    // the Simd4 block should be completely optimized away if not active
    if (simd4BlockWidth != 0)
    {
        // stores stripes of the row vector v
        simd4_float_type vSimd4(0);
        // stores stripes of column vectors of a
        simd4_float_type aSimd4(0);
        // stores stripes of the result column vectors w
        simd4_float_type wSimd4(0);

        for (index_type i = 0; i < n; ++i)
        {
            // scalar product of vector v (row 1 of the 1 x m matrix) and column vector j of matrix a
            for (index_type j = simdBlockWidth; j < (simdBlockWidth + simd4BlockWidth); j += simd4Width)
            {
                vSimd4 = gmx::load4(v.data() + j);
                aSimd4 = gmx::load4(aT(i) + j);
                wSimd4 = gmx::fma(vSimd4, aSimd4, wSimd4);
            }
            w[i]  += gmx::reduce(wSimd4);
            wSimd4 = gmx::setZero();
        }
    }

    // compute the remaining parts with scalar math
    for (index_type i = 0; i < n; ++i)
    {
        // c[1][j] = v[1][:] b[j][:], since b = a^T
        // c[1][j] = v[1][:] a[:][j]
        // scalar product of vector v (row 1 of the 1 x m matrix) and column vector j of matrix a
        for (index_type j = (simdBlockWidth + simd4BlockWidth); j < m; ++j)
        {
            w[i] += v[j] * aT(i, j);
        }
    }
}

/*! \brief vector-matrix product \f$ \mathbf{w} = \mathbf{v}^{\mathrm{T}}\,\mathbf{A} \f$
           SIMD implementation, only enabled for float or double, if the respective data
           type is supported by the GROMACS simd module on the current architecture.

           \verbatim
            row vector    n column vectors    row vector
            of length n   of of length m      of length m
            ----                              -
            vvvv ... vv       a | ...         wwww ... ww
                              a | ...
                              a | ...
                         x    a | ...        =
                              :   ...
                              a   ...
                              a   ...
           \endverbatim

           use SIMD vector operations for the marked parts that match the SIMD width

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    allocator type for T used by the function arguments

    \param[in]   v        vector \f$ \mathbf{v} \f$
    \param[in]   aT       matrix \f$ \mathbf{A}^{\mathrm{T}} \f$
    \param[out]  w        the result vector \f$ \mathbf{w} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline
typename std::enable_if<gmx::SimdSetup<T>::doSimd_ && !(gmx::SimdSetup<T>::doSimd4_), void>::type
matrixProductVectorMatrixTSimd(const std::vector<T, Alloc> &v,
                               const IrregArray2D<T, Alloc, false> &aT,
                               std::vector<T, Alloc> &w)
{
    typedef typename IrregArray2D<T, Alloc, false>::value_type value_type;
    typedef typename IrregArray2D<T, Alloc, false>::index_type index_type;

    // set up at compile time
    static constexpr index_type simdWidth  = SimdSetup<T>::simdWidth_;
    typedef typename   SimdSetup<T>::simd_type     simd_float_type;

    // number of rows    of matrix A, columns of A^T
    const index_type m = aT.length2(0);
    // number of columns of matrix A, rows    of A^T
    const index_type n = aT.length1();
    std::fill(w.begin(), w.end(), static_cast<value_type>(0));

    // stores stripes of the row vector v
    simd_float_type  vSimd(0);
    // stores stripes of column vectors of a
    simd_float_type  aSimd(0);
    // stores stripes of the result column vectors w
    simd_float_type  wSimd(0);

    const index_type simdBlockWidth = m - m % simdWidth;
    if (simdBlockWidth != 0)
    {
        for (index_type i = 0; i < n; ++i)
        {
            // scalar product of vector v (row 1 of the 1 x m matrix) and column vector j of matrix a
            for (index_type j = 0; j < simdBlockWidth; j += simdWidth)
            {
                vSimd = gmx::load<simd_float_type>(v.data() + j);
                aSimd = gmx::load<simd_float_type>(aT(i) + j);
                wSimd = gmx::fma(vSimd, aSimd, wSimd);
            }
            w[i] += gmx::reduce(wSimd);
            wSimd = gmx::setZero();
        }
    }

    // compute the remaining parts with scalar math
    for (index_type i = 0; i < n; ++i)
    {
        // c[1][j] = v[1][:] b[j][:], since b = a^T
        // c[1][j] = v[1][:] a[:][j]
        // scalar product of vector v (row 1 of the 1 x m matrix) and column vector j of matrix a
        for (index_type j = simdBlockWidth; j < m; ++j)
        {
            w[i] += v[j] * aT(i, j);
        }
    }
}

#endif // GMX_SIMD

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
#if GMX_SIMD
    return matrixProductVectorMatrixTSimd(v, aT, w);
#else
    return matrixProductVectorMatrixTRef(v, aT, w);
#endif
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------

/*! \brief vector-matrix product \f$ \mathbf{w} = \mathbf{v}^{\mathrm{T}}\,\mathbf{A} \f$
           Currently, SIMD support is limited to float and/or double depending on the architecture.
           This template function is only active for unsupported data types and returns the reference
           implementation.

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   v        vector \f$ \mathbf{v} \f$
    \param[in]   a        matrix \f$ \mathbf{A} \f$
    \param[out]  w        the result vector \f$ \mathbf{w} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline
typename std::enable_if<!(gmx::SimdSetup<T>::doSimd_), void>::type
matrixProductVectorMatrixSimd(const std::vector<T, Alloc> &v,
                              const IrregArray2D<T, Alloc, false> &a,
                              std::vector<T, Alloc> &w)
{
    return matrixProductVectorMatrixRef(v, a, w);
}

#if GMX_SIMD // can't get rid of these guards because the functions load store etc are undefined if !GMX_SIMD

/*! \brief vector-matrix product \f$ \mathbf{w} = \mathbf{v}^{\mathrm{T}}\,\mathbf{A} \f$,
           SIMD implementation, only enabled for float or double, if the respective data
           type is supported by the GROMACS SIMD module on the current architecture.
           This variant also requires reduced-SIMD-width SIMD4 support.

           The data is expected to be aligned. In principle the alignment requirements
           could be matched against values provided by the allocator at compile time, but,
           since there is no standardized way for allocators to provide this information,
           the template function would not work for allocators that provide this information
           differently or not at all.

           \verbatim
            row vector    n column vectors    row vector
            of length n   of of length m      of length m
            -             - - - -
            vvvv vv       a a a a ...         wwww ww
                          a a a a ...
                          a a a a ...
                     x    a a a a ...    =
                          a a a a ...
                          a a a a ...
           \endverbatim

           use SIMD vector operations for the marked parts that match the SIMD width

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    allocator type for T used by the function arguments

    \param[in]   v        vector \f$ \mathbf{v} \f$
    \param[in]   a        matrix \f$ \mathbf{A} \f$
    \param[out]  w        the result vector \f$ \mathbf{w} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline
typename std::enable_if<gmx::SimdSetup<T>::doSimd_ &&gmx::SimdSetup<T>::doSimd4_, void>::type
matrixProductVectorMatrixSimd(const std::vector<T, Alloc> &v,
                              const IrregArray2D<T, Alloc, false> &a,
                              std::vector<T, Alloc> &w)
{
    typedef typename IrregArray2D<T, Alloc, false>::value_type value_type;
    typedef typename IrregArray2D<T, Alloc, false>::index_type index_type;

    // set up at compile time
    static constexpr index_type simdWidth  = SimdSetup<T>::simdWidth_;
    static constexpr index_type simd4Width = SimdSetup<T>::simd4Width_;
    typedef typename   SimdSetup<T>::simd_type     simd_float_type;
    typedef typename   SimdSetup<T>::simd4_type   simd4_float_type;

    // number of rows    of matrix A
    const index_type m = a.length1();
    // number of columns of matrix A
    const index_type n = a.length2(0);
    std::fill(w.begin(), w.end(), static_cast<value_type>(0));

    const index_type simdBlockWidth = n - n % simdWidth;
    if (simdBlockWidth != 0)
    {
        // stores stripes of the row vector v
        simd_float_type vSimd(0);
        // stores stripes of column vectors of a
        simd_float_type aSimd(0);
        // stores stripes of the result column vectors w
        simd_float_type wSimd(0);

        for (index_type j = 0; j < m; ++j)
        {
            // stores a replicated entry value of vector v
            vSimd = v[j];
            // scalar product of vector v (row 1 of the 1 x m matrix) and column vector j of matrix a
            for (index_type i = 0; i < simdBlockWidth; i += simdWidth)
            {
                wSimd = gmx::load<simd_float_type>(w.data() + i);
                aSimd = gmx::load<simd_float_type>(a(j) + i);
                wSimd = gmx::fma(vSimd, aSimd, wSimd);
                gmx::store(w.data() + i, wSimd);
            }
        }
    }

    const index_type rest            = n - simdBlockWidth;
    const index_type simd4BlockWidth = rest - rest % simd4Width;
    if (simd4BlockWidth != 0)
    {
        // stores stripes of the row vector v
        simd4_float_type vSimd4(0);
        // stores stripes of column vectors of a
        simd4_float_type aSimd4(0);
        // stores stripes of the result column vectors w
        simd4_float_type wSimd4(0);
        for (index_type j = 0; j < m; ++j)
        {
            // stores a replicated entry value of vector v
            vSimd4 = v[j];
            // scalar product of vector v (row 1 of the 1 x m matrix) and column vector j of matrix a
            for (index_type i = simdBlockWidth; i < (simdBlockWidth + simd4BlockWidth); i += simd4Width)
            {
                wSimd4 = gmx::load4(w.data() + i);
                aSimd4 = gmx::load4(a(j) + i);
                wSimd4 = gmx::fma(vSimd4, aSimd4, wSimd4);
                gmx::store4(w.data() + i, wSimd4);
            }
        }
    }

    // compute the remaining parts with scalar math
    for (index_type j = 0; j < m; ++j)
    {
        // w[1][i] = v[1][:] a[:][i]
        // scalar product of vector v (row 1 of the 1 x m matrix) and column vector j of matrix a
        for (index_type i = (simdBlockWidth + simd4BlockWidth); i < n; ++i)
        {
            w[i] += v[j] * a(j, i);
        }
    }
}

/*! \brief vector-matrix product \f$ \mathbf{w} = \mathbf{v}^{\mathrm{T}}\,\mathbf{A} \f$,
           SIMD implementation, only enabled for float or double, if the respective data
           type is supported by the GROMACS simd module on the current architecture.

           \verbatim
            row vector    n column vectors    row vector
            of length n   of of length m      of length m
            -             - - - -
            vvvv vv       a a a a ...         wwww ww
                          a a a a ...
                          a a a a ...
                     x    a a a a ...    =
                          a a a a ...
                          a a a a ...
           \endverbatim

           use SIMD vector operations for the marked parts that match the SIMD width

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    allocator type for T used by the function arguments

    \param[in]   v        vector \f$ \mathbf{v} \f$
    \param[in]   a        matrix \f$ \mathbf{A} \f$
    \param[out]  w        the result vector \f$ \mathbf{w} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline
typename std::enable_if<gmx::SimdSetup<T>::doSimd_ && !(gmx::SimdSetup<T>::doSimd4_), void>::type
matrixProductVectorMatrixSimd(const std::vector<T, Alloc> &v,
                              const IrregArray2D<T, Alloc, false> &a,
                              std::vector<T, Alloc> &w)
{
    typedef typename IrregArray2D<T, Alloc, false>::value_type value_type;
    typedef typename IrregArray2D<T, Alloc, false>::index_type index_type;

    // set up at compile time
    static constexpr index_type simdWidth  = SimdSetup<T>::simdWidth_;
    typedef typename   SimdSetup<T>::simd_type     simd_float_type;

    // number of rows    of matrix A
    const index_type m = a.length1();
    // number of columns of matrix A
    const index_type n = a.length2(0);
    std::fill(w.begin(), w.end(), static_cast<value_type>(0));

    const index_type simdBlockWidth = n - n % simdWidth;
    if (simdBlockWidth != 0)
    {
        // stores stripes of the row vector v
        simd_float_type vSimd(0);
        // stores stripes of column vectors of a
        simd_float_type aSimd(0);
        // stores stripes of the result column vectors w
        simd_float_type wSimd(0);
        for (index_type j = 0; j < m; ++j)
        {
            // stores a replicated entry value of vector v
            vSimd = v[j];
            // scalar product of vector v (row 1 of the 1 x m matrix) and column vector j of matrix a
            for (index_type i = 0; i < simdBlockWidth; i += simdWidth)
            {
                wSimd = gmx::load<simd_float_type>(w.data() + i);
                aSimd = gmx::load<simd_float_type>(a(j) + i);
                wSimd = gmx::fma(vSimd, aSimd, wSimd);
                gmx::store(w.data() + i, wSimd);
            }
        }
    }

    // compute the remaining parts with scalar math
    for (index_type j = 0; j < m; ++j)
    {
        // w[1][i] = v[1][:] a[:][i]
        // scalar product of vector v (row 1 of the 1 x m matrix) and column vector j of matrix a
        for (index_type i = simdBlockWidth; i < n; ++i)
        {
            w[i] += v[j] * a(j, i);
        }
    }
}

#endif // GMX_SIMD

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
#if GMX_SIMD
    return matrixProductVectorMatrixSimd(v, a, w);
#else
    return matrixProductVectorMatrixRef(v, a, w);
#endif
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------

/*! \brief vector-matrix product \f$ \mathbf{w} = \mathbf{A}\,\mathbf{v} \f$
           Currently, SIMD support is limited to float and/or double depending on the architecture.
           This template function is only active for unsupported data types and returns the reference
           implementation.

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   aT       matrix \f$ \mathbf{A}^{\mathrm{T}} \f$
    \param[in]   v        vector \f$ \mathbf{v} \f$
    \param[out]  w        the result vector \f$ \mathbf{w} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline
typename std::enable_if<!(gmx::SimdSetup<T>::doSimd_), void>::type
matrixProductMatrixTVectorSimd(const IrregArray2D<T, Alloc, false> &aT,
                               const std::vector<T, Alloc> &v,
                               std::vector<T, Alloc> &w)
{
    return matrixProductMatrixTVectorRef(aT, v, w);
}

#if GMX_SIMD // can't get rid of these guards because the functions load store etc are undefined if !GMX_SIMD

/*! \brief vector-matrix product \f$ \mathbf{w} = \mathbf{v}^{\mathrm{T}}\,\mathbf{A} \f$,
           SIMD implementation, only enabled for float or double, if the respective data
           type is supported by the GROMACS SIMD module on the current architecture.
           This variant also requires reduced-SIMD-width SIMD4 support.

           \verbatim
            row vector    n column vectors    row vector
            of length n   of of length m      of length m
            -             - - - -             ----
            vvvv vv       a a a a ...         wwww ww
                          a a a a ...
                          a a a a ...
                     x    a a a a ...    =
                          a a a a ...
                          a a a a ...
           \endverbatim

           use SIMD vector operations for the marked parts that match the SIMD width

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    allocator type for T used by the function arguments

    \param[in]   aT       matrix \f$ \mathbf{A}^{\mathrm{T}} \f$
    \param[in]   v        vector \f$ \mathbf{v} \f$
    \param[out]  w        the result vector \f$ \mathbf{w} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline
typename std::enable_if<gmx::SimdSetup<T>::doSimd_ &&gmx::SimdSetup<T>::doSimd4_, void>::type
matrixProductMatrixTVectorSimd(const IrregArray2D<T, Alloc, false> &aT,
                               const std::vector<T, Alloc> &v,
                               std::vector<T, Alloc> &w)
{
    typedef typename IrregArray2D<T, Alloc, false>::value_type value_type;
    typedef typename IrregArray2D<T, Alloc, false>::index_type index_type;

    // set up at compile time
    static constexpr index_type simdWidth  = SimdSetup<T>::simdWidth_;
    static constexpr index_type simd4Width = SimdSetup<T>::simd4Width_;
    typedef typename   SimdSetup<T>::simd_type     simd_float_type;
    typedef typename   SimdSetup<T>::simd4_type   simd4_float_type;

    // number of rows    of matrix A, columns of A^T
    const index_type m = aT.length2(0);
    // number of columns of matrix A, rows    of A^T
    const index_type n = aT.length1();

    std::fill(w.begin(), w.end(), static_cast<value_type>(0));

    const index_type simdBlockWidth = m - m % simdWidth;
    if (simdBlockWidth != 0)
    {
        // stores stripes of the row vector v
        simd_float_type vSimd(0);
        // stores stripes of column vectors of a
        simd_float_type aSimd(0);
        // stores stripes of the result column vectors w
        simd_float_type wSimd(0);
        for (index_type j = 0; j < n; ++j)
        {
            // stores a replicated entry value of vector v
            vSimd = v[j];
            // scalar product of vector v (row 1 of the 1 x m matrix) and column vector j of matrix a
            for (index_type i = 0; i < simdBlockWidth; i += simdWidth)
            {
                wSimd = gmx::load<simd_float_type>(w.data() + i);
                aSimd = gmx::load<simd_float_type>(aT(j) + i);
                wSimd = gmx::fma(vSimd, aSimd, wSimd);
                gmx::store(w.data() + i, wSimd);
            }
        }
    }

    const index_type rest            = m - simdBlockWidth;
    const index_type simd4BlockWidth = rest - rest % simd4Width;
    if (simd4BlockWidth != 0)
    {
        // stores stripes of the row vector v
        simd4_float_type vSimd4(0);
        // stores stripes of column vectors of a
        simd4_float_type aSimd4(0);
        // stores stripes of the result column vectors w
        simd4_float_type wSimd4(0);
        for (index_type j = 0; j < n; ++j)
        {
            // stores a replicated entry value of vector v
            vSimd4 = v[j];
            // scalar product of vector v (row 1 of the 1 x m matrix) and column vector j of matrix a
            for (index_type i = simdBlockWidth; i < (simdBlockWidth + simd4BlockWidth); i += simd4Width)
            {
                wSimd4 = gmx::load4(w.data() + i);
                aSimd4 = gmx::load4(aT(j) + i);
                wSimd4 = gmx::fma(vSimd4, aSimd4, wSimd4);
                gmx::store4(w.data() + i, wSimd4);
            }
        }
    }

    // compute the remaining parts with scalar math
    for (index_type i = (simdBlockWidth + simd4BlockWidth); i < m; ++i)
    {
        // w[1][i] = v[1][:] a[:][i]
        // scalar product of vector v (row 1 of the 1 x m matrix) and column vector j of matrix a
        for (index_type j = 0; j < n; ++j)
        {
            w[i] += aT(j, i) * v[j];
        }
    }
}

/*! \brief vector-matrix product \f$ \mathbf{w} = \mathbf{A}\,\mathbf{v} \f$,
           SIMD implementation, only enabled for float or double, if the respective data
           type is supported by the GROMACS simd module on the current architecture.

           \verbatim
            row vector    n column vectors    row vector
            of length n   of of length m      of length m
            -             - - - -             ----
            vvvv vv       a a a a ...         wwww ww
                          a a a a ...
                          a a a a ...
                     x    a a a a ...    =
                          a a a a ...
                          a a a a ...
           \endverbatim

           use SIMD vector operations for the marked parts that match the SIMD width

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    allocator type for T used by the function arguments

    \param[in]   aT       matrix \f$ \mathbf{A}^{\mathrm{T}} \f$
    \param[in]   v        vector \f$ \mathbf{v} \f$
    \param[out]  w        the result vector \f$ \mathbf{w} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline
typename std::enable_if<gmx::SimdSetup<T>::doSimd_ && !(gmx::SimdSetup<T>::doSimd4_), void>::type
matrixProductMatrixTVectorSimd(const IrregArray2D<T, Alloc, false> &aT,
                               const std::vector<T, Alloc> &v,
                               std::vector<T, Alloc> &w)
{
    typedef typename IrregArray2D<T, Alloc, false>::value_type value_type;
    typedef typename IrregArray2D<T, Alloc, false>::index_type index_type;

    // set up at compile time
    static constexpr index_type simdWidth  = SimdSetup<T>::simdWidth_;
    typedef typename   SimdSetup<T>::simd_type     simd_float_type;

    // number of rows    of matrix A, columns of A^T
    const index_type m = aT.length2(0);
    // number of columns of matrix A, rows    of A^T
    const index_type n = aT.length1();

    std::fill(w.begin(), w.end(), static_cast<value_type>(0));

    const index_type simdBlockWidth = m - m % simdWidth;
    if (simdBlockWidth != 0)
    {
        // stores stripes of the row vector v
        simd_float_type vSimd(0);
        // stores stripes of column vectors of a
        simd_float_type aSimd(0);
        // stores stripes of the result column vectors w
        simd_float_type wSimd(0);
        for (index_type j = 0; j < n; ++j)
        {
            // stores a replicated entry value of vector v
            vSimd = v[j];
            // scalar product of vector v (row 1 of the 1 x m matrix) and column vector j of matrix a
            for (index_type i = 0; i < simdBlockWidth; i += simdWidth)
            {
                wSimd = gmx::load<simd_float_type>(w.data() + i);
                aSimd = gmx::load<simd_float_type>(aT(j) + i);
                wSimd = gmx::fma(vSimd, aSimd, wSimd);
                gmx::store(w.data() + i, wSimd);
            }
        }
    }

    // compute the remaining parts with scalar math
    for (index_type i = simdBlockWidth; i < m; ++i)
    {
        // w[1][i] = v[1][:] a[:][i]
        // scalar product of vector v (row 1 of the 1 x m matrix) and column vector j of matrix a
        for (index_type j = 0; j < n; ++j)
        {
            w[i] += aT(j, i) * v[j];
        }
    }
}

#endif // GMX_SIMD

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
#if GMX_SIMD
    return matrixProductMatrixTVectorSimd(aT, v, w);
#else
    return matrixProductMatrixTVectorRef(aT, v, w);
#endif
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------

/*! \brief matrix product \f$ \mathbf{C} = \mathbf{A}\,\mathbf{B} \f$
           Currently, SIMD support is limited to float and/or double depending on the architecture.
           This template function is only active for unsupported data types and returns the reference
           implementation.

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   a        \f$ m \times l \f$ matrix \f$ \mathbf{A} \f$
    \param[in]   b        \f$ l \times n \f$ matrix \f$ \mathbf{B} \f$
    \param[out]  c        the result: \f$ m \times  n \f$ matrix \f$ \mathbf{C} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline
typename std::enable_if<!(gmx::SimdSetup<T>::doSimd_), void>::type
matrixProductSimd(const IrregArray2D<T, Alloc, false> &a,
                  const IrregArray2D<T, Alloc, false> &b,
                  IrregArray2D<T, Alloc, false> &c)
{
    return matrixProductRef(a, b, c);
}

#if GMX_SIMD // can't get rid of these guards because the functions load store etc are undefined if !GMX_SIMD

/*! \brief matrix product \f$ \mathbf{C} = \mathbf{A}\,\mathbf{B} \f$,
           SIMD implementation with SIMD4 support

           The function is only enabled for float or double, if the respective data
           type is supported by the GROMACS SIMD module on the current architecture.
           This variant also requires reduced-SIMD-width SIMD4 support.
           The data is expected to be aligned.

           \verbatim
            m row vectors    n column vectors    m rows x n columns
            of length l      of length l         m row    vectors of length n
                                                 n column vectors of length m

            -                    - - - -             - - - -
            a a a a a a          b b b b ...         c c c c . .
            a a a a a a          b b b b ...         c .
            a a a a a a          b b b b ...         c   .
            a a a a a a    x     b b b b ...    =    c     .
            . . . . . .          b b b b ...         c       .
            . . . . . .          b b b b ...         c . . . . .
           \endverbatim

           use SIMD vector operations for the marked parts that match the SIMD width

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   a        \f$ m \times l \f$ matrix \f$ \mathbf{A} \f$
    \param[in]   b        \f$ l \times n \f$ matrix \f$ \mathbf{B} \f$
    \param[out]  c        the result: \f$ m \times  n \f$ matrix \f$ \mathbf{C} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline
typename std::enable_if<gmx::SimdSetup<T>::doSimd_ &&gmx::SimdSetup<T>::doSimd4_, void>::type
matrixProductSimd(const IrregArray2D<T, Alloc, false> &a,
                  const IrregArray2D<T, Alloc, false> &b,
                  IrregArray2D<T, Alloc, false> &c)
{
    typedef typename IrregArray2D<T, Alloc, false>::value_type value_type;
    typedef typename IrregArray2D<T, Alloc, false>::index_type index_type;

    // set up at compile time
    static constexpr index_type simdWidth  = SimdSetup<T>::simdWidth_;
    static constexpr index_type simd4Width = SimdSetup<T>::simd4Width_;
    typedef typename   SimdSetup<T>::simd_type     simd_float_type;
    typedef typename   SimdSetup<T>::simd4_type   simd4_float_type;

    // number of rows    of matrix A
    const index_type m = a.length1();
    // number of columns of matrix A/rows of matrix B
    const index_type l = a.length2(0);
    // number of columns of matrix B
    const index_type n = b.length2(0);

    // initialize matrix c one row vector at a time
    for (index_type i = 0; i < m; ++i)
    {
        std::fill(c(i), c(i) + static_cast<ptrdiff_t>(n), static_cast<value_type>(0));
    }

    const index_type simdBlockWidth = n - n % simdWidth;
    if (simdBlockWidth != 0)
    {
        // stores replicated elements of the row vectors of a
        simd_float_type aSimd(0);
        // stores stripes of column vectors of b
        simd_float_type bSimd(0);
        // stores stripes of the result column vectors c
        simd_float_type cSimd(0);

        // Ideally, we would directly compute c[i][j] from the scalar product of
        // row vector a[i][:] and column vector b[j][:]. But SIMD loads/stores
        // of partial vectors is only possible along rows. Therefore, we compute
        // parts of multiple elements along rows of c in each vector operation.
        // c[i][j0..j1] += a[i][k] * b[k][j0..j1] only sequential reading
        for (index_type i = 0; i < m; ++i)
        {
            for (index_type k = 0; k < l; ++k)
            {
                // stores a replicated entry value of matrix a
                aSimd = a(i, k);
                for (index_type j = 0; j < simdBlockWidth; j += simdWidth)
                {
                    // stores a replicated entry value of the column vector vb
                    bSimd = gmx::load<simd_float_type>(b(k) + j);
                    cSimd = gmx::load<simd_float_type>(c(i) + j);
                    cSimd = gmx::fma(aSimd, bSimd, cSimd);
                    gmx::store(c(i) + j, cSimd);
                }
            }
        }
    }

    const index_type rest               = n - simdBlockWidth;
    const index_type simd4BlockWidth    = rest - rest % simd4Width;
    const index_type simdFullBlockWidth = simdBlockWidth + simd4BlockWidth;
    if (simd4BlockWidth != 0)
    {
        // stores replicated elements of the row vectors of a
        simd4_float_type aSimd4(0);
        // stores stripes of column vectors of b
        simd4_float_type bSimd4(0);
        // stores stripes of the result column vectors c
        simd4_float_type cSimd4(0);

        for (index_type i = 0; i < m; ++i)
        {
            for (index_type k = 0; k < l; ++k)
            {
                // stores a replicated entry value of matrix a
                aSimd4 = a(i, k);
                for (index_type j = simdBlockWidth; j < simdFullBlockWidth; j += simd4Width)
                {
                    bSimd4 = gmx::load4(b(k) + j);
                    cSimd4 = gmx::load4(c(i) + j);
                    cSimd4 = gmx::fma(aSimd4, bSimd4, cSimd4);
                    gmx::store4(c(i) + j, cSimd4);
                }
            }
        }
    }

    // compute the remaining parts with scalar math
    for (index_type i = 0; i < m; ++i)
    {
        for (index_type k = 0; k < l; ++k)
        {
            for (index_type j = simdFullBlockWidth; j < n; ++j)
            {
                c(i, j) += a(i, k) * b(k, j);
            }
        }
    }
}

/*! \brief matrix product \f$ \mathbf{C} = \mathbf{A}\,\mathbf{B} \f$,
           SIMD implementation

           The function is only enabled for float or double, if the respective data
           type is supported by the GROMACS SIMD module on the current architecture.
           The data is expected to be aligned.

           \verbatim
            m row vectors    n column vectors    m rows x n columns
            of length l      of length l         m row    vectors of length n
                                                 n column vectors of length m

            -                    - - - -             - - - -
            a a a a a a          b b b b ...         c c c c . .
            a a a a a a          b b b b ...         c .
            a a a a a a          b b b b ...         c   .
            a a a a a a    x     b b b b ...    =    c     .
            . . . . . .          b b b b ...         c       .
            . . . . . .          b b b b ...         c . . . . .
           \endverbatim

           use SIMD vector operations for the marked parts that match the SIMD width
           Other variants of using SIMD operations for computing the matrix product
           are conceivable, but this turned out to work particularly well.

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   a        \f$ m \times l \f$ matrix \f$ \mathbf{A} \f$
    \param[in]   b        \f$ l \times n \f$ matrix \f$ \mathbf{B} \f$
    \param[out]  c        the result: \f$ m \times  n \f$ matrix \f$ \mathbf{C} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline
typename std::enable_if<gmx::SimdSetup<T>::doSimd_ && !(gmx::SimdSetup<T>::doSimd4_), void>::type
matrixProductSimd(const IrregArray2D<T, Alloc, false> &a,
                  const IrregArray2D<T, Alloc, false> &b,
                  IrregArray2D<T, Alloc, false> &c)
{
    typedef typename IrregArray2D<T, Alloc, false>::value_type value_type;
    typedef typename IrregArray2D<T, Alloc, false>::index_type index_type;

    // set up at compile time
    static constexpr index_type simdWidth  = SimdSetup<T>::simdWidth_;
    typedef typename   SimdSetup<T>::simd_type     simd_float_type;

    // number of rows    of matrix A
    const index_type m = a.length1();
    // number of columns of matrix A/rows of matrix B
    const index_type l = a.length2(0);
    // number of columns of matrix B
    const index_type n = b.length2(0);

    // initialize matrix c one row vector at a time
    for (index_type i = 0; i < m; ++i)
    {
        std::fill(c(i), c(i) + static_cast<ptrdiff_t>(n), static_cast<value_type>(0));
    }

    const index_type simdBlockWidth = n - n % simdWidth;

    if (simdBlockWidth != 0)
    {
        // stores replicated elements of the row vectors of a
        simd_float_type aSimd(0);
        // stores stripes of column vectors of b
        simd_float_type bSimd(0);
        // stores stripes of the result column vectors c
        simd_float_type cSimd(0);

        // Ideally, we would directly compute c[i][j] from the scalar product of
        // row vector a[i][:] and column vector b[j][:]. But SIMD loads/stores
        // of partial vectors is only possible along rows. Therefore, we compute
        // parts of multiple elements along rows of c in each vector operation.
        // c[i][j0..j1] += a[i][k] * b[k][j0..j1] only sequential reading
        for (index_type i = 0; i < m; ++i)
        {
            for (index_type k = 0; k < l; ++k)
            {
                // stores a replicated entry value of matrix a
                aSimd = a(i, k);
                for (index_type j = 0; j < simdBlockWidth; j += simdWidth)
                {
                    // stores a replicated entry value of the column vector vb
                    bSimd = gmx::load<simd_float_type>(b(k) + j);
                    cSimd = gmx::load<simd_float_type>(c(i) + j);
                    cSimd = gmx::fma(aSimd, bSimd, cSimd);
                    gmx::store(c(i) + j, cSimd);
                }
            }
        }
    }

    // compute the remaining parts with scalar math
    for (index_type i = 0; i < m; ++i)
    {
        for (index_type k = 0; k < l; ++k)
        {
            for (index_type j = simdBlockWidth; j < n; ++j)
            {
                c(i, j) += a(i, k) * b(k, j);
            }
        }
    }
}

#endif // GMX_SIMD

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
#if GMX_SIMD
    return matrixProductSimd(a, b, c);
#else
    return matrixProductRef(a, b, c);
#endif
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

/*! \brief matrix product \f$ \mathbf{C} = \mathbf{A}\,\mathbf{B} \f$
           Currently, SIMD support is limited to float and/or double depending on the architecture.
           This template function is only active for unsupported data types and returns the reference
           implementation.

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   a        \f$ m \times l \f$ matrix \f$ \mathbf{A} \f$
    \param[in]   bT       \f$ n \times l \f$ matrix \f$ \mathbf{B}^{\mathrm{T}} \f$
    \param[out]  c        the result: \f$ m \times  n \f$ matrix \f$ \mathbf{C} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline
typename std::enable_if<!(gmx::SimdSetup<T>::doSimd_), void>::type
matrixProductMatrixMatrixTSimd(const IrregArray2D<T, Alloc, false> &a,
                               const IrregArray2D<T, Alloc, false> &bT,
                               IrregArray2D<T, Alloc, false> &c)
{
    return matrixProductMatrixMatrixTRef(a, bT, c);
}

#if GMX_SIMD // can't get rid of these guards because the functions load store etc are undefined if !GMX_SIMD

/*! \brief matrix product \f$ \mathbf{C} = \mathbf{A}\,\mathbf{B} \f$,
           SIMD implementation with SIMD4 support

           The function is only enabled for float or double, if the respective data
           type is supported by the GROMACS SIMD module on the current architecture.
           This variant also requires reduced-SIMD-width SIMD4 support.
           The data is expected to be aligned.

           \verbatim
            m row vectors    n column vectors    m rows x n columns
            of length l      of length l         m row    vectors of length n
                                                 n column vectors of length m

            - - - -              - - - -             - - - -
            a a a a . .          b b b b . .         c c c c . .
            a a a a . .          b b b b . .         c .
            a a a a . .          b b b b . .         c   .
            a a a a . .    x     b b b b . .    =    c     .
            . . . . . .          . . . . . .         c       .
            . . . . . .          . . . . . .         c . . . . .
           \endverbatim

           use SIMD vector operations for the marked parts that match the SIMD width

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   a        \f$ m \times l \f$ matrix \f$ \mathbf{A} \f$
    \param[in]   bT       \f$ n \times l \f$ matrix \f$ \mathbf{B}^{\mathrm{T}} \f$
    \param[out]  c        the result: \f$ m \times  n \f$ matrix \f$ \mathbf{C} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline
typename std::enable_if<gmx::SimdSetup<T>::doSimd_ &&gmx::SimdSetup<T>::doSimd4_, void>::type
matrixProductMatrixMatrixTSimd(const IrregArray2D<T, Alloc, false> &a,
                               const IrregArray2D<T, Alloc, false> &bT,
                               IrregArray2D<T, Alloc, false> &c)
{
    typedef typename IrregArray2D<T, Alloc, false>::value_type value_type;
    typedef typename IrregArray2D<T, Alloc, false>::index_type index_type;

    // set up at compile time
    static constexpr index_type simdWidth  = SimdSetup<T>::simdWidth_;
    static constexpr index_type simd4Width = SimdSetup<T>::simd4Width_;
    typedef typename   SimdSetup<T>::simd_type     simd_float_type;
    typedef typename   SimdSetup<T>::simd4_type   simd4_float_type;

    // number of rows    of matrix A
    const index_type m = a.length1();
    // number of columns of matrix A/rows of matrix B
    const index_type l = a.length2(0);
    // number of columns of matrix B
    const index_type n = bT.length1();

    // initialize matrix c one row vector at a time
    for (index_type i = 0; i < m; ++i)
    {
        std::fill(c(i), c(i) + static_cast<ptrdiff_t>(n), static_cast<value_type>(0));
    }

    const index_type simdBlockWidth = l - l % simdWidth;
    if (simdBlockWidth != 0)
    {
        // stores replicated elements of the row vectors of a
        simd_float_type aSimd(0);
        // stores stripes of column vectors of b
        simd_float_type bSimd(0);
        // stores stripes of the result column vectors c
        simd_float_type cSimd(0);

        // Since we have the transpose of B, reading along rows of A and columns of B
        // is efficient. Thus contributions to a single element of C can be summed up.
        for (index_type i = 0; i < m; ++i)
        {
            for (index_type j = 0; j < n; ++j)
            {
                // stores a replicated entry value of matrix a
                for (index_type k = 0; k < simdBlockWidth; k += simdWidth)
                {
                    // store matching entry stretches of the column vector vb and row vector va
                    aSimd = gmx::load<simd_float_type>(a(i) + k);
                    bSimd = gmx::load<simd_float_type>(bT(j) + k);
                    cSimd = gmx::fma(aSimd, bSimd, cSimd);
                }
                c(i, j) = gmx::reduce(cSimd);
                cSimd   = gmx::setZero();
            }
        }
    }

    const index_type rest               = l - simdBlockWidth;
    const index_type simd4BlockWidth    = rest - rest % simd4Width;
    const index_type simdFullBlockWidth = simdBlockWidth + simd4BlockWidth;
    if (simd4BlockWidth != 0)
    {
        // stores replicated elements of the row vectors of a
        simd4_float_type aSimd4(0);
        // stores stripes of column vectors of b
        simd4_float_type bSimd4(0);
        // stores stripes of the result column vectors c
        simd4_float_type cSimd4(0);

        // Since we have the transpose of B, reading along rows of A and columns of B
        // is efficient. Thus contributions to a single element of C can be summed up.
        for (index_type i = 0; i < m; ++i)
        {
            for (index_type j = 0; j < n; ++j)
            {
                for (index_type k = simdBlockWidth; k < simdFullBlockWidth; k += simd4Width)
                {
                    // store matching entry stretches of the column vector vb and row vector va
                    aSimd4 = gmx::load4(a(i) + k);
                    bSimd4 = gmx::load4(bT(j) + k);
                    cSimd4 = gmx::fma(aSimd4, bSimd4, cSimd4);
                }
                c(i, j) += gmx::reduce(cSimd4);
                cSimd4   = gmx::setZero();
            }
        }
    }

    // compute the remaining parts with scalar math
    for (index_type i = 0; i < m; ++i)
    {
        for (index_type j = 0; j < n; ++j)
        {
            for (index_type k = simdFullBlockWidth; k < l; ++k)
            {
                c(i, j) += a(i, k) * bT(j, k);
            }
        }
    }
}

/*! \brief matrix product \f$ \mathbf{C} = \mathbf{A}\,\mathbf{B} \f$,
           SIMD implementation

           The function is only enabled for float or double, if the respective data
           type is supported by the GROMACS SIMD module on the current architecture.
           The data is expected to be aligned.

           \verbatim
            m row vectors    n column vectors    m rows x n columns
            of length l      of length l         m row    vectors of length n
                                                 n column vectors of length m

            - - - -              - - - -             - - - -
            a a a a . .          b b b b . .         c c c c . .
            a a a a . .          b b b b . .         c .
            a a a a . .          b b b b . .         c   .
            a a a a . .    x     b b b b . .    =    c     .
            . . . . . .          . . . . . .         c       .
            . . . . . .          . . . . . .         c . . . . .
           \endverbatim

           use SIMD vector operations for the marked parts that match the SIMD width
           Other variants of using SIMD operations for computing the matrix product
           are conceivable, but this turned out to work well.

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   a        \f$ m \times l \f$ matrix \f$ \mathbf{A} \f$
    \param[in]   bT       \f$ n \times l \f$ matrix \f$ \mathbf{B}^{\mathrm{T}} \f$
    \param[out]  c        the result: \f$ m \times  n \f$ matrix \f$ \mathbf{C} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline
typename std::enable_if<gmx::SimdSetup<T>::doSimd_ && !(gmx::SimdSetup<T>::doSimd4_), void>::type
matrixProductMatrixMatrixTSimd(const IrregArray2D<T, Alloc, false> &a,
                               const IrregArray2D<T, Alloc, false> &bT,
                               IrregArray2D<T, Alloc, false> &c)
{
    typedef typename IrregArray2D<T, Alloc, false>::value_type value_type;
    typedef typename IrregArray2D<T, Alloc, false>::index_type index_type;

    // set up at compile time
    static constexpr index_type simdWidth  = SimdSetup<T>::simdWidth_;
    typedef typename   SimdSetup<T>::simd_type     simd_float_type;

    // number of rows    of matrix A
    const index_type m = a.length1();
    // number of columns of matrix A/rows of matrix B
    const index_type l = a.length2(0);
    // number of columns of matrix B
    const index_type n = bT.length1();

    // initialize matrix c one row vector at a time
    for (index_type i = 0; i < m; ++i)
    {
        std::fill(c(i), c(i) + static_cast<ptrdiff_t>(n), static_cast<value_type>(0));
    }

    const index_type simdBlockWidth = l - l % simdWidth;
    if (simdBlockWidth != 0)
    {
        // stores replicated elements of the row vectors of a
        simd_float_type aSimd(0);
        // stores stripes of column vectors of b
        simd_float_type bSimd(0);
        // stores stripes of the result column vectors c
        simd_float_type cSimd(0);

        // Since we have the transpose of B, reading along rows of A and columns of B
        // is efficient. Thus contributions to a single element of C can be summed up.
        for (index_type i = 0; i < m; ++i)
        {
            for (index_type j = 0; j < n; ++j)
            {
                for (index_type k = 0; k < simdBlockWidth; k += simdWidth)
                {
                    // store matching entry stretches of the column vector vb and row vector va
                    aSimd = gmx::load<simd_float_type>(a(i) + k);
                    bSimd = gmx::load<simd_float_type>(bT(j) + k);
                    cSimd = gmx::fma(aSimd, bSimd, cSimd);
                }
                c(i, j) = gmx::reduce(cSimd);
                cSimd   = gmx::setZero();
            }
        }
    }

    // compute the remaining parts with scalar math
    for (index_type i = 0; i < m; ++i)
    {
        for (index_type j = 0; j < n; ++j)
        {
            for (index_type k = simdBlockWidth; k < l; ++k)
            {
                c(i, j) += a(i, k) * bT(j, k);
            }
        }
    }
}

#endif // GMX_SIMD

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
#if GMX_SIMD
    return matrixProductMatrixMatrixTSimd(a, bT, c);
#else
    return matrixProductMatrixMatrixTRef(a, bT, c);
#endif
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------

/*! \brief matrix product \f$ \mathbf{C} = \mathbf{A}\,\mathbf{B} \f$
           Currently, SIMD support is limited to float and/or double depending on the architecture.
           This template function is only active for unsupported data types and returns the reference
           implementation.

    \note  Is there a good way of SIMDzation for this data layout, i.e., C = A^T B^T? If so please add.
           Temporarily transposing A^T requires transient memory allocation, but allows to return
           C = A B^T.

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   aT       \f$ l \times m \f$ matrix \f$ \mathbf{A}^{\mathrm{T}} \f$
    \param[in]   bT       \f$ n \times l \f$ matrix \f$ \mathbf{B}^{\mathrm{T}} \f$
    \param[out]  c        the result: \f$ m \times  n \f$ matrix \f$ \mathbf{C} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline void
//typename std::enable_if<!(gmx::SimdSetup<T>::doSimd_), void>::type
matrixProductMatrixTMatrixTSimd(const IrregArray2D<T, Alloc, false> &aT,
                                const IrregArray2D<T, Alloc, false> &bT,
                                IrregArray2D<T, Alloc, false> &c)
{
    typedef IrregArray2D<T, Alloc, false> array_type;
    typedef typename array_type::index_type index_type;

    // number of rows    of matrix A^T
    const index_type     l = aT.length1();
    // number of columns of matrix A^T
    const index_type     m = aT.length2(0);

    constexpr index_type cutNElements = 1000;
    if (l * m < cutNElements)
    {
        return matrixProductMatrixTMatrixTRef(aT, bT, c);
    }
    else
    {
        array_type a(m, l);

        for (index_type i = 0; i < m; ++i)
        {
            for (index_type j = 0; j < l; ++j)
            {
                a(i, j) = aT(j, i);
            }
        }

        return matrixProductMatrixMatrixTSimd(a, bT, c);
    }
}


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
#if GMX_SIMD
    return matrixProductMatrixTMatrixTSimd(aT, bT, c);
#else
    return matrixProductMatrixTMatrixTRef(aT, bT, c);
#endif
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------

/*! \brief matrix product \f$ \mathbf{C} = \mathbf{A}\,\mathbf{B} \f$
           Currently, SIMD support is limited to float and/or double depending on the architecture.
           This template function is only active for unsupported data types and returns the reference
           implementation.

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   aT       \f$ l \times m \f$ matrix \f$ \mathbf{A}^{\mathrm{T}} \f$
    \param[in]   b        \f$ l \times n \f$ matrix \f$ \mathbf{B} \f$
    \param[out]  c        the result: \f$ m \times  n \f$ matrix \f$ \mathbf{C} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline
typename std::enable_if<!(gmx::SimdSetup<T>::doSimd_), void>::type
matrixProductMatrixTMatrixSimd(const IrregArray2D<T, Alloc, false> &aT,
                               const IrregArray2D<T, Alloc, false> &b,
                               IrregArray2D<T, Alloc, false> &c)
{
    return matrixProductMatrixTMatrixRef(aT, b, c);
}

#if GMX_SIMD // can't get rid of these guards because the functions load store etc are undefined if !GMX_SIMD

/*! \brief matrix product \f$ \mathbf{C} = \mathbf{A}\,\mathbf{B} \f$,
           SIMD implementation with SIMD4 support

           The function is only enabled for float or double, if the respective data
           type is supported by the GROMACS SIMD module on the current architecture.
           This variant also requires reduced-SIMD-width SIMD4 support.
           The data is expected to be aligned.

           \verbatim
            m row vectors    n column vectors    m rows x n columns
            of length l      of length l         m row    vectors of length n
                                                 n column vectors of length m

                                 - - - -             - - - -
            a a a a a a |        b b b b ...         c c c c . .
            a a a a a a          b b b b ...         c .
            a a a a a a          b b b b ...         c   .
            a a a a a a    x     b b b b ...    =    c     .
            . . . . . .          b b b b ...         c       .
            . . . . . .          b b b b ...         c . . . . .
           \endverbatim

           use SIMD vector operations for the marked parts that match the SIMD width

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   aT       \f$ l \times m \f$ matrix \f$ \mathbf{A}^{\mathrm{T}} \f$
    \param[in]   b        \f$ l \times n \f$ matrix \f$ \mathbf{B} \f$
    \param[out]  c        the result: \f$ m \times  n \f$ matrix \f$ \mathbf{C} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline
typename std::enable_if<gmx::SimdSetup<T>::doSimd_ &&gmx::SimdSetup<T>::doSimd4_, void>::type
matrixProductMatrixTMatrixSimd(const IrregArray2D<T, Alloc, false> &aT,
                               const IrregArray2D<T, Alloc, false> &b,
                               IrregArray2D<T, Alloc, false> &c)
{
    typedef typename IrregArray2D<T, Alloc, false>::value_type value_type;
    typedef typename IrregArray2D<T, Alloc, false>::index_type index_type;

    // set up at compile time
    static constexpr index_type simdWidth  = SimdSetup<T>::simdWidth_;
    static constexpr index_type simd4Width = SimdSetup<T>::simd4Width_;
    typedef typename   SimdSetup<T>::simd_type     simd_float_type;
    typedef typename   SimdSetup<T>::simd4_type   simd4_float_type;

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

    const index_type simdBlockWidth = n - n % simdWidth;
    if (simdBlockWidth != 0)
    {
        // stores replicated elements of the row vectors of a
        simd_float_type aSimd(0);
        // stores stripes of column vectors of b
        simd_float_type bSimd(0);
        // stores stripes of the result column vectors c
        simd_float_type cSimd(0);

        // Ideally, we would directly compute c[i][j] from the scalar product of
        // row vector a[i][:] and column vector b[j][:]. But SIMD loads/stores
        // of partial vectors is only possible along rows. Therefore, we compute
        // parts of multiple elements along rows of c in each vector operation.
        // c[i][j0..j1] += a[i][k] * b[k][j0..j1] only sequential reading
        for (index_type k = 0; k < l; ++k)
        {
            for (index_type i = 0; i < m; ++i)
            {
                // stores a replicated entry value of matrix a
                aSimd = aT(k, i);
                for (index_type j = 0; j < simdBlockWidth; j += simdWidth)
                {
                    // stores a replicated entry value of the column vector vb
                    bSimd = gmx::load<simd_float_type>(b(k) + j);
                    cSimd = gmx::load<simd_float_type>(c(i) + j);
                    cSimd = gmx::fma(aSimd, bSimd, cSimd);
                    gmx::store(c(i) + j, cSimd);
                }
            }
        }
    }

    const index_type rest               = n - simdBlockWidth;
    const index_type simd4BlockWidth    = rest - rest % simd4Width;
    const index_type simdFullBlockWidth = simdBlockWidth + simd4BlockWidth;
    if (simd4BlockWidth != 0)
    {
        // stores replicated elements of the row vectors of a
        simd4_float_type aSimd4(0);
        // stores stripes of column vectors of b
        simd4_float_type bSimd4(0);
        // stores stripes of the result column vectors c
        simd4_float_type cSimd4(0);

        for (index_type k = 0; k < l; ++k)
        {
            for (index_type i = 0; i < m; ++i)
            {
                // stores a replicated entry value of matrix a
                aSimd4 = aT(k, i);
                for (index_type j = simdBlockWidth; j < simdFullBlockWidth; j += simd4Width)
                {
                    bSimd4 = gmx::load4(b(k) + j);
                    cSimd4 = gmx::load4(c(i) + j);
                    cSimd4 = gmx::fma(aSimd4, bSimd4, cSimd4);
                    gmx::store4(c(i) + j, cSimd4);
                }
            }
        }
    }

    // compute the remaining parts with scalar math
    for (index_type k = 0; k < l; ++k)
    {
        for (index_type i = 0; i < m; ++i)
        {
            for (index_type j = simdFullBlockWidth; j < n; ++j)
            {
                c(i, j) += aT(k, i) * b(k, j);
            }
        }
    }
}

/*! \brief matrix product \f$ \mathbf{C} = \mathbf{A}\,\mathbf{B} \f$,
           SIMD implementation

           The function is only enabled for float or double, if the respective data
           type is supported by the GROMACS SIMD module on the current architecture.
           This variant also requires reduced-SIMD-width SIMD4 support.
           The data is expected to be aligned.

           \verbatim
            m row vectors    n column vectors    m rows x n columns
            of length l      of length l         m row    vectors of length n
                                                 n column vectors of length m

            -                    - - - -             - - - -
            a a a a a a          b b b b ...         c c c c . .
            a a a a a a          b b b b ...         c .
            a a a a a a          b b b b ...         c   .
            a a a a a a    x     b b b b ...    =    c     .
            . . . . . .          b b b b ...         c       .
            . . . . . .          b b b b ...         c . . . . .
           \endverbatim

           use SIMD vector operations for the marked parts that match the SIMD width
           Other variants of using SIMD operations for computing the matrix product
           are conceivable, but this turned out to work particularly well.

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   aT       \f$ l \times m \f$ matrix \f$ \mathbf{A}^{\mathrm{T}} \f$
    \param[in]   b        \f$ l \times n \f$ matrix \f$ \mathbf{B} \f$
    \param[out]  c        the result: \f$ m \times  n \f$ matrix \f$ \mathbf{C} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline
typename std::enable_if<gmx::SimdSetup<T>::doSimd_ && !(gmx::SimdSetup<T>::doSimd4_), void>::type
matrixProductMatrixTMatrixSimd(const IrregArray2D<T, Alloc, false> &aT,
                               const IrregArray2D<T, Alloc, false> &b,
                               IrregArray2D<T, Alloc, false> &c)
{
    typedef typename IrregArray2D<T, Alloc, false>::value_type value_type;
    typedef typename IrregArray2D<T, Alloc, false>::index_type index_type;

    // set up at compile time
    static constexpr index_type simdWidth  = SimdSetup<T>::simdWidth_;
    typedef typename   SimdSetup<T>::simd_type     simd_float_type;

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

    const index_type simdBlockWidth = n - n % simdWidth;

    if (simdBlockWidth != 0)
    {
        // stores replicated elements of the row vectors of a
        simd_float_type aSimd(0);
        // stores stripes of column vectors of b
        simd_float_type bSimd(0);
        // stores stripes of the result column vectors c
        simd_float_type cSimd(0);

        // Ideally, we would directly compute c[i][j] from the scalar product of
        // row vector a[i][:] and column vector b[j][:]. But SIMD loads/stores
        // of partial vectors is only possible along rows. Therefore, we compute
        // parts of multiple elements along rows of c in each vector operation.
        // c[i][j0..j1] += a[i][k] * b[k][j0..j1] only sequential reading
        for (index_type k = 0; k < l; ++k)
        {
            for (index_type i = 0; i < m; ++i)
            {
                // stores a replicated entry value of matrix a
                aSimd = aT(k, i);
                for (index_type j = 0; j < simdBlockWidth; j += simdWidth)
                {
                    // stores a replicated entry value of the column vector vb
                    bSimd = gmx::load<simd_float_type>(b(k) + j);
                    cSimd = gmx::load<simd_float_type>(c(i) + j);
                    cSimd = gmx::fma(aSimd, bSimd, cSimd);
                    gmx::store(c(i) + j, cSimd);
                }
            }
        }
    }

    // compute the remaining parts with scalar math
    for (index_type k = 0; k < l; ++k)
    {
        for (index_type i = 0; i < m; ++i)
        {
            for (index_type j = simdBlockWidth; j < n; ++j)
            {
                c(i, j) += aT(k, i) * b(k, j);
            }
        }
    }
}

#endif // GMX_SIMD

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
#if GMX_SIMD
    return matrixProductMatrixTMatrixSimd(aT, b, c);
#else
    return matrixProductMatrixTMatrixRef(aT, b, c);
#endif
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------

} // end namespace gmx

#endif
