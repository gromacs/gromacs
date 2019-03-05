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
/*! \libinternal \file

    \brief   matrix inverse and determinant based on LUP decomposition utilizing IrregArray2D

    The functions are intentionally declared as free functions to allow
    storage of non-arithmetic data types in IrregArray2D, too.

   \note The organization of the functions is the same as explained
         in header irreg_array_2_matrix_products.h except for omission
         of the SIMD4 variant of the functions.

   \author R. Thomas Ullmann <tullman@gwdg.de>
   \inlibraryapi
   \ingroup module_math
 */
#ifndef GMX_MATH_IRREG_ARRAY_MATH_IRREG_ARRAY_2D_LU_DECOMPOSITION_H
#define GMX_MATH_IRREG_ARRAY_MATH_IRREG_ARRAY_2D_LU_DECOMPOSITION_H
#include <numeric>
#include <type_traits>

#include "gromacs/simd/simd_math.h"
#include "gromacs/simd/simd_setup.h"
#include "gromacs/utility/data_structures/irreg_array_2d.h"

namespace gmx
{

// --------------------------------------------------------------------------------------------------------------------------------------------------------

/*! \brief LU decomposition with row pivoting \f$ \mathbf{A} = \mathbf{P}^{-1}\,\mathbf{L}\,\mathbf{U} \f$
           Currently, SIMD support is limited to float and/or double depending on the architecture.
           This template function is only active for unsupported data types and returns the reference
           implementation.

    \ingroup module_math
    \libinternal

    \tparam      T           data type to be used in the calculation
    \tparam      Alloc       data type to be used in the calculation
    \tparam      TInt        integer data type to be used for indexing array elements
    \tparam      AllocInt    allocator for TInt

    \param[in]   a        \f$ n \times n\f$ matrix \f$ \mathbf{A} \f$ to be decomposed
    \param[out]  p        vector representation of the \f$ n \times n\f$ row pivot matrix \f$ \mathbf{P} \f$, entry \f$ i \f$ indicates the new row index of the current row \f$ i \f$
    \param[out]  lu       the lower and upper triangular matrices \f$ \mathbf{L} \f$ and \f$ \mathbf{U} \f$  stored in the upper right and lower left omitting the diagonal entries of \f$ \mathbf{L} \f$, which are all equal to \f$ 1 \f$
 */
template<typename T, class Alloc = std::allocator<T>, typename TInt, class AllocInt = std::allocator<TInt>  >
inline
typename std::enable_if<!(gmx::SimdSetup<T>::doSimd_), bool>::type
LUPDecompositionSimd(const IrregArray2D<T, Alloc, false> &a, std::vector<TInt, AllocInt> &p, IrregArray2D<T, Alloc, false> &lu)
{
    return LUPDecompositionRef(a, p, lu);
}

/*! \brief LU decomposition with row pivoting \f$ \mathbf{A} = \mathbf{P}^{-1}\,\mathbf{L}\,\mathbf{U} \f$,
           SIMD implementation, returns true if successful and false otherwise

    \ingroup module_math
    \libinternal

    \tparam      T           data type to be used in the calculation
    \tparam      Alloc       data type to be used in the calculation
    \tparam      TInt        integer data type to be used for indexing array elements
    \tparam      AllocInt    allocator for TInt

    \param[in]   a        \f$ n \times n\f$ matrix \f$ \mathbf{A} \f$ to be decomposed
    \param[out]  p        vector representation of the \f$ n \times n\f$ row pivot matrix \f$ \mathbf{P} \f$, entry \f$ i \f$ indicates the new row index of the current row \f$ i \f$
    \param[out]  lu       the lower and upper triangular matrices \f$ \mathbf{L} \f$ and \f$ \mathbf{U} \f$  stored in the upper right and lower left omitting the diagonal entries of \f$ \mathbf{L} \f$, which are all equal to \f$ 1 \f$
 */
template<typename T, class Alloc = std::allocator<T>, typename TInt, class AllocInt = std::allocator<TInt> >
inline
typename std::enable_if<gmx::SimdSetup<T>::doSimd_, bool>::type
LUPDecompositionSimd(const IrregArray2D<T, Alloc, false> &a, std::vector<TInt, AllocInt> &p, IrregArray2D<T, Alloc, false> &lu)
{
    typedef typename IrregArray2D<T, Alloc>::value_type value_type;
    typedef typename IrregArray2D<T, Alloc>::index_type index_type;
    typedef typename   SimdSetup<T>::simd_type     simd_float_type;

    // set up at compile time
    static constexpr index_type simdWidth  = SimdSetup<T>::simdWidth_;

    simd_float_type             aSimd(0);
    simd_float_type             bSimd(0);

    const index_type            n               = a.length1();
    const index_type            maxSimdBlockEnd = n - n % simdWidth;

    // initialize P such that each row keeps its original index
    std::iota(p.begin(), p.end(), static_cast<index_type>(0));

    // initialize LU as copy of A
    for (index_type i = 0; i < n; ++i)
    {
        std::copy(a(i), a(i) + n, lu(i));
    }

    // determine the pivot matrix by searching for the largest entry in each column
    // After row pivoting (reordering of rows), the first row has the largest entry
    // in column 1, the second row has the largest value in column 2 and so on.
    // The last element of a row will remain in each row during the decomposition.
    for (index_type j = 0; j < (n - 1); ++j)
    {
        value_type maxColVal         = 0;
        index_type maxColValRowIndex = 0;

        for (index_type i = j; i < n; ++i)
        {
            const value_type tmpf = std::abs(lu(i, j));
            if (tmpf > maxColVal)
            {
                maxColVal         = tmpf;
                maxColValRowIndex = i;
            }
        }

        if (maxColVal == 0)
        {
            std::fprintf(stderr, "Warning in gmx::LUPDecomposition(): input matrix A is singular (all values in column %lui are zero), decomposition failed.", static_cast<unsigned long int>(j));
            return false;
        }

        // swap rows according to the pivot determined above
        std::swap(p[j], p[maxColValRowIndex]);
        std::swap_ranges(lu(j), lu(j) + n, lu(maxColValRowIndex));

        // Gaussian elimination
        // Eliminate the entries of column j (all entries below the diagonal entry of row j).
        // The remaining non-zero entries form the upper right triangular matrix
        // In place of the zero-entries created by the Gaussian elimination step,
        // the multiplier used in subtracting C times row j from row i is recorded
        const index_type j_plus_1       = j + 1;
        const index_type l              = (j_plus_1) % simdWidth;
        const index_type m              = (simdWidth - l) * static_cast<index_type>(l != 0);
        const index_type simdBlockBegin = std::min(j_plus_1 + m, maxSimdBlockEnd);
        for (index_type i = j_plus_1; i < n; ++i)
        {
            // the multiplier, forms an entry of L, the pivoting step guarantees lu(j, j) != 0
            const value_type c = lu(i, j) / lu(j, j);
            lu(i, j) = c;
            simd_float_type  cSimd(c);

            // subtract C times row j from row i to eliminate all entries in column j
            if (simdBlockBegin > j && simdBlockBegin < maxSimdBlockEnd)
            {
                std::transform(lu(i) + j_plus_1, lu(i) + static_cast<ptrdiff_t>(simdBlockBegin),
                               lu(j) + j_plus_1,
                               lu(i) + j_plus_1,
                               [&c](value_type &elementRowI, value_type &elementRowJ) { return elementRowI - c * elementRowJ; });

                for (index_type k = simdBlockBegin; k < maxSimdBlockEnd; k += simdWidth)
                {
                    aSimd = gmx::load<simd_float_type>(lu(j) + k);
                    bSimd = gmx::load<simd_float_type>(lu(i) + k);
                    aSimd = bSimd - c * aSimd;
                    gmx::store(lu(i) + k, aSimd);
                }

                std::transform(lu(i) + maxSimdBlockEnd, lu(i) + n,
                               lu(j) + maxSimdBlockEnd,
                               lu(i) + maxSimdBlockEnd,
                               [&c](value_type &elementRowI, value_type &elementRowJ) { return elementRowI - c * elementRowJ; });
            }
            else
            {
                std::transform(lu(i) + j_plus_1, lu(i) + n,
                               lu(j) + j_plus_1,
                               lu(i) + j_plus_1,
                               [&c](value_type &elementRowI, value_type &elementRowJ) { return elementRowI - c * elementRowJ; });
            }
        }
    }

    return true;
}

/*! \brief LU decomposition with row pivoting \f$ \mathbf{A} = \mathbf{P}^{-1}\,\mathbf{L}\,\mathbf{U} \f$
           reference implementation, returns true if successful and false otherwise

    \ingroup module_math
    \libinternal

    \tparam      T           data type to be used in the calculation
    \tparam      Alloc       data type to be used in the calculation
    \tparam      TInt        integer data type to be used for indexing array elements
    \tparam      AllocInt    allocator for TInt

    \param[in]   a     \f$ n \times n\f$ matrix \f$ \mathbf{A} \f$ to be decomposed
    \param[out]  p     vector representation of the \f$ n \times n\f$ row pivot matrix \f$ \mathbf{P} \f$, entry \f$ i \f$ indicates the new row index of the current row \f$ i \f$
    \param[out]  lu    the lower and upper triangular matrices \f$ \mathbf{L} \f$ and \f$ \mathbf{U} \f$  stored in the upper right and lower left omitting the diagonal entries of \f$ \mathbf{L} \f$, which are all equal to \f$ 1 \f$
 */
template<typename T, class Alloc = std::allocator<T>, typename TInt, class AllocInt = std::allocator<TInt> >
inline bool LUPDecompositionRef(const IrregArray2D<T, Alloc, false> &a, std::vector<TInt, AllocInt> &p, IrregArray2D<T, Alloc, false> &lu)
{
    typedef typename IrregArray2D<T, Alloc, false>::value_type value_type;
    typedef typename IrregArray2D<T, Alloc, false>::index_type index_type;

    const index_type n = a.length1();

    // initialize P such that each row keeps its original index
    std::iota(p.begin(), p.end(), static_cast<index_type>(0));

    // initialize LU as copy of A
    for (index_type i = 0; i < n; ++i)
    {
        std::copy(a(i), a(i) + static_cast<ptrdiff_t>(n), lu(i));
    }

    // determine the pivot matrix by searching for the largest entry in each column
    // After row pivoting (reordering of rows), the first row has the largest entry
    // in column 1, the second row has the largest value in column 2 and so on.
    // The last element of a row will remain in each row during the decomposition.
    for (index_type j = 0; j < (n - 1); ++j)
    {
        value_type maxColVal         = 0;
        index_type maxColValRowIndex = 0;

        for (index_type i = j; i < n; ++i)
        {
            const value_type tmpf = std::abs(lu(i, j));
            if (tmpf > maxColVal)
            {
                maxColVal         = tmpf;
                maxColValRowIndex = i;
            }
        }

        if (maxColVal == 0)
        {
            std::fprintf(stderr, "Warning in gmx::LUPDecomposition(): input matrix A is singular (all values in column %lui are zero), decomposition failed.", static_cast<unsigned long int>(j));
            return false;
        }

        // swap rows according to the pivot determined above
        std::swap(p[j], p[maxColValRowIndex]);
        std::swap_ranges(lu(j), lu(j) + static_cast<ptrdiff_t>(n), lu(maxColValRowIndex));

        // Gaussian elimination
        // Eliminate the entries of column j (all entries below the diagonal entry of row j).
        // The remaining non-zero entries form the upper right triangular matrix
        // In place of the zero-entries created by the Gaussian elimination step,
        // the multiplier used in subtracting C times row j from row i is recorded
        const index_type j_plus_1 = j + 1;
        for (index_type i = j_plus_1; i < n; ++i)
        {
            // the multiplier, forms an entry of L, the pivoting step guarantees lu(j, j) != 0
            const value_type c = lu(i, j) / lu(j, j);
            lu(i, j) = c;
            // subtract the C times row j from row i to eliminate all entries in column j
            std::transform(lu(i) + j_plus_1, lu(i) + n,
                           lu(j) + j_plus_1,
                           lu(i) + j_plus_1,
                           [&c](value_type &elementRowI, value_type &elementRowJ) { return elementRowI - c * elementRowJ; });
        }
    }

    return true;
}

/*! \brief LU decomposition with row pivoting for numerical stability \f$ \mathbf{A} = \mathbf{P}^{-1}\,\mathbf{L}\,\mathbf{U} \f$,
           Returns true if successful and false otherwise.

    \f$
        \mathrm{LU} =
        \begin{bmatrix}
                         \mathrm{U}_{1,1}  &              \mathrm{U}_{1,2}  &              \mathrm{U}_{1,3}  &              \mathrm{U}_{1,4}  &              \cdots              &              \mathrm{U}_{  1,n-1}  & \mathrm{U}_{1,n} \\
            \color{cyan}{\mathrm{L}_{2,1}} &              \mathrm{U}_{2,2}  &              \mathrm{U}_{2,3}  &              \mathrm{U}_{2,4}  &              \cdots              &              \mathrm{U}_{  2,n-1}  & \mathrm{U}_{2,n} \\
            \color{cyan}{\mathrm{L}_{3,1}} & \color{cyan}{\mathrm{L}_{3,2}} &              \mathrm{U}_{3,3}  &              \mathrm{U}_{3,4}  &              \cdots              &              \mathrm{U}_{  3,n-1}  & \mathrm{U}_{3,n} \\
            \color{cyan}{\mathrm{L}_{4,1}} & \color{cyan}{\mathrm{L}_{4,2}} & \color{cyan}{\mathrm{L}_{4,3}} &              \mathrm{U}_{4,4}  &              \cdots              &              \mathrm{U}_{  4,n-1}  & \mathrm{U}_{4,n} \\
            \color{cyan}{\mathrm{L}_{5,1}} & \color{cyan}{\mathrm{L}_{5,2}} & \color{cyan}{\mathrm{L}_{5,3}} & \color{cyan}{\mathrm{L}_{5,4}} &              \ddots              &              \vdots                & \vdots           \\
            \color{cyan}{\vdots          } & \color{cyan}{\vdots          } & \color{cyan}{\vdots          } & \color{cyan}{\ddots          } & \color{cyan}{\ddots            } &              \ddots                & \vdots           \\
            \color{cyan}{\mathrm{L}_{n,1}} & \color{cyan}{\mathrm{L}_{n,2}} & \color{cyan}{\mathrm{L}_{n,3}} & \color{cyan}{\cdots          } & \color{cyan}{\mathrm{L}_{n,n-2}} & \color{cyan}{\mathrm{L}_{n  ,n-1}} & \mathrm{U}_{n,n} \\
        \end{bmatrix}
    \f$

    the upper and lower triangular matrices \f$ \mathbf{U} \f$ and \f$ \mathbf{L} \f$ are
    stored in the upper right and lower left. The diagonal entries of \f$ \mathbf{L} \f$
    are omitted. These entries are all set equal to \f$ 1 \f$, which is a common convention
    to obtain a unique decomposition.

    \ingroup module_math
    \libinternal

    \tparam      T           data type to be used in the calculation
    \tparam      Alloc       data type to be used in the calculation
    \tparam      TInt        integer data type to be used for indexing array elements
    \tparam      AllocInt    allocator for TInt

    \param[in]   a    \f$ n \times n\f$ matrix \f$ \mathbf{A} \f$ to be decomposed
    \param[out]  p    vector representation of the \f$ n \times n\f$ row pivot matrix \f$ \mathbf{P} \f$, entry \f$ i \f$ indicates the new row index of the current row \f$ i \f$
    \param[out]  lu   the lower and upper triangular matrices \f$ \mathbf{L} \f$ and \f$ \mathbf{U} \f$  stored in the upper right and lower left, omitting the diagonal entries of \f$ \mathbf{L} \f$, which are all equal to \f$ 1 \f$
 */
template<typename T, class Alloc = std::allocator<T>, typename TInt, class AllocInt = std::allocator<TInt> >
inline bool LUPDecomposition(const IrregArray2D<T, Alloc, false> &a, std::vector<TInt, AllocInt> &p, IrregArray2D<T, Alloc, false> &lu)
{
#if GMX_SIMD
    return LUPDecompositionSimd(a, p, lu);
#else
    return LUPDecompositionRef(a, p, lu);
#endif
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------

/*! \brief matrix inversion based on LUP decomposition, returns true if successful and false otherwise
           Currently, SIMD support is limited to float and/or double depending on the architecture.
           This template function is only active for unsupported data types and returns the reference
           implementation.

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   a        \f$ n \times n\f$ matrix \f$ \mathbf{A} \f$ to be decomposed
    \param[out]  aInv     the inverse \f$ \mathbf{A}^{-1} \f$ of matrix \f$ \mathbf{A} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline
typename std::enable_if<!(gmx::SimdSetup<T>::doSimd_), bool>::type
matrixInverseLUPSimd(const IrregArray2D<T, Alloc, false> &a, IrregArray2D<T, Alloc, false> &aInv)
{
    return matrixInverseLUPRef(a, aInv);
}

/*! \brief matrix inversion based on LUP decomposition, SIMD implemention, returns true if successful and false otherwise,
           SIMD implementation

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   a        \f$ n \times n\f$ matrix \f$ \mathbf{A} \f$ to be decomposed
    \param[out]  aInv     the inverse \f$ \mathbf{A}^{-1} \f$ of matrix \f$ A \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline
typename std::enable_if<gmx::SimdSetup<T>::doSimd_, bool>::type
matrixInverseLUPSimd(const IrregArray2D<T, Alloc, false> &a, IrregArray2D<T, Alloc, false> &aInv)
{
    typedef           IrregArray2D<T, Alloc, false>            array_type;
    typedef typename  array_type::value_type                   value_type;
    typedef typename  array_type::allocator_type               allocator_type;
    typedef typename  array_type::index_type                   index_type;
    typedef typename  std::vector<value_type, allocator_type>  vector_type;

    typedef typename  gmx::SimdSetup<index_type>::allocator_type     index_allocator_type;
    typedef typename  std::vector<index_type, index_allocator_type>     index_vector_type;

    // set up at compile time
    static constexpr index_type simdWidth  = SimdSetup<T>::simdWidth_;
    typedef typename   SimdSetup<T>::simd_type     simd_float_type;

    const index_type n                = a.length1();
    const index_type maxSimdBlockEnd  = (n > simdWidth ? n - n % simdWidth : 0);

    // vector representation of the permutation matrix P in LUPDecomposition()
    // original rows of matrix A are mapped to new rows after pivoting row i -> row pVec[i]
    index_vector_type p(n, static_cast<value_type>(0));
    // packed upper and lower triangular matrices, diagonal entries of L omitted (all equal to 1)
    array_type        lu(n, n, static_cast<value_type>(0));

    // temporary work space
    vector_type  x(n, static_cast<value_type>(0));
    vector_type  y(n, static_cast<value_type>(0));
    array_type   b(n, n, static_cast<value_type>(0));

    // call the reference version here too for a clean performance comparison
    const bool lupSuccess = gmx::LUPDecompositionSimd(a, p, lu);
    if (!lupSuccess)
    {
        std::fprintf(stderr, "Warning in gmx::matrixInverseLUP(): matrix inversion failed because the LUP decomposition of input matrix A failed.");
        return false;
    }

    simd_float_type  aSimd(0);
    simd_float_type  bSimd(0);
    simd_float_type  cSimd(0);

    // solve L U x = P e, where e is a column vector of the identity matrix I
    for (index_type i = 0; i < n; ++i)
    {
        // store the column vector i of I in row i of B,
        // B should be value initialized all-zero entries
        b(i, i) = 1;

        // solve L y = P B
        for (index_type j = 0; j < n; ++j)
        {
            const index_type simdBlockWidth = j - j % simdWidth;
            cSimd = gmx::setZero();
            for (index_type k = 0; k < simdBlockWidth; k += simdWidth)
            {
                aSimd = gmx::load<simd_float_type>(lu(j)    + k);
                bSimd = gmx::load<simd_float_type>(y.data() + k);
                cSimd = gmx::fma(aSimd, bSimd, cSimd);
            }
            value_type tmpf = gmx::reduce(cSimd);
            // accumulate the rest in scalar math
            tmpf += std::inner_product(lu(j)    + simdBlockWidth, lu(j) + j,
                                       y.data() + simdBlockWidth, static_cast<value_type>(0));
            y[j] = b(i, p[j]) - tmpf;
        }

        // solve U x = y
        for (index_type j = n - 1; j >= 0; --j)
        {
            const index_type j_plus_1         = j + 1;
            const index_type l                = j_plus_1 % simdWidth;
            const index_type m                = (simdWidth - l) * static_cast<index_type>(l != 0);
            const index_type simdBlockBegin   = std::min(j_plus_1 + m, maxSimdBlockEnd);
            value_type       tmpf             = 0;

            if (simdBlockBegin > j && simdBlockBegin < maxSimdBlockEnd)
            {
                tmpf = std::inner_product(lu(j)    + j_plus_1, lu(j) + simdBlockBegin,
                                          x.data() + j_plus_1, static_cast<value_type>(0));

                cSimd = gmx::setZero();
                for (index_type k = simdBlockBegin; k < maxSimdBlockEnd; k += simdWidth)
                {
                    aSimd = gmx::load<simd_float_type>(lu(j)    + k);
                    bSimd = gmx::load<simd_float_type>(x.data() + k);
                    cSimd = gmx::fma(aSimd, bSimd, cSimd);
                }
                tmpf += gmx::reduce(cSimd);

                tmpf += std::inner_product(lu(j)    + maxSimdBlockEnd, lu(j) + n,
                                           x.data() + maxSimdBlockEnd, static_cast<value_type>(0));
            }
            else
            {
                tmpf = std::inner_product(lu(j)    + j_plus_1, lu(j) + n,
                                          x.data() + j_plus_1, static_cast<value_type>(0));
            }
            x[j] = (y[j] - tmpf) / lu(j, j);
        }

        // x contains a column vector of A^{-1}, which is copied to the output matrix now
        for (index_type j = 0; j < n; ++j)
        {
            aInv(j, i) = x[j];
        }
    }

    return true;
}

/*! \brief matrix inversion based on LUP decomposition, reference implemention, returns true if successful and false otherwise,
           reference implementation

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   a        \f$ n \times n\f$ matrix \f$ \mathbf{A} \f$ to be decomposed
    \param[out]  aInv     the inverse \f$ \mathbf{A}^{-1} \f$ of matrix \f$ \mathbf{A} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline bool matrixInverseLUPRef(const IrregArray2D<T, Alloc, false> &a, IrregArray2D<T, Alloc, false> &aInv)
{
    typedef           IrregArray2D<T, Alloc, false>            array_type;
    typedef typename  array_type::value_type                   value_type;
    typedef typename  array_type::allocator_type               allocator_type;
    typedef typename  array_type::index_type                   index_type;
    typedef typename  std::vector<value_type, allocator_type>  vector_type;

    typedef typename  gmx::SimdSetup<index_type>::allocator_type     index_allocator_type;
    typedef typename  std::vector<index_type, index_allocator_type>     index_vector_type;

    const index_type n = a.length1();

    // vector representation of the permutation matrix P in LUPDecomposition()
    // original rows of matrix A are mapped to new rows after pivoting row i -> row pVec[i]
    index_vector_type p(n, static_cast<value_type>(0));
    // packed upper and lower triangular matrices, diagonal entries of L omitted (all equal to 1)
    array_type        lu(n, n, static_cast<value_type>(0));

    // temporary work space
    vector_type  x(n, static_cast<value_type>(0));
    vector_type  y(n, static_cast<value_type>(0));
    array_type   b(n, n, static_cast<value_type>(0));

    // call the reference version here too for a clean performance comparison
    const bool lupSuccess = gmx::LUPDecompositionRef(a, p, lu);
    if (!lupSuccess)
    {
        std::fprintf(stderr, "Warning in gmx::matrixInverseLUP(): matrix inversion failed because the LUP decomposition of input matrix A failed.");
        return false;
    }

    // solve L U x = P e, where e is a column vector of the identity matrix I
    for (index_type i = 0; i < n; ++i)
    {
        // store the column vector i of I in row i of B,
        // B should be value initialized all-zero entries
        b(i, i) = 1;

        // solve L y = P B
        for (index_type j = 0; j < n; ++j)
        {
            const value_type tmpf = std::inner_product(lu(j), lu(j) + j, y.data(), static_cast<value_type>(0));
            y[j] = b(i, p[j]) - tmpf;
        }

        // solve U x = y
        for (index_type j = n - 1; j >= 0; --j)
        {
            const value_type tmpf = std::inner_product(lu(j)    + j + 1, lu(j) + n,
                                                       x.data() + j + 1, static_cast<value_type>(0));
            x[j] = (y[j] - tmpf) / lu(j, j);
        }

        // x contains a column vector of A^{-1}, which is copied to the output matrix now
        for (index_type j = 0; j < n; ++j)
        {
            aInv(j, i) = x[j];
        }
    }

    return true;
}

/*! \brief matrix inversion based on LUP decomposition, returns true if successful and false otherwise

    Compute \f$ \mathbf{A}^{-1} \f$, such that \f$ \mathbf{A}\,\mathbf{A}^{-1} = \mathbf{I} \f$,
    where I is the identity matrix. This matrix inversion variant is based on LU decomposition
    with row pivoting (LUP decomposition).

    The inverse matrix has to fulfill the equation

    \f$
     \mathbf{A} \mathbf{A}^{-1} = \mathbf{I}
    \f$

    where \f$ \mathbf{I} \f$ is the identity matrix.
    The LUP decomposition of matrix \f$ \mathbf{A} \f$ allows us to write (see LUPDecomposition())

    \f$
     \mathbf{P},\mathbf{A} = \mathbf{L}\,\mathbf{U}
    \f$

    For each column vector \f$ \mathbf{x} \f$ of the identity matrix

    \f$
     \mathbf{L} \left(\mathbf{U}\,\mathbf{A}^{-1}\right) = \mathbf{P} \mathbf{I}
    \f$

    and defining \f$ \mathbf{U}\,\mathbf{A}^{-1} = \mathbf{Y}\f$, one obtains
    a system of two equations

    \f$
     \left(1\right)\quad \mathbf{L}\,\mathbf{Y}      = \mathbf{P}\,\mathbf{I}
    \f$

    \f$
     \left(2\right)\quad \mathbf{U}\,\mathbf{A}^{-1} = \mathbf{X}
    \f$

    This equation system can be solved for each column vector \f$ \mathbf{x},\mathbf{y}, \mathbf{e} \f$
    of the matrices \f$ \mathbf{X},\mathbf{Y},\mathbf{I} \f$, respectively.
    The triangular form of the matrices \f$ \mathbf{L}\,\mathbf{U} \f$ allows
    for an easy solution via forward-substitution.

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \ingroup module_math
    \libinternal

    \param[in]   a        \f$ n \times n\f$ matrix \f$ \mathbf{A} \f$ to be decomposed
    \param[out]  aInv     the inverse \f$ \mathbf{A}^{-1} \f$ of matrix \f$ \mathbf{A} \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline bool matrixInverseLUP(const IrregArray2D<T, Alloc, false> &a, IrregArray2D<T, Alloc, false> &aInv)
{
#if GMX_SIMD
    return matrixInverseLUPSimd(a, aInv);
#else
    return matrixInverseLUPRef(a, aInv);
#endif
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------

/*! \brief compute the determinant of a matrix based on LUP decomposition, SIMD implemention,
           returns true if successful and false otherwise

           Currently, SIMD support is limited to float and/or double depending on the architecture.
           This template function is only active for unsupported data types and returns the reference
           implementation.

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   a        \f$ n \times n\f$ matrix \f$ \mathbf{A} \f$ to be decomposed
    \param[out]  det      the determinant of the input matrix \f$ \mathrm{det} \left|\mathbf{A}\right| \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline
typename std::enable_if<!(gmx::SimdSetup<T>::doSimd_), bool>::type
matrixDeterminantLUPSimd(const IrregArray2D<T, Alloc, false> &a, T *det)
{
    return matrixDeterminantLUPRef(a, det);
}

/*! \brief compute the determinant of a matrix based on LUP decomposition,
           SIMD implemention, returns true if successful and false otherwise

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   a        \f$ n \times n\f$ matrix \f$ \mathbf{A} \f$ to be decomposed
    \param[out]  det      the determinant of the input matrix \f$ \mathrm{det} \left|\mathbf{A}\right| \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline
typename std::enable_if<gmx::SimdSetup<T>::doSimd_, bool>::type
matrixDeterminantLUPSimd(const IrregArray2D<T, Alloc, false> &a, T *det)
{
    typedef           IrregArray2D<T, Alloc, false>            array_type;
    typedef typename  array_type::value_type                   value_type;
    typedef typename  array_type::allocator_type               allocator_type;
    typedef typename  array_type::index_type                   index_type;
    typedef typename  std::vector<value_type, allocator_type>  vector_type;

    typedef typename  gmx::SimdSetup<index_type>::allocator_type     index_allocator_type;
    typedef typename  std::vector<index_type, index_allocator_type>     index_vector_type;

    // set up at compile time
    static constexpr index_type simdWidth  = SimdSetup<T>::simdWidth_;
    typedef typename   SimdSetup<T>::simd_type     simd_float_type;

    const index_type n = a.length1();

    // vector representation of the permutation matrix P in LUPDecomposition()
    // original rows of matrix A are mapped to new rows after pivoting row i -> row pVec[i]
    index_vector_type p(n, static_cast<value_type>(0));
    // packed upper and lower triangular matrices, diagonal entries of L omitted (all equal to 1)
    array_type        lu(n, n, static_cast<value_type>(0));

    // temporary work space for the diagonal entries of U
    vector_type udiag(n, static_cast<value_type>(0));

    // call the reference version here too for a clean performance comparison
    const bool lupSuccess = gmx::LUPDecompositionSimd(a, p, lu);
    if (!lupSuccess)
    {
        std::fprintf(stderr, "Warning in gmx::matrixDeterminantLUP(): determinant computation failed because the LUP decomposition of input matrix A failed.");
        return false;
    }

    // first write all diagonal entries to the vector and then accumulate the product in SIMD
    for (index_type i = 0; i < n; ++i)
    {
        udiag[i] = lu(i, i);
    }
    *det = static_cast<value_type>(1);
    // accumulate the product of the elements outside simdBlockWidth with scalar math
    const index_type simdBlockWidth  = n - n % simdWidth;
    for (index_type j = simdBlockWidth; j < n; ++j)
    {
        *det *= udiag[j];
    }
    if (simdBlockWidth != 0)
    {
        // accumulate the product in SIMD vector operations
        simd_float_type  aSimd(1);
        simd_float_type  bSimd(0);
        for (index_type j = 0; j < simdBlockWidth; j += simdWidth)
        {
            bSimd = gmx::load<simd_float_type>(udiag.data() + j);
            aSimd = aSimd * bSimd;
        }
        // extract the vector elements and accumulate their product with scalar math
        gmx::store(udiag.data(), aSimd);
        for (index_type j = 0; j < simdWidth; ++j)
        {
            *det *= udiag[j];
        }
    }

    // the signum of the permutation, anything to be gained with SIMD here?, probably only for large n
    index_type s = 0;
    for (index_type i = 0; i < n; ++i)
    {
        for (index_type j = 0; j < i; ++j)
        {
            s += static_cast<index_type>(p[i] < p[j]);
        }
    }

    *det *= ((s % 2) == 0 ? static_cast<value_type>(1) : static_cast<value_type>(-1));

    return true;
}

/*! \brief compute the determinant of a matrix based on LUP decomposition,
           reference implemention, returns true if successful and false otherwise

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   a        \f$ n \times n\f$ matrix \f$ \mathbf{A} \f$ to be decomposed
    \param[out]  det      the determinant of the input matrix \f$ \mathrm{det} \left|\mathbf{A}\right| \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline bool
matrixDeterminantLUPRef(const IrregArray2D<T, Alloc, false> &a, T *det)
{
    typedef           IrregArray2D<T, Alloc, false>            array_type;
    typedef typename  array_type::value_type                   value_type;
    typedef typename  array_type::allocator_type               allocator_type;
    typedef typename  array_type::index_type                   index_type;
    typedef typename  std::vector<value_type, allocator_type>  vector_type;

    typedef typename  gmx::SimdSetup<index_type>::allocator_type     index_allocator_type;
    typedef typename  std::vector<index_type, index_allocator_type>     index_vector_type;

    const index_type n = a.length1();

    // vector representation of the permutation matrix P in LUPDecomposition()
    // original rows of matrix A are mapped to new rows after pivoting row i -> row pVec[i]
    index_vector_type p(n, static_cast<value_type>(0));
    // packed upper and lower triangular matrices, diagonal entries of L omitted (all equal to 1)
    array_type        lu(n, n, static_cast<value_type>(0));

    // temporary work space for the diagonal entries of U
    vector_type udiag(n, static_cast<value_type>(0));

    // call the reference version here too for a clean performance comparison
    const bool lupSuccess = gmx::LUPDecompositionRef(a, p, lu);
    if (!lupSuccess)
    {
        std::fprintf(stderr, "Warning in gmx::matrixDeterminantLUP(): determinant computation failed because the LUP decomposition of input matrix A failed.");
        return false;
    }

    // first write all diagonal entries to the vector and then accumulate the product
    for (index_type i = 0; i < n; ++i)
    {
        udiag[i] = lu(i, i);
    }
    *det = std::accumulate(udiag.begin(), udiag.end(), static_cast<value_type>(1), std::multiplies<value_type>());

    // the signum of the permutation
    index_type s = 0;
    for (index_type i = 0; i < n; ++i)
    {
        for (index_type j = 0; j < i; ++j)
        {
            s += static_cast<index_type>(p[i] < p[j]);
        }
    }

    *det *= ((s % 2) == 0 ? static_cast<value_type>(1) : static_cast<value_type>(-1));

    return true;
}

/*! \brief compute the determinant of a matrix based on LUP decomposition,
           returns true if successful and false otherwise

    \f$
        \mathrm{det}\left|\mathbf{A}\right| = \mathrm{det}\left|\mathbf{P}^{-1}\right|\,\mathrm{det}\left|\mathbf{L}\right|\,\mathrm{det}\left|\mathbf{U}\right|
    \f$

    where \f$ \mathrm{det}\left|\mathbf{P}^{-1}\right| = \left(-1\right)^{s} \f$ is the signum of
    the permutation. The signum is computed using the vector representation \f$ \mathbf{p} \f$ of
    the permutation matrix \f$ \mathbf{P} \f$, whose entries \f$ p_{i}, p_{j} \f$ specify the new
    index of rows \f$ i, j \f$ after permutation/pivoting. \f$ s \f$ can be set equal to the number
    of exchanges between rows \f$ j > i \f$ with \f$ p_{j} - p_{i} < 0 \f$. Formally, the determinant
    of the permutation matrix \f$ \mathrm{det}\left|\mathbf{P}^{-1}\right| \f$ can be written as

    \f$
       \left(-1\right)^{s} = \prod\limits_{i}^{n}\prod\limits_{j > i}^{n} \frac{ p_{j} - p_{i} }{ j - i }
    \f$

    The entries \f$ p_{i}, p_{j} \f$ of the vector representation \f$ \mathbf{p} \f$ of the permutation
    matrix specify the new index of rows \f$ i, j \f$ after permutation/pivoting.

    \f$ \mathrm{det}\left|\mathbf{L}\right| = 1\f$ because \f$ \mathbf{L} \f$ is triangular and all its
    diagonal entries are equal to \f$ 1 \f$ by convention (see LUPDecomposition()). Since \f$ \mathbf{U} \f$
    is also a triangular matrix, its determinant is equal to its trace, i.e., to the product of its diagonal
    entries. Thus,

    \f$
        \mathrm{det}\left|\mathbf{A}\right| = \left(-1\right)^{s}\,\prod\limits_{i}^{n} U_{i, i}
    \f$

    \ingroup module_math
    \libinternal

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

    \param[in]   a        \f$ n \times n\f$ matrix \f$ \mathbf{A} \f$ to be decomposed
    \param[out]  det      the determinant of the input matrix \f$ \mathrm{det} \left|\mathbf{A}\right| \f$
 */
template<typename T, class Alloc = std::allocator<T> >
inline bool
matrixDeterminantLUP(const IrregArray2D<T, Alloc, false> &a, T *det)
{
#if GMX_SIMD
    return matrixDeterminantLUPSimd(a, det);
#else
    return matrixDeterminantLUPRef(a, det);
#endif
}


} // end namespace gmx

#endif
