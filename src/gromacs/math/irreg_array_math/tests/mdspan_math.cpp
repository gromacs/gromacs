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
/*! \internal \file
   \brief
   Tests for math operations with IrregArray classes

   For development, the tests can be run with a '-stdout' command-line option
   to print out the help to stdout instead of using the XML reference
   framework.

   \author R. Thomas Ullmann <tullman@gwdg.de>

   \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/simd/simd_setup.h"
#include "gromacs/mdspan/mdspan.h"

#include "irreg_array_math_test_commons.h"

#define DO_TIMING
#ifdef DO_TIMING
#include <chrono>
#endif

namespace gmx
{

namespace irreg_array_math_test
{

// add the parameter interface for automated testing for multiple data types,
// the template parameter for ::testing::WithParamInterface<>, is a data type
// for value parametrized testing
template<typename TReal = double>
class ParamMdspanMathTest : public IrregArrayMathTest,
                                public ::testing::WithParamInterface<TReal>
{
    public:
        //! Unfortunately, for typed test fixtures one still has to explicitly
        //! state this->member within the test case instead of just member to
        //! access members of the base class. The using statement has no effect
        //! there. Base class members can directly be used if the derived fixture
        //! class is not a class template, which is not possible for a typed test.
        using IrregArrayMathTest::debug;
        using IrregArrayMathTest::error;
};

namespace
{

using std::endl;

using namespace std::experimental;

/*! \brief  output stream operator to insert a string representation of
            a vector into an output stream, e.g., cout

    \tparam   TReal        floating point data type

    \param[in]   output   ostream in which the content is to be inserted
    \param[in]   vec      the object whose contents are to be inserted
 */
template<typename TReal = double,
         class Extents = extents<dynamic_extent, dynamic_extent>,
         class LayoutPolicy,
         class AccessorPolicy>
std::ostream &operator<<(std::ostream &output,
                         const basic_mdspan<TReal, Extents, LayoutPolicy, AccessorPolicy> &mat)
{
    typedef basic_mdspan<TReal, Extents, LayoutPolicy, AccessorPolicy> matrix_type;
    typedef typename matrix_type::index_type index_type;

    if (mat.extent(0) == 0 || mat.extent(1) == 0)
    {
        return output;
    }
    constexpr size_t extraFloatWidth = 4;
    for (index_type i = 0; i < mat.extent(0); ++i)
    {
        for (index_type j = 0; j < mat.extent(1); ++j)
        {
            output << std::setw(extraFloatWidth + std::numeric_limits<TReal>::digits)
                   << std::fixed << std::setprecision(std::numeric_limits<TReal>::digits) << mat(i, j);
            if (j < (mat.extent(1) - 1))
            {
                output << "  ";
            }
            else
            {
                output << "\n";
            }
        }
    }
    return output;
}

/*! \brief  output stream operator to insert a string representation of
            a vector into an output stream, e.g., cout

    \tparam   TReal        floating point data type

    \param[in]   output   ostream in which the content is to be inserted
    \param[in]   vec      the object whose contents are to be inserted
 */
template<typename TReal = double, class Alloc = std::allocator<TReal> >
std::ostream &operator<<(std::ostream &output, const std::vector<TReal, Alloc> &vec)
{
    if (vec.empty())
    {
        return output;
    }
    constexpr size_t extraFloatWidth = 4;
    for (typename std::vector<TReal, Alloc>::size_type i = 0; i < vec.size(); ++i)
    {
        output << std::setw(extraFloatWidth + std::numeric_limits<TReal>::digits) << std::fixed << std::setprecision(std::numeric_limits<TReal>::digits) << vec[i];
        if (i < (vec.size() - 1))
        {
            output << "  ";
        }
    }
    return output;
}

/*********************************************************************/

/*! \brief helper function for computing the matrix product

    \tparam   TMat    matrix type

    \param[in]    a   input  matrix A, M x L elements
    \param[in]    b   input  matrix B, L x N elements
    \param[out]   c   output matrix C, M x N elements, c = a b
 */
template<class TMat>
void computeMatrixProduct(const TMat a, const TMat &b, TMat &c)
{
    typedef typename TMat::value_type float_type;
    typedef typename TMat::index_type index_type;

    // number of columns    of matrix A
    const index_type m = a.extent(0);
    // number of columns    of matrix A
    const index_type l = a.extent(1);
    // number of rows       of matrix B
    const index_type n = b.extent(1);

    // initialize matrix c filling the whole buffer
    std::fill(c.data(), c.data() + static_cast<ptrdiff_t>(m * n), static_cast<float_type>(0));

    for (index_type i = 0; i < m; ++i)
    {
        for (index_type k = 0; k < l; ++k)
        {
            // c[i][j] = v[i][:] a[:][j]
            // scalar product of row vector i of m x l matrix A and column vector j of l x n matrix B
            for (index_type j = 0; j < n; ++j)
            {
                c(i, j) += a(i, k) * b(k, j);
            }
        }
    }
}

/*! \brief helper function for computing the matrix product

    \tparam   TMat    matrix type

    \param[in]    a   input  matrix A^T, L x M elements
    \param[in]    b   input  matrix B  , L x N elements
    \param[out]   c   output matrix C  , M x N elements, c = a b
 */
template<class TMat>
void computeMatrixProductMatrixTMatrix(const TMat a, const TMat &b, TMat &c)
{
    typedef typename TMat::value_type float_type;
    typedef typename TMat::index_type index_type;

    // number of columns    of matrix A
    const index_type m = a.extent(1);
    // number of columns    of matrix A
    const index_type l = a.extent(0);
    // number of rows       of matrix B
    const index_type n = b.extent(1);

    // initialize matrix c filling the whole buffer
    std::fill(c.data(), c.data() + static_cast<ptrdiff_t>(m * n), static_cast<float_type>(0));

    for (index_type k = 0; k < l; ++k)
    {
        for (index_type i = 0; i < m; ++i)
        {
            // c[i][j] = v[i][:] a[:][j]
            // scalar product of row vector i of m x l matrix A and column vector j of l x n matrix B
            for (index_type j = 0; j < n; ++j)
            {
                c(i, j) += a(k, i) * b(k, j);
            }
        }
    }
}

/*! \brief helper function for computing the matrix product

    \tparam   TMat    matrix type

    \param[in]    a   input  matrix A^T, L x M elements
    \param[in]    b   input  matrix B^T, N x L elements
    \param[out]   c   output matrix C  , M x N elements, c = a b
 */
template<class TMat>
void computeMatrixProductMatrixTMatrixT(const TMat a, const TMat &b, TMat &c)
{
    typedef typename TMat::value_type float_type;
    typedef typename TMat::index_type index_type;

    // number of columns    of matrix A
    const index_type m = a.extent(1);
    // number of columns    of matrix A
    const index_type l = a.extent(0);
    // number of rows       of matrix B
    const index_type n = b.extent(0);

    // initialize matrix c filling the whole buffer
    std::fill(c.data(), c.data() + static_cast<ptrdiff_t>(m * n), static_cast<float_type>(0));

    for (index_type j = 0; j < n; ++j)
    {
        for (index_type i = 0; i < m; ++i)
        {
            for (index_type k = 0; k < l; ++k)
            {
                // c[i][j] = v[i][:] a[:][j]
                // scalar product of row vector i of m x l matrix A and column vector j of l x n matrix B
                c(i, j) += a(k, i) * b(j, k);
            }
        }
    }
}

/*! \brief helper function for computing the matrix product

    \tparam   TMat    matrix type

    \param[in]    a   input  matrix A,   M x L elements
    \param[in]    b   input  matrix B^T, N x L elements
    \param[out]   c   output matrix,     M x N elements, c = a b
 */
template<class TMat>
void computeMatrixProductMatrixMatrixT(const TMat a, const TMat &b, TMat &c)
{
    typedef typename TMat::value_type float_type;
    typedef typename TMat::index_type index_type;

    // number of columns    of matrix A
    const index_type m = a.extent(0);
    // number of columns    of matrix A
    const index_type l = a.extent(1);
    // number of rows       of matrix B
    const index_type n = b.extent(0);

    for (index_type i = 0; i < m; ++i)
    {
        for (index_type j = 0; j < n; ++j)
        {
            // c[i][j] = v[i][:] a[:][j]
            // scalar product of row vector i of m x l matrix A and column vector j of l x n matrix B
            float_type tmpf = 0;
            for (index_type k = 0; k < l; ++k)
            {
                tmpf += a(i, k) * b(j, k);
            }
            c(i, j) = tmpf;
        }
    }
}

/*! \brief helper function for computing the vector matrix product w = v A^T

    \tparam   TVec    vector type for storing TVal
    \tparam   TMat    matrix type for storing TVal

    \param[in]    v   input  vector, 1 x M elements
    \param[in]    aT  input  matrix, N x M elements
    \param[out]   c   output vector, 1 x N elements, c = v * a^T
 */
template<class TVec, class TMat = gmx::IrregArray2D<typename TVec::value_type> >
void computeMatrixProductVecMatTranspose(const TVec &v, const TMat &aT, TVec &c)
{
    typedef typename TVec::value_type float_type;

    // number of columns    of the matrix
    const typename TMat::index_type M = aT.extent(1);
    // number of rows       of the matrix
    const typename TMat::index_type N = aT.extent(0);

    std::fill(c.begin(), c.end(), static_cast<float_type>(0));

    for (typename TMat::index_type cj = 0; cj < N; ++cj)
    {
        // c[1][cj] = v[1][:] b[cj][:], since b = a^T
        // c[1][cj] = v[1][:] a[:][cj]
        // scalar product of vector v (row 1 of the 1xM matrix) and column vector cj of matrix a
        for (typename TMat::index_type vi = 0; vi < M; ++vi)
        {
            c[cj] += v[vi] * aT(cj, vi);
        }
    }
}

/*! \brief helper function for computing the vector matrix product w = v A

    \tparam   TVec    vector type for storing TVal
    \tparam   TMat    matrix type for storing TVal

    \param[in]    v   input  vector, 1 x M elements
    \param[in]    a   input  matrix, M x N elements
    \param[out]   c   output vector, 1 x N elements, c = v * a
 */
template<class TVec, class TMat = gmx::IrregArray2D<typename TVec::value_type> >
void computeMatrixProductVecMat(const TVec &v, const TMat &a, TVec &c)
{
    typedef typename TVec::value_type float_type;

    // number of rows       of the matrix
    const typename TMat::index_type M = a.extent(0);
    // number of columns    of the matrix
    const typename TMat::index_type N = a.extent(1);

    std::fill(c.begin(), c.end(), static_cast<float_type>(0));

    for (typename TMat::index_type cj = 0; cj < N; ++cj)
    {
        // c[1][cj] = v[1][:] b[cj][:], since b = a^T
        // c[1][cj] = v[1][:] a[:][cj]
        // scalar product of vector v (row 1 of the 1xM matrix) and column vector cj of matrix a
        for (typename TMat::index_type vi = 0; vi < M; ++vi)
        {
            c[cj] += v[vi] * a(vi, cj);
        }
    }
}

/*! \brief helper function for computing the matrix vector product w = A v

    \tparam   TVec    vector type for storing TVal
    \tparam   TMat    matrix type for storing TVal

    \param[in]    a   input  matrix, M x N elements
    \param[in]    v   input  vector, N x 1 elements
    \param[out]   c   output vector, N x 1 elements, c = a * v
 */
template<class TVec, class TMat = gmx::IrregArray2D<typename TVec::value_type> >
void computeMatrixProductMatVec(const TMat &a, const TVec &v, TVec &c)
{
    typedef typename TVec::value_type float_type;

    // number of rows       of the matrix
    const typename TMat::index_type M = a.extent(0);
    // number of columns    of the matrix
    const typename TMat::index_type N = a.extent(1);

    std::fill(c.begin(), c.end(), static_cast<float_type>(0));

    for (typename TMat::index_type cj = 0; cj < M; ++cj)
    {
        // c[1][cj] = v[1][:] b[cj][:], since b = a^T
        // c[1][cj] = v[1][:] a[:][cj]
        // scalar product of vector v (row 1 of the 1xM matrix) and column vector cj of matrix a
        for (typename TMat::index_type vi = 0; vi < N; ++vi)
        {
            c[cj] += a(cj, vi) * v[vi];
        }
    }
}

/*! \brief helper function for computing the matrix vector product w = A^T v

    \tparam   TVec    vector type for storing TVal
    \tparam   TMat    matrix type for storing TVal

    \param[in]    aT  input  matrix, N x M elements
    \param[in]    v   input  vector, N x 1 elements
    \param[out]   c   output vector, N x 1 elements, c = a * v
 */
template<class TVec, class TMat = gmx::IrregArray2D<typename TVec::value_type> >
void computeMatrixProductMatTVec(const TMat &aT, const TVec &v, TVec &c)
{
    typedef typename TVec::value_type float_type;

    // number of rows       of the matrix A
    const typename TMat::index_type N = aT.extent(0);
    // number of columns    of the matrix A
    const typename TMat::index_type M = aT.extent(1);

    std::fill(c.begin(), c.end(), static_cast<float_type>(0));

    for (typename TMat::index_type cj = 0; cj < M; ++cj)
    {
        for (typename TMat::index_type vi = 0; vi < N; ++vi)
        {
            c[cj] += aT(vi, cj) * v[vi];
        }
    }
}

/*! \brief helper function for element-wise comparison of two vectors, returns the maximum absolute deviation

    \tparam   TVec    vector type

    \param[in]    a   input  vector, N elements
    \param[in]    b   input  vector, N elements
 */
template<class TVec>
typename TVec::value_type maxAbsDeviationVectorElements(const TVec &a, const TVec &b)
{
    typedef typename TVec::value_type float_type;

    float_type maxDev = 0;
    for (size_t i = 0; i < a.size(); ++i)
    {
        const float_type tmpd = (a[i] > b[i] ? (a[i] - b[i]) : (b[i] - a[i]));
        if (tmpd > maxDev)
        {
            maxDev = tmpd;
        }
    }

    return maxDev;
}

/*! \brief helper function for element-wise comparison of two matrices, returns the maximum absolute deviation

    \tparam   TMat    matrix type

    \param[in]    a   input  vector, N elements
    \param[in]    b   input  vector, N elements
 */
template<class TMat>
typename TMat::value_type maxAbsDeviationMatrixElements(const TMat &a, const TMat &b)
{
    typedef typename TMat::value_type float_type;

    float_type maxDev = 0;
    for (size_t i = 0; i < a.extent(0); ++i)
    {
        for (size_t j = 0; j < a.extent(1); ++j)
        {
            const float_type tmpd = (a(i, j) > b(i, j) ? (a(i, j) - b(i, j)) : (b(i, j) - a(i, j)));
            if (tmpd > maxDev)
            {
                maxDev = tmpd;
            }
        }
    }

    return maxDev;
}


/*********************************************************************/

// associate a list of types with the test case, which will be repeated for each type in the list
// list all tests names as macro parameters following the test case name
TYPED_TEST_CASE(ParamMdspanMathTest, TestFloatTypes);


TYPED_TEST(ParamMdspanMathTest, VectorMatrixProduct)
{
    typedef TypeParam float_type;

    // tolerance needs to be that high for larger matrices
    const float_type tolFac  = 100;
    const float_type tolMax  = 0.1;
    float_type       tolDiff = static_cast<float_type>(tolFac) / static_cast<float_type>(pow(std::numeric_limits<float_type>::digits, 2));
    if (tolDiff > tolMax)
    {
        tolDiff = tolMax;
    }
    typedef typename gmx::SimdSetup<float_type>::allocator_type allocator_type;

    typedef std::vector<float_type, allocator_type>                 vector_type;
    typedef basic_mdspan<float_type, extents<dynamic_extent, dynamic_extent>,
                         layout_left, accessor_basic<float_type> >   array_type;

    // testing the vector matrix product cVec^T  = vVec^T  aMat
    //                                   (1 x m) = (1 x n) (n x m)

    // n x m matrix sizes
    const std::vector<size_t> n = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 60, 67, 1200, 8000};
    const std::vector<size_t> m = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 60, 73, 1200, 8000};

    const size_t              nruns = 10;

    // testing the vector matrix product cVec = vVec aMat
    for (size_t i = 0; i < m.size(); ++i)
    {
        vector_type vVec(m[i], static_cast<float_type>(0));
        for (size_t j = 0; j < n.size(); ++j)
        {

            this->debug << "\nTesting vector-matrix product " << " (float type " << typeid(float_type).name() << ") for a " << m[i] << " x " << n[j] << " matrix ..." << endl;

            vector_type  cVec(    n[j], static_cast<float_type>(0));
            vector_type  cVecRef( n[j], static_cast<float_type>(0));
            vector_type  cVecSimd(n[j], static_cast<float_type>(0));

            vector_type aBuf(m[i] * n[j], static_cast<float_type>(0));
            vector_type aBufT(m[i] * n[j], static_cast<float_type>(0));
            array_type   aMat(aBuf.data(), m[i],  n[j]);
            array_type   aMatT(aBufT.data(), n[j], m[i]);

            // make sure that the result fits by filling the input vector and array with non-uniform values
            // vVec[i] = i / C
            // C = series y_i = sum_i^n x_i with x_i = x_{i-1} + 1,  n * (n + 1) / 2
            // aMat(i, j) = vVec[i] * j / C'
            // C = series y_i = sum_i^n x_i with x_i = x_{i-1} + 1,  m * (m + 1) / 2
            const float_type normfac = 10.0 * std::sqrt(static_cast<float_type>(m[i] * n[j])) / static_cast<float_type>(n[j] * (n[j] + 1));
            for (size_t k = 0; k < m[i]; ++k)
            {
                vVec[k] = static_cast<float_type>(k + 1) * normfac;
            }

            for (size_t k = 0; k < m[i]; ++k)
            {
                const float_type normfacRow = 2.0 / static_cast<float_type>(m[i] * (m[i] + 1));
                for (size_t l = 0; l < n[j]; ++l)
                {
                    aMat(k, l)  = vVec[k] * static_cast<float_type>(l + 1) * normfacRow;
                    aMatT(l, k) = aMat(k, l);
                }
            }


            this->debug << "Test input data:"        << endl;
            this->debug << "v        ="   << vVec    << endl;
            this->debug << "A        =\n" << aMat    << endl;
            this->debug << "A^T      =\n" << aMatT   << endl;

            // ==============================================
            //         using the matrix directly
            // ==============================================

#ifdef DO_TIMING
            auto start1 = std::chrono::high_resolution_clock::now();
#endif
            for (size_t r = 0; r < nruns; ++r)
            {
                // the reference solution
                computeMatrixProductVecMat(vVec, aMat, cVecRef);
            }

#ifdef DO_TIMING
            auto stop1 = std::chrono::high_resolution_clock::now();
            this->debug << "Runtime of the reference vector-matrix product using A  : "
            << std::chrono::duration<double, std::milli>(stop1 - start1).count() << " ms"
            << endl;
#endif

/*
#ifdef DO_TIMING
            auto start2 = std::chrono::high_resolution_clock::now();
#endif
            for (size_t r = 0; r < nruns; ++r)
            {
                // the data structure-provided reference solution
                gmx::matrixProductVectorMatrixRef(vVec, aMat, cVec);
            }
#ifdef DO_TIMING
            auto stop2 = std::chrono::high_resolution_clock::now();
            this->debug << "Runtime of the scalar    vector-matrix product using A  : "
            << std::chrono::duration<double, std::milli>(stop2 - start2).count() << " ms"
            << endl;
#endif

#ifdef DO_TIMING
            auto start3 = std::chrono::high_resolution_clock::now();
#endif
            for (size_t r = 0; r < nruns; ++r)
            {
                // the data structure-provided reference solution
                gmx::matrixProductVectorMatrixSimd(vVec, aMat, cVecSimd);
            }
#ifdef DO_TIMING
            auto stop3 = std::chrono::high_resolution_clock::now();
            this->debug << "Runtime of the SIMD      vector-matrix product using A  : "
            << std::chrono::duration<double, std::milli>(stop3 - start3).count() << " ms"
            << endl;
#endif
*/
            this->debug << "vector-matrix product using A:" << endl;
            this->debug << "c (scalar) ="   << cVec         << endl;
            this->debug << "c (SIMD)   ="   << cVecSimd     << endl;
            this->debug << "c (ref.)   ="   << cVecRef      << endl;

/*
            EXPECT_NEAR(maxAbsDeviationVectorElements(cVec,     cVecRef), 0, tolDiff)
            << "Scalar vector matrix product directly using the matrix as input differs from the reference solution." << std::endl;
            EXPECT_NEAR(maxAbsDeviationVectorElements(cVecSimd, cVecRef), 0, tolDiff)
            << "SIMD   vector matrix product directly using the matrix as input differs from the reference solution." << std::endl;
*/
            std::fill(    cVec.begin(),     cVec.end(), static_cast<float_type>(0));
            std::fill( cVecRef.begin(),  cVecRef.end(), static_cast<float_type>(0));
            std::fill(cVecSimd.begin(), cVecSimd.end(), static_cast<float_type>(0));

            // ==============================================
            //         using the transpose matrix
            // ==============================================

#ifdef DO_TIMING
            auto start4 = std::chrono::high_resolution_clock::now();
#endif
            for (size_t r = 0; r < nruns; ++r)
            {
                // the reference solution
                computeMatrixProductVecMatTranspose(vVec, aMatT, cVecRef);
            }
#ifdef DO_TIMING
            auto stop4 = std::chrono::high_resolution_clock::now();
            this->debug << "Runtime of the reference vector-matrix product using A^T: "
            << std::chrono::duration<double, std::milli>(stop4 - start4).count() << " ms"
            << endl;
#endif

/*
#ifdef DO_TIMING
            auto start5 = std::chrono::high_resolution_clock::now();
#endif
            for (size_t r = 0; r < nruns; ++r)
            {
                // the data structure-provided reference solution
                gmx::matrixProductVectorMatrixTRef(vVec, aMatT, cVec);
            }
#ifdef DO_TIMING
            auto stop5 = std::chrono::high_resolution_clock::now();
            this->debug << "Runtime of the scalar    vector-matrix product using A^T: "
            << std::chrono::duration<double, std::milli>(stop5 - start5).count() << " ms"
            << endl;
#endif

#ifdef DO_TIMING
            auto start6 = std::chrono::high_resolution_clock::now();
#endif
            for (size_t r = 0; r < nruns; ++r)
            {
                // the reference solution
                gmx::matrixProductVectorMatrixTSimd(vVec, aMatT, cVecSimd);
            }
#ifdef DO_TIMING
            auto stop6 = std::chrono::high_resolution_clock::now();
            this->debug << "Runtime of the SIMD      vector-matrix product using A^T: "
            << std::chrono::duration<double, std::milli>(stop6 - start6).count() << " ms"
            << endl;
#endif
*/

            this->debug << "vector-matrix product using A^T:" << endl;
            this->debug << "c (scalar) ="   << cVec           << endl;
            this->debug << "c (SIMD)   ="   << cVecSimd       << endl;
            this->debug << "c (ref.)   ="   << cVecRef        << endl;

/*
            EXPECT_NEAR(maxAbsDeviationVectorElements(cVec,     cVecRef), 0, tolDiff)
            << "Scalar vector matrix product using the transpose matrix as input deviates from the reference solution." << std::endl;
            EXPECT_NEAR(maxAbsDeviationVectorElements(cVecSimd, cVecRef), 0, tolDiff)
            << "SIMD   vector matrix product using the transpose matrix as input deviates from the reference solution." << std::endl;
*/
        }
    }

}

TYPED_TEST(ParamMdspanMathTest, MatrixVectorProduct)
{
    typedef TypeParam float_type;

    // tolerance needs to be that high for larger matrices
    const float_type tolFac  = 100;
    const float_type tolMax  = 0.1;
    float_type       tolDiff = static_cast<float_type>(tolFac) / static_cast<float_type>(pow(std::numeric_limits<float_type>::digits, 2));
    if (tolDiff > tolMax)
    {
        tolDiff = tolMax;
    }

    typedef typename gmx::SimdSetup<float_type>::allocator_type allocator_type;

    typedef std::vector<float_type, allocator_type>                 vector_type;
    typedef basic_mdspan<float_type, extents<dynamic_extent, dynamic_extent>,
                         layout_left, accessor_basic<float_type> >   array_type;

    // testing the vector matrix product cVec^T  =  aMat    vVec^T
    //                                   (m x 1) = (m x n) (n x 1)

    // m x n matrix sizes
    const std::vector<size_t> m = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 60, 67, 1200, 8000};
    const std::vector<size_t> n = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 60, 73, 1200, 8000};

    const size_t              nruns = 10;

    // testing the vector matrix product cVec = vVec aMat
    for (size_t i = 0; i < m.size(); ++i)
    {
        for (size_t j = 0; j < n.size(); ++j)
        {
            this->debug << "\nTesting matrix-vector product " << " (float type " << typeid(float_type).name() << ") for a " << m[i] << " x " << n[j] << " matrix ..." << endl;

            vector_type  vVec(n[j], static_cast<float_type>(0));

            vector_type  cVec(    m[i], static_cast<float_type>(0));
            vector_type  cVecRef( m[i], static_cast<float_type>(0));
            vector_type  cVecSimd(m[i], static_cast<float_type>(0));

            vector_type aBuf(m[i] * n[j], static_cast<float_type>(0));
            vector_type aBufT(m[i] * n[j], static_cast<float_type>(0));
            array_type   aMat(aBuf.data(), m[i],  n[j]);
            array_type   aMatT(aBufT.data(), n[j], m[i]);

            // make sure that the result fits by filling the input vector and array with non-uniform values
            // vVec[i] = i / C
            // C = series y_i = sum_i^n x_i with x_i = x_{i-1} + 1,  n * (n + 1) / 2
            // aMat(i, j) = vVec[i] * j / C'
            // C = series y_i = sum_i^n x_i with x_i = x_{i-1} + 1,  m * (m + 1) / 2
            const float_type normfac = 10.0 * std::sqrt(static_cast<float_type>(m[i] * n[j])) / static_cast<float_type>(n[j] * (n[j] + 1));
            for (size_t k = 0; k < n[j]; ++k)
            {
                vVec[k] = static_cast<float_type>(k + 1) * normfac;
            }

            for (size_t k = 0; k < n[j]; ++k)
            {
                const float_type normfacRow = 2.0 / static_cast<float_type>(m[i] * (m[i] + 1));
                for (size_t l = 0; l < m[i]; ++l)
                {
                    aMat(l, k)  = vVec[k] * static_cast<float_type>(l + 1) * normfacRow;
                    aMatT(k, l) = aMat(l, k);
                }
            }


            this->debug << "Test input data:"        << endl;
            this->debug << "v        ="   << vVec    << endl;
            this->debug << "A        =\n" << aMat    << endl;
            this->debug << "A^T      =\n" << aMatT   << endl;

            // ==============================================
            //         using the matrix directly
            // ==============================================

#ifdef DO_TIMING
            auto start1 = std::chrono::high_resolution_clock::now();
#endif
            for (size_t r = 0; r < nruns; ++r)
            {
                // the reference solution
                computeMatrixProductMatVec(aMat, vVec, cVecRef);
            }

#ifdef DO_TIMING
            auto stop1 = std::chrono::high_resolution_clock::now();
            this->debug << "Runtime of the reference matrix-vector product using A  : "
            << std::chrono::duration<double, std::milli>(stop1 - start1).count() << " ms"
            << endl;
#endif

/*
#ifdef DO_TIMING
            auto start2 = std::chrono::high_resolution_clock::now();
#endif
            for (size_t r = 0; r < nruns; ++r)
            {
                // the data structure-provided reference solution
                gmx::matrixProductMatrixVectorRef(aMat, vVec, cVec);
            }
#ifdef DO_TIMING
            auto stop2 = std::chrono::high_resolution_clock::now();
            this->debug << "Runtime of the scalar    matrix-vector product using A  : "
            << std::chrono::duration<double, std::milli>(stop2 - start2).count() << " ms"
            << endl;
#endif

#ifdef DO_TIMING
            auto start3 = std::chrono::high_resolution_clock::now();
#endif
            for (size_t r = 0; r < nruns; ++r)
            {
                // the data structure-provided reference solution
                gmx::matrixProductMatrixVectorSimd(aMat, vVec, cVecSimd);
            }
#ifdef DO_TIMING
            auto stop3 = std::chrono::high_resolution_clock::now();
            this->debug << "Runtime of the SIMD      matrix-vector product using A  : "
            << std::chrono::duration<double, std::milli>(stop3 - start3).count() << " ms"
            << endl;
#endif
*/
            this->debug << "matrix-vector product using A:" << endl;
            this->debug << "c (scalar) ="   << cVec         << endl;
            this->debug << "c (SIMD)   ="   << cVecSimd     << endl;
            this->debug << "c (ref.)   ="   << cVecRef      << endl;

/*
            EXPECT_NEAR(maxAbsDeviationVectorElements(cVec,     cVecRef), 0, tolDiff)
            << "Scalar vector matrix product directly using the matrix as input differs from the reference solution." << std::endl;
            EXPECT_NEAR(maxAbsDeviationVectorElements(cVecSimd, cVecRef), 0, tolDiff)
            << "SIMD   vector matrix product directly using the matrix as input differs from the reference solution." << std::endl;
*/

            std::fill(    cVec.begin(),     cVec.end(), static_cast<float_type>(0));
            std::fill( cVecRef.begin(),  cVecRef.end(), static_cast<float_type>(0));
            std::fill(cVecSimd.begin(), cVecSimd.end(), static_cast<float_type>(0));

            // ==============================================
            //         using the transpose matrix
            // ==============================================

#ifdef DO_TIMING
            auto start4 = std::chrono::high_resolution_clock::now();
#endif
            for (size_t r = 0; r < nruns; ++r)
            {
                // the reference solution
                computeMatrixProductMatTVec(aMatT, vVec, cVecRef);
            }
#ifdef DO_TIMING
            auto stop4 = std::chrono::high_resolution_clock::now();
            this->debug << "Runtime of the reference matrix-vector product using A^T: "
            << std::chrono::duration<double, std::milli>(stop4 - start4).count() << " ms"
            << endl;
#endif

/*
#ifdef DO_TIMING
            auto start5 = std::chrono::high_resolution_clock::now();
#endif
            for (size_t r = 0; r < nruns; ++r)
            {
                // the data structure-provided reference solution
                gmx::matrixProductMatrixTVectorRef(aMatT, vVec, cVec);
            }
#ifdef DO_TIMING
            auto stop5 = std::chrono::high_resolution_clock::now();
            this->debug << "Runtime of the scalar    matrix-vector product using A^T: "
            << std::chrono::duration<double, std::milli>(stop5 - start5).count() << " ms"
            << endl;
#endif

#ifdef DO_TIMING
            auto start6 = std::chrono::high_resolution_clock::now();
#endif
            for (size_t r = 0; r < nruns; ++r)
            {
                // the reference solution
                gmx::matrixProductMatrixTVectorSimd(aMatT, vVec, cVecSimd);
            }
#ifdef DO_TIMING
            auto stop6 = std::chrono::high_resolution_clock::now();
            this->debug << "Runtime of the SIMD      matrix-vector product using A^T: "
            << std::chrono::duration<double, std::milli>(stop6 - start6).count() << " ms"
            << endl;
#endif
*/
            this->debug << "matrix-vector product using A^T:" << endl;
            this->debug << "c (scalar) ="   << cVec           << endl;
            this->debug << "c (SIMD)   ="   << cVecSimd       << endl;
            this->debug << "c (ref.)   ="   << cVecRef        << endl;

/*
            EXPECT_NEAR(maxAbsDeviationVectorElements(cVec,     cVecRef), 0, tolDiff)
            << "Scalar matrix vector product using the transpose matrix as input deviates from the reference solution." << std::endl;
            EXPECT_NEAR(maxAbsDeviationVectorElements(cVecSimd, cVecRef), 0, tolDiff)
            << "SIMD   matrix vector product using the transpose matrix as input deviates from the reference solution." << std::endl;
*/

        }
    }

}

TYPED_TEST(ParamMdspanMathTest, MatrixProduct)
{
    typedef TypeParam float_type;

    // tolerance needs to be that high for larger matrices
    const float_type tolFac  = 100;
    const float_type tolMax  = 0.1;
    float_type       tolDiff = static_cast<float_type>(tolFac) / static_cast<float_type>(pow(std::numeric_limits<float_type>::digits, 2));
    if (tolDiff > tolMax)
    {
        tolDiff = tolMax;
    }

    typedef typename gmx::SimdSetup<float_type>::allocator_type allocator_type;

    typedef std::vector<float_type, allocator_type>                 vector_type;
    typedef basic_mdspan<float_type, extents<dynamic_extent, dynamic_extent>,
                         layout_left, accessor_basic<float_type> >   array_type;

    // testing the matrix product C  = A B or A B^T or A^T B or A^T B^T

    // n x m matrix sizes
    const std::vector<size_t> n = {1, 2, 3, 4, 5, 8, 10, 20, 40, 43, 240, 400};
    const std::vector<size_t> m = {1, 2, 3, 4, 5, 8, 10, 20, 40, 43, 240, 400};

    const size_t              nruns = 100;

    // testing the vector matrix product cVec = vVec aMat
    for (size_t i = 0; i < m.size(); ++i)
    {
        for (size_t j = 0; j < n.size(); ++j)
        {

            // length of the contracted dimension l
//            const std::vector<size_t> l = {std::max(minnm - 1, (size_t)1), m[i], (maxnm + 1), (maxnm * 2 + 1)};
            const std::vector<size_t> l = {m[i]};

            for (size_t k = 0; k < l.size(); ++k)
            {

                this->debug << "\nTesting matrix product " << " (float type " << typeid(float_type).name() << ") multiplying "
                << m[i] << " by " << l[k] << " matrix A and "
                << l[k] << " by " << n[j] << " matrix B resulting in "
                << m[i] << " by " << n[j] << " matrix C ..." << endl;

                vector_type  aBuf(m[i] * l[k], static_cast<float_type>(0));
                vector_type  aBufT(l[k] * m[i], static_cast<float_type>(0));
                array_type   aMat(aBuf.data(), m[i],  l[k]);
                array_type   aMatT(aBufT.data(), l[k], m[i]);

                vector_type  bBuf(l[k] * n[j], static_cast<float_type>(0));
                vector_type  bBufT(n[j] * l[k], static_cast<float_type>(0));
                array_type   bMat(bBuf.data(), l[k],  n[j]);
                array_type   bMatT(bBufT.data(), n[j], l[k]);

                vector_type  cBuf(m[i] * n[j], static_cast<float_type>(0));
                vector_type  cBufRef(m[i] * n[j], static_cast<float_type>(0));
                vector_type  cBufSimd(m[i] * n[j], static_cast<float_type>(0));
                array_type   cMat(cBuf.data(), m[i],  n[j]);
                array_type   cMatRef(cBufRef.data(), m[i], n[j]);
                array_type   cMatSimd(cBufSimd.data(), m[i], n[j]);

                // make sure that the result fits by filling the input vector and array with non-uniform values
                // a[i][j] = i / C + j / C
                // C = series y_i = sum_i^n x_i with x_i = x_{i-1} + 1,  n * (n + 1) / 2
                // aMat(i, j) = vVec[i] * j / C'
                // C = series y_i = sum_i^n x_i with x_i = x_{i-1} + 1,  m * (m + 1) / 2
                const float_type normfaci = 2.0 * std::sqrt(static_cast<float_type>(l[k] * (l[k] + 1))) / static_cast<float_type>(m[i] * (m[i] + 1));
                const float_type normfacj = 2.0 * std::sqrt(static_cast<float_type>(l[k] * (l[k] + 1))) / static_cast<float_type>(n[j] * (n[j] + 1));

                for (size_t ii = 0; ii < m[i]; ++ii)
                {
                    for (size_t jj = 0; jj < n[j]; ++jj)
                    {
                        for (size_t kk = 0; kk < l[k]; ++kk)
                        {
                            const float_type normfack = 100.0 / static_cast<float_type>(kk + jj + ii + 1);
                            aMat(ii, kk)  = normfaci * normfack;
                            aMatT(kk, ii) = normfaci * normfack;
                            bMat(kk, jj)  = normfacj * normfack;
                            bMatT(jj, kk) = normfacj * normfack;
                        }
                    }
                }

                this->debug << "Test input data:"     << endl;
                this->debug << "A   (address " << &aMat  << ") =\n" << aMat  << endl;
                this->debug << "A^T (address " << &aMatT << ") =\n" << aMatT << endl;
                this->debug << "B   (address " << &bMat  << ") =\n" << bMat  << endl;
                this->debug << "B^T (address " << &bMatT << ") =\n" << bMatT << endl;

                // ==============================================
                //         C = A B
                // ==============================================

#ifdef DO_TIMING
                auto start1 = std::chrono::high_resolution_clock::now();
#endif
                for (size_t r = 0; r < nruns; ++r)
                {
                    // the reference solution
                    computeMatrixProduct(aMat, bMat, cMatRef);
                }

#ifdef DO_TIMING
                auto stop1 = std::chrono::high_resolution_clock::now();
                this->debug << "Runtime of the reference matrix product A   B  : "
                << std::chrono::duration<double, std::milli>(stop1 - start1).count() << " ms"
                << endl;
#endif

/*
#ifdef DO_TIMING
                auto start2 = std::chrono::high_resolution_clock::now();
#endif
                for (size_t r = 0; r < nruns; ++r)
                {
                    // the data structure-provided reference solution
                    gmx::matrixProductRef(aMat, bMat, cMat);
                }
#ifdef DO_TIMING
                auto stop2 = std::chrono::high_resolution_clock::now();
                this->debug << "Runtime of the scalar    matrix product A   B  : "
                << std::chrono::duration<double, std::milli>(stop2 - start2).count() << " ms"
                << endl;
#endif

#ifdef DO_TIMING
                auto start3 = std::chrono::high_resolution_clock::now();
#endif
                for (size_t r = 0; r < nruns; ++r)
                {
                    // the data structure-provided SIMD solution
                    gmx::matrixProductSimd(aMat, bMat, cMatSimd);
                }
#ifdef DO_TIMING
                auto stop3 = std::chrono::high_resolution_clock::now();
                this->debug << "Runtime of the SIMD      matrix product A   B  : "
                << std::chrono::duration<double, std::milli>(stop3 - start3).count() << " ms"
                << endl;
#endif
*/

                this->debug << "matrix product result:"       << endl;
                this->debug << "c (scalar) =\n"   << cMat     << endl;
                this->debug << "c (SIMD)   =\n"   << cMatSimd << endl;
                this->debug << "c (ref.)   =\n"   << cMatRef  << endl;

/*
                EXPECT_NEAR(maxAbsDeviationMatrixElements(cMat,     cMatRef), 0, tolDiff)
                << "Scalar matrix product differs from the reference solution." << std::endl;
                EXPECT_NEAR(maxAbsDeviationMatrixElements(cMatSimd, cMatRef), 0, tolDiff)
                << "SIMD   matrix product differs from the reference solution." << std::endl;
*/
                // ==============================================
                //         C = A B^T
                // ==============================================

#ifdef DO_TIMING
                auto start4 = std::chrono::high_resolution_clock::now();
#endif
                for (size_t r = 0; r < nruns; ++r)
                {
                    // the reference solution
                    computeMatrixProductMatrixMatrixT(aMat, bMatT, cMatRef);
                }

#ifdef DO_TIMING
                auto stop4 = std::chrono::high_resolution_clock::now();
                this->debug << "Runtime of the reference matrix product A   B^T: "
                << std::chrono::duration<double, std::milli>(stop4 - start4).count() << " ms"
                << endl;
#endif

/*
#ifdef DO_TIMING
                auto start5 = std::chrono::high_resolution_clock::now();
#endif
                for (size_t r = 0; r < nruns; ++r)
                {
                    // the data structure-provided reference solution
                    gmx::matrixProductMatrixMatrixTRef(aMat, bMatT, cMat);
                }
#ifdef DO_TIMING
                auto stop5 = std::chrono::high_resolution_clock::now();
                this->debug << "Runtime of the scalar    matrix product A   B^T: "
                << std::chrono::duration<double, std::milli>(stop5 - start5).count() << " ms"
                << endl;
#endif

#ifdef DO_TIMING
                auto start6 = std::chrono::high_resolution_clock::now();
#endif
                for (size_t r = 0; r < nruns; ++r)
                {
                    // the data structure-provided reference solution
                    gmx::matrixProductMatrixMatrixTSimd(aMat, bMatT, cMatSimd);
                }
#ifdef DO_TIMING
                auto stop6 = std::chrono::high_resolution_clock::now();
                this->debug << "Runtime of the SIMD      matrix product A   B^T: "
                << std::chrono::duration<double, std::milli>(stop6 - start6).count() << " ms"
                << endl;
#endif
*/

                this->debug << "matrix product result:"       << endl;
                this->debug << "c (scalar) =\n"   << cMat     << endl;
                this->debug << "c (SIMD)   =\n"   << cMatSimd << endl;
                this->debug << "c (ref.)   =\n"   << cMatRef  << endl;

/*
                EXPECT_NEAR(maxAbsDeviationMatrixElements(cMat,     cMatRef), 0, tolDiff)
                << "Scalar matrix product differs from the reference solution." << std::endl;
                EXPECT_NEAR(maxAbsDeviationMatrixElements(cMatSimd, cMatRef), 0, tolDiff)
                << "SIMD   matrix product differs from the reference solution." << std::endl;
*/
                // ==============================================
                //         C = A^T B
                // ==============================================

#ifdef DO_TIMING
                auto start7 = std::chrono::high_resolution_clock::now();
#endif
                for (size_t r = 0; r < nruns; ++r)
                {
                    // the reference solution
                    computeMatrixProductMatrixTMatrix(aMatT, bMat, cMatRef);
                }

#ifdef DO_TIMING
                auto stop7 = std::chrono::high_resolution_clock::now();
                this->debug << "Runtime of the reference matrix product A^T B  : "
                << std::chrono::duration<double, std::milli>(stop7 - start7).count() << " ms"
                << endl;
#endif

/*
#ifdef DO_TIMING
                auto start8 = std::chrono::high_resolution_clock::now();
#endif
                for (size_t r = 0; r < nruns; ++r)
                {
                    // the data structure-provided reference solution
                    gmx::matrixProductMatrixTMatrixRef(aMatT, bMat, cMat);
                }
#ifdef DO_TIMING
                auto stop8 = std::chrono::high_resolution_clock::now();
                this->debug << "Runtime of the scalar    matrix product A^T B  : "
                << std::chrono::duration<double, std::milli>(stop8 - start8).count() << " ms"
                << endl;
#endif

#ifdef DO_TIMING
                auto start9 = std::chrono::high_resolution_clock::now();
#endif
                for (size_t r = 0; r < nruns; ++r)
                {
                    // the data structure-provided reference solution
                    gmx::matrixProductMatrixTMatrixSimd(aMatT, bMat, cMatSimd);
                }
#ifdef DO_TIMING
                auto stop9 = std::chrono::high_resolution_clock::now();
                this->debug << "Runtime of the SIMD      matrix product A^T B  : "
                << std::chrono::duration<double, std::milli>(stop9 - start9).count() << " ms"
                << endl;
#endif
*/
                this->debug << "matrix product result:"       << endl;
                this->debug << "c (scalar) =\n"   << cMat     << endl;
                this->debug << "c (SIMD)   =\n"   << cMatSimd << endl;
                this->debug << "c (ref.)   =\n"   << cMatRef  << endl;

/*
                EXPECT_NEAR(maxAbsDeviationMatrixElements(cMat,     cMatRef), 0, tolDiff)
                << "Scalar matrix product differs from the reference solution." << std::endl;
                EXPECT_NEAR(maxAbsDeviationMatrixElements(cMatSimd, cMatRef), 0, tolDiff)
                << "SIMD   matrix product differs from the reference solution." << std::endl;
*/
                // ==============================================
                //         C = A^T B^T
                // ==============================================

#ifdef DO_TIMING
                auto start10 = std::chrono::high_resolution_clock::now();
#endif
                for (size_t r = 0; r < nruns; ++r)
                {
                    // the reference solution
                    computeMatrixProductMatrixTMatrixT(aMatT, bMatT, cMatRef);
                }

#ifdef DO_TIMING
                auto stop10 = std::chrono::high_resolution_clock::now();
                this->debug << "Runtime of the reference matrix product A^T B^T: "
                << std::chrono::duration<double, std::milli>(stop10 - start10).count() << " ms"
                << endl;
#endif

/*
#ifdef DO_TIMING
                auto start11 = std::chrono::high_resolution_clock::now();
#endif
                for (size_t r = 0; r < nruns; ++r)
                {
                    // the data structure-provided reference solution
                    gmx::matrixProductMatrixTMatrixTRef(aMatT, bMatT, cMat);
                }
#ifdef DO_TIMING
                auto stop11 = std::chrono::high_resolution_clock::now();
                this->debug << "Runtime of the scalar    matrix product A^T B^T: "
                << std::chrono::duration<double, std::milli>(stop11 - start11).count() << " ms"
                << endl;
#endif

#ifdef DO_TIMING
                auto start12 = std::chrono::high_resolution_clock::now();
#endif
                for (size_t r = 0; r < nruns; ++r)
                {
                    // the data structure-provided reference solution
                    gmx::matrixProductMatrixTMatrixTSimd(aMatT, bMatT, cMatSimd);
                }
#ifdef DO_TIMING
                auto stop12 = std::chrono::high_resolution_clock::now();
                this->debug << "Runtime of the SIMD      matrix product A^T B^T: "
                << std::chrono::duration<double, std::milli>(stop12 - start12).count() << " ms"
                << endl;
#endif
*/
                this->debug << "matrix product result:"       << endl;
                this->debug << "c (scalar) =\n"   << cMat     << endl;
                this->debug << "c (SIMD)   =\n"   << cMatSimd << endl;
                this->debug << "c (ref.)   =\n"   << cMatRef  << endl;

/*
                EXPECT_NEAR(maxAbsDeviationMatrixElements(cMat,     cMatRef), 0, tolDiff)
                << "Scalar matrix product differs from the reference solution." << std::endl;
                EXPECT_NEAR(maxAbsDeviationMatrixElements(cMatSimd, cMatRef), 0, tolDiff)
                << "SIMD   matrix product differs from the reference solution." << std::endl;
*/
            }
        }
    }

}

} // namespace

} // namespace irreg_array_math_test

} // namespace gmx
