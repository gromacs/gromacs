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
#include "gromacs/utility/data_structures/flat_irreg_array_4d.h"
#include "gromacs/utility/data_structures/irreg_array_4d.h"

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
class ParamIrregArrayMathTest : public IrregArrayMathTest,
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
    index_type m = a.getLength1();
    // number of columns    of matrix A
    index_type l = a.getLength2(0);
    // number of rows       of matrix B
    index_type n = b.getLength2(0);

    // initialize matrix c one row vector at a time
    for (index_type i = 0; i < m; ++i)
    {
        std::fill(c(i), c(i) + static_cast<ptrdiff_t>(n), (float_type)0);
    }

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
    index_type m = a.getLength2(0);
    // number of columns    of matrix A
    index_type l = a.getLength1();
    // number of rows       of matrix B
    index_type n = b.getLength2(0);

    // initialize matrix c one row vector at a time
    for (index_type i = 0; i < m; ++i)
    {
        std::fill(c(i), c(i) + static_cast<ptrdiff_t>(n), (float_type)0);
    }

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
    index_type m = a.getLength2(0);
    // number of columns    of matrix A
    index_type l = a.getLength1();
    // number of rows       of matrix B
    index_type n = b.getLength1();

    // initialize matrix c one row vector at a time
    for (index_type i = 0; i < m; ++i)
    {
        std::fill(c(i), c(i) + static_cast<ptrdiff_t>(n), (float_type)0);
    }

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
    index_type m = a.getLength1();
    // number of columns    of matrix A
    index_type l = a.getLength2(0);
    // number of rows       of matrix B
    index_type n = b.getLength1();

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
    const typename TMat::index_type M = aT.getLength2(0);
    // number of rows       of the matrix
    const typename TMat::index_type N = aT.getLength1();

    std::fill(c.begin(), c.end(), (float_type)0);

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
    const typename TMat::index_type M = a.getLength1();
    // number of columns    of the matrix
    const typename TMat::index_type N = a.getLength2(0);

    std::fill(c.begin(), c.end(), (float_type)0);

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
    const typename TMat::index_type M = a.getLength1();
    // number of columns    of the matrix
    const typename TMat::index_type N = a.getLength2(0);

    std::fill(c.begin(), c.end(), (float_type)0);

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
    const typename TMat::index_type N = aT.getLength1();
    // number of columns    of the matrix A
    const typename TMat::index_type M = aT.getLength2(0);

    std::fill(c.begin(), c.end(), (float_type)0);

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
    for (size_t i = 0; i < a.getLength1(); ++i)
    {
        for (size_t j = 0; j < a.getLength2(i); ++j)
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
TYPED_TEST_CASE(ParamIrregArrayMathTest, TestFloatTypes);


TYPED_TEST(ParamIrregArrayMathTest, VectorMatrixProduct)
{
    typedef TypeParam float_type;

    // tolerance needs to be that high for larger matrices
    float_type tolDiff = (float_type)100 / static_cast<float_type>(pow(std::numeric_limits<float_type>::digits, 2));
    if (tolDiff > 0.1)
    {
        tolDiff = 0.1;
    }

    typedef typename gmx::SimdSetup<float_type>::allocator_type allocator_type;
//    typedef gmx::AlignedAllocator<float_type> allocator_type;
//    typedef std::allocator<float_type>                             allocator_type;

    typedef std::vector<float_type, allocator_type>                 vector_type;
    typedef gmx::IrregArray2D<float_type, allocator_type>            array_type;

    // testing the vector matrix product cVec^T  = vVec^T  aMat
    //                                   (1 x m) = (1 x n) (n x m)

    // n x m matrix sizes
    const std::vector<size_t> n = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 60, 67};
    const std::vector<size_t> m = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 60, 73};

//    const std::vector<size_t> m = {5};
//    const std::vector<size_t> n = {4};

    const size_t nruns = 1;

    // testing the vector matrix product cVec = vVec aMat
    for (size_t i = 0; i < m.size(); ++i)
    {
        vector_type vVec(m[i], (float_type)0);
        for (size_t j = 0; j < n.size(); ++j)
        {

            this->debug << "\nTesting vector-matrix product " << " (float type " << typeid(float_type).name() << ") for a " << m[i] << " x " << n[j] << " matrix ..." << endl;

            vector_type  cVec(    n[j], (float_type)0);
            vector_type  cVecRef( n[j], (float_type)0);
            vector_type  cVecSimd(n[j], (float_type)0);

            array_type   aMat(0, m[i] - 1, 0, n[j] - 1, (float_type)0);
            array_type   aMatT(0, n[j] - 1, 0, m[i] - 1, (float_type)0);

            // make sure that the result fits by filling the input vector and array with non-uniform values
            // vVec[i] = i / C
            // C = series y_i = sum_i^n x_i with x_i = x_{i-1} + 1,  n * (n + 1) / 2
            // aMat(i, j) = vVec[i] * j / C'
            // C = series y_i = sum_i^n x_i with x_i = x_{i-1} + 1,  m * (m + 1) / 2
            float_type normfac = 10.0 * std::sqrt(static_cast<float_type>(m[i] * n[j])) / static_cast<float_type>(n[j] * (n[j] + 1));
            for (size_t k = 0; k < m[i]; ++k)
            {
                vVec[k] = static_cast<float_type>(k + 1) * normfac;
            }

            for (size_t k = 0; k < m[i]; ++k)
            {
                float_type normfacRow = 2.0 / static_cast<float_type>(m[i] * (m[i] + 1));
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

            this->debug << "vector-matrix product using A:" << endl;
            this->debug << "c (scalar) ="   << cVec         << endl;
            this->debug << "c (SIMD)   ="   << cVecSimd     << endl;
            this->debug << "c (ref.)   ="   << cVecRef      << endl;

            EXPECT_NEAR(maxAbsDeviationVectorElements(cVec,     cVecRef), 0, tolDiff)
            << "Scalar vector matrix product directly using the matrix as input differs from the reference solution." << std::endl;
            EXPECT_NEAR(maxAbsDeviationVectorElements(cVecSimd, cVecRef), 0, tolDiff)
            << "SIMD   vector matrix product directly using the matrix as input differs from the reference solution." << std::endl;

            std::fill(    cVec.begin(),     cVec.end(), (float_type)0);
            std::fill( cVecRef.begin(),  cVecRef.end(), (float_type)0);
            std::fill(cVecSimd.begin(), cVecSimd.end(), (float_type)0);

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

            this->debug << "vector-matrix product using A^T:" << endl;
            this->debug << "c (scalar) ="   << cVec           << endl;
            this->debug << "c (SIMD)   ="   << cVecSimd       << endl;
            this->debug << "c (ref.)   ="   << cVecRef        << endl;

            EXPECT_NEAR(maxAbsDeviationVectorElements(cVec,     cVecRef), 0, tolDiff)
            << "Scalar vector matrix product using the transpose matrix as input deviates from the reference solution." << std::endl;
            EXPECT_NEAR(maxAbsDeviationVectorElements(cVecSimd, cVecRef), 0, tolDiff)
            << "SIMD   vector matrix product using the transpose matrix as input deviates from the reference solution." << std::endl;

        }
    }

}

TYPED_TEST(ParamIrregArrayMathTest, MatrixVectorProduct)
{
    typedef TypeParam float_type;

    // tolerance needs to be that high for larger matrices
    float_type tolDiff = (float_type)100 / static_cast<float_type>(pow(std::numeric_limits<float_type>::digits, 2));
    if (tolDiff > 0.1)
    {
        tolDiff = 0.1;
    }

    typedef typename gmx::SimdSetup<float_type>::allocator_type allocator_type;
//    typedef gmx::AlignedAllocator<float_type> allocator_type;
//    typedef std::allocator<float_type>                             allocator_type;

    typedef std::vector<float_type, allocator_type>                 vector_type;
    typedef gmx::IrregArray2D<float_type, allocator_type>            array_type;

    // testing the vector matrix product cVec^T  =  aMat    vVec^T
    //                                   (m x 1) = (m x n) (n x 1)

    // m x n matrix sizes
    const std::vector<size_t> m = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 60, 67};
    const std::vector<size_t> n = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 60, 73};

/*
    const std::vector<size_t> m = {9};
    const std::vector<size_t> n = {5};
 */

    const size_t nruns = 1;

    // testing the vector matrix product cVec = vVec aMat
    for (size_t i = 0; i < m.size(); ++i)
    {
        for (size_t j = 0; j < n.size(); ++j)
        {
            this->debug << "\nTesting matrix-vector product " << " (float type " << typeid(float_type).name() << ") for a " << m[i] << " x " << n[j] << " matrix ..." << endl;

            vector_type  vVec(n[j], (float_type)0);

            vector_type  cVec(    m[i], (float_type)0);
            vector_type  cVecRef( m[i], (float_type)0);
            vector_type  cVecSimd(m[i], (float_type)0);

            array_type   aMat( 0, m[i] - 1, 0, n[j] - 1, (float_type)0);
            array_type   aMatT(0, n[j] - 1, 0, m[i] - 1, (float_type)0);

            // make sure that the result fits by filling the input vector and array with non-uniform values
            // vVec[i] = i / C
            // C = series y_i = sum_i^n x_i with x_i = x_{i-1} + 1,  n * (n + 1) / 2
            // aMat(i, j) = vVec[i] * j / C'
            // C = series y_i = sum_i^n x_i with x_i = x_{i-1} + 1,  m * (m + 1) / 2
            float_type normfac = 10.0 * std::sqrt(static_cast<float_type>(m[i] * n[j])) / static_cast<float_type>(n[j] * (n[j] + 1));
            for (size_t k = 0; k < n[j]; ++k)
            {
                vVec[k] = static_cast<float_type>(k + 1) * normfac;
            }

            for (size_t k = 0; k < n[j]; ++k)
            {
                float_type normfacRow = 2.0 / static_cast<float_type>(m[i] * (m[i] + 1));
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

            this->debug << "matrix-vector product using A:" << endl;
            this->debug << "c (scalar) ="   << cVec         << endl;
            this->debug << "c (SIMD)   ="   << cVecSimd     << endl;
            this->debug << "c (ref.)   ="   << cVecRef      << endl;

            EXPECT_NEAR(maxAbsDeviationVectorElements(cVec,     cVecRef), 0, tolDiff)
            << "Scalar vector matrix product directly using the matrix as input differs from the reference solution." << std::endl;
            EXPECT_NEAR(maxAbsDeviationVectorElements(cVecSimd, cVecRef), 0, tolDiff)
            << "SIMD   vector matrix product directly using the matrix as input differs from the reference solution." << std::endl;

            std::fill(    cVec.begin(),     cVec.end(), (float_type)0);
            std::fill( cVecRef.begin(),  cVecRef.end(), (float_type)0);
            std::fill(cVecSimd.begin(), cVecSimd.end(), (float_type)0);

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

            this->debug << "matrix-vector product using A^T:" << endl;
            this->debug << "c (scalar) ="   << cVec           << endl;
            this->debug << "c (SIMD)   ="   << cVecSimd       << endl;
            this->debug << "c (ref.)   ="   << cVecRef        << endl;

            EXPECT_NEAR(maxAbsDeviationVectorElements(cVec,     cVecRef), 0, tolDiff)
            << "Scalar matrix vector product using the transpose matrix as input deviates from the reference solution." << std::endl;
            EXPECT_NEAR(maxAbsDeviationVectorElements(cVecSimd, cVecRef), 0, tolDiff)
            << "SIMD   matrix vector product using the transpose matrix as input deviates from the reference solution." << std::endl;

        }
    }

}

TYPED_TEST(ParamIrregArrayMathTest, MatrixProduct)
{
    typedef TypeParam float_type;

    // tolerance needs to be that high for larger matrices
    float_type tolDiff = (float_type)100 / static_cast<float_type>(pow(std::numeric_limits<float_type>::digits, 2));
    if (tolDiff > 0.1)
    {
        tolDiff = 0.1;
    }

    typedef typename gmx::SimdSetup<float_type>::allocator_type allocator_type;
//    typedef gmx::AlignedAllocator<float_type> allocator_type;
//    typedef std::allocator<float_type>                             allocator_type;

    typedef gmx::IrregArray2D<float_type, allocator_type>            array_type;

    // testing the matrix product C  = A B or A B^T or A^T B or A^T B^T

    // n x m matrix sizes
    const std::vector<size_t> n = {1, 2, 3, 4, 5, 8, 10, 20, 40, 43};
    const std::vector<size_t> m = {1, 2, 3, 4, 5, 8, 10, 20, 40, 43};

    const size_t              nruns = 1;

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

                array_type     aMat(    0, m[i] - 1, 0, l[k] - 1, (float_type)0);
                array_type     aMatT(   0, l[k] - 1, 0, m[i] - 1, (float_type)0);
                array_type     bMat(    0, l[k] - 1, 0, n[j] - 1, (float_type)0);
                array_type     bMatT(   0, n[j] - 1, 0, l[k] - 1, (float_type)0);
                array_type     cMat(    0, m[i] - 1, 0, n[j] - 1, (float_type)0);
                array_type     cMatRef( 0, m[i] - 1, 0, n[j] - 1, (float_type)0);
                array_type     cMatSimd(0, m[i] - 1, 0, n[j] - 1, (float_type)0);

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

                this->debug << "matrix product result:"       << endl;
                this->debug << "c (scalar) =\n"   << cMat     << endl;
                this->debug << "c (SIMD)   =\n"   << cMatSimd << endl;
                this->debug << "c (ref.)   =\n"   << cMatRef  << endl;

                EXPECT_NEAR(maxAbsDeviationMatrixElements(cMat,     cMatRef), 0, tolDiff)
                << "Scalar matrix product differs from the reference solution." << std::endl;
                EXPECT_NEAR(maxAbsDeviationMatrixElements(cMatSimd, cMatRef), 0, tolDiff)
                << "SIMD   matrix product differs from the reference solution." << std::endl;

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

                this->debug << "matrix product result:"       << endl;
                this->debug << "c (scalar) =\n"   << cMat     << endl;
                this->debug << "c (SIMD)   =\n"   << cMatSimd << endl;
                this->debug << "c (ref.)   =\n"   << cMatRef  << endl;

                EXPECT_NEAR(maxAbsDeviationMatrixElements(cMat,     cMatRef), 0, tolDiff)
                << "Scalar matrix product differs from the reference solution." << std::endl;
                EXPECT_NEAR(maxAbsDeviationMatrixElements(cMatSimd, cMatRef), 0, tolDiff)
                << "SIMD   matrix product differs from the reference solution." << std::endl;

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

                this->debug << "matrix product result:"       << endl;
                this->debug << "c (scalar) =\n"   << cMat     << endl;
                this->debug << "c (SIMD)   =\n"   << cMatSimd << endl;
                this->debug << "c (ref.)   =\n"   << cMatRef  << endl;

                EXPECT_NEAR(maxAbsDeviationMatrixElements(cMat,     cMatRef), 0, tolDiff)
                << "Scalar matrix product differs from the reference solution." << std::endl;
                EXPECT_NEAR(maxAbsDeviationMatrixElements(cMatSimd, cMatRef), 0, tolDiff)
                << "SIMD   matrix product differs from the reference solution." << std::endl;

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

                this->debug << "matrix product result:"       << endl;
                this->debug << "c (scalar) =\n"   << cMat     << endl;
                this->debug << "c (SIMD)   =\n"   << cMatSimd << endl;
                this->debug << "c (ref.)   =\n"   << cMatRef  << endl;

                EXPECT_NEAR(maxAbsDeviationMatrixElements(cMat,     cMatRef), 0, tolDiff)
                << "Scalar matrix product differs from the reference solution." << std::endl;
                EXPECT_NEAR(maxAbsDeviationMatrixElements(cMatSimd, cMatRef), 0, tolDiff)
                << "SIMD   matrix product differs from the reference solution." << std::endl;
            }
        }
    }

}


/*! \brief check the results of LU decomposition with row pivoting \f$ \mathbf{A} = \mathbf{P}^{-1}\,\mathbf{L}\,\mathbf{U} \f$

    \tparam      T        data type to be used in the calculation
    \tparam      Alloc    data type to be used in the calculation

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

    The upper and lower triangular matrices \f$ \mathbf{L} \f$ and \f$ \mathbf{P} \f$ stored
    in the upper right and lower left. The diagonal entries of \f$ \mathbf{L} \f$ are omitted.
    These entries are all set equal to \f$ 1 \f$, which is a common convention to obtain a
    unique decomposition.

    \param[in]       aMat          \f$ n \times n\f$ matrix \f$ \mathbf{A} \f$ to be decomposed
    \param[in]       pVec          vector representation of the \f$ n \times n\f$ row pivot matrix \f$ \mathbf{P} \f$, entry \f$ i \f$ indicates the new row index of the current row \f$ i \f$
    \param[in]       luMatPacked   the upper and lower triangular matrices \f$ \mathbf{L} \f$ and \f$ \mathbf{P} \f$  stored in the upper right and lower left, omitting the diagonal entries of \f$ \mathbf{L} \f$, which are all equal to \f$ 1 \f$
    \param[in,out]   debug         ostrema for debug output
 */
template<typename TMat, typename TVec>
inline bool checkLUPDecomposition(const TMat &aMat, const TVec &pVec, const TMat &luMatPacked, std::ostream &debug)
{
    typedef           TMat                    array_type;
    typedef typename  array_type::index_type  index_type;
    typedef typename  array_type::value_type  float_type;

    // tolerance needs to be that high for larger matrices
    float_type tolDiff = (float_type)10.0 / static_cast<float_type>(pow(std::numeric_limits<float_type>::digits, 2));
    if (tolDiff > 0.01)
    {
        tolDiff = 0.01;
    }

    const index_type n = aMat.getLength1();

    // create separate matrices for verifying A = B = P^{-1] L U
    array_type  bMat( 0, n - 1, 0, n - 1, (float_type)0);
    array_type  ipMat(0, n - 1, 0, n - 1, (float_type)0);
    array_type  lMat( 0, n - 1, 0, n - 1, (float_type)0);
    array_type  uMat( 0, n - 1, 0, n - 1, (float_type)0);
    array_type  luMat(0, n - 1, 0, n - 1, (float_type)0);

    // expand the vector representation of the permutation matrix P
    // the inverse of  a permutation matrix is equal to its transpose
    for (index_type i = 0; i < n; ++i)
    {
        ipMat(pVec[i], i) = (float_type)1.0;
    }
    // unpack the lower and upper triangular matrices L and U
    for (index_type i = 0; i < n; ++i)
    {
        lMat(i, i) = (float_type)1.0;
        for (index_type j = 0; j < i; ++j)
        {
            lMat(i, j) = luMatPacked(i, j);
        }
        for (index_type j = i; j < n; ++j)
        {
            uMat(i, j) = luMatPacked(i, j);
        }
    }

    gmx::matrixProduct(lMat, uMat, luMat);
    gmx::matrixProduct(ipMat, luMat, bMat);

    debug << "LUP decomposition, input matrix A should equal B = P^{-1] L U:" << endl;
    debug << "P (vector representation) = "          << pVec        << endl;
    debug << "P^{-1} = P^{T} (expanded matrix)  =\n" << ipMat       << endl;
    debug << "LU packed =\n"                         << luMatPacked << endl;
    debug << "L         =\n"                         << lMat        << endl;
    debug << "U         =\n"                         << uMat        << endl;
    debug << "L U       =\n"                         << luMat       << endl;
    debug << "A         =\n"                         << aMat        << endl;
    debug << "B         =\n"                         << bMat        << endl;

    const float_type maxDiff = maxAbsDeviationMatrixElements(aMat, bMat);
    EXPECT_NEAR(maxDiff, 0, tolDiff)
    << "One or more elements of matrix A and its reconstruction from the LUP decomposition B = P^{-1} L U differ." << std::endl;

    // true if all entries of the original and reconstructed matrices match within tolDiff, false otherwise
    return (maxDiff < tolDiff);
}

TYPED_TEST(ParamIrregArrayMathTest, LUPDecomposition)
{
    typedef TypeParam float_type;

    typedef  typename gmx::SimdSetup<float_type>::allocator_type  allocator_type;
    typedef  gmx::IrregArray2D<float_type, allocator_type>            array_type;

    typedef  typename  array_type::index_type                                   index_type;
    typedef  typename  gmx::SimdSetup<index_type>::allocator_type     index_allocator_type;
    typedef  typename  std::vector<index_type, index_allocator_type>     index_vector_type;

    // n x n matrix sizes
    const std::vector<size_t> n = {1, 2, 3, 4, 5, 8, 10, 20, 40, 67};

    const size_t              nruns = 1;

    // testing the vector matrix product cVec = vVec aMat
    for (size_t i = 0; i < n.size(); ++i)
    {
        this->debug << "\nTesting LUP decomposition " << " (float type " << typeid(float_type).name() << ") of "
        << n[i] << " by " << n[i] << " matrix A" << endl;

        // the matrix to be decomposed
        array_type aMat(     0,  n[i] - 1,              0, n[i] - 1,              (float_type)0);
        // packed upper and lower triangular matrices, diagonal entries of L omitted (all equal to 1)
        array_type luMatRef( 0,  aMat.getLength1() - 1, 0, aMat.getLength1() - 1, (float_type)0);
        array_type luMatSimd(0,  aMat.getLength1() - 1, 0, aMat.getLength1() - 1, (float_type)0);
        // matrix for storing the reconstructed matrix A = P^{-1} L U
        array_type bMat( 0,  aMat.getLength1() - 1, 0, aMat.getLength1() - 1, (float_type)0);
        // vector representation of the permutation matrix P, original rows
        // of matrix A are mapped to new rows after pivoting row i -> row pVec[i]
        index_vector_type pVecRef( aMat.getLength1(), (index_type)0);
        index_vector_type pVecSimd(aMat.getLength1(), (index_type)0);

        // make sure that the result fits by filling the array with non-uniform values
        const float_type normfac1 = 2.0 * std::sqrt(static_cast<float_type>(n[i] * (n[i] + 1)) * static_cast<float_type>(i + 1));
        for (size_t ii = 0; ii < n[i]; ++ii)
        {
            for (size_t jj = 0; jj < n[i]; ++jj)
            {
                const float_type normfac2 = 10.0 / static_cast<float_type>(ii + jj + 1);
                aMat(ii, jj)  = normfac1 * normfac2;
            }
        }

        this->debug << "Test input data:" << endl;
        this->debug << "A   (address " << &aMat  << ") =\n" << aMat  << endl;

        bool successRef = false;
#ifdef DO_TIMING
        auto start1 = std::chrono::high_resolution_clock::now();
#endif
        for (size_t r = 0; r < nruns; ++r)
        {
            successRef = gmx::LUPDecompositionRef(aMat, pVecRef, luMatRef);
        }
#ifdef DO_TIMING
        auto stop1 = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of the scalar    LUP decomposition: "
        << std::chrono::duration<double, std::milli>(stop1 - start1).count() << " ms"
        << endl;
#endif

        bool successSimd = false;
#ifdef DO_TIMING
        auto start2 = std::chrono::high_resolution_clock::now();
#endif
        for (size_t r = 0; r < nruns; ++r)
        {
            successSimd = gmx::LUPDecompositionSimd(aMat, pVecSimd, luMatSimd);
        }
#ifdef DO_TIMING
        auto stop2 = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of the SIMD      LUP decomposition: "
        << std::chrono::duration<double, std::milli>(stop2 - start2).count() << " ms"
        << endl;
#endif

        EXPECT_TRUE(successRef) << "Scalar LUP decomposition returned failure." << std::endl;
        EXPECT_TRUE(checkLUPDecomposition(aMat, pVecRef,  luMatRef,  this->debug))
        << "Scalar LUP decomposition B = P^{-1} L U differs from the " << n[i] << "x" << n[i] << " input matrix A." << std::endl;

        EXPECT_TRUE(successSimd) << "SIMD LUP decomposition returned failure." << std::endl;
        EXPECT_TRUE(checkLUPDecomposition(aMat, pVecSimd,  luMatSimd,  this->debug))
        << "SIMD LUP decomposition B = P^{-1} L U differs from the " << n[i] << "x" << n[i] << " input matrix A." << std::endl;
    }

}

TYPED_TEST(ParamIrregArrayMathTest, MatrixInverse)
{
    typedef TypeParam float_type;

    typedef  typename gmx::SimdSetup<float_type>::allocator_type  allocator_type;
    typedef  gmx::IrregArray2D<float_type, allocator_type>            array_type;

    // tolerance needs to be that high for larger matrices
    float_type tolDiff = (float_type)10.0 / static_cast<float_type>(pow(std::numeric_limits<float_type>::digits, 2));
    if (tolDiff > 0.01)
    {
        tolDiff = 0.01;
    }

    // n x n matrix sizes
    const std::vector<size_t> n = {1, 2, 3, 4, 5, 8, 10, 20, 40, 67};

    const size_t              nruns = 1;

    // testing the vector matrix product cVec = vVec aMat
    for (size_t i = 0; i < n.size(); ++i)
    {
        // single precision doesn't suffice for inverting larger matrices, for less well behaved matrices, even n = 8 failed sometimes
        if (std::is_same<float_type, float>::value && n[i] > 40)
        {
            continue;
        }
        this->debug << "\nTesting matrix inversion via LUP decomposition " << " (float type " << typeid(float_type).name() << ") of "
        << n[i] << " by " << n[i] << " matrix A" << endl;

        // the matrix to be decomposed
        array_type aMat(             0,  n[i] - 1,              0, n[i] - 1,              (float_type)0);
        // A^{-1}
        array_type aMatInvScalar(    0,  aMat.getLength1() - 1, 0, aMat.getLength1() - 1, (float_type)0);
        array_type aMatInvSimd(      0,  aMat.getLength1() - 1, 0, aMat.getLength1() - 1, (float_type)0);
        // matrix for storing the A A^{-1}
        array_type bMat(             0,  aMat.getLength1() - 1, 0, aMat.getLength1() - 1, (float_type)0);
        // the identity matrix
        array_type iMat(             0,  aMat.getLength1() - 1, 0, aMat.getLength1() - 1, (float_type)0);
        for (size_t ii = 0; ii < n[i]; ++ii)
        {
            iMat(ii, ii) = 1;
        }

        // Make sure that the result fits by filling the array with non-uniform values.
        // The series used in other tests lead to matrices whose inversion fails in single precision.
        for (size_t ii = 0; ii < n[i]; ++ii)
        {
            for (size_t jj = 0; jj < n[i]; ++jj)
            {
                aMat(ii, jj) = std::cos(static_cast<float_type>(ii * jj * jj - ii * ii * jj + jj))
                    + static_cast<float_type>(ii)/static_cast<float_type>(n[i]);
            }
        }

        this->debug << "Test input data:" << endl;
        this->debug << "A   (address " << &aMat  << ") =\n" << aMat  << endl;

        bool successRef = false;
#ifdef DO_TIMING
        auto start1 = std::chrono::high_resolution_clock::now();
#endif
        for (size_t r = 0; r < nruns; ++r)
        {
            successRef = gmx::matrixInverseLUPRef(aMat, aMatInvScalar);
        }
#ifdef DO_TIMING
        auto stop1 = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of the scalar    matrix inversion: "
        << std::chrono::duration<double, std::milli>(stop1 - start1).count() << " ms"
        << endl;
#endif

        bool successSimd = false;
#ifdef DO_TIMING
        auto start2 = std::chrono::high_resolution_clock::now();
#endif
        for (size_t r = 0; r < nruns; ++r)
        {
            successSimd = gmx::matrixInverseLUPSimd(aMat, aMatInvSimd);
        }
#ifdef DO_TIMING
        auto stop2 = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of the SIMD      matrix inversion: "
        << std::chrono::duration<double, std::milli>(stop2 - start2).count() << " ms"
        << endl;
#endif

        gmx::matrixProduct(aMat, aMatInvScalar, bMat);
        this->debug << "A                 =\n" << aMat          << std::endl;
        this->debug << "A^{-1}   (Scalar) =\n" << aMatInvScalar << std::endl;
        this->debug << "A A^{-1} (Scalar) =\n" << bMat          << std::endl;

        EXPECT_TRUE(successRef) << "Scalar implementation of matrix inverse via LUP decomposition returned failure." << std::endl;
        EXPECT_NEAR(maxAbsDeviationMatrixElements(bMat, iMat), 0, tolDiff)
        << "Product A A^{-1} != I for the SIMD matrix inversion of the " << n[i] << "x" << n[i] << " input matrix A." << std::endl;


        gmx::matrixProduct(aMat, aMatInvSimd, bMat);
        this->debug << "A                 =\n" << aMat        << std::endl;
        this->debug << "A^{-1}   (SIMD)   =\n" << aMatInvSimd << std::endl;
        this->debug << "A A^{-1} (SIMD)   =\n" << bMat        << std::endl;

        EXPECT_TRUE(successSimd) << "SIMD implementation of matrix inverse via LUP decomposition returned failure." << std::endl;
        EXPECT_NEAR(maxAbsDeviationMatrixElements(bMat, iMat), 0, tolDiff)
        << "Product A A^{-1} != I for the SIMD matrix inversion of the " << n[i] << "x" << n[i] << " input matrix A." << std::endl;
    }

}


/*! \brief helper function for computing the determinant of a matrix via Leibniz' rule
           \f$ \left| \mathbf{A} \right| = \mathrm{det}\,\mathbf{A} \f$

    \tparam   TReal   floating point data type
    \tparam   TVec    vector type for storing TReal

    \param[in]   matrix   square matrix whose determinant is to be computed
 */
template<class TMat>
inline typename TMat::value_type computeDeterminant(const TMat &matrix)
{
    typedef typename TMat::value_type         float_type;
    typedef typename TMat::allocator_type allocator_type;

    EXPECT_TRUE(matrix.getLength1() == matrix.getLength2())
    << "array dimensions differ, need a square matrix for computing its determinant." << std::endl;

    // side length of the N x N
    const typename TMat::index_type N = matrix.getLength1();

    // helper array for storing permutations of the array indices
    gmx::IrregArray1D<size_t, allocator_type> permutation(0, N - 1, (size_t)0);
    for (typename TMat::index_type i = 0; i < N; ++i)
    {
        permutation[i] = i + 1;
    }

    long double result = 0;
    do
    {
        bool sign = true;
        for (typename TMat::index_type si = 1; si < N; ++si)
        {
            for (typename TMat::index_type sj = 0; sj < si; ++sj)
            {
                if (permutation[si] < permutation[sj])
                {
                    sign = (!sign);
                }
            }
        }

        long double prod = (sign ? (long double)1 : (long double)-1);
        for (typename TMat::index_type i = 0; i < N; ++i)
        {
            prod *= matrix(i, (permutation[i] - 1));
        }

        result += prod;
    }
    while (std::next_permutation(permutation.getArray(), permutation.getArray() + N));

    return static_cast<float_type>(result);
}

TYPED_TEST(ParamIrregArrayMathTest, MatrixDeterminant)
{
    typedef TypeParam float_type;

    typedef  typename gmx::SimdSetup<float_type>::allocator_type  allocator_type;
    typedef  gmx::IrregArray2D<float_type, allocator_type>            array_type;

    // tolerance needs to be that high for larger matrices
    float_type tolDiff = (float_type)10.0 / static_cast<float_type>(pow(std::numeric_limits<float_type>::digits, 2));
    if (tolDiff > 0.01)
    {
        tolDiff = 0.01;
    }

    // n x n matrix sizes
    // The reference solution using the Leibniz rule is very slow for large matrices since it scales with O(n!)
    const std::vector<size_t> n = {1, 2, 3, 4, 5, 6, 7, 8};

    const size_t              nruns = 1;

    // testing the vector matrix product cVec = vVec aMat
    for (size_t i = 0; i < n.size(); ++i)
    {
        // single precision doesn't suffice for inverting larger matrices, for less well behaved matrices, even n = 8 failed sometimes
        if (std::is_same<float_type, float>::value && n[i] > 40)
        {
            continue;
        }
        this->debug << "\nTesting matrix inversion via LUP decomposition " << " (float type " << typeid(float_type).name() << ") of "
        << n[i] << " by " << n[i] << " matrix A" << endl;

        // the matrix to be decomposed
        array_type aMat( 0,  n[i] - 1, 0, n[i] - 1, (float_type)0);
        float_type detAScalar = 0;
        float_type detASimd   = 0;
        float_type detAExpect = 0;

        // Make sure that the result fits by filling the array with non-uniform values.
        // The series used in other tests lead to matrices whose inversion fails in single precision.
        for (size_t ii = 0; ii < n[i]; ++ii)
        {
            for (size_t jj = 0; jj < n[i]; ++jj)
            {
                aMat(ii, jj) = static_cast<float_type>(std::cos(static_cast<double>(ii * jj * jj - ii * ii * jj + jj))
                                                       + static_cast<double>(ii)/static_cast<double>(n[i]));
            }
        }

        this->debug << "Test input data:" << endl;
        this->debug << "A   (address " << &aMat  << ") =\n" << aMat  << endl;

        bool successRef = false;
#ifdef DO_TIMING
        auto start1 = std::chrono::high_resolution_clock::now();
#endif
        for (size_t r = 0; r < nruns; ++r)
        {
            successRef = gmx::matrixDeterminantLUPRef(aMat, detAScalar);
        }
#ifdef DO_TIMING
        auto stop1 = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of the scalar        matrix determinant computation: "
        << std::chrono::duration<double, std::milli>(stop1 - start1).count() << " ms"
        << endl;
#endif

        bool successSimd = false;
#ifdef DO_TIMING
        auto start2 = std::chrono::high_resolution_clock::now();
#endif
        for (size_t r = 0; r < nruns; ++r)
        {
            successSimd = gmx::matrixDeterminantLUPSimd(aMat, detASimd);
        }
#ifdef DO_TIMING
        auto stop2 = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of the SIMD          matrix determinant computation: "
        << std::chrono::duration<double, std::milli>(stop2 - start2).count() << " ms"
        << endl;
#endif

#ifdef DO_TIMING
        auto start3 = std::chrono::high_resolution_clock::now();
#endif
        detAExpect = computeDeterminant(aMat);
#ifdef DO_TIMING
        auto stop3 = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of the Leibniz rule matrix determinant computation (single run only): "
        << std::chrono::duration<double, std::milli>(stop3 - start3).count() << " ms"
        << endl;
#endif

        this->debug << "det|A| expected =\n" << detAExpect << std::endl;
        this->debug << "det|A| (Scalar) =\n" << detAScalar << std::endl;
        this->debug << "det|A| (SIMD)   =\n" << detASimd   << std::endl;

        EXPECT_TRUE(successRef) << "Scalar implementation of matrix determinant via LUP decomposition returned failure." << std::endl;
        EXPECT_NEAR(detAExpect, detAScalar, tolDiff)
        << "Matrix determinant (scalar) of the " << n[i] << "x" << n[i] << " input matrix A differs from the reference value." << std::endl;

        EXPECT_TRUE(successSimd) << "SIMD implementation of matrix determinant via LUP decomposition returned failure." << std::endl;
        EXPECT_NEAR(detAExpect, detASimd, tolDiff)
        << "Matrix determinant (SIMD) of the " << n[i] << "x" << n[i] << " input matrix A differs from the reference value." << std::endl;
    }
}


} // namespace

} // namespace irreg_array_math_test

} // namespace gmx
