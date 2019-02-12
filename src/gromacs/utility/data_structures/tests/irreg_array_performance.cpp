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
/*! \internal \file
   \brief
   Tests of different access variants of the irregular array data structures

   For development, the tests can be run with a '-stdout' command-line option
   to print out the help to stdout instead of using the XML reference
   framework.

   The tests in this file are not implemented as type parametrized test
   becausethe underlying array differs for flat and non-flat array variants.

   \author R. Thomas Ullmann <tullman@gwdg.de>

   \ingroup module_utility
 */
#include "gmxpre.h"

#include "data_struct_test_commons.h"

#define DO_TIMING
#ifdef DO_TIMING
#include <chrono>
#endif

namespace gmx
{

namespace data_struct_test
{

namespace
{

using std::endl;

class IrregArrayPerformanceTest : public gmx::data_struct_test::DataStructTest
{
};

/*********************************************************************/

// NOTE:
// precomputing the start and end indices of the array stripe in the last dimension is much faster
// but here I want to compare the relative cost of the different access variants
// results: operator[] is nearly twofold slower than operator() -- is there a better implementation?
//          the ratio between book-keeping and storage data is heavily dependent on the length of the
//           array stripes in the highest dimension,
//          For realistic examples, 3.2 bytes book-keeping data per stored array element was
//           the worst case tested (for a 1000 x 5 x 1000 x 5 IrregArray4D).

/*! \brief compute the sum of all array elements for testing access times

    \tparam   T       data type stored by the array
    \tparam   Alloc   allocator use by the array

    \param[in]    test_arr    array to be tested
    \param[in]    ncycles     number of times to sum over the array elements

    \returns      ncycles times the sum over all array elements
 */
template<typename T, class Alloc>
inline long double sumArray(const gmx::IrregArray2D<T, Alloc, false> &test_arr, const size_t ncycles)
{
    long double sum = 0;
    if (!test_arr.isIrreg())
    {
        for (size_t c = 0; c < ncycles; ++c)
        {
            for (size_t i = 0; i < test_arr.length1(); ++i)
            {
                for (size_t j = 0; j < test_arr.length2(); ++j)
                {
                    sum += test_arr.data()[i][j];
                }
            }
        }
    }
    else
    {
        for (size_t c = 0; c < ncycles; ++c)
        {
            for (size_t i = 0; i < test_arr.length1(); ++i)
            {
                for (size_t j = 0; j < test_arr.length2(i); ++j)
                {
                    sum += test_arr.data()[i][j];
                }
            }
        }
    }
    return sum;
}

/*! \brief compute the sum of all array elements for testing access times

    \tparam   T       data type stored by the array
    \tparam   Alloc   allocator use by the array

    \param[in]    test_arr    array to be tested
    \param[in]    ncycles     number of times to sum over the array elements

    \returns      cycles times the sum over all array elements
 */
template<typename T, class Alloc>
inline long double sumArray(const gmx::IrregArray2D<T, Alloc, true> &test_arr, const size_t ncycles)
{
    long double sum = 0;
    if (!test_arr.isIrreg())
    {

        for (size_t c = 0; c < ncycles; ++c)
        {
            for (size_t i = 0; i < test_arr.length1(); ++i)
            {
                for (size_t j = test_arr.arrayIndex(i, 0); j < test_arr.arrayIndex(i, 0) + test_arr.length2(); ++j)
                {
                    sum += test_arr.data()[j];
                }
            }
        }
    }
    else
    {
        for (size_t c = 0; c < ncycles; ++c)
        {
            for (size_t i = 0; i < test_arr.length1(); ++i)
            {
                for (size_t j = test_arr.arrayIndex(i, 0); j < test_arr.arrayIndex(i, 0) + test_arr.length2(i); ++j)
                {
                    sum += test_arr.data()[j];
                }
            }
        }
    }
    return sum;
}

TYPED_TEST(ParamIrregArray2DTest, 2DArrayAccessOperatorPerformance)
{
    typedef TypeParam                                        array_type;
    typedef typename array_type::size_type                    size_type;
    typedef typename array_type::size_1d_array_type  size_1d_array_type;

    constexpr size_type      ncycles = 1;
    constexpr size_type      length1 = 100;
    constexpr size_type      l2      = 10;
    const size_1d_array_type length2(length1, l2);
    array_type               test_arr(length1, l2, 1);
    array_type               test_arr2(length2, 1);

    {
#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        const long double sum = sumArray(test_arr, ncycles);

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of bare array access test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        this->debug << "sum1 " << sum << endl;
    }

    {
        long double sum = 0;

#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        for (size_type c = 0; c < ncycles; ++c)
        {
            for (size_type i = 0; i < test_arr.length1(); ++i)
            {
                for (size_type j = 0; j < test_arr.length2(); ++j)
                {
                    sum += test_arr(i, j);
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of operator() test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        this->debug << "sum1 " << sum << endl;
    }

    {
        long double sum = 0;

#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        for (size_type c = 0; c < ncycles; ++c)
        {
            for (size_type i = 0; i < test_arr.length1(); ++i)
            {
                for (size_type j = 0; j < test_arr.length2(); ++j)
                {
                    sum += test_arr[i][j];
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of operator[] test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        this->debug << "sum2 " << sum << endl;
    }

    {
#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        const long double sum = sumArray(test_arr2, ncycles);

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of bare array access test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        this->debug << "sum1 " << sum << endl;
    }

    {
        long double sum = 0;

#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        for (size_type c = 0; c < ncycles; ++c)
        {
            for (size_type i = 0; i < test_arr2.length1(); ++i)
            {
                for (size_type j = 0; j < test_arr2.length2(i); ++j)
                {
                    sum += test_arr2(i, j);
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of operator() test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        this->debug << "sum1 " << sum << endl;
    }

    {
        long double sum = 0;

#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        for (size_type c = 0; c < ncycles; ++c)
        {
            for (size_type i = 0; i < test_arr2.length1(); ++i)
            {
                for (size_type j = 0; j < test_arr2.length2(i); ++j)
                {
                    sum += test_arr2[i][j];
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of operator[] test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        this->debug << "sum2 " << sum << endl;
    }

}


/*! \brief compute the sum of all array elements for testing access times

    \tparam   T       data type stored by the array
    \tparam   Alloc   allocator use by the array

    \param[in]    test_arr    array to be tested
    \param[in]    ncycles     number of times to sum over the array elements

    \returns      ncycles times the sum over all array elements
 */
template<typename T, class Alloc>
inline long double sumArray(const gmx::IrregArray3D<T, Alloc, false> &test_arr, const size_t ncycles)
{
    long double sum = 0;
    if (!test_arr.isIrreg())
    {
        for (size_t c = 0; c < ncycles; ++c)
        {
            for (size_t i = 0; i < test_arr.length1(); ++i)
            {
                for (size_t j = 0; j < test_arr.length2(); ++j)
                {
                    for (size_t k = 0; k < test_arr.length3(); ++k)
                    {
                        sum += test_arr.data()[i][j][k];
                    }
                }
            }
        }
    }
    else
    {
        for (size_t c = 0; c < ncycles; ++c)
        {
            for (size_t i = 0; i < test_arr.length1(); ++i)
            {
                for (size_t j = 0; j < test_arr.length2(i); ++j)
                {
                    for (size_t k = 0; k < test_arr.length3(i, j); ++k)
                    {
                        sum += test_arr.data()[i][j][k];
                    }
                }
            }
        }
    }
    return sum;
}

/*! \brief compute the sum of all array elements for testing access times

    \tparam   T       data type stored by the array
    \tparam   Alloc   allocator use by the array

    \param[in]    test_arr    array to be tested
    \param[in]    ncycles     number of times to sum over the array elements

    \returns      ncycles times the sum over all array elements
 */
template<typename T, class Alloc>
inline long double sumArray(const gmx::IrregArray3D<T, Alloc, true> &test_arr, const size_t ncycles)
{
    long double sum = 0;
    if (!test_arr.isIrreg())
    {
        for (size_t c = 0; c < ncycles; ++c)
        {
            for (size_t i = 0; i < test_arr.length1(); ++i)
            {
                for (size_t j = 0; j < test_arr.length2(); ++j)
                {
                    for (size_t k = test_arr.arrayIndex(i, j, 0); k < test_arr.arrayIndex(i, j, 0) + test_arr.length3(); ++k)
                    {
                        sum += test_arr.data()[k];
                    }
                }
            }
        }
    }
    else
    {
        for (size_t c = 0; c < ncycles; ++c)
        {
            for (size_t i = 0; i < test_arr.length1(); ++i)
            {
                for (size_t j = 0; j < test_arr.length2(i); ++j)
                {
                    for (size_t k = test_arr.arrayIndex(i, j, 0); k < test_arr.arrayIndex(i, j, 0) + test_arr.length3(i, j); ++k)
                    {
                        sum += test_arr.data()[k];
                    }
                }
            }
        }
    }
    return sum;
}

TYPED_TEST(ParamIrregArray3DTest, 3DArrayAccessOperatorPerformance)
{
// The Clang static analyzer doesn't track the value of the isIrreg_ member variable
// resulting in a false positive warning for a nullptr dereference.
// https://bugs.llvm.org/show_bug.cgi?id=39064
#ifndef __clang_analyzer__
    typedef TypeParam                                        array_type;
    typedef typename array_type::size_type                    size_type;
    typedef typename array_type::size_1d_array_type  size_1d_array_type;
    typedef typename array_type::size_2d_array_type  size_2d_array_type;
    const size_t                    ncycles = 1;
    constexpr size_type             l1      = 30;
    constexpr size_type             l2      = 5;
    const size_t                    l3      = 30;
    array_type test_arr(l1, l2, l3, 1);
    const size_1d_array_type        length2(l1, l2);
    const size_2d_array_type        length3(length2, l3);
    array_type test_arr2(length3, 1);

    {
#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        const long double sum = sumArray(test_arr, ncycles);

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of bare array access test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        this->debug << "sum " << sum << endl;
    }

    {
        long double sum = 0;

#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        for (size_t c = 0; c < ncycles; ++c)
        {
            for (size_type i = 0; i < test_arr.length1(); ++i)
            {
                for (size_type j = 0; j < test_arr.length2(); ++j)
                {
                    for (size_type k = 0; k < test_arr.length3(); ++k)
                    {
                        sum += test_arr(i, j, k);
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of operator() test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        this->debug << "sum " << sum << endl;
    }

    {
        long double sum = 0;

#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        for (size_type c = 0; c < ncycles; ++c)
        {
            for (size_type i = 0; i < test_arr.length1(); ++i)
            {
                for (size_type j = 0; j < test_arr.length2(); ++j)
                {
                    for (size_type k = 0; k < test_arr.length3(); ++k)
                    {
                        sum += test_arr[i][j][k];
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of operator[] test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        this->debug << "sum " << sum << endl;
    }

    {
#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        const long double sum = sumArray(test_arr2, ncycles);

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of bare array access test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        this->debug << "sum1 " << sum << endl;
    }

    {
        long double sum = 0;

#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        for (size_type c = 0; c < ncycles; ++c)
        {
            for (size_type i = 0; i < test_arr2.length1(); ++i)
            {
                for (size_type j = 0; j < test_arr2.length2(i); ++j)
                {
                    for (size_type k = 0; k < test_arr2.length3(i, j); ++k)
                    {
                        sum += test_arr2(i, j, k);
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of operator() test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        this->debug << "sum1 " << sum << endl;
    }

    {
        long double sum = 0;

#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        for (size_type c = 0; c < ncycles; ++c)
        {
            for (size_type i = 0; i < test_arr2.length1(); ++i)
            {
                for (size_type j = 0; j < test_arr2.length2(i); ++j)
                {
                    for (size_type k = 0; k < test_arr2.length3(i, j); ++k)
                    {
                        sum += test_arr2[i][j][k];
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of operator[] test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        this->debug << "sum2 " << sum << endl;
    }

#endif          // end __clang_analyzer__
}

/*! \brief compute the sum of all array elements for testing access times

    \tparam   T       data type stored by the array
    \tparam   Alloc   allocator use by the array

    \param[in]    test_arr    array to be tested
    \param[in]    ncycles     number of times to sum over the array elements

    \returns      ncycles times the sum over all array elements
 */
template<typename T, class Alloc>
inline long double sumArray(const gmx::IrregArray4D<T, Alloc, false> &test_arr, const size_t ncycles)
{
    long double sum = 0;
    if (!test_arr.isIrreg())
    {
        for (size_t c = 0; c < ncycles; ++c)
        {
            for (size_t i = 0; i < test_arr.length1(); ++i)
            {
                for (size_t j = 0; j < test_arr.length2(); ++j)
                {
                    for (size_t k = 0; k < test_arr.length3(); ++k)
                    {
                        for (size_t l = 0; l < test_arr.length4(); ++l)
                        {
                            sum += test_arr.data()[i][j][k][l];
                        }
                    }
                }
            }
        }
    }
    else
    {
        for (size_t c = 0; c < ncycles; ++c)
        {
            for (size_t i = 0; i < test_arr.length1(); ++i)
            {
                for (size_t j = 0; j < test_arr.length2(i); ++j)
                {
                    for (size_t k = 0; k < test_arr.length3(i, j); ++k)
                    {
                        for (size_t l = 0; l < test_arr.length4(i, j, k); ++l)
                        {
                            sum += test_arr.data()[i][j][k][l];
                        }
                    }
                }
            }
        }
    }
    return sum;
}

/*! \brief compute the sum of all array elements for testing access times

    \tparam   T       data type stored by the array
    \tparam   Alloc   allocator use by the array

    \param[in]    test_arr    array to be tested
    \param[in]    ncycles     number of times to sum over the array elements

    \returns      ncycles times the sum over all array elements
 */
template<typename T, class Alloc>
inline long double sumArray(const gmx::IrregArray4D<T, Alloc, true> &test_arr, const size_t ncycles)
{
    long double sum = 0;
    if (!test_arr.isIrreg())
    {
        for (size_t c = 0; c < ncycles; ++c)
        {
            for (size_t i = 0; i < test_arr.length1(); ++i)
            {
                for (size_t j = 0; j < test_arr.length2(); ++j)
                {
                    for (size_t k = 0; k < test_arr.length3(); ++k)
                    {
                        for (size_t l = test_arr.arrayIndex(i, j, k, 0); l < test_arr.arrayIndex(i, j, k, 0) + test_arr.length4(); ++l)
                        {
                            sum += test_arr.data()[l];
                        }
                    }
                }
            }
        }
    }
    else
    {
        for (size_t c = 0; c < ncycles; ++c)
        {
            for (size_t i = 0; i < test_arr.length1(); ++i)
            {
                for (size_t j = 0; j < test_arr.length2(i); ++j)
                {
                    for (size_t k = 0; k < test_arr.length3(i, j); ++k)
                    {
                        for (size_t l = test_arr.arrayIndex(i, j, k, 0); l < test_arr.arrayIndex(i, j, k, 0) + test_arr.length4(i, j, k); ++l)
                        {
                            sum += test_arr.data()[l];
                        }
                    }
                }
            }
        }
    }
    return sum;
}

/*! \brief compute the sum of all array elements for testing access times

    \tparam   T       data type stored by the array
    \tparam   Alloc   allocator use by the array

    \param[in]    test_arr    array to be tested
    \param[in]    ncycles     number of times to sum over the array elements

    \returns      ncycles times the sum over all array elements
 */
template<typename T, class Alloc>
inline long double sumArrayCachedIndices(const gmx::IrregArray4D<T, Alloc, false> &test_arr, const size_t ncycles)
{
    long double sum = 0;
    if (!test_arr.isIrreg())
    {
        for (size_t c = 0; c < ncycles; ++c)
        {
            const size_t l1 = test_arr.length1();
            const size_t l2 = test_arr.length2();
            const size_t l3 = test_arr.length3();
            const size_t l4 = test_arr.length4();
            for (size_t i = 0; i < l1; ++i)
            {
                for (size_t j = 0; j < l2; ++j)
                {
                    for (size_t k = 0; k < l3; ++k)
                    {
                        for (size_t l = 0; l < l4; ++l)
                        {
                            sum += test_arr.data()[i][j][k][l];
                        }
                    }
                }
            }
        }
    }
    else
    {
        for (size_t c = 0; c < ncycles; ++c)
        {
            const size_t l1 = test_arr.length1();
            for (size_t i = 0; i < l1; ++i)
            {
                const size_t l2 = test_arr.length2(i);
                for (size_t j = 0; j < l2; ++j)
                {
                    const size_t l3 = test_arr.length3(i, j);
                    for (size_t k = 0; k < l3; ++k)
                    {
                        const size_t l4  = test_arr.length4(i, j, k);
                        for (size_t l = 0; l < l4; ++l)
                        {
                            sum += test_arr.data()[i][j][k][l];
                        }
                    }
                }
            }
        }
    }
    return sum;
}

/*! \brief compute the sum of all array elements for testing access times

    \tparam   T       data type stored by the array
    \tparam   Alloc   allocator use by the array

    \param[in]    test_arr    array to be tested
    \param[in]    ncycles     number of times to sum over the array elements

    \returns      ncycles times the sum over all array elements
 */
template<typename T, class Alloc>
inline long double sumArrayCachedIndices(const gmx::IrregArray4D<T, Alloc, true> &test_arr, const size_t ncycles)
{
    long double sum = 0;
    if (!test_arr.isIrreg())
    {
        for (size_t c = 0; c < ncycles; ++c)
        {
            const size_t l1 = test_arr.length1();
            const size_t l2 = test_arr.length2();
            const size_t l3 = test_arr.length3();
            for (size_t i = 0; i < l1; ++i)
            {
                for (size_t j = 0; j < l2; ++j)
                {
                    for (size_t k = 0; k < l3; ++k)
                    {
                        const size_t l4 = test_arr.length4();
                        for (size_t l = 0; l < l4; ++l)
                        {
                            sum += test_arr.data()[l];
                        }
                    }
                }
            }
        }
    }
    else
    {
        for (size_t c = 0; c < ncycles; ++c)
        {
            const size_t l1 = test_arr.length1();
            for (size_t i = 0; i < l1; ++i)
            {
                const size_t l2 = test_arr.length2(i);
                for (size_t j = 0; j < l2; ++j)
                {
                    const size_t l3 = test_arr.length3(i, j);
                    for (size_t k = 0; k < l3; ++k)
                    {
                        const size_t l4 = test_arr.length4(i, j, k) - 1;
                        for (size_t l = 0; l < l4; ++l)
                        {
                            sum += test_arr.data()[l];
                        }
                    }
                }
            }
        }
    }
    return sum;
}

TYPED_TEST(ParamIrregArray4DTest, 4DArrayAccessOperatorPerformance)
{
    typedef TypeParam                                        array_type;
    typedef typename array_type::size_type                    size_type;
    typedef typename array_type::size_1d_array_type  size_1d_array_type;
    typedef typename array_type::size_2d_array_type  size_2d_array_type;
    typedef typename array_type::size_3d_array_type  size_3d_array_type;
    constexpr size_type      ncycles = 1;
    constexpr size_type      len1    = 30;
    constexpr size_type      len2    = 5;
    constexpr size_type      len3    = 30;
    constexpr size_type      len4    = 5;
    const size_1d_array_type length2(len1, len2);
    const size_2d_array_type length3(length2, len3);
    const size_3d_array_type length4(length3, len4);
    array_type               test_arr(len1, len2, len3, len4, 1);
    array_type               test_arr2(length4, 1);

    {
#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        const long double sum = sumArray(test_arr, ncycles);

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of bare array access test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        this->debug << "sum1 " << sum << endl;
    }

    {
        long double sum = 0;

#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        for (size_type c = 0; c < ncycles; ++c)
        {
            for (size_type i = 0; i < test_arr.length1(); ++i)
            {
                for (size_type j = 0; j < test_arr.length2(); ++j)
                {
                    for (size_type k = 0; k < test_arr.length3(); ++k)
                    {
                        for (size_type l = 0; l < test_arr.length4(); ++l)
                        {
                            sum += test_arr(i, j, k, l);
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of operator() test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        this->debug << "sum1 " << sum << endl;
    }

    {
        long double sum = 0;

#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        for (size_type c = 0; c < ncycles; ++c)
        {
            for (size_type i = 0; i < test_arr.length1(); ++i)
            {
                for (size_type j = 0; j < test_arr.length2(); ++j)
                {
                    for (size_type k = 0; k < test_arr.length3(); ++k)
                    {
                        for (size_type l = 0; l < test_arr.length4(); ++l)
                        {
                            sum += test_arr[i][j][k][l];
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of operator[] test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        this->debug << "sum2 " << sum << endl;
    }

    {
#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        const long double sum = sumArray(test_arr2, ncycles);

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of bare array access test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        this->debug << "sum1 " << sum << endl;
    }

    {
        long double sum = 0;

#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        for (size_type c = 0; c < ncycles; ++c)
        {
            for (size_type i = 0; i < test_arr2.length1(); ++i)
            {
                for (size_type j = 0; j < test_arr2.length2(i); ++j)
                {
                    for (size_type k = 0; k < test_arr2.length3(i, j); ++k)
                    {
                        for (size_type l = 0; l < test_arr2.length4(i, j, k); ++l)
                        {
                            sum += test_arr2(i, j, k, l);
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of operator() test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        this->debug << "sum1 " << sum << endl;
    }

    {
        long double sum = 0;

#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        for (size_type c = 0; c < ncycles; ++c)
        {
            for (size_type i = 0; i < test_arr2.length1(); ++i)
            {
                for (size_type j = 0; j < test_arr2.length2(i); ++j)
                {
                    for (size_type k = 0; k < test_arr2.length3(i, j); ++k)
                    {
                        for (size_type l = 0; l < test_arr2.length4(i, j, k); ++l)
                        {
                            sum += test_arr2[i][j][k][l];
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of operator[] test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        this->debug << "sum2 " << sum << endl;
    }
}


TYPED_TEST(ParamIrregArray4DTest, 4DArrayAccessOperatorPerformanceCachedIndices)
{
    typedef TypeParam                                        array_type;
    typedef typename array_type::size_type                    size_type;
    typedef typename array_type::index_type                  index_type;
    typedef typename array_type::size_1d_array_type  size_1d_array_type;
    typedef typename array_type::size_2d_array_type  size_2d_array_type;
    typedef typename array_type::size_3d_array_type  size_3d_array_type;
    const size_type          ncycles = 1;
    const size_type          len1    = 30;
    const size_type          len2    = 5;
    const size_type          len3    = 30;
    const size_type          len4    = 5;
    const size_1d_array_type length2(len1, len2);
    const size_2d_array_type length3(length2, len3);
    const size_3d_array_type length4(length3, len4);
    array_type               test_arr(len1, len2, len3, len4, 1);
    array_type               test_arr2(length4, 1);

    {
#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        const long double sum = sumArrayCachedIndices(test_arr, ncycles);

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of bare array access test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        this->debug << "sum1 " << sum << endl;
    }

    {
        long double sum = 0;

#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        for (size_type c = 0; c < ncycles; ++c)
        {
            const index_type l1 = test_arr.length1();
            const index_type l2 = test_arr.length2();
            const index_type l3 = test_arr.length3();
            const index_type l4 = test_arr.length4();
            for (index_type i = 0; i < l1; ++i)
            {
                for (index_type j = 0; j < l2; ++j)
                {
                    for (index_type k = 0; k < l3; ++k)
                    {
                        for (ssize_t l = 0; l < l4; ++l)
                        {
                            sum += test_arr(i, j, k, l);
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of operator() test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        this->debug << "sum1 " << sum << endl;
    }

    {
        long double sum = 0;

#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        for (size_type c = 0; c < ncycles; ++c)
        {
            const index_type l1 = test_arr.length1();
            const index_type l2 = test_arr.length2();
            const index_type l3 = test_arr.length3();
            const index_type l4 = test_arr.length4();
            for (index_type i = 0; i < l1; ++i)
            {
                for (index_type j = 0; j < l2; ++j)
                {
                    for (index_type k = 0; k < l3; ++k)
                    {
                        for (ssize_t l = 0; l < l4; ++l)
                        {
                            sum += test_arr[i][j][k][l];
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of operator[] test for the regular special case of an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        this->debug << "sum2 " << sum << endl;
    }

    {
#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        const long double sum = sumArrayCachedIndices(test_arr2, ncycles);

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of bare array access test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        this->debug << "sum1 " << sum << endl;
    }

    {
        long double sum = 0;

#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        for (size_type c = 0; c < ncycles; ++c)
        {
            const index_type l1  = test_arr2.length1();
            for (index_type i = 0; i < l1; ++i)
            {
                const index_type l2  = test_arr2.length2(i);
                for (index_type j = 0; j < l2; ++j)
                {
                    const index_type l3  = test_arr2.length3(i, j);
                    for (index_type k = 0; k < l3; ++k)
                    {
                        const index_type l4  = test_arr2.length4(i, j, k);
                        for (index_type l = 0; l < l4; ++l)
                        {
                            sum += test_arr2(i, j, k, l);
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of operator() test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        this->debug << "sum1 " << sum << endl;
    }

    {
        long double sum = 0;

#ifdef DO_TIMING
        // compare performance of the access operator variants
        auto start = std::chrono::high_resolution_clock::now();
#endif

        for (size_type c = 0; c < ncycles; ++c)
        {
            const index_type l1  = test_arr2.length1();
            for (index_type i = 0; i < l1; ++i)
            {
                const index_type l2  = test_arr2.length2(i);
                for (index_type j = 0; j < l2; ++j)
                {
                    const index_type l3  = test_arr2.length3(i, j);
                    for (index_type k = 0; k < l3; ++k)
                    {
                        const index_type l4  = test_arr2.length4(i, j, k);
                        for (index_type l = 0; l < l4; ++l)
                        {
                            sum += test_arr2[i][j][k][l];
                        }
                    }
                }
            }
        }

#ifdef DO_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        this->debug << "Runtime of operator[] test for an irregular array according to chrono: "
        << std::chrono::duration<double, std::milli>(stop - start).count() << " ms"
        << endl;
#endif
        this->debug << "sum2 " << sum << endl;
    }

}

} // namespace

} // namespace data_struct_test

} // namespace gmx

#ifdef DO_TIMING
#undef DO_TIMING
#endif
