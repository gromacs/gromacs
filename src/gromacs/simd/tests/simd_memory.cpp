/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017, by the GROMACS development team, led by
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
 * \brief Tests for gmx::ArrayRef for SIMD types.
 *
 * \author Roland Schulz<roland.schulz@intel.com>
 * \ingroup module_simd
 */
#include "gmxpre.h"

#include <numeric>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/simd/simd.h"

#include "simd.h"

namespace gmx
{

namespace
{

#if GMX_SIMD_HAVE_REAL

TEST(EmptyArrayRefTest, IsEmpty)
{
    EmptyArrayRef      emptyArray = {};
    ArrayRef<SimdReal> empty(emptyArray);

    EXPECT_EQ(0U, empty.size());
    EXPECT_TRUE(empty.empty());
}

#ifdef GTEST_HAS_TYPED_TEST

/*! \brief Permit all the tests to run on all kinds of ArrayRefs
 *
 * The main objective is to verify that all the different kinds of
 * construction lead to the expected result. */
template <typename TypeParam>
class ArrayRefTest : public test::SimdTest
{
    public:
        using ArrayRefType = TypeParam;
        using PointerType  = typename ArrayRefType::pointer;
        using ValueType    = typename ArrayRefType::value_type;
        using ElementType  = typename std::remove_const<typename SimdTraits<ValueType>::type>::type;
        static constexpr int width = SimdTraits<ValueType>::width;

        /*! \brief Run the same tests all the time
         *
         * Note that test cases must call this->runTests(), because
         * that's how the derived-class templates that implement
         * type-parameterized tests actually work. */
        void runReadOnlyTests(PointerType   a,
                              size_t        aSize,
                              ArrayRefType &arrayRef)
        {
            ASSERT_EQ(aSize, arrayRef.size());
            ASSERT_FALSE(arrayRef.empty());

            GMX_EXPECT_SIMD_EQ(load<ValueType>(a), arrayRef.front());
            GMX_EXPECT_SIMD_EQ(load<ValueType>(a+(aSize-1)*width), arrayRef.back());

            auto it = arrayRef.begin();
            for (size_t i = 0; i != aSize; ++i, ++it)
            {
                GMX_EXPECT_SIMD_EQ(load<ValueType>(a+i*width), arrayRef[i]);
                GMX_EXPECT_SIMD_EQ(load<ValueType>(a+i*width), *it);
            }

            EXPECT_EQ(aSize, arrayRef.end()-arrayRef.begin());
            EXPECT_EQ(arrayRef.begin()+1, arrayRef.end()-(aSize-1));
            EXPECT_LT(arrayRef.end()-1, arrayRef.end());
            EXPECT_GT(arrayRef.begin()+1, arrayRef.begin());
            EXPECT_LE(arrayRef.end()-1, arrayRef.end());
            EXPECT_GE(arrayRef.begin()+1, arrayRef.begin());

            it = arrayRef.begin();
            it++;
            ASSERT_EQ(arrayRef.begin()+1, it);
            it--;
            ASSERT_EQ(arrayRef.begin(), it);
            it += 1;
            ASSERT_EQ(arrayRef.begin()+1, it);
            it -= 1;
            ASSERT_EQ(arrayRef.begin(), it);
            ++it;
            ASSERT_EQ(arrayRef.begin()+1, it);
            --it;
            ASSERT_EQ(arrayRef.begin(), it);
        }
};

using ArrayRefTypes = ::testing::Types<ArrayRef<SimdReal>, ArrayRef<const SimdReal>,
                                       ArrayRef<SimdInt32>, ArrayRef<const SimdInt32> >;
TYPED_TEST_CASE(ArrayRefTest, ArrayRefTypes);

TYPED_TEST(ArrayRefTest, ConstructFromPointersWorks)
{
    alignas(TestFixture::width*sizeof(typename TestFixture::ElementType))
    std::array<typename TestFixture::ElementType, TestFixture::width*3> a;

    std::iota(a.begin(), a.end(), 0);
    typename TestFixture::ArrayRefType arrayRef(a.data(), a.data()+a.size());
    this->runReadOnlyTests(a.data(), 3, arrayRef);
}

TYPED_TEST(ArrayRefTest, ConstructFromArrayRefWorks)
{
    alignas(TestFixture::width*sizeof(typename TestFixture::ElementType))
    std::array<typename TestFixture::ElementType, TestFixture::width*3> a;

    std::iota(a.begin(), a.end(), 0);
    ArrayRef<typename std::remove_const<typename TestFixture::ValueType>::type>
    ref(a.data(), a.data()+a.size());
    typename TestFixture::ArrayRefType        arrayRef(ref);
    this->runReadOnlyTests(a.data(), 3, arrayRef);
}

TYPED_TEST(ArrayRefTest, ConstructFromArrayWorks)
{
    alignas(TestFixture::width*sizeof(typename TestFixture::ElementType))
    std::array<typename TestFixture::ElementType, TestFixture::width*3> a;

    std::iota(a.begin(), a.end(), 0);
    typename TestFixture::ArrayRefType arrayRef(a);
    this->runReadOnlyTests(a.data(), 3, arrayRef);
}

template <typename TypeParam>
using ArrayRefReadWriteTest = ArrayRefTest<TypeParam>;

using ArrayRefReadWriteTypes = ::testing::Types< ArrayRef<SimdReal>, ArrayRef<SimdInt32> >;
TYPED_TEST_CASE(ArrayRefReadWriteTest, ArrayRefReadWriteTypes);

TYPED_TEST(ArrayRefReadWriteTest, Assignment)
{
    constexpr int width = TestFixture::width;

    alignas(width*sizeof(typename TestFixture::ElementType))
    std::array<typename TestFixture::ElementType, width*4> a;

    typename TestFixture::ArrayRefType arrayRef(a.data(), a.data()+a.size());

    arrayRef.front() = 1;
    EXPECT_EQ(1, a[0*width]);
    (arrayRef.front() = 1) = 2;
    EXPECT_EQ(2, a[0*width]);

    arrayRef[1] = 2;
    EXPECT_EQ(2, a[1*width]);
    *(arrayRef.begin()+2) = 3;
    EXPECT_EQ(3, a[2*width]);
    arrayRef.back() = 4;
    EXPECT_EQ(4, a[3*width]);
}

template <typename TypeParam>
using ArrayRefArithmeticTest = ArrayRefTest<TypeParam>;

using ArrayRefArithmeticTypes = ::testing::Types< ArrayRef<SimdReal>
#if GMX_SIMD_HAVE_INT32_ARITHMETICS
                                                  , ArrayRef<SimdInt32>
#endif
                                                  >;
TYPED_TEST_CASE(ArrayRefArithmeticTest, ArrayRefArithmeticTypes);

TYPED_TEST(ArrayRefArithmeticTest, Basic)
{
    constexpr int width = TestFixture::width;

    alignas(width*sizeof(typename TestFixture::ElementType))
    std::array<typename TestFixture::ElementType, width> a;

    typename TestFixture::ArrayRefType arrayRef(a.data(), a.data()+a.size());

    arrayRef.front() = 1;
    ASSERT_EQ(1, a[0*width]);
    arrayRef.front() += 1;
    ASSERT_EQ(2, a[0*width]);
    arrayRef.front() *= 2;
    ASSERT_EQ(4, a[0*width]);
    arrayRef.front() -= 3;
    ASSERT_EQ(1, a[0*width]);
}

#endif // GTEST_HAS_TYPED_TEST

#endif // GMX_HAVE_SIMD_REAL

}      // namespace

}      // namespace
