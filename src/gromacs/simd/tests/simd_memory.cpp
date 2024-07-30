/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief Tests for gmx::ArrayRef for SIMD types.
 *
 * \author Roland Schulz<roland.schulz@intel.com>
 * \ingroup module_simd
 */
#include "gmxpre.h"

#include <cstddef>

#include <array>
#include <numeric>
#include <string>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/simd/simd.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

#include "simd.h"

namespace gmx
{
#if GMX_SIMD_HAVE_REAL

/* SimdInt32 is a strange type which would never belong in an interface,
 * because its properties are peculiar to int-to-float conversions within
 * SIMD types, so there is no need to support it as a specialization of
 * SimdArrayRef. But it is useful here for better test coverage of
 * SimdArrayRef. */

template<>
class ArrayRef<SimdInt32> : public internal::SimdArrayRef<SimdInt32>
{
    using Base = internal::SimdArrayRef<SimdInt32>;
    using Base::Base;
};
template<>
class ArrayRef<const SimdInt32> : public internal::SimdArrayRef<const SimdInt32>
{
    using Base = internal::SimdArrayRef<const SimdInt32>;
    using Base::Base;
};

namespace test
{
namespace
{

TEST(EmptyArrayRefTest, IsEmpty)
{
    ArrayRef<SimdReal> empty = ArrayRef<real>();

    EXPECT_EQ(0U, empty.size());
    EXPECT_TRUE(empty.empty());
}

#    ifdef GTEST_HAS_TYPED_TEST

/*! \brief Permit all the tests to run on all kinds of ArrayRefs
 *
 * The main objective is to verify that all the different kinds of
 * construction lead to the expected result. */
template<typename TypeParam>
class ArrayRefTest : public test::SimdTest
{
public:
    using ArrayRefType         = TypeParam;
    using PointerType          = typename ArrayRefType::pointer;
    using ValueType            = typename ArrayRefType::value_type;
    using ElementType          = std::remove_const_t<gmx::internal::SimdTraitsT<ValueType>>;
    static constexpr int width = gmx::internal::SimdTraits<ValueType>::width;

    /*! \brief Run the same tests all the time
     *
     * Note that test cases must call this->runTests(), because
     * that's how the derived-class templates that implement
     * type-parameterized tests actually work. */
    void runReadOnlyTests(PointerType a, size_t aSize, ArrayRefType& arrayRef)
    {
        ASSERT_EQ(aSize, arrayRef.size());
        ASSERT_FALSE(arrayRef.empty());

        GMX_EXPECT_SIMD_EQ(load<ValueType>(a), arrayRef.front());
        GMX_EXPECT_SIMD_EQ(load<ValueType>(a + (aSize - 1) * width), arrayRef.back());

        auto it = arrayRef.begin();
        for (size_t i = 0; i != aSize; ++i, ++it)
        {
            GMX_EXPECT_SIMD_EQ(load<ValueType>(a + i * width), arrayRef[i]);
            GMX_EXPECT_SIMD_EQ(load<ValueType>(a + i * width), *it);
        }

        EXPECT_EQ(aSize, arrayRef.end() - arrayRef.begin());
        EXPECT_EQ(arrayRef.begin() + 1, arrayRef.end() - (aSize - 1));
        EXPECT_LT(arrayRef.end() - 1, arrayRef.end());
        EXPECT_GT(arrayRef.begin() + 1, arrayRef.begin());
        EXPECT_LE(arrayRef.end() - 1, arrayRef.end());
        EXPECT_GE(arrayRef.begin() + 1, arrayRef.begin());

        it = arrayRef.begin();
        it++;
        ASSERT_EQ(arrayRef.begin() + 1, it);
        it--;
        ASSERT_EQ(arrayRef.begin(), it);
        it += 1;
        ASSERT_EQ(arrayRef.begin() + 1, it);
        it -= 1;
        ASSERT_EQ(arrayRef.begin(), it);
        ++it;
        ASSERT_EQ(arrayRef.begin() + 1, it);
        --it;
        ASSERT_EQ(arrayRef.begin(), it);
    }
};

using ArrayRefTypes =
        ::testing::Types<ArrayRef<SimdReal>, ArrayRef<const SimdReal>, ArrayRef<SimdInt32>, ArrayRef<const SimdInt32>>;
TYPED_TEST_SUITE(ArrayRefTest, ArrayRefTypes);

TYPED_TEST(ArrayRefTest, ConstructFromPointersWorks)
{
    alignas(TestFixture::width * sizeof(typename TestFixture::ElementType))
            std::array<typename TestFixture::ElementType, TestFixture::width * 3>
                    a;

    std::iota(a.begin(), a.end(), 0);
    typename TestFixture::ArrayRefType arrayRef(a.data(), a.data() + a.size());
    this->runReadOnlyTests(a.data(), 3, arrayRef);
}

TYPED_TEST(ArrayRefTest, ConstructFromArrayRefWorks)
{
    alignas(TestFixture::width * sizeof(typename TestFixture::ElementType))
            std::array<typename TestFixture::ElementType, TestFixture::width * 3>
                    a;

    std::iota(a.begin(), a.end(), 0);
    ArrayRef<std::remove_const_t<typename TestFixture::ValueType>> ref(a.data(), a.data() + a.size());
    typename TestFixture::ArrayRefType                             arrayRef(ref);
    this->runReadOnlyTests(a.data(), 3, arrayRef);
}

TYPED_TEST(ArrayRefTest, ConstructFromArrayWorks)
{
    alignas(TestFixture::width * sizeof(typename TestFixture::ElementType))
            std::array<typename TestFixture::ElementType, TestFixture::width * 3>
                    a;

    std::iota(a.begin(), a.end(), 0);
    typename TestFixture::ArrayRefType arrayRef(a);
    this->runReadOnlyTests(a.data(), 3, arrayRef);
}

template<typename TypeParam>
using ArrayRefReadWriteTest = ArrayRefTest<TypeParam>;

using ArrayRefReadWriteTypes = ::testing::Types<ArrayRef<SimdReal>, ArrayRef<SimdInt32>>;
TYPED_TEST_SUITE(ArrayRefReadWriteTest, ArrayRefReadWriteTypes);

TYPED_TEST(ArrayRefReadWriteTest, Assignment)
{
    constexpr int width = TestFixture::width;

    alignas(width * sizeof(typename TestFixture::ElementType)) std::array<typename TestFixture::ElementType, width * 4> a;

    typename TestFixture::ArrayRefType arrayRef(a.data(), a.data() + a.size());

    arrayRef.front() = 1;
    EXPECT_EQ(1, a[0 * width]);
    (arrayRef.front() = 1) = 2;
    EXPECT_EQ(2, a[0 * width]);

    arrayRef[1] = 2;
    EXPECT_EQ(2, a[1 * width]);
    *(arrayRef.begin() + 2) = 3;
    EXPECT_EQ(3, a[2 * width]);
    arrayRef.back() = 4;
    EXPECT_EQ(4, a[3 * width]);
}

template<typename TypeParam>
using ArrayRefArithmeticTest = ArrayRefTest<TypeParam>;

using ArrayRefArithmeticTypes = ::testing::Types<ArrayRef<SimdReal>
#        if GMX_SIMD_HAVE_INT32_ARITHMETICS
                                                 ,
                                                 ArrayRef<SimdInt32>
#        endif
                                                 >;
TYPED_TEST_SUITE(ArrayRefArithmeticTest, ArrayRefArithmeticTypes);

TYPED_TEST(ArrayRefArithmeticTest, Basic)
{
    constexpr int width = TestFixture::width;

    alignas(width * sizeof(typename TestFixture::ElementType)) std::array<typename TestFixture::ElementType, width> a;

    typename TestFixture::ArrayRefType arrayRef(a.data(), a.data() + a.size());

    arrayRef.front() = 1;
    ASSERT_EQ(1, a[0 * width]);
    arrayRef.front() += 1;
    ASSERT_EQ(2, a[0 * width]);
    arrayRef.front() *= 2;
    ASSERT_EQ(4, a[0 * width]);
    arrayRef.front() -= 3;
    ASSERT_EQ(1, a[0 * width]);
}

#    endif // GTEST_HAS_TYPED_TEST

} // namespace
} // namespace test

#endif // GMX_HAVE_SIMD_REAL

} // namespace gmx
