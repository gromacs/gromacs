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
 * \brief Tests for gmx::ArrayRef.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/arrayref.h"

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

namespace gmx
{

namespace
{

TEST(EmptyArrayRefTest, IsEmpty)
{
    EmptyArrayRef  emptyArray = {};
    ArrayRef<real> empty(emptyArray);

    EXPECT_EQ(0U, empty.size());
    EXPECT_TRUE(empty.empty());
}

TEST(EmptyConstArrayRefTest, IsEmpty)
{
    EmptyArrayRef        emptyArray = {};
    ArrayRef<const real> empty(emptyArray);

    EXPECT_EQ(0U, empty.size());
    EXPECT_TRUE(empty.empty());
}

#ifdef GTEST_HAS_TYPED_TEST

//! Define the types that end up being available as TypeParam in the test cases for both kinds of ArrayRef
typedef ::testing::Types<
        ArrayRef<char>,
        ArrayRef<unsigned char>,
        ArrayRef<int>,
        ArrayRef<unsigned int>,
        ArrayRef<long>,
        ArrayRef<unsigned long>,
        ArrayRef<gmx_int64_t>,
        ArrayRef<gmx_uint64_t>,
        ArrayRef<float>,
        ArrayRef<double>,
        ArrayRef<const char>,
        ArrayRef<const unsigned char>,
        ArrayRef<const int>,
        ArrayRef<const unsigned int>,
        ArrayRef<const long>,
        ArrayRef<const unsigned long>,
        ArrayRef<const gmx_int64_t>,
        ArrayRef<const gmx_uint64_t>,
        ArrayRef<const float>,
        ArrayRef<const double>
        > ArrayRefTypes;

/*! \brief Permit all the tests to run on all kinds of ArrayRefs
 *
 * The main objective is to verify that all the different kinds of
 * construction lead to the expected result. */
template <typename TypeParam>
class ArrayRefTest : public ::testing::Test
{
    public:
        typedef TypeParam ArrayRefType;
        typedef typename ArrayRefType::value_type ValueType;
        typedef typename std::remove_const<ValueType>::type NonConstValueType;

        /*! \brief Run the same tests all the time
         *
         * Note that test cases must call this->runTests(), because
         * that's how the derived-class templates that implement
         * type-parameterized tests actually work. */
        void runTests(ValueType     *a,
                      size_t         aSize,
                      ValueType     *aData,
                      ArrayRefType  &arrayRef)
        {
            ASSERT_EQ(aSize, arrayRef.size());
            ASSERT_FALSE(arrayRef.empty());
            EXPECT_EQ(aData, arrayRef.data());
            EXPECT_EQ(a[0], arrayRef.front());
            EXPECT_EQ(a[aSize-1], arrayRef.back());
            for (size_t i = 0; i != aSize; ++i)
            {
                EXPECT_EQ(a[i], arrayRef[i]);
            }
        }
};

TYPED_TEST_CASE(ArrayRefTest, ArrayRefTypes);

/* Welcome back to the past. While you can declare a static array[] of
   templated type in a class, in C++98, you have to define it outside
   the class, and when you do, the compiler knows the declaration is
   incomplete and can't match the types to actual functions. So,
   declaring locals is the only choice available, so we need macros to
   avoid duplication. Lovely. */
#define DEFINE_ARRAY(a, aSize)                                  \
    typename TestFixture::ValueType (a)[] = {                   \
        static_cast<typename TestFixture::ValueType>(1.2),      \
        static_cast<typename TestFixture::ValueType>(2.4),      \
        static_cast<typename TestFixture::ValueType>(3.1)       \
    };                                                          \
    size_t (aSize) = sizeof((a)) / sizeof(typename TestFixture::ValueType);

TYPED_TEST(ArrayRefTest, MakeWithAssignmentWorks)
{
    DEFINE_ARRAY(a, aSize);
    typename TestFixture::ArrayRefType arrayRef = a;
    this->runTests(a, aSize, a, arrayRef);
}

TYPED_TEST(ArrayRefTest, ConstructWithTemplateConstructorWorks)
{
    DEFINE_ARRAY(a, aSize);
    typename TestFixture::ArrayRefType arrayRef(a);
    this->runTests(a, aSize, a, arrayRef);
}

TYPED_TEST(ArrayRefTest, ConstructFromPointersWorks)
{
    DEFINE_ARRAY(a, aSize);
    typename TestFixture::ArrayRefType arrayRef(a, a + aSize);
    this->runTests(a, aSize, a, arrayRef);
}

TYPED_TEST(ArrayRefTest, MakeFromPointersWorks)
{
    DEFINE_ARRAY(a, aSize);
    typename TestFixture::ArrayRefType arrayRef
        = TestFixture::ArrayRefType::fromPointers(a, a + aSize);
    this->runTests(a, aSize, a, arrayRef);
}

TYPED_TEST(ArrayRefTest, MakeFromArrayWorks)
{
    DEFINE_ARRAY(a, aSize);
    typename TestFixture::ArrayRefType arrayRef
        = TestFixture::ArrayRefType::fromArray(a, aSize);
    this->runTests(a, aSize, a, arrayRef);
}

TYPED_TEST(ArrayRefTest, ConstructFromVectorWorks)
{
    DEFINE_ARRAY(a, aSize);
    std::vector<typename TestFixture::NonConstValueType> v(a, a + aSize);
    typename TestFixture::ArrayRefType                   arrayRef(v);
    this->runTests(a, v.size(), v.data(), arrayRef);
}

TYPED_TEST(ArrayRefTest, MakeFromVectorWorks)
{
    DEFINE_ARRAY(a, aSize);
    std::vector<typename TestFixture::NonConstValueType> v(a, a + aSize);
    typename TestFixture::ArrayRefType arrayRef
        = TestFixture::ArrayRefType::fromVector(v.begin(), v.end());
    this->runTests(a, v.size(), v.data(), arrayRef);
}

//! Helper struct for the case actually used in mdrun signalling
template <typename T>
struct Helper
{
    public:
        T   a[3];
        int size;
};

/*! \brief Test of the case actually used in mdrun signalling
 *
 * There, we take a non-const struct-field array of static length and
 * make an ArrayRef to it using the template constructor that is
 * supposed to infer the length from the static size. This has
 * been a problem (for a compiler that we no longer support),
 * so we test it.
 */

TYPED_TEST(ArrayRefTest, ConstructFromStructFieldWithTemplateConstructorWorks)
{
    DEFINE_ARRAY(a, aSize);
    Helper<typename TestFixture::NonConstValueType> h;
    h.size = aSize;
    for (int i = 0; i != h.size; ++i)
    {
        h.a[i] = a[i];
    }
    typename TestFixture::ArrayRefType arrayRef(h.a);
    this->runTests(h.a, h.size, h.a, arrayRef);
}

#else   // GTEST_HAS_TYPED_TEST

/* A dummy test that at least signals that something is missing if one runs the
 * unit test executable itself.
 */
TEST(DISABLED_ArrayRefTest, GenericTests)
{
    ADD_FAILURE()
    << "Tests for generic ArrayRef functionality require support for "
    << "Google Test typed tests, which was not available when the tests "
    << "were compiled.";
}

#endif // GTEST_HAS_TYPED_TEST

}      // namespace

}      // namespace
