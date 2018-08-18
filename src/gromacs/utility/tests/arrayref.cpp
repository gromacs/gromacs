/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017,2018, by the GROMACS development team, led by
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
        ArrayRef<int64_t>,
        ArrayRef<uint64_t>,
        ArrayRef<float>,
        ArrayRef<double>,
        ArrayRef<const char>,
        ArrayRef<const unsigned char>,
        ArrayRef<const int>,
        ArrayRef<const unsigned int>,
        ArrayRef<const long>,
        ArrayRef<const unsigned long>,
        ArrayRef<const int64_t>,
        ArrayRef<const uint64_t>,
        ArrayRef<const float>,
        ArrayRef<const double>
        > ArrayRefTypes;

constexpr index aSize = 3;

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
        void runTests(ValueType     *aData,
                      ArrayRefType  &arrayRef)
        {
            ASSERT_EQ(aSize, arrayRef.size());
            ASSERT_FALSE(arrayRef.empty());
            EXPECT_EQ(aData, arrayRef.data());
            EXPECT_EQ(a[0], arrayRef.front());
            EXPECT_EQ(a[aSize-1], arrayRef.back());
            for (index i = 0; i != aSize; ++i)
            {
                EXPECT_EQ(a[i], arrayRef[i]);
            }
        }

        ValueType         a[aSize]  = { ValueType(1.2), ValueType(2.4), ValueType(3.1) };
        NonConstValueType ma[aSize] = { ValueType(1.2), ValueType(2.4), ValueType(3.1) };
};

TYPED_TEST_CASE(ArrayRefTest, ArrayRefTypes);


TYPED_TEST(ArrayRefTest, MakeWithAssignmentWorks)
{
    typename TestFixture::ArrayRefType arrayRef = this->a;
    this->runTests(this->a, arrayRef);
}

TYPED_TEST(ArrayRefTest, MakeWithNonConstAssignmentWorks)
{
    typename TestFixture::ArrayRefType arrayRef = this->ma;
    this->runTests(this->ma, arrayRef);
}

TYPED_TEST(ArrayRefTest, ConstructWithTemplateConstructorWorks)
{
    typename TestFixture::ArrayRefType arrayRef(this->a);
    this->runTests(this->a, arrayRef);
}

TYPED_TEST(ArrayRefTest, ConstructWithNonConstTemplateConstructorWorks)
{
    typename TestFixture::ArrayRefType arrayRef(this->ma);
    this->runTests(this->ma, arrayRef);
}

TYPED_TEST(ArrayRefTest, ConstructFromPointersWorks)
{
    typename TestFixture::ArrayRefType arrayRef(this->a, this->a + aSize);
    this->runTests(this->a, arrayRef);
}

TYPED_TEST(ArrayRefTest, ConstructFromNonConstPointersWorks)
{
    typename TestFixture::ArrayRefType arrayRef(this->ma, this->ma + aSize);
    this->runTests(this->ma, arrayRef);
}

template<bool c, typename T>
using makeConstIf_t = typename std::conditional<c, const T, T>::type;

TYPED_TEST(ArrayRefTest, ConstructFromVectorWorks)
{
    makeConstIf_t<std::is_const<typename TestFixture::ValueType>::value,
                  std::vector<typename TestFixture::NonConstValueType> > v(this->a, this->a + aSize);
    typename TestFixture::ArrayRefType                                   arrayRef(v);
    this->runTests(v.data(), arrayRef);
}

TYPED_TEST(ArrayRefTest, ConstructFromNonConstVectorWorks)
{
    std::vector<typename TestFixture::NonConstValueType> v(this->a, this->a + aSize);
    typename TestFixture::ArrayRefType                   arrayRef(v);
    this->runTests(v.data(), arrayRef);
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
    Helper<typename TestFixture::NonConstValueType> h;
    h.size = aSize;
    for (int i = 0; i != h.size; ++i)
    {
        h.a[i] = this->a[i];
    }
    typename TestFixture::ArrayRefType arrayRef(h.a);
    this->runTests(h.a, arrayRef);
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

}      // namespace gmx
