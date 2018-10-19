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
 * \brief Tests for gmx::ArrayRefWithPadding.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/arrayrefwithpadding.h"

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/paddedvector.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

namespace gmx
{

namespace
{

TEST(EmptyArrayRefWithPaddingTest, IsEmpty)
{
    EmptyArrayRef             emptyArray = {};
    ArrayRefWithPadding<real> empty(emptyArray);

    EXPECT_EQ(0U, empty.size());
    EXPECT_TRUE(empty.empty());
}

TEST(EmptyConstArrayRefWithPaddingTest, IsEmpty)
{
    EmptyArrayRef                   emptyArray = {};
    ArrayRefWithPadding<const real> empty(emptyArray);

    EXPECT_EQ(0U, empty.size());
    EXPECT_TRUE(empty.empty());
}

#ifdef GTEST_HAS_TYPED_TEST

//! Define the types that end up being available as TypeParam in the test cases for both kinds of ArrayRefWithPadding
typedef ::testing::Types<
        ArrayRefWithPadding<int32_t>,
        ArrayRefWithPadding<float>,
        ArrayRefWithPadding<double>
        > ArrayRefTypes;

//! Helper constant used in the text fixture.
constexpr index aSize = 3;

/*! \brief Permit all the tests to run on all kinds of ArrayRefWithPadding types
 *
 * The main objective is to verify that all the different kinds of
 * construction lead to the expected result. */
template <typename TypeParam>
class ArrayRefWithPaddingTest : public ::testing::Test
{
    public:
        typedef TypeParam ArrayRefType;
        typedef typename ArrayRefType::value_type ValueType;
        typedef PaddedVector<ValueType> PaddedVectorType;

        ValueType        a[aSize] = { ValueType(1.2), ValueType(2.4), ValueType(3.1) };
        PaddedVectorType v;

        //! Constructor, which prepares a padded vector to take array ref of.
        ArrayRefWithPaddingTest() : v(aSize)
        {
            std::copy(a, a + aSize, v.begin());
        }

        /*! \brief Run the same tests all the time
         *
         * Note that test cases must call this->runTests(), because
         * that's how the derived-class templates that implement
         * type-parameterized tests actually work. */
        void runTests(ArrayRefType  &arrayRefWithPadding)
        {
            ASSERT_LE(aSize, arrayRefWithPadding.size());
            ASSERT_FALSE(arrayRefWithPadding.empty());
            auto unpaddedArrayRef = arrayRefWithPadding.unpaddedArrayRef();
            auto paddedArrayRef   = arrayRefWithPadding.paddedArrayRef();
            EXPECT_LE(unpaddedArrayRef.size(), paddedArrayRef.size());
            for (index i = 0; i != aSize; ++i)
            {
                EXPECT_EQ(paddedArrayRef[i], unpaddedArrayRef[i]);
                EXPECT_EQ(a[i], unpaddedArrayRef[i]);
            }
        }
};

TYPED_TEST_CASE(ArrayRefWithPaddingTest, ArrayRefTypes);


TYPED_TEST(ArrayRefWithPaddingTest, AssignFromPaddedVectorWorks)
{
    typename TestFixture::ArrayRefType arrayRef = this->v.arrayRefWithPadding();
    this->runTests(arrayRef);
}

TYPED_TEST(ArrayRefWithPaddingTest, ConstructFromPointersWorks)
{
    typename TestFixture::ArrayRefType arrayRef(this->v.data(),
                                                this->v.data() + this->v.size(),
                                                this->v.data() + this->v.paddedSize());
    this->runTests(arrayRef);
}

template<bool c, typename T>
using makeConstIf_t = typename std::conditional<c, const T, T>::type;

//! Helper struct for the case actually used in mdrun signalling
template <typename T>
struct Helper
{
    public:
        T   a[3];
        int size;
};

#else   // GTEST_HAS_TYPED_TEST

/* A dummy test that at least signals that something is missing if one runs the
 * unit test executable itself.
 */
TEST(DISABLED_ArrayRefWithPaddingTest, GenericTests)
{
    ADD_FAILURE()
    << "Tests for generic ArrayRefWithPadding functionality require support for "
    << "Google Test typed tests, which was not available when the tests "
    << "were compiled.";
}

#endif // GTEST_HAS_TYPED_TEST

}      // namespace

}      // namespace gmx
