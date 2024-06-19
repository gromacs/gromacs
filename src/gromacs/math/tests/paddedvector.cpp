/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * \brief
 * Tests implementation of gmx::PaddedVector
 *
 * In particular we need to check that the padding works, and that
 * unpadded_end() is maintained correctly.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/paddedvector.h"

#include <cstdint>

#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/tests/testarrayrefs.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"

namespace gmx
{
namespace test
{

//! Typed test fixture
template<typename T>
class PaddedVectorTest : public ::testing::Test
{
public:
};

//! The types used in testing
using Implementations = ::testing::Types<std::allocator<int32_t>,
                                         std::allocator<float>,
                                         std::allocator<double>,
                                         std::allocator<BasicVector<float>>,
                                         std::allocator<BasicVector<double>>,
                                         AlignedAllocator<int32_t>,
                                         AlignedAllocator<float>,
                                         AlignedAllocator<double>,
                                         AlignedAllocator<BasicVector<float>>,
                                         AlignedAllocator<BasicVector<double>>>;
TYPED_TEST_SUITE(PaddedVectorTest, Implementations);

TYPED_TEST(PaddedVectorTest, DefaultConstructorWorks)
{
    using VectorType = PaddedVector<typename TypeParam::value_type, TypeParam>;

    VectorType v;
    EXPECT_EQ(v.size(), 0);
    EXPECT_EQ(v.paddedSize(), 0);
    EXPECT_TRUE(v.empty());
    EXPECT_EQ(v.begin(), v.end());
    EXPECT_EQ(v.cbegin(), v.cend());
}

TYPED_TEST(PaddedVectorTest, ResizeWorks)
{
    using VectorType = PaddedVector<typename TypeParam::value_type, TypeParam>;

    VectorType v;

    resizeAndFillInput(&v, 3, 1);
    EXPECT_EQ(v.size(), 3);
    EXPECT_GE(v.paddedSize(), 3);
    EXPECT_EQ(v.paddedSize(), v.arrayRefWithPadding().size());
    EXPECT_LE(v.size(), v.arrayRefWithPadding().size());
}

TYPED_TEST(PaddedVectorTest, ReserveWorks)
{
    using VectorType = PaddedVector<typename TypeParam::value_type, TypeParam>;

    VectorType vReserved;
    vReserved.reserveWithPadding(5);
    resizeAndFillInput(&vReserved, 3, 1);

    EXPECT_EQ(vReserved.paddedSize(), vReserved.arrayRefWithPadding().size());
    EXPECT_LE(vReserved.size(), vReserved.arrayRefWithPadding().size());
}

TYPED_TEST(PaddedVectorTest, ReserveWorksTheSameAsNoReserve)
{
    using VectorType = PaddedVector<typename TypeParam::value_type, TypeParam>;

    VectorType v;
    resizeAndFillInput(&v, 3, 1);

    {
        SCOPED_TRACE("Test when the reservation is larger than needed");
        VectorType vReserved;
        vReserved.reserveWithPadding(5);
        resizeAndFillInput(&vReserved, 3, 1);

        EXPECT_EQ(v.size(), vReserved.size());
        EXPECT_LE(v.paddedSize(), vReserved.paddedSize());
    }
    {
        SCOPED_TRACE("Test when the reservation is smaller than needed");
        VectorType vReserved;
        vReserved.reserveWithPadding(1);
        resizeAndFillInput(&vReserved, 3, 1);

        EXPECT_EQ(v.size(), vReserved.size());
        EXPECT_GE(v.paddedSize(), vReserved.paddedSize());
    }
}

TYPED_TEST(PaddedVectorTest, MoveConstructorWorks)
{
    using VectorType = PaddedVector<typename TypeParam::value_type, TypeParam>;

    VectorType vOriginal;
    resizeAndFillInput(&vOriginal, 3, 1);

    VectorType v(std::move(vOriginal));
    EXPECT_EQ(v.size(), 3);
    EXPECT_GE(v.paddedSize(), 3);
    EXPECT_EQ(v.paddedSize(), v.arrayRefWithPadding().size());
    EXPECT_LE(v.size(), v.arrayRefWithPadding().size());
}

TYPED_TEST(PaddedVectorTest, MoveConstructorWithAllocatorWorks)
{
    using VectorType = PaddedVector<typename TypeParam::value_type, TypeParam>;

    VectorType vOriginal;
    resizeAndFillInput(&vOriginal, 3, 1);

    TypeParam  allocatorToTest;
    VectorType v(std::move(vOriginal), allocatorToTest);
    EXPECT_EQ(v.size(), 3);
    EXPECT_GE(v.paddedSize(), 3);
    EXPECT_EQ(v.paddedSize(), v.arrayRefWithPadding().size());
    EXPECT_LE(v.size(), v.arrayRefWithPadding().size());
}

TYPED_TEST(PaddedVectorTest, MoveAssignmentWorks)
{
    using VectorType = PaddedVector<typename TypeParam::value_type, TypeParam>;

    VectorType vOriginal;
    resizeAndFillInput(&vOriginal, 3, 1);

    VectorType v;
    v = std::move(vOriginal);
    EXPECT_EQ(v.size(), 3);
    EXPECT_GE(v.paddedSize(), 3);
    EXPECT_EQ(v.paddedSize(), v.arrayRefWithPadding().size());
    EXPECT_LE(v.size(), v.arrayRefWithPadding().size());
}

TYPED_TEST(PaddedVectorTest, ArrayRefConversionsAreIdentical)
{
    using VectorType = PaddedVector<typename TypeParam::value_type, TypeParam>;

    VectorType v;
    resizeAndFillInput(&v, 3, 1);

    SCOPED_TRACE("Comparing different paths to create identical unpadded views");
    compareViews(makeArrayRef(v), v.arrayRefWithPadding().unpaddedArrayRef());
    compareViews(makeConstArrayRef(v), v.constArrayRefWithPadding().unpaddedConstArrayRef());
    compareViews(makeConstArrayRef(v), v.constArrayRefWithPadding().unpaddedArrayRef());
    compareViews(makeConstArrayRef(v), v.arrayRefWithPadding().unpaddedConstArrayRef());

    SCOPED_TRACE("Comparing const to non-const unpadded views");
    compareViewsIgnoreConst(makeArrayRef(v), makeConstArrayRef(v));
}

TYPED_TEST(PaddedVectorTest, CanCopyAssign)
{
    using VectorType = PaddedVector<typename TypeParam::value_type, TypeParam>;

    VectorType v, w;
    resizeAndFillInput(&v, 3, 1);
    resizeAndFillInput(&w, 3, 2);

    w = v;
    compareViews(v.arrayRefWithPadding().unpaddedArrayRef(), w.arrayRefWithPadding().unpaddedArrayRef());
    compareViews(makeArrayRef(v), makeArrayRef(w));
}

TYPED_TEST(PaddedVectorTest, CanMoveAssign)
{
    using VectorType = PaddedVector<typename TypeParam::value_type, TypeParam>;

    VectorType v, w, x;
    resizeAndFillInput(&v, 3, 1);
    resizeAndFillInput(&w, 3, 2);
    resizeAndFillInput(&x, 3, 1);

    SCOPED_TRACE("Comparing padded views before move");
    compareViews(v.arrayRefWithPadding().unpaddedArrayRef(), x.arrayRefWithPadding().unpaddedArrayRef());
    SCOPED_TRACE("Comparing unpadded views before move");
    compareViews(makeArrayRef(v), makeArrayRef(x));

    w = std::move(x);

    SCOPED_TRACE("Comparing padded views");
    compareViews(v.arrayRefWithPadding().unpaddedArrayRef(), w.arrayRefWithPadding().unpaddedArrayRef());
    SCOPED_TRACE("Comparing unpadded views");
    compareViews(makeArrayRef(v), makeArrayRef(w));
}

TYPED_TEST(PaddedVectorTest, CanSwap)
{
    using VectorType = PaddedVector<typename TypeParam::value_type, TypeParam>;

    VectorType v, w, x;
    resizeAndFillInput(&v, 3, 1);
    resizeAndFillInput(&w, 3, 2);
    resizeAndFillInput(&x, 3, 1);

    std::swap(w, x);
    compareViews(v.arrayRefWithPadding().unpaddedArrayRef(), w.arrayRefWithPadding().unpaddedArrayRef());
    compareViews(makeArrayRef(v), makeArrayRef(w));
}

} // namespace test
} // namespace gmx
