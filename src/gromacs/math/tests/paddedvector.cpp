/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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

#include <memory>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/basedefinitions.h"

#include "gromacs/math/tests/testarrayrefs.h"

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
TYPED_TEST_CASE(PaddedVectorTest, Implementations);

TYPED_TEST(PaddedVectorTest, ConstructsResizesAndReserves)
{
    using VectorType = PaddedVector<typename TypeParam::value_type, TypeParam>;

    VectorType v;
    fillInput(&v, 1);

    EXPECT_EQ(v.size(), v.size());
    EXPECT_EQ(v.paddedSize(), v.arrayRefWithPadding().size());
    EXPECT_LE(v.size(), v.arrayRefWithPadding().size());

    VectorType vReserved;
    vReserved.reserveWithPadding(5);
    fillInput(&vReserved, 1);

    EXPECT_EQ(vReserved.size(), vReserved.size());
    EXPECT_EQ(vReserved.paddedSize(), vReserved.arrayRefWithPadding().size());
    EXPECT_LE(vReserved.size(), vReserved.arrayRefWithPadding().size());

    EXPECT_LE(v.paddedSize(), vReserved.paddedSize());
}

TYPED_TEST(PaddedVectorTest, CanCopyAssign)
{
    using VectorType = PaddedVector<typename TypeParam::value_type, TypeParam>;

    VectorType v, w;
    fillInput(&v, 1);
    fillInput(&w, 2);

    w = v;
    compareViews(v.arrayRefWithPadding().unpaddedArrayRef(), w.arrayRefWithPadding().unpaddedArrayRef());
    compareViews(makeArrayRef(v), makeArrayRef(v));
}

TYPED_TEST(PaddedVectorTest, CanMoveAssign)
{
    using VectorType = PaddedVector<typename TypeParam::value_type, TypeParam>;

    VectorType v, w, x;
    fillInput(&v, 1);
    fillInput(&w, 2);
    fillInput(&x, 1);

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
    fillInput(&v, 1);
    fillInput(&w, 2);
    fillInput(&x, 1);

    std::swap(w, x);
    compareViews(v.arrayRefWithPadding().unpaddedArrayRef(), w.arrayRefWithPadding().unpaddedArrayRef());
    compareViews(makeArrayRef(v), makeArrayRef(w));
}

} // namespace test
} // namespace gmx
