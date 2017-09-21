/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * Tests implementation of gmx::SimdMemory
 *
 * In particular we need to check that the padding works, and that
 * unpadded_end_ is maintained correctly.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_simd
 */
#include "gmxpre.h"

#include "gromacs/simd/simdmemory.h"

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"

namespace gmx
{
namespace
{

//! Templated test fixture.
template <typename T>
class SimdMemoryTest : public ::testing::Test
{
    public:
};

//! Declare allocator types to test.
using SimdMemoryTypesToTest = ::testing::Types<SimdMemory<real>
                                               /*, SimdMemory<RVec> TODO*/
                                               >;
TYPED_TEST_CASE(SimdMemoryTest, SimdMemoryTypesToTest);

// TODO Test other constructors and setPaddingRegionTo

TYPED_TEST(SimdMemoryTest, ConstructsAndResizesCorrectly)
{
    for (size_t sizeOfVector = 0; sizeOfVector <= 8; sizeOfVector++)
    {
        TypeParam v;
        EXPECT_EQ(v.size(), 0);
        EXPECT_TRUE(v.empty());

        v.resizeWithPadding(sizeOfVector);
        EXPECT_EQ(v.size(), sizeOfVector);

        // The iterable SIMD range should not exceed the actual range
        auto simdRef = v.getSimdArrayRef();
        EXPECT_LE(sizeOfVector, simdRef.size() * TypeParam::s_simdWidth) << "size of Vector is " << sizeOfVector;

        // The padding should be at least as large as 3 SIMD widths
        ArrayRef<real> paddingRef = v.getPaddingArrayRef();
        if (sizeOfVector == simdRef.size() * TypeParam::s_simdWidth)
        {
            EXPECT_EQ(paddingRef.size(), 3 * TypeParam::s_simdWidth);
        }
        else
        {
            EXPECT_GT(paddingRef.size(), 3 * TypeParam::s_simdWidth);
        }

        // The iterable range in T should be the size that was
        // requested.
        auto rvecs = v.getUnpaddedArrayRef();
        EXPECT_EQ(sizeOfVector, rvecs.size()) << "size of Vector is " << sizeOfVector;

        v.clear();
        EXPECT_EQ(0, v.size());
        EXPECT_TRUE(v.empty());
    }
}

TYPED_TEST(SimdMemoryTest, ConstructsReservesAndResizesCorrectly)
{
    for (int sizeOfVector = 0; sizeOfVector <= 8; sizeOfVector++)
    {
        TypeParam v;
        EXPECT_EQ(v.size(), 0);
        EXPECT_TRUE(v.empty());

        // Reservation should leave the memory empty
        v.reserveWithPadding(sizeOfVector);
        auto oldData = v.getUnpaddedArrayRef().data();
        EXPECT_EQ(v.size(), 0);

        // No reallocations should take place
        v.resizeWithPadding(sizeOfVector);
        EXPECT_EQ(v.size(), sizeOfVector);
        {
            auto newData = v.getUnpaddedArrayRef().data();
            EXPECT_EQ(oldData, newData);
        }

        if (sizeOfVector > 0)
        {
            // No reallocations should take place
            v.resizeWithPadding(sizeOfVector - 1);
            EXPECT_EQ(v.size(), sizeOfVector - 1);
            {
                auto newData = v.getUnpaddedArrayRef().data();
                EXPECT_EQ(oldData, newData);
            }
        }
    }
}

} // namespace
} // namespace
