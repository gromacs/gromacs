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
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/simdmemory.h"

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"

namespace gmx
{
namespace
{

// TODO lots of other behaviours to test

//! Templated test fixture.
template <typename T>
class SimdMemoryTest : public ::testing::Test
{
    public:
};

//! Declare allocator types to test.
using SimdMemoryTypesToTest = ::testing::Types<SimdMemory<real>,
                                               SimdMemory<RVec>
                                               >;
TYPED_TEST_CASE(SimdMemoryTest, SimdMemoryTypesToTest);

// TODO lots of other behaviours to test

TYPED_TEST(SimdMemoryTest, ConstructsAndResizes)
{
    for (int sizeOfVector = 0; sizeOfVector <= 8; sizeOfVector++)
    {
        TypeParam v;
        v.resizeWithPadding(sizeOfVector);

        ArrayRef<SimdReal> paddedRef = v.getPaddedArrayRef();
        EXPECT_LE(sizeOfVector * TypeParam::s_packSize, paddedRef.size() * TypeParam::s_simdWidth) << "size of Vector is " << sizeOfVector;

        ArrayRef<real> paddingRef = v.getPaddingArrayRef();
        EXPECT_EQ(sizeOfVector * TypeParam::s_packSize + paddingRef.size(), paddedRef.size() * TypeParam::s_simdWidth) << "size of Vector is " << sizeOfVector;

        //ArrayRef<RVec>
        auto rvecs = v.getUnpaddedArrayRef();
        EXPECT_EQ(sizeOfVector, rvecs.size()) << "size of Vector is " << sizeOfVector;
    }
}

} // namespace
} // namespace
