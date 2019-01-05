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
 * \brief
 * Tests multidimensional arrays
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/multidimarray.h"

#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

namespace
{

class MultiDimArrayTest : public ::testing::Test
{
    public:
        MultiDimArrayTest()
        {
            std::fill(std::begin(twoDArray_.data()), std::end(twoDArray_.data()), testNumber_ - 1);
        }
    protected:
        using container_type = std::vector<float>;
        using extents_type   = extents<3, 3>;
        using array_type     = MultiDimArray < container_type, extents_type>;

        array_type twoDArray_  = array_type::allocate();
        float      testNumber_ = 42;
};


TEST_F(MultiDimArrayTest, canConstructAndFill)
{
    for (const auto &x : twoDArray_.data())
    {
        EXPECT_EQ(testNumber_ - 1, x);
    }

}

TEST_F(MultiDimArrayTest, canSetValues)
{
    twoDArray_(1, 1) = testNumber_;
    EXPECT_EQ(testNumber_, twoDArray_(1, 1) );

}

TEST_F(MultiDimArrayTest, canCopy)
{
    const auto other       = twoDArray_;
    auto       twoDArrayIt = twoDArray_.data().begin();
    for (const auto &x : other.data())
    {
        EXPECT_EQ(*twoDArrayIt, x);
        ++twoDArrayIt;
    }
}


} // namespace

} // namespace test

} // namespace gmx
