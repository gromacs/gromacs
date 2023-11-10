/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * This implements basic util::array tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "nblib/util/array.hpp"

#include <vector>

#include "testutils/testasserts.h"


namespace nblib
{
namespace test
{
namespace
{

TEST(Array, construct)
{
    util::array<int, 3> a{ 0, 1, 2 };

    EXPECT_EQ(a[0], 0);
    EXPECT_EQ(a[1], 1);
    EXPECT_EQ(a[2], 2);
}

TEST(Array, plusEqual)
{
    util::array<int, 3> a{ 0, 1, 2 };
    util::array<int, 3> b{ 1, 2, 3 };

    a += b;

    EXPECT_EQ(a[0], 1);
    EXPECT_EQ(a[1], 3);
    EXPECT_EQ(a[2], 5);
}

TEST(Array, minusEqual)
{
    util::array<int, 3> a{ 0, 1, 2 };
    util::array<int, 3> b{ 1, 0, 3 };

    a -= b;

    EXPECT_EQ(a[0], -1);
    EXPECT_EQ(a[1], 1);
    EXPECT_EQ(a[2], -1);
}

TEST(Array, equal)
{
    util::array<int, 3> a{ 0, 1, 2 };
    util::array<int, 3> b{ 1, 2, 3 };

    util::array<int, 3> c{ 1, 2, 3 };

    EXPECT_FALSE(a == b);
    EXPECT_TRUE(b == c);
}

TEST(Array, unequal)
{
    util::array<int, 3> a{ 0, 1, 2 };
    util::array<int, 3> b{ 1, 2, 3 };
    util::array<int, 3> c{ 1, 2, 3 };

    EXPECT_TRUE(a != b);
    EXPECT_FALSE(b != c);
}

TEST(Array, smaller)
{
    {
        util::array<int, 3> a{ 0, 0, 0 };
        util::array<int, 3> b{ 0, 0, 0 };
        EXPECT_FALSE(a < b);
    }
    {
        util::array<int, 3> a{ 0, 0, 0 };
        util::array<int, 3> b{ 1, 0, 0 };
        EXPECT_TRUE(a < b);
        EXPECT_FALSE(b < a);
    }
    {
        util::array<int, 3> a{ 0, 0, 0 };
        util::array<int, 3> b{ 0, 1, 0 };
        EXPECT_TRUE(a < b);
        EXPECT_FALSE(b < a);
    }
    {
        util::array<int, 3> a{ 0, 0, 0 };
        util::array<int, 3> b{ 0, 0, 1 };
        EXPECT_TRUE(a < b);
        EXPECT_FALSE(b < a);
    }
    {
        util::array<int, 3> a{ 2, 4, 0 };
        util::array<int, 3> b{ 3, 0, 0 };
        EXPECT_TRUE(a < b);
        EXPECT_FALSE(b < a);
    }
}

TEST(Array, scalarMultiply)
{
    util::array<int, 3> a{ 2, 4, 6 };
    a *= 2;
    EXPECT_EQ(a[0], 4);
    EXPECT_EQ(a[1], 8);
    EXPECT_EQ(a[2], 12);
}

TEST(Array, scalarDivide)
{
    util::array<int, 3> a{ 2, 4, 6 };
    a /= 2;
    EXPECT_EQ(a[0], 1);
    EXPECT_EQ(a[1], 2);
    EXPECT_EQ(a[2], 3);
}

TEST(Array, binaryAdd)
{
    util::array<int, 3> a{ 2, 4, 0 };
    util::array<int, 3> b{ 3, 0, 1 };

    util::array<int, 3> s{ 5, 4, 1 };
    EXPECT_EQ(s, a + b);
}

TEST(Array, binarySub)
{
    util::array<int, 3> a{ 2, 4, 0 };
    util::array<int, 3> b{ 3, 0, 1 };

    util::array<int, 3> d{ -1, 4, -1 };
    EXPECT_EQ(d, a - b);
}

TEST(Array, freeScalarMultiply)
{
    util::array<int, 3> a{ 2, 4, 0 };

    util::array<int, 3> p{ 4, 8, 0 };
    EXPECT_EQ(p, a * 2);
}

TEST(Array, dot)
{
    util::array<int, 3> a{ 2, 4, 2 };
    util::array<int, 3> b{ 4, 8, 1 };
    EXPECT_EQ(dot(a, b), 42);
}

TEST(Array, negate)
{
    util::array<int, 3> a{ 2, 4, 2 };
    util::array<int, 3> b = -a;

    util::array<int, 3> ref{ -2, -4, -2 };
    EXPECT_EQ(b, ref);
}

} // namespace
} // namespace test
} // namespace nblib
