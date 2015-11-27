/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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
 * Tests various corners of gmx::RVec implementation.
 *
 * The main point of these tests is to check that all different constructs
 * using gmx::RVec compile, and that some of the non-trivial conversions
 * to/from rvec work as intended.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/vectypes.h"

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/vec.h"

namespace
{

using gmx::RVec;

TEST(RVecTest, CanBeStoredInVector)
{
    std::vector<RVec> v;
    v.push_back(RVec(1, 2, 3));
    v.resize(2);
    EXPECT_EQ(1, v[0][XX]);
    EXPECT_EQ(2, v[0][YY]);
    EXPECT_EQ(3, v[0][ZZ]);
}

TEST(RVecTest, ConvertsImplicitlyFrom_rvec)
{
    std::vector<RVec> v;
    rvec              x = { 1, 2, 3 };
    v.push_back(x);
    EXPECT_EQ(1, v[0][XX]);
    EXPECT_EQ(2, v[0][YY]);
    EXPECT_EQ(3, v[0][ZZ]);
}

TEST(RVecTest, ConvertsImplicitlyTo_rvec)
{
    std::vector<RVec> v;
    v.push_back(RVec(1, 2, 3));
    rvec              x;
    copy_rvec(v[0], x);
    EXPECT_EQ(1, x[XX]);
    EXPECT_EQ(2, x[YY]);
    EXPECT_EQ(3, x[ZZ]);
}

TEST(RVecTest, WorksAsMutable_rvec)
{
    std::vector<RVec> v;
    v.push_back(RVec(1, 2, 3));
    rvec              x = {2, 3, 4};
    copy_rvec(x, v[0]);
    EXPECT_EQ(2, v[0][XX]);
    EXPECT_EQ(3, v[0][YY]);
    EXPECT_EQ(4, v[0][ZZ]);
}

TEST(RVecTest, WorksAs_rvec_Array)
{
    std::vector<RVec> v;
    v.push_back(RVec(1, 2, 3));
    v.push_back(RVec(2, 3, 4));
    const rvec *r = as_rvec_array(v.data());
    EXPECT_EQ(1, r[0][XX]);
    EXPECT_EQ(2, r[0][YY]);
    EXPECT_EQ(3, r[0][ZZ]);
    EXPECT_EQ(2, r[1][XX]);
    EXPECT_EQ(3, r[1][YY]);
    EXPECT_EQ(4, r[1][ZZ]);
}

/*! \brief
 * Test overloaded math operations
 */

TEST(RVecTest, CanAssignRVecToRvec)
{
    RVec a(1, 2, 3);
    RVec b;
    b = a;
    EXPECT_EQ(1, b[XX]);
    EXPECT_EQ(2, b[YY]);
    EXPECT_EQ(3, b[ZZ]);
}

TEST(RVecTest, CanAddRVecToRvec)
{
    RVec a(1, 2, 3);
    RVec b(3, 2, 1);
    RVec c;
    c = a + b;
    EXPECT_EQ(4, c[XX]);
    EXPECT_EQ(4, c[YY]);
    EXPECT_EQ(4, c[ZZ]);
}

TEST(RVecTest, CanAddAssignRVecToRvec)
{
    RVec a(1, 2, 3);
    RVec b(3, 2, 1);
    a += b;
    EXPECT_EQ(4, a[XX]);
    EXPECT_EQ(4, a[YY]);
    EXPECT_EQ(4, a[ZZ]);
}


TEST(RVecTest, CanSubstractRVecFromRvec)
{
    RVec a(1, 2, 3);
    RVec b(3, 2, 1);
    RVec c;
    c = b - a;
    EXPECT_EQ(2, c[XX]);
    EXPECT_EQ(0, c[YY]);
    EXPECT_EQ(-2, c[ZZ]);
}

TEST(RVecTest, CanSubstractAssignRVecFromRvec)
{
    RVec a(1, 2, 3);
    RVec b(3, 2, 1);
    RVec c;
    b -= a;
    EXPECT_EQ(2, b[XX]);
    EXPECT_EQ(0, b[YY]);
    EXPECT_EQ(-2, b[ZZ]);
}

TEST(RVecTest, CanDotProductRVecByRvec)
{
    RVec  a(1, 2, 3);
    RVec  b(3, 2, 1);
    float c;
    c = a.dot(b);
    EXPECT_EQ(10, c);
}

TEST(RVecTest, CanCrossProductRVecByRvec)
{
    RVec a(1, 2, 3);
    RVec b(3, 2, 1);
    RVec c;
    c = a.cross(b);
    EXPECT_EQ(-4, c[XX]);
    EXPECT_EQ(8, c[YY]);
    EXPECT_EQ(-4, c[ZZ]);
}

TEST(RVecTest, CanScaleRVec)
{
    RVec a(1, 2, 3);
    real b = 2.0;
    RVec c;
    c = a.scale(b);
    EXPECT_EQ(2, c[XX]);
    EXPECT_EQ(4, c[YY]);
    EXPECT_EQ(6, c[ZZ]);
}

/*! \brief
 * Test operations imported from vec.h
 */

TEST(RVecTest, CanClearRVec)
{
    RVec a(1, 2, 3);
    RVec b;
    b = a.clear();
    EXPECT_EQ(0, b[XX]);
    EXPECT_EQ(0, b[YY]);
    EXPECT_EQ(0, b[ZZ]);
}

TEST(RVecTest, CanNorm2OfRVec)
{
    RVec a(1, 2, 2);
    real b;
    b = a.norm2();
    EXPECT_EQ(9, b);
}

TEST(RVecTest, CanNormOfRVec)
{
    RVec a(1, 2, 2);
    real b;
    b = a.norm();
    EXPECT_EQ(3, b);
}


/*! \brief
 * Helper function for testing RVec to rvec conversions.
 */
const rvec *testFunction(const rvec &x)
{
    return &x;
}

TEST(RVecTest, WorksAs_rvec_Reference)
{
    RVec        v(1, 2, 3);
    const rvec *r = testFunction(v.as_vec());
    EXPECT_EQ(1, r[0][XX]);
    EXPECT_EQ(2, r[0][YY]);
    EXPECT_EQ(3, r[0][ZZ]);
}

} // namespace
