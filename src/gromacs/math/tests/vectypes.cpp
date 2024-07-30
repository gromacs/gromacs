/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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

#include <array>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/vec.h"
#include "gromacs/utility/real.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

using gmx::DVec;
using gmx::IVec;
using gmx::RVec;

TEST(RVecTest, CanBeStoredInVector)
{
    std::vector<RVec> v;
    v.emplace_back(1, 2, 3);
    v.resize(2);
    EXPECT_EQ(1, v[0][XX]);
    EXPECT_EQ(2, v[0][YY]);
    EXPECT_EQ(3, v[0][ZZ]);
}

TEST(RVecTest, ConvertsImplicitlyFrom_rvec)
{
    std::vector<RVec> v;
    rvec              x = { 1, 2, 3 };
    v.emplace_back(x);
    EXPECT_EQ(1, v[0][XX]);
    EXPECT_EQ(2, v[0][YY]);
    EXPECT_EQ(3, v[0][ZZ]);
}

TEST(RVecTest, ConvertsImplicitlyTo_rvec)
{
    std::vector<RVec> v;
    v.emplace_back(1, 2, 3);
    rvec x;
    copy_rvec(v[0], x);
    EXPECT_EQ(1, x[XX]);
    EXPECT_EQ(2, x[YY]);
    EXPECT_EQ(3, x[ZZ]);
}

TEST(RVecTest, WorksAsMutable_rvec)
{
    std::vector<RVec> v;
    v.emplace_back(1, 2, 3);
    rvec x = { 2, 3, 4 };
    copy_rvec(x, v[0]);
    EXPECT_EQ(2, v[0][XX]);
    EXPECT_EQ(3, v[0][YY]);
    EXPECT_EQ(4, v[0][ZZ]);
}

TEST(RVecTest, WorksAs_rvec_Array)
{
    std::vector<RVec> v;
    v.emplace_back(1, 2, 3);
    v.emplace_back(2, 3, 4);
    const rvec* r = as_rvec_array(v.data());
    EXPECT_EQ(1, r[0][XX]);
    EXPECT_EQ(2, r[0][YY]);
    EXPECT_EQ(3, r[0][ZZ]);
    EXPECT_EQ(2, r[1][XX]);
    EXPECT_EQ(3, r[1][YY]);
    EXPECT_EQ(4, r[1][ZZ]);
}

TEST(RVecTest, ComparesEqual)
{
    RVec a(1, 2, 3);
    RVec b(1, 2, 3);
    RVec c(1, 2, 2);
    EXPECT_EQ(a == b, true);
    EXPECT_EQ(a == c, false);
}

TEST(RVecTest, ComparesUnequal)
{
    RVec a(1, 2, 3);
    RVec b(1, 2, 3);
    RVec c(1, 2, 2);
    EXPECT_EQ(a != b, false);
    EXPECT_EQ(a != c, true);
}

/*! \brief
 * Test overloaded math operations
 */

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


TEST(RVecTest, CanSubtractRVecFromRvec)
{
    RVec a(1, 2, 3);
    RVec b(3, 2, 1);
    RVec c;
    c = b - a;
    EXPECT_EQ(2, c[XX]);
    EXPECT_EQ(0, c[YY]);
    EXPECT_EQ(-2, c[ZZ]);
}

TEST(RVecTest, CanSubtractAssignRVecFromRvec)
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

/*! \brief
 * Test for inplace operations imported from vec.h
 */

TEST(RVecTest, CanDivideRVecInplace)
{
    RVec a(1, 2, 3);
    real b = 0.5;
    RVec c;
    c = a / b;
    EXPECT_EQ(2, c[XX]);
    EXPECT_EQ(4, c[YY]);
    EXPECT_EQ(6, c[ZZ]);
}

TEST(RVecTest, CanScaleRVec)
{
    RVec a(1, 2, 3);
    real b = 2.0;
    a *= b;
    EXPECT_EQ(2, a[XX]);
    EXPECT_EQ(4, a[YY]);
    EXPECT_EQ(6, a[ZZ]);
}

TEST(RVecTest, CanDivideRVec)
{
    RVec a(1, 2, 3);
    real b = 0.5;
    a /= b;
    EXPECT_EQ(2, a[XX]);
    EXPECT_EQ(4, a[YY]);
    EXPECT_EQ(6, a[ZZ]);
}

TEST(RVecTest, CanDoUnitvFromRVec)
{
    RVec a(3, 0, 0);
    RVec b;
    b = a.unitVector();
    EXPECT_REAL_EQ(1, b[XX]);
    EXPECT_REAL_EQ(0, b[YY]);
    EXPECT_REAL_EQ(0, b[ZZ]);
}

TEST(RVecTest, CanSqLengthOfRVec)
{
    RVec a(1, 2, 2);
    real b;
    b = a.norm2();
    EXPECT_REAL_EQ(9, b);
}

TEST(RVecTest, CanLengthOfRVec)
{
    RVec a(1, 2, 2);
    real b;
    b = a.norm();
    EXPECT_REAL_EQ(3, b);
}

TEST(RVecTest, CanCastToRVec)
{
    DVec a(1, 2, 2);
    RVec b;
    b = a.toRVec();
    EXPECT_EQ(1, b[XX]);
    EXPECT_EQ(2, b[YY]);
    EXPECT_EQ(2, b[ZZ]);
}

TEST(RVecTest, CanCastToDVec)
{
    RVec a(1, 2, 2);
    DVec b;
    b = a.toDVec();
    EXPECT_EQ(1, b[XX]);
    EXPECT_EQ(2, b[YY]);
    EXPECT_EQ(2, b[ZZ]);
}


/*! \brief
 * Tests for out of class functions
 */
TEST(RVecTest, CanLeftScalarMultiply)
{
    RVec a(1, 2, 3);
    real b = 2.0;
    RVec c;
    c = b * a;
    EXPECT_EQ(2, c[XX]);
    EXPECT_EQ(4, c[YY]);
    EXPECT_EQ(6, c[ZZ]);
}

TEST(RVecTest, CanRightScalarMultiply)
{
    RVec a(1, 2, 3);
    real b = 2.0;
    RVec c;
    c = a * b;
    EXPECT_EQ(2, c[XX]);
    EXPECT_EQ(4, c[YY]);
    EXPECT_EQ(6, c[ZZ]);
}

TEST(RVecTest, CanGetUnitvFromRVec)
{
    RVec a(3, 0, 0);
    RVec b;
    b = gmx::unitVector(a);
    EXPECT_REAL_EQ(1, b[XX]);
    EXPECT_REAL_EQ(0, b[YY]);
    EXPECT_REAL_EQ(0, b[ZZ]);
}

TEST(RVecTest, CanGetSqLengthOfRVec)
{
    RVec a(1, 2, 2);
    real b;
    b = gmx::norm2(a);
    EXPECT_REAL_EQ(9, b);
}

TEST(RVecTest, CanGetLengthOfRVec)
{
    RVec a(1, 2, 2);
    real b;
    b = gmx::norm(a);
    EXPECT_REAL_EQ(3, b);
}

TEST(RVecTest, CanDoCrossProductOfRVec)
{
    RVec a(1, 2, 3);
    RVec b(3, 2, 1);
    RVec c;
    c = gmx::cross(a, b);
    EXPECT_EQ(-4, c[XX]);
    EXPECT_EQ(8, c[YY]);
    EXPECT_EQ(-4, c[ZZ]);
}

TEST(RVecTest, CanDoDotProductOfRVec)
{
    RVec  a(1, 2, 3);
    RVec  b(3, 2, 1);
    float c;
    c = gmx::dot(a, b);
    EXPECT_EQ(10, c);
}

TEST(RVecTest, CanScaleByVector)
{
    RVec a(1, 2, 3);
    RVec b(3, 2, 1);
    RVec scaled = scaleByVector(a, b);
    EXPECT_REAL_EQ(3, scaled[XX]);
    EXPECT_REAL_EQ(4, scaled[YY]);
    EXPECT_REAL_EQ(3, scaled[ZZ]);
}

TEST(RVecTest, CanNegate)
{
    RVec a(1, -2, 3);
    RVec b = -a;
    EXPECT_REAL_EQ(-1, b[XX]);
    EXPECT_REAL_EQ(2, b[YY]);
    EXPECT_REAL_EQ(-3, b[ZZ]);
}

TEST(RVecTest, asIVec)
{
    RVec a(1.2, 2.7, -3e3);
    auto asIvec = a.toIVec();

    EXPECT_REAL_EQ(1, asIvec[XX]);
    EXPECT_REAL_EQ(2, asIvec[YY]);
    EXPECT_REAL_EQ(-3000, asIvec[ZZ]);
}

TEST(RVecTest, elementWiseMin)
{
    RVec a(1, 2, 3);
    RVec b(3, 2, 1);
    auto minAB = elementWiseMin(a, b);

    EXPECT_REAL_EQ(1, minAB[XX]);
    EXPECT_REAL_EQ(2, minAB[YY]);
    EXPECT_REAL_EQ(1, minAB[ZZ]);
}

TEST(RVecTest, elementWiseMax)
{
    RVec a(1, 2, 3);
    RVec b(3, 2, 1);
    auto maxAB = elementWiseMax(a, b);

    EXPECT_REAL_EQ(3, maxAB[XX]);
    EXPECT_REAL_EQ(2, maxAB[YY]);
    EXPECT_REAL_EQ(3, maxAB[ZZ]);
}

/*! \brief
 * Helper function for testing DVec to dvec conversions.
 */
const dvec* testFunction(const dvec& x)
{
    return &x;
}

TEST(RVecTest, WorksAs_dvec_Reference)
{
    DVec        v(1, 2, 3);
    const dvec* r = testFunction(v.as_vec());
    EXPECT_EQ(1, r[0][XX]);
    EXPECT_EQ(2, r[0][YY]);
    EXPECT_EQ(3, r[0][ZZ]);
}

/*! \brief
 * Helper function for testing IVec to ivec conversions.
 */
const ivec* testFunction(const ivec& x)
{
    return &x;
}

TEST(RVecTest, WorksAs_ivec_Reference)
{
    IVec        v(1, 2, 3);
    const ivec* r = testFunction(v.as_vec());
    EXPECT_EQ(1, r[0][XX]);
    EXPECT_EQ(2, r[0][YY]);
    EXPECT_EQ(3, r[0][ZZ]);
}

/*! \brief
 * Helper function for testing RVec to rvec conversions.
 */
#if !GMX_DOUBLE // otherwise rvec==dvec
const rvec* testFunction(const rvec& x)
{
    return &x;
}
#endif

TEST(RVecTest, WorksAs_rvec_Reference)
{
    RVec        v(1, 2, 3);
    const rvec* r = testFunction(v);
    EXPECT_EQ(1, r[0][XX]);
    EXPECT_EQ(2, r[0][YY]);
    EXPECT_EQ(3, r[0][ZZ]);
}

TEST(RVecTest, CopyConstructorWorks)
{
    RVec v(1, 2, 3);
    RVec copy(v);
    EXPECT_EQ(1, copy[XX]);
    EXPECT_EQ(2, copy[YY]);
    EXPECT_EQ(3, copy[ZZ]);
}

TEST(RVecTest, CopyAssignmentWorks)
{
    RVec v(1, 2, 3);
    RVec copy;
    copy = v;
    EXPECT_EQ(1, copy[XX]);
    EXPECT_EQ(2, copy[YY]);
    EXPECT_EQ(3, copy[ZZ]);
}

TEST(RVecTest, MoveConstructorWorks)
{
    RVec v(1, 2, 3);
    // We are testing for correctness and don't care about performance
    // NOLINTNEXTLINE(performance-move-const-arg)
    RVec copy(std::move(v));
    EXPECT_EQ(1, copy[XX]);
    EXPECT_EQ(2, copy[YY]);
    EXPECT_EQ(3, copy[ZZ]);
}

TEST(RVecTest, MoveAssignmentWorks)
{
    RVec v(1, 2, 3);
    RVec copy;
    // We are testing for correctness and don't care about performance
    // NOLINTNEXTLINE(performance-move-const-arg)
    copy = std::move(v);
    EXPECT_EQ(1, copy[XX]);
    EXPECT_EQ(2, copy[YY]);
    EXPECT_EQ(3, copy[ZZ]);
}

TEST(RVecTest, UsableInConstexpr)
{
    // Check that we can use gmx::RVec as constexpr and common operations work
    constexpr std::array<RVec, 2> a{ RVec{ 0, 1, 2 }, RVec{ -1, -2, -3.3 } };
    static_assert(a[0][0] == 0);
    constexpr RVec b = [](RVec v) {
        v *= 2;
        return v;
    }(a[0]);
    static_assert(b == RVec{ 0, 2, 4 });
    static_assert(b.toIVec() == IVec{ 0, 2, 4 });
    static_assert(b.toDVec() == DVec{ 0, 2, 4 });
    static_assert(scaleByVector(a[0], a[0]) == RVec{ 0, 1, 4 });
    static_assert(a[0].norm2() == 5);
    static_assert(a[0].dot(b) == 10);
    static_assert(elementWiseMax(a[0], a[1]) == a[0]);
    static_assert(elementWiseMin(a[0], a[1]) == a[1]);
}

} // namespace
} // namespace test
} // namespace gmx
