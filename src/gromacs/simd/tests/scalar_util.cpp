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
#include "gmxpre.h"

#include <cmath>
#include <cstdint>

#include <string>

#include <gtest/gtest.h>

#include "gromacs/simd/simd.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#include "testutils/testasserts.h"

#include "data.h"

namespace gmx
{
namespace test
{

namespace
{

/*! \cond internal */
/*! \addtogroup module_simd */
/*! \{ */

TEST(SimdScalarUtilTest, gatherLoadTranspose)
{
    real         data[8] = { c0, c1, c2, c3, c4, c5, c6, c7 };
    std::int32_t offset  = 1;
    real         v0, v1, v2, v3;

    gatherLoadTranspose<4>(data, &offset, &v0, &v1, &v2, &v3);

    EXPECT_EQ(data[4], v0);
    EXPECT_EQ(data[5], v1);
    EXPECT_EQ(data[6], v2);
    EXPECT_EQ(data[7], v3);

    gatherLoadTranspose<2>(data, &offset, &v0, &v1);

    EXPECT_EQ(data[2], v0);
    EXPECT_EQ(data[3], v1);
}

TEST(SimdScalarUtilTest, gatherLoadUTranspose)
{
    real         data[6] = { c0, c1, c2, c3, c4, c5 };
    std::int32_t offset  = 1;
    real         v0, v1, v2;

    gatherLoadUTranspose<3>(data, &offset, &v0, &v1, &v2);

    EXPECT_EQ(data[3], v0);
    EXPECT_EQ(data[4], v1);
    EXPECT_EQ(data[5], v2);
}

TEST(SimdScalarUtilTest, transposeScatterStoreU)
{
    real         data[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    std::int32_t offset  = 1;
    real         v0      = 1;
    real         v1      = 2;
    real         v2      = 3;

    transposeScatterStoreU<3>(data, &offset, v0, v1, v2);

    EXPECT_EQ(czero, data[0]);
    EXPECT_EQ(czero, data[1]);
    EXPECT_EQ(czero, data[2]);
    EXPECT_EQ(v0, data[3]);
    EXPECT_EQ(v1, data[4]);
    EXPECT_EQ(v2, data[5]);
    EXPECT_EQ(czero, data[6]);
    EXPECT_EQ(czero, data[7]);
    EXPECT_EQ(czero, data[8]);
}

TEST(SimdScalarUtilTest, transposeScatterIncrU)
{
    real         data[9] = { 10, 20, 30, 40, 50, 60, 70, 80, 90 };
    std::int32_t offset  = 1;
    real         v0      = c1;
    real         v1      = c2;
    real         v2      = c3;

    transposeScatterIncrU<3>(data, &offset, v0, v1, v2);

    EXPECT_EQ(real(10), data[0]);
    EXPECT_EQ(real(20), data[1]);
    EXPECT_EQ(real(30), data[2]);
    EXPECT_EQ(real(40 + c1), data[3]);
    EXPECT_EQ(real(50 + c2), data[4]);
    EXPECT_EQ(real(60 + c3), data[5]);
    EXPECT_EQ(real(70), data[6]);
    EXPECT_EQ(real(80), data[7]);
    EXPECT_EQ(real(90), data[8]);
}

TEST(SimdScalarUtilTest, transposeScatterDecrU)
{
    real         data[9] = { 10, 20, 30, 40, 50, 60, 70, 80, 90 };
    std::int32_t offset  = 1;
    real         v0      = c1;
    real         v1      = c2;
    real         v2      = c3;

    transposeScatterDecrU<3>(data, &offset, v0, v1, v2);

    EXPECT_EQ(real(10), data[0]);
    EXPECT_EQ(real(20), data[1]);
    EXPECT_EQ(real(30), data[2]);
    EXPECT_EQ(real(40 - c1), data[3]);
    EXPECT_EQ(real(50 - c2), data[4]);
    EXPECT_EQ(real(60 - c3), data[5]);
    EXPECT_EQ(real(70), data[6]);
    EXPECT_EQ(real(80), data[7]);
    EXPECT_EQ(real(90), data[8]);
}


TEST(SimdScalarTest, expandScalarsToTriplets)
{
    real scalar = c1;
    real t0, t1, t2;

    expandScalarsToTriplets(scalar, &t0, &t1, &t2);

    EXPECT_EQ(scalar, t0);
    EXPECT_EQ(scalar, t1);
    EXPECT_EQ(scalar, t2);
}

TEST(SimdScalarUtilTest, gatherLoadBySimdIntTranspose)
{
    real         data[8] = { c0, c1, c2, c3, c4, c5, c6, c7 };
    std::int32_t offset  = 1;
    real         v0, v1, v2, v3;

    gatherLoadBySimdIntTranspose<4>(data, offset, &v0, &v1, &v2, &v3);

    EXPECT_EQ(data[4], v0);
    EXPECT_EQ(data[5], v1);
    EXPECT_EQ(data[6], v2);
    EXPECT_EQ(data[7], v3);

    gatherLoadBySimdIntTranspose<2>(data, offset, &v0, &v1);

    EXPECT_EQ(data[2], v0);
    EXPECT_EQ(data[3], v1);
}

TEST(SimdScalarUtilTest, gatherLoadUBySimdIntTranspose)
{
    real         data[8] = { c0, c1, c2, c3, c4, c5, c6, c7 };
    std::int32_t offset  = 1;
    real         v0, v1;

    gatherLoadUBySimdIntTranspose<4>(data, offset, &v0, &v1);

    EXPECT_EQ(data[4], v0);
    EXPECT_EQ(data[5], v1);
}

TEST(SimdScalarUtilTest, reduceIncr4ReturnSum)
{
    real data[6] = { 0, 0, 0, 0, 0, 0 };
    real v0      = c1;
    real v1      = c2;
    real v2      = c3;
    real v3      = c4;
    real sum;

    sum = reduceIncr4ReturnSum(data + 1, v0, v1, v2, v3);

    EXPECT_EQ(czero, data[0]);
    EXPECT_EQ(v0, data[1]);
    EXPECT_EQ(v1, data[2]);
    EXPECT_EQ(v2, data[3]);
    EXPECT_EQ(v3, data[4]);
    EXPECT_EQ(czero, data[5]);

    EXPECT_REAL_EQ_TOL(real(v0 + v1 + v2 + v3), sum, defaultRealTolerance());
}

/*! \} */
/*! \endcond internal */

} // namespace
} // namespace test
} // namespace gmx
