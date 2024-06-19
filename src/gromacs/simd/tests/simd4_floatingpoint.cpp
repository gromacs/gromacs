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

#include <limits>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#include "testutils/testasserts.h"

#include "data.h"
#include "simd4.h"

#if GMX_SIMD

namespace gmx
{
namespace test
{
namespace
{

/*! \cond internal */
/*! \addtogroup module_simd */
/*! \{ */

#    if GMX_SIMD4_HAVE_REAL

/*! \brief Test fixture for SIMD4 floating-point operations (identical to the SIMD4 \ref Simd4Test) */
typedef Simd4Test Simd4FloatingpointTest;

TEST_F(Simd4FloatingpointTest, setZero)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom1R(0.0), setZero());
}

TEST_F(Simd4FloatingpointTest, set)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom1R(c1), Simd4Real(c1));
}

TEST_F(Simd4FloatingpointTest, add)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(c0 + c3, c1 + c4, c2 + c5), rSimd4_c0c1c2 + rSimd4_c3c4c5);
}

TEST_F(Simd4FloatingpointTest, sub)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(c0 - c3, c1 - c4, c2 - c5), rSimd4_c0c1c2 - rSimd4_c3c4c5);
}

TEST_F(Simd4FloatingpointTest, mul)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(c0 * c3, c1 * c4, c2 * c5), rSimd4_c0c1c2 * rSimd4_c3c4c5);
}

TEST_F(Simd4FloatingpointTest, fma)
{
    // The last bit of FMA operations depends on hardware, so we don't require exact match
    GMX_EXPECT_SIMD4_REAL_NEAR(setSimd4RealFrom3R(c0 * c3 + c6, c1 * c4 + c7, c2 * c5 + c8),
                               fma(rSimd4_c0c1c2, rSimd4_c3c4c5, rSimd4_c6c7c8));
}

TEST_F(Simd4FloatingpointTest, fms)
{
    // The last bit of FMA operations depends on hardware, so we don't require exact match
    GMX_EXPECT_SIMD4_REAL_NEAR(setSimd4RealFrom3R(c0 * c3 - c6, c1 * c4 - c7, c2 * c5 - c8),
                               fms(rSimd4_c0c1c2, rSimd4_c3c4c5, rSimd4_c6c7c8));
}

TEST_F(Simd4FloatingpointTest, fnma)
{
    // The last bit of FMA operations depends on hardware, so we don't require exact match
    GMX_EXPECT_SIMD4_REAL_NEAR(setSimd4RealFrom3R(c6 - c0 * c3, c7 - c1 * c4, c8 - c2 * c5),
                               fnma(rSimd4_c0c1c2, rSimd4_c3c4c5, rSimd4_c6c7c8));
}

TEST_F(Simd4FloatingpointTest, fnms)
{
    // The last bit of FMA operations depends on hardware, so we don't require exact match
    GMX_EXPECT_SIMD4_REAL_NEAR(setSimd4RealFrom3R(-c0 * c3 - c6, -c1 * c4 - c7, -c2 * c5 - c8),
                               fnms(rSimd4_c0c1c2, rSimd4_c3c4c5, rSimd4_c6c7c8));
}

TEST_F(Simd4FloatingpointTest, abs)
{
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd4_c0c1c2, abs(rSimd4_c0c1c2)); // fabs(x)=x
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd4_c0c1c2, abs(rSimd4_m0m1m2)); // fabs(-x)=x
}

TEST_F(Simd4FloatingpointTest, neg)
{
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd4_m0m1m2, -(rSimd4_c0c1c2)); // fneg(x)=-x
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd4_c0c1c2, -(rSimd4_m0m1m2)); // fneg(-x)=x
}

#        if GMX_SIMD_HAVE_LOGICAL
TEST_F(Simd4FloatingpointTest, and)
{
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd4_logicalResultAnd, (rSimd4_logicalA & rSimd4_logicalB));
}

TEST_F(Simd4FloatingpointTest, or)
{
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd4_logicalResultOr, (rSimd4_logicalA | rSimd4_logicalB));
}

TEST_F(Simd4FloatingpointTest, xor)
{
    /* Test xor by taking xor with a number and its negative. This should result
     * in only the sign bit being set. We then use this bit change the sign of
     * different numbers.
     */
    Simd4Real signbit = Simd4Real(c1) ^ Simd4Real(-c1);
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(-c2, c3, -c4), signbit ^ setSimd4RealFrom3R(c2, -c3, c4));
}

TEST_F(Simd4FloatingpointTest, andNot)
{
    /* Use xor (which we already tested, so fix that first if both tests fail)
     * to extract the sign bit, and then use andnot to take absolute values.
     */
    Simd4Real signbit = Simd4Real(c1) ^ Simd4Real(-c1);
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(c2, c3, c4),
                             andNot(signbit, setSimd4RealFrom3R(-c2, c3, -c4)));
}

#        endif

TEST_F(Simd4FloatingpointTest, max)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(c3, c1, c4), max(rSimd4_c0c1c2, rSimd4_c3c0c4));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(c3, c1, c4), max(rSimd4_c3c0c4, rSimd4_c0c1c2));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(-c0, -c0, -c2), max(rSimd4_m0m1m2, rSimd4_m3m0m4));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(-c0, -c0, -c2), max(rSimd4_m3m0m4, rSimd4_m0m1m2));
}

TEST_F(Simd4FloatingpointTest, min)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(c0, c0, c2), min(rSimd4_c0c1c2, rSimd4_c3c0c4));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(c0, c0, c2), min(rSimd4_c3c0c4, rSimd4_c0c1c2));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(-c3, -c1, -c4), min(rSimd4_m0m1m2, rSimd4_m3m0m4));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(-c3, -c1, -c4), min(rSimd4_m3m0m4, rSimd4_m0m1m2));
}

TEST_F(Simd4FloatingpointTest, round)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom1R(2), round(Simd4Real(2.25)));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom1R(4), round(Simd4Real(3.75)));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom1R(-2), round(Simd4Real(-2.25)));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom1R(-4), round(Simd4Real(-3.75)));
}

TEST_F(Simd4FloatingpointTest, trunc)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom1R(2), trunc(rSimd4_2p25));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom1R(3), trunc(rSimd4_3p75));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom1R(-2), trunc(rSimd4_m2p25));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom1R(-3), trunc(rSimd4_m3p75));
}

/* We do extensive 1/sqrt(x) and 1/x accuracy testing in the tests for
 * the SIMD math functions, so we just make sure the lookup instructions
 * appear to work for a few values here.
 */
TEST_F(Simd4FloatingpointTest, gmxSimd4RsqrtR)
{
    Simd4Real x   = setSimd4RealFrom3R(4.0, M_PI, 1234567890.0);
    Simd4Real ref = setSimd4RealFrom3R(0.5, 1.0 / std::sqrt(M_PI), 1.0 / std::sqrt(1234567890.0));
    int       shiftbits = std::numeric_limits<real>::digits - GMX_SIMD_RSQRT_BITS;

    if (shiftbits < 0)
    {
        shiftbits = 0;
    }

    // The allowed Ulp deviation is 2 to the power of the number of mantissa
    // digits, minus the number of bits provided by the table lookup
    setUlpTol(1LL << shiftbits);
    GMX_EXPECT_SIMD4_REAL_NEAR(ref, rsqrt(x));
}

TEST_F(Simd4FloatingpointTest, cmpEqAndSelectByMask)
{
    Simd4Bool eq = rSimd4_c4c6c8 == rSimd4_c6c7c8;
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(0, 0, c2), selectByMask(rSimd4_c0c1c2, eq));
}

TEST_F(Simd4FloatingpointTest, selectByNotMask)
{
    Simd4Bool eq = rSimd4_c4c6c8 == rSimd4_c6c7c8;
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(c0, c1, 0), selectByNotMask(rSimd4_c0c1c2, eq));
}

TEST_F(Simd4FloatingpointTest, cmpNe)
{
    Simd4Bool eq = rSimd4_c4c6c8 != rSimd4_c6c7c8;
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(c0, c1, 0), selectByMask(rSimd4_c0c1c2, eq));
}

TEST_F(Simd4FloatingpointTest, cmpLe)
{
    Simd4Bool le = rSimd4_c4c6c8 <= rSimd4_c6c7c8;
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd4_c0c1c2, selectByMask(rSimd4_c0c1c2, le));
}

TEST_F(Simd4FloatingpointTest, cmpLt)
{
    Simd4Bool lt = rSimd4_c4c6c8 < rSimd4_c6c7c8;
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(c0, c1, 0), selectByMask(rSimd4_c0c1c2, lt));
}

TEST_F(Simd4FloatingpointTest, andB)
{
    Simd4Bool eq = rSimd4_c4c6c8 == rSimd4_c6c7c8;
    Simd4Bool le = rSimd4_c4c6c8 <= rSimd4_c6c7c8;
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(0, 0, c2), selectByMask(rSimd4_c0c1c2, (eq && le)));
}

TEST_F(Simd4FloatingpointTest, orB)
{
    Simd4Bool eq = rSimd4_c4c6c8 == rSimd4_c6c7c8;
    Simd4Bool lt = rSimd4_c4c6c8 < rSimd4_c6c7c8;
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd4_c0c1c2, selectByMask(rSimd4_c0c1c2, (eq || lt)));
}

TEST_F(Simd4FloatingpointTest, anyTrue)
{
    Simd4Bool eq;

    /* this test is a bit tricky since we don't know the simd width.
     * We cannot check for truth values for "any" element beyond the first,
     * since that part of the data will not be used if simd width is 1.
     */
    eq = (rSimd4_c4c6c8 == setSimd4RealFrom3R(c4, 0, 0));
    EXPECT_TRUE(anyTrue(eq));

    eq = (rSimd4_c0c1c2 == rSimd4_c3c4c5);
    EXPECT_FALSE(anyTrue(eq));
}

TEST_F(Simd4FloatingpointTest, blend)
{
    Simd4Bool lt = rSimd4_c4c6c8 < rSimd4_c6c7c8;
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(c3, c4, c2), blend(rSimd4_c0c1c2, rSimd4_c3c4c5, lt));
}

TEST_F(Simd4FloatingpointTest, reduce)
{
    // The horizontal sum of the SIMD variable depends on the width, so
    // simply store it an extra time and calculate what the sum should be
    std::vector<real> v   = simd4Real2Vector(rSimd4_c3c4c5);
    real              sum = 0.0;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        sum += v[i];
    }

    EXPECT_REAL_EQ_TOL(sum, reduce(rSimd4_c3c4c5), defaultRealTolerance());
}


TEST_F(Simd4FloatingpointTest, dotProduct)
{
    real res = c0 * c3 + c1 * c4 + c2 * c5;

    EXPECT_REAL_EQ_TOL(res, dotProduct(rSimd4_c0c1c2, rSimd4_c3c4c5), defaultRealTolerance());
}

TEST_F(Simd4FloatingpointTest, transpose)
{
    Simd4Real v0, v1, v2, v3;
    int       i;
    // aligned pointers
    alignas(GMX_SIMD_ALIGNMENT) real p0[4 * GMX_SIMD4_WIDTH];
    real*                            p1 = p0 + GMX_SIMD4_WIDTH;
    real*                            p2 = p0 + 2 * GMX_SIMD4_WIDTH;
    real*                            p3 = p0 + 3 * GMX_SIMD4_WIDTH;

    // Assign data with tens as row, single-digit as column
    for (i = 0; i < 4; i++)
    {
        // Scale by 1+100*eps to use low bits tii
        p0[i] = (0 * 10 + i * 1) * (1.0 + 100 * GMX_REAL_EPS);
        p1[i] = (1 * 10 + i * 1) * (1.0 + 100 * GMX_REAL_EPS);
        p2[i] = (2 * 10 + i * 1) * (1.0 + 100 * GMX_REAL_EPS);
        p3[i] = (3 * 10 + i * 1) * (1.0 + 100 * GMX_REAL_EPS);
    }

    v0 = load4(p0);
    v1 = load4(p1);
    v2 = load4(p2);
    v3 = load4(p3);

    transpose(&v0, &v1, &v2, &v3);

    store4(p0, v0);
    store4(p1, v1);
    store4(p2, v2);
    store4(p3, v3);

    for (i = 0; i < 4; i++)
    {
        EXPECT_REAL_EQ_TOL((i * 10 + 0) * (1.0 + 100 * GMX_REAL_EPS), p0[i], defaultRealTolerance());
        EXPECT_REAL_EQ_TOL((i * 10 + 1) * (1.0 + 100 * GMX_REAL_EPS), p1[i], defaultRealTolerance());
        EXPECT_REAL_EQ_TOL((i * 10 + 2) * (1.0 + 100 * GMX_REAL_EPS), p2[i], defaultRealTolerance());
        EXPECT_REAL_EQ_TOL((i * 10 + 3) * (1.0 + 100 * GMX_REAL_EPS), p3[i], defaultRealTolerance());
    }
}

#    endif // GMX_SIMD4_HAVE_REAL

/*! \} */
/*! \endcond */

} // namespace
} // namespace test
} // namespace gmx

#endif // GMX_SIMD
