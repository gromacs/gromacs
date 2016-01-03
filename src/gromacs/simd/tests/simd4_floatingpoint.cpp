/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include <cmath>

#include "gromacs/math/utilities.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/basedefinitions.h"

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

#if GMX_SIMD4_HAVE_REAL

/*! \brief Test fixture for SIMD4 floating-point operations (identical to the SIMD4 \ref Simd4Test) */
typedef Simd4Test Simd4FloatingpointTest;

TEST_F(Simd4FloatingpointTest, setZero)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom1R(0.0), setZero());
}

TEST_F(Simd4FloatingpointTest, set)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom1R(1.0), Simd4Real(1.0));
}

TEST_F(Simd4FloatingpointTest, add)
{
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd4_5_7_9, rSimd4_1_2_3 + rSimd4_4_5_6); // 1+4=5, 2+5=7, 3+6=9
}

TEST_F(Simd4FloatingpointTest, sub)
{
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd4_4_5_6, rSimd4_5_7_9 - rSimd4_1_2_3); // 5-1=4, 7-2=5, 9-3=6
}

TEST_F(Simd4FloatingpointTest, mul)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(4, 10, 18), rSimd4_1_2_3 * rSimd4_4_5_6);
}

TEST_F(Simd4FloatingpointTest, fma)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(11, 18, 27), fma(rSimd4_1_2_3, rSimd4_4_5_6, rSimd4_7_8_9)); // 1*4+7, etc.
}

TEST_F(Simd4FloatingpointTest, fms)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(-3, 2, 9), fms(rSimd4_1_2_3, rSimd4_4_5_6, rSimd4_7_8_9)); // 1*4-7, etc.
}

TEST_F(Simd4FloatingpointTest, fnma)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(3, -2, -9), fnma(rSimd4_1_2_3, rSimd4_4_5_6, rSimd4_7_8_9)); // -1*4+7, etc.
}

TEST_F(Simd4FloatingpointTest, fnms)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(-11, -18, -27), fnms(rSimd4_1_2_3, rSimd4_4_5_6, rSimd4_7_8_9)); // -1*4-7, etc.
}

TEST_F(Simd4FloatingpointTest, abs)
{
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd4_1_2_3, abs(rSimd4_1_2_3));    // fabs(x)=x
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd4_1_2_3, abs(rSimd4_m1_m2_m3)); // fabs(-x)=x
}

TEST_F(Simd4FloatingpointTest, neg)
{
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd4_m1_m2_m3, -rSimd4_1_2_3);   // fneg(x)=-x
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd4_1_2_3,   -rSimd4_m1_m2_m3); // fneg(-x)=x
}

#if GMX_SIMD_HAVE_LOGICAL
/* 1.3333282470703125 has mantissa 0101010101010101 (followed by zeros)
 * 1.79998779296875   has mantissa 1100110011001100 (followed by zeros)
 * 1.26666259765625   has mantissa 0100010001000100 (followed by zeros)
 * 1.8666534423828125 has mantissa 1101110111011101 (followed by zeros)
 *
 * Since all of them have the same exponent (2^0), the exponent will
 * not change with AND or OR operations.
 */
TEST_F(Simd4FloatingpointTest, and)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom1R(1.26666259765625),
                             Simd4Real(1.3333282470703125) & Simd4Real(1.79998779296875));
}

TEST_F(Simd4FloatingpointTest, or)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom1R(1.8666534423828125),
                             Simd4Real(1.3333282470703125) | Simd4Real(1.79998779296875));
}

TEST_F(Simd4FloatingpointTest, xor)
{
    /* Test xor by taking xor with a number and its negative. This should result
     * in only the sign bit being set. We then use this bit change the sign of
     * different numbers.
     */
    Simd4Real signbit = Simd4Real(1.5) ^ Simd4Real(-1.5);
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(-1, 2, -3), signbit ^ setSimd4RealFrom3R(1, -2, 3));
}

TEST_F(Simd4FloatingpointTest, andNot)
{
    /* Use xor (which we already tested, so fix that first if both tests fail)
     * to extract the sign bit, and then use andnot to take absolute values.
     */
    Simd4Real signbit = Simd4Real(1.5) ^ Simd4Real(-1.5);
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(1, 2, 3), andNot(signbit, setSimd4RealFrom3R(-1, 2, -3)));
}

#endif

TEST_F(Simd4FloatingpointTest, max)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(3, 2, 4), max(rSimd4_1_2_3, rSimd4_3_1_4));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(3, 2, 4), max(rSimd4_3_1_4, rSimd4_1_2_3));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(-1, -1, -3), max(rSimd4_m1_m2_m3, rSimd4_m3_m1_m4));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(-1, -1, -3), max(rSimd4_m3_m1_m4, rSimd4_m1_m2_m3));
}

TEST_F(Simd4FloatingpointTest, min)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(1, 1, 3), min(rSimd4_1_2_3, rSimd4_3_1_4));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(1, 1, 3), min(rSimd4_3_1_4, rSimd4_1_2_3));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(-3, -2, -4), min(rSimd4_m1_m2_m3, rSimd4_m3_m1_m4));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(-3, -2, -4), min(rSimd4_m3_m1_m4, rSimd4_m1_m2_m3));
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
    Simd4Real        x                   = setSimd4RealFrom3R(4.0, M_PI, 1234567890.0);
    Simd4Real        ref                 = setSimd4RealFrom3R(0.5, 1.0/std::sqrt(M_PI), 1.0/std::sqrt(1234567890.0));
    int              shiftbits           = std::numeric_limits<real>::digits-GMX_SIMD_RSQRT_BITS;

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
    Simd4Bool eq   = (rSimd4_5_7_9 == rSimd4_7_8_9);
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(0, 0, 3), selectByMask(rSimd4_1_2_3, eq));
}

TEST_F(Simd4FloatingpointTest, selectByNotMask)
{
    Simd4Bool eq   = (rSimd4_5_7_9 == rSimd4_7_8_9);
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(1, 2, 0), selectByNotMask(rSimd4_1_2_3, eq));
}

TEST_F(Simd4FloatingpointTest, cmpNe)
{
    Simd4Bool eq   = (rSimd4_5_7_9 != rSimd4_7_8_9);
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(1, 2, 0), selectByMask(rSimd4_1_2_3, eq));
}

TEST_F(Simd4FloatingpointTest, cmpLe)
{
    Simd4Bool le   = (rSimd4_5_7_9 <= rSimd4_7_8_9);
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd4_1_2_3, selectByMask(rSimd4_1_2_3, le));
}

TEST_F(Simd4FloatingpointTest, cmpLt)
{
    Simd4Bool lt   = (rSimd4_5_7_9 < rSimd4_7_8_9);
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(1, 2, 0), selectByMask(rSimd4_1_2_3, lt));
}

TEST_F(Simd4FloatingpointTest, andB)
{
    Simd4Bool eq   = (rSimd4_5_7_9 == rSimd4_7_8_9);
    Simd4Bool le   = (rSimd4_5_7_9 <= rSimd4_7_8_9);
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(0, 0, 3), selectByMask(rSimd4_1_2_3, eq && le));
}

TEST_F(Simd4FloatingpointTest, orB)
{
    Simd4Bool eq   = (rSimd4_5_7_9 == rSimd4_7_8_9);
    Simd4Bool lt   = (rSimd4_5_7_9 < rSimd4_7_8_9);
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(1, 2, 3), selectByMask(rSimd4_1_2_3, eq || lt));
}

TEST_F(Simd4FloatingpointTest, anyTrue)
{
    Simd4Bool eq;

    /* this test is a bit tricky since we don't know the simd width.
     * We cannot check for truth values for "any" element beyond the first,
     * since that part of the data will not be used if simd width is 1.
     */
    eq = (rSimd4_5_7_9 == setSimd4RealFrom3R(5, 0, 0));
    EXPECT_TRUE(anyTrue(eq));

    eq = (rSimd4_1_2_3 == rSimd4_4_5_6);
    EXPECT_FALSE(anyTrue(eq));
}

TEST_F(Simd4FloatingpointTest, blend)
{
    Simd4Bool lt   = (rSimd4_5_7_9 < rSimd4_7_8_9);
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(4, 5, 3), blend(rSimd4_1_2_3, rSimd4_4_5_6, lt));
}

TEST_F(Simd4FloatingpointTest, reduce)
{
    // The horizontal sum of the SIMD variable depends on the width, so
    // simply store it an extra time and calculate what the sum should be
    std::vector<real> v   = simd4Real2Vector(rSimd4_1_2_3);
    real              sum = 0.0;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        sum += v[i];
    }

    EXPECT_EQ(sum, reduce(rSimd4_1_2_3));
}


TEST_F(Simd4FloatingpointTest, dotProduct)
{
    Simd4Real v1 = setSimd4RealFrom3R(1, 4, 5);
    Simd4Real v2 = setSimd4RealFrom3R(3, 8, 2);
#    if GMX_DOUBLE
    EXPECT_DOUBLE_EQ(45.0, dotProduct(v1, v2));
#    else
    EXPECT_FLOAT_EQ(45.0, dotProduct(v1, v2));
#    endif
}

TEST_F(Simd4FloatingpointTest, transpose)
{
    Simd4Real        v0, v1, v2, v3;
    int              i;
    // aligned pointers
    GMX_ALIGNED(real, GMX_SIMD4_WIDTH) p0[4*GMX_SIMD4_WIDTH];
    real          *  p1 = p0 + GMX_SIMD4_WIDTH;
    real          *  p2 = p0 + 2*GMX_SIMD4_WIDTH;
    real          *  p3 = p0 + 3*GMX_SIMD4_WIDTH;

    // Assign data with tens as row, single-digit as column
    for (i = 0; i < 4; i++)
    {
        p0[i] = 0*10 + i*1;
        p1[i] = 1*10 + i*1;
        p2[i] = 2*10 + i*1;
        p3[i] = 3*10 + i*1;
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
        EXPECT_EQ(i*10+0, p0[i]);
        EXPECT_EQ(i*10+1, p1[i]);
        EXPECT_EQ(i*10+2, p2[i]);
        EXPECT_EQ(i*10+3, p3[i]);
    }
}

#endif      // GMX_SIMD4_HAVE_REAL

/*! \} */
/*! \endcond */

}      // namespace
}      // namespace
}      // namespace

#endif // GMX_SIMD
