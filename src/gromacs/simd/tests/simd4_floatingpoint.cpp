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

TEST_F(Simd4FloatingpointTest, gmxSimd4SetZeroR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom1R(0.0), simd4SetZero());
}

TEST_F(Simd4FloatingpointTest, gmxSimd4Set1R)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom1R(1.0), simd4Set1(1.0));
}

TEST_F(Simd4FloatingpointTest, gmxSimd4Load1R)
{
    real r = 2.0;
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom1R(r), simd4Load1(&r));
}

TEST_F(Simd4FloatingpointTest, gmxSimd4AddR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd4_5_7_9, simd4Add(rSimd4_1_2_3, rSimd4_4_5_6)); // 1+4=5, 2+5=7, 3+6=9
}

TEST_F(Simd4FloatingpointTest, gmxSimd4SubR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd4_4_5_6, simd4Sub(rSimd4_5_7_9, rSimd4_1_2_3)); // 5-1=4, 7-2=5, 9-3=6
}

TEST_F(Simd4FloatingpointTest, gmxSimd4MulR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(4, 10, 18), simd4Mul(rSimd4_1_2_3, rSimd4_4_5_6));
}

TEST_F(Simd4FloatingpointTest, gmxSimd4FmaddR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(11, 18, 27), simd4Fmadd(rSimd4_1_2_3, rSimd4_4_5_6, rSimd4_7_8_9)); // 1*4+7, etc.
}

TEST_F(Simd4FloatingpointTest, gmxSimd4FmsubR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(-3, 2, 9), simd4Fmsub(rSimd4_1_2_3, rSimd4_4_5_6, rSimd4_7_8_9)); // 1*4-7, etc.
}

TEST_F(Simd4FloatingpointTest, gmxSimd4FnmaddR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(3, -2, -9), simd4Fnmadd(rSimd4_1_2_3, rSimd4_4_5_6, rSimd4_7_8_9)); // -1*4+7, etc.
}

TEST_F(Simd4FloatingpointTest, gmxSimd4FnmsubR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(-11, -18, -27), simd4Fnmsub(rSimd4_1_2_3, rSimd4_4_5_6, rSimd4_7_8_9)); // -1*4-7, etc.
}

TEST_F(Simd4FloatingpointTest, gmxSimd4FabsR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd4_1_2_3, simd4Abs(rSimd4_1_2_3));    // fabs(x)=x
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd4_1_2_3, simd4Abs(rSimd4_m1_m2_m3)); // fabs(-x)=x
}

TEST_F(Simd4FloatingpointTest, gmxSimd4FnegR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd4_m1_m2_m3, simd4Neg(rSimd4_1_2_3));   // fneg(x)=-x
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd4_1_2_3,   simd4Neg(rSimd4_m1_m2_m3)); // fneg(-x)=x
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
TEST_F(Simd4FloatingpointTest, gmxSimd4AndR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom1R(1.26666259765625),
                             simd4And(simd4Set1(1.3333282470703125),
                                      simd4Set1(1.79998779296875)));
}

TEST_F(Simd4FloatingpointTest, gmxSimd4OrR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom1R(1.8666534423828125),
                             simd4Or(simd4Set1(1.3333282470703125),
                                     simd4Set1(1.79998779296875)));
}

TEST_F(Simd4FloatingpointTest, gmxSimd4XorR)
{
    /* Test xor by taking xor with a number and its negative. This should result
     * in only the sign bit being set. We then use this bit change the sign of
     * different numbers.
     */
    Simd4Real signbit = simd4Xor(simd4Set1(1.5), simd4Set1(-1.5));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(-1, 2, -3), simd4Xor(signbit, setSimd4RealFrom3R(1, -2, 3)));
}

TEST_F(Simd4FloatingpointTest, gmxSimd4AndnotR)
{
    /* Use xor (which we already tested, so fix that first if both tests fail)
     * to extract the sign bit, and then use andnot to take absolute values.
     */
    Simd4Real signbit = simd4Xor(simd4Set1(1.5), simd4Set1(-1.5));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(1, 2, 3), simd4AndNot(signbit, setSimd4RealFrom3R(-1, 2, -3)));
}

#endif

TEST_F(Simd4FloatingpointTest, gmxSimd4MaxR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(3, 2, 4), simd4Max(rSimd4_1_2_3, rSimd4_3_1_4));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(3, 2, 4), simd4Max(rSimd4_3_1_4, rSimd4_1_2_3));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(-1, -1, -3), simd4Max(rSimd4_m1_m2_m3, rSimd4_m3_m1_m4));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(-1, -1, -3), simd4Max(rSimd4_m3_m1_m4, rSimd4_m1_m2_m3));
}

TEST_F(Simd4FloatingpointTest, gmxSimd4MinR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(1, 1, 3), simd4Min(rSimd4_1_2_3, rSimd4_3_1_4));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(1, 1, 3), simd4Min(rSimd4_3_1_4, rSimd4_1_2_3));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(-3, -2, -4), simd4Min(rSimd4_m1_m2_m3, rSimd4_m3_m1_m4));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(-3, -2, -4), simd4Min(rSimd4_m3_m1_m4, rSimd4_m1_m2_m3));
}

TEST_F(Simd4FloatingpointTest, gmxSimd4RoundR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom1R(2), simd4Round(simd4Set1(2.25)));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom1R(4), simd4Round(simd4Set1(3.75)));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom1R(-2), simd4Round(simd4Set1(-2.25)));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom1R(-4), simd4Round(simd4Set1(-3.75)));
}

TEST_F(Simd4FloatingpointTest, gmxSimd4TruncR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom1R(2), simd4Trunc(rSimd4_2p25));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom1R(3), simd4Trunc(rSimd4_3p75));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom1R(-2), simd4Trunc(rSimd4_m2p25));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom1R(-3), simd4Trunc(rSimd4_m3p75));
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
    GMX_EXPECT_SIMD4_REAL_NEAR(ref, simd4Rsqrt(x));
}

TEST_F(Simd4FloatingpointTest, gmxSimd4BoolCmpEqAndBlendZeroR)
{
    Simd4Bool eq   = simd4CmpEq(rSimd4_5_7_9, rSimd4_7_8_9);
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(0, 0, 3), simd4Mask(rSimd4_1_2_3, eq));
}

TEST_F(Simd4FloatingpointTest, gmxSimd4BlendNotZeroR)
{
    Simd4Bool eq   = simd4CmpEq(rSimd4_5_7_9, rSimd4_7_8_9);
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(1, 2, 0), simd4MaskNot(rSimd4_1_2_3, eq));
}

TEST_F(Simd4FloatingpointTest, gmxSimd4BoolCmpLER)
{
    Simd4Bool le   = simd4CmpLe(rSimd4_5_7_9, rSimd4_7_8_9);
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd4_1_2_3, simd4Mask(rSimd4_1_2_3, le));
}

TEST_F(Simd4FloatingpointTest, gmxSimd4BoolCmpLTR)
{
    Simd4Bool lt   = simd4CmpLt(rSimd4_5_7_9, rSimd4_7_8_9);
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(1, 2, 0), simd4Mask(rSimd4_1_2_3, lt));
}

TEST_F(Simd4FloatingpointTest, gmxSimd4BoolAndB)
{
    Simd4Bool eq   = simd4CmpEq(rSimd4_5_7_9, rSimd4_7_8_9);
    Simd4Bool le   = simd4CmpLe(rSimd4_5_7_9, rSimd4_7_8_9);
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(0, 0, 3), simd4Mask(rSimd4_1_2_3, simd4AndB(eq, le)));
}

TEST_F(Simd4FloatingpointTest, gmxSimd4BoolOrB)
{
    Simd4Bool eq   = simd4CmpEq(rSimd4_5_7_9, rSimd4_7_8_9);
    Simd4Bool lt   = simd4CmpLt(rSimd4_5_7_9, rSimd4_7_8_9);
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(1, 2, 3), simd4Mask(rSimd4_1_2_3, simd4OrB(eq, lt)));
}

TEST_F(Simd4FloatingpointTest, gmxSimd4AnytrueB)
{
    Simd4Bool eq;

    /* this test is a bit tricky since we don't know the simd width.
     * We cannot check for truth values for "any" element beyond the first,
     * since that part of the data will not be used if simd width is 1.
     */
    eq = simd4CmpEq(rSimd4_5_7_9, setSimd4RealFrom3R(5, 0, 0));
    EXPECT_NE(0, simd4AnyTrueB(eq));

    eq = simd4CmpEq(rSimd4_1_2_3, rSimd4_4_5_6);
    EXPECT_EQ(0, simd4AnyTrueB(eq));
}

TEST_F(Simd4FloatingpointTest, gmxSimd4BlendvR)
{
    Simd4Bool lt   = simd4CmpLt(rSimd4_5_7_9, rSimd4_7_8_9);
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd4RealFrom3R(4, 5, 3), simd4Blend(rSimd4_1_2_3, rSimd4_4_5_6, lt));
}

TEST_F(Simd4FloatingpointTest, gmxSimd4ReduceR)
{
    // The horizontal sum of the SIMD variable depends on the width, so
    // simply store it an extra time and calculate what the sum should be
    std::vector<real> v   = simd4Real2Vector(rSimd4_1_2_3);
    real              sum = 0.0;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        sum += v[i];
    }

    EXPECT_EQ(sum, simd4Reduce(rSimd4_1_2_3));
}


TEST_F(Simd4FloatingpointTest, gmxSimd4DotProduct3R)
{
    Simd4Real v1 = setSimd4RealFrom3R(1, 4, 5);
    Simd4Real v2 = setSimd4RealFrom3R(3, 8, 2);
#    ifdef GMX_DOUBLE
    EXPECT_DOUBLE_EQ(45.0, simd4DotProduct(v1, v2));
#    else
    EXPECT_FLOAT_EQ(45.0, simd4DotProduct(v1, v2));
#    endif
}

#if GMX_SIMD == 2
TEST_F(Simd4FloatingpointTest, gmxSimd4TransposeR)
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

    v0 = simd4Load(p0);
    v1 = simd4Load(p1);
    v2 = simd4Load(p2);
    v3 = simd4Load(p3);

    simd4Transpose(v0, v1, v2, v3);

    simd4Store(p0, v0);
    simd4Store(p1, v1);
    simd4Store(p2, v2);
    simd4Store(p3, v3);

    for (i = 0; i < 4; i++)
    {
        EXPECT_EQ(i*10+0, p0[i]);
        EXPECT_EQ(i*10+1, p1[i]);
        EXPECT_EQ(i*10+2, p2[i]);
        EXPECT_EQ(i*10+3, p3[i]);
    }
}
#endif      // GMX_SIMD==2

#endif      // GMX_SIMD4_HAVE_REAL

/*! \} */
/*! \endcond */

}      // namespace
}      // namespace
}      // namespace

#endif // GMX_SIMD
