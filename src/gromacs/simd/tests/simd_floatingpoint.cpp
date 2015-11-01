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

#include <math.h>

#include "gromacs/math/utilities.h"
#include "gromacs/simd/simd.h"

#include "simd.h"

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

#if GMX_SIMD_HAVE_REAL

/*! \brief Test fixture for floating-point tests (identical to the generic \ref SimdTest) */
typedef SimdTest SimdFloatingpointTest;

TEST_F(SimdFloatingpointTest, gmxSimdSetZeroR)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(0.0), simdSetZero());
}

TEST_F(SimdFloatingpointTest, gmxSimdSet1R)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(1.0), simdSet1(1.0));
}

TEST_F(SimdFloatingpointTest, gmxSimdLoad1R)
{
    real r = 2.0;
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(r), simdLoad1(&r));
}

TEST_F(SimdFloatingpointTest, gmxSimdAddR)
{
    GMX_EXPECT_SIMD_REAL_EQ(rSimd_5_7_9,
                            simdAdd(rSimd_1_2_3, rSimd_4_5_6)); // 1+4=5, 2+5=7, 3+6=9
}

TEST_F(SimdFloatingpointTest, gmxSimdSubR)
{
    GMX_EXPECT_SIMD_REAL_EQ(rSimd_4_5_6,
                            simdSub(rSimd_5_7_9, rSimd_1_2_3)); // 5-1=4, 7-2=5, 9-3=6
}

TEST_F(SimdFloatingpointTest, gmxSimdMulR)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(4, 10, 18),
                            simdMul(rSimd_1_2_3, rSimd_4_5_6));
}

TEST_F(SimdFloatingpointTest, gmxSimdFmaddR)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(11, 18, 27),
                            simdFmadd(rSimd_1_2_3, rSimd_4_5_6, rSimd_7_8_9)); // 1*4+7, etc.
}

TEST_F(SimdFloatingpointTest, gmxSimdFmsubR)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(-3, 2, 9),
                            simdFmsub(rSimd_1_2_3, rSimd_4_5_6, rSimd_7_8_9)); // 1*4-7, etc.
}

TEST_F(SimdFloatingpointTest, gmxSimdFnmaddR)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(3, -2, -9),
                            simdFnmadd(rSimd_1_2_3, rSimd_4_5_6, rSimd_7_8_9)); // -1*4+7, etc.
}

TEST_F(SimdFloatingpointTest, gmxSimdFnmsubR)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(-11, -18, -27),
                            simdFnmsub(rSimd_1_2_3, rSimd_4_5_6, rSimd_7_8_9)); // -1*4-7, etc.
}

TEST_F(SimdFloatingpointTest, gmxSimdFabsR)
{
    GMX_EXPECT_SIMD_REAL_EQ(rSimd_1_2_3, simdAbs(rSimd_1_2_3));    // fabs(x)=x
    GMX_EXPECT_SIMD_REAL_EQ(rSimd_1_2_3, simdAbs(rSimd_m1_m2_m3)); // fabs(-x)=x
}

TEST_F(SimdFloatingpointTest, gmxSimdFnegR)
{
    GMX_EXPECT_SIMD_REAL_EQ(rSimd_m1_m2_m3, simdNeg(rSimd_1_2_3));    // fneg(x)=-x
    GMX_EXPECT_SIMD_REAL_EQ(rSimd_1_2_3,    simdNeg(rSimd_m1_m2_m3)); // fneg(-x)=x
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
TEST_F(SimdFloatingpointTest, gmxSimdAndR)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(1.26666259765625),
                            simdAnd(simdSet1(1.3333282470703125),
                                    simdSet1(1.79998779296875)));
}

TEST_F(SimdFloatingpointTest, gmxSimdOrR)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(1.8666534423828125),
                            simdOr(simdSet1(1.3333282470703125),
                                   simdSet1(1.79998779296875)));
}

TEST_F(SimdFloatingpointTest, gmxSimdXorR)
{
    /* Test xor by taking xor with a number and its negative. This should result
     * in only the sign bit being set. We then use this bit change the sign of
     * different numbers.
     */
    SimdReal signbit = simdXor(simdSet1(1.5), simdSet1(-1.5));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(-1, 2, -3), simdXor(signbit, setSimdRealFrom3R(1, -2, 3)));
}

TEST_F(SimdFloatingpointTest, gmxSimdAndnotR)
{
    /* Use xor (which we already tested, so fix that first if both tests fail)
     * to extract the sign bit, and then use andnot to take absolute values.
     */
    SimdReal signbit = simdXor(simdSet1(1.5), simdSet1(-1.5));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(1, 2, 3), simdAndNot(signbit, setSimdRealFrom3R(-1, 2, -3)));
}

#endif

TEST_F(SimdFloatingpointTest, gmxSimdMaxR)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(3, 2, 4), simdMax(rSimd_1_2_3, rSimd_3_1_4));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(3, 2, 4), simdMax(rSimd_3_1_4, rSimd_1_2_3));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(-1, -1, -3), simdMax(rSimd_m1_m2_m3, rSimd_m3_m1_m4));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(-1, -1, -3), simdMax(rSimd_m3_m1_m4, rSimd_m1_m2_m3));
}

TEST_F(SimdFloatingpointTest, gmxSimdMinR)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(1, 1, 3), simdMin(rSimd_1_2_3, rSimd_3_1_4));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(1, 1, 3), simdMin(rSimd_3_1_4, rSimd_1_2_3));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(-3, -2, -4), simdMin(rSimd_m1_m2_m3, rSimd_m3_m1_m4));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(-3, -2, -4), simdMin(rSimd_m3_m1_m4, rSimd_m1_m2_m3));
}

TEST_F(SimdFloatingpointTest, gmxSimdRoundR)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(2), simdRound(simdSet1(2.25)));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(4), simdRound(simdSet1(3.75)));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(-2), simdRound(simdSet1(-2.25)));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(-4), simdRound(simdSet1(-3.75)));
}

TEST_F(SimdFloatingpointTest, gmxSimdTruncR)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(2), simdTrunc(rSimd_2p25));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(3), simdTrunc(rSimd_3p75));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(-2), simdTrunc(rSimd_m2p25));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(-3), simdTrunc(rSimd_m3p75));
}

TEST_F(SimdFloatingpointTest, gmxSimdFractionR)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(0.25), simdFraction(rSimd_2p25));   // fract(2.25)=0.25
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(0.75), simdFraction(rSimd_3p75));   // fract(3.75)=0.75
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(-0.25), simdFraction(rSimd_m2p25)); // fract(-2.25)=-0.25
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(-0.75), simdFraction(rSimd_m3p75)); // fract(-3.75)=-0.75
}

// We explicitly test the exponent/mantissa routines with double precision data,
// since these usually rely on direct manipulation and shift of the SIMD registers,
// where it is easy to make mistakes with single vs double precision.

TEST_F(SimdFloatingpointTest, gmxSimdGetExponentR)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(60.0, -41.0, 54.0), simdGetExponent(rSimd_Exp));
#if GMX_SIMD_HAVE_DOUBLE && defined GMX_DOUBLE
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(587.0, -462.0, 672.0), simdGetExponent(rSimd_ExpDouble));
#endif
}

TEST_F(SimdFloatingpointTest, gmxSimdGetMantissaR)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(1.219097320577810839026256,
                                              1.166738027848349235071623,
                                              1.168904015004464724825084), simdGetMantissa(rSimd_Exp));
#if GMX_SIMD_HAVE_DOUBLE && defined GMX_DOUBLE
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(1.241261238952345623563251,
                                              1.047294723759123852359232,
                                              1.856066204750275957395734), simdGetMantissa(rSimd_ExpDouble));
#endif
}

TEST_F(SimdFloatingpointTest, gmxSimdSetExponentR)
{
    SimdReal x0 = setSimdRealFrom3R(0.5, 11.5, 99.5);
    SimdReal x1 = setSimdRealFrom3R(-0.5, -11.5, -99.5);

    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(pow(2.0, 60.0), pow(2.0, -41.0), pow(2.0, 54.0)),
                            simdSetExponent(setSimdRealFrom3R(60.0, -41.0, 54.0)));
#if GMX_SIMD_HAVE_DOUBLE && defined GMX_DOUBLE
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(pow(2.0, 587.0), pow(2.0, -462.0), pow(2.0, 672.0)),
                            simdSetExponent(setSimdRealFrom3R(587.0, -462.0, 672.0)));
#endif
    /* Rounding mode in simdSetExponent() must be consistent with simdRound() */
    GMX_EXPECT_SIMD_REAL_EQ(simdSetExponent(simdRound(x0)), simdSetExponent(x0));
    GMX_EXPECT_SIMD_REAL_EQ(simdSetExponent(simdRound(x1)), simdSetExponent(x1));
}

/*
 * We do extensive 1/sqrt(x) and 1/x accuracy testing in the math module, so
 * we just make sure the lookup instructions appear to work here
 */

TEST_F(SimdFloatingpointTest, gmxSimdRsqrtR)
{
    SimdReal        x                  = setSimdRealFrom3R(4.0, M_PI, 1234567890.0);
    SimdReal        ref                = setSimdRealFrom3R(0.5, 1.0/sqrt(M_PI), 1.0/sqrt(1234567890.0));
    int             shiftbits          = std::numeric_limits<real>::digits-GMX_SIMD_RSQRT_BITS;

    if (shiftbits < 0)
    {
        shiftbits = 0;
    }

    /* Set the allowed ulp error as 2 to the power of the number of bits in
     * the mantissa that do not have to be correct after the table lookup.
     */
    setUlpTol(1LL << shiftbits);
    GMX_EXPECT_SIMD_REAL_NEAR(ref, simdRsqrt(x));
}

TEST_F(SimdFloatingpointTest, gmxSimdRcpR)
{
    SimdReal        x                  = setSimdRealFrom3R(4.0, M_PI, 1234567890.0);
    SimdReal        ref                = setSimdRealFrom3R(0.25, 1.0/M_PI, 1.0/1234567890.0);
    int             shiftbits          = std::numeric_limits<real>::digits-GMX_SIMD_RCP_BITS;

    if (shiftbits < 0)
    {
        shiftbits = 0;
    }

    /* Set the allowed ulp error as 2 to the power of the number of bits in
     * the mantissa that do not have to be correct after the table lookup.
     */
    setUlpTol(1LL << shiftbits);
    GMX_EXPECT_SIMD_REAL_NEAR(ref, simdRcp(x));
}

TEST_F(SimdFloatingpointTest, gmxSimdBoolCmpEqAndBlendZeroR)
{
    SimdBool eq   = simdCmpEq(rSimd_5_7_9, rSimd_7_8_9);
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(0, 0, 3), simdMask(rSimd_1_2_3, eq));
}

TEST_F(SimdFloatingpointTest, gmxSimdBlendNotZeroR)
{
    SimdBool eq   = simdCmpEq(rSimd_5_7_9, rSimd_7_8_9);
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(1, 2, 0), simdMaskNot(rSimd_1_2_3, eq));
}

TEST_F(SimdFloatingpointTest, gmxSimdBoolCmpLER)
{
    SimdBool le   = simdCmpLe(rSimd_5_7_9, rSimd_7_8_9);
    GMX_EXPECT_SIMD_REAL_EQ(rSimd_1_2_3, simdMask(rSimd_1_2_3, le));
}

TEST_F(SimdFloatingpointTest, gmxSimdBoolCmpLTR)
{
    SimdBool lt   = simdCmpLt(rSimd_5_7_9, rSimd_7_8_9);
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(1, 2, 0), simdMask(rSimd_1_2_3, lt));
}

TEST_F(SimdFloatingpointTest, gmxSimdBoolAndB)
{
    SimdBool eq   = simdCmpEq(rSimd_5_7_9, rSimd_7_8_9);
    SimdBool le   = simdCmpLe(rSimd_5_7_9, rSimd_7_8_9);
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(0, 0, 3), simdMask(rSimd_1_2_3, simdAndB(eq, le)));
}

TEST_F(SimdFloatingpointTest, gmxSimdBoolOrB)
{
    SimdBool eq   = simdCmpEq(rSimd_5_7_9, rSimd_7_8_9);
    SimdBool lt   = simdCmpLt(rSimd_5_7_9, rSimd_7_8_9);
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(1, 2, 3), simdMask(rSimd_1_2_3, simdOrB(eq, lt)));
}

TEST_F(SimdFloatingpointTest, gmxSimdAnytrueB)
{
    SimdBool eq;

    /* this test is a bit tricky since we don't know the simd width.
     * We cannot check for truth values for "any" element beyond the first,
     * since that part of the data will not be used if simd width is 1.
     */
    eq = simdCmpEq(rSimd_5_7_9, setSimdRealFrom3R(5, 0, 0));
    EXPECT_NE(0, simdAnyTrueB(eq));

    eq = simdCmpEq(rSimd_1_2_3, rSimd_4_5_6);
    EXPECT_EQ(0, simdAnyTrueB(eq));
}

TEST_F(SimdFloatingpointTest, gmxSimdBlendvR)
{
    SimdBool lt   = simdCmpLt(rSimd_5_7_9, rSimd_7_8_9);
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(4, 5, 3), simdBlend(rSimd_1_2_3, rSimd_4_5_6, lt));
}

TEST_F(SimdFloatingpointTest, gmxSimdReduceR)
{
    // The horizontal sum of the SIMD variable depends on the width, so
    // simply store it an extra time and calculate what the sum should be
    std::vector<real> v   = simdReal2Vector(rSimd_4_5_6);
    real              sum = 0.0;

    for (int i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        sum += v[i];
    }

    EXPECT_EQ(sum, simdReduce(rSimd_4_5_6));
}

#endif      // GMX_SIMD_HAVE_REAL

/*! \} */
/*! \endcond */

}      // namespace
}      // namespace
}      // namespace

#endif // GMX_SIMD
