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

#include "testutils/testasserts.h"

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

TEST_F(SimdFloatingpointTest, setZero)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(0.0), setZero());
}

TEST_F(SimdFloatingpointTest, set)
{
    real  r  = 2.0;
    real *p  = &r;
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(1.0), SimdReal(1.0));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(r), SimdReal(*p));
}

TEST_F(SimdFloatingpointTest, add)
{
    GMX_EXPECT_SIMD_REAL_EQ(rSimd_5_7_9,
                            (rSimd_1_2_3 + rSimd_4_5_6)); // 1+4=5, 2+5=7, 3+6=9
}

TEST_F(SimdFloatingpointTest, maskAdd)
{
    SimdBool m = setSimdRealFrom3R(1, 0, 1) != setZero();
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(5, 2, 9),
                            maskAdd(rSimd_1_2_3, rSimd_4_5_6, m));
}

TEST_F(SimdFloatingpointTest, sub)
{
    GMX_EXPECT_SIMD_REAL_EQ(rSimd_4_5_6,
                            (rSimd_5_7_9 - rSimd_1_2_3)); // 5-1=4, 7-2=5, 9-3=6
}

TEST_F(SimdFloatingpointTest, mul)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(4, 10, 18),
                            (rSimd_1_2_3 * rSimd_4_5_6));
}

TEST_F(SimdFloatingpointTest, maskzMul)
{
    SimdBool m = setSimdRealFrom3R(1, 0, 1) != setZero();
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(4, 0, 18),
                            maskzMul(rSimd_1_2_3, rSimd_4_5_6, m));
}

TEST_F(SimdFloatingpointTest, fma)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(11, 18, 27),
                            fma(rSimd_1_2_3, rSimd_4_5_6, rSimd_7_8_9)); // 1*4+7, etc.
}

TEST_F(SimdFloatingpointTest, maskzFma)
{
    SimdBool m = setSimdRealFrom3R(1, 0, 1) != setZero();
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(11, 0, 27),
                            maskzFma(rSimd_1_2_3, rSimd_4_5_6, rSimd_7_8_9, m)); // 1*4+7, etc.
}

TEST_F(SimdFloatingpointTest, fms)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(-3, 2, 9),
                            fms(rSimd_1_2_3, rSimd_4_5_6, rSimd_7_8_9)); // 1*4-7, etc.
}

TEST_F(SimdFloatingpointTest, fnma)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(3, -2, -9),
                            fnma(rSimd_1_2_3, rSimd_4_5_6, rSimd_7_8_9)); // -1*4+7, etc.
}

TEST_F(SimdFloatingpointTest, fnms)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(-11, -18, -27),
                            fnms(rSimd_1_2_3, rSimd_4_5_6, rSimd_7_8_9)); // -1*4-7, etc.
}

TEST_F(SimdFloatingpointTest, abs)
{
    GMX_EXPECT_SIMD_REAL_EQ(rSimd_1_2_3, abs(rSimd_1_2_3));    // fabs(x)=x
    GMX_EXPECT_SIMD_REAL_EQ(rSimd_1_2_3, abs(rSimd_m1_m2_m3)); // fabs(-x)=x
}

TEST_F(SimdFloatingpointTest, neg)
{
    GMX_EXPECT_SIMD_REAL_EQ(rSimd_m1_m2_m3, -(rSimd_1_2_3));    // fneg(x)=-x
    GMX_EXPECT_SIMD_REAL_EQ(rSimd_1_2_3,    -(rSimd_m1_m2_m3)); // fneg(-x)=x
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
TEST_F(SimdFloatingpointTest, and)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(1.26666259765625),
                            (SimdReal(1.3333282470703125) & SimdReal(1.79998779296875)));
}

TEST_F(SimdFloatingpointTest, or)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(1.8666534423828125),
                            (SimdReal(1.3333282470703125) | SimdReal(1.79998779296875)));
}

TEST_F(SimdFloatingpointTest, xor)
{
    /* Test xor by taking xor with a number and its negative. This should result
     * in only the sign bit being set. We then use this bit change the sign of
     * different numbers.
     */
    SimdReal signbit = SimdReal(1.5) ^ SimdReal(-1.5);
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(-1, 2, -3), (signbit ^ setSimdRealFrom3R(1, -2, 3)));
}

TEST_F(SimdFloatingpointTest, andNot)
{
    /* Use xor (which we already tested, so fix that first if both tests fail)
     * to extract the sign bit, and then use andnot to take absolute values.
     */
    SimdReal signbit = SimdReal(1.5) ^ SimdReal(-1.5);
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(1, 2, 3), andNot(signbit, setSimdRealFrom3R(-1, 2, -3)));
}

#endif

TEST_F(SimdFloatingpointTest, max)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(3, 2, 4), max(rSimd_1_2_3, rSimd_3_1_4));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(3, 2, 4), max(rSimd_3_1_4, rSimd_1_2_3));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(-1, -1, -3), max(rSimd_m1_m2_m3, rSimd_m3_m1_m4));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(-1, -1, -3), max(rSimd_m3_m1_m4, rSimd_m1_m2_m3));
}

TEST_F(SimdFloatingpointTest, min)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(1, 1, 3), min(rSimd_1_2_3, rSimd_3_1_4));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(1, 1, 3), min(rSimd_3_1_4, rSimd_1_2_3));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(-3, -2, -4), min(rSimd_m1_m2_m3, rSimd_m3_m1_m4));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(-3, -2, -4), min(rSimd_m3_m1_m4, rSimd_m1_m2_m3));
}

TEST_F(SimdFloatingpointTest, round)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(2), round(SimdReal(2.25)));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(4), round(SimdReal(3.75)));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(-2), round(SimdReal(-2.25)));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(-4), round(SimdReal(-3.75)));
}

TEST_F(SimdFloatingpointTest, trunc)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(2), trunc(rSimd_2p25));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(3), trunc(rSimd_3p75));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(-2), trunc(rSimd_m2p25));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(-3), trunc(rSimd_m3p75));
}

// We explicitly test the exponent/mantissa routines with double precision data,
// since these usually rely on direct manipulation and shift of the SIMD registers,
// where it is easy to make mistakes with single vs double precision.

TEST_F(SimdFloatingpointTest, frexp)
{
    SimdReal  fraction;
    SimdInt32 exponent;

    fraction = frexp(rSimd_Exp, &exponent);

    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(0.609548660288905419513128,
                                              0.5833690139241746175358116,
                                              -0.584452007502232362412542),
                            fraction);
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom3I(61, -40, 55), exponent);


#if GMX_SIMD_HAVE_DOUBLE && GMX_DOUBLE
    fraction = frexp(rSimd_ExpDouble, &exponent);

    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(0.6206306194761728178832527,
                                              0.5236473618795619566768096,
                                              -0.9280331023751380303821179),
                            fraction);
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom3I(588, -461, 673), exponent);
#endif
}

TEST_F(SimdFloatingpointTest, ldexp)
{
    SimdReal x0  = setSimdRealFrom3R(0.5, 11.5, 99.5);
    SimdReal x1  = setSimdRealFrom3R(-0.5, -11.5, -99.5);
    SimdReal one = setSimdRealFrom1R(1.0);

    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(pow(2.0, 60.0), pow(2.0, -41.0), pow(2.0, 54.0)),
                            ldexp(one, setSimdIntFrom3I(60, -41, 54)));
#if GMX_SIMD_HAVE_DOUBLE && GMX_DOUBLE
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(pow(2.0, 587.0), pow(2.0, -462.0), pow(2.0, 672.0)),
                            ldexp(one, setSimdIntFrom3I(587, -462, 672)));
#endif
    /* Rounding mode in conversions must be consistent with simdRound() for SetExponent() to work */
    GMX_EXPECT_SIMD_REAL_EQ(ldexp(one, cvtR2I(round(x0))), ldexp(one, cvtR2I(x0)));
    GMX_EXPECT_SIMD_REAL_EQ(ldexp(one, cvtR2I(round(x1))), ldexp(one, cvtR2I(x1)));
}

/*
 * We do extensive 1/sqrt(x) and 1/x accuracy testing in the math module, so
 * we just make sure the lookup instructions appear to work here
 */

TEST_F(SimdFloatingpointTest, rsqrt)
{
    SimdReal        x                  = setSimdRealFrom3R(4.0, M_PI, 1234567890.0);
    SimdReal        ref                = setSimdRealFrom3R(0.5, 1.0/std::sqrt(M_PI), 1.0/std::sqrt(1234567890.0));
    int             shiftbits          = std::numeric_limits<real>::digits-GMX_SIMD_RSQRT_BITS;

    if (shiftbits < 0)
    {
        shiftbits = 0;
    }

    /* Set the allowed ulp error as 2 to the power of the number of bits in
     * the mantissa that do not have to be correct after the table lookup.
     */
    setUlpTol(1LL << shiftbits);
    GMX_EXPECT_SIMD_REAL_NEAR(ref, rsqrt(x));
}

TEST_F(SimdFloatingpointTest, maskzRsqrt)
{
    SimdReal        x                  = setSimdRealFrom3R(M_PI, -4.0, 0.0);
    // simdCmpLe is tested separately further down
    SimdBool        m                  = setZero() < x;
    SimdReal        ref                = setSimdRealFrom3R(1.0/std::sqrt(M_PI), 0.0, 0.0);
    int             shiftbits          = std::numeric_limits<real>::digits-GMX_SIMD_RSQRT_BITS;

    if (shiftbits < 0)
    {
        shiftbits = 0;
    }

    /* Set the allowed ulp error as 2 to the power of the number of bits in
     * the mantissa that do not have to be correct after the table lookup.
     */
    setUlpTol(1LL << shiftbits);
    GMX_EXPECT_SIMD_REAL_NEAR(ref, maskzRsqrt(x, m));
}

TEST_F(SimdFloatingpointTest, rcp)
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
    GMX_EXPECT_SIMD_REAL_NEAR(ref, rcp(x));
}

TEST_F(SimdFloatingpointTest, maskzRcp)
{
    SimdReal        x                  = setSimdRealFrom3R(M_PI, 0.0, -1234567890.0);
    SimdBool        m                  = (x != setZero());
    SimdReal        ref                = setSimdRealFrom3R(1.0/M_PI, 0.0, -1.0/1234567890.0);
    int             shiftbits          = std::numeric_limits<real>::digits-GMX_SIMD_RCP_BITS;

    if (shiftbits < 0)
    {
        shiftbits = 0;
    }

    /* Set the allowed ulp error as 2 to the power of the number of bits in
     * the mantissa that do not have to be correct after the table lookup.
     */
    setUlpTol(1LL << shiftbits);
    GMX_EXPECT_SIMD_REAL_NEAR(ref, maskzRcp(x, m));
}

TEST_F(SimdFloatingpointTest, cmpEqAndSelectByMask)
{
    SimdBool eq   = rSimd_5_7_9 == rSimd_7_8_9;
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(0, 0, 3), selectByMask(rSimd_1_2_3, eq));
}

TEST_F(SimdFloatingpointTest, selectByNotMask)
{
    SimdBool eq   = rSimd_5_7_9 == rSimd_7_8_9;
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(1, 2, 0), selectByNotMask(rSimd_1_2_3, eq));
}

TEST_F(SimdFloatingpointTest, cmpNe)
{
    SimdBool eq   = rSimd_5_7_9 != rSimd_7_8_9;
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(1, 2, 0), selectByMask(rSimd_1_2_3, eq));
}

TEST_F(SimdFloatingpointTest, cmpLe)
{
    SimdBool le   = rSimd_5_7_9 <= rSimd_7_8_9;
    GMX_EXPECT_SIMD_REAL_EQ(rSimd_1_2_3, selectByMask(rSimd_1_2_3, le));
}

TEST_F(SimdFloatingpointTest, cmpLt)
{
    SimdBool lt   = rSimd_5_7_9 < rSimd_7_8_9;
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(1, 2, 0), selectByMask(rSimd_1_2_3, lt));
}

#if GMX_SIMD_HAVE_INT32_LOGICAL || GMX_SIMD_HAVE_LOGICAL
TEST_F(SimdFloatingpointTest, testBits)
{
    SimdBool eq   = testBits(setSimdRealFrom3R(1, 0, 2));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(1, 0, 3), selectByMask(rSimd_1_2_3, eq));

    // Test if we detect only the sign bit being set
    eq            = testBits(setSimdRealFrom1R(GMX_REAL_NEGZERO));
    GMX_EXPECT_SIMD_REAL_EQ(rSimd_1_2_3, selectByMask(rSimd_1_2_3, eq));
}
#endif

TEST_F(SimdFloatingpointTest, andB)
{
    SimdBool eq   = rSimd_5_7_9 == rSimd_7_8_9;
    SimdBool le   = rSimd_5_7_9 <= rSimd_7_8_9;
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(0, 0, 3), selectByMask(rSimd_1_2_3, (eq && le)));
}

TEST_F(SimdFloatingpointTest, orB)
{
    SimdBool eq   = rSimd_5_7_9 == rSimd_7_8_9;
    SimdBool lt   = rSimd_5_7_9 < rSimd_7_8_9;
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(1, 2, 3), selectByMask(rSimd_1_2_3, (eq || lt)));
}

TEST_F(SimdFloatingpointTest, anyTrueB)
{
    SimdBool eq;

    /* this test is a bit tricky since we don't know the simd width.
     * We cannot check for truth values for "any" element beyond the first,
     * since that part of the data will not be used if simd width is 1.
     */
    eq = rSimd_5_7_9 == setSimdRealFrom3R(5, 0, 0);
    EXPECT_TRUE(anyTrue(eq));

    eq = rSimd_1_2_3 == rSimd_4_5_6;
    EXPECT_FALSE(anyTrue(eq));
}

TEST_F(SimdFloatingpointTest, blend)
{
    SimdBool lt   = rSimd_5_7_9 < rSimd_7_8_9;
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(4, 5, 3), blend(rSimd_1_2_3, rSimd_4_5_6, lt));
}

TEST_F(SimdFloatingpointTest, reduce)
{
    // The horizontal sum of the SIMD variable depends on the width, so
    // simply store it an extra time and calculate what the sum should be
    std::vector<real> v   = simdReal2Vector(rSimd_4_5_6);
    real              sum = 0.0;

    for (int i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        sum += v[i];
    }

    EXPECT_EQ(sum, reduce(rSimd_4_5_6));
}

#endif      // GMX_SIMD_HAVE_REAL

#if GMX_SIMD_HAVE_FLOAT && GMX_SIMD_HAVE_DOUBLE
TEST_F(SimdFloatingpointTest, cvtFloat2Double)
{
    GMX_ALIGNED(float, GMX_SIMD_FLOAT_WIDTH)   f[GMX_SIMD_FLOAT_WIDTH];
    GMX_ALIGNED(double, GMX_SIMD_DOUBLE_WIDTH) d[GMX_SIMD_FLOAT_WIDTH];  // Yes, double array length should be same as float

    int                               i;
    SimdFloat                         vf;
    SimdDouble                        vd0;
    FloatingPointTolerance            tolerance(defaultRealTolerance());

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        f[i] = i;
    }

    vf = load(f);
#if (GMX_SIMD_FLOAT_WIDTH == 2*GMX_SIMD_DOUBLE_WIDTH)
    SimdDouble vd1;
    cvtF2DD(vf, &vd0, &vd1);
    store(d + GMX_SIMD_DOUBLE_WIDTH, vd1); // Store upper part halfway through array
#elif (GMX_SIMD_FLOAT_WIDTH == GMX_SIMD_DOUBLE_WIDTH)
    vd0 = cvtF2D(vf);
#else
#    error Width of float SIMD must either be identical to double, or twice the width.
#endif
    store(d, vd0); // store lower (or whole) part from start of vector

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        EXPECT_REAL_EQ_TOL(f[i], d[i], tolerance);
    }
}

TEST_F(SimdFloatingpointTest, cvtDouble2Float)
{
    GMX_ALIGNED(float, GMX_SIMD_FLOAT_WIDTH)   f[GMX_SIMD_FLOAT_WIDTH];
    GMX_ALIGNED(double, GMX_SIMD_DOUBLE_WIDTH) d[GMX_SIMD_FLOAT_WIDTH];  // Yes, double array length should be same as float
    int                               i;
    SimdFloat                         vf;
    SimdDouble                        vd0;
    FloatingPointTolerance            tolerance(defaultRealTolerance());

    // This fills elements for pd1 too when double width is 2*single width
    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        d[i] = i;
    }

    vd0 = load(d);
#if (GMX_SIMD_FLOAT_WIDTH == 2*GMX_SIMD_DOUBLE_WIDTH)
    SimdDouble vd1 = load(d + GMX_SIMD_DOUBLE_WIDTH); // load upper half of data
    vf = cvtDD2F(vd0, vd1);
#elif (GMX_SIMD_FLOAT_WIDTH == GMX_SIMD_DOUBLE_WIDTH)
    vf = cvtD2F(vd0);
#else
#    error Width of float SIMD must either be identical to double, or twice the width.
#endif
    store(f, vf);

    // This will check elements in pd1 too when double width is 2*single width
    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        EXPECT_REAL_EQ_TOL(d[i], f[i], tolerance);
    }
}
#endif      // GMX_SIMD_HAVE_FLOAT && GMX_SIMD_HAVE_DOUBLE

/*! \} */
/*! \endcond */

}      // namespace
}      // namespace
}      // namespace

#endif // GMX_SIMD
