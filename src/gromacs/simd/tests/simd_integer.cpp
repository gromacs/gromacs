/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016,2017, by the GROMACS development team, led by
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

#include "gromacs/simd/simd.h"
#include "gromacs/utility/basedefinitions.h"

#include "simd.h"

/* Some notes on the setup of these tests:
 *
 * It might seem strange to mix different instructions for "setting" SIMD
 * registers, but the difference is that the routines like setSimdIntFrom1I()
 * only use the load/store operations that we already test separately in
 * bootstrap_loadstore.cpp. Since these are "known good" if the bootstrap
 * tests pass, we use them to test the normal SIMD implementation instructions.
 */

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

/*! \brief Test fixture for integer tests (identical to the generic \ref SimdTest) */
typedef SimdTest SimdIntegerTest;

/* Yes, Virginia. We test for real even for integers. This is because we use
 * the floating-point type when no real integer SIMD type exists (which in turn
 * is because the results of real-to-integer conversions end up there). This
 * means the basic integer SIMD type is available whenever the real one is,
 * but depending on the precision selected that might not be the case.
 *
 * The second we have default-precision floating-point SIMD, we also have
 * the integer SIMD dataype and the most fundamental load/store ops.
 */
#if GMX_SIMD_HAVE_REAL

TEST_F(SimdIntegerTest, setZero)
{
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(0), setZero());
}
TEST_F(SimdIntegerTest, set)
{
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(1), SimdInt32(1));
}
#endif      // GMX_SIMD_HAVE_REAL

#if GMX_SIMD_HAVE_INT32_ARITHMETICS
TEST_F(SimdIntegerTest, add)
{
    GMX_EXPECT_SIMD_INT_EQ(iSimd_5_7_9, iSimd_1_2_3 + iSimd_4_5_6 );         // short add
    GMX_EXPECT_SIMD_INT_EQ(iSimd_5M_7M_9M, iSimd_1M_2M_3M + iSimd_4M_5M_6M); // 32 bit add
}

TEST_F(SimdIntegerTest, sub)
{
    GMX_EXPECT_SIMD_INT_EQ(iSimd_1_2_3, iSimd_5_7_9 - iSimd_4_5_6 );          // short sub
    GMX_EXPECT_SIMD_INT_EQ(iSimd_1M_2M_3M, iSimd_5M_7M_9M - iSimd_4M_5M_6M ); // 32 bit sub
}

TEST_F(SimdIntegerTest, mul)
{
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom3I(4, 10, 18), iSimd_1_2_3 * iSimd_4_5_6);            // 2*3=6 (short mul)
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(268435456), SimdInt32(16384) * SimdInt32(16384) ); // 16384*16384 = 268435456 (long mul)
}

#endif                                                                     // GMX_SIMD_HAVE_INT32_ARITHMETICS

#if GMX_SIMD_HAVE_INT32_LOGICAL
TEST_F(SimdIntegerTest, and)
{
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(0xC0C0C0C0), iSimd_0xF0F0F0F0 & iSimd_0xCCCCCCCC);
}

TEST_F(SimdIntegerTest, andNot)
{
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(0x0C0C0C0C), andNot(iSimd_0xF0F0F0F0, iSimd_0xCCCCCCCC));
}

TEST_F(SimdIntegerTest, or)
{
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(0xFCFCFCFC), iSimd_0xF0F0F0F0 | iSimd_0xCCCCCCCC);
}

TEST_F(SimdIntegerTest, xor)
{
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(0x3C3C3C3C), iSimd_0xF0F0F0F0 ^iSimd_0xCCCCCCCC);
}
#endif      // GMX_SIMD_HAVE_INT32_LOGICAL

#if GMX_SIMD_HAVE_INT32_EXTRACT
TEST_F(SimdIntegerTest, extract)
{
    GMX_ALIGNED(int, GMX_SIMD_REAL_WIDTH)  idata[GMX_SIMD_REAL_WIDTH];
    SimdInt32 simd;

    for (int i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        idata[i] = i+1;
    }
    simd = load<SimdInt32>(idata);

    /* We cannot do a loop here, since
     * - C++ gets confused about signed/unsigned if SSE macros are used in EXPECT_EQ()
     * - Extract macros can only take immediates (not variables) on some archs,
     *   and some compilers are not smart enough to expand the for loop.
     *
     * To solve this we use a few values manually instead of a for-loop.
     */
    int extracted_int;
    extracted_int = extract<0>(simd);
    EXPECT_EQ(1, extracted_int);
#if GMX_SIMD_REAL_WIDTH >= 2
    extracted_int = extract<1>(simd);
    EXPECT_EQ(2, extracted_int);
#endif
#if GMX_SIMD_REAL_WIDTH >= 4
    extracted_int = extract<3>(simd);
    EXPECT_EQ(4, extracted_int);
#endif
#if GMX_SIMD_REAL_WIDTH >= 6
    extracted_int = extract<5>(simd);
    EXPECT_EQ(6, extracted_int);
#endif
#if GMX_SIMD_REAL_WIDTH >= 8
    extracted_int = extract<7>(simd);
    EXPECT_EQ(8, extracted_int);
#endif
}
#endif      // GMX_SIMD_HAVE_INT32_EXTRACT

#if GMX_SIMD_HAVE_REAL
TEST_F(SimdIntegerTest, cvtR2I)
{
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(2), cvtR2I(rSimd_2p25));
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(-2), cvtR2I(rSimd_m2p25));
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(4), cvtR2I(rSimd_3p75));
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(-4), cvtR2I(rSimd_m3p75));
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(3), cvtR2I(rSimd_3p25));
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(-3), cvtR2I(rSimd_m3p25));

    // Test multi-byte numbers
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(123457), cvtR2I(setSimdRealFrom1R(123456.7)));
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(-123457), cvtR2I(setSimdRealFrom1R(-123456.7)));
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(123456), cvtR2I(setSimdRealFrom1R(123456.3)));
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(-123456), cvtR2I(setSimdRealFrom1R(-123456.3)));

#if GMX_DOUBLE
    // Test number with more digits than we can represent in single.
    // Note that our SIMD integers are only 32 bits, so we cannot go beyond that.
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(12345679), cvtR2I(setSimdRealFrom1R(12345678.6)));
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(-12345679), cvtR2I(setSimdRealFrom1R(-12345678.6)));
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(12345678), cvtR2I(setSimdRealFrom1R(12345678.3)));
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(-12345678), cvtR2I(setSimdRealFrom1R(-12345678.3)));
#endif
}

TEST_F(SimdIntegerTest, cvttR2I)
{
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(2), cvttR2I(rSimd_2p25));
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(-2), cvttR2I(rSimd_m2p25));
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(3), cvttR2I(rSimd_3p75));
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(-3), cvttR2I(rSimd_m3p75));
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(3), cvttR2I(rSimd_3p25));
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(-3), cvttR2I(rSimd_m3p25));

    // Test multi-byte numbers
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(123456), cvttR2I(setSimdRealFrom1R(123456.7)));
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(-123456), cvttR2I(setSimdRealFrom1R(-123456.7)));
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(123456), cvttR2I(setSimdRealFrom1R(123456.3)));
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(-123456), cvttR2I(setSimdRealFrom1R(-123456.3)));

#if GMX_DOUBLE
    // Test number with more digits than we can represent in single.
    // Note that our SIMD integers are only 32 bits, so we cannot go beyond that.
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(12345678), cvttR2I(setSimdRealFrom1R(12345678.6)));
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(-12345678), cvttR2I(setSimdRealFrom1R(-12345678.6)));
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(12345678), cvttR2I(setSimdRealFrom1R(12345678.3)));
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(-12345678), cvttR2I(setSimdRealFrom1R(-12345678.3)));
#endif
}

TEST_F(SimdIntegerTest, cvtI2R)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(2.0), cvtI2R(SimdInt32(2)));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(-2.0), cvtI2R(SimdInt32(-2)));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(102448689), cvtI2R(SimdInt32(102448689)));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(-102448689), cvtI2R(SimdInt32(-102448689)));
}
#endif      // GMX_SIMD_HAVE_REAL

#if GMX_SIMD_HAVE_INT32_ARITHMETICS
TEST_F(SimdIntegerTest, cmpEqAndSelectMask)
{
    SimdIBool eq   = (iSimd_5_7_9 == iSimd_7_8_9);
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom3I(0, 0, 3), selectByMask(iSimd_1_2_3, eq));
}

TEST_F(SimdIntegerTest, cmpEqAndSelectNotMask)
{
    SimdIBool eq   = (iSimd_5_7_9 == iSimd_7_8_9);
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom3I(1, 2, 0), selectByNotMask(iSimd_1_2_3, eq));
}

TEST_F(SimdIntegerTest, cmpLt)
{
    SimdIBool lt   = (iSimd_5_7_9 < iSimd_7_8_9);
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom3I(1, 2, 0), selectByMask(iSimd_1_2_3, lt));
}

TEST_F(SimdIntegerTest, testBits)
{
    SimdIBool eq   = testBits(setSimdIntFrom3I(1, 0, 2));
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom3I(1, 0, 3), selectByMask(iSimd_1_2_3, eq));

    // Test if we detect only the sign bit being set
    eq            = testBits(setSimdIntFrom1I(0x80000000));
    GMX_EXPECT_SIMD_INT_EQ(iSimd_1_2_3, selectByMask(iSimd_1_2_3, eq));
}

TEST_F(SimdIntegerTest, andB)
{
    SimdIBool eq1  = (iSimd_5_7_9 == iSimd_7_8_9);
    SimdIBool eq2  = (iSimd_5_7_9 == iSimd_5_7_9);
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom3I(0, 0, 3), selectByMask(iSimd_1_2_3, eq1 && eq2));
}

TEST_F(SimdIntegerTest, orB)
{
    SimdIBool eq1  = (iSimd_5_7_9 == iSimd_7_8_9);
    SimdIBool eq2  = (iSimd_5_7_9 == setSimdIntFrom3I(5, 0, 0));
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom3I(1, 0, 3), selectByMask(iSimd_1_2_3, eq1 || eq2));
}

TEST_F(SimdIntegerTest, anyTrue)
{
    SimdIBool eq;

    /* See comment in floatingpoint.cpp. We should only check the first element here,
     * since the SIMD width could be 1 as a special case.
     */
    eq = (iSimd_5_7_9 == setSimdIntFrom3I(5, 0, 0));
    EXPECT_TRUE(anyTrue(eq));

    eq = (iSimd_1_2_3 == iSimd_4_5_6);
    EXPECT_FALSE(anyTrue(eq));
}

TEST_F(SimdIntegerTest, blend)
{
    SimdIBool lt   = (iSimd_5_7_9 < iSimd_7_8_9);
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom3I(4, 5, 3), blend(iSimd_1_2_3, iSimd_4_5_6, lt));
}
#endif      // GMX_SIMD_HAVE_INT32_ARITHMETICS

#if GMX_SIMD_HAVE_REAL && GMX_SIMD_HAVE_INT32_ARITHMETICS
TEST_F(SimdIntegerTest, cvtB2IB)
{
    SimdBool  eq   = (rSimd_c3c4c5 == rSimd_c3c0c4);  // eq should be T,F,F
    SimdIBool eqi  = cvtB2IB(eq);
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom3I(1, 0, 0), selectByMask(iSimd_1_2_3, eqi));

}

TEST_F(SimdIntegerTest, cvtIB2B)
{
    SimdIBool eqi  = (iSimd_5_7_9 == setSimdIntFrom3I(5, 0, 0));  // eq should be T,F,F
    SimdBool  eq   = cvtIB2B(eqi);
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(c0, 0, 0), selectByMask(rSimd_c0c1c2, eq));
}
#endif      // GMX_SIMD_HAVE_REAL && GMX_SIMD_HAVE_INT32_ARITHMETICS

/*! \} */
/*! \endcond */

}      // namespace
}      // namespace
}      // namespace

#endif // GMX_SIMD
