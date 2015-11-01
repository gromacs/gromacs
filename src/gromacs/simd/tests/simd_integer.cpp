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

#include "gromacs/simd/simd.h"
#include "gromacs/utility/basedefinitions.h"

#include "simd.h"

/* Some notes on the setup of these tests:
 *
 * It might seem strange to mix different instructions for "setting" SIMD
 * registers, but the difference is that the routines like setSimdIntFrom1I()
 * only use the load/store operations that we already test separately in
 * bootstrap_loadstore.cpp. Since these are "known good" if the bootstrap
 * tests pass, we use them to test the normal SIMD implementation instructions
 * that all have gmx_simd_ prefixes.
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

#if GMX_SIMD_HAVE_INT32

/*! \brief Test fixture for integer tests (identical to the generic \ref SimdTest) */
typedef SimdTest SimdIntegerTest;

TEST_F(SimdIntegerTest, gmxSimdSetZeroI)
{
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(0), simdSetZeroI());
}

TEST_F(SimdIntegerTest, gmxSimdSet1I)
{
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(1), simdSet1I(1));
}

#if GMX_SIMD_HAVE_FINT32_ARITHMETICS
TEST_F(SimdIntegerTest, gmxSimdAddI)
{
    GMX_EXPECT_SIMD_INT_EQ(iSimd_5_7_9,   simdAddI(iSimd_1_2_3, iSimd_4_5_6)    );    // short add
    GMX_EXPECT_SIMD_INT_EQ(iSimd_5M_7M_9M, simdAddI(iSimd_1M_2M_3M, iSimd_4M_5M_6M)); // 32 bit add
}

TEST_F(SimdIntegerTest, gmxSimdSubI)
{
    GMX_EXPECT_SIMD_INT_EQ(iSimd_1_2_3,   simdSubI(iSimd_5_7_9, iSimd_4_5_6)    );    // short sub
    GMX_EXPECT_SIMD_INT_EQ(iSimd_1M_2M_3M, simdSubI(iSimd_5M_7M_9M, iSimd_4M_5M_6M)); // 32 bit sub
}

TEST_F(SimdIntegerTest, gmxSimdMulI)
{
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom3I(4, 10, 18), simdMulI(iSimd_1_2_3, iSimd_4_5_6));           // 2*3=6 (short mul)
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(268435456), simdMulI(simdSet1I(16384), simdSet1I(16384))); // 16384*16384 = 268435456 (long mul)
}
#endif

#if GMX_SIMD_HAVE_FINT32_LOGICAL
TEST_F(SimdIntegerTest, gmxSimdSlliI)
{
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(4194304), simdSlliI(simdSet1I(2), 21)); // 2 << 21 = 4194304
}

TEST_F(SimdIntegerTest, gmxSimdSrliI)
{
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(4), simdSrliI(simdSet1I(4194304), 20)); // 4194304 >> 20 = 4
}

TEST_F(SimdIntegerTest, gmxSimdAndI)
{
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(0xC0C0C0C0), simdAndI(iSimd_0xF0F0F0F0, iSimd_0xCCCCCCCC));
}

TEST_F(SimdIntegerTest, gmxSimdAndnotI)
{
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(0x0C0C0C0C), simdAndNotI(iSimd_0xF0F0F0F0, iSimd_0xCCCCCCCC));
}

TEST_F(SimdIntegerTest, gmxSimdOrI)
{
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(0xFCFCFCFC), simdOrI(iSimd_0xF0F0F0F0, iSimd_0xCCCCCCCC));
}

TEST_F(SimdIntegerTest, gmxSimdXorI)
{
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(0x3C3C3C3C), simdXorI(iSimd_0xF0F0F0F0, iSimd_0xCCCCCCCC));
}
#endif

#if GMX_SIMD_HAVE_INT32_EXTRACT
TEST_F(SimdIntegerTest, gmxSimdExtractI)
{
    GMX_ALIGNED(int, GMX_SIMD_INT32_WIDTH)  idata[GMX_SIMD_INT32_WIDTH];
    SimdInt32 simd;

    for (int i = 0; i < GMX_SIMD_INT32_WIDTH; i++)
    {
        idata[i] = i+1;
    }
    simd = simdLoadI(idata);

    /* We cannot do a loop here, since
     * - C++ gets confused about signed/unsigned if SSE macros are used in EXPECT_EQ()
     * - Extract macros can only take immediates (not variables) on some archs,
     *   and some compilers are not smart enough to expand the for loop.
     *
     * To solve this we use a few values manually instead of a for-loop.
     */
    int extracted_int;
    extracted_int = simdExtractI<0>(simd);
    EXPECT_EQ(1, extracted_int);
#if GMX_SIMD_INT32_WIDTH >= 2
    extracted_int = simdExtractI<1>(simd);
    EXPECT_EQ(2, extracted_int);
#endif
#if GMX_SIMD_INT32_WIDTH >= 4
    extracted_int = simdExtractI<3>(simd);
    EXPECT_EQ(4, extracted_int);
#endif
#if GMX_SIMD_INT32_WIDTH >= 6
    extracted_int = simdExtractI<5>(simd);
    EXPECT_EQ(6, extracted_int);
#endif
#if GMX_SIMD_INT32_WIDTH >= 8
    extracted_int = simdExtractI<7>(simd);
    EXPECT_EQ(8, extracted_int);
#endif
}
#endif

#if GMX_SIMD_HAVE_REAL
TEST_F(SimdIntegerTest, gmxSimdCvtR2I)
{
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(4), simdCvtR2I(rSimd_3p75));
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(-4), simdCvtR2I(rSimd_m3p75));
}

TEST_F(SimdIntegerTest, gmxSimdCvttR2I)
{
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(3), simdCvttR2I(rSimd_3p75));
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom1I(-3), simdCvttR2I(rSimd_m3p75));
}

TEST_F(SimdIntegerTest, gmxSimdCvtI2R)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(2.0), simdCvtI2R(simdSet1I(2)));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom1R(-2.0), simdCvtI2R(simdSet1I(-2)));
}
#endif

#if GMX_SIMD_HAVE_FINT32_ARITHMETICS
TEST_F(SimdIntegerTest, gmxSimdBoolCmpEqAndBlendZeroI)
{
    SimdIBool eq   = simdCmpEqI(iSimd_5_7_9, iSimd_7_8_9);
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom3I(0, 0, 3), simdMaskI(iSimd_1_2_3, eq));
}

TEST_F(SimdIntegerTest, gmxSimdBlendNotZeroI)
{
    SimdIBool eq   = simdCmpEqI(iSimd_5_7_9, iSimd_7_8_9);
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom3I(1, 2, 0), simdMaskNotI(iSimd_1_2_3, eq));
}

TEST_F(SimdIntegerTest, gmxSimdBoolCmpLTI)
{
    SimdIBool lt   = simdCmpLtI(iSimd_5_7_9, iSimd_7_8_9);
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom3I(1, 2, 0), simdMaskI(iSimd_1_2_3, lt));
}

TEST_F(SimdIntegerTest, gmxSimdBoolAndIB)
{
    SimdIBool eq1  = simdCmpEqI(iSimd_5_7_9, iSimd_7_8_9);
    SimdIBool eq2  = simdCmpEqI(iSimd_5_7_9, iSimd_5_7_9);
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom3I(0, 0, 3), simdMaskI(iSimd_1_2_3, simdAndIB(eq1, eq2)));
}

TEST_F(SimdIntegerTest, gmxSimdBoolOrIB)
{
    SimdIBool eq1  = simdCmpEqI(iSimd_5_7_9, iSimd_7_8_9);
    SimdIBool eq2  = simdCmpEqI(iSimd_5_7_9, setSimdIntFrom3I(5, 0, 0));
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom3I(1, 0, 3), simdMaskI(iSimd_1_2_3, simdOrIB(eq1, eq2)));
}

TEST_F(SimdIntegerTest, gmxSimdAnytrueIB)
{
    SimdIBool eq;

    /* See comment in floatingpoint.cpp. We should only check the first element here,
     * since the SIMD width could be 1 as a special case.
     */
    eq = simdCmpEqI(iSimd_5_7_9, setSimdIntFrom3I(5, 0, 0));
    EXPECT_NE(0, simdAnyTrueIB(eq));

    eq = simdCmpEqI(iSimd_1_2_3, iSimd_4_5_6);
    EXPECT_EQ(0, simdAnyTrueIB(eq));
}

TEST_F(SimdIntegerTest, gmxSimdBlendvI)
{
    SimdIBool lt   = simdCmpLtI(iSimd_5_7_9, iSimd_7_8_9);
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom3I(4, 5, 3), simdBlendI(iSimd_1_2_3, iSimd_4_5_6, lt));
}
#endif

#if GMX_SIMD_HAVE_REAL && GMX_SIMD_HAVE_FINT32_ARITHMETICS
TEST_F(SimdIntegerTest, gmxSimdCvtB2IB)
{
    SimdBool  eq   = simdCmpEq(rSimd_5_7_9, setSimdRealFrom3R(5, 0, 0));  // eq should be T,F,F
    SimdIBool eqi  = simdCvtB2IB(eq);
    GMX_EXPECT_SIMD_INT_EQ(setSimdIntFrom3I(1, 0, 0), simdMaskI(iSimd_1_2_3, eqi));

}

TEST_F(SimdIntegerTest, gmxSimdCvtIB2B)
{
    SimdIBool eqi  = simdCmpEqI(iSimd_5_7_9, setSimdIntFrom3I(5, 0, 0));  // eq should be T,F,F
    SimdBool  eq   = simdCvtIB2B(eqi);
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(1.0, 0, 0), simdMask(rSimd_1_2_3, eq));
}
#endif

#endif      // GMX_SIMD_HAVE_INT32

/*! \} */
/*! \endcond */

}      // namespace
}      // namespace
}      // namespace

#endif // GMX_SIMD
