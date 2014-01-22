/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include "gromacs/math/utilities.h"

#include "util.h"

namespace simdTest
{

// Unfortunately I had to remove the fixture class here after introducing it,
// and replace it with static const variables in the simdTest namespace.
// The problem is that SIMD memory need to be aligned, and in particular
// this applies to automatic storage of variables. For SSE registers this means
// 16-byte alignment (which seems to work), but AVX requires 32-bit alignment.
// At least both gcc-4.7.3 and Apple clang-5.0 (OS X 10.9) fail to align these
// variables when they are stored as data in a class.

#ifdef GMX_SIMD_HAVE_REAL

TEST_F(SimdTest, gmxSimdSetZeroR)
{
    GMX_EXPECT_SIMD_REAL_EQ(0.0, gmx_simd_setzero_r());
}

TEST_F(SimdTest, gmxSimdSet1R)
{
    GMX_EXPECT_SIMD_REAL_EQ(1.0, gmx_simd_set1_r(1.0));
}

TEST_F(SimdTest, gmxSimdLoad1R)
{
    real r = 2.0;
    GMX_EXPECT_SIMD_REAL_EQ(2.0, gmx_simd_load1_r(&r));
}

TEST_F(SimdTest, gmxSimdAddR)
{
    GMX_EXPECT_SIMD_REAL_EQ(rSimd_5_7_9, gmx_simd_add_r(rSimd_1_2_3, rSimd_4_5_6)); // 1+4=5, 2+5=7, 3+6=9
}

TEST_F(SimdTest, gmxSimdSubR)
{
    GMX_EXPECT_SIMD_REAL_EQ(rSimd_4_5_6, gmx_simd_sub_r(rSimd_5_7_9, rSimd_1_2_3)); // 5-1=4, 7-2=5, 9-3=6
}

TEST_F(SimdTest, gmxSimdMulR)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimd3R(4, 10, 18), gmx_simd_mul_r(rSimd_1_2_3, rSimd_4_5_6));
}

TEST_F(SimdTest, gmxSimdFmaddR)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimd3R(11, 18, 27), gmx_simd_fmadd_r(rSimd_1_2_3, rSimd_4_5_6, rSimd_7_8_9)); // 1*4+7, etc.
}

TEST_F(SimdTest, gmxSimdFmsubR)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimd3R(-3, 2, 9), gmx_simd_fmsub_r(rSimd_1_2_3, rSimd_4_5_6, rSimd_7_8_9)); // 1*4-7, etc.
}

TEST_F(SimdTest, gmxSimdFnmaddR)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimd3R(3, -2, -9), gmx_simd_fnmadd_r(rSimd_1_2_3, rSimd_4_5_6, rSimd_7_8_9)); // -1*4+7, etc.
}

TEST_F(SimdTest, gmxSimdFnmsubR)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimd3R(-11, -18, -27), gmx_simd_fnmsub_r(rSimd_1_2_3, rSimd_4_5_6, rSimd_7_8_9)); // -1*4-7, etc.
}

TEST_F(SimdTest, gmxSimdFabsR)
{
    GMX_EXPECT_SIMD_REAL_EQ(rSimd_1_2_3, gmx_simd_fabs_r(rSimd_1_2_3));    // fabs(x)=x
    GMX_EXPECT_SIMD_REAL_EQ(rSimd_1_2_3, gmx_simd_fabs_r(rSimd_m1_m2_m3)); // fabs(-x)=x
}

TEST_F(SimdTest, gmxSimdFnegR)
{
    GMX_EXPECT_SIMD_REAL_EQ(rSimd_m1_m2_m3, gmx_simd_fneg_r(rSimd_1_2_3));   // fneg(x)=-x
    GMX_EXPECT_SIMD_REAL_EQ(rSimd_1_2_3,   gmx_simd_fneg_r(rSimd_m1_m2_m3)); // fneg(-x)=x
}

#ifdef GMX_SIMD_HAVE_LOGICAL
TEST_F(SimdTest, gmxSimdAndR)
{
    GMX_EXPECT_SIMD_REAL_EQ(rSimd_Bits3, gmx_simd_and_r(rSimd_Bits1, rSimd_Bits2)); // Bits1 & Bits2 = Bits3
}

TEST_F(SimdTest, gmxSimdAndnotR)
{
    GMX_EXPECT_SIMD_REAL_EQ(rSimd_Bits4, gmx_simd_andnot_r(rSimd_Bits1, rSimd_Bits2)); // (~Bits1) & Bits2 = Bits3
}

TEST_F(SimdTest, gmxSimdOrR)
{
    GMX_EXPECT_SIMD_REAL_EQ(rSimd_Bits5, gmx_simd_or_r(rSimd_Bits1, rSimd_Bits2)); // Bits1 | Bits2 = Bits3
}

TEST_F(SimdTest, gmxSimdXorR)
{
    GMX_EXPECT_SIMD_REAL_EQ(rSimd_Bits6, gmx_simd_xor_r(rSimd_Bits1, rSimd_Bits2)); // Bits1 ^ Bits2 = Bits3
}
#endif

TEST_F(SimdTest, gmxSimdMaxR)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimd3R(3, 2, 4), gmx_simd_max_r(rSimd_1_2_3, rSimd_3_1_4));
    GMX_EXPECT_SIMD_REAL_EQ(setSimd3R(3, 2, 4), gmx_simd_max_r(rSimd_3_1_4, rSimd_1_2_3));
    GMX_EXPECT_SIMD_REAL_EQ(setSimd3R(-1, -1, -3), gmx_simd_max_r(rSimd_m1_m2_m3, rSimd_m3_m1_m4));
    GMX_EXPECT_SIMD_REAL_EQ(setSimd3R(-1, -1, -3), gmx_simd_max_r(rSimd_m3_m1_m4, rSimd_m1_m2_m3));
}

TEST_F(SimdTest, gmxSimdMinR)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimd3R(1, 1, 3), gmx_simd_min_r(rSimd_1_2_3, rSimd_3_1_4));
    GMX_EXPECT_SIMD_REAL_EQ(setSimd3R(1, 1, 3), gmx_simd_min_r(rSimd_3_1_4, rSimd_1_2_3));
    GMX_EXPECT_SIMD_REAL_EQ(setSimd3R(-3, -2, -4), gmx_simd_min_r(rSimd_m1_m2_m3, rSimd_m3_m1_m4));
    GMX_EXPECT_SIMD_REAL_EQ(setSimd3R(-3, -2, -4), gmx_simd_min_r(rSimd_m3_m1_m4, rSimd_m1_m2_m3));
}

TEST_F(SimdTest, gmxSimdRoundR)
{
    GMX_EXPECT_SIMD_REAL_EQ( 2, gmx_simd_round_r(gmx_simd_set1_r(2.25)));
    GMX_EXPECT_SIMD_REAL_EQ( 4, gmx_simd_round_r(gmx_simd_set1_r(3.75)));
    GMX_EXPECT_SIMD_REAL_EQ(-2, gmx_simd_round_r(gmx_simd_set1_r(-2.25)));
    GMX_EXPECT_SIMD_REAL_EQ(-4, gmx_simd_round_r(gmx_simd_set1_r(-3.75)));
}

TEST_F(SimdTest, gmxSimdTruncR)
{
    GMX_EXPECT_SIMD_REAL_EQ( 2, gmx_simd_trunc_r(rSimd_2p25));
    GMX_EXPECT_SIMD_REAL_EQ( 3, gmx_simd_trunc_r(rSimd_3p75));
    GMX_EXPECT_SIMD_REAL_EQ(-2, gmx_simd_trunc_r(rSimd_m2p25));
    GMX_EXPECT_SIMD_REAL_EQ(-3, gmx_simd_trunc_r(rSimd_m3p75));
}

TEST_F(SimdTest, gmxSimdFractionR)
{
    GMX_EXPECT_SIMD_REAL_EQ( 0.25, gmx_simd_fraction_r(rSimd_2p25));  // fract(2.25)=0.25
    GMX_EXPECT_SIMD_REAL_EQ( 0.75, gmx_simd_fraction_r(rSimd_3p75));  // fract(3.75)=0.75
    GMX_EXPECT_SIMD_REAL_EQ(-0.25, gmx_simd_fraction_r(rSimd_m2p25)); // fract(-2.25)=-0.25
    GMX_EXPECT_SIMD_REAL_EQ(-0.75, gmx_simd_fraction_r(rSimd_m3p75)); // fract(-3.75)=-0.75
}

TEST_F(SimdTest, gmxSimdGetExponentR)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimd3R(60.0, -41.0, 54.0), gmx_simd_get_exponent_r(rSimd_Exp));
}

TEST_F(SimdTest, gmxSimdGetMantissaR)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimd3R(1.219097320577810839026256,
                                      1.166738027848349235071623,
                                      1.168904015004464724825084), gmx_simd_get_mantissa_r(rSimd_Exp));
}

TEST_F(SimdTest, gmxSimdSetExponentR)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimd3R(pow(2.0, 60.0), pow(2.0, -41.0), pow(2.0, 54.0)),
                            gmx_simd_set_exponent_r(setSimd3R(60.0, -41.0, 54.0)));
}

// We do extensive 1/sqrt(x) and 1/x accuracy testing in the math module, so
// we just make sure the lookup instructions appear to work here

TEST_F(SimdTest, gmxSimdRsqrtR)
{
    gmx_simd_real_t x      = setSimd3R(4.0, M_PI, 1234567890.0);
    gmx_simd_real_t ref    = setSimd3R(0.5, 1.0/sqrt(M_PI), 1.0/sqrt(1234567890.0));

    setUlpTol(1LL << (mantissaBits-GMX_SIMD_RSQRT_BITS));

    GMX_EXPECT_SIMD_REAL_NEAR(ref, gmx_simd_rsqrt_r(x));
}

TEST_F(SimdTest, gmxSimdRcpR)
{
    gmx_simd_real_t x      = setSimd3R(4.0, M_PI, 1234567890.0);
    gmx_simd_real_t ref    = setSimd3R(0.25, 1.0/M_PI, 1.0/1234567890.0);

    setUlpTol(1LL << (mantissaBits-GMX_SIMD_RCP_BITS));

    GMX_EXPECT_SIMD_REAL_NEAR(ref, gmx_simd_rcp_r(x));
}

TEST_F(SimdTest, gmxSimdBoolCmpEqAndBlendZeroR)
{
    gmx_simd_bool_t eq   = gmx_simd_cmpeq_r(rSimd_5_7_9, rSimd_7_8_9);
    GMX_EXPECT_SIMD_REAL_EQ(setSimd3R(0, 0, 3), gmx_simd_blendzero_r(rSimd_1_2_3, eq));
}

TEST_F(SimdTest, gmxSimdBlendNotZeroR)
{
    gmx_simd_bool_t eq   = gmx_simd_cmpeq_r(rSimd_5_7_9, rSimd_7_8_9);
    GMX_EXPECT_SIMD_REAL_EQ(setSimd3R(1, 2, 0), gmx_simd_blendnotzero_r(rSimd_1_2_3, eq));
}

TEST_F(SimdTest, gmxSimdBoolCmpLER)
{
    gmx_simd_bool_t le   = gmx_simd_cmple_r(rSimd_5_7_9, rSimd_7_8_9);
    GMX_EXPECT_SIMD_REAL_EQ(rSimd_1_2_3, gmx_simd_blendzero_r(rSimd_1_2_3, le));
}

TEST_F(SimdTest, gmxSimdBoolCmpLTR)
{
    gmx_simd_bool_t lt   = gmx_simd_cmplt_r(rSimd_5_7_9, rSimd_7_8_9);
    GMX_EXPECT_SIMD_REAL_EQ(setSimd3R(1, 2, 0), gmx_simd_blendzero_r(rSimd_1_2_3, lt));
}

TEST_F(SimdTest, gmxSimdBoolAndB)
{
    gmx_simd_bool_t eq   = gmx_simd_cmpeq_r(rSimd_5_7_9, rSimd_7_8_9);
    gmx_simd_bool_t le   = gmx_simd_cmple_r(rSimd_5_7_9, rSimd_7_8_9);
    GMX_EXPECT_SIMD_REAL_EQ(setSimd3R(0, 0, 3), gmx_simd_blendzero_r(rSimd_1_2_3, gmx_simd_and_b(eq, le)));
}

TEST_F(SimdTest, gmxSimdBoolOrB)
{
    gmx_simd_bool_t eq   = gmx_simd_cmpeq_r(rSimd_5_7_9, rSimd_7_8_9);
    gmx_simd_bool_t lt   = gmx_simd_cmplt_r(rSimd_5_7_9, rSimd_7_8_9);
    GMX_EXPECT_SIMD_REAL_EQ(setSimd3R(1, 2, 3), gmx_simd_blendzero_r(rSimd_1_2_3, gmx_simd_or_b(eq, lt)));
}

TEST_F(SimdTest, gmxSimdAnytrueB)
{
    gmx_simd_bool_t eq;

    // this test is a bit tricky since we don't know the simd width.
    // We cannot check for truth values for "any" variable with index>1,
    // since that part of the data will not be used if simd width is 1.
    eq = gmx_simd_cmpeq_r(rSimd_5_7_9, setSimd3R(5, 0, 0));
    EXPECT_NE(0, gmx_simd_anytrue_b(eq));

    eq = gmx_simd_cmpeq_r(rSimd_1_2_3, rSimd_4_5_6);
    EXPECT_EQ(0, gmx_simd_anytrue_b(eq));
}

TEST_F(SimdTest, gmxSimdBlendvR)
{
    gmx_simd_bool_t lt   = gmx_simd_cmplt_r(rSimd_5_7_9, rSimd_7_8_9);
    GMX_EXPECT_SIMD_REAL_EQ(setSimd3R(4, 5, 3), gmx_simd_blendv_r(rSimd_1_2_3, rSimd_4_5_6, lt));
}

#endif // GMX_SIMD_HAVE_REAL

}
