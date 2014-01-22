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

#include "simd4_util.h"

namespace simd4Test
{

#ifdef GMX_SIMD4_HAVE_REAL

TEST(Simd4TestBase, gmxSimd4SetZeroR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(0.0, gmx_simd4_setzero_r());
}

TEST(Simd4TestBase, gmxSimd4Set1R)
{
    GMX_EXPECT_SIMD4_REAL_EQ(1.0, gmx_simd4_set1_r(1.0));
}

TEST(Simd4TestBase, gmxSimd4Load1R)
{
    real r = 2.0;
    GMX_EXPECT_SIMD4_REAL_EQ(2.0, gmx_simd4_load1_r(&r));
}

TEST(Simd4TestBase, gmxSimd4AddR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd_5_7_9, gmx_simd4_add_r(rSimd_1_2_3, rSimd_4_5_6)); // 1+4=5, 2+5=7, 3+6=9
}

TEST(Simd4TestBase, gmxSimd4SubR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd_4_5_6, gmx_simd4_sub_r(rSimd_5_7_9, rSimd_1_2_3)); // 5-1=4, 7-2=5, 9-3=6
}

TEST(Simd4TestBase, gmxSimd4MulR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd3R(4, 10, 18), gmx_simd4_mul_r(rSimd_1_2_3, rSimd_4_5_6));
}

TEST(Simd4TestBase, gmxSimd4FmaddR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd3R(11, 18, 27), gmx_simd4_fmadd_r(rSimd_1_2_3, rSimd_4_5_6, rSimd_7_8_9)); // 1*4+7, etc.
}

TEST(Simd4TestBase, gmxSimd4FmsubR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd3R(-3, 2, 9), gmx_simd4_fmsub_r(rSimd_1_2_3, rSimd_4_5_6, rSimd_7_8_9)); // 1*4-7, etc.
}

TEST(Simd4TestBase, gmxSimd4FnmaddR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd3R(3, -2, -9), gmx_simd4_fnmadd_r(rSimd_1_2_3, rSimd_4_5_6, rSimd_7_8_9)); // -1*4+7, etc.
}

TEST(Simd4TestBase, gmxSimd4FnmsubR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd3R(-11, -18, -27), gmx_simd4_fnmsub_r(rSimd_1_2_3, rSimd_4_5_6, rSimd_7_8_9)); // -1*4-7, etc.
}

TEST(Simd4TestBase, gmxSimd4FabsR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd_1_2_3, gmx_simd4_fabs_r(rSimd_1_2_3));    // fabs(x)=x
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd_1_2_3, gmx_simd4_fabs_r(rSimd_m1_m2_m3)); // fabs(-x)=x
}

TEST(Simd4TestBase, gmxSimd4FnegR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd_m1_m2_m3, gmx_simd4_fneg_r(rSimd_1_2_3));   // fneg(x)=-x
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd_1_2_3,   gmx_simd4_fneg_r(rSimd_m1_m2_m3)); // fneg(-x)=x
}

#ifdef GMX_SIMD4_HAVE_LOGICAL
TEST(Simd4TestBase, gmxSimd4AndR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd_Bits3, gmx_simd4_and_r(rSimd_Bits1, rSimd_Bits2)); // Bits1 & Bits2 = Bits3
}

TEST(Simd4TestBase, gmxSimd4AndnotR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd_Bits4, gmx_simd4_andnot_r(rSimd_Bits1, rSimd_Bits2)); // (~Bits1) & Bits2 = Bits3
}

TEST(Simd4TestBase, gmxSimd4OrR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd_Bits5, gmx_simd4_or_r(rSimd_Bits1, rSimd_Bits2)); // Bits1 | Bits2 = Bits3
}

TEST(Simd4TestBase, gmxSimd4XorR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd_Bits6, gmx_simd4_xor_r(rSimd_Bits1, rSimd_Bits2)); // Bits1 ^ Bits2 = Bits3
}
#endif

TEST(Simd4TestBase, gmxSimd4MaxR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd3R(3, 2, 4), gmx_simd4_max_r(rSimd_1_2_3, rSimd_3_1_4));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd3R(3, 2, 4), gmx_simd4_max_r(rSimd_3_1_4, rSimd_1_2_3));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd3R(-1, -1, -3), gmx_simd4_max_r(rSimd_m1_m2_m3, rSimd_m3_m1_m4));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd3R(-1, -1, -3), gmx_simd4_max_r(rSimd_m3_m1_m4, rSimd_m1_m2_m3));
}

TEST(Simd4TestBase, gmxSimd4MinR)
{
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd3R(1, 1, 3), gmx_simd4_min_r(rSimd_1_2_3, rSimd_3_1_4));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd3R(1, 1, 3), gmx_simd4_min_r(rSimd_3_1_4, rSimd_1_2_3));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd3R(-3, -2, -4), gmx_simd4_min_r(rSimd_m1_m2_m3, rSimd_m3_m1_m4));
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd3R(-3, -2, -4), gmx_simd4_min_r(rSimd_m3_m1_m4, rSimd_m1_m2_m3));
}

TEST(Simd4TestBase, gmxSimd4RoundR)
{
    GMX_EXPECT_SIMD4_REAL_EQ( 2, gmx_simd4_round_r(gmx_simd4_set1_r(2.25)));
    GMX_EXPECT_SIMD4_REAL_EQ( 4, gmx_simd4_round_r(gmx_simd4_set1_r(3.75)));
    GMX_EXPECT_SIMD4_REAL_EQ(-2, gmx_simd4_round_r(gmx_simd4_set1_r(-2.25)));
    GMX_EXPECT_SIMD4_REAL_EQ(-4, gmx_simd4_round_r(gmx_simd4_set1_r(-3.75)));
}

TEST(Simd4TestBase, gmxSimd4TruncR)
{
    GMX_EXPECT_SIMD4_REAL_EQ( 2, gmx_simd4_trunc_r(rSimd_2p25));
    GMX_EXPECT_SIMD4_REAL_EQ( 3, gmx_simd4_trunc_r(rSimd_3p75));
    GMX_EXPECT_SIMD4_REAL_EQ(-2, gmx_simd4_trunc_r(rSimd_m2p25));
    GMX_EXPECT_SIMD4_REAL_EQ(-3, gmx_simd4_trunc_r(rSimd_m3p75));
}

// We do extensive 1/sqrt(x) and 1/x accuracy testing in the math module, so
// we just make sure the lookup instructions appear to work here

TEST(Simd4TestBase, gmxSimd4RsqrtR)
{
    gmx_int64_t      ulpTol = 1LL << (mantissaBits-GMX_SIMD_RSQRT_BITS);
    gmx_simd4_real_t x      = setSimd3R(4.0, M_PI, 1234567890.0);
    gmx_simd4_real_t ref    = setSimd3R(0.5, 1.0/sqrt(M_PI), 1.0/sqrt(1234567890.0));

    GMX_EXPECT_SIMD4_REAL_NEAR(ref, gmx_simd4_rsqrt_r(x), ulpTol, 0);
}

TEST(Simd4TestBase, gmxSimd4BoolCmpEqAndBlendZeroR)
{
    gmx_simd4_bool_t eq   = gmx_simd4_cmpeq_r(rSimd_5_7_9, rSimd_7_8_9);
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd3R(0, 0, 3), gmx_simd4_blendzero_r(rSimd_1_2_3, eq));
}

TEST(Simd4TestBase, gmxSimd4BlendNotZeroR)
{
    gmx_simd4_bool_t eq   = gmx_simd4_cmpeq_r(rSimd_5_7_9, rSimd_7_8_9);
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd3R(1, 2, 0), gmx_simd4_blendnotzero_r(rSimd_1_2_3, eq));
}

TEST(Simd4TestBase, gmxSimd4BoolCmpLER)
{
    gmx_simd4_bool_t le   = gmx_simd4_cmple_r(rSimd_5_7_9, rSimd_7_8_9);
    GMX_EXPECT_SIMD4_REAL_EQ(rSimd_1_2_3, gmx_simd4_blendzero_r(rSimd_1_2_3, le));
}

TEST(Simd4TestBase, gmxSimd4BoolCmpLTR)
{
    gmx_simd4_bool_t lt   = gmx_simd4_cmplt_r(rSimd_5_7_9, rSimd_7_8_9);
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd3R(1, 2, 0), gmx_simd4_blendzero_r(rSimd_1_2_3, lt));
}

TEST(Simd4TestBase, gmxSimd4BoolAndB)
{
    gmx_simd4_bool_t eq   = gmx_simd4_cmpeq_r(rSimd_5_7_9, rSimd_7_8_9);
    gmx_simd4_bool_t le   = gmx_simd4_cmple_r(rSimd_5_7_9, rSimd_7_8_9);
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd3R(0, 0, 3), gmx_simd4_blendzero_r(rSimd_1_2_3, gmx_simd4_and_b(eq, le)));
}

TEST(Simd4TestBase, gmxSimd4BoolOrB)
{
    gmx_simd4_bool_t eq   = gmx_simd4_cmpeq_r(rSimd_5_7_9, rSimd_7_8_9);
    gmx_simd4_bool_t lt   = gmx_simd4_cmplt_r(rSimd_5_7_9, rSimd_7_8_9);
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd3R(1, 2, 3), gmx_simd4_blendzero_r(rSimd_1_2_3, gmx_simd4_or_b(eq, lt)));
}

TEST(Simd4TestBase, gmxSimd4AnytrueB)
{
    gmx_simd4_bool_t eq;

    // this test is a bit tricky since we don't know the simd width.
    // We cannot check for truth values for "any" variable with index>1,
    // since that part of the data will not be used if simd width is 1.
    eq = gmx_simd4_cmpeq_r(rSimd_5_7_9, setSimd3R(5, 0, 0));
    ASSERT_NE(0, gmx_simd4_anytrue_b(eq));

    eq = gmx_simd4_cmpeq_r(rSimd_1_2_3, rSimd_4_5_6);
    ASSERT_EQ(0, gmx_simd4_anytrue_b(eq));
}

TEST(Simd4TestBase, gmxSimd4BlendvR)
{
    gmx_simd4_bool_t lt   = gmx_simd4_cmplt_r(rSimd_5_7_9, rSimd_7_8_9);
    GMX_EXPECT_SIMD4_REAL_EQ(setSimd3R(4, 5, 3), gmx_simd4_blendv_r(rSimd_1_2_3, rSimd_4_5_6, lt));
}

TEST(Simd4TestBase, gmxSimd4Dotproduct3R)
{
    gmx_simd4_real_t v1 = setSimd3R(1, 4, 5);
    gmx_simd4_real_t v2 = setSimd3R(3, 8, 2);
#    ifdef GMX_DOUBLE
    ASSERT_DOUBLE_EQ(45.0, gmx_simd4_dotproduct3_r(v1, v2));
#    else
    ASSERT_FLOAT_EQ(45.0, gmx_simd4_dotproduct3_r(v1, v2));
#    endif
}

#endif // GMX_SIMD4_HAVE_REAL

}
