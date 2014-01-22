/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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

#ifdef GMX_SIMD_HAVE_INT32

TEST(SimdTestBase, gmxSimdSetZeroI)
{
    GMX_ASSERT_SIMD_INT_EQ(0, gmx_simd_setzero_i());
}

TEST(SimdTestBase, gmxSimdSet1I)
{
    GMX_ASSERT_SIMD_INT_EQ(1, gmx_simd_set1_i(1));
}

#ifdef GMX_SIMD_HAVE_FINT32_ARITHMETICS
TEST(SimdTestBase, gmxSimdAddI)
{
    GMX_ASSERT_SIMD_INT_EQ(iSimd_5_7_9,   gmx_simd_add_i(iSimd_1_2_3, iSimd_4_5_6)    );    // short add
    GMX_ASSERT_SIMD_INT_EQ(iSimd_5M_7M_9M, gmx_simd_add_i(iSimd_1M_2M_3M, iSimd_4M_5M_6M)); // 32 bit add
}

TEST(SimdTestBase, gmxSimdSubI)
{
    GMX_ASSERT_SIMD_INT_EQ(iSimd_1_2_3,   gmx_simd_sub_i(iSimd_5_7_9, iSimd_4_5_6)    );    // short sub
    GMX_ASSERT_SIMD_INT_EQ(iSimd_1M_2M_3M, gmx_simd_sub_i(iSimd_5M_7M_9M, iSimd_4M_5M_6M)); // 32 bit sub
}

TEST(SimdTestBase, gmxSimdMulI)
{
    GMX_ASSERT_SIMD_INT_EQ(setSimd3I(4, 10, 18), gmx_simd_mul_i(iSimd_1_2_3, iSimd_4_5_6));            // 2*3=6 (short mul)
    GMX_ASSERT_SIMD_INT_EQ(268435456, gmx_simd_mul_i(gmx_simd_set1_i(16384), gmx_simd_set1_i(16384))); // 16384*16384 = 268435456 (long mul)
}
#endif

#ifdef GMX_SIMD_HAVE_FINT32_LOGICAL
TEST(SimdTestBase, gmxSimdSlliI)
{
    GMX_ASSERT_SIMD_INT_EQ(4194304, gmx_simd_slli_i(gmx_simd_set1_i(2), 21)); // 2 << 21 = 4194304
}

TEST(SimdTestBase, gmxSimdSrliI)
{
    GMX_ASSERT_SIMD_INT_EQ(4, gmx_simd_srli_i(gmx_simd_set1_i(4194304), 20)); // 4194304 >> 20 = 4
}

TEST(SimdTestBase, gmxSimdAndI)
{
    GMX_ASSERT_SIMD_INT_EQ((int)0xC0C0C0C0, gmx_simd_and_i(iSimd_0xF0F0F0F0, iSimd_0xCCCCCCCC));
}

TEST(SimdTestBase, gmxSimdAndnotI)
{
    GMX_ASSERT_SIMD_INT_EQ((int)0x0C0C0C0C, gmx_simd_andnot_i(iSimd_0xF0F0F0F0, iSimd_0xCCCCCCCC));
}

TEST(SimdTestBase, gmxSimdOrI)
{
    GMX_ASSERT_SIMD_INT_EQ((int)0xFCFCFCFC, gmx_simd_or_i(iSimd_0xF0F0F0F0, iSimd_0xCCCCCCCC));
}

TEST(SimdTestBase, gmxSimdXorI)
{
    GMX_ASSERT_SIMD_INT_EQ((int)0x3C3C3C3C, gmx_simd_xor_i(iSimd_0xF0F0F0F0, iSimd_0xCCCCCCCC));
}
#endif

#ifdef GMX_SIMD_HAVE_INT32_EXTRACT
TEST(SimdTestBase, gmxSimdExtractI)
{
    int              idata[GMX_SIMD_INT32_WIDTH*2];
    int *            p = gmx_simd_align_i(idata);
    gmx_simd_int32_t simd;
    int              i, extracted_int;

    for (i = 0; i < GMX_SIMD_INT32_WIDTH; i++)
    {
        p[i] = i+1;
    }
    simd = gmx_simd_load_i(p);

    // 1) C++ gets confused about signed/unsigned if SSE macros is used in ASSERT_EQ()
    // 2) Extract macros can only take immediates (not variables) on some archs,
    //    so we try a few values manually instead of using a loop
    extracted_int = gmx_simd_extract_i(simd, 0);
    ASSERT_EQ(1, extracted_int);
    if (GMX_SIMD_INT32_WIDTH >= 2)
    {
        extracted_int = gmx_simd_extract_i(simd, 1);
        ASSERT_EQ(2, extracted_int);
    }
    if (GMX_SIMD_INT32_WIDTH >= 4)
    {
        extracted_int = gmx_simd_extract_i(simd, 3);
        ASSERT_EQ(4, extracted_int);
    }
    if (GMX_SIMD_INT32_WIDTH >= 6)
    {
        extracted_int = gmx_simd_extract_i(simd, 5);
        ASSERT_EQ(6, extracted_int);
    }
    if (GMX_SIMD_INT32_WIDTH >= 8)
    {
        extracted_int = gmx_simd_extract_i(simd, 7);
        ASSERT_EQ(8, extracted_int);
    }
}
#endif

#ifdef GMX_SIMD_HAVE_REAL
TEST(SimdTestBase, gmxSimdCvtR2I)
{
    GMX_ASSERT_SIMD_INT_EQ( 4, gmx_simd_cvt_r2i(rSimd_3p75));
    GMX_ASSERT_SIMD_INT_EQ(-4, gmx_simd_cvt_r2i(rSimd_m3p75));
}

TEST(SimdTestBase, gmxSimdCvttR2I)
{
    GMX_ASSERT_SIMD_INT_EQ( 3, gmx_simd_cvtt_r2i(rSimd_3p75));
    GMX_ASSERT_SIMD_INT_EQ(-3, gmx_simd_cvtt_r2i(rSimd_m3p75));
}

TEST(SimdTestBase, gmxSimdCvtI2R)
{
    GMX_ASSERT_SIMD_REAL_EQ( 2.0, gmx_simd_cvt_i2r(gmx_simd_set1_i(2)));
    GMX_ASSERT_SIMD_REAL_EQ(-2.0, gmx_simd_cvt_i2r(gmx_simd_set1_i(-2)));
}
#endif

#ifdef GMX_SIMD_HAVE_FINT32_ARITHMETICS
TEST(SimdTestBase, gmxSimdBoolCmpEqAndBlendZeroI)
{
    gmx_simd_ibool_t eq   = gmx_simd_cmpeq_i(iSimd_5_7_9, iSimd_7_8_9);
    GMX_ASSERT_SIMD_INT_EQ(setSimd3I(0, 0, 3), gmx_simd_blendzero_i(iSimd_1_2_3, eq));
}

TEST(SimdTestBase, gmxSimdBoolCmpLTI)
{
    gmx_simd_ibool_t lt   = gmx_simd_cmplt_i(iSimd_5_7_9, iSimd_7_8_9);
    GMX_ASSERT_SIMD_INT_EQ(setSimd3I(1, 2, 0), gmx_simd_blendzero_i(iSimd_1_2_3, lt));
}

TEST(SimdTestBase, gmxSimdBoolAndIB)
{
    gmx_simd_ibool_t eq1  = gmx_simd_cmpeq_i(iSimd_5_7_9, iSimd_7_8_9);
    gmx_simd_ibool_t eq2  = gmx_simd_cmpeq_i(iSimd_5_7_9, iSimd_5_7_9);
    GMX_ASSERT_SIMD_INT_EQ(setSimd3I(0, 0, 3), gmx_simd_blendzero_i(iSimd_1_2_3, gmx_simd_and_ib(eq1, eq2)));
}

TEST(SimdTestBase, gmxSimdBoolOrIB)
{
    gmx_simd_ibool_t eq1  = gmx_simd_cmpeq_i(iSimd_5_7_9, iSimd_7_8_9);
    gmx_simd_ibool_t eq2  = gmx_simd_cmpeq_i(iSimd_5_7_9, setSimd3I(5, 0, 0));
    GMX_ASSERT_SIMD_INT_EQ(setSimd3I(1, 0, 3), gmx_simd_blendzero_i(iSimd_1_2_3, gmx_simd_or_ib(eq1, eq2)));
}

TEST(SimdTestBase, gmxSimdAnytrueIB)
{
    gmx_simd_ibool_t eq;

    // See comment in floatingpoint.cpp. We should only check the first element here.
    eq = gmx_simd_cmpeq_i(iSimd_5_7_9, setSimd3I(5, 0, 0));
    ASSERT_NE(0, gmx_simd_anytrue_ib(eq));

    eq = gmx_simd_cmpeq_i(iSimd_1_2_3, iSimd_4_5_6);
    ASSERT_EQ(0, gmx_simd_anytrue_ib(eq));
}

TEST(SimdTestBase, gmxSimdBlendvI)
{
    gmx_simd_ibool_t lt   = gmx_simd_cmplt_i(iSimd_5_7_9, iSimd_7_8_9);
    GMX_ASSERT_SIMD_INT_EQ(setSimd3I(4, 5, 3), gmx_simd_blendv_i(iSimd_1_2_3, iSimd_4_5_6, lt));
}
#endif

#if (defined GMX_SIMD_HAVE_REAL) && (defined GMX_SIMD_HAVE_FINT32_ARITHMETICS)
TEST(SimdTestBase, gmxSimdCvtB2IB)
{
    gmx_simd_bool_t  eq   = gmx_simd_cmpeq_r(rSimd_5_7_9, setSimd3R(5, 0, 0));  // eq should be T,F,F
    gmx_simd_ibool_t eqi  = gmx_simd_cvt_b2ib(eq);
    GMX_ASSERT_SIMD_INT_EQ(setSimd3I(1, 0, 0), gmx_simd_blendzero_i(iSimd_1_2_3, eqi));

}

TEST(SimdTestBase, gmxSimdCvtIB2B)
{
    gmx_simd_ibool_t eqi  = gmx_simd_cmpeq_i(iSimd_5_7_9, setSimd3I(5, 0, 0));  // eq should be T,F,F
    gmx_simd_bool_t  eq   = gmx_simd_cvt_ib2b(eqi);
    GMX_ASSERT_SIMD_REAL_EQ(setSimd3R(1, 0, 0), gmx_simd_blendzero_r(rSimd_1_2_3, eq));
}
#endif

#endif // GMX_SIMD_HAVE_INT32

}
