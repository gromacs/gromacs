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

#include "util.h"
#include "smalloc.h"

namespace simdTest
{

/* Unfortunately I had to remove the fixture class here after introducing it,
 * and replace it with static const variables in the simdTest namespace.
 * The problem is that SIMD memory need to be aligned, and in particular
 * this applies to automatic storage of variables. For SSE registers this means
 * 16-byte alignment (which seems to work), but AVX requires 32-bit alignment.
 * At least both gcc-4.7.3 and Apple clang-5.0 (OS X 10.9) fail to align these
 * variables when they are stored as data in a class.
 */
#ifdef GMX_SIMD_HAVE_REAL
const gmx_simd_real_t rSimd_1_2_3    = setSimd3R(1, 2, 3);
const gmx_simd_real_t rSimd_4_5_6    = setSimd3R(4, 5, 6);
const gmx_simd_real_t rSimd_7_8_9    = setSimd3R(7, 8, 9);
const gmx_simd_real_t rSimd_5_7_9    = setSimd3R(5, 7, 9);
const gmx_simd_real_t rSimd_m1_m2_m3 = setSimd3R(-1, -2, -3);
const gmx_simd_real_t rSimd_3_1_4    = setSimd3R(3, 1, 4);
const gmx_simd_real_t rSimd_m3_m1_m4 = setSimd3R(-3, -1, -4);
const gmx_simd_real_t rSimd_2p25     = gmx_simd_set1_r(2.25);
const gmx_simd_real_t rSimd_3p75     = gmx_simd_set1_r(3.75);
const gmx_simd_real_t rSimd_m2p25    = gmx_simd_set1_r(-2.25);
const gmx_simd_real_t rSimd_m3p75    = gmx_simd_set1_r(-3.75);
const gmx_simd_real_t rSimd_Exp      = setSimd3R( 1.4055235171027452623914516e+18,
                                                  5.3057102734253445623914516e-13,
                                                  -2.1057102745623934534514516e+16);
#    ifdef GMX_DOUBLE
// Magic FP numbers corresponding to specific bit patterns
const gmx_simd_real_t rSimd_Bits1    = gmx_simd_set1_r(-1.07730874267432137e+236);
const gmx_simd_real_t rSimd_Bits2    = gmx_simd_set1_r(-9.25596313493178307e+061);
const gmx_simd_real_t rSimd_Bits3    = gmx_simd_set1_r(-8.57750588235293981e+003);
const gmx_simd_real_t rSimd_Bits4    = gmx_simd_set1_r( 1.22416778341839096e-250);
const gmx_simd_real_t rSimd_Bits5    = gmx_simd_set1_r(-1.15711777004554095e+294);
const gmx_simd_real_t rSimd_Bits6    = gmx_simd_set1_r( 1.53063836115600621e-018);
#    else
const gmx_simd_real_t rSimd_Bits1    = gmx_simd_set1_r(-5.9654142337e+29);
const gmx_simd_real_t rSimd_Bits2    = gmx_simd_set1_r(-1.0737417600e+08);
const gmx_simd_real_t rSimd_Bits3    = gmx_simd_set1_r(-6.0235290527e+00);
const gmx_simd_real_t rSimd_Bits4    = gmx_simd_set1_r( 1.0788832913e-31);
const gmx_simd_real_t rSimd_Bits5    = gmx_simd_set1_r(-1.0508719529e+37);
const gmx_simd_real_t rSimd_Bits6    = gmx_simd_set1_r( 1.1488970369e-02);
#    endif
#endif // GMX_SIMD_HAVE_REAL
#ifdef GMX_SIMD_HAVE_INT32
const gmx_simd_int32_t iSimd_1_2_3      = setSimd3I(1, 2, 3);
const gmx_simd_int32_t iSimd_4_5_6      = setSimd3I(4, 5, 6);
const gmx_simd_int32_t iSimd_7_8_9      = setSimd3I(7, 8, 9);
const gmx_simd_int32_t iSimd_5_7_9      = setSimd3I(5, 7, 9);
const gmx_simd_int32_t iSimd_1M_2M_3M   = setSimd3I(1000000, 2000000, 3000000);
const gmx_simd_int32_t iSimd_4M_5M_6M   = setSimd3I(4000000, 5000000, 6000000);
const gmx_simd_int32_t iSimd_5M_7M_9M   = setSimd3I(5000000, 7000000, 9000000);
const gmx_simd_int32_t iSimd_0xF0F0F0F0 = gmx_simd_set1_i(0xF0F0F0F0);
const gmx_simd_int32_t iSimd_0xCCCCCCCC = gmx_simd_set1_i(0xCCCCCCCC);
#endif // GMX_SIMD_HAVE_INT32

#ifdef GMX_SIMD_HAVE_REAL
::std::vector<real>
simdReal2Vector(const gmx_simd_real_t simd)
{
    real                mem[GMX_SIMD_REAL_WIDTH*2];
    real *              p = gmx_simd_align_r(mem);

    gmx_simd_store_r(p, simd);
    std::vector<real>   v(p, p+GMX_SIMD_REAL_WIDTH);

    return v;
}

gmx_simd_real_t
vector2SimdReal(const std::vector<real> &v)
{
    real                mem[GMX_SIMD_REAL_WIDTH*2];
    real *              p = gmx_simd_align_r(mem);

    for (int i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        p[i] = v[i % v.size()];  // repeat vector contents to fill simd width
    }
    return gmx_simd_load_r(p);
}

gmx_simd_real_t
setSimd3R(real r0, real r1, real r2)
{
    std::vector<real> v(3);
    v[0] = r0;
    v[1] = r1;
    v[2] = r2;
    return vector2SimdReal(v);
}

std::string
printSimdReal(const gmx_simd_real_t simd)
{
    std::stringstream   buf;
    std::vector<real>   v = simdTest::simdReal2Vector(simd);
    buf << "[ ";
    buf << std::setprecision(20);
    copy(v.begin(), v.end(), std::ostream_iterator<real>(buf, " "));
    buf << "]";
    return buf.str();
}

testing::AssertionResult
SimdTest::compareSimdRealUlp(const char *  refExpr,     const char *  tstExpr,
                             const gmx_simd_real_t ref, const gmx_simd_real_t tst)
{
    std::vector<real>             vref = simdTest::simdReal2Vector(ref);
    std::vector<real>             vtst = simdTest::simdReal2Vector(tst);
    std::vector<gmx_int64_t>      ulpDiff(GMX_SIMD_REAL_WIDTH);
    bool                          eq     = true;
    bool                          signOk = true;
    int                           i;
#    ifdef GMX_DOUBLE
    union {
        double r; gmx_int64_t i;
    } conv0, conv1;
#    else
    union {
        float  r; int i;
    } conv0, conv1;
#    endif

    for (i = 0; i < GMX_SIMD_REAL_WIDTH && eq == true; i++)
    {
        eq     = eq && ( fabs(vref[i]-vtst[i]) < absTol_ );
        signOk = signOk && ( vref[i]*vtst[i] >= 0 );
    }
    if (eq == true)
    {
        // Pass test if all values were within absolute tolerance
        return ::testing::AssertionSuccess();
    }
    else if (signOk == false)
    {
        return ::testing::AssertionFailure()
               << "Failing comparison between " << refExpr << " and " << tstExpr << std::endl
               << "SIMD contents differs in sign from reference: " << std::endl
               << "Ref values:  " << printSimdReal(ref) << std::endl
               << "SIMD values: " << printSimdReal(tst) << std::endl;
    }

    for (i = 0, eq = true; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        conv0.r    = vref[i];
        conv1.r    = vtst[i];
        ulpDiff[i] = llabs(conv0.i-conv1.i);
        eq         = eq && (ulpDiff[i] <= ulpTol_);
    }

    if (eq == true)
    {
        return ::testing::AssertionSuccess();
    }
    else
    {
        return ::testing::AssertionFailure()
               << "Failing ulp comparison between " << refExpr << " and " << tstExpr << std::endl
               << "Requested ulp tolerance: " << ulpTol_ << std::endl
               << "Requested abs tolerance: " << absTol_ << std::endl
               << "Ref  values: " << printSimdReal(ref) << std::endl
               << "SIMD values: " << printSimdReal(tst) << std::endl
               << "Ulp diff.:   " << ::testing::PrintToString(ulpDiff) << std::endl;
    }
}

testing::AssertionResult
SimdTest::compareSimdRealUlp(const char * refExpr, const char * tstExpr,
                             real ref,             const gmx_simd_real_t tst)
{
    return compareSimdRealUlp(refExpr, tstExpr, gmx_simd_set1_r(ref), tst);
}
#endif // GMX_SIMD_HAVE_REAL

#ifdef GMX_SIMD_HAVE_INT32
std::vector<int>
simdInt2Vector(const gmx_simd_int32_t simd)
{
    int *               p;
    snew_aligned(p, GMX_SIMD_INT32_WIDTH, GMX_SIMD_INT32_WIDTH*sizeof(int));

    gmx_simd_store_i(p, simd);
    std::vector<int>    v(p, p+GMX_SIMD_INT32_WIDTH);
    sfree_aligned(p);

    return v;
}

gmx_simd_int32_t
vector2SimdInt(const std::vector<int> &v)
{
    int                 mem[GMX_SIMD_INT32_WIDTH*2];
    int *               p = gmx_simd_align_i(mem);

    for (int i = 0; i < GMX_SIMD_INT32_WIDTH; i++)
    {
        p[i] = v[i % v.size()];  // repeat vector contents to fill simd width
    }
    return gmx_simd_load_i(p);
}

gmx_simd_int32_t
setSimd3I(int i0, int i1, int i2)
{
    std::vector<int> v(3);
    v[0] = i0;
    v[1] = i1;
    v[2] = i2;
    return vector2SimdInt(v);
}

std::string
printSimdInt32(const gmx_simd_int32_t simd)
{
    std::stringstream   buf;
    std::vector<int>    v = simdTest::simdInt2Vector(simd);
    buf << "[ ";
    copy(v.begin(), v.end(), std::ostream_iterator<int>(buf, " "));
    buf << "]";
    return buf.str();
}

testing::AssertionResult
SimdTest::compareSimdInt32(const char *  refExpr,      const char *  tstExpr,
                           const gmx_simd_int32_t ref, const gmx_simd_int32_t tst)
{
    std::vector<int>    vref = simdTest::simdInt2Vector(ref);
    std::vector<int>    vtst = simdTest::simdInt2Vector(tst);
    bool                eq   = true;
    int                 i;

    for (i = 0; i < GMX_SIMD_INT32_WIDTH && eq == true; i++)
    {
        eq     = eq && ( vref[i] == vtst[i] );
    }
    if (eq == true)
    {
        return ::testing::AssertionSuccess();
    }
    else
    {
        return ::testing::AssertionFailure()
               << "Failing SIMD comparison between " << refExpr << " and " << tstExpr << std::endl
               << "Ref. values: " << printSimdInt32(ref) << std::endl
               << "Test values: " << printSimdInt32(tst) << std::endl;
    }
}

testing::AssertionResult
SimdTest::compareSimdInt32(const char * refExpr, const char * tstExpr,
                           int ref,              const gmx_simd_int32_t tst)
{
    return compareSimdInt32(refExpr, tstExpr, gmx_simd_set1_i(ref), tst);
}
#endif // GMX_SIMD_HAVE_INT32

}      // namespace simdTest
