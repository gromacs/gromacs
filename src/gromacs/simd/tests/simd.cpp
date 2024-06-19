/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#include "gmxpre.h"

#include "simd.h"

#include <string>

#include "gromacs/simd/simd.h"
#include "gromacs/simd/tests/data.h"
#include "gromacs/utility/basedefinitions.h"

#if GMX_SIMD

namespace gmx
{
namespace test
{

/*! \cond internal */
/*! \addtogroup module_simd */
/*! \{ */

/* Unfortunately we cannot keep static SIMD constants in the test fixture class.
 * The problem is that SIMD memory need to be aligned, and in particular
 * this applies to automatic storage of variables in classes. For SSE registers
 * this means 16-byte alignment (which seems to work), but AVX requires 32-bit
 * alignment. At least both gcc-4.7.3 and Apple clang-5.0 (OS X 10.9) fail to
 * align these variables when they are stored as data in a class.
 *
 * In theory we could set some of these on-the-fly e.g. with setSimdRealFrom3R()
 * instead (although that would mean repeating code between tests), but many of
 * the constants depend on the current precision not to mention they
 * occasionally have many digits that need to be exactly right, and keeping
 * them in a single place makes sure they are consistent.
 */
#    if GMX_SIMD_HAVE_REAL
const SimdReal rSimd_c0c1c2 = setSimdRealFrom3R(c0, c1, c2);
const SimdReal rSimd_c3c4c5 = setSimdRealFrom3R(c3, c4, c5);
const SimdReal rSimd_c6c7c8 = setSimdRealFrom3R(c6, c7, c8);
const SimdReal rSimd_c3c0c4 = setSimdRealFrom3R(c3, c0, c4);
const SimdReal rSimd_c4c6c8 = setSimdRealFrom3R(c4, c6, c8);
const SimdReal rSimd_c7c2c3 = setSimdRealFrom3R(c7, c2, c3);
const SimdReal rSimd_m0m1m2 = setSimdRealFrom3R(-c0, -c1, -c2);
const SimdReal rSimd_m3m0m4 = setSimdRealFrom3R(-c3, -c0, -c4);

const SimdReal rSimd_2p25  = setSimdRealFrom1R(2.25);
const SimdReal rSimd_3p25  = setSimdRealFrom1R(3.25);
const SimdReal rSimd_3p75  = setSimdRealFrom1R(3.75);
const SimdReal rSimd_m2p25 = setSimdRealFrom1R(-2.25);
const SimdReal rSimd_m3p25 = setSimdRealFrom1R(-3.25);
const SimdReal rSimd_m3p75 = setSimdRealFrom1R(-3.75);
const SimdReal rSimd_Exp   = setSimdRealFrom3R(1.4055235171027452623914516e+18,
                                             5.3057102734253445623914516e-13,
                                             -2.1057102745623934534514516e+16);

#        if GMX_SIMD_HAVE_DOUBLE && GMX_DOUBLE
// Make sure we also test exponents outside single precision when we use double
const SimdReal rSimd_ExpDouble1 =
        setSimdRealFrom3R(0.0, 8.794495252903116023030553e-140, -3.637060701570496477655022e+202);
const SimdReal rSimd_ExpDouble2 =
        setSimdRealFrom3R(6.287393598732017379054414e+176, 0.0, -3.637060701570496477655022e+202);
#        endif // GMX_SIMD_HAVE_DOUBLE && GMX_DOUBLE

#        if GMX_SIMD_HAVE_LOGICAL
// The numbers below all have exponent (2^0), which will not change with AND/OR operations.
// We also leave the last part of the mantissa as zeros, to avoid rounding issues in the compiler
#            if GMX_DOUBLE
const SimdReal rSimd_logicalA =
        setSimdRealFrom1R(1.3333333332557231188); // mantissa 01010101010101010101010101010101
const SimdReal rSimd_logicalB =
        setSimdRealFrom1R(1.7999999998137354851); // mantissa 11001100110011001100110011001100
const SimdReal rSimd_logicalResultAnd =
        setSimdRealFrom1R(1.266666666604578495); // mantissa 01000100010001000100010001000100
const SimdReal rSimd_logicalResultOr =
        setSimdRealFrom1R(1.8666666664648801088); // mantissa 11011101110111011101110111011101
#            else                                 // GMX_DOUBLE
const SimdReal rSimd_logicalA = setSimdRealFrom1R(1.3333282470703125); // mantissa 0101010101010101
const SimdReal rSimd_logicalB = setSimdRealFrom1R(1.79998779296875);   // mantissa 1100110011001100
const SimdReal rSimd_logicalResultAnd = setSimdRealFrom1R(1.26666259765625); // mantissa 0100010001000100
const SimdReal rSimd_logicalResultOr = setSimdRealFrom1R(1.8666534423828125); // mantissa 1101110111011101
#            endif                                // GMX_DOUBLE
#        endif                                    // GMX_SIMD_HAVE_LOGICAL

#    endif // GMX_SIMD_HAVE_REAL
#    if GMX_SIMD_HAVE_INT32_ARITHMETICS
const SimdInt32 iSimd_1_2_3    = setSimdIntFrom3I(1, 2, 3);
const SimdInt32 iSimd_4_5_6    = setSimdIntFrom3I(4, 5, 6);
const SimdInt32 iSimd_7_8_9    = setSimdIntFrom3I(7, 8, 9);
const SimdInt32 iSimd_5_7_9    = setSimdIntFrom3I(5, 7, 9);
const SimdInt32 iSimd_1M_2M_3M = setSimdIntFrom3I(1000000, 2000000, 3000000);
const SimdInt32 iSimd_4M_5M_6M = setSimdIntFrom3I(4000000, 5000000, 6000000);
const SimdInt32 iSimd_5M_7M_9M = setSimdIntFrom3I(5000000, 7000000, 9000000);
#    endif
#    if GMX_SIMD_HAVE_INT32_LOGICAL
const SimdInt32 iSimd_0xF0F0F0F0 = setSimdIntFrom1I(0xF0F0F0F0);
const SimdInt32 iSimd_0xCCCCCCCC = setSimdIntFrom1I(0xCCCCCCCC);
#    endif

#    if GMX_SIMD_HAVE_REAL
TEST(SimdTest, GmxAligned)
{
    // Test alignment with two variables that must be aligned, and one that
    // doesn't have to be. The order of variables is up to the compiler, but
    // if it ignores alignment it is highly unlikely that both r1/r3 still end
    // up being aligned by mistake.
    alignas(GMX_SIMD_ALIGNMENT) real r1;
    real                             r2;
    alignas(GMX_SIMD_ALIGNMENT) real r3;

    std::uint64_t addr1 = reinterpret_cast<std::uint64_t>(&r1);
    std::uint64_t addr2 = reinterpret_cast<std::uint64_t>(&r2);
    std::uint64_t addr3 = reinterpret_cast<std::uint64_t>(&r3);

    EXPECT_EQ(0, addr1 % GMX_SIMD_ALIGNMENT);
    EXPECT_NE(0, addr2); // Just so r2 is not optimized away
    EXPECT_EQ(0, addr3 % GMX_SIMD_ALIGNMENT);

    alignas(GMX_SIMD_ALIGNMENT) std::int32_t i1;
    std::int32_t                             i2;
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t i3;

    addr1 = reinterpret_cast<std::uint64_t>(&i1);
    addr2 = reinterpret_cast<std::uint64_t>(&i2);
    addr3 = reinterpret_cast<std::uint64_t>(&i3);

    EXPECT_EQ(0, addr1 % GMX_SIMD_ALIGNMENT);
    EXPECT_NE(0, addr2); // Just so i2 is not optimized away
    EXPECT_EQ(0, addr3 % GMX_SIMD_ALIGNMENT);
}


::std::vector<real> simdReal2Vector(const SimdReal simd)
{
    alignas(GMX_SIMD_ALIGNMENT) real mem[GMX_SIMD_REAL_WIDTH];

    store(mem, simd);
    std::vector<real> v(mem, mem + GMX_SIMD_REAL_WIDTH);

    return v;
}

SimdReal vector2SimdReal(const std::vector<real>& v)
{
    alignas(GMX_SIMD_ALIGNMENT) real mem[GMX_SIMD_REAL_WIDTH];

    for (int i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        mem[i] = v[i % v.size()]; // repeat vector contents to fill simd width
    }
    return load<SimdReal>(mem);
}

SimdReal setSimdRealFrom3R(real r0, real r1, real r2)
{
    std::vector<real> v(3);
    v[0] = r0;
    v[1] = r1;
    v[2] = r2;
    return vector2SimdReal(v);
}

SimdReal setSimdRealFrom1R(real value)
{
    std::vector<real> v(GMX_SIMD_REAL_WIDTH);
    for (int i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        v[i] = value;
    }
    return vector2SimdReal(v);
}

testing::AssertionResult SimdTest::compareSimdRealUlp(const char*    refExpr,
                                                      const char*    tstExpr,
                                                      const SimdReal ref,
                                                      const SimdReal tst)
{
    return compareVectorRealUlp(refExpr, tstExpr, simdReal2Vector(ref), simdReal2Vector(tst));
}

testing::AssertionResult SimdTest::compareSimdEq(const char*    refExpr,
                                                 const char*    tstExpr,
                                                 const SimdReal ref,
                                                 const SimdReal tst)
{
    return compareVectorEq(refExpr, tstExpr, simdReal2Vector(ref), simdReal2Vector(tst));
}

std::vector<std::int32_t> simdInt2Vector(const SimdInt32 simd)
{
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t mem[GMX_SIMD_REAL_WIDTH];

    store(mem, simd);
    std::vector<std::int32_t> v(mem, mem + GMX_SIMD_REAL_WIDTH);

    return v;
}

SimdInt32 vector2SimdInt(const std::vector<std::int32_t>& v)
{
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t mem[GMX_SIMD_REAL_WIDTH];

    for (int i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        mem[i] = v[i % v.size()]; // repeat vector contents to fill simd width
    }
    return load<SimdInt32>(mem);
}

SimdInt32 setSimdIntFrom3I(int i0, int i1, int i2)
{
    std::vector<int> v(3);
    v[0] = i0;
    v[1] = i1;
    v[2] = i2;
    return vector2SimdInt(v);
}

SimdInt32 setSimdIntFrom1I(int value)
{
    std::vector<int> v(GMX_SIMD_REAL_WIDTH);
    for (int i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        v[i] = value;
    }
    return vector2SimdInt(v);
}

::testing::AssertionResult SimdTest::compareSimdEq(const char*     refExpr,
                                                   const char*     tstExpr,
                                                   const SimdInt32 ref,
                                                   const SimdInt32 tst)
{
    return compareVectorEq(refExpr, tstExpr, simdInt2Vector(ref), simdInt2Vector(tst));
}

#    endif // GMX_SIMD_HAVE_REAL

/*! \} */
/*! \endcond */

} // namespace test
} // namespace gmx

#endif // GMX_SIMD
