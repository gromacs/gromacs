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

#include "simd.h"

#include "gromacs/simd/simd.h"
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
#if GMX_SIMD_HAVE_REAL
const SimdReal rSimd_1_2_3    = setSimdRealFrom3R(1, 2, 3);
const SimdReal rSimd_4_5_6    = setSimdRealFrom3R(4, 5, 6);
const SimdReal rSimd_7_8_9    = setSimdRealFrom3R(7, 8, 9);
const SimdReal rSimd_5_7_9    = setSimdRealFrom3R(5, 7, 9);
const SimdReal rSimd_m1_m2_m3 = setSimdRealFrom3R(-1, -2, -3);
const SimdReal rSimd_3_1_4    = setSimdRealFrom3R(3, 1, 4);
const SimdReal rSimd_m3_m1_m4 = setSimdRealFrom3R(-3, -1, -4);
const SimdReal rSimd_2p25     = setSimdRealFrom1R(2.25);
const SimdReal rSimd_3p25     = setSimdRealFrom1R(3.25);
const SimdReal rSimd_3p75     = setSimdRealFrom1R(3.75);
const SimdReal rSimd_m2p25    = setSimdRealFrom1R(-2.25);
const SimdReal rSimd_m3p25    = setSimdRealFrom1R(-3.25);
const SimdReal rSimd_m3p75    = setSimdRealFrom1R(-3.75);
const SimdReal rSimd_Exp      = setSimdRealFrom3R( 1.4055235171027452623914516e+18,
                                                   5.3057102734253445623914516e-13,
                                                   -2.1057102745623934534514516e+16);
#    if GMX_SIMD_HAVE_DOUBLE && GMX_DOUBLE
// Make sure we also test exponents outside single precision when we use double
const SimdReal rSimd_ExpDouble = setSimdRealFrom3R( 6.287393598732017379054414e+176,
                                                    8.794495252903116023030553e-140,
                                                    -3.637060701570496477655022e+202);
#    endif
#endif  // GMX_SIMD_HAVE_REAL
#if GMX_SIMD_HAVE_INT32_ARITHMETICS
const SimdInt32 iSimd_1_2_3      = setSimdIntFrom3I(1, 2, 3);
const SimdInt32 iSimd_4_5_6      = setSimdIntFrom3I(4, 5, 6);
const SimdInt32 iSimd_7_8_9      = setSimdIntFrom3I(7, 8, 9);
const SimdInt32 iSimd_5_7_9      = setSimdIntFrom3I(5, 7, 9);
const SimdInt32 iSimd_1M_2M_3M   = setSimdIntFrom3I(1000000, 2000000, 3000000);
const SimdInt32 iSimd_4M_5M_6M   = setSimdIntFrom3I(4000000, 5000000, 6000000);
const SimdInt32 iSimd_5M_7M_9M   = setSimdIntFrom3I(5000000, 7000000, 9000000);
#endif
#if GMX_SIMD_HAVE_INT32_LOGICAL
const SimdInt32 iSimd_0xF0F0F0F0 = setSimdIntFrom1I(0xF0F0F0F0);
const SimdInt32 iSimd_0xCCCCCCCC = setSimdIntFrom1I(0xCCCCCCCC);
#endif

#if GMX_SIMD_HAVE_REAL
::std::vector<real>
simdReal2Vector(const SimdReal simd)
{
    GMX_ALIGNED(real, GMX_SIMD_REAL_WIDTH)  mem[GMX_SIMD_REAL_WIDTH];

    store(mem, simd);
    std::vector<real>   v(mem, mem+GMX_SIMD_REAL_WIDTH);

    return v;
}

SimdReal
vector2SimdReal(const std::vector<real> &v)
{
    GMX_ALIGNED(real, GMX_SIMD_REAL_WIDTH)  mem[GMX_SIMD_REAL_WIDTH];

    for (int i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        mem[i] = v[i % v.size()];  // repeat vector contents to fill simd width
    }
    return load(mem);
}

SimdReal
setSimdRealFrom3R(real r0, real r1, real r2)
{
    std::vector<real> v(3);
    v[0] = r0;
    v[1] = r1;
    v[2] = r2;
    return vector2SimdReal(v);
}

SimdReal
setSimdRealFrom1R(real value)
{
    std::vector<real> v(GMX_SIMD_REAL_WIDTH);
    for (int i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        v[i] = value;
    }
    return vector2SimdReal(v);
}

testing::AssertionResult
SimdTest::compareSimdRealUlp(const char *  refExpr,     const char *  tstExpr,
                             const SimdReal ref, const SimdReal tst)
{
    return compareVectorRealUlp(refExpr, tstExpr, simdReal2Vector(ref), simdReal2Vector(tst));
}

testing::AssertionResult
SimdTest::compareSimdRealEq(const char * refExpr, const char * tstExpr,
                            const SimdReal ref, const SimdReal tst)
{
    return compareVectorEq(refExpr, tstExpr, simdReal2Vector(ref), simdReal2Vector(tst));
}

std::vector<int>
simdInt2Vector(const SimdInt32 simd)
{
    GMX_ALIGNED(int, GMX_SIMD_REAL_WIDTH)  mem[GMX_SIMD_REAL_WIDTH];

    store(mem, simd);
    std::vector<int>    v(mem, mem+GMX_SIMD_REAL_WIDTH);

    return v;
}

SimdInt32
vector2SimdInt(const std::vector<int> &v)
{
    GMX_ALIGNED(int, GMX_SIMD_REAL_WIDTH)  mem[GMX_SIMD_REAL_WIDTH];

    for (int i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        mem[i] = v[i % v.size()];  // repeat vector contents to fill simd width
    }
    return load(mem);
}

SimdInt32
setSimdIntFrom3I(int i0, int i1, int i2)
{
    std::vector<int> v(3);
    v[0] = i0;
    v[1] = i1;
    v[2] = i2;
    return vector2SimdInt(v);
}

SimdInt32
setSimdIntFrom1I(int value)
{
    std::vector<int> v(GMX_SIMD_REAL_WIDTH);
    for (int i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        v[i] = value;
    }
    return vector2SimdInt(v);
}

::testing::AssertionResult
SimdTest::compareSimdInt32(const char *  refExpr,      const char *  tstExpr,
                           const SimdInt32 ref, const SimdInt32 tst)
{
    return compareVectorEq(refExpr, tstExpr, simdInt2Vector(ref), simdInt2Vector(tst));
}

#endif  // GMX_SIMD_HAVE_REAL

/*! \} */
/*! \endcond */

}      // namespace
}      // namespace

#endif // GMX_SIMD
