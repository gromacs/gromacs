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

/*! \internal \file
 * \brief
 * Implements test fixture for SimdTests
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 * \ingroup module_simd
 */

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
const gmx_simd_real_t rSimd_1_2_3    = setSimd3R(1, 2, 3);    //!< Generic (different) fp values
const gmx_simd_real_t rSimd_4_5_6    = setSimd3R(4, 5, 6);    //!< Generic (different) fp values
const gmx_simd_real_t rSimd_7_8_9    = setSimd3R(7, 8, 9);    //!< Generic (different) fp values
const gmx_simd_real_t rSimd_5_7_9    = setSimd3R(5, 7, 9);    //!< rSimd_1_2_3 + rSimd_4_5_6
const gmx_simd_real_t rSimd_m1_m2_m3 = setSimd3R(-1, -2, -3); //!< Generic negative fp values
const gmx_simd_real_t rSimd_3_1_4    = setSimd3R(3, 1, 4);    //!< Used to test min/max
const gmx_simd_real_t rSimd_m3_m1_m4 = setSimd3R(-3, -1, -4); //!< negative rSimd_3_1_4
const gmx_simd_real_t rSimd_2p25     = gmx_simd_set1_r(2.25); //!< Rounds down
const gmx_simd_real_t rSimd_3p75     = gmx_simd_set1_r(3.75); //!< Rounds up
const gmx_simd_real_t rSimd_m2p25    = gmx_simd_set1_r(-2.25);//!< Negative rounds up
const gmx_simd_real_t rSimd_m3p75    = gmx_simd_set1_r(-3.75);//!< Negative rounds down
//! Three large floating-point values whose exponents are >32.
const gmx_simd_real_t rSimd_Exp      = setSimd3R( 1.4055235171027452623914516e+18,
                                                  5.3057102734253445623914516e-13,
                                                  -2.1057102745623934534514516e+16);
#    ifdef GMX_DOUBLE
// Magic FP numbers corresponding to specific bit patterns
const gmx_simd_real_t rSimd_Bits1    = gmx_simd_set1_r(-1.07730874267432137e+236); //!< 0xF0F0F0F0F0F0F0F0
const gmx_simd_real_t rSimd_Bits2    = gmx_simd_set1_r(-9.25596313493178307e+061); //!< 0xCCCCCCCCCCCCCCCC
const gmx_simd_real_t rSimd_Bits3    = gmx_simd_set1_r(-8.57750588235293981e+003); //!< 0xC0C0C0C0C0C0C0C0
const gmx_simd_real_t rSimd_Bits4    = gmx_simd_set1_r( 1.22416778341839096e-250); //!< 0x0C0C0C0C0C0C0C0C
const gmx_simd_real_t rSimd_Bits5    = gmx_simd_set1_r(-1.15711777004554095e+294); //!< 0xFCFCFCFCFCFCFCFC
const gmx_simd_real_t rSimd_Bits6    = gmx_simd_set1_r( 1.53063836115600621e-018); //!< 0x3C3C3C3C3C3C3C3C
#    else
const gmx_simd_real_t rSimd_Bits1    = gmx_simd_set1_r(-5.9654142337e+29);         //!< 0xF0F0F0F0
const gmx_simd_real_t rSimd_Bits2    = gmx_simd_set1_r(-1.0737417600e+08);         //!< 0xCCCCCCCC
const gmx_simd_real_t rSimd_Bits3    = gmx_simd_set1_r(-6.0235290527e+00);         //!< 0xC0C0C0C0
const gmx_simd_real_t rSimd_Bits4    = gmx_simd_set1_r( 1.0788832913e-31);         //!< 0x0C0C0C0C
const gmx_simd_real_t rSimd_Bits5    = gmx_simd_set1_r(-1.0508719529e+37);         //!< 0xFCFCFCFC
const gmx_simd_real_t rSimd_Bits6    = gmx_simd_set1_r( 1.1488970369e-02);         //!< 0x3C3C3C3C
#    endif
#    ifdef GMX_DOUBLE
const int mantissaBits = 52;  //!< Size of double mantissa
const int defaultTol   = 255; //!< Aim for ~twice the precision we have in single
#    else
const int mantissaBits = 23;  //!< Size of double mantissa
const int defaultTol   = 10;  //!< Be a bit liberal so compiler optimization doesn't bite us
#    endif
#endif
#ifdef GMX_SIMD_HAVE_INT32
const gmx_simd_int32_t iSimd_1_2_3      = setSimd3I(1, 2, 3); //!< Three generic ints
const gmx_simd_int32_t iSimd_4_5_6      = setSimd3I(4, 5, 6); //!< Three generic ints
const gmx_simd_int32_t iSimd_7_8_9      = setSimd3I(7, 8, 9); //!< Three generic ints
const gmx_simd_int32_t iSimd_5_7_9      = setSimd3I(5, 7, 9); //!< iSimd_1_2_3 + iSimd_4_5_6
const gmx_simd_int32_t iSimd_1M_2M_3M   = setSimd3I(1000000, 2000000, 3000000); //!< Term for 32bit add/sub
const gmx_simd_int32_t iSimd_4M_5M_6M   = setSimd3I(4000000, 5000000, 6000000); //!< Term for 32bit add/sub
const gmx_simd_int32_t iSimd_5M_7M_9M   = setSimd3I(5000000, 7000000, 9000000); //!< Sum for 32bit add/sub
const gmx_simd_int32_t iSimd_0xF0F0F0F0 = gmx_simd_set1_i(0xF0F0F0F0); //!< Integer bitpattern
const gmx_simd_int32_t iSimd_0xCCCCCCCC = gmx_simd_set1_i(0xCCCCCCCC); //!< Integer bitpattern
#endif

#ifdef GMX_SIMD_HAVE_REAL
/// Return stl vector from floating-point simd. The vector will have the same
/// length as the SIMD width.
::std::vector<real>
simd2Vector(const gmx_simd_real_t simd)
{
    real                mem[GMX_SIMD_REAL_WIDTH*2];
    real *              p = gmx_simd_align_r(mem);

    gmx_simd_store_r(p, simd);
    std::vector<real>   v(p, p+GMX_SIMD_REAL_WIDTH);

    return v;
}

/// Return floating-point simd value from stl vector. If the vector is longer
/// than SIMD width, only the first elements will be used. If it is shorter,
/// the contents will be repeated to fill the SIMD register.
gmx_simd_real_t
vector2Simd(const std::vector<real> &v)
{
    real                mem[GMX_SIMD_REAL_WIDTH*2];
    real *              p = gmx_simd_align_r(mem);

    for (int i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        p[i] = v[i % v.size()];  // repeat vector contents to fill simd width
    }
    return gmx_simd_load_r(p);
}

/// Set SIMD register contents from three real values. Our reason for using
/// three values is that 3 is not a factor in any known SIMD width, so this way
/// there will not be any simple repeated patterns e.g. between the low/high
/// 64/128/256 bits in the SIMD register, which could hide subtle bugs.
gmx_simd_real_t
setSimd3R(real r0, real r1, real r2)
{
    std::vector<real> v(3);
    v[0] = r0;
    v[1] = r1;
    v[2] = r2;
    return vector2Simd(v);
}

/// Simd datatypes are typically implemented as quite low-level
/// compiler-specific constructs, and at least for intel we have not
/// managed to apply the operator <<, but use a routine instead.
std::string
printSimdReal(const gmx_simd_real_t simd)
{
    std::stringstream   buf;
    std::vector<real>   v = simdTest::simd2Vector(simd);
    buf << "[ ";
    buf << std::setprecision(20);
    copy(v.begin(), v.end(), std::ostream_iterator<real>(buf, " "));
    buf << "]";
    return buf.str();
}
#endif


#ifdef GMX_SIMD_HAVE_REAL
/// Internal implementation - use the macro GMX_ASSERT_SIMD_REAL_EQ() instead.
/// This routine is designed according to the Google test specs, so the char
/// strings will describe the arguments to the macro, while the SIMD and
/// tolerance arguments are used to decide if the values are approximately equal.
testing::AssertionResult
compareSimdRealUlp(const char *  refExpr,    const char *  tstExpr,
                   const char *  ulpTolExpr, const char *  absTolExpr,
                   const gmx_simd_real_t ref, const gmx_simd_real_t tst, gmx_int64_t ulpTol, real absTol)
{
    std::vector<real>             vref = simdTest::simd2Vector(ref);
    std::vector<real>             vtst = simdTest::simd2Vector(tst);
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
        eq     = eq && ( fabs(vref[i]-vtst[i]) < absTol );
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
        eq         = eq && (ulpDiff[i] <= ulpTol);
    }

    if (eq == true)
    {
        return ::testing::AssertionSuccess();
    }
    else
    {
        return ::testing::AssertionFailure()
               << "Failing ulp comparison between " << refExpr << " and " << tstExpr << std::endl
               << "Requested ulp tolerance: " << ulpTolExpr << " (" << ulpTol << ")" << std::endl
               << "Requested abs tolerance: " << absTolExpr << " (" << absTol << ")" << std::endl
               << "Ref  values: " << printSimdReal(ref) << std::endl
               << "SIMD values: " << printSimdReal(tst) << std::endl
               << "Ulp diff.:   " << ::testing::PrintToString(ulpDiff) << std::endl;
    }
}

/// Internal implementation - use the macro GMX_ASSERT_SIMD_REAL_EQ() instead.
/// This routine expands the scalar reference to a SIMD register and performs
/// the comparison on all elements.
testing::AssertionResult
compareSimdRealUlp(const char * refExpr,    const char * tstExpr,
                   const char * ulpTolExpr, const char *absTolExpr,
                   real ref, const gmx_simd_real_t tst, gmx_int64_t ulpTol, real absTol)
{
    return compareSimdRealUlp(refExpr, tstExpr, ulpTolExpr, absTolExpr, gmx_simd_set1_r(ref), tst, ulpTol, absTol);
}

/// Internal implementation - use the macro GMX_ASSERT_SIMD_FUNC_NEAR instead.
/// The reference and SIMD functions are called for 10,000 points over the
/// specified range. They are considered equal if they are either within
/// the absolute tolerance (with same sign), or within the ULP tolerance.
::testing::AssertionResult
compareSimdMathFunction(const char * refFuncExpr, const char *simdFuncExpr,
                        const char * refRange,
                        const char * ulpTolExpr,  const char *  absTolExpr,
                        real refFunc(real x),     gmx_simd_real_t simdFunc(gmx_simd_real_t x),
                        std::pair<real, real> range,
                        gmx_int64_t ulpTol, real absTol)
{
    std::vector<real>            vsimd;
    const int                    npoints = 10000;
    real                         x, dx, refval;
    gmx_simd_real_t              simdx;
    gmx_int64_t                  ulpDiff, maxUlpDiff;
    real                         maxUlpDiffPos;
    real                         refValMaxUlpDiff, simdValMaxUlpDiff;
    bool                         eq, signOk;
    int                          i;
#    ifdef GMX_DOUBLE
    union {
        double r; gmx_int64_t i;
    } conv0, conv1;
#    else
    union {
        float  r; int i;
    } conv0, conv1;
#    endif

    maxUlpDiff = 0;
    dx         = (range.second-range.first)/npoints;

    for (x = range.first; x < range.second; x += dx)
    {
        simdx   = gmx_simd_set1_r(x);
        refval  = refFunc(x);
        vsimd   = simd2Vector(simdFunc(simdx));

        for (i = 0, eq = true, signOk = true; i < GMX_SIMD_REAL_WIDTH && eq == true; i++)
        {
            eq     = eq && ( fabs(refval-vsimd[i]) < absTol );
            signOk = signOk && ( refval*vsimd[i] >= 0 );
        }
        if (eq == true)
        {
            // Go to next point if everything within absolute tolerance
            continue;
        }
        else if (signOk == false)
        {
            return ::testing::AssertionFailure()
                   << "Failing math function comparison due to sign differences." << std::endl
                   << "Reference function: " << refFuncExpr << std::endl
                   << "Simd function:      " << simdFuncExpr << std::endl
                   << "Test range is " << refRange << " ( " << range.first << " , " << range.second << " ) " << std::endl
                   << "First sign difference for x=" << std::setprecision(20) << x << std::endl
                   << "Ref value:  " << std::setprecision(20) << refval << std::endl
                   << "SIMD value: " << std::setprecision(20) << vsimd[0] << std::endl;
        }

        for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
        {
            conv0.r = refval;
            conv1.r = vsimd[i];
            ulpDiff = llabs(conv0.i-conv1.i);
            if (ulpDiff > maxUlpDiff)
            {
                maxUlpDiff        = ulpDiff;
                maxUlpDiffPos     = x;
                refValMaxUlpDiff  = refval;
                simdValMaxUlpDiff = vsimd[0];
            }
        }
    }

    if (maxUlpDiff <= ulpTol)
    {
        return ::testing::AssertionSuccess();
    }
    else
    {
        return ::testing::AssertionFailure()
               << "Failing ulp comparison between " << refFuncExpr << " and " << simdFuncExpr << std::endl
               << "Requested ulp tolerance: " << ulpTolExpr << " (" << ulpTol << ")" << std::endl
               << "Requested abs tolerance: " << absTolExpr << " (" << absTol << ")" << std::endl
               << "Largest Ulp difference occurs for x=" << std::setprecision(20) << maxUlpDiffPos << std::endl
               << "Ref  value: " << std::setprecision(20) << refValMaxUlpDiff << std::endl
               << "SIMD value: " << std::setprecision(20) << simdValMaxUlpDiff << std::endl
               << "Ulp diff.:   " << std::setprecision(20) << maxUlpDiff << std::endl;
    }
}
#endif

#ifdef GMX_SIMD_HAVE_INT32
//! Convert integer SIMD to STL vector
std::vector<int>
simd2Vector(const gmx_simd_int32_t simd)
{
    int                 mem[GMX_SIMD_INT32_WIDTH*2];
    int *               p = gmx_simd_align_i(mem);

    gmx_simd_store_i(p, simd);
    std::vector<int>    v(p, p+GMX_SIMD_INT32_WIDTH);

    return v;
}
    
/// Return floating-point simd value from stl vector. If the vector is longer
/// than SIMD width, only the first elements will be used. If it is shorter,
/// the contents will be repeated to fill the SIMD register.
gmx_simd_int32_t
vector2Simd(const std::vector<int> &v)
{
    int                 mem[GMX_SIMD_INT32_WIDTH*2];
    int *               p = gmx_simd_align_i(mem);

    for (int i = 0; i < GMX_SIMD_INT32_WIDTH; i++)
    {
        p[i] = v[i % v.size()];  // repeat vector contents to fill simd width
    }
    return gmx_simd_load_i(p);
}

/// Set SIMD register contents from three int values. Our reason for using
/// three values is that 3 is not a factor in any known SIMD width, so this way
/// there will not be any simple repeated patterns e.g. between the low/high
/// 64/128/256 bits in the SIMD register, which could hide subtle bugs.
gmx_simd_int32_t
setSimd3I(int i0, int i1, int i2)
{
    std::vector<int> v(3);
    v[0] = i0;
    v[1] = i1;
    v[2] = i2;
    return vector2Simd(v);
}

//! Returns a string with a printout of SIMD contents
std::string
printSimdInt32(const gmx_simd_int32_t simd)
{
    std::stringstream   buf;
    std::vector<int>    v = simdTest::simd2Vector(simd);
    buf << "[ ";
    copy(v.begin(), v.end(), std::ostream_iterator<int>(buf, " "));
    buf << "]";
    return buf.str();
}

/// Internal routine. Use the macro GMX_ASSERT_SIMD_INT_EQ instead.
/// This function checks if each element in the SIMD test value is identical
/// to the corresponding element in ref.
testing::AssertionResult
compareSimdInt32(const char *  refExpr,      const char *  tstExpr,
                 const gmx_simd_int32_t ref, const gmx_simd_int32_t tst)
{
    std::vector<int>    vref = simdTest::simd2Vector(ref);
    std::vector<int>    vtst = simdTest::simd2Vector(tst);
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

/// Internal routine. Use the macro GMX_ASSERT_SIMD_INT_EQ instead.
/// This function expands the scalar reference into a SIMD integer and checks
/// that each element of test is equal to this value.
testing::AssertionResult
compareSimdInt32(const char * refExpr,    const char * tstExpr,
                 int ref, const gmx_simd_int32_t tst)
{
    return compareSimdInt32(refExpr, tstExpr, gmx_simd_set1_i(ref), tst);
}
#endif

} // namespace simdTest
