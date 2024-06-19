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

#include "gromacs/simd/simd_math.h"

#include "config.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>

#include <algorithm>
#include <iomanip>
#include <limits>
#include <map>
#include <ostream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/tests/data.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

#include "simd.h"

#if GMX_SIMD

namespace gmx
{
namespace test
{

/*! \cond internal */
/*! \addtogroup module_simd */
/*! \{ */

#    if GMX_SIMD_HAVE_REAL

class SimdMathTest : public SimdTest
{
public:
    /*! \brief Type for half-open intervals specifying test ranges */
    typedef std::pair<real, real> Range;

    /*! \brief Control what is considered matching values
     *
     * Normal simply means that we request the values to be equal
     * to within the specified tolerance.
     * However, there are also two more cases that are special:
     *
     * - Even if we only care about normal (i.e., not denormal) values, some math
     *   libraries might clamp the value to zero, which means our SIMD output
     *   might not match their values. By using MatchRule::Dtz, we will consider
     *   all values both from the reference and test functions that are within the
     *   requested ulp tolerance of a denormal number to be equivalent to 0.0.
     * - For some older architectures without fused multiply-add units (e.g. x86 SSE2),
     *   we might end up clamping the results to zero just before reaching
     *   denormal output, since the intermediate results e.g. in polynomial
     *   approximations can be smaller than the final one. We often simply don't
     *   care about those values, and then one can use
     *   MatchRule::ReferenceOrZero to allow the test value to either match
     *   the reference or be zero.
     */
    enum class MatchRule
    {
        Normal, //!< Match function values
        Dtz, //!< Match function values after setting denormals to zero both in test and reference
        ReferenceOrZero, //!< Test values can either match reference or be zero
    };

    const std::map<MatchRule, std::string> matchRuleNames_ = {
        { MatchRule::Normal, "Test should match reference." },
        { MatchRule::Dtz, "Test should match reference, with denormals treated as 0.0." },
        { MatchRule::ReferenceOrZero, "Test should match reference or 0.0." }
    };

    /*! \brief Settings used for simd math function comparisons */
    struct CompareSettings
    {
        Range        range;     //!< Range over which to test function
        std::int64_t ulpTol;    //!< Ulp tolerance
        real         absTol;    //!< Absolute tolerance
        MatchRule    matchRule; //!< Decide what we consider a match
    };

    ::testing::AssertionResult compareSimdMathFunction(const char*            refFuncExpr,
                                                       const char*            simdFuncExpr,
                                                       const char*            compareSettingsExpr,
                                                       real                   refFunc(real x),
                                                       SimdReal gmx_simdcall  simdFunc(SimdReal x),
                                                       const CompareSettings& compareSettings);

    /*! \brief Generate test point vector
     *
     *  \param range  The test interval, half open. Upper limit is not included.
     *                Pass by value, since we need to modify in method anyway.
     *  \param points Number of points to generate. This might be increased
     *                slightly to account both for extra special values like 0.0
     *                and the SIMD width.
     *
     * This routine generates a vector with test points separated by constant
     * multiplicative factors, based on the range and number of points in the
     * class. If the range includes both negative and positive values, points
     * will be generated separately for the negative/positive intervals down
     * to the smallest real number that can be represented, and we also include
     * 0.0 explicitly.
     *
     * This is highly useful for large test ranges. For example, with a linear
     * 1000-point division of the range (1,1e10) the first three values to test
     * would be 1, 10000000.999, and 20000000.998, etc. For large values we would
     * commonly hit the point where adding the small delta has no effect due to
     * limited numerical precision.
     * When we instead use this routine, the values will be 1, 1.0239, 1.0471, etc.
     * This will spread the entropy over all bits in the IEEE754 representation,
     * and be a much better test of all potential input values.
     *
     *  \note We do not use the static variable s_nPoints in the parent class
     *        to avoid altering any value the user has set on the command line; since
     *        it's a static member, changing it would have permanent effect.
     */
    static std::vector<real> generateTestPoints(Range range, std::size_t points);

    /*! \brief Test routine for the test point vector generation
     */
    static void generateTestPointsTest();
};

/*! \brief Test approximate equality of SIMD vs reference version of a function.
 *
 * This macro takes vanilla C and SIMD flavors of a function and tests it with
 * the number of points, range, and tolerances specified by the test fixture class.
 *
 * The third option controls the range, tolerances, and match settings.
 */
#        define GMX_EXPECT_SIMD_FUNC_NEAR(refFunc, tstFunc, compareSettings) \
            EXPECT_PRED_FORMAT3(compareSimdMathFunction, refFunc, tstFunc, compareSettings)

std::vector<real> SimdMathTest::generateTestPoints(Range inputRange, std::size_t inputPoints)
{

    std::vector<real> testPoints;
    testPoints.reserve(inputPoints);

    GMX_RELEASE_ASSERT(inputRange.first < inputRange.second,
                       "The start of the interval must come before the end");

    std::vector<Range> testRanges;

    if (inputRange.first < 0 && inputRange.second > 0)
    {
        testRanges.emplace_back(Range({ inputRange.first, -std::numeric_limits<real>::min() }));
        testRanges.emplace_back(Range({ 0.0, inputRange.second }));
    }
    else
    {
        if (inputRange.second == 0)
        {
            inputRange.second = -std::numeric_limits<real>::min();
            inputRange.first  = std::min(inputRange.first, inputRange.second);
        }
        testRanges.push_back(inputRange);
    }

    for (Range& range : testRanges)
    {
        std::size_t points = inputPoints / testRanges.size();

        // The value 0 is special, and can only occur at the start of
        // the interval after the corrections outside this loop.
        // Add it explicitly, and adjust the interval to continue
        // at the first valid non-zero positive number.
        if (range.first == 0)
        {
            testPoints.push_back(0.0);
            range.first = std::numeric_limits<real>::min();
            points--; // Used one point
        }

        union
        {
            real                                                                               r;
            std::conditional<sizeof(real) == sizeof(double), std::int64_t, std::int32_t>::type i;
        } low, high, x;

        low.r  = range.first;
        high.r = range.second;

        // IEEE754 floating-point numbers have the cool property that for any range of
        // constant sign, for all non-zero numbers a constant (i.e., linear) difference
        // in the bitwise representation corresponds to a constant multiplicative factor.
        //
        // Divide the ulp difference evenly
        std::int64_t ulpDiff = high.i - low.i;
        // dividend and divisor must both be signed types
        std::int64_t ulpDelta    = ulpDiff / static_cast<std::int64_t>(points);
        std::int64_t minUlpDelta = (ulpDiff > 0) ? 1 : -1;

        if (ulpDelta == 0)
        {
            // Very short interval or very many points caused round-to-zero.
            // Select the smallest possible change, which is one ulp (with correct sign)
            ulpDelta = minUlpDelta;
            points   = std::abs(ulpDiff);
        }

        x.r = low.r;
        // Use an index-based loop to avoid floating-point comparisons with
        // values that might have overflowed. Save one point for the very last
        // bitwise value that is part of the interval
        for (std::size_t i = 0; i < points - 1; i++)
        {
            testPoints.push_back(x.r);
            x.i += ulpDelta;
        }

        // Make sure we test the very last point that is inside the interval
        x.r = high.r;
        x.i -= minUlpDelta;
        testPoints.push_back(x.r);
    }
    return testPoints;
}

/*! \brief Implementation routine to compare SIMD vs reference functions.
 *
 * \param refFuncExpr         Description of reference function expression
 * \param simdFuncExpr        Description of SIMD function expression
 * \param compareSettingsExpr Description of compareSettings
 * \param refFunc             Reference math function pointer
 * \param simdFunc            SIMD math function pointer
 * \param compareSettings     Structure with the range, tolerances, and
 *                            matching rules to use for the comparison.
 *
 * \note You should not never call this function directly,  but use the
 *       macro GMX_EXPECT_SIMD_FUNC_NEAR(refFunc,tstFunc,matchRule) instead.
 */
::testing::AssertionResult SimdMathTest::compareSimdMathFunction(const char* refFuncExpr,
                                                                 const char* simdFuncExpr,
                                                                 const char gmx_unused* compareSettingsExpr,
                                                                 real refFunc(real x),
                                                                 SimdReal gmx_simdcall simdFunc(SimdReal x),
                                                                 const CompareSettings& compareSettings)
{
    std::vector<real> vx(GMX_SIMD_REAL_WIDTH);
    std::vector<real> vref(GMX_SIMD_REAL_WIDTH);
    std::vector<real> vtst(GMX_SIMD_REAL_WIDTH);
    real              absDiff;
    std::int64_t      ulpDiff;
    std::int64_t      maxUlpDiff = 0;
    real              maxUlpDiffPos;
    real              refValMaxUlpDiff, simdValMaxUlpDiff;
    const int         niter = s_nPoints / GMX_SIMD_REAL_WIDTH;

    union
    {
        real                                                                               r;
        std::conditional<sizeof(real) == sizeof(double), std::int64_t, std::int32_t>::type i;
    } conv0, conv1;

    // Allow zero-size intervals - nothing to test means we succeeded at it
    if (compareSettings.range.first == compareSettings.range.second)
    {
        ::testing::AssertionSuccess();
    }

    // Calculate the tolerance limit to use for denormals - we want
    // values that are within the ulp tolerance of denormals to be considered matching
    conv0.r = std::numeric_limits<real>::min();
    conv0.i += compareSettings.ulpTol - 1; // min() itself is not denormal, but one ulp larger
    const real denormalLimit = conv0.r;

    // We want to test as many diverse bit combinations as possible over the range requested,
    // and in particular do it evenly spaced in bit-space.
    // Due to the way IEEE754 floating-point is represented, that means we should have a
    // constant multiplicative factor between adjacent values. This gets a bit complicated
    // when we have both positive and negative values, so we offload the generation of the
    // specific testing values to a separate routine
    std::vector<real> testPoints = generateTestPoints(compareSettings.range, s_nPoints);

    size_t pointIndex = 0;

    for (int iter = 0; iter < niter; iter++)
    {
        for (int i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
        {
            vx[i]   = testPoints[pointIndex];
            vref[i] = refFunc(vx[i]);
            // If we reach the end of the points, stop increasing index so we pad with
            // extra copies of the last element up to the SIMD width
            if (pointIndex + 1 < testPoints.size())
            {
                pointIndex++;
            }
        }
        vtst = simdReal2Vector(simdFunc(vector2SimdReal(vx)));

        bool absOk = true, signOk = true;
        for (int i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
        {
            if (compareSettings.matchRule == MatchRule::Dtz && std::abs(vref[i]) <= denormalLimit
                && std::abs(vtst[i]) <= denormalLimit)
            {
                continue;
            }

            if (compareSettings.matchRule == MatchRule::ReferenceOrZero && vtst[i] == 0.0)
            {
                // If we accept 0.0 for the test function, we can continue to the next loop iteration.
                continue;
            }

            absDiff = std::abs(vref[i] - vtst[i]);
            absOk   = absOk && (absDiff < compareSettings.absTol);
            signOk  = signOk && ((vref[i] >= 0 && vtst[i] >= 0) || (vref[i] <= 0 && vtst[i] <= 0));

            if (absDiff >= compareSettings.absTol)
            {
                /* We replicate the trivial ulp differences comparison here rather than
                 * calling the lower-level routine for comparing them, since this enables
                 * us to run through the entire test range and report the largest deviation
                 * without lots of extra glue routines.
                 */
                conv0.r = vref[i];
                conv1.r = vtst[i];
                ulpDiff = std::llabs(conv0.i - conv1.i);
                if (ulpDiff > maxUlpDiff)
                {
                    maxUlpDiff        = ulpDiff;
                    maxUlpDiffPos     = vx[i];
                    refValMaxUlpDiff  = vref[i];
                    simdValMaxUlpDiff = vtst[i];
                }
            }
        }
        if ((!absOk) && (!signOk))
        {
            return ::testing::AssertionFailure()
                   << "Failing SIMD math function comparison due to sign differences." << std::endl
                   << "Reference function: " << refFuncExpr << std::endl
                   << "Simd function:      " << simdFuncExpr << std::endl
                   << "Test range is ( " << compareSettings.range.first << " , "
                   << compareSettings.range.second << " ) " << std::endl
                   << "Match rule: " << matchRuleNames_.at(compareSettings.matchRule) << std::endl
                   << "First sign difference around x=" << std::setprecision(20)
                   << ::testing::PrintToString(vx) << std::endl
                   << "Ref values:  " << std::setprecision(20) << ::testing::PrintToString(vref)
                   << std::endl
                   << "SIMD values: " << std::setprecision(20) << ::testing::PrintToString(vtst)
                   << std::endl;
        }
    }

    GMX_RELEASE_ASSERT(compareSettings.ulpTol >= 0, "Invalid ulp value.");
    if (maxUlpDiff <= compareSettings.ulpTol)
    {
        return ::testing::AssertionSuccess();
    }
    else
    {
        return ::testing::AssertionFailure()
               << "Failing SIMD math function ulp comparison between " << refFuncExpr << " and "
               << simdFuncExpr << std::endl
               << "Requested ulp tolerance: " << compareSettings.ulpTol << std::endl
               << "Requested abs tolerance: " << compareSettings.absTol << std::endl
               << "Match rule: " << matchRuleNames_.at(compareSettings.matchRule) << std::endl
               << "Largest Ulp difference occurs for x=" << std::setprecision(20) << maxUlpDiffPos
               << std::endl
               << "Ref  values: " << std::setprecision(20) << refValMaxUlpDiff << std::endl
               << "SIMD values: " << std::setprecision(20) << simdValMaxUlpDiff << std::endl
               << "Ulp diff.:   " << std::setprecision(20) << maxUlpDiff << std::endl;
    }
}

// Actual routine to generate a small set of test points in current precision. This will
// be called by either the double or single precision test fixture, since we need different
// test names to compare to the right reference data.
void SimdMathTest::generateTestPointsTest()
{
    int                             points(10);
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    std::vector<real> result;

    result = generateTestPoints(Range(-1e10, -1), points);
    checker.checkSequence(result.begin(), result.end(), "Test points for interval [-1e10,-1[");

    result = generateTestPoints(Range(-1e10, -1e-10), points);
    checker.checkSequence(result.begin(), result.end(), "Test points for interval [-1e10, -1e-10[");

    result = generateTestPoints(Range(1, 1e10), points);
    checker.checkSequence(result.begin(), result.end(), "Test points for interval [1, 1e10[");

    result = generateTestPoints(Range(1e-10, 1e10), points);
    checker.checkSequence(result.begin(), result.end(), "Test points for interval [1e-10, 1e10[");

    result = generateTestPoints(Range(-1e10, 1e-10), points);
    checker.checkSequence(result.begin(), result.end(), "Test points for interval [-1e10, 1e-10[");

    result = generateTestPoints(Range(-1e-10, 1e-10), points);
    checker.checkSequence(result.begin(), result.end(), "Test points for interval [-1e-10, 1e-10[");

    result = generateTestPoints(Range(-1e-10, 1e10), points);
    checker.checkSequence(result.begin(), result.end(), "Test points for interval [-1e-10, 1e10[");

    result = generateTestPoints(Range(-1e10, 1e10), points);
    checker.checkSequence(result.begin(), result.end(), "Test points for interval [-1e10, 1e10[");

    result = generateTestPoints(Range(-1000, 0), points);
    checker.checkSequence(result.begin(), result.end(), "Test points for interval [-1000, 0[");

    result = generateTestPoints(Range(0, 1000), points);
    checker.checkSequence(result.begin(), result.end(), "Test points for interval [0, 1000[");
}

/*! \} */
/*! \endcond */


// Actual math function tests below

/*! \cond internal */
/*! \addtogroup module_simd */
/*! \{ */

namespace
{

// Reference data is selected based on test name, so make the test name precision-dependent
#        if GMX_DOUBLE
TEST_F(SimdMathTest, generateTestPointsDouble)
{
    generateTestPointsTest();
}
#        else
TEST_F(SimdMathTest, generateTestPointsFloat)
{
    generateTestPointsTest();
}
#        endif

TEST_F(SimdMathTest, copysign)
{
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(-c0, c1, c2),
                            copysign(setSimdRealFrom3R(c0, c1, c2), setSimdRealFrom3R(-c3, c4, 0)));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(-c0, c1, c2),
                            copysign(setSimdRealFrom3R(-c0, -c1, -c2), setSimdRealFrom3R(-c3, c4, 0)));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(c0, -c1, c2),
                            copysign(setSimdRealFrom3R(c0, c1, c2), setSimdRealFrom3R(c3, -c4, 0)));
    GMX_EXPECT_SIMD_REAL_EQ(setSimdRealFrom3R(c0, -c1, c2),
                            copysign(setSimdRealFrom3R(-c0, -c1, -c2), setSimdRealFrom3R(c3, -c4, 0)));
}

/*! \brief Function wrapper to evaluate reference 1/sqrt(x) */
real refInvsqrt(real x)
{
    return 1.0 / std::sqrt(x);
}

TEST_F(SimdMathTest, invsqrt)
{
    const real      low  = std::numeric_limits<float>::min();
    const real      high = std::numeric_limits<float>::max();
    CompareSettings settings{ Range(low, high), ulpTol_, absTol_, MatchRule::Normal };

    GMX_EXPECT_SIMD_FUNC_NEAR(refInvsqrt, invsqrt, settings);
}

TEST_F(SimdMathTest, maskzInvsqrt)
{
    SimdReal x   = setSimdRealFrom3R(c1, 0.0, c2);
    SimdBool m   = (setZero() < x);
    SimdReal ref = setSimdRealFrom3R(1.0 / std::sqrt(c1), 0.0, 1.0 / std::sqrt(c2));
    GMX_EXPECT_SIMD_REAL_NEAR(ref, maskzInvsqrt(x, m));
}

/*! \brief Function wrapper to return first result when testing \ref invsqrtPair */
SimdReal gmx_simdcall tstInvsqrtPair0(SimdReal x)
{
    SimdReal r0, r1;
    invsqrtPair(x, x, &r0, &r1);
    return r0;
}

/*! \brief Function wrapper to return second result when testing \ref invsqrtPair */
SimdReal gmx_simdcall tstInvsqrtPair1(SimdReal x)
{
    SimdReal r0, r1;
    invsqrtPair(x, x, &r0, &r1);
    return r1;
}

TEST_F(SimdMathTest, invsqrtPair)
{
    const real low  = std::numeric_limits<float>::min();
    const real high = std::numeric_limits<float>::max();

    // Accuracy conversions lose a bit of accuracy compared to all-double,
    // so increase the tolerance to 4*ulpTol_
    CompareSettings settings{ Range(low, high), 4 * ulpTol_, absTol_, MatchRule::Normal };

    GMX_EXPECT_SIMD_FUNC_NEAR(refInvsqrt, tstInvsqrtPair0, settings);
    GMX_EXPECT_SIMD_FUNC_NEAR(refInvsqrt, tstInvsqrtPair1, settings);
}

/*! \brief Function wrapper to evaluate reference sqrt(x) */
real refSqrt(real x)
{
    return std::sqrt(x);
}

TEST_F(SimdMathTest, sqrt)
{
    // Since the first lookup step is sometimes performed in single precision,
    // our SIMD sqrt can only handle single-precision input values, even when
    // compiled in double precision.

    const real      minFloat     = std::numeric_limits<float>::min();
    const real      minSafeFloat = minFloat * 10;
    const real      maxSafeFloat = std::numeric_limits<float>::max() * 0.1;
    CompareSettings settings;
    // The accuracy conversions lose a bit of extra accuracy compared to
    // doing the iterations in all-double.
    setUlpTol(4 * ulpTol_);

    // First test that 0.0 and a few other values work
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(0, std::sqrt(c1), std::sqrt(c2)),
                              sqrt(setSimdRealFrom3R(0, c1, c2)));

#        if GMX_DOUBLE
    // As mentioned above, we cannot guarantee that very small double precision
    // input values (below std::numeric_limits<float>::min()) are handled correctly,
    // so our implementation will clamp it to zero. In this range we allow either
    // the correct value or zero, but it's important that it does not result in NaN or Inf values.
    //
    // This test range must not be called for single precision, since if we try to divide
    // the interval (0.0, low( in npoints we will try to multiply by factors so small that
    // they end up being flushed to zero, and the loop would never end.
    settings = { Range(0.0, minFloat), ulpTol_, absTol_, MatchRule::ReferenceOrZero };
    GMX_EXPECT_SIMD_FUNC_NEAR(refSqrt, sqrt, settings);
#        endif

    // Next range: Just about minFloat the lookup should always work, but the results
    // might be a bit fragile due to issues with the N-R iterations being flushed to zero
    // for denormals. We can probably relax the latter in double precision, but since we
    // anyway cannot handle numbers that cannot be represented in single it's not worth
    // worrying too much about whether we have zero or an exact values around 10^-38....
    settings = { Range(minFloat, minSafeFloat), ulpTol_, absTol_, MatchRule::ReferenceOrZero };
    GMX_EXPECT_SIMD_FUNC_NEAR(refSqrt, sqrt, settings);

    settings = { Range(minSafeFloat, maxSafeFloat), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(refSqrt, sqrt, settings);
}

TEST_F(SimdMathTest, sqrtUnsafe)
{
    const real minSafeFloat = std::numeric_limits<float>::min() * 10;
    const real maxSafeFloat = std::numeric_limits<float>::max() * 0.1;

    // The accuracy conversions lose a bit of extra accuracy compared to
    // doing the iterations in all-double, so we use 4*ulpTol_
    setUlpTol(4 * ulpTol_);

    CompareSettings settings{ Range(minSafeFloat, maxSafeFloat), 4 * ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(refSqrt, sqrt<MathOptimization::Unsafe>, settings);
}

/*! \brief Function wrapper to evaluate reference 1/x */
real refInv(real x)
{
    return 1.0 / x;
}

TEST_F(SimdMathTest, inv)
{
    // Since the first lookup step is sometimes performed in single precision,
    // our SIMD 1/x can only handle single-precision input values, even when
    // compiled in double precision.

    // Relevant threshold points
    const real minSafeFloat = std::numeric_limits<float>::min()
                              * 10; // X value guaranteed not to result in Inf intermediates for 1/x calc.
    const real maxSafeFloat = std::numeric_limits<float>::max()
                              * 0.1; // X value guaranteed not to result in DTZ intermediates for 1/x calc.
    // Scale highest value by 1-eps, since we will do some arithmetics on this value
    const real maxFloat =
            std::numeric_limits<float>::max() * (1.0 - std::numeric_limits<float>::epsilon());
    CompareSettings settings;

    // Danger zone where intermediates might be flushed to zero and produce 1/x==0.0
    settings = { Range(-maxFloat, -maxSafeFloat), ulpTol_, absTol_, MatchRule::ReferenceOrZero };
    GMX_EXPECT_SIMD_FUNC_NEAR(refInv, inv, settings);

    // Normal checks for x < 0
    settings = { Range(-maxSafeFloat, -minSafeFloat), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(refInv, inv, settings);

    // We do not care about the small range -minSafeFloat < x < +minSafeFloat where the result can be +/- Inf, since we don't require strict IEEE754.

    // Normal checks for x > 0
    settings = { Range(minSafeFloat, maxSafeFloat), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(refInv, inv, settings);

    // Danger zone where intermediates might be flushed to zero and produce 1/x==0.0
    settings = { Range(maxSafeFloat, maxFloat), ulpTol_, absTol_, MatchRule::ReferenceOrZero };
    GMX_EXPECT_SIMD_FUNC_NEAR(refInv, inv, settings);
}

TEST_F(SimdMathTest, maskzInv)
{
    SimdReal x   = setSimdRealFrom3R(c1, 0.0, c2);
    SimdBool m   = (setZero() < x);
    SimdReal ref = setSimdRealFrom3R(1.0 / c1, 0.0, 1.0 / c2);
    GMX_EXPECT_SIMD_REAL_NEAR(ref, maskzInv(x, m));
}

TEST_F(SimdMathTest, cbrt)
{
    const real low  = -std::numeric_limits<real>::max();
    const real high = std::numeric_limits<real>::max();

    CompareSettings settings{ Range(low, high), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::cbrt, cbrt, settings);
}

/*! \brief Function wrapper to evaluate reference 1/cbrt(x) */
real refInvCbrt(real x)
{
    return 1.0 / std::cbrt(x);
}

TEST_F(SimdMathTest, invcbrt)
{
    // Negative values first
    real low  = -std::numeric_limits<real>::max();
    real high = -std::numeric_limits<real>::min();

    CompareSettings settings{ Range(low, high), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(refInvCbrt, invcbrt, settings);

    // Positive values
    low      = std::numeric_limits<real>::min();
    high     = std::numeric_limits<real>::max();
    settings = { Range(low, high), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(refInvCbrt, invcbrt, settings);
}

TEST_F(SimdMathTest, log2)
{
    const real low  = std::numeric_limits<real>::min();
    const real high = std::numeric_limits<real>::max();

    CompareSettings settings{ Range(low, high), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::log2, log2, settings);
}

TEST_F(SimdMathTest, log)
{
    const real low  = std::numeric_limits<real>::min();
    const real high = std::numeric_limits<real>::max();

    CompareSettings settings{ Range(low, high), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::log, log, settings);
}

TEST_F(SimdMathTest, exp2)
{
    // Relevant threshold points
    constexpr real lowestReal = -std::numeric_limits<real>::max();
    constexpr real lowestRealThatProducesNormal =
            std::numeric_limits<real>::min_exponent
            - 1; // adding the significant corresponds to one more unit in exponent
    constexpr real lowestRealThatProducesDenormal =
            lowestRealThatProducesNormal
            - std::numeric_limits<real>::digits; // digits refer to bits in significand, so 24/53 for float/double
    constexpr real highestRealThatProducesNormal =
            std::numeric_limits<real>::max_exponent
            - 1; // adding the significant corresponds to one more unit in exponent
    CompareSettings settings;

    // Below subnormal range all results should be zero (so, match the reference)
    settings = { Range(lowestReal, lowestRealThatProducesDenormal), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp2, exp2, settings);

    // Subnormal range, require matching, but DTZ is fine
    settings = { Range(lowestRealThatProducesDenormal, lowestRealThatProducesNormal),
                 ulpTol_,
                 absTol_,
                 MatchRule::Dtz };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp2, exp2, settings);

    // Normal range, standard result expected
    settings = { Range(lowestRealThatProducesNormal, highestRealThatProducesNormal),
                 ulpTol_,
                 absTol_,
                 MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp2, exp2, settings);
}

TEST_F(SimdMathTest, exp2Unsafe)
{
    // The unsafe version is only defined in the normal range
    constexpr real lowestRealThatProducesNormal =
            std::numeric_limits<real>::min_exponent
            - 1; // adding the significant corresponds to one more unit in exponent
    constexpr real highestRealThatProducesNormal =
            std::numeric_limits<real>::max_exponent
            - 1; // adding the significant corresponds to one more unit in exponent

    CompareSettings settings{ Range(lowestRealThatProducesNormal, highestRealThatProducesNormal),
                              ulpTol_,
                              absTol_,
                              MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp2, exp2<MathOptimization::Unsafe>, settings);
}

TEST_F(SimdMathTest, exp)
{
    // Relevant threshold points. See the exp2 test for more details about the values; these are
    // simply scaled by log(2) due to the difference between exp2 and exp.
    const real lowestReal = -std::numeric_limits<real>::max();
    // In theory the smallest value should be (min_exponent-1)*log(2), but rounding after the multiplication will cause this
    // value to be a single ulp too low. This might cause failed tests on CPUs that use different DTZ modes for SIMD vs.
    // non-SIMD arithmetics (e.g. ARM v7), so multiply by (1.0-eps) to increase it by a single ulp.
    const real lowestRealThatProducesNormal = (std::numeric_limits<real>::min_exponent - 1)
                                              * std::log(2.0)
                                              * (1 - std::numeric_limits<real>::epsilon());
    const real lowestRealThatProducesDenormal =
            lowestRealThatProducesNormal - std::numeric_limits<real>::digits * std::log(2.0);
    const real highestRealThatProducesNormal =
            (std::numeric_limits<real>::max_exponent - 1) * std::log(2.0);
    CompareSettings settings;

    // Below subnormal range all results should be zero (so, match the reference)
    settings = { Range(lowestReal, lowestRealThatProducesDenormal), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp, exp, settings);

    // Subnormal range, require matching, but DTZ is fine
    settings = { Range(lowestRealThatProducesDenormal, lowestRealThatProducesNormal),
                 ulpTol_,
                 absTol_,
                 MatchRule::Dtz };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp, exp, settings);

    // Normal range, standard result expected
    settings = { Range(lowestRealThatProducesNormal, highestRealThatProducesNormal),
                 ulpTol_,
                 absTol_,
                 MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp, exp, settings);
}

TEST_F(SimdMathTest, expUnsafe)
{
    // See test of exp() for comments about test ranges
    const real lowestRealThatProducesNormal = (std::numeric_limits<real>::min_exponent - 1)
                                              * std::log(2.0)
                                              * (1 - std::numeric_limits<real>::epsilon());
    const real highestRealThatProducesNormal =
            (std::numeric_limits<real>::max_exponent - 1) * std::log(2.0);

    CompareSettings settings{ Range(lowestRealThatProducesNormal, highestRealThatProducesNormal),
                              ulpTol_,
                              absTol_,
                              MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp, exp<MathOptimization::Unsafe>, settings);
}

TEST_F(SimdMathTest, pow)
{
    // We already test the log2/exp2 components of pow() extensively above, and it's a very
    // simple single-line function, so here we just test a handful of values to catch typos
    // and then some special values.

    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::pow(c0, c3), std::pow(c1, c4), std::pow(c2, c5)),
                              pow(rSimd_c0c1c2, rSimd_c3c4c5));

    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::pow(c0, -c3), std::pow(c1, -c0), std::pow(c2, -c4)),
                              pow(rSimd_c0c1c2, rSimd_m3m0m4));

    // 0^0 = 1 , 0^c1=0, -c1^0=1
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(1.0, 0.0, 1.0),
                              pow(setSimdRealFrom3R(0, 0.0, -c1), setSimdRealFrom3R(0.0, c1, 0.0)));
}

TEST_F(SimdMathTest, powUnsafe)
{
    // We already test the log2/exp2 components of pow() extensively above, and it's a very
    // simple single-line function, so here we just test a handful of values to catch typos
    // and then some special values.

    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::pow(c0, c3), std::pow(c1, c4), std::pow(c2, c5)),
                              pow<MathOptimization::Unsafe>(rSimd_c0c1c2, rSimd_c3c4c5));

    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::pow(c0, -c3), std::pow(c1, -c0), std::pow(c2, -c4)),
                              pow<MathOptimization::Unsafe>(rSimd_c0c1c2, rSimd_m3m0m4));
}

/*! \brief Function wrapper for erf(x), with argument/return in default Gromacs precision.
 *
 * \note Single-precision erf() in some libraries can be slightly lower precision
 * than the SIMD flavor, so we use a cast to force double precision for reference.
 */
real refErf(real x)
{
    return std::erf(static_cast<double>(x));
}

TEST_F(SimdMathTest, erf)
{
    CompareSettings settings{ Range(-9, 9), ulpTol_, std::numeric_limits<real>::min(), MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(refErf, erf, settings);
}

/*! \brief Function wrapper for erfc(x), with argument/return in default Gromacs precision.
 *
 * \note Single-precision erfc() in some libraries can be slightly lower precision
 * than the SIMD flavor, so we use a cast to force double precision for reference.
 */
real refErfc(real x)
{
    return std::erfc(static_cast<double>(x));
}

TEST_F(SimdMathTest, erfc)
{
    // Our erfc algorithm has 4 ulp accuracy, so relax tolerance a bit to 4*ulpTol
    CompareSettings settings{ Range(-9, 9), 4 * ulpTol_, std::numeric_limits<real>::min(), MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(refErfc, erfc, settings);
}

TEST_F(SimdMathTest, sin)
{
    CompareSettings settings{ Range(-8 * M_PI, 8 * M_PI), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::sin, sin, settings);

    // Range reduction leads to accuracy loss, so we might want higher tolerance here
    settings = { Range(-10000, 10000), 2 * ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::sin, sin, settings);
}

TEST_F(SimdMathTest, cos)
{
    CompareSettings settings{ Range(-8 * M_PI, 8 * M_PI), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::cos, cos, settings);

    // Range reduction leads to accuracy loss, so we might want higher tolerance here
    settings = { Range(-10000, 10000), 2 * ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::cos, cos, settings);
}

TEST_F(SimdMathTest, tan)
{
    // Tan(x) is a little sensitive due to the division in the algorithm.
    // Rather than using lots of extra FP operations, we accept the algorithm
    // presently only achieves a ~3 ulp error and use the medium tolerance.
    CompareSettings settings{ Range(-8 * M_PI, 8 * M_PI), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::tan, tan, settings);

    // Range reduction leads to accuracy loss, so we might want higher tolerance here
    settings = { Range(-10000, 10000), 2 * ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::tan, tan, settings);
}

TEST_F(SimdMathTest, asin)
{
    // Our present asin(x) algorithm achieves 2-3 ulp accuracy
    CompareSettings settings{ Range(-1, 1), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::asin, asin, settings);
}

TEST_F(SimdMathTest, acos)
{
    // Our present acos(x) algorithm achieves 2-3 ulp accuracy
    CompareSettings settings{ Range(-1, 1), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::acos, acos, settings);
}

TEST_F(SimdMathTest, atan)
{
    // Our present atan(x) algorithm achieves 1 ulp accuracy
    CompareSettings settings{ Range(-10000, 10000), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::atan, atan, settings);
}

TEST_F(SimdMathTest, atan2)
{
    // test each quadrant
    GMX_EXPECT_SIMD_REAL_NEAR(
            setSimdRealFrom3R(std::atan2(c0, c3), std::atan2(c1, c4), std::atan2(c2, c5)),
            atan2(rSimd_c0c1c2, rSimd_c3c4c5));
    GMX_EXPECT_SIMD_REAL_NEAR(
            setSimdRealFrom3R(std::atan2(-c0, c3), std::atan2(-c1, c4), std::atan2(-c2, c5)),
            atan2(rSimd_m0m1m2, rSimd_c3c4c5));
    GMX_EXPECT_SIMD_REAL_NEAR(
            setSimdRealFrom3R(std::atan2(-c0, -c3), std::atan2(-c1, -c0), std::atan2(-c2, -c4)),
            atan2(rSimd_m0m1m2, rSimd_m3m0m4));
    GMX_EXPECT_SIMD_REAL_NEAR(
            setSimdRealFrom3R(std::atan2(c0, -c3), std::atan2(c1, -c0), std::atan2(c2, -c4)),
            atan2(rSimd_c0c1c2, rSimd_m3m0m4));

    // cases important for calculating angles
    // values on coordinate axes
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::atan2(0, c0), std::atan2(0, c1), std::atan2(0, c2)),
                              atan2(setZero(), rSimd_c0c1c2));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::atan2(c0, 0), std::atan2(c1, 0), std::atan2(c2, 0)),
                              atan2(rSimd_c0c1c2, setZero()));
    GMX_EXPECT_SIMD_REAL_NEAR(
            setSimdRealFrom3R(std::atan2(0, -c0), std::atan2(0, -c1), std::atan2(0, -c2)),
            atan2(setZero(), rSimd_m0m1m2));
    GMX_EXPECT_SIMD_REAL_NEAR(
            setSimdRealFrom3R(std::atan2(-c0, 0), std::atan2(-c1, 0), std::atan2(-c2, 0)),
            atan2(rSimd_m0m1m2, setZero()));
    // degenerate value (origin) should return 0.0. At least IBM xlc 13.1.5 gets the reference
    // value wrong (-nan) at -O3 optimization, so we compare to the correct value (0.0) instead.
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(0.0), atan2(setSimdRealFrom3R(0.0, 0.0, 0.0), setZero()));
}

/*! \brief Evaluate reference version of PME force correction. */
real refPmeForceCorrection(real x)
{
    if (x != 0)
    {
        real y = std::sqrt(x);
        return 2 * std::exp(-x) / (std::sqrt(M_PI) * x) - std::erf(static_cast<double>(y)) / (x * y);
    }
    else
    {
        return -4 / (3 * std::sqrt(M_PI));
    }
}

// The PME corrections will be added to ~1/r2, so absolute tolerance of EPS is fine.
TEST_F(SimdMathTest, pmeForceCorrection)
{
    // Pme correction relative accuracy only needs to be ~1e-6 accuracy single, 1e-10 double
    const std::int64_t ulpTol = (GMX_DOUBLE ? 5e-10 : 5e-6) / GMX_REAL_EPS;

    CompareSettings settings{ Range(0.15, 4), ulpTol, GMX_REAL_EPS, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(refPmeForceCorrection, pmeForceCorrection, settings);
}

/*! \brief Evaluate reference version of PME potential correction. */
real refPmePotentialCorrection(real x)
{
    real y = std::sqrt(x);
    return std::erf(static_cast<double>(y)) / y;
}

// The PME corrections will be added to ~1/r, so absolute tolerance of EPS is fine.
TEST_F(SimdMathTest, pmePotentialCorrection)
{
    // Pme correction relative accuracy only needs to be ~1e-6 accuracy single, 1e-10 double
    const std::int64_t ulpTol = (GMX_DOUBLE ? 5e-10 : 5e-6) / GMX_REAL_EPS;

    CompareSettings settings{ Range(0.15, 4), ulpTol, GMX_REAL_EPS, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(refPmePotentialCorrection, pmePotentialCorrection, settings);
}

// Functions that only target single accuracy, even for double SIMD data

TEST_F(SimdMathTest, invsqrtSingleAccuracy)
{
    // Here we always use float limits, since the lookup is not defined for numbers that
    // cannot be represented in single precision.
    const real low  = std::numeric_limits<float>::min();
    const real high = std::numeric_limits<float>::max();
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTolSingleAccuracy(ulpTol_);

    CompareSettings settings{ Range(low, high), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(refInvsqrt, invsqrtSingleAccuracy, settings);
}

/*! \brief Function wrapper to return first result when testing \ref invsqrtPairSingleAccuracy */
SimdReal gmx_simdcall tst_invsqrt_SingleAccuracy_pair0(SimdReal x)
{
    SimdReal r0, r1;
    invsqrtPairSingleAccuracy(x, x, &r0, &r1);
    return r0;
}

/*! \brief Function wrapper to return second result when testing \ref invsqrtPairSingleAccuracy */
SimdReal gmx_simdcall tst_invsqrt_SingleAccuracy_pair1(SimdReal x)
{
    SimdReal r0, r1;
    invsqrtPairSingleAccuracy(x, x, &r0, &r1);
    return r1;
}

TEST_F(SimdMathTest, invsqrtPairSingleAccuracy)
{
    // Float limits since lookup is always performed in single
    const real low  = std::numeric_limits<float>::min();
    const real high = std::numeric_limits<float>::max();
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTolSingleAccuracy(ulpTol_);

    CompareSettings settings{ Range(low, high), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(refInvsqrt, tst_invsqrt_SingleAccuracy_pair0, settings);
    GMX_EXPECT_SIMD_FUNC_NEAR(refInvsqrt, tst_invsqrt_SingleAccuracy_pair1, settings);
}

TEST_F(SimdMathTest, sqrtSingleAccuracy)
{
    // Since the first lookup step is sometimes performed in single precision,
    // our SIMD sqrt can only handle single-precision input values, even when
    // compiled in double precision - thus we use single precision limits here.

    // Scale lowest value by 1+eps, since we will do some arithmetics on this value
    const real low = std::numeric_limits<float>::min() * (1.0 + std::numeric_limits<float>::epsilon());
    const real      high = std::numeric_limits<float>::max();
    CompareSettings settings;

    // Increase the allowed error by the difference between the actual precision and single
    setUlpTolSingleAccuracy(ulpTol_);

    // First test that 0.0 and a few other values works
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(0, std::sqrt(c0), std::sqrt(c1)),
                              sqrtSingleAccuracy(setSimdRealFrom3R(0, c0, c1)));

#        if GMX_DOUBLE
    // As mentioned above, we cannot guarantee that very small double precision
    // input values (below std::numeric_limits<float>::min()) are handled correctly,
    // so our implementation will clamp it to zero. In this range we allow either
    // the correct value or zero, but it's important that it does not result in NaN or Inf values.
    //
    // This test range must not be called for single precision, since if we try to divide
    // the interval (0.0, low( in npoints we will try to multiply by factors so small that
    // they end up being flushed to zero, and the loop would never end.
    settings = { Range(0.0, low), ulpTol_, absTol_, MatchRule::ReferenceOrZero };
    GMX_EXPECT_SIMD_FUNC_NEAR(refSqrt, sqrtSingleAccuracy, settings);
#        endif

    settings = { Range(low, high), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(refSqrt, sqrtSingleAccuracy, settings);
}

TEST_F(SimdMathTest, sqrtSingleAccuracyUnsafe)
{
    // Test the full range, but stick to float limits since lookup is done in single.
    const real low  = std::numeric_limits<float>::min();
    const real high = std::numeric_limits<float>::max();

    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTolSingleAccuracy(ulpTol_);

    CompareSettings settings{ Range(low, high), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(refSqrt, sqrtSingleAccuracy<MathOptimization::Unsafe>, settings);
}

TEST_F(SimdMathTest, invSingleAccuracy)
{
    // Since the first lookup step is sometimes performed in single precision,
    // our SIMD 1/x can only handle single-precision input values, even when
    // compiled in double precision.

    // Relevant threshold points
    const real minSafeFloat = std::numeric_limits<float>::min()
                              * 10; // X value guaranteed not to result in Inf intermediates for 1/x calc.
    const real maxSafeFloat = std::numeric_limits<float>::max()
                              * 0.1; // X value guaranteed not to result in DTZ intermediates for 1/x calc.
    // Scale highest value by 1-eps, since we will do some arithmetics on this value
    const real maxFloat =
            std::numeric_limits<float>::max() * (1.0 - std::numeric_limits<float>::epsilon());
    CompareSettings settings;

    // Increase the allowed error by the difference between the actual precision and single
    setUlpTolSingleAccuracy(ulpTol_);

    // Danger zone where intermediates might be flushed to zero and produce 1/x==0.0
    settings = { Range(-maxFloat, -maxSafeFloat), ulpTol_, absTol_, MatchRule::ReferenceOrZero };
    GMX_EXPECT_SIMD_FUNC_NEAR(refInv, inv, settings);

    // Normal checks for x < 0
    settings = { Range(-maxSafeFloat, -minSafeFloat), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(refInv, inv, settings);

    // We do not care about the small range -minSafeFloat < x < +minSafeFloat where the result can be +/- Inf, since we don't require strict IEEE754.

    // Normal checks for x > 0
    settings = { Range(minSafeFloat, maxSafeFloat), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(refInv, inv, settings);

    // Danger zone where intermediates might be flushed to zero and produce 1/x==0.0
    settings = { Range(maxSafeFloat, maxFloat), ulpTol_, absTol_, MatchRule::ReferenceOrZero };
    GMX_EXPECT_SIMD_FUNC_NEAR(refInv, inv, settings);
}

TEST_F(SimdMathTest, cbrtSingleAccuracy)
{
    const real low  = -std::numeric_limits<real>::max();
    const real high = std::numeric_limits<real>::max();

    // Increase the allowed error by the difference between the actual precision and single
    setUlpTolSingleAccuracy(ulpTol_);

    CompareSettings settings{ Range(low, high), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::cbrt, cbrtSingleAccuracy, settings);
}

TEST_F(SimdMathTest, invcbrtSingleAccuracy)
{
    // Increase the allowed error by the difference between the actual precision and single
    setUlpTolSingleAccuracy(ulpTol_);

    // Negative values first
    real low  = -std::numeric_limits<real>::max();
    real high = -std::numeric_limits<real>::min();

    CompareSettings settings{ Range(low, high), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(refInvCbrt, invcbrtSingleAccuracy, settings);

    // Positive values
    low      = std::numeric_limits<real>::min();
    high     = std::numeric_limits<real>::max();
    settings = { Range(low, high), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(refInvCbrt, invcbrtSingleAccuracy, settings);
}

TEST_F(SimdMathTest, log2SingleAccuracy)
{
    const real low  = std::numeric_limits<real>::min();
    const real high = std::numeric_limits<real>::max();

    // Increase the allowed error by the difference between the actual precision and single
    setUlpTolSingleAccuracy(ulpTol_);

    CompareSettings settings{ Range(low, high), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::log2, log2SingleAccuracy, settings);
}

TEST_F(SimdMathTest, logSingleAccuracy)
{
    const real low  = std::numeric_limits<real>::min();
    const real high = std::numeric_limits<real>::max();

    // Increase the allowed error by the difference between the actual precision and single
    setUlpTolSingleAccuracy(ulpTol_);

    CompareSettings settings{ Range(low, high), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::log, logSingleAccuracy, settings);
}

TEST_F(SimdMathTest, exp2SingleAccuracy)
{
    // Relevant threshold points - float limits since we only target single accuracy
    constexpr real lowestReal = -std::numeric_limits<real>::max();
    constexpr real lowestRealThatProducesNormal =
            std::numeric_limits<real>::min_exponent
            - 1; // adding the significant corresponds to one more unit in exponent
    constexpr real lowestRealThatProducesDenormal =
            lowestRealThatProducesNormal
            - std::numeric_limits<real>::digits; // digits refer to bits in significand, so 24/53 for float/double
    constexpr real highestRealThatProducesNormal =
            std::numeric_limits<real>::max_exponent
            - 1; // adding the significant corresponds to one more unit in exponent
    CompareSettings settings;

    // Increase the allowed error by the difference between the actual precision and single
    setUlpTolSingleAccuracy(ulpTol_);

    // Below subnormal range all results should be zero
    settings = { Range(lowestReal, lowestRealThatProducesDenormal), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp2, exp2SingleAccuracy, settings);

    // Subnormal range, require matching, but DTZ is fine
    settings = { Range(lowestRealThatProducesDenormal, lowestRealThatProducesNormal),
                 ulpTol_,
                 absTol_,
                 MatchRule::Dtz };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp2, exp2SingleAccuracy, settings);

    // Normal range, standard result expected
    settings = { Range(lowestRealThatProducesNormal, highestRealThatProducesNormal),
                 ulpTol_,
                 absTol_,
                 MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp2, exp2SingleAccuracy, settings);
}

TEST_F(SimdMathTest, exp2SingleAccuracyUnsafe)
{
    // The unsafe version is only defined in the normal range
    constexpr real lowestRealThatProducesNormal =
            std::numeric_limits<real>::min_exponent
            - 1; // adding the significant corresponds to one more unit in exponent
    constexpr real highestRealThatProducesNormal =
            std::numeric_limits<real>::max_exponent
            - 1; // adding the significant corresponds to one more unit in exponent

    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTolSingleAccuracy(ulpTol_);

    CompareSettings settings{ Range(lowestRealThatProducesNormal, highestRealThatProducesNormal),
                              ulpTol_,
                              absTol_,
                              MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp2, exp2SingleAccuracy<MathOptimization::Unsafe>, settings);
}

TEST_F(SimdMathTest, expSingleAccuracy)
{
    // See threshold point comments in normal exp() test
    const real lowestReal = -std::numeric_limits<real>::max();
    // In theory the smallest value should be (min_exponent-1)*log(2), but rounding after the multiplication will cause this
    // value to be a single ulp too low. This might cause failed tests on CPUs that use different DTZ modes for SIMD vs.
    // non-SIMD arithmetics (e.g. ARM v7), so multiply by (1.0-eps) to increase it by a single ulp.
    const real lowestRealThatProducesNormal = (std::numeric_limits<real>::min_exponent - 1)
                                              * std::log(2.0)
                                              * (1.0 - std::numeric_limits<real>::epsilon());
    const real lowestRealThatProducesDenormal =
            lowestRealThatProducesNormal - std::numeric_limits<real>::digits * std::log(2.0);
    const real highestRealThatProducesNormal =
            (std::numeric_limits<real>::max_exponent - 1) * std::log(2.0);
    CompareSettings settings;

    // Increase the allowed error by the difference between the actual precision and single
    setUlpTolSingleAccuracy(ulpTol_);

    // Below subnormal range all results should be zero (so, match the reference)
    settings = { Range(lowestReal, lowestRealThatProducesDenormal), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp, expSingleAccuracy, settings);

    // Subnormal range, require matching, but DTZ is fine
    settings = { Range(lowestRealThatProducesDenormal, lowestRealThatProducesNormal),
                 ulpTol_,
                 absTol_,
                 MatchRule::Dtz };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp, expSingleAccuracy, settings);

    // Normal range, standard result expected
    settings = { Range(lowestRealThatProducesNormal, highestRealThatProducesNormal),
                 ulpTol_,
                 absTol_,
                 MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp, expSingleAccuracy, settings);
}

TEST_F(SimdMathTest, expSingleAccuracyUnsafe)
{
    // See test of exp() for comments about test ranges
    const real lowestRealThatProducesNormal = (std::numeric_limits<real>::min_exponent - 1)
                                              * std::log(2.0)
                                              * (1 - std::numeric_limits<real>::epsilon());
    const real highestRealThatProducesNormal =
            (std::numeric_limits<real>::max_exponent - 1) * std::log(2.0);

    // Increase the allowed error by the difference between the actual precision and single
    setUlpTolSingleAccuracy(ulpTol_);

    CompareSettings settings{ Range(lowestRealThatProducesNormal, highestRealThatProducesNormal),
                              ulpTol_,
                              absTol_,
                              MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::exp, expSingleAccuracy<MathOptimization::Unsafe>, settings);
}

TEST_F(SimdMathTest, powSingleAccuracy)
{
    // We already test the log2/exp2 components of pow() extensively above, and it's a very
    // simple single-line function, so here we just test a handful of values to catch typos
    // and then some special values.

    // Increase the allowed error by the difference between the actual precision and single
    setUlpTolSingleAccuracy(ulpTol_);

    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::pow(c0, c3), std::pow(c1, c4), std::pow(c2, c5)),
                              powSingleAccuracy(rSimd_c0c1c2, rSimd_c3c4c5));

    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::pow(c0, -c3), std::pow(c1, -c0), std::pow(c2, -c4)),
                              powSingleAccuracy(rSimd_c0c1c2, rSimd_m3m0m4));

    // 0^0 = 1 , 0^c1=0, -c1^0=1
    GMX_EXPECT_SIMD_REAL_NEAR(
            setSimdRealFrom3R(1.0, 0.0, 1.0),
            powSingleAccuracy(setSimdRealFrom3R(0, 0.0, -c1), setSimdRealFrom3R(0.0, c1, 0.0)));
}

TEST_F(SimdMathTest, powSingleAccuracyUnsafe)
{
    // We already test the log2/exp2 components of pow() extensively above, and it's a very
    // simple single-line function, so here we just test a handful of values to catch typos
    // and then some special values.

    // Increase the allowed error by the difference between the actual precision and single
    setUlpTolSingleAccuracy(ulpTol_);

    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::pow(c0, c3), std::pow(c1, c4), std::pow(c2, c5)),
                              powSingleAccuracy<MathOptimization::Unsafe>(rSimd_c0c1c2, rSimd_c3c4c5));

    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::pow(c0, -c3), std::pow(c1, -c0), std::pow(c2, -c4)),
                              powSingleAccuracy<MathOptimization::Unsafe>(rSimd_c0c1c2, rSimd_m3m0m4));
}

TEST_F(SimdMathTest, erfSingleAccuracy)
{
    // Increase the allowed error by the difference between the actual precision and single
    setUlpTolSingleAccuracy(ulpTol_);

    CompareSettings settings{ Range(-9, 9), ulpTol_, GMX_REAL_MIN, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(refErf, erfSingleAccuracy, settings);
}

TEST_F(SimdMathTest, erfcSingleAccuracy)
{
    // Increase the allowed error by the difference between the actual precision and single
    setUlpTolSingleAccuracy(ulpTol_);

    // Our erfc algorithm has 4 ulp accuracy, so relax tolerance a bit
    CompareSettings settings{ Range(-9, 9), 4 * ulpTol_, GMX_REAL_MIN, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(refErfc, erfcSingleAccuracy, settings);
}


TEST_F(SimdMathTest, sinSingleAccuracy)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTolSingleAccuracy(ulpTol_);

    CompareSettings settings{ Range(-8 * M_PI, 8 * M_PI), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::sin, sinSingleAccuracy, settings);

    // Range reduction leads to accuracy loss, so we might want higher tolerance here
    settings = { Range(-10000, 10000), 2 * ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::sin, sinSingleAccuracy, settings);
}

TEST_F(SimdMathTest, cosSingleAccuracy)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTolSingleAccuracy(ulpTol_);

    CompareSettings settings{ Range(-8 * M_PI, 8 * M_PI), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::cos, cosSingleAccuracy, settings);

    // Range reduction leads to accuracy loss, so we might want higher tolerance here
    settings = { Range(-10000, 10000), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::cos, cosSingleAccuracy, settings);
}

TEST_F(SimdMathTest, tanSingleAccuracy)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTolSingleAccuracy(ulpTol_);

    // Tan(x) is a little sensitive due to the division in the algorithm.
    // Rather than using lots of extra FP operations, we accept the algorithm
    // presently only achieves a ~3 ulp error and use the medium tolerance.
    CompareSettings settings{ Range(-8 * M_PI, 8 * M_PI), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::tan, tanSingleAccuracy, settings);

    // Range reduction leads to accuracy loss, so we might want higher tolerance here
    settings = { Range(-10000, 10000), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::tan, tanSingleAccuracy, settings);
}

TEST_F(SimdMathTest, asinSingleAccuracy)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTolSingleAccuracy(ulpTol_);

    // Our present asin(x) algorithm achieves 2-3 ulp accuracy
    CompareSettings settings{ Range(-1, 1), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::asin, asinSingleAccuracy, settings);
}

TEST_F(SimdMathTest, acosSingleAccuracy)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTolSingleAccuracy(ulpTol_);

    // Our present acos(x) algorithm achieves 2-3 ulp accuracy
    CompareSettings settings{ Range(-1, 1), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::acos, acosSingleAccuracy, settings);
}

TEST_F(SimdMathTest, atanSingleAccuracy)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTolSingleAccuracy(ulpTol_);

    // Our present atan(x) algorithm achieves 1 ulp accuracy
    CompareSettings settings{ Range(-10000, 10000), ulpTol_, absTol_, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(std::atan, atanSingleAccuracy, settings);
}

TEST_F(SimdMathTest, atan2SingleAccuracy)
{
    /* Increase the allowed error by the difference between the actual precision and single */
    setUlpTolSingleAccuracy(ulpTol_);

    // test each quadrant
    GMX_EXPECT_SIMD_REAL_NEAR(
            setSimdRealFrom3R(std::atan2(c0, c3), std::atan2(c1, c4), std::atan2(c2, c5)),
            atan2SingleAccuracy(rSimd_c0c1c2, rSimd_c3c4c5));
    GMX_EXPECT_SIMD_REAL_NEAR(
            setSimdRealFrom3R(std::atan2(-c0, c3), std::atan2(-c1, c4), std::atan2(-c2, c5)),
            atan2SingleAccuracy(rSimd_m0m1m2, rSimd_c3c4c5));
    GMX_EXPECT_SIMD_REAL_NEAR(
            setSimdRealFrom3R(std::atan2(-c0, -c3), std::atan2(-c1, -c0), std::atan2(-c2, -c4)),
            atan2SingleAccuracy(rSimd_m0m1m2, rSimd_m3m0m4));
    GMX_EXPECT_SIMD_REAL_NEAR(
            setSimdRealFrom3R(std::atan2(c0, -c3), std::atan2(c1, -c0), std::atan2(c2, -c4)),
            atan2SingleAccuracy(rSimd_c0c1c2, rSimd_m3m0m4));
    // cases important for calculating angles
    // values on coordinate axes
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::atan2(0, c0), std::atan2(0, c1), std::atan2(0, c2)),
                              atan2SingleAccuracy(setZero(), rSimd_c0c1c2));
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom3R(std::atan2(c0, 0), std::atan2(c1, 0), std::atan2(c2, 0)),
                              atan2SingleAccuracy(rSimd_c0c1c2, setZero()));
    GMX_EXPECT_SIMD_REAL_NEAR(
            setSimdRealFrom3R(std::atan2(0, -c0), std::atan2(0, -c1), std::atan2(0, -c2)),
            atan2SingleAccuracy(setZero(), rSimd_m0m1m2));
    GMX_EXPECT_SIMD_REAL_NEAR(
            setSimdRealFrom3R(std::atan2(-c0, 0), std::atan2(-c1, 0), std::atan2(-c2, 0)),
            atan2SingleAccuracy(rSimd_m0m1m2, setZero()));

    // degenerate value (origin) should return 0.0. At least IBM xlc 13.1.5 gets the reference
    // value wrong (-nan) at -O3 optimization, so we compare to the correct value (0.0) instead.
    GMX_EXPECT_SIMD_REAL_NEAR(setSimdRealFrom1R(0.0),
                              atan2SingleAccuracy(setSimdRealFrom3R(0.0, 0.0, 0.0), setZero()));
}

TEST_F(SimdMathTest, pmeForceCorrectionSingleAccuracy)
{
    // The PME corrections will be added to ~1/r2, so absolute tolerance of EPS is fine.
    // Pme correction only needs to be ~1e-6 accuracy single.
    // Then increase the allowed error by the difference between the actual precision and single.
    setUlpTolSingleAccuracy(std::int64_t(5e-6 / GMX_FLOAT_EPS));

    CompareSettings settings{ Range(0.15, 4), ulpTol_, GMX_FLOAT_EPS, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(refPmeForceCorrection, pmeForceCorrectionSingleAccuracy, settings);
}

TEST_F(SimdMathTest, pmePotentialCorrectionSingleAccuracy)
{
    // The PME corrections will be added to ~1/r, so absolute tolerance of EPS is fine.
    // Pme correction only needs to be ~1e-6 accuracy single.
    // Then increase the allowed error by the difference between the actual precision and single.
    setUlpTolSingleAccuracy(std::int64_t(5e-6 / GMX_FLOAT_EPS));

    CompareSettings settings{ Range(0.15, 4), ulpTol_, GMX_FLOAT_EPS, MatchRule::Normal };
    GMX_EXPECT_SIMD_FUNC_NEAR(refPmePotentialCorrection, pmePotentialCorrectionSingleAccuracy, settings);
}

} // namespace

#    endif // GMX_SIMD_HAVE_REAL

/*! \} */
/*! \endcond */

} // namespace test
} // namespace gmx

#endif // GMX_SIMD
